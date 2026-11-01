//! Parts that fold into themselves. Every flat cell is cut into prisms (a
//! triangle through the sheet, in thin slices across a bend), each taken to
//! 3D by its zone's fold map. Two prisms of different zones that overlap
//! deeper than `TOUCH` clash; the overlap is mapped back to the flat, on both
//! features, so it can be shown where the part runs into itself.

use crate::flat::{Arc, Cell, Zone, P2};
use crate::fold::{Model, P3};
use crate::spec::PartSpec;
use crate::SheetError;
use schemars::JsonSchema;
use serde::Serialize;
use std::cell::OnceCell;
use std::collections::BTreeMap;

/// Overlap up to this deep is contact, not a clash (mm).
const TOUCH: f64 = 0.01;
/// Largest gap between an arc of the outline and its chords (mm).
const CHORD: f64 = TOUCH / 4.0;
/// Largest gap between a bend's faces and the slices across it (mm).
const SLICE: f64 = TOUCH;

/// Two features whose material overlaps once folded.
#[derive(Clone, Debug, PartialEq, Serialize, JsonSchema)]
pub struct Clash {
    /// Flange, contour, fold or hem ids, "base" for the base; the same id
    /// twice when a feature folds into itself.
    pub features: [String; 2],
    /// How deep they overlap, at the deepest (mm).
    pub depth: f64,
    /// Where, on the flat pattern: convex pieces on each feature.
    pub areas: [Vec<Vec<P2>>; 2],
}

impl std::fmt::Display for Clash {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [a, b] = &self.features;
        let at = |k: usize| {
            let c = centroid(&self.areas[k]);
            format!("({:.1}, {:.1})", c[0], c[1])
        };
        if a == b {
            write!(f, "{a} folds into itself")?;
        } else {
            write!(f, "{a} and {b} fold into each other")?;
        }
        write!(
            f,
            ", {:.2} mm deep: around {} and {} on the flat pattern",
            self.depth,
            at(0),
            at(1)
        )
    }
}

/// Every pair of features that folds into each other, deepest first.
pub fn clashes(spec: &PartSpec) -> Result<Vec<Clash>, SheetError> {
    let spec = &*spec.resolved()?;
    let layout = crate::flat::layout(spec)?;
    Ok(Model::new(spec, &layout)
        .overlaps()
        .into_iter()
        .map(|o| o.clash)
        .collect())
}

/// The folded part as loose triangles (the sheet's two faces, no walls),
/// built even where it folds into itself, and its clashing areas on both
/// faces. For showing a part the kernel can't build.
pub struct Preview {
    pub sheet: Vec<[P3; 3]>,
    pub clashes: Vec<[P3; 3]>,
}

pub fn preview(spec: &PartSpec) -> Result<Preview, SheetError> {
    let spec = &*spec.resolved()?;
    let layout = crate::flat::layout(spec)?;
    let model = Model::new(spec, &layout);
    let mut sheet = Vec::new();
    for cell in &layout.cells {
        for tri in model.triangles(cell) {
            for z in [0.0, model.t] {
                sheet.push(tri.map(|p| model.map(cell.zone, p, z)));
            }
        }
    }
    let mut clashes = Vec::new();
    for o in model.overlaps() {
        for (zone, poly) in o.zoned {
            for k in 1..poly.len().saturating_sub(1) {
                for z in [0.0, model.t] {
                    clashes.push([poly[0], poly[k], poly[k + 1]].map(|p| model.map(zone, p, z)));
                }
            }
        }
    }
    Ok(Preview { sheet, clashes })
}

struct Overlap {
    clash: Clash,
    /// Each area with the zone it lies in.
    zoned: Vec<(Zone, Vec<P2>)>,
}

/// A triangle of a cell through the sheet, in 3D: points 0..3 on the z = 0
/// face, 3..6 above them on z = T.
struct Prism {
    zone: Zone,
    /// The triangle in the flat.
    flat: [P2; 3],
    pts: [P3; 6],
    lo: P3,
    hi: P3,
    /// How far the hull's faces may lie off the sheet's faces.
    slack: f64,
    /// Worked out when first tested.
    hull: OnceCell<Hull>,
}

/// A prism's convex hull.
struct Hull {
    /// Outward unit normals and offsets.
    planes: Vec<(P3, f64)>,
    /// Edge directions, unit, no two parallel.
    edges: Vec<P3>,
}

impl Model<'_> {
    /// Refuses a part that folds into itself, naming the deepest clash.
    pub(crate) fn refuse_clashes(&self) -> Result<(), SheetError> {
        match self.overlaps().into_iter().next() {
            Some(o) => Err(SheetError::Invalid(o.clash.to_string())),
            None => Ok(()),
        }
    }

    fn feature(&self, zone: Zone) -> String {
        match zone {
            Zone::Base => "base".into(),
            Zone::Bend(i) | Zone::Flange(i) => self.layout.folds[i].flange.clone(),
        }
    }

    fn overlaps(&self) -> Vec<Overlap> {
        let mut prisms: Vec<Prism> = Vec::new();
        for cell in &self.layout.cells {
            for tri in self.triangles(cell) {
                let pts: [P3; 6] = std::array::from_fn(|k| {
                    self.map(cell.zone, tri[k % 3], if k < 3 { 0.0 } else { self.t })
                });
                prisms.push(Prism::new(cell.zone, tri, pts));
            }
        }
        prisms.sort_by(|a, b| a.lo[0].total_cmp(&b.lo[0]));
        // A bend can't reach the two zones it joins short of a full half turn.
        let joined = |a: Zone, b: Zone| match (a, b) {
            (Zone::Bend(i), other) | (other, Zone::Bend(i)) => {
                let f = &self.layout.folds[i];
                let parent = f.parent.map_or(Zone::Base, Zone::Flange);
                f.angle <= std::f64::consts::PI + 1e-9
                    && (other == Zone::Flange(i) || other == parent)
            }
            _ => false,
        };

        // Per pair of zones: the deepest clash, and the areas on each of
        // every overlap (also those too shallow to count on their own, so
        // the areas don't miss pieces along a triangle's edge).
        type Found = (f64, Vec<Vec<P2>>, Vec<Vec<P2>>);
        let mut found: BTreeMap<(usize, usize), (Zone, Zone, Found)> = BTreeMap::new();
        let index = |z: Zone| match z {
            Zone::Base => 0,
            Zone::Bend(i) => 2 * i + 1,
            Zone::Flange(i) => 2 * i + 2,
        };
        for (i, p) in prisms.iter().enumerate() {
            for q in &prisms[i + 1..] {
                if q.lo[0] > p.hi[0] - TOUCH {
                    break;
                }
                if p.zone == q.zone
                    || (1..3).any(|j| q.lo[j] > p.hi[j] - TOUCH || p.lo[j] > q.hi[j] - TOUCH)
                    || joined(p.zone, q.zone)
                    || touch(&p.flat, &q.flat)
                {
                    continue;
                }
                let Some(depth) = depth(p, q, 0.0) else {
                    continue;
                };
                let clash = if depth > TOUCH + p.slack + q.slack {
                    depth
                } else {
                    0.0
                };
                let (p, q) = if index(p.zone) < index(q.zone) {
                    (p, q)
                } else {
                    (q, p)
                };
                let corners = intersection(p, q);
                let area = |prism: &Prism| {
                    hull(corners.iter().map(|&x| self.unmap(prism.zone, x)).collect())
                };
                let e = found.entry((index(p.zone), index(q.zone))).or_insert((
                    p.zone,
                    q.zone,
                    (0.0, Vec::new(), Vec::new()),
                ));
                e.2 .0 = e.2 .0.max(clash);
                e.2 .1.push(area(p));
                e.2 .2.push(area(q));
            }
        }

        // Per pair of features.
        let mut out: Vec<Overlap> = Vec::new();
        for (za, zb, (depth, aa, ab)) in found.into_values().filter(|f| f.2 .0 > 0.0) {
            let features = [self.feature(za), self.feature(zb)];
            let zoned = aa
                .iter()
                .map(|a| (za, a.clone()))
                .chain(ab.iter().map(|a| (zb, a.clone())));
            match out.iter_mut().find(|o| o.clash.features == features) {
                Some(o) => {
                    o.clash.depth = o.clash.depth.max(depth);
                    o.clash.areas[0].extend(aa.iter().cloned());
                    o.clash.areas[1].extend(ab.iter().cloned());
                    o.zoned.extend(zoned);
                }
                None => out.push(Overlap {
                    zoned: zoned.collect(),
                    clash: Clash {
                        features,
                        depth,
                        areas: [aa, ab],
                    },
                }),
            }
        }
        out.sort_by(|a, b| b.clash.depth.total_cmp(&a.clash.depth));
        out
    }

    /// A cell in triangles; in a bend, none spans more than a slice.
    fn triangles(&self, cell: &Cell) -> Vec<[P2; 3]> {
        let mut data: Vec<f64> = Vec::new();
        let mut starts: Vec<usize> = Vec::new();
        let mut push = |ring: Vec<P2>, hole: bool| {
            if hole {
                starts.push(data.len() / 2);
            }
            data.extend(ring.iter().flat_map(|p| *p));
        };
        push(polyline(&cell.polygon, &cell.arcs), false);
        for r in &cell.inner {
            push(polyline(&r.points, &r.arcs), true);
        }
        for h in &cell.holes {
            let (c, r) = (h.center, h.diameter / 2.0);
            let n = pieces(r, std::f64::consts::TAU);
            let ring = (0..n)
                .map(|k| {
                    let (s, co) = (std::f64::consts::TAU * k as f64 / n as f64).sin_cos();
                    [c[0] + r * co, c[1] + r * s]
                })
                .collect();
            push(ring, true);
        }
        let point = |k: usize| [data[2 * k], data[2 * k + 1]];
        let mut tris: Vec<[P2; 3]> = earcutr::earcut(&data, &starts, 2)
            .as_chunks::<3>()
            .0
            .iter()
            .map(|t| t.map(point))
            .collect();
        let Zone::Bend(i) = cell.zone else {
            return tris;
        };
        // A slice spans at most this much flat length along n.
        let f = &self.folds[i].fold;
        let outer = f.radius + self.t;
        let step = 2.0 * (1.0 - (SLICE / outer).min(1.0)).acos() * (f.radius + f.k_factor * self.t);
        let s = |p: P2| (p[0] - f.edge_point[0]) * f.n[0] + (p[1] - f.edge_point[1]) * f.n[1];
        let mut out = Vec::new();
        while let Some(t) = tris.pop() {
            let (k, span) = (0..3)
                .map(|k| (k, (s(t[k]) - s(t[(k + 1) % 3])).abs()))
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();
            if span <= step {
                out.push(t);
                continue;
            }
            let (a, b, c) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
            let m = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
            tris.push([a, m, c]);
            tris.push([m, b, c]);
        }
        out
    }
}

impl Prism {
    fn new(zone: Zone, flat: [P2; 3], pts: [P3; 6]) -> Prism {
        // A slice's faces are chords of the bend's. Its side faces twist
        // too, so its hull bulges past them; only into a neighbour that
        // shares them, which never counts (`touch`).
        let slack = match zone {
            Zone::Bend(_) => SLICE,
            _ => CHORD,
        };
        let (mut lo, mut hi) = (pts[0], pts[0]);
        for p in &pts {
            for j in 0..3 {
                lo[j] = lo[j].min(p[j]);
                hi[j] = hi[j].max(p[j]);
            }
        }
        Prism {
            zone,
            flat,
            pts,
            lo,
            hi,
            slack,
            hull: OnceCell::new(),
        }
    }

    fn hull(&self) -> &Hull {
        self.hull.get_or_init(|| Hull::new(&self.pts))
    }

    fn span(&self, axis: P3) -> (f64, f64) {
        self.pts
            .iter()
            .map(|&p| dot(axis, p))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), x| {
                (a.min(x), b.max(x))
            })
    }
}

impl Hull {
    fn new(pts: &[P3; 6]) -> Hull {
        let scale = pts.iter().flatten().fold(1.0f64, |m, x| m.max(x.abs()));
        let eps = 1e-9 * scale;
        let mut planes: Vec<(P3, f64)> = Vec::new();
        for i in 0..6 {
            for j in i + 1..6 {
                for k in j + 1..6 {
                    let n = cross(sub(pts[j], pts[i]), sub(pts[k], pts[i]));
                    let l = norm(n);
                    if l < 1e-12 * scale * scale {
                        continue;
                    }
                    let n = scale3(n, 1.0 / l);
                    let d = dot(n, pts[i]);
                    let plane = if pts.iter().all(|&p| dot(n, p) - d <= eps) {
                        (n, d)
                    } else if pts.iter().all(|&p| dot(n, p) - d >= -eps) {
                        (scale3(n, -1.0), -d)
                    } else {
                        continue;
                    };
                    if !planes.iter().any(|q| dot(q.0, plane.0) > 1.0 - 1e-12) {
                        planes.push(plane);
                    }
                }
            }
        }
        let mut edges: Vec<P3> = Vec::new();
        for i in 0..6 {
            for j in i + 1..6 {
                let on = |p: &&(P3, f64)| {
                    (dot(p.0, pts[i]) - p.1).abs() < eps && (dot(p.0, pts[j]) - p.1).abs() < eps
                };
                if planes.iter().filter(on).count() < 2 {
                    continue;
                }
                let d = sub(pts[j], pts[i]);
                let l = norm(d);
                if l < eps {
                    continue;
                }
                let d = scale3(d, 1.0 / l);
                if !edges.iter().any(|e| dot(*e, d).abs() > 1.0 - 1e-12) {
                    edges.push(d);
                }
            }
        }
        Hull { planes, edges }
    }
}

/// How deep two convex prisms overlap (separating axes), if deeper than `limit`.
fn depth(p: &Prism, q: &Prism, limit: f64) -> Option<f64> {
    let (hp, hq) = (p.hull(), q.hull());
    if hp.planes.len() < 4 || hq.planes.len() < 4 {
        // A sliver with no volume.
        return None;
    }
    let mut best = f64::INFINITY;
    let mut along = |axis: P3| {
        let (a0, a1) = p.span(axis);
        let (b0, b1) = q.span(axis);
        best = best.min((a1 - b0).min(b1 - a0));
        best > limit
    };
    for &(n, _) in hp.planes.iter().chain(&hq.planes) {
        if !along(n) {
            return None;
        }
    }
    for &e in &hp.edges {
        for &f in &hq.edges {
            let c = cross(e, f);
            let l = norm(c);
            if l > 1e-9 && !along(scale3(c, 1.0 / l)) {
                return None;
            }
        }
    }
    Some(best)
}

/// Whether two flat triangles meet: neighbours across a cell's edge, whose
/// prisms share a face and so always overlap a little.
fn touch(a: &[P2; 3], b: &[P2; 3]) -> bool {
    let apart = |t: &[P2; 3]| {
        (0..3).any(|k| {
            let (p, q) = (t[k], t[(k + 1) % 3]);
            let n = [q[1] - p[1], p[0] - q[0]];
            let l = (n[0] * n[0] + n[1] * n[1]).sqrt();
            if l < 1e-12 {
                return false;
            }
            let f = |x: &P2| (n[0] * x[0] + n[1] * x[1]) / l;
            let (a0, a1) = a
                .iter()
                .map(f)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                    (lo.min(x), hi.max(x))
                });
            let (b0, b1) = b
                .iter()
                .map(f)
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                    (lo.min(x), hi.max(x))
                });
            b0 > a1 + 1e-6 || a0 > b1 + 1e-6
        })
    };
    !apart(a) && !apart(b)
}

/// The corners of two convex prisms' common part.
fn intersection(p: &Prism, q: &Prism) -> Vec<P3> {
    let planes: Vec<(P3, f64)> = p
        .hull()
        .planes
        .iter()
        .chain(&q.hull().planes)
        .copied()
        .collect();
    let n = planes.len();
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let (a, b, c) = (planes[i], planes[j], planes[k]);
                let bc = cross(b.0, c.0);
                let det = dot(a.0, bc);
                if det.abs() < 1e-9 {
                    continue;
                }
                let x = scale3(
                    add(
                        add(scale3(bc, a.1), scale3(cross(c.0, a.0), b.1)),
                        scale3(cross(a.0, b.0), c.1),
                    ),
                    1.0 / det,
                );
                if planes.iter().all(|pl| dot(pl.0, x) - pl.1 < 1e-7) {
                    out.push(x);
                }
            }
        }
    }
    out
}

/// Counter-clockwise convex hull.
fn hull(pts: Vec<P2>) -> Vec<P2> {
    // On a grid, so points in line (an edge along x = 100 that comes out as
    // 99.99999999999999 and 100.00000000000001) sort and turn as such.
    let mut pts: Vec<P2> = pts
        .iter()
        .map(|p| p.map(|x| (x * 1e9).round() / 1e9))
        .collect();
    pts.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
    pts.dedup();
    if pts.len() < 3 {
        return pts;
    }
    let turn = |o: P2, a: P2, b: P2| (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    let mut h: Vec<P2> = Vec::new();
    for pass in 0..2 {
        let start = h.len();
        let order: Vec<P2> = if pass == 0 {
            pts.clone()
        } else {
            pts.iter().rev().copied().collect()
        };
        for p in order {
            while h.len() >= start + 2 && turn(h[h.len() - 2], h[h.len() - 1], p) <= 1e-12 {
                h.pop();
            }
            h.push(p);
        }
        h.pop();
    }
    h
}

/// Area-weighted centre of convex pieces.
fn centroid(areas: &[Vec<P2>]) -> P2 {
    let (mut sum, mut at) = (0.0, [0.0; 2]);
    for poly in areas {
        for k in 1..poly.len().saturating_sub(1) {
            let (a, b, c) = (poly[0], poly[k], poly[k + 1]);
            let w = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.0;
            sum += w;
            for j in 0..2 {
                at[j] += w * (a[j] + b[j] + c[j]) / 3.0;
            }
        }
    }
    if sum.abs() < 1e-12 {
        return areas.iter().flatten().next().copied().unwrap_or([0.0; 2]);
    }
    at.map(|x| x / sum)
}

/// A closed boundary as points, arcs in chords within `CHORD`.
fn polyline(points: &[P2], arcs: &[Option<Arc>]) -> Vec<P2> {
    let n = points.len();
    let mut out = Vec::new();
    for k in 0..n {
        let (p, q) = (points[k], points[(k + 1) % n]);
        out.push(p);
        if let Some(a) = arcs.get(k).copied().flatten() {
            let u = [p[0] - a.center[0], p[1] - a.center[1]];
            let r = (u[0] * u[0] + u[1] * u[1]).sqrt();
            let sweep = a.sweep(p, q);
            let m = pieces(r, sweep.abs());
            for i in 1..m {
                let (s, c) = (sweep * i as f64 / m as f64).sin_cos();
                out.push([
                    a.center[0] + u[0] * c - u[1] * s,
                    a.center[1] + u[0] * s + u[1] * c,
                ]);
            }
        }
    }
    out
}

/// Chords for an arc of radius r and the given sweep.
fn pieces(r: f64, sweep: f64) -> usize {
    let most = 2.0 * (1.0 - (CHORD / r).min(1.0)).acos();
    ((sweep / most.max(1e-3)).ceil() as usize).max(1)
}

fn add(a: P3, b: P3) -> P3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale3(a: P3, s: f64) -> P3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn folded(angle: f64) -> PartSpec {
        // A 100 x 60 blank folded along y = 25 and y = 35: two flaps of 25
        // over a strip of 10.
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "holes": [] }},
                 "folds": [{{ "id": "a", "line": [[0,25],[100,25]], "angle": {angle},
                              "location": "centerline", "side": "up" }},
                           {{ "id": "b", "line": [[100,35],[0,35]], "angle": {angle},
                              "location": "centerline", "side": "up" }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn flaps_folded_past_each_other_clash_where_they_cross() {
        assert_eq!(clashes(&folded(90.0)).unwrap(), vec![]);

        // At 120 degrees the flaps meet over the middle (y = 30), about 10
        // out from their bend lines, and cross on from there.
        let c = clashes(&folded(120.0)).unwrap();
        assert_eq!(c.len(), 1);
        assert_eq!(c[0].features, ["a".to_string(), "b".to_string()]);
        assert!(c[0].depth > 5.0, "{}", c[0].depth);
        // Each area on its own flap, about 10 + half the bend from its fold line.
        let ys = |k: usize| {
            c[0].areas[k]
                .iter()
                .flatten()
                .map(|p| p[1])
                .collect::<Vec<_>>()
        };
        assert!(ys(0).iter().all(|&y| y < 24.0), "{:?}", ys(0));
        assert!(ys(1).iter().all(|&y| y > 36.0), "{:?}", ys(1));
        let (ca, cb) = (centroid(&c[0].areas[0]), centroid(&c[0].areas[1]));
        assert!(
            (ca[0] - 50.0).abs() < 1.0 && (cb[0] - 50.0).abs() < 1.0,
            "{ca:?} {cb:?}"
        );
        assert!(
            (ca[1] - 12.5).abs() < 2.0 && (cb[1] - 47.5).abs() < 2.0,
            "{ca:?} {cb:?}"
        );

        match crate::fold::fold(&folded(120.0)) {
            Err(SheetError::Invalid(m)) => {
                assert!(m.starts_with("a and b fold into each other, 9."), "{m}")
            }
            other => panic!("{:?}", other.err()),
        }
        assert!(matches!(
            crate::fold::exact_mass(&folded(120.0)),
            Err(SheetError::Invalid(_))
        ));
    }

    #[test]
    fn a_hull_keeps_its_corner_on_an_edge_that_comes_out_a_hair_off_line() {
        let h = hull(vec![
            [99.99999999999999, 11.0],
            [100.00000000000001, 13.3],
            [100.0, 12.1],
            [100.0, 14.5],
            [65.5, 14.5],
            [50.1, 11.0],
        ]);
        assert_eq!(
            h,
            vec![[50.1, 11.0], [100.0, 11.0], [100.0, 14.5], [65.5, 14.5]]
        );
    }

    #[test]
    fn the_preview_shows_the_clash_on_both_faces_of_both_flaps() {
        let p = preview(&folded(120.0)).unwrap();
        assert!(!p.sheet.is_empty() && !p.clashes.is_empty());
        // The clash triangles lie where the flaps cross, over the middle
        // strip (the bends lift the flaps, so higher than 8.7).
        for t in &p.clashes {
            for v in t {
                assert!(
                    (v[1] - 30.0).abs() < 5.0 && v[2] > 8.0 && v[2] < 25.0,
                    "{v:?}"
                );
            }
        }
    }
}
