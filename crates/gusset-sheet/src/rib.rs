//! Gussets (`gussets`): a rib pressed across a bend into its inside. In a
//! section square to the bend axis the sheet runs straight across the corner,
//! from the base side's inner face to the flange's (the rib's inner face, a
//! plane parallel to the axis), with the groove behind it on the outside (a
//! parallel plane the rib's wall further out). Its side walls are planes
//! square to the axis: where the rib's side walls stand, the sheet fills the
//! corner from its outer face up to the rib's inner face; between them only
//! the wall is left, the groove open to the outside.
//!
//! Laid out as the plain part with the rib's footprint on the inner face (a
//! rectangle across the bend: the base side's part, the bend, the flange's)
//! and the groove's on the outer face cut out of the cells. The outer face
//! goes on under the rib's side walls and up to the groove, and with a wall
//! thinner than the sheet the inner face goes on over the groove's ends:
//! those pieces (`skins`) are cells on one level only. The fold model builds the rib's own faces on their edges: its inner
//! face, the groove's, and the side walls of each. All planes; their edges
//! lines and the bend's arcs.

use crate::flat::{path_area, Cell, Curve, FlatForm, Layout, Ring, Zone, P2};
use crate::fold::{Mass, Rigid, GAUSS};
use crate::spec::{Cutout, GussetShape, PartSpec};
use crate::SheetError;
use std::f64::consts::{FRAC_PI_2, PI};

const EPS: f64 = 1e-9;

/// A gusset, worked out.
#[derive(Clone, Debug)]
pub(crate) struct Formed {
    pub id: String,
    /// Its fold, and that fold's flat frame: (u, s) -> edge_point + u u + s n.
    pub fold: usize,
    edge_point: P2,
    u: P2,
    n: P2,
    /// The footprint on the inner face (u and s ranges) and the groove's on
    /// the outer face.
    pub u_in: [f64; 2],
    pub s_in: [f64; 2],
    pub u_groove: [f64; 2],
    pub s_groove: [f64; 2],
    /// The bend strip's s range.
    pub s_bend: [f64; 2],
    /// The inner face's level (1: the top face, for a bend up).
    pub inner: usize,
    /// The faces' pieces inside the cut-out, each with its level.
    pub skins: Vec<(Cell, usize)>,
    /// In the parent's section (s along n, z): the rib's inner face's
    /// normal, into its wall.
    normal: P2,
    /// Volume, area (the rib's faces and the skins) and the centre of
    /// volume in the parent's section (s, z) at u.
    volume: f64,
    area: f64,
    centre: [f64; 3],
}

impl Formed {
    /// The flat point at (u, s).
    pub fn flat(&self, u: f64, s: f64) -> P2 {
        [
            self.edge_point[0] + u * self.u[0] + s * self.n[0],
            self.edge_point[1] + u * self.u[1] + s * self.n[1],
        ]
    }

    /// The direction along the bend (flat u) and the rib's inner face's
    /// normal into the wall, in the parent's flat frame (x, y, z).
    pub fn directions(&self) -> ([f64; 3], [f64; 3]) {
        let (nu, n) = (self.normal, self.n);
        (
            [self.u[0], self.u[1], 0.0],
            [nu[0] * n[0], nu[0] * n[1], nu[1]],
        )
    }

    /// Mass of the rib (and the skins' area), placed by the parent's map.
    pub fn mass(&self, parent: &Rigid) -> Mass {
        let [u, s, z] = self.centre;
        let p = self.flat(u, s);
        Mass {
            volume: self.volume,
            area: self.area,
            center_of_mass: parent.apply([p[0], p[1], z]),
        }
    }
}

/// The plain layout, cut and given its gussets.
pub(crate) fn layout(spec: &PartSpec, plain: Layout) -> Result<Layout, SheetError> {
    let unsupported = |m: &str| Err(SheetError::Unsupported(format!("gussets: {m}")));
    if spec.base.cone.is_some() || spec.base.loft.is_some() {
        return unsupported("on a cone or a loft");
    }
    if !spec.tab_slots.is_empty() {
        return unsupported("with tabs and slots");
    }
    if !spec.faces.is_empty() {
        return unsupported("on a part of faces joined by bends");
    }
    let mut cut = spec.clone();
    cut.gussets.clear();
    let mut formed = Vec::new();
    let mut forms = Vec::new();
    for g in &spec.gussets {
        let f = new(spec, &plain, g)?;
        let corners = [
            f.flat(f.u_in[0], f.s_in[0]),
            f.flat(f.u_in[1], f.s_in[0]),
            f.flat(f.u_in[1], f.s_in[1]),
            f.flat(f.u_in[0], f.s_in[1]),
        ];
        let mut outline = corners.to_vec();
        if path_area(&outline, &[None; 4]) < 0.0 {
            outline.reverse();
        }
        // Both footprints: the groove's reaches past the rib's ends with a
        // wall thinner than the sheet.
        let (ui, ug, si, sg) = (f.u_in, f.u_groove, f.s_in, f.s_groove);
        let mut both = vec![[ui[0], si[0]]];
        if sg[0] < si[0] {
            both.extend([
                [ug[0], si[0]],
                [ug[0], sg[0]],
                [ug[1], sg[0]],
                [ug[1], si[0]],
            ]);
        }
        both.extend([[ui[1], si[0]], [ui[1], si[1]]]);
        if sg[1] > si[1] {
            both.extend([
                [ug[1], si[1]],
                [ug[1], sg[1]],
                [ug[0], sg[1]],
                [ug[0], si[1]],
            ]);
        }
        both.push([ui[0], si[1]]);
        let mut cutout: Vec<P2> = both.iter().map(|p| f.flat(p[0], p[1])).collect();
        if path_area(&cutout, &vec![None; cutout.len()]) < 0.0 {
            cutout.reverse();
        }
        cut.base.cutouts.push(Cutout {
            outline: cutout,
            bulges: Vec::new(),
        });
        let mid = |j: usize| (outline[0][j] + outline[2][j]) / 2.0;
        forms.push(FlatForm {
            id: g.id.clone(),
            center: [mid(0), mid(1)],
            outline: Ring {
                points: outline,
                arcs: vec![None; 4],
            },
        });
        formed.push(f);
    }
    let mut layout = crate::flat::layout_cells(&cut)?;
    if layout.folds.len() != plain.folds.len() {
        return Err(SheetError::Kernel("gussets: the folds changed".into()));
    }
    layout.flat = plain.flat;
    layout.flat.forms.extend(forms);
    layout.ribs = formed;
    Ok(layout)
}

/// A segment of a section: a line, or an arc about c of radius rho from
/// angle a0 to a1 (counter-clockwise from +s).
#[derive(Clone, Copy)]
enum Seg {
    Line(P2, P2),
    Arc { c: P2, rho: f64, a0: f64, a1: f64 },
}

/// Signed area and its first moments (of s, of z) of a closed path.
fn section(path: &[Seg]) -> (f64, f64, f64) {
    let (mut a, mut ms, mut mz) = (0.0, 0.0, 0.0);
    for seg in path {
        let pieces: Vec<Seg> = match *seg {
            Seg::Line(..) => vec![*seg],
            Seg::Arc { c, rho, a0, a1 } => {
                let n = ((a1 - a0).abs() / (PI / 8.0)).ceil().max(1.0) as usize;
                (0..n)
                    .map(|i| {
                        let f = |k: usize| a0 + (a1 - a0) * k as f64 / n as f64;
                        Seg::Arc {
                            c,
                            rho,
                            a0: f(i),
                            a1: f(i + 1),
                        }
                    })
                    .collect()
            }
        };
        for piece in pieces {
            for (x, w) in GAUSS {
                let (p, d) = match piece {
                    Seg::Line(p, q) => (
                        [p[0] + x * (q[0] - p[0]), p[1] + x * (q[1] - p[1])],
                        [q[0] - p[0], q[1] - p[1]],
                    ),
                    Seg::Arc { c, rho, a0, a1 } => {
                        let t = a0 + x * (a1 - a0);
                        let (sn, cs) = t.sin_cos();
                        (
                            [c[0] + rho * cs, c[1] + rho * sn],
                            [-rho * sn * (a1 - a0), rho * cs * (a1 - a0)],
                        )
                    }
                };
                a += (p[0] * d[1] - p[1] * d[0]) / 2.0 * w;
                ms += p[0] * p[0] / 2.0 * d[1] * w;
                mz -= p[1] * p[1] / 2.0 * d[0] * w;
            }
        }
    }
    (a, ms, mz)
}

fn new(spec: &PartSpec, plain: &Layout, g: &crate::spec::BendGusset) -> Result<Formed, SheetError> {
    let invalid = |m: &str| Err(SheetError::Invalid(format!("gusset {}: {m}", g.id)));
    let unsupported = |m: &str| Err(SheetError::Unsupported(format!("gusset {}: {m}", g.id)));
    if g.shape != GussetShape::Flat {
        return unsupported("a round one");
    }
    if g.draft != 0.0 {
        return unsupported("drafted side walls");
    }
    let t = spec.thickness;
    let tg = g.thickness.unwrap_or(t);
    if !(g.width > 0.0 && tg > 0.0) {
        return invalid("it needs a width and a wall");
    }
    if !(g.width > 2.0 * tg + EPS) {
        return invalid("its width leaves no groove between its walls");
    }
    let ks: Vec<usize> = (0..plain.folds.len())
        .filter(|&k| plain.folds[k].flange == g.on)
        .collect();
    let &[k] = ks.as_slice() else {
        return if ks.is_empty() {
            invalid(&format!("no flange {}", g.on))
        } else {
            unsupported(&format!("{} has more than one bend", g.on))
        };
    };
    let fold = &plain.folds[k];
    let theta = fold.angle;
    if !(theta > EPS && theta < PI - 1e-6) {
        return unsupported("on a bend of 180 degrees or more");
    }
    let (r, ro) = (fold.radius, fold.radius + t);
    let sigma = if fold.up { 1.0 } else { -1.0 };
    let (z_in, z_out) = if fold.up { (t, 0.0) } else { (0.0, t) };
    let (s0, s1) = (fold.s0, fold.s0 + fold.bend_allowance);
    let c = [s0, if fold.up { t + r } else { -r }];
    // On the bend at angle phi, radius rho; the flange's direction there.
    let on = |phi: f64, rho: f64| [c[0] + rho * phi.sin(), c[1] - sigma * rho * phi.cos()];
    let along = [theta.cos(), sigma * theta.sin()];
    let half = theta / 2.0;
    let sharp = [s0 + r * half.tan(), z_in];
    // The rib's inner face: from the base side's inner face to the flange's.
    let (l1, l2) = match (g.depth, g.legs) {
        (Some(d), None) => {
            if !(d > 0.0) {
                return invalid("its depth must be above 0");
            }
            let l = r * half.tan() - (r - d - r * half.cos()) / half.sin();
            (l, l)
        }
        (None, Some([l1, l2])) => (l1, l2),
        _ => return invalid("give its depth or its legs"),
    };
    let a = [sharp[0] - l1, z_in];
    let b = [sharp[0] + l2 * along[0], sharp[1] + l2 * along[1]];
    // Along the flange from the bend's end, on the inner face.
    let tangent_in = on(theta, r);
    let lb = (b[0] - tangent_in[0]) * along[0] + (b[1] - tangent_in[1]) * along[1];
    if !(a[0] < s0 - EPS && lb > EPS) {
        return invalid("its inner face must meet the inner faces past the bend");
    }
    // Into the wall: toward the sharp corner.
    let ab = [b[0] - a[0], b[1] - a[1]];
    let len = (ab[0] * ab[0] + ab[1] * ab[1]).sqrt();
    let mut normal = [-ab[1] / len, ab[0] / len];
    if (sharp[0] - a[0]) * normal[0] + (sharp[1] - a[1]) * normal[1] < 0.0 {
        normal = [-normal[0], -normal[1]];
    }
    // The groove's face, the wall further out, on the outer faces.
    let g0 = [a[0] + tg * normal[0], a[1] + tg * normal[1]];
    let a2 = [g0[0] + (z_out - g0[1]) / ab[1] * ab[0], z_out];
    if ab[1].abs() < EPS {
        return invalid("its inner face runs along the base");
    }
    let tangent_out = on(theta, ro);
    // Where the groove's line meets the flange's outer face: solve
    // g0 + x ab = tangent_out + y along.
    let det = ab[0] * (-along[1]) - ab[1] * (-along[0]);
    let rhs = [tangent_out[0] - g0[0], tangent_out[1] - g0[1]];
    let y = (ab[0] * rhs[1] - ab[1] * rhs[0]) / det;
    let b2 = [tangent_out[0] + y * along[0], tangent_out[1] + y * along[1]];
    let lb2 = y;
    if !(a2[0] < s0 - EPS && lb2 > EPS) {
        return invalid("its groove must meet the outer faces past the bend");
    }
    let (s_lo, s_hi) = (a[0].min(a2[0]), s1 + lb.max(lb2));

    // Along the bend: the middle of the bend strip by default.
    let u_of =
        |p: P2| (p[0] - fold.edge_point[0]) * fold.u[0] + (p[1] - fold.edge_point[1]) * fold.u[1];
    let s_of =
        |p: P2| (p[0] - fold.edge_point[0]) * fold.n[0] + (p[1] - fold.edge_point[1]) * fold.n[1];
    let bend_u = plain
        .cells
        .iter()
        .filter(|cell| cell.zone == Zone::Bend(k))
        .flat_map(|cell| cell.polygon.iter().map(|&p| u_of(p)))
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
            (lo.min(x), hi.max(x))
        });
    let u0 = g.at.unwrap_or((bend_u.0 + bend_u.1) / 2.0);
    let w = g.width;
    let u_in = [u0 - w / 2.0, u0 + w / 2.0];
    let u_groove = [u_in[0] + tg, u_in[1] - tg];
    let s_in = [a[0], s1 + lb];
    let s_groove = [a2[0], s1 + lb2];

    // Nothing in its way: no cell edge, cut-out or hole inside the
    // footprint but the bend's own lines, and the sheet under it of the
    // parent, the bend and the flange.
    let lo = [u_in[0] + 1e-7, s_lo + 1e-7];
    let hi = [u_in[1] - 1e-7, s_hi - 1e-7];
    let crosses = |p: P2, q: P2| -> bool {
        let (p, q) = ([u_of(p), s_of(p)], [u_of(q), s_of(q)]);
        let on_bend_line = |x: f64| (x - s0).abs() < 1e-7 || (x - s1).abs() < 1e-7;
        if on_bend_line(p[1]) && on_bend_line(q[1]) && (p[1] - q[1]).abs() < 1e-7 {
            return false;
        }
        // Liang-Barsky against the open box.
        let (mut t0, mut t1) = (0.0f64, 1.0f64);
        for j in 0..2 {
            let d = q[j] - p[j];
            for (bound, sign) in [(lo[j], -1.0), (hi[j], 1.0)] {
                let (num, den) = (sign * (bound - p[j]), sign * d);
                if den.abs() < 1e-15 {
                    if num < 0.0 {
                        return false;
                    }
                } else {
                    let x = num / den;
                    if den < 0.0 {
                        t0 = t0.max(x);
                    } else {
                        t1 = t1.min(x);
                    }
                }
            }
        }
        t0 < t1
    };
    let ring_crosses = |pts: &[P2], arcs: &[Option<Curve>]| {
        (0..pts.len()).any(|j| {
            let (p, q) = (pts[j], pts[(j + 1) % pts.len()]);
            match arcs[j] {
                None => crosses(p, q),
                // An arc: by its chords.
                Some(Curve::Arc(arc)) => {
                    let rr =
                        ((p[0] - arc.center[0]).powi(2) + (p[1] - arc.center[1]).powi(2)).sqrt();
                    let a0 = (p[1] - arc.center[1]).atan2(p[0] - arc.center[0]);
                    let sweep = arc.sweep(p, q);
                    let at = |x: f64| {
                        let t = a0 + sweep * x;
                        [arc.center[0] + rr * t.cos(), arc.center[1] + rr * t.sin()]
                    };
                    (0..16).any(|i| crosses(at(i as f64 / 16.0), at((i + 1) as f64 / 16.0)))
                }
                Some(Curve::Cut(_)) => true,
            }
        })
    };
    let in_way = plain.cells.iter().any(|cell| {
        ring_crosses(&cell.polygon, &cell.arcs)
            || cell
                .inner
                .iter()
                .any(|ring| ring_crosses(&ring.points, &ring.arcs))
            || cell.holes.iter().any(|h| {
                let (cu, cs) = (u_of(h.center), s_of(h.center));
                let du = (cu - cu.clamp(lo[0], hi[0])).abs();
                let ds = (cs - cs.clamp(lo[1], hi[1])).abs();
                du * du + ds * ds < (h.diameter / 2.0).powi(2)
            })
    });
    if in_way {
        return invalid("something is in its way (an edge, a hole or a cut-out)");
    }
    let flat_at = |u: f64, s: f64| {
        [
            fold.edge_point[0] + u * fold.u[0] + s * fold.n[0],
            fold.edge_point[1] + u * fold.u[1] + s * fold.n[1],
        ]
    };
    let parent_zone = fold.parent.map_or(Zone::Base, Zone::Flange);
    for (s, zone) in [
        ((s_lo + s0) / 2.0, parent_zone),
        ((s0 + s1) / 2.0, Zone::Bend(k)),
        ((s1 + s_hi) / 2.0, Zone::Flange(k)),
    ] {
        let p = flat_at(u0, s);
        if !plain
            .cells
            .iter()
            .any(|cell| cell.zone == zone && crate::cut::inside(p, &cell.polygon, &cell.arcs))
        {
            return invalid("it doesn't fit on its bend and the faces beside it");
        }
    }

    // The outer face under the walls: on the base side and the flange a U
    // round the groove, across the bend a strip under each wall.
    let cell = |name: String, zone: Zone, pts: Vec<[f64; 2]>| {
        let mut polygon: Vec<P2> = pts.iter().map(|p| flat_at(p[0], p[1])).collect();
        if path_area(&polygon, &vec![None; polygon.len()]) < 0.0 {
            polygon.reverse();
        }
        Cell {
            name,
            arcs: vec![None; polygon.len()],
            polygon,
            holes: Vec::new(),
            inner: Vec::new(),
            zone,
        }
    };
    let (ui, ug, si, sg) = (u_in, u_groove, s_in, s_groove);
    let (inner, outer) = (if fold.up { 1 } else { 0 }, if fold.up { 0 } else { 1 });
    let mut skins = Vec::new();
    let mut skin = |name: &str, zone: Zone, level: usize, pts: Vec<[f64; 2]>| {
        skins.push((
            cell(format!("gusset.{}.skin.{name}", g.id), zone, pts),
            level,
        ));
    };
    // The base side and the flange: the outer face round the groove (a U,
    // or the two sides when the groove reaches past), the inner face over
    // the groove's end past the rib's.
    for (name, zone, root, end_in, end_g) in [
        ("base", parent_zone, s0, si[0], sg[0]),
        ("flange", Zone::Flange(k), s1, si[1], sg[1]),
    ] {
        if (end_g - root).abs() > (end_in - root).abs() {
            skin(
                &format!("{name}0"),
                zone,
                outer,
                vec![
                    [ui[0], end_in],
                    [ug[0], end_in],
                    [ug[0], root],
                    [ui[0], root],
                ],
            );
            skin(
                &format!("{name}1"),
                zone,
                outer,
                vec![
                    [ug[1], end_in],
                    [ui[1], end_in],
                    [ui[1], root],
                    [ug[1], root],
                ],
            );
            skin(
                &format!("{name}.over"),
                zone,
                inner,
                vec![
                    [ug[0], end_g],
                    [ug[1], end_g],
                    [ug[1], end_in],
                    [ug[0], end_in],
                ],
            );
        } else {
            skin(
                name,
                zone,
                outer,
                vec![
                    [ui[0], end_in],
                    [ui[1], end_in],
                    [ui[1], root],
                    [ug[1], root],
                    [ug[1], end_g],
                    [ug[0], end_g],
                    [ug[0], root],
                    [ui[0], root],
                ],
            );
        }
    }
    skin(
        "bend0",
        Zone::Bend(k),
        outer,
        vec![[ui[0], s0], [ug[0], s0], [ug[0], s1], [ui[0], s1]],
    );
    skin(
        "bend1",
        Zone::Bend(k),
        outer,
        vec![[ug[1], s0], [ui[1], s0], [ui[1], s1], [ug[1], s1]],
    );

    // The sections: where the walls stand (the sheet from its outer face up
    // to the rib's inner face), the groove (outer face to the groove's face),
    // and the corner under the rib (inner face to the rib's inner face).
    let arc = |rho: f64| Seg::Arc {
        c,
        rho,
        a0: sigma * -FRAC_PI_2,
        a1: sigma * (theta - FRAC_PI_2),
    };
    let (t1_in, t1_out) = ([s0, z_in], [s0, z_out]);
    let b_out = [
        tangent_out[0] + lb * along[0],
        tangent_out[1] + lb * along[1],
    ];
    let a_out = [a[0], z_out];
    let walls = section(&[
        Seg::Line(a_out, t1_out),
        arc(ro),
        Seg::Line(tangent_out, b_out),
        Seg::Line(b_out, b),
        Seg::Line(b, a),
        Seg::Line(a, a_out),
    ]);
    // Between the walls, up to the groove's ends where they reach past.
    let far = lb.max(lb2);
    let a_lo = [s_lo, z_in];
    let b_hi = [
        tangent_in[0] + far * along[0],
        tangent_in[1] + far * along[1],
    ];
    let b_hi_out = [
        tangent_out[0] + far * along[0],
        tangent_out[1] + far * along[1],
    ];
    let between = section(&[
        Seg::Line([s_lo, z_out], t1_out),
        arc(ro),
        Seg::Line(tangent_out, b_hi_out),
        Seg::Line(b_hi_out, b_hi),
        Seg::Line(b_hi, b),
        Seg::Line(b, a),
        Seg::Line(a, a_lo),
        Seg::Line(a_lo, [s_lo, z_out]),
    ]);
    let groove = section(&[
        Seg::Line(a2, t1_out),
        arc(ro),
        Seg::Line(tangent_out, b2),
        Seg::Line(b2, a2),
    ]);
    let corner = section(&[
        Seg::Line(a, t1_in),
        arc(r),
        Seg::Line(tangent_in, b),
        Seg::Line(b, a),
    ]);
    // Signs: every path runs the same way round, so the signs agree.
    let (aw, ab2, ag) = (walls.0, between.0, groove.0);
    let side = 2.0 * tg;
    let middle = w - 2.0 * tg;
    let signed = aw * side + (ab2 - ag) * middle;
    let volume = signed.abs();
    let ms = walls.1 * side + (between.1 - groove.1) * middle;
    let mz = walls.2 * side + (between.2 - groove.2) * middle;
    let sign = signed.signum();
    let centre = [u0, ms * sign / volume, mz * sign / volume];
    let dist = |p: P2, q: P2| ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2)).sqrt();
    let skin_area: f64 = skins
        .iter()
        .map(|(cell, _)| match cell.zone {
            Zone::Bend(_) => tg * ro * theta,
            _ => path_area(&cell.polygon, &cell.arcs).abs(),
        })
        .sum();
    let area =
        w * dist(a, b) + middle * dist(a2, b2) + 2.0 * corner.0.abs() + 2.0 * ag.abs() + skin_area;
    Ok(Formed {
        id: g.id.clone(),
        fold: k,
        edge_point: fold.edge_point,
        u: fold.u,
        n: fold.n,
        u_in,
        s_in,
        u_groove,
        s_groove,
        s_bend: [s0, s1],
        inner,
        skins,
        normal,
        volume,
        area,
        centre,
    })
}

#[cfg(test)]
mod tests {
    use crate::spec::PartSpec;

    fn part(gussets: &str, angle: f64) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "holes": [] }},
                "flanges": [{{ "id": "f1", "on": {{ "base": 2 }}, "angle": {angle}, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }}],
                "gussets": [{gussets}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_gusset_across_a_90_degree_bend() {
        let spec = part(
            r#"{ "id": "g", "on": "f1", "at": 50, "width": 10, "depth": 10 }"#,
            90.0,
        );
        let exact = crate::fold::exact_mass(&spec).unwrap();
        // The reference CAD: 18111.750333.
        assert!(
            (exact.volume - 18111.750333).abs() < 1e-5,
            "{}",
            exact.volume
        );
        let brep = crate::fold::fold(&spec).unwrap();
        let m = crate::fold::mass(&brep, 0.002).unwrap();
        assert!(
            (m.volume - exact.volume).abs() / exact.volume < 1e-4,
            "{} {}",
            m.volume,
            exact.volume
        );
        assert!(
            (m.area - exact.area).abs() / exact.area < 1e-4,
            "{} {}",
            m.area,
            exact.area
        );
        for j in 0..3 {
            assert!((m.center_of_mass[j] - exact.center_of_mass[j]).abs() < 1e-3);
        }
        // The flat: the plain sheet, the footprint marked.
        let flat = crate::flat::flat_pattern(&spec).unwrap();
        let plain = crate::flat::flat_pattern(&part("", 90.0)).unwrap();
        assert_eq!(flat.outline, plain.outline);
        assert!(flat.cutouts.is_empty());
        assert_eq!(flat.forms.len(), 1);
    }

    #[test]
    fn a_wall_thinner_than_the_sheet() {
        // The groove reaches past the rib's footprint: the inner face goes on over it.
        let spec = part(
            r#"{ "id": "g", "on": "f1", "at": 50, "width": 10, "depth": 10, "thickness": 1 }"#,
            90.0,
        );
        let exact = crate::fold::exact_mass(&spec).unwrap();
        // The reference CAD: 17652.343422.
        assert!(
            (exact.volume - 17652.343422).abs() < 1e-5,
            "{}",
            exact.volume
        );
        let brep = crate::fold::fold(&spec).unwrap();
        let m = crate::fold::mass(&brep, 0.002).unwrap();
        assert!(
            (m.volume - exact.volume).abs() / exact.volume < 1e-4,
            "{} {}",
            m.volume,
            exact.volume
        );
        assert!(
            (m.area - exact.area).abs() / exact.area < 1e-4,
            "{} {}",
            m.area,
            exact.area
        );
    }

    #[test]
    fn gussets_across_45_and_135_degree_bends() {
        // The reference CAD's volumes, depth 5.
        for (angle, volume) in [(45.0, 18207.934521), (135.0, 16540.217192)] {
            let spec = part(
                r#"{ "id": "g", "on": "f1", "width": 10, "depth": 5 }"#,
                angle,
            );
            let exact = crate::fold::exact_mass(&spec).unwrap();
            assert!(
                (exact.volume - volume).abs() < 1e-5,
                "{angle}: {}",
                exact.volume
            );
            let brep = crate::fold::fold(&spec).unwrap();
            let m = crate::fold::mass(&brep, 0.002).unwrap();
            assert!(
                (m.volume - volume).abs() / volume < 1e-4,
                "{angle}: {}",
                m.volume
            );
        }
    }
}
