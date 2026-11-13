//! Holes and corner reliefs cut into the flat layout. A hole is a circle strictly inside the flat
//! outline, anywhere: on the base, across a bend, on a flange. A cell it
//! crosses loses the disc (its boundary follows the circle there, and it may
//! split into pieces); a cell that holds it whole gets it as an inner loop.

use crate::flat::{path_area, Arc, Cell, Curve, FlatPattern, Layout, Ring, P2};
use crate::spec::Hole;
use crate::SheetError;

const EPS: f64 = 1e-9;

pub(crate) fn cut_holes(layout: &mut Layout, holes: &[Hole]) -> Result<(), SheetError> {
    for (i, h) in holes.iter().enumerate() {
        check(&layout.flat, h)?;
        for g in &holes[..i] {
            if dist(h.center, g.center) <= (h.diameter + g.diameter) / 2.0 + EPS {
                return Err(SheetError::Invalid(format!(
                    "holes at {:?} and {:?} overlap",
                    g.center, h.center
                )));
            }
        }
    }
    for h in holes {
        let mut cells = Vec::new();
        for cell in layout.cells.drain(..) {
            let clip = Clip::Disc {
                c: h.center,
                r: h.diameter / 2.0,
            };
            cells.extend(cut(cell, &clip, false, &format!("hole at {:?}", h.center))?);
        }
        layout.cells = cells;
    }
    Ok(())
}

/// Whether a hole's circle crosses the flat outline (a round notch).
pub(crate) fn crosses_outline(flat: &FlatPattern, h: &Hole) -> bool {
    let (pts, arcs) = (&flat.outline, &flat.outline_arcs);
    (0..pts.len()).any(|k| {
        let (p, q) = (pts[k], pts[(k + 1) % pts.len()]);
        let d = match arc(arcs[k]) {
            None => segment_distance(h.center, p, q),
            Some(arc) => arc_distance(h.center, p, q, arc),
        };
        d < h.diameter / 2.0 - EPS
    })
}

/// A hole must lie inside the flat outline without touching it.
fn check(flat: &FlatPattern, h: &Hole) -> Result<(), SheetError> {
    let r = h.diameter / 2.0;
    let invalid = |what: &str| {
        Err(SheetError::Invalid(format!(
            "hole at {:?}: {what}",
            h.center
        )))
    };
    if !(r > EPS) {
        return invalid("needs a positive diameter");
    }
    let (pts, arcs) = (&flat.outline, &flat.outline_arcs);
    for k in 0..pts.len() {
        let (p, q) = (pts[k], pts[(k + 1) % pts.len()]);
        let d = match arc(arcs[k]) {
            None => segment_distance(h.center, p, q),
            Some(arc) => arc_distance(h.center, p, q, arc),
        };
        if d <= r + EPS {
            return invalid("touches or crosses the outline");
        }
    }
    if !inside(h.center, pts, arcs) {
        return invalid("lies outside the part");
    }
    Ok(())
}

/// One crossing of the clip's boundary with a cell's boundary.
#[derive(Clone, Copy)]
struct Crossing {
    /// Boundary segment and position along it (0..1).
    segment: usize,
    t: f64,
    p: P2,
    /// The boundary runs into the clip here.
    entering: bool,
}

/// A region cut out of the flat: a disc, or a profile of lines and arcs
/// (counter-clockwise, convex or not).
#[derive(Clone, Debug)]
pub(crate) enum Clip {
    Disc { c: P2, r: f64 },
    Profile(Ring),
}

impl Clip {
    pub(crate) fn contains(&self, p: P2) -> bool {
        match self {
            Clip::Disc { c, r } => dist(p, *c) < *r,
            Clip::Profile(ring) => inside(p, &ring.points, &ring.arcs),
        }
    }

    /// Whether p lies on the boundary (to EPS).
    pub(crate) fn on_boundary(&self, p: P2) -> bool {
        match self {
            Clip::Disc { c, r } => (dist(p, *c) - r).abs() <= EPS,
            Clip::Profile(ring) => (0..ring.points.len()).any(|k| {
                let (a, b) = segment(ring, k);
                match arc(ring.arcs[k]) {
                    None => segment_distance(p, a, b) <= EPS,
                    Some(arc) => arc_distance(p, a, b, arc) <= EPS,
                }
            }),
        }
    }

    /// Parameters t in (0, 1) where a -> b crosses the boundary, in order;
    /// Err when it only touches it, or passes through a corner of it.
    pub(crate) fn crossings(&self, a: P2, b: P2) -> Result<Vec<f64>, ()> {
        match self {
            Clip::Disc { c, r } => line_circle(a, b, *c, *r)
                .map(|ts| ts.into_iter().filter(|t| *t > 0.0 && *t < 1.0).collect())
                .map_err(|_| ()),
            Clip::Profile(ring) => {
                let d = sub(b, a);
                let mut ts = Vec::new();
                for k in 0..ring.points.len() {
                    let (c0, c1) = segment(ring, k);
                    if segment_distance(c0, a, b) <= EPS {
                        // Through a corner of the profile.
                        return Err(());
                    }
                    match arc(ring.arcs[k]) {
                        None => {
                            let e = sub(c1, c0);
                            let den = cross(d, e);
                            if den.abs() <= EPS * dot(d, d).sqrt() * dot(e, e).sqrt() {
                                if segment_distance(a, c0, c1) <= EPS
                                    || segment_distance(b, c0, c1) <= EPS
                                {
                                    return Err(());
                                }
                                continue;
                            }
                            let w = sub(c0, a);
                            let (t, u) = (cross(w, e) / den, cross(w, d) / den);
                            if t > 0.0 && t < 1.0 && u > 0.0 && u < 1.0 {
                                ts.push(t);
                            }
                        }
                        Some(arc) => {
                            let r = dist(c0, arc.center);
                            // Touching the circle matters only on the arc.
                            let hits = match line_circle(a, b, arc.center, r) {
                                Ok(ts) => ts,
                                Err(t) if on_arc(add(a, scale(d, t)), c0, c1, arc) => {
                                    return Err(())
                                }
                                Err(_) => Vec::new(),
                            };
                            for t in hits {
                                if t > 0.0 && t < 1.0 && on_arc(add(a, scale(d, t)), c0, c1, arc) {
                                    ts.push(t);
                                }
                            }
                        }
                    }
                }
                ts.sort_by(f64::total_cmp);
                Ok(ts)
            }
        }
    }

    /// Position of a boundary point counter-clockwise round the clip, and the
    /// full turn.
    fn param(&self, p: P2) -> (f64, f64) {
        match self {
            Clip::Disc { c, .. } => {
                let tau = std::f64::consts::TAU;
                ((p[1] - c[1]).atan2(p[0] - c[0]).rem_euclid(tau), tau)
            }
            Clip::Profile(ring) => {
                let (k, local) = locate(ring, p);
                (starts(ring)[k] + local, starts(ring)[ring.points.len()])
            }
        }
    }

    /// The boundary clockwise from e to l, as path pieces: e first, each
    /// point with the segment that leaves it.
    fn clockwise(&self, e: P2, l: P2, pts: &mut Vec<P2>, arcs: &mut Vec<Option<Curve>>) {
        match self {
            Clip::Disc { c, .. } => {
                pts.push(e);
                arcs.push(Some(Curve::Arc(Arc {
                    center: *c,
                    ccw: false,
                })));
            }
            Clip::Profile(ring) => {
                let n = ring.points.len();
                let reverse = |k: usize| {
                    arc(ring.arcs[k]).map(|a| {
                        Curve::Arc(Arc {
                            center: a.center,
                            ccw: !a.ccw,
                        })
                    })
                };
                let (mut k, le) = locate(ring, e);
                let (kl, ll) = locate(ring, l);
                pts.push(e);
                arcs.push(reverse(k));
                // Back to the start of each segment until l's segment, where
                // l lies behind the current point.
                let mut first = true;
                loop {
                    if k == kl && (!first || ll < le) {
                        break;
                    }
                    first = false;
                    let prev = (k + n - 1) % n;
                    pts.push(ring.points[k]);
                    arcs.push(reverse(prev));
                    k = prev;
                }
            }
        }
    }
}

/// Segment k of a ring: its end points.
fn segment(ring: &Ring, k: usize) -> (P2, P2) {
    (ring.points[k], ring.points[(k + 1) % ring.points.len()])
}

fn segment_length(ring: &Ring, k: usize) -> f64 {
    let (a, b) = segment(ring, k);
    match arc(ring.arcs[k]) {
        None => dist(a, b),
        Some(arc) => dist(a, arc.center) * arc.sweep(a, b).abs(),
    }
}

/// Where each segment starts along the ring, and (last) its full length.
fn starts(ring: &Ring) -> Vec<f64> {
    let mut out = vec![0.0];
    for k in 0..ring.points.len() {
        out.push(out[k] + segment_length(ring, k));
    }
    out
}

/// The segment a boundary point lies on and how far along it.
fn locate(ring: &Ring, p: P2) -> (usize, f64) {
    let mut best = (f64::MAX, 0, 0.0);
    for k in 0..ring.points.len() {
        let (a, b) = segment(ring, k);
        let (d, along) = match arc(ring.arcs[k]) {
            None => {
                let v = sub(b, a);
                let t = (dot(sub(p, a), v) / dot(v, v)).clamp(0.0, 1.0);
                (dist(p, add(a, scale(v, t))), t * dist(a, b))
            }
            Some(arc) => {
                let r = dist(a, arc.center);
                let part = arc.sweep(a, p);
                let total = arc.sweep(a, b);
                let along = if (total > 0.0) == (part > 0.0) && part.abs() <= total.abs() + 1e-12 {
                    part.abs()
                } else {
                    0.0
                };
                (arc_distance(p, a, b, arc), along * r)
            }
        };
        if d < best.0 {
            best = (d, k, along);
        }
    }
    (best.1, best.2)
}

/// Parameters where the line a + t (b - a) meets a circle; Err(t) when it
/// only touches it, at t within the segment.
fn line_circle(a: P2, b: P2, c: P2, r: f64) -> Result<Vec<f64>, f64> {
    let d = sub(b, a);
    let f = sub(a, c);
    let (qa, qb, qc) = (dot(d, d), 2.0 * dot(f, d), dot(f, f) - r * r);
    let disc = qb * qb - 4.0 * qa * qc;
    if disc.abs() <= EPS * qa {
        let t = -qb / (2.0 * qa);
        return if (0.0..=1.0).contains(&t) {
            Err(t)
        } else {
            Ok(Vec::new())
        };
    }
    if disc < 0.0 {
        return Ok(Vec::new());
    }
    Ok([-1.0, 1.0]
        .map(|s| (-qb + s * disc.sqrt()) / (2.0 * qa))
        .to_vec())
}

/// The cell minus the clip: zero, one or several cells. A clip wholly
/// inside the cell becomes a hole (discs only); a cell wholly inside the clip
/// goes when `swallow`, and is refused otherwise.
pub(crate) fn cut(
    cell: Cell,
    clip: &Clip,
    swallow: bool,
    what: &str,
) -> Result<Vec<Cell>, SheetError> {
    let pts = &cell.polygon;
    let n = pts.len();
    let unsupported = |why: &str| Err(SheetError::Unsupported(format!("{what}: {why}")));
    let mut crossings = Vec::new();
    let on: Vec<bool> = pts.iter().map(|&p| clip.on_boundary(p)).collect();
    // Just off p towards q.
    let near = |p: P2, q: P2| add(p, scale(sub(q, p), 1e-6 / dist(p, q).max(1e-300)));
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        if on[k] {
            // A corner on the clip's edge: a crossing where the boundary runs
            // from one side to the other there, a touch otherwise.
            let prev = pts[(k + n - 1) % n];
            let (before, after) = (clip.contains(near(p, prev)), clip.contains(near(p, q)));
            if before != after {
                crossings.push(Crossing {
                    segment: k,
                    t: 0.0,
                    p,
                    entering: after,
                });
            }
        }
        if let Some(Curve::Cut(c)) = cell.arcs[k] {
            // A mitered bend's end: where the clip's boundary crosses it.
            let (sp, sq) = (c.s_of(p), c.s_of(q));
            for s in cut_crossings(&c, sp, sq, clip) {
                let point = c.point(s);
                if (on[k] && dist(point, p) <= EPS) || (on[(k + 1) % n] && dist(point, q) <= EPS) {
                    continue;
                }
                let ahead = c.point(s + (sq - sp).signum() * 1e-6);
                crossings.push(Crossing {
                    segment: k,
                    t: (s - sp) / (sq - sp),
                    p: point,
                    entering: clip.contains(ahead),
                });
            }
            continue;
        }
        if cell.arcs[k].is_some() {
            // Arcs are outline (relief) segments, which cuts keep clear of.
            continue;
        }
        // Solve in a canonical direction so neighbouring cells get the very
        // same points for their shared edge.
        let flip = (p[0], p[1]) > (q[0], q[1]);
        let (a, b) = if flip { (q, p) } else { (p, q) };
        let Ok(ts) = clip.crossings(a, b) else {
            return unsupported("its edge touches an edge of the part");
        };
        for t in ts {
            let point = add(a, scale(sub(b, a), t));
            let t = if flip { 1.0 - t } else { t };
            // At a corner on the clip's edge: handled with the corner.
            if (on[k] && dist(point, p) <= EPS) || (on[(k + 1) % n] && dist(point, q) <= EPS) {
                continue;
            }
            // Into the clip when just past the crossing is inside.
            let ahead = add(point, scale(sub(q, p), 1e-6 / dist(p, q).max(1e-300)));
            crossings.push(Crossing {
                segment: k,
                t,
                p: point,
                entering: clip.contains(ahead),
            });
        }
    }

    if crossings.is_empty() {
        if pts.iter().all(|&p| clip.contains(p)) && swallow {
            return Ok(Vec::new());
        }
        if pts.iter().any(|&p| clip.contains(p)) {
            return unsupported("it swallows a corner of a cell");
        }
        let mut cell = cell;
        match clip {
            Clip::Disc { c, r } => {
                if inside(*c, &cell.polygon, &cell.arcs) {
                    cell.holes.push(Hole {
                        center: *c,
                        diameter: 2.0 * r,
                    });
                }
            }
            Clip::Profile(ring) => {
                if inside(ring.points[0], &cell.polygon, &cell.arcs) {
                    cell.inner.push(ring.reversed());
                }
            }
        }
        return Ok(vec![cell]);
    }
    crossings.sort_by(|x, y| (x.segment, x.t).partial_cmp(&(y.segment, y.t)).unwrap());

    // Walk: from a crossing where the boundary leaves the clip, along the
    // boundary to the next crossing (it enters there), then clockwise round
    // the clip to the next crossing (the boundary leaves there), and so on.
    let mut visited = vec![false; crossings.len()];
    let mut pieces = Vec::new();
    while let Some(start) = (0..crossings.len()).find(|&i| !visited[i] && !crossings[i].entering) {
        let (mut pts_out, mut arcs_out): (Vec<P2>, Vec<Option<Curve>>) = (Vec::new(), Vec::new());
        let mut cur = start;
        loop {
            visited[cur] = true;
            let x = crossings[cur];
            pts_out.push(x.p);
            // The rest of x's segment: a line, or the same curve.
            arcs_out.push(cell.arcs[x.segment]);
            // Along the boundary: the next crossing in boundary order.
            let next = (cur + 1) % crossings.len();
            let e = crossings[next];
            if !e.entering {
                return unsupported("its crossings with a cell don't alternate");
            }
            // The vertices in between: none when e is further along the same
            // segment, all the way round when it's behind x on it.
            if !(e.segment == x.segment && e.t > x.t) {
                let mut k = x.segment;
                loop {
                    k = (k + 1) % n;
                    pts_out.push(pts[k]);
                    arcs_out.push(cell.arcs[k]);
                    if k == e.segment {
                        break;
                    }
                }
            }
            visited[next] = true;
            // Round the clip, clockwise: the nearest crossing.
            let (from, turn) = clip.param(e.p);
            let after = (0..crossings.len())
                .filter(|&i| i != next)
                .min_by(|&i, &j| {
                    let di = (from - clip.param(crossings[i].p).0).rem_euclid(turn);
                    let dj = (from - clip.param(crossings[j].p).0).rem_euclid(turn);
                    di.total_cmp(&dj)
                })
                .unwrap();
            if crossings[after].entering {
                return unsupported("its crossings with a cell don't alternate");
            }
            clip.clockwise(e.p, crossings[after].p, &mut pts_out, &mut arcs_out);
            if after == start {
                break;
            }
            cur = after;
        }
        if path_area(&pts_out, &arcs_out) <= EPS {
            return unsupported("it leaves a sliver of a cell");
        }
        pieces.push((pts_out, arcs_out));
    }
    if visited.iter().any(|v| !v) {
        return unsupported("its crossings with a cell don't pair up");
    }
    let many = pieces.len() > 1;
    Ok(pieces
        .into_iter()
        .enumerate()
        .map(|(i, (polygon, arcs))| Cell {
            // `#`, not `.`: a contour's later steps are named `c.2`, `c.3`.
            name: if many {
                format!("{}#{i}", cell.name)
            } else {
                cell.name.clone()
            },
            // Earlier holes held whole stay with the piece around them.
            holes: cell
                .holes
                .iter()
                .filter(|g| inside(g.center, &polygon, &arcs))
                .cloned()
                .collect(),
            inner: cell
                .inner
                .iter()
                .filter(|g| inside(g.points[0], &polygon, &arcs))
                .cloned()
                .collect(),
            polygon,
            arcs,
            zone: cell.zone,
        })
        .collect())
}

/// Cuts the flat pattern's cut-outs out of the cells, after checking each
/// lies inside the outline, clear of it, the holes and the other cut-outs.
pub(crate) fn cut_cutouts(layout: &mut Layout) -> Result<(), SheetError> {
    let rings = layout.flat.cutouts.clone();
    let clips: Vec<Clip> = rings.iter().cloned().map(Clip::Profile).collect();
    let flat = &layout.flat;
    for (k, ring) in rings.iter().enumerate() {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("cut-out {k}: {what}")));
        let n = ring.points.len();
        if n < 2 || ring.arcs.len() != n || ring.area() <= EPS {
            return invalid("needs a counter-clockwise outline");
        }
        let outline = Ring {
            points: flat.outline.clone(),
            arcs: flat.outline_arcs.clone(),
        };
        // One that crosses the outline is a notch; one wholly outside is wrong.
        if !touches(&clips[k], &outline) && !inside(ring.points[0], &outline.points, &outline.arcs)
        {
            return invalid("lies outside the part");
        }
        for h in &flat.holes {
            let disc = Clip::Disc {
                c: h.center,
                r: h.diameter / 2.0,
            };
            if touches(&disc, ring) || clips[k].contains(h.center) {
                return invalid("overlaps a hole");
            }
        }
        for other in &rings[..k] {
            if touches(&clips[k], other) || clips[k].contains(other.points[0]) {
                return invalid("overlaps another cut-out");
            }
        }
    }
    let mut kept = Vec::new();
    for (k, clip) in clips.iter().enumerate() {
        let what = format!("cut-out {k}");
        let outline = Ring {
            points: layout.flat.outline.clone(),
            arcs: layout.flat.outline_arcs.clone(),
        };
        if touches(clip, &outline) {
            // A notch: it becomes part of the outline.
            cut_across(layout, clip, true, false, &what)?;
            continue;
        }
        let mut cells = Vec::new();
        for cell in layout.cells.drain(..) {
            cells.extend(cut(cell, clip, false, &what)?);
        }
        layout.cells = cells;
        kept.push(rings[k].clone());
    }
    layout.flat.cutouts = kept;
    Ok(())
}

/// Cuts a closed ring out of the cells only (an emboss's footprint: the flat
/// pattern stays whole), clear of the outline, holes, cut-outs, `discs` and
/// `rings`.
pub(crate) fn cut_inner(
    layout: &mut Layout,
    ring: &Ring,
    discs: &[Hole],
    rings: &[Ring],
    what: &str,
) -> Result<(), SheetError> {
    let invalid = |m: &str| Err(SheetError::Invalid(format!("{what}: {m}")));
    let clip = Clip::Profile(ring.clone());
    let outline = Ring {
        points: layout.flat.outline.clone(),
        arcs: layout.flat.outline_arcs.clone(),
    };
    if touches(&clip, &outline) {
        return invalid("crosses the outline");
    }
    if !inside(ring.points[0], &outline.points, &outline.arcs) {
        return invalid("lies outside the part");
    }
    for h in layout.flat.holes.iter().chain(discs) {
        let disc = Clip::Disc {
            c: h.center,
            r: h.diameter / 2.0,
        };
        if touches(&disc, ring) || clip.contains(h.center) {
            return invalid("overlaps a hole or a round emboss");
        }
    }
    for other in layout.flat.cutouts.iter().chain(rings) {
        if touches(&clip, other)
            || clip.contains(other.points[0])
            || Clip::Profile(other.clone()).contains(ring.points[0])
        {
            return invalid("overlaps a cut-out or another emboss");
        }
    }
    let mut cells = Vec::new();
    for cell in layout.cells.drain(..) {
        cells.extend(cut(cell, &clip, false, what)?);
    }
    layout.cells = cells;
    Ok(())
}

/// Whether a ring meets the clip's boundary, or starts inside it (arcs of
/// the ring checked at points along them).
fn touches(clip: &Clip, ring: &Ring) -> bool {
    // A profile partly inside the ring and partly outside crosses it (a
    // long arc of the ring can pass between the samples below).
    if let Clip::Profile(profile) = clip {
        let mut sides = samples(profile).map(|x| inside(x, &ring.points, &ring.arcs));
        let first = sides.next();
        if sides.any(|s| Some(s) != first) {
            return true;
        }
    }
    let n = ring.points.len();
    (0..n).any(|k| {
        let (a, b) = segment(ring, k);
        if clip.contains(a) || clip.on_boundary(a) {
            return true;
        }
        match arc(ring.arcs[k]) {
            None => !matches!(clip.crossings(a, b), Ok(ts) if ts.is_empty()),
            Some(arc) => {
                let (r, a0, sweep) = (
                    dist(a, arc.center),
                    (a[1] - arc.center[1]).atan2(a[0] - arc.center[0]),
                    arc.sweep(a, b),
                );
                (1..32).any(|i| {
                    let t = a0 + sweep * i as f64 / 32.0;
                    let x = [arc.center[0] + r * t.cos(), arc.center[1] + r * t.sin()];
                    clip.contains(x) || clip.on_boundary(x)
                })
            }
        }
    })
}

/// A ring's corners, and 31 points along each of its arcs.
fn samples(ring: &Ring) -> impl Iterator<Item = P2> + '_ {
    let n = ring.points.len();
    (0..n).flat_map(move |k| {
        let (a, b) = segment(ring, k);
        let along = arc(ring.arcs[k]).map(|arc| {
            let (r, a0, sweep) = (
                dist(a, arc.center),
                (a[1] - arc.center[1]).atan2(a[0] - arc.center[0]),
                arc.sweep(a, b),
            );
            (1..32).map(move |i| {
                let t = a0 + sweep * i as f64 / 32.0;
                [arc.center[0] + r * t.cos(), arc.center[1] + r * t.sin()]
            })
        });
        std::iter::once(a).chain(along.into_iter().flatten())
    })
}

/// Cuts a corner relief out of the layout.
/// A corner relief; `cells` false cuts the flat pattern only (outline and
/// bend lines), for a small one.
pub(crate) fn cut_corner(layout: &mut Layout, clip: &Clip, cells: bool) -> Result<(), SheetError> {
    cut_across(layout, clip, cells, true, "corner relief")
}

/// Cuts a region across the outline (a corner relief, a notch): from the
/// cells it crosses (cells inside it go), the flat outline and the bend
/// lines. The part must stay in one piece. At a corner a bend line that
/// passes through near its end loses the stub beyond.
fn cut_across(
    layout: &mut Layout,
    clip: &Clip,
    cut_cells: bool,
    corner: bool,
    what: &str,
) -> Result<(), SheetError> {
    if cut_cells {
        let mut cells = Vec::new();
        for cell in layout.cells.drain(..) {
            cells.extend(cut(cell, clip, true, what)?);
        }
        layout.cells = cells;
    }
    let outline = Cell {
        name: "outline".into(),
        polygon: std::mem::take(&mut layout.flat.outline),
        arcs: std::mem::take(&mut layout.flat.outline_arcs),
        holes: Vec::new(),
        inner: Vec::new(),
        zone: crate::flat::Zone::Base,
    };
    let mut pieces = cut(outline, clip, false, what)?;
    if pieces.len() != 1 {
        return Err(SheetError::Unsupported(format!(
            "{what}: it splits the part"
        )));
    }
    let piece = pieces.remove(0);
    layout.flat.outline = piece.polygon;
    layout.flat.outline_arcs = piece.arcs;
    // Bend lines end where the relief cuts them.
    for bend in &mut layout.flat.bends {
        let [a, b] = bend.line;
        let Ok(ts) = clip.crossings(a, b) else {
            continue;
        };
        if clip.contains(a) {
            if let Some(&t) = ts.first() {
                bend.line[0] = add(a, scale(sub(b, a), t));
            }
        }
        if clip.contains(b) {
            if let Some(&t) = ts.last() {
                bend.line[1] = add(a, scale(sub(b, a), t));
            }
        }
        // Through it near one end (a bend run on to a sloped flange end, past
        // a small relief): the stub beyond goes.
        if corner && !clip.contains(a) && !clip.contains(b) && ts.len() >= 2 {
            let (t0, t1) = (ts[0], ts[ts.len() - 1]);
            if t0 > 1.0 - t1 {
                bend.line[1] = add(a, scale(sub(b, a), t0));
            } else {
                bend.line[0] = add(a, scale(sub(b, a), t1));
            }
        }
    }
    Ok(())
}

/// Flat s where a mitered bend's end from s = a to b crosses the clip's
/// boundary, in order from a. Found in s ascending whichever way it runs, so
/// a cell and the outline sharing the curve get the same points.
fn cut_crossings(c: &crate::flat::Cut, a: f64, b: f64, clip: &Clip) -> Vec<f64> {
    let (lo, hi) = (a.min(b), a.max(b));
    let pieces = 256;
    let s = |j: usize| lo + (hi - lo) * j as f64 / pieces as f64;
    let mut out = Vec::new();
    for j in 0..pieces {
        let (mut s0, mut s1) = (s(j), s(j + 1));
        let inside = clip.contains(c.point(s0));
        if inside == clip.contains(c.point(s1)) {
            continue;
        }
        for _ in 0..80 {
            let m = (s0 + s1) / 2.0;
            if clip.contains(c.point(m)) == inside {
                s0 = m;
            } else {
                s1 = m;
            }
        }
        out.push((s0 + s1) / 2.0);
    }
    if a > b {
        out.reverse();
    }
    out
}

/// Whether p lies inside a closed boundary of lines and arcs: crossings of a
/// ray in an odd direction (so it misses vertices and tangencies).
pub(crate) fn inside(p: P2, pts: &[P2], arcs: &[Option<Curve>]) -> bool {
    if arcs.iter().any(|c| matches!(c, Some(Curve::Cut(_)))) {
        // A mitered bend's end as short lines (it never comes near a hole's
        // centre or a cut-out's first point).
        let (mut q, mut qa) = (Vec::new(), Vec::new());
        for k in 0..pts.len() {
            let (a, b) = (pts[k], pts[(k + 1) % pts.len()]);
            match arcs[k] {
                Some(Curve::Cut(c)) => {
                    let line = c.points(a, b, 64);
                    qa.extend(std::iter::repeat_n(None, line.len() - 1));
                    q.extend(&line[..line.len() - 1]);
                }
                other => {
                    q.push(a);
                    qa.push(other);
                }
            }
        }
        return inside(p, &q, &qa);
    }
    let d = [1.0, 0.371_293_417];
    let mut count = 0;
    for k in 0..pts.len() {
        let (a, b) = (pts[k], pts[(k + 1) % pts.len()]);
        match arc(arcs[k]) {
            None => {
                // p + s d = a + t (b - a), s > 0, 0 <= t < 1.
                let e = sub(b, a);
                let den = cross(d, e);
                if den.abs() < 1e-15 {
                    continue;
                }
                let w = sub(a, p);
                let (s, t) = (cross(w, e) / den, cross(w, d) / den);
                if s > 0.0 && (0.0..1.0).contains(&t) {
                    count += 1;
                }
            }
            Some(arc) => {
                let f = sub(p, arc.center);
                let rr = dot(sub(a, arc.center), sub(a, arc.center));
                let (qa, qb, qc) = (dot(d, d), 2.0 * dot(f, d), dot(f, f) - rr);
                let disc = qb * qb - 4.0 * qa * qc;
                if disc <= 0.0 {
                    continue;
                }
                for s in [-1.0, 1.0] {
                    let t = (-qb + s * disc.sqrt()) / (2.0 * qa);
                    if t > 0.0 && on_arc(add(p, scale(d, t)), a, b, arc) {
                        count += 1;
                    }
                }
            }
        }
    }
    count % 2 == 1
}

/// Whether point x (on the arc's circle) lies on the arc from a to b.
fn on_arc(x: P2, a: P2, b: P2, arc: Arc) -> bool {
    let sweep = arc.sweep(a, b);
    let part = arc.sweep(a, x);
    if sweep > 0.0 {
        part > 0.0 && part < sweep
    } else {
        part < 0.0 && part > sweep
    }
}

fn segment_distance(x: P2, a: P2, b: P2) -> f64 {
    let d = sub(b, a);
    let t = (dot(sub(x, a), d) / dot(d, d)).clamp(0.0, 1.0);
    dist(x, add(a, scale(d, t)))
}

fn arc_distance(x: P2, a: P2, b: P2, arc: Arc) -> f64 {
    let r = dist(a, arc.center);
    let to = sub(x, arc.center);
    let l = dot(to, to).sqrt();
    if l > EPS {
        let foot = add(arc.center, scale(to, r / l));
        if on_arc(foot, a, b, arc) {
            return (l - r).abs();
        }
    }
    dist(x, a).min(dist(x, b))
}

fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale(a: P2, s: f64) -> P2 {
    [a[0] * s, a[1] * s]
}
fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn dist(a: P2, b: P2) -> f64 {
    dot(sub(a, b), sub(a, b)).sqrt()
}

/// The arc of a segment that is one (a mitered bend's end is not).
fn arc(c: Option<Curve>) -> Option<Arc> {
    c.and_then(Curve::arc)
}
