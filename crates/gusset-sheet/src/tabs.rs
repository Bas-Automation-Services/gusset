//! Tabs on a flange's tip and slots where that tip lands on a face of the same
//! part. Laid out twice: once without them, to fold the tip onto the face and
//! find each tab's footprint there (the slots, in the face's flat frame);
//! then with the slots as cut-outs, the tabs added to the tip as bumps.

use crate::flat::{layout_cells, path_area, Curve, Layout, Zone, P2};
use crate::spec::{Cutout, PartSpec};
use crate::SheetError;

const EPS: f64 = 1e-9;

pub(crate) fn layout(spec: &PartSpec) -> Result<Layout, SheetError> {
    let invalid = |m: String| Err(SheetError::Invalid(format!("tab_slots: {m}")));
    let unsupported = |m: &str| Err(SheetError::Unsupported(format!("tab_slots: {m}")));
    if spec.base.cone.is_some() || spec.base.loft.is_some() {
        return unsupported("on a cone or a loft");
    }
    let mut plain = spec.clone();
    plain.tab_slots.clear();
    let first = layout_cells(&plain)?;
    let model = crate::fold::Model::new(&plain, &first);
    let t = spec.thickness;
    let mut cut = plain.clone();
    // Per tip: the fold, the tip's s, the tabs' t ranges, their height.
    let mut bumps: Vec<(usize, f64, Vec<(f64, f64)>, f64)> = Vec::new();
    for ts in &spec.tab_slots {
        let fold_of = |id: &str| first.folds.iter().rposition(|f| f.flange == id);
        let (Some(kl), Some(kf)) = (fold_of(&ts.tip), fold_of(&ts.face)) else {
            return invalid(format!("no flange {} or {}", ts.tip, ts.face));
        };
        if !(ts.length > 0.0) || ts.count == 0 {
            return invalid("tabs need a length and a count".into());
        }
        let lip = &first.folds[kl];
        let (s_of, t_of) = (
            |p: P2| dot(sub(p, lip.edge_point), lip.n),
            |p: P2| dot(sub(p, lip.edge_point), lip.u),
        );
        let at = |tt: f64, s: f64| add(lip.edge_point, add(scale(lip.u, tt), scale(lip.n, s)));
        // The tip: the tip flange's furthest edge, square to its edge.
        let cell = first
            .cells
            .iter()
            .filter(|c| c.zone == Zone::Flange(kl))
            .max_by(|a, b| {
                path_area(&a.polygon, &a.arcs).total_cmp(&path_area(&b.polygon, &b.arcs))
            })
            .ok_or_else(|| {
                SheetError::Invalid(format!("tab_slots: {} has no flat face", ts.tip))
            })?;
        let s_tip = cell
            .polygon
            .iter()
            .map(|&p| s_of(p))
            .fold(f64::MIN, f64::max);
        let on_tip: Vec<f64> = cell
            .polygon
            .iter()
            .filter(|&&p| (s_of(p) - s_tip).abs() < 1e-7)
            .map(|&p| t_of(p))
            .collect();
        let (t_lo, t_hi) = (
            on_tip.iter().copied().fold(f64::MAX, f64::min),
            on_tip.iter().copied().fold(f64::MIN, f64::max),
        );
        let (o0, o1) = (ts.offset_start.unwrap_or(0.0), ts.offset_end.unwrap_or(0.0));
        let span = t_hi - t_lo - o0 - o1;
        if span < ts.length * ts.count as f64 - EPS {
            return invalid(format!(
                "{} tabs of {} don't fit on {}'s tip",
                ts.count, ts.length, ts.tip
            ));
        }
        let ranges: Vec<(f64, f64)> = (0..ts.count)
            .map(|k| {
                let a = if ts.count == 1 {
                    t_lo + o0
                } else {
                    t_lo + o0 + k as f64 * (span - ts.length) / (ts.count - 1) as f64
                };
                (a, a + ts.length)
            })
            .collect();
        // Folded: the tip's way out and the face's normal (from its z = 0
        // side to its z = T side).
        let corner = |tt: f64, z: f64| model.map(Zone::Flange(kl), at(tt, s_tip), z);
        let out = unit3(sub3(
            model.map(Zone::Flange(kl), at(t_lo, s_tip + 1.0), 0.0),
            corner(t_lo, 0.0),
        ));
        let face_cell = first
            .cells
            .iter()
            .find(|c| c.zone == Zone::Flange(kf))
            .ok_or_else(|| {
                SheetError::Invalid(format!("tab_slots: {} has no flat face", ts.face))
            })?;
        let p0 = face_cell.polygon[0];
        let normal = sub3(
            model.map(Zone::Flange(kf), p0, 1.0),
            model.map(Zone::Flange(kf), p0, 0.0),
        );
        let toward = dot3(out, normal);
        if (toward.abs() - 1.0).abs() > 1e-9 {
            return unsupported("a tip that doesn't meet its face square on");
        }
        // The face it meets first, and the gap to it.
        let z_face = if toward > 0.0 { 0.0 } else { t };
        let level = |p: [f64; 3]| dot3(p, normal) * toward;
        let face_level = level(model.map(Zone::Flange(kf), p0, z_face));
        let gap = face_level - level(corner(t_lo, 0.0));
        if gap <= EPS {
            return invalid(format!(
                "{}'s tip must stand off {}'s face (folded, touching, they'd be one; gap {gap})",
                ts.tip, ts.face
            ));
        }
        let height = ts.height.unwrap_or(gap + t);
        if !(height > gap) {
            return invalid(format!(
                "tabs of height {height} don't reach the face {gap} away"
            ));
        }
        let (clearance, across_clearance) = (
            ts.clearance.unwrap_or(0.1),
            ts.clearance_across.unwrap_or(0.1),
        );
        // Each slot: the tab's footprint on the face, in the face's flat
        // frame, grown by the clearance.
        for &(a, b) in &ranges {
            let foot: Vec<P2> = [(a, 0.0), (b, 0.0), (b, t), (a, t)]
                .iter()
                .map(|&(tt, z)| {
                    let q = add3(corner(tt, z), scale3(out, gap));
                    let f = model.unmap_flange(kf, q);
                    [f[0], f[1]]
                })
                .collect();
            let centre = scale(add(add(foot[0], foot[1]), add(foot[2], foot[3])), 0.25);
            let (along, across) = (unit(sub(foot[1], foot[0])), unit(sub(foot[3], foot[0])));
            let grown: Vec<P2> = foot
                .iter()
                .map(|&p| {
                    let d = sub(p, centre);
                    let (x, y) = (dot(d, along), dot(d, across));
                    add(
                        centre,
                        add(
                            scale(along, x + clearance * x.signum()),
                            scale(across, y + across_clearance * y.signum()),
                        ),
                    )
                })
                .collect();
            let mut outline = grown;
            if crate::flat::path_area(&outline, &[None, None, None, None]) < 0.0 {
                outline.reverse();
            }
            cut.base.cutouts.push(Cutout {
                outline,
                bulges: Vec::new(),
            });
        }
        bumps.push((kl, s_tip, ranges, height));
    }
    let mut layout = layout_cells(&cut)?;
    for (kl, s_tip, ranges, height) in bumps {
        let fold = layout.folds[kl].clone();
        let at = |tt: f64, s: f64| add(fold.edge_point, add(scale(fold.u, tt), scale(fold.n, s)));
        let s_of = |p: P2| dot(sub(p, fold.edge_point), fold.n);
        let t_of = |p: P2| dot(sub(p, fold.edge_point), fold.u);
        let add_tabs = |pts: &mut Vec<P2>, arcs: &mut Vec<Option<Curve>>| -> bool {
            let mut done = 0;
            for &(a, b) in &ranges {
                let n = pts.len();
                let Some(k) = (0..n).find(|&k| {
                    let (p, q) = (pts[k], pts[(k + 1) % n]);
                    arcs[k].is_none()
                        && (s_of(p) - s_tip).abs() < 1e-7
                        && (s_of(q) - s_tip).abs() < 1e-7
                        && t_of(p).min(t_of(q)) < a + EPS
                        && t_of(p).max(t_of(q)) > b - EPS
                }) else {
                    continue;
                };
                let forward = t_of(pts[(k + 1) % n]) > t_of(pts[k]);
                let (u, v) = if forward { (a, b) } else { (b, a) };
                let bump = [
                    at(u, s_tip),
                    at(u, s_tip + height),
                    at(v, s_tip + height),
                    at(v, s_tip),
                ];
                for (j, p) in bump.into_iter().enumerate() {
                    pts.insert(k + 1 + j, p);
                    arcs.insert(k + 1 + j, None);
                }
                done += 1;
            }
            done == ranges.len()
        };
        let mut placed = false;
        for cell in layout
            .cells
            .iter_mut()
            .filter(|c| c.zone == Zone::Flange(kl))
        {
            placed |= add_tabs(&mut cell.polygon, &mut cell.arcs);
        }
        if !placed || !add_tabs(&mut layout.flat.outline, &mut layout.flat.outline_arcs) {
            return Err(SheetError::Unsupported(
                "tab_slots: tabs on a tip that isn't one straight edge".into(),
            ));
        }
    }
    Ok(layout)
}

fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn scale(a: P2, s: f64) -> P2 {
    [a[0] * s, a[1] * s]
}
fn unit(a: P2) -> P2 {
    scale(a, 1.0 / dot(a, a).sqrt())
}
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn unit3(a: [f64; 3]) -> [f64; 3] {
    scale3(a, 1.0 / dot3(a, a).sqrt())
}

#[cfg(test)]
mod tests {
    use crate::spec::{PartSpec, TabSlot};

    /// A box section: f2 up 40 on y = 0, f1 up 30 on y = 60, and a lip on f1's
    /// tip back across the part, its end `gap` short of f2's inner face.
    fn closed(gap: f64, tabs: Vec<TabSlot>) -> PartSpec {
        let mut part: PartSpec = serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0, 0], [100, 0], [100, 60], [0, 60]] }},
                 "flanges": [
                   {{ "id": "f1", "on": {{ "base": 2 }}, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }},
                   {{ "id": "f2", "on": {{ "base": 0 }}, "angle": 90, "length": 40, "datum": "outer", "position": "outside", "side": "up" }},
                   {{ "id": "lip1", "on": {{ "tip": "f1" }}, "angle": 90, "length": {}, "datum": "outer", "position": "outside", "side": "up" }}] }}"#,
            62.0 - gap
        ))
        .unwrap();
        part.tab_slots = tabs;
        part
    }

    fn tabs(count: usize, length: f64, offsets: (f64, f64)) -> TabSlot {
        TabSlot {
            tip: "lip1".into(),
            face: "f2".into(),
            count,
            length,
            height: None,
            clearance: Some(0.2),
            clearance_across: None,
            offset_start: Some(offsets.0),
            offset_end: Some(offsets.1),
        }
    }

    #[test]
    fn tabs_add_and_slots_take_their_volume() {
        let volume = |p: &PartSpec| crate::fold::exact_mass(p).unwrap().volume;
        let plain = volume(&closed(0.1, Vec::new()));
        let tabbed = closed(0.1, vec![tabs(3, 10.0, (10.0, 10.0))]);
        crate::fold::fold(&tabbed).unwrap();
        // Through by default: height gap + T; slots 0.2 longer each end,
        // 0.1 wider across the sheet each side.
        let expected = plain + 3.0 * 10.0 * 2.1 * 2.0 - 3.0 * 10.4 * 2.2 * 2.0;
        assert!(
            (volume(&tabbed) - expected).abs() < 1e-6,
            "{} {expected}",
            volume(&tabbed)
        );
        let flat = crate::flat::flat_pattern(&tabbed).unwrap();
        assert_eq!(flat.cutouts.len(), 3);
        assert_eq!(flat.outline.len(), 4 + 3 * 4);
    }

    #[test]
    fn one_tab_starts_at_the_start_offset() {
        // The lip's tip runs along f1's base edge, (100, 60) -> (0, 60): its
        // start is at x = 100.
        let flat =
            crate::flat::flat_pattern(&closed(0.1, vec![tabs(1, 20.0, (30.0, 0.0))])).unwrap();
        let top = flat.outline.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
        let xs: Vec<f64> = flat
            .outline
            .iter()
            .filter(|p| (p[1] - top).abs() < 1e-9)
            .map(|p| p[0])
            .collect();
        assert_eq!(xs.len(), 2);
        let (lo, hi) = (xs[0].min(xs[1]), xs[0].max(xs[1]));
        assert!(
            (lo - 50.0).abs() < 1e-9 && (hi - 70.0).abs() < 1e-9,
            "{xs:?}"
        );
    }

    #[test]
    fn a_tip_touching_its_face_or_too_many_tabs_are_refused() {
        let invalid = |p: PartSpec| {
            matches!(
                crate::flat::flat_pattern(&p),
                Err(crate::SheetError::Invalid(_))
            )
        };
        assert!(invalid(closed(0.0, vec![tabs(2, 10.0, (0.0, 0.0))])));
        assert!(invalid(closed(0.1, vec![tabs(6, 20.0, (0.0, 0.0))])));
    }
}
