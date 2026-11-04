//! Folds along a line across the base or a flange (a fold feature). The flat
//! pattern doesn't change: each fold splits the region its bend crosses into
//! what stays, a bend strip and what turns. Flanges and hems wholly on the
//! turning side turn with the fold; a fold line may not run through them.

use crate::flat::{Arc, Cell, Curve, FlatBend, Fold, Layout, Zone, P2};
use crate::spec::{BendLocation, PartSpec, Side};
use crate::SheetError;

const EPS: f64 = 1e-9;

type Piece = (Vec<P2>, Vec<Option<Curve>>);

/// Adds the spec's folds to a layout of the base with its flanges and hems
/// (before holes are cut).
pub(crate) fn apply_folds(layout: &mut Layout, spec: &PartSpec) -> Result<(), SheetError> {
    let unsupported = |what: String| Err(SheetError::Unsupported(what));

    // Each fold as a virtual edge: s along the turning side's normal from the line.
    let mut lines = Vec::new();
    for f in &spec.folds {
        let (a, b) = (f.line[0], f.line[1]);
        let d = sub(b, a);
        let len = dot(d, d).sqrt();
        if len < EPS || !(f.angle > 0.0 && f.angle < 180.0) {
            return Err(SheetError::Invalid(format!(
                "fold {}: needs a line and an angle between 0 and 180",
                f.id
            )));
        }
        let u = [d[0] / len, d[1] / len];
        let r = f.bend_radius.unwrap_or(spec.bend_radius);
        let (bend_allowance, k) = spec.bend_allowance(f.angle, r, f.k_factor)?;
        let s0 = match f.location {
            BendLocation::Centerline => -bend_allowance / 2.0,
            BendLocation::Start => 0.0,
            BendLocation::End => -bend_allowance,
        };
        lines.push(Fold {
            flange: f.id.clone(),
            parent: None,
            edge_point: a,
            u,
            n: [u[1], -u[0]],
            s0,
            bend_allowance,
            angle: f.angle.to_radians(),
            radius: r,
            k_factor: k,
            up: f.side == Side::Up,
        });
    }
    let nf = lines.len();
    let s = |f: &Fold, p: P2| dot(sub(p, f.edge_point), f.n);

    // The flanges' folds go after the fold lines (parents come first).
    for f in &mut layout.folds {
        f.parent = f.parent.map(|p| p + nf);
    }
    for c in &mut layout.cells {
        c.zone = match c.zone {
            Zone::Bend(j) => Zone::Bend(j + nf),
            Zone::Flange(j) => Zone::Flange(j + nf),
            zone => zone,
        };
    }
    let mut folds = lines.clone();
    folds.append(&mut layout.folds);

    // Each fold line splits the one region its bend crosses: the base, a
    // flange (or hem, contour) strip, or the part an earlier fold turns (then
    // it is nested in that fold). Below the bend stays, the strip bends, the
    // rest turns.
    let mut bends = Vec::new();
    for i in 0..nf {
        let mut region = None;
        for c in &layout.cells {
            let foldable = match c.zone {
                Zone::Base => true,
                Zone::Flange(j) => j < i || j >= nf,
                Zone::Bend(_) => false,
            };
            let f = &lines[i];
            let (lo, hi) = extent(c, |p| s(f, p));
            if foldable && lo < f.s0 + f.bend_allowance - EPS && hi > f.s0 + EPS {
                if region.is_some_and(|z| z != c.zone) {
                    return unsupported(format!(
                        "fold {} crosses more than one region of the part",
                        f.flange
                    ));
                }
                region = Some(c.zone);
            }
        }
        let Some(region) = region else {
            return Err(SheetError::Invalid(format!(
                "fold {}: its line must cross the part",
                lines[i].flange
            )));
        };
        if let Zone::Flange(j) = region {
            folds[i].parent = Some(j);
            // On a strip, the side towards the strip's own bend stays: a fold
            // the other way round is the same shape with the other side turning.
            let toward_bend = layout.cells.iter().any(|c| {
                c.zone == Zone::Bend(j) && extent(c, |p| s(&lines[i], p)).1 > lines[i].s0 + EPS
            });
            if toward_bend {
                for f in [&mut lines[i], &mut folds[i]] {
                    f.u = [-f.u[0], -f.u[1]];
                    f.n = [-f.n[0], -f.n[1]];
                    f.s0 = -f.bend_allowance - f.s0;
                }
            }
        }
        let f = &lines[i];
        let mut cells = Vec::new();
        let (mut strips, mut turns) = (Vec::new(), Vec::new());
        for cell in layout.cells.drain(..) {
            if cell.zone != region {
                cells.push(cell);
                continue;
            }
            let (below, rest) = split(&cell.polygon, &cell.arcs, f, f.s0)?;
            for piece in below {
                cells.push(piece_cell(&cell.name, piece, region));
            }
            for piece in rest {
                let (strip, above) = split(&piece.0, &piece.1, f, f.s0 + f.bend_allowance)?;
                strips.extend(strip);
                turns.extend(above);
            }
        }
        if strips.is_empty() || turns.is_empty() {
            return Err(SheetError::Invalid(format!(
                "fold {}: its line must cross the part",
                f.flange
            )));
        }
        if strips.len() > 1 {
            return unsupported(format!(
                "fold {}: a line that crosses the base more than once",
                f.flange
            ));
        }
        bends.push(FlatBend {
            flange: f.flange.clone(),
            line: chord(&strips[0].0, f, f.s0 + f.bend_allowance / 2.0),
            angle: f.angle.to_degrees(),
            inner_radius: f.radius,
            k_factor: f.k_factor,
            up: f.up,
        });
        for (n, piece) in strips.into_iter().enumerate() {
            let label = if n == 0 {
                format!("{}.bend", f.flange)
            } else {
                format!("{}.bend.{n}", f.flange)
            };
            cells.push(piece_cell(&label, piece, Zone::Bend(i)));
        }
        let many = turns.len() > 1;
        for (n, piece) in turns.into_iter().enumerate() {
            let label = if many {
                format!("{}.{n}", f.flange)
            } else {
                f.flange.clone()
            };
            cells.push(piece_cell(&label, piece, Zone::Flange(i)));
        }
        layout.cells = cells;
    }

    // A flange or hem hangs on the deepest fold line whose turning side holds
    // its bend, among those in the region it sits on (the base, or its
    // parent's strip): a fold line's host is its first ancestor that isn't one.
    let depth = |folds: &[Fold], mut j: usize| {
        let mut d = 0;
        while let Some(p) = folds[j].parent {
            d += 1;
            j = p;
        }
        d
    };
    let host = |folds: &[Fold], mut j: usize| loop {
        match folds[j].parent {
            Some(p) if p < nf => j = p,
            other => return other,
        }
    };
    for k in nf..folds.len() {
        let on = folds[k].parent;
        let deepest = layout
            .cells
            .iter()
            .filter(|c| c.zone == Zone::Bend(k))
            .flat_map(|c| {
                lines
                    .iter()
                    .enumerate()
                    .filter(|&(i, f)| {
                        host(&folds, i) == on
                            && extent(c, |p| s(f, p)).0 >= f.s0 + f.bend_allowance - EPS
                    })
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>()
            })
            .max_by_key(|&i| depth(&folds, i));
        if deepest.is_some() {
            folds[k].parent = deepest;
        }
    }

    // Every cell in a fold line's host region lies on the turning side of the
    // folds it hangs on (through its chain) and on the staying side of every
    // other fold.
    let hangs_on = |mut j: usize, i: usize| loop {
        if j == i {
            return true;
        }
        match folds[j].parent {
            Some(p) => j = p,
            None => return false,
        }
    };
    for (i, f) in lines.iter().enumerate() {
        let h = host(&folds, i);
        for c in &layout.cells {
            let owner = match c.zone {
                Zone::Base => None,
                Zone::Bend(j) | Zone::Flange(j) => Some(j),
            };
            if let Some(h) = h {
                if !owner.is_some_and(|j| hangs_on(j, h)) {
                    continue;
                }
            }
            let (lo, hi) = extent(c, |p| s(f, p));
            let ok = match c.zone {
                Zone::Bend(j) if j == i => true,
                Zone::Bend(j) | Zone::Flange(j) if hangs_on(j, i) => {
                    lo >= f.s0 + f.bend_allowance - EPS
                }
                _ => hi <= f.s0 + EPS,
            };
            if !ok {
                return unsupported(format!("fold {} runs through {}", f.flange, c.name));
            }
        }
    }

    bends.append(&mut layout.flat.bends);
    layout.flat.bends = bends;
    // Parents before children (a fold line on a flange comes after it).
    let mut order = Vec::with_capacity(folds.len());
    let mut placed = vec![false; folds.len()];
    while order.len() < folds.len() {
        for j in 0..folds.len() {
            if !placed[j] && folds[j].parent.is_none_or(|p| placed[p]) {
                placed[j] = true;
                order.push(j);
            }
        }
    }
    let mut new = vec![0; folds.len()];
    for (k, &j) in order.iter().enumerate() {
        new[j] = k;
    }
    layout.folds = order
        .iter()
        .map(|&j| Fold {
            parent: folds[j].parent.map(|p| new[p]),
            ..folds[j].clone()
        })
        .collect();
    for c in &mut layout.cells {
        c.zone = match c.zone {
            Zone::Bend(j) => Zone::Bend(new[j]),
            Zone::Flange(j) => Zone::Flange(new[j]),
            zone => zone,
        };
    }
    Ok(())
}

fn piece_cell(name: &str, (polygon, arcs): Piece, zone: Zone) -> Cell {
    Cell {
        name: name.into(),
        polygon,
        arcs,
        holes: Vec::new(),
        inner: Vec::new(),
        zone,
    }
}

/// Lowest and highest value of `f` over a cell's boundary (arcs sampled).
fn extent(c: &Cell, f: impl Fn(P2) -> f64) -> (f64, f64) {
    let mut values = Vec::new();
    let n = c.polygon.len();
    for k in 0..n {
        let (p, q) = (c.polygon[k], c.polygon[(k + 1) % n]);
        values.push(f(p));
        if let Some(arc) = c.arcs[k] {
            for x in curve_points(p, q, arc, 16) {
                values.push(f(x));
            }
        }
    }
    values
        .iter()
        .fold((f64::MAX, f64::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)))
}

fn arc_points(p: P2, q: P2, arc: Arc, n: usize) -> Vec<P2> {
    let c = arc.center;
    let r = dot(sub(p, c), sub(p, c)).sqrt();
    let a0 = (p[1] - c[1]).atan2(p[0] - c[0]);
    let sweep = arc.sweep(p, q);
    (1..n)
        .map(|i| {
            let a = a0 + sweep * i as f64 / n as f64;
            [c[0] + r * a.cos(), c[1] + r * a.sin()]
        })
        .collect()
}

/// Splits a counter-clockwise boundary of lines and arcs by the line s = c of
/// fold f into its pieces below and above. The line must cross straight
/// segments, through no vertex.
fn split(
    pts: &[P2],
    arcs: &[Option<Curve>],
    f: &Fold,
    c: f64,
) -> Result<(Vec<Piece>, Vec<Piece>), SheetError> {
    let s = |p: P2| dot(sub(p, f.edge_point), f.n) - c;
    let unsupported = |what: &str| {
        Err(SheetError::Unsupported(format!(
            "fold {}: {what}",
            f.flange
        )))
    };
    let n = pts.len();
    // Crossings: (segment, position along it, point, going up, along the line).
    let mut crossings: Vec<(usize, f64, P2, bool, f64)> = Vec::new();
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        let (sp, sq) = (s(p), s(q));
        if sp.abs() < EPS {
            return unsupported("its line runs through a corner of the part");
        }
        if let Some(arc) = arcs[k] {
            if curve_points(p, q, arc, 32)
                .into_iter()
                .any(|x| (s(x) > 0.0) != (sp > 0.0) || (s(x) > 0.0) != (sq > 0.0))
            {
                return unsupported("its line runs through a curved edge");
            }
            continue;
        }
        if (sp > 0.0) != (sq > 0.0) {
            let t = sp / (sp - sq);
            let x = [p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])];
            crossings.push((k, t, x, sq > 0.0, dot(x, f.u)));
        }
    }
    if crossings.is_empty() {
        let whole = (pts.to_vec(), arcs.to_vec());
        return Ok(if s(pts[0]) < 0.0 {
            (vec![whole], Vec::new())
        } else {
            (Vec::new(), vec![whole])
        });
    }
    // Along the line, the inside of the boundary is between crossings 0-1,
    // 2-3, ...: each crossing's partner.
    let mut order: Vec<usize> = (0..crossings.len()).collect();
    order.sort_by(|&a, &b| crossings[a].4.total_cmp(&crossings[b].4));
    let mut partner = vec![0; crossings.len()];
    for pair in order.chunks(2) {
        if pair.len() < 2 {
            return unsupported("its crossings with the part don't pair up");
        }
        partner[pair[0]] = pair[1];
        partner[pair[1]] = pair[0];
    }
    // Pieces on one side: from a crossing where the boundary enters that
    // side, along the boundary to where it leaves, then along the line to
    // that crossing's partner, and so on.
    let pieces = |up: bool| -> Result<Vec<Piece>, SheetError> {
        let mut visited = vec![false; crossings.len()];
        let mut out = Vec::new();
        while let Some(start) = (0..crossings.len()).find(|&i| !visited[i] && crossings[i].3 == up)
        {
            let (mut ps, mut ars) = (Vec::new(), Vec::new());
            let mut cur = start;
            loop {
                visited[cur] = true;
                let x = crossings[cur];
                ps.push(x.2);
                ars.push(None);
                let next = (cur + 1) % crossings.len();
                let e = crossings[next];
                if !(e.0 == x.0 && e.1 > x.1) {
                    let mut k = x.0;
                    loop {
                        k = (k + 1) % n;
                        ps.push(pts[k]);
                        ars.push(arcs[k]);
                        if k == e.0 {
                            break;
                        }
                    }
                }
                visited[next] = true;
                ps.push(e.2);
                ars.push(None);
                cur = partner[next];
                if cur == start {
                    break;
                }
                if visited[cur] {
                    return Err(SheetError::Unsupported(format!(
                        "fold {}: its line splits the part ambiguously",
                        f.flange
                    )));
                }
            }
            out.push((ps, ars));
        }
        Ok(out)
    };
    Ok((pieces(false)?, pieces(true)?))
}

/// The bend line: where s = c crosses the strip, end to end.
fn chord(poly: &[P2], f: &Fold, c: f64) -> [P2; 2] {
    let s = |p: P2| dot(sub(p, f.edge_point), f.n) - c;
    let mut hits = Vec::new();
    for i in 0..poly.len() {
        let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
        let (sp, sq) = (s(p), s(q));
        if (sp <= 0.0) != (sq <= 0.0) {
            let t = sp / (sp - sq);
            hits.push([p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])]);
        }
    }
    let along = |p: &P2| dot(*p, f.u);
    let lo = hits
        .iter()
        .min_by(|a, b| along(a).total_cmp(&along(b)))
        .unwrap();
    let hi = hits
        .iter()
        .max_by(|a, b| along(a).total_cmp(&along(b)))
        .unwrap();
    [*lo, *hi]
}

fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[cfg(test)]
mod tests {
    use crate::flat::layout;
    use crate::fold::{fold, mass};
    use crate::spec::PartSpec;

    fn part(location: &str, side: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "folds": [{{ "id": "f", "line": [[100,40],[0,40]], "angle": 90,
                             "location": "{location}", "side": "{side}" }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_fold_keeps_the_flat_and_bends_on_its_line() {
        // Line y = 40, the +y side turns up: bend axis at
        // y 37.80 / 40 / 35.60 for centreline / start / end; volume 12062.8319.
        for (location, axis) in [("centerline", 37.8009), ("start", 40.0), ("end", 35.6018)] {
            let spec = part(location, "up");
            let l = layout(&spec).unwrap();
            assert_eq!(l.flat.outline, spec.base.outline);
            // n points to +y here, so the bend starts at y = 40 + s0.
            assert!((40.0 + l.folds[0].s0 - axis).abs() < 1e-4, "{location}");
            let m = mass(&fold(&spec).unwrap(), 1e-5).unwrap();
            assert!(
                (m.volume - 12062.8319).abs() < 1e-2,
                "{location}: {}",
                m.volume
            );
        }
    }

    #[test]
    fn a_fold_line_may_cross_the_outline_at_an_angle() {
        // The bend's ends are slanted cuts: helices and wrapped planes.
        for (line, side, holes) in [
            ("[[100,30],[0,50]]", "up", ""),
            ("[[100,30],[0,50]]", "down", ""),
            (
                "[[0,10],[100,45]]",
                "up",
                r#"{ "center": [50, 27.5], "diameter": 7 }"#,
            ),
        ] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "holes": [{holes}] }},
                     "folds": [{{ "id": "f", "line": {line}, "angle": 90,
                                 "location": "centerline", "side": "{side}" }}] }}"#
            ))
            .unwrap();
            let brep = fold(&spec).unwrap_or_else(|e| panic!("{line} {side}: {e}"));
            let m = mass(&brep, 1e-3).unwrap();
            let e = crate::fold::exact_mass(&spec).unwrap();
            assert!(
                ((m.volume - e.volume) / e.volume).abs() < 2e-5,
                "{line} {side}: mesh {} exact {}",
                m.volume,
                e.volume
            );
            assert!(((m.area - e.area) / e.area).abs() < 2e-4, "{line} {side}");
        }
    }

    fn with_flange(fold_line: &str, edge: usize, outline: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": {outline} }},
                 "flanges": [{{ "id": "f1", "on": {{ "base": {edge} }}, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }}],
                 "folds": [{{ "id": "f", "line": {fold_line}, "angle": 90,
                             "location": "centerline", "side": "up" }}] }}"#
        ))
        .unwrap()
    }

    const RECT: &str = "[[0,0],[100,0],[100,60],[0,60]]";

    #[test]
    fn a_flange_on_the_turning_side_turns_with_the_fold() {
        // Expected: 17805.3096 mm3 with the flange on either side.
        for edge in [2, 0] {
            let spec = with_flange("[[100,40],[0,40]]", edge, RECT);
            let l = layout(&spec).unwrap();
            // The fold line comes first; the flange hangs on it only when it turns.
            assert_eq!(l.folds[0].flange, "f");
            assert_eq!(l.folds[1].parent, (edge == 2).then_some(0));
            let m = crate::fold::exact_mass(&spec).unwrap();
            assert!((m.volume - 17805.3096).abs() < 1e-3, "{}", m.volume);
            fold(&spec).unwrap();
        }
    }

    #[test]
    fn a_fold_line_on_a_flange_folds_the_flange() {
        // Across the 30 mm flange on edge 0, 20 mm up (flat y = -18.3982;
        // 15 mm up for the fold down): whichever way the line runs, the tip
        // side turns. The flat stays; a lip on the tip turns along.
        let (y, y2) = (-18.398229715, -13.398229715);
        let lip = r#"{ "id": "l", "on": { "tip": "f1" }, "angle": 90, "length": 10,
                       "datum": "outer", "position": "outside", "side": "up" }"#;
        for (line, location, side, with_lip, volume) in [
            (
                format!("[[0,{y}],[100,{y}]]"),
                "centerline",
                "up",
                false,
                17805.3096,
            ),
            (
                format!("[[100,{y}],[0,{y}]]"),
                "centerline",
                "up",
                false,
                17805.3096,
            ),
            (
                format!("[[100,{y2}],[0,{y2}]]"),
                "start",
                "down",
                false,
                17805.3096,
            ),
            (
                format!("[[0,{y}],[100,{y}]]"),
                "centerline",
                "up",
                true,
                19547.7874,
            ),
            (
                "[[70,-28.398229715],[100,-10]]".into(),
                "centerline",
                "up",
                false,
                17764.5897,
            ),
        ] {
            let what = format!("{line} {location} {side} lip {with_lip}");
            let mut spec = with_flange(&line, 0, RECT);
            spec.folds[0].location = serde_json::from_str(&format!("\"{location}\"")).unwrap();
            spec.folds[0].side = serde_json::from_str(&format!("\"{side}\"")).unwrap();
            if with_lip {
                spec.flanges.push(serde_json::from_str(lip).unwrap());
            }
            let mut plain = spec.clone();
            plain.folds.clear();
            let l = layout(&spec).unwrap_or_else(|e| panic!("{what}: {e}"));
            assert_eq!(
                l.flat.outline,
                layout(&plain).unwrap().flat.outline,
                "{what}"
            );
            // The fold line hangs on the flange, after it; the lip on the fold line.
            let at = |id: &str| l.folds.iter().position(|f| f.flange == id).unwrap();
            assert_eq!(l.folds[at("f")].parent, Some(at("f1")), "{what}");
            assert!(at("f1") < at("f"), "{what}");
            if with_lip {
                assert_eq!(l.folds[at("l")].parent, Some(at("f")), "{what}");
            }
            let e = crate::fold::exact_mass(&spec).unwrap();
            assert!((e.volume - volume).abs() < 1e-3, "{what}: {}", e.volume);
            // The mesh falls short by ~1e-5 per bend at this deflection.
            let m = mass(&fold(&spec).unwrap(), 1e-3).unwrap();
            assert!(
                ((m.volume - e.volume) / e.volume).abs() < 5e-5,
                "{what}: mesh {}",
                m.volume
            );
        }
    }

    #[test]
    fn a_fold_line_may_not_run_through_a_flange() {
        // Edge 1 (x = 100) carries the flange; the line crosses it. The flange
        // would have to fold along the line too, which a sheet can't.
        let spec = with_flange("[[100,40],[0,40]]", 1, RECT);
        assert!(matches!(
            layout(&spec),
            Err(crate::SheetError::Unsupported(_))
        ));
    }

    #[test]
    fn a_fold_splits_a_non_convex_base() {
        // An L; the line crosses its upper leg only. Expected: 7825.1327 mm3.
        let mut spec = with_flange(
            "[[40,40],[0,40]]",
            0,
            "[[0,0],[100,0],[100,25],[40,25],[40,60],[0,60]]",
        );
        spec.flanges.clear();
        let m = crate::fold::exact_mass(&spec).unwrap();
        assert!((m.volume - 7825.1327).abs() < 1e-3, "{}", m.volume);
        fold(&spec).unwrap();
    }

    #[test]
    fn a_fold_in_the_part_another_fold_turns_hangs_on_it() {
        // Expected: 12125.6637 (C and Z), 14868.1415 with a flange beyond both.
        for (side, flange, volume) in [
            ("up", "", 12125.6637),
            ("down", "", 12125.6637),
            (
                "up",
                r#"{ "id": "f3", "on": { "base": 2 }, "angle": 90, "length": 15,
                     "datum": "outer", "position": "outside", "side": "up" }"#,
                14868.1415,
            ),
        ] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "flanges": [{flange}],
                     "folds": [
                       {{ "id": "f1", "line": [[100,20],[0,20]], "angle": 90, "location": "centerline", "side": "up" }},
                       {{ "id": "f2", "line": [[100,45],[0,45]], "angle": 90, "location": "centerline", "side": "{side}" }} ] }}"#
            ))
            .unwrap();
            let l = layout(&spec).unwrap();
            assert_eq!(l.folds[1].parent, Some(0));
            if !flange.is_empty() {
                assert_eq!(
                    l.folds[2].parent,
                    Some(1),
                    "the flange hangs on the inner fold"
                );
            }
            let m = crate::fold::exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{side}: {}", m.volume);
            fold(&spec).unwrap();
        }
    }
}

/// Points strictly inside a curved segment.
fn curve_points(p: P2, q: P2, curve: Curve, pieces: usize) -> Vec<P2> {
    match curve {
        Curve::Arc(arc) => arc_points(p, q, arc, pieces),
        Curve::Cut(c) => {
            let line = c.points(p, q, pieces);
            line[1..line.len() - 1].to_vec()
        }
    }
}
