//! Flat pattern as an ASCII DXF (R12): the blank outline, holes and bend
//! lines, on the layers sheet metal flat patterns commonly use
//! (Outer, Inner, BendUp, BendDown). Millimetres, base frame.

use gusset_sheet::flat::{lines_and_arcs, FlatPattern};
use std::fmt::Write;

pub fn flat_dxf(flat: &FlatPattern) -> String {
    let mut out = String::new();
    let mut pair = |code: u32, value: &str| writeln!(out, "{code}\n{value}").unwrap();
    pair(0, "SECTION");
    pair(2, "HEADER");
    pair(9, "$INSUNITS");
    pair(70, "4"); // millimetres
    pair(0, "ENDSEC");
    pair(0, "SECTION");
    pair(2, "ENTITIES");

    let (outline, arcs) = lines_and_arcs(&flat.outline, &flat.outline_arcs);
    ring(&mut pair, "Outer", &outline, &arcs);
    for cutout in &flat.cutouts {
        let (points, arcs) = lines_and_arcs(&cutout.points, &cutout.arcs);
        ring(&mut pair, "Inner", &points, &arcs);
    }
    for hole in &flat.holes {
        pair(0, "CIRCLE");
        pair(8, "Inner");
        pair(10, &num(hole.center[0]));
        pair(20, &num(hole.center[1]));
        pair(30, "0.0");
        pair(40, &num(hole.diameter / 2.0));
    }
    for bend in &flat.bends {
        let layer = if bend.up { "BendUp" } else { "BendDown" };
        line(&mut pair, layer, bend.line[0], bend.line[1]);
    }
    pair(0, "ENDSEC");
    pair(0, "EOF");
    out
}

/// A closed boundary of lines and arcs.
fn ring(
    pair: &mut impl FnMut(u32, &str),
    layer: &str,
    points: &[[f64; 2]],
    arcs: &[Option<gusset_sheet::flat::Arc>],
) {
    let n = points.len();
    for i in 0..n {
        let (p, q) = (points[i], points[(i + 1) % n]);
        match arcs[i] {
            None => line(pair, layer, p, q),
            Some(arc) => {
                // DXF arcs run counter-clockwise from start to end angle.
                let (s, e) = if arc.ccw { (p, q) } else { (q, p) };
                let c = arc.center;
                let angle = |v: [f64; 2]| (v[1] - c[1]).atan2(v[0] - c[0]).to_degrees();
                pair(0, "ARC");
                pair(8, layer);
                pair(10, &num(c[0]));
                pair(20, &num(c[1]));
                pair(30, "0.0");
                pair(
                    40,
                    &num(((p[0] - c[0]).powi(2) + (p[1] - c[1]).powi(2)).sqrt()),
                );
                pair(50, &num(angle(s)));
                pair(51, &num(angle(e)));
            }
        }
    }
}

fn line(pair: &mut impl FnMut(u32, &str), layer: &str, a: [f64; 2], b: [f64; 2]) {
    pair(0, "LINE");
    pair(8, layer);
    pair(10, &num(a[0]));
    pair(20, &num(a[1]));
    pair(30, "0.0");
    pair(11, &num(b[0]));
    pair(21, &num(b[1]));
    pair(31, "0.0");
}

/// Shortest text that reads back as exactly the same f64.
fn num(v: f64) -> String {
    format!("{v:?}")
}

#[cfg(test)]
mod tests {
    use super::*;
    use gusset_sheet::flat::flat_pattern;
    use gusset_sheet::spec::PartSpec;

    #[test]
    fn a_flanged_plate_reads_back_with_its_outline_holes_and_bend() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                           "holes": [{ "center": [20, 30], "diameter": 8 }] },
                 "flanges": [{ "id": "f1", "on": { "base": 2 }, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }] }"#,
        )
        .unwrap();
        let flat = flat_pattern(&spec).unwrap();
        let text = flat_dxf(&flat);
        let drawing = ::dxf::Drawing::load(&mut text.as_bytes()).unwrap();

        let mut layers: Vec<(String, String)> = drawing
            .entities()
            .map(|e| {
                let kind = match &e.specific {
                    ::dxf::entities::EntityType::Line(_) => "line",
                    ::dxf::entities::EntityType::Circle(_) => "circle",
                    _ => "other",
                };
                (e.common.layer.clone(), kind.to_string())
            })
            .collect();
        layers.sort();
        let count = |layer: &str, kind: &str| {
            layers
                .iter()
                .filter(|(l, k)| l == layer && k == kind)
                .count()
        };
        assert_eq!(count("Outer", "line"), 4);
        assert_eq!(count("Inner", "circle"), 1);
        assert_eq!(count("BendUp", "line"), 1);

        // Coordinates survive exactly: the far edge sits at the flat width.
        let far = drawing.entities().find_map(|e| match &e.specific {
            ::dxf::entities::EntityType::Line(l) if e.common.layer == "Outer" && l.p1.y > 80.0 => {
                Some(l.p1.y)
            }
            _ => None,
        });
        assert_eq!(far, Some(flat.outline[2][1]));
        assert!((flat.outline[2][1] - 88.3982297150257).abs() < 1e-9);
    }
}

#[cfg(test)]
mod arc_tests {
    use super::*;
    use gusset_sheet::flat::flat_pattern;
    use gusset_sheet::spec::PartSpec;

    #[test]
    fn round_reliefs_are_written_as_arcs() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "relief": { "shape": "round", "width": 2, "depth": 1, "minimum_remnant": 4 },
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
                 "flanges": [{ "id": "f1", "on": { "base": 2 }, "width": { "offset_start": 30, "offset_end": 20 },
                               "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }] }"#,
        )
        .unwrap();
        let text = flat_dxf(&flat_pattern(&spec).unwrap());
        let drawing = ::dxf::Drawing::load(&mut text.as_bytes()).unwrap();
        let arcs: Vec<_> = drawing
            .entities()
            .filter_map(|e| match &e.specific {
                ::dxf::entities::EntityType::Arc(a) => Some(a.clone()),
                _ => None,
            })
            .collect();
        assert_eq!(arcs.len(), 2);
        for a in &arcs {
            assert!((a.radius - 1.0).abs() < 1e-12);
            // Half circles that bulge into the base (toward -y): 180 -> 360.
            let sweep = (a.end_angle - a.start_angle).rem_euclid(360.0);
            assert!(
                (sweep - 180.0).abs() < 1e-9
                    && (a.start_angle.rem_euclid(360.0) - 180.0).abs() < 1e-9,
                "{a:?}"
            );
        }
    }
}
