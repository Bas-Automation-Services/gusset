//! Flat pattern as SVG, for display: the blank outline and cut-outs as paths
//! (lines and arcs), holes as circles, bend lines solid (up) or dashed
//! (down). Millimetres; y is flipped so the base frame's +y points up.
//! Classes for styling: `cut`, `hole`, `bend up` / `bend down`, `clash`.

use gusset_sheet::flat::{Arc, FlatPattern};
use std::fmt::Write;

pub fn flat_svg(flat: &FlatPattern) -> String {
    flat_svg_marked(flat, &[])
}

/// With areas marked in red under the bend lines: where a part folds into
/// itself (`gusset_sheet::clash`).
pub fn flat_svg_marked(flat: &FlatPattern, marks: &[Vec<[f64; 2]>]) -> String {
    let mut pts: Vec<[f64; 2]> = flat.outline.clone();
    // Arcs can bulge past their end points: take their full circles' boxes.
    for (k, arc) in flat.outline_arcs.iter().enumerate() {
        if let Some(a) = arc {
            let r = dist(flat.outline[k], a.center);
            pts.push([a.center[0] - r, a.center[1] - r]);
            pts.push([a.center[0] + r, a.center[1] + r]);
        }
    }
    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
    for p in &pts {
        for j in 0..2 {
            lo[j] = lo[j].min(p[j]);
            hi[j] = hi[j].max(p[j]);
        }
    }
    let margin = 0.05 * (hi[0] - lo[0]).max(hi[1] - lo[1]);
    let mut out = String::new();
    writeln!(
        out,
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="{} {} {} {}">"#,
        num(lo[0] - margin),
        num(-hi[1] - margin),
        num(hi[0] - lo[0] + 2.0 * margin),
        num(hi[1] - lo[1] + 2.0 * margin)
    )
    .unwrap();
    let mut d = ring(&flat.outline, &flat.outline_arcs);
    for c in &flat.cutouts {
        d += &ring(&c.points, &c.arcs);
    }
    writeln!(
        out,
        r#"<path class="cut" d="{d}" fill-rule="evenodd" vector-effect="non-scaling-stroke"/>"#
    )
    .unwrap();
    for h in &flat.holes {
        writeln!(
            out,
            r#"<circle class="hole" cx="{}" cy="{}" r="{}" vector-effect="non-scaling-stroke"/>"#,
            num(h.center[0]),
            num(-h.center[1]),
            num(h.diameter / 2.0)
        )
        .unwrap();
    }
    if !marks.is_empty() {
        let d: String = marks
            .iter()
            .filter(|m| m.len() > 2)
            .map(|m| ring(m, &vec![None; m.len()]))
            .collect();
        writeln!(
            out,
            r##"<path class="clash" d="{d}" fill="#d22" fill-opacity="0.75" stroke="#d22" stroke-width="2" stroke-linejoin="round" vector-effect="non-scaling-stroke"/>"##
        )
        .unwrap();
    }
    for b in &flat.bends {
        writeln!(
            out,
            r#"<line class="bend {}" x1="{}" y1="{}" x2="{}" y2="{}" vector-effect="non-scaling-stroke"/>"#,
            if b.up { "up" } else { "down" },
            num(b.line[0][0]),
            num(-b.line[0][1]),
            num(b.line[1][0]),
            num(-b.line[1][1])
        )
        .unwrap();
    }
    out.push_str("</svg>\n");
    out
}

/// A closed boundary as path data, y flipped.
fn ring(points: &[[f64; 2]], arcs: &[Option<Arc>]) -> String {
    let n = points.len();
    let mut d = format!("M{} {}", num(points[0][0]), num(-points[0][1]));
    for k in 0..n {
        let (p, q) = (points[k], points[(k + 1) % n]);
        match arcs[k] {
            None => write!(d, "L{} {}", num(q[0]), num(-q[1])).unwrap(),
            Some(a) => {
                // The sweep from p to q around the centre, in the arc's sense.
                let (u, v) = (sub(p, a.center), sub(q, a.center));
                let mut sweep = (u[0] * v[1] - u[1] * v[0]).atan2(u[0] * v[0] + u[1] * v[1]);
                if !a.ccw {
                    sweep = -sweep;
                }
                if sweep < 0.0 {
                    sweep += std::f64::consts::TAU;
                }
                let r = dist(p, a.center);
                // Counter-clockwise with y up is clockwise on screen: sweep flag 1.
                write!(
                    d,
                    "A{} {} 0 {} {} {} {}",
                    num(r),
                    num(r),
                    u8::from(sweep > std::f64::consts::PI),
                    u8::from(a.ccw),
                    num(q[0]),
                    num(-q[1])
                )
                .unwrap();
            }
        }
    }
    d.push('Z');
    d
}

fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

fn dist(a: [f64; 2], b: [f64; 2]) -> f64 {
    let d = sub(a, b);
    (d[0] * d[0] + d[1] * d[1]).sqrt()
}

fn num(x: f64) -> String {
    let s = format!("{x:.4}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" {
        "0".into()
    } else {
        s.into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use gusset_sheet::{flat::flat_pattern, spec::PartSpec};

    #[test]
    fn a_rounded_plate_with_a_flange_draws_arcs_holes_and_a_bend() {
        let q = (std::f64::consts::PI / 8.0).tan();
        let spec: PartSpec = serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[10,0],[90,0],[100,10],[100,60],[0,60],[0,10]],
                           "bulges": [0,{q},0,0,0,{q}],
                           "holes": [{{ "center": [50, 30], "diameter": 10 }}] }},
                 "flanges": [{{ "id": "f", "on": {{ "base": 3 }}, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "down" }}] }}"#
        ))
        .unwrap();
        let svg = flat_svg(&flat_pattern(&spec).unwrap());
        // Two R10 quarter arcs, counter-clockwise (sweep flag 1), not large.
        assert_eq!(svg.matches("A10 10 0 0 1").count(), 2, "{svg}");
        assert!(svg.contains(r#"<circle class="hole" cx="50" cy="-30" r="5""#));
        assert!(svg.contains(r#"class="bend down""#));
        // The flange unfolds past y = 60: the view reaches above it.
        assert!(svg.contains("viewBox="));
    }

    #[test]
    fn clashing_areas_are_marked_under_the_bend_lines() {
        let flat = flat_pattern(
            &serde_json::from_str::<PartSpec>(
                r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] } }"#,
            )
            .unwrap(),
        )
        .unwrap();
        assert!(!flat_svg(&flat).contains("clash"));
        let svg = flat_svg_marked(&flat, &[vec![[10.0, 10.0], [20.0, 10.0], [20.0, 15.0]]]);
        assert!(
            svg.contains(r#"<path class="clash" d="M10 -10L20 -10L20 -15L10 -10Z""#),
            "{svg}"
        );
    }
}
