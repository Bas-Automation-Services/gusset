//! Manufacturability checks: what a press brake and a punch or laser shop
//! usually ask of a part, as rules of thumb in multiples of the thickness.
//! Shops differ; every limit is a field of `Rules`.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::flat::{layout, Curve, Zone};
use crate::spec::PartSpec;
use crate::SheetError;

type P2 = [f64; 2];

/// Limits, in multiples of the sheet thickness.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Rules {
    /// Smallest inner bend radius (hems excepted). Default 1.
    pub min_radius: f64,
    /// Shortest straight leg past a bend, for the press brake's die. Default 3.
    pub min_flange: f64,
    /// Closest a hole or cut-out may come to a bend zone in the flat, before
    /// it distorts. Default 2.
    pub hole_to_bend: f64,
    /// Closest a hole or cut-out may come to the blank's outline. Default 1.
    pub hole_to_edge: f64,
    /// Smallest hole diameter, for punching. Default 1.
    pub min_hole: f64,
}

impl Default for Rules {
    fn default() -> Self {
        Rules {
            min_radius: 1.0,
            min_flange: 3.0,
            hole_to_bend: 2.0,
            hole_to_edge: 1.0,
            min_hole: 1.0,
        }
    }
}

/// One rule a part breaks.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Finding {
    pub rule: Rule,
    /// The feature: a flange, contour, fold or hem id, "hole 1" or "cut-out 2".
    pub feature: String,
    /// What the part has, and the limit (mm).
    pub value: f64,
    pub limit: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Rule {
    BendRadius,
    FlangeLength,
    HoleToBend,
    HoleToEdge,
    HoleSize,
}

impl std::fmt::Display for Finding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let what = match self.rule {
            Rule::BendRadius => "inner bend radius",
            Rule::FlangeLength => "straight leg past the bend",
            Rule::HoleToBend => "distance to a bend zone",
            Rule::HoleToEdge => "distance to the outline",
            Rule::HoleSize => "diameter",
        };
        write!(
            f,
            "{}: {what} {:.2} mm, at least {:.2} mm",
            self.feature, self.value, self.limit
        )
    }
}

/// Every rule the part breaks, in the order above.
pub fn check(spec: &PartSpec, rules: &Rules) -> Result<Vec<Finding>, SheetError> {
    let spec = &*spec.resolved()?;
    let t = spec.thickness;
    let mut out = Vec::new();
    let mut find = |rule: Rule, feature: String, value: f64, limit: f64| {
        if value < limit - 1e-9 {
            out.push(Finding {
                rule,
                feature,
                // No "-0.00" for a hole touching a bend.
                value: value + 0.0,
                limit,
            });
        }
    };

    // Bends: flanges, contours (hems are tight on purpose) and folds.
    let hem = |id: &str| {
        spec.hems
            .iter()
            .any(|h| id == h.id || id.starts_with(&format!("{}.", h.id)))
    };
    let steps = crate::hem::steps(spec)?;
    for s in steps.iter().filter(|s| !hem(&s.id)) {
        find(
            Rule::BendRadius,
            s.id.clone(),
            s.radius,
            rules.min_radius * t,
        );
    }
    for f in &spec.folds {
        let r = f.bend_radius.unwrap_or(spec.bend_radius);
        find(Rule::BendRadius, f.id.clone(), r, rules.min_radius * t);
    }
    // A leg's straight length; one with a bend right after it is a short leg too.
    for s in steps.iter().filter(|s| !hem(&s.id)) {
        let tip_has_hem = spec
            .hems
            .iter()
            .any(|h| matches!(&h.on, crate::spec::EdgeRef::Tip(p) if *p == s.id));
        if !tip_has_hem {
            find(
                Rule::FlangeLength,
                s.id.clone(),
                s.straight,
                rules.min_flange * t,
            );
        }
    }

    // Holes and cut-outs, in the flat.
    let l = layout(spec)?;
    let outline = polyline(&l.flat.outline, &l.flat.outline_arcs);
    let bends: Vec<(String, Vec<P2>)> = l
        .cells
        .iter()
        .filter_map(|c| match c.zone {
            Zone::Bend(k) => Some((l.folds[k].flange.clone(), polyline(&c.polygon, &c.arcs))),
            _ => None,
        })
        .collect();
    for (n, h) in l.holes.iter().enumerate() {
        let name = format!("hole {}", n + 1);
        let r = h.diameter / 2.0;
        find(Rule::HoleSize, name.clone(), h.diameter, rules.min_hole * t);
        let to_bend = bends
            .iter()
            .map(|(_, b)| signed(h.center, b) - r)
            .fold(f64::INFINITY, f64::min);
        find(
            Rule::HoleToBend,
            name.clone(),
            to_bend,
            rules.hole_to_bend * t,
        );
        find(
            Rule::HoleToEdge,
            name,
            distance(h.center, &outline) - r,
            rules.hole_to_edge * t,
        );
    }
    for (n, c) in l.flat.cutouts.iter().enumerate() {
        let name = format!("cut-out {}", n + 1);
        let ring = polyline(&c.points, &c.arcs);
        let to_bend = bends
            .iter()
            .map(|(_, b)| between(&ring, b))
            .fold(f64::INFINITY, f64::min);
        find(
            Rule::HoleToBend,
            name.clone(),
            to_bend,
            rules.hole_to_bend * t,
        );
        find(
            Rule::HoleToEdge,
            name,
            clearance(&ring, &outline),
            rules.hole_to_edge * t,
        );
    }
    Ok(out)
}

/// A closed boundary as points, arcs in pieces of at most 3 degrees.
fn polyline(points: &[P2], arcs: &[Option<Curve>]) -> Vec<P2> {
    let n = points.len();
    let mut out = Vec::new();
    for k in 0..n {
        let (p, q) = (points[k], points[(k + 1) % n]);
        out.push(p);
        if let Some(Curve::Cut(c)) = arcs.get(k).copied().flatten() {
            let line = c.points(p, q, 64);
            out.extend(&line[1..line.len() - 1]);
            continue;
        }
        if let Some(a) = arcs.get(k).copied().flatten().and_then(Curve::arc) {
            let (u, v) = (sub(p, a.center), sub(q, a.center));
            let mut sweep = (u[0] * v[1] - u[1] * v[0]).atan2(u[0] * v[0] + u[1] * v[1]);
            if !a.ccw {
                sweep = -sweep;
            }
            if sweep < 0.0 {
                sweep += std::f64::consts::TAU;
            }
            let sweep = if a.ccw { sweep } else { -sweep };
            let pieces = (sweep.abs() / 3f64.to_radians()).ceil().max(1.0) as usize;
            for i in 1..pieces {
                let (s, c) = (sweep * i as f64 / pieces as f64).sin_cos();
                out.push([
                    a.center[0] + u[0] * c - u[1] * s,
                    a.center[1] + u[0] * s + u[1] * c,
                ]);
            }
        }
    }
    out
}

fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}

/// From x to the closed polyline.
fn distance(x: P2, poly: &[P2]) -> f64 {
    let n = poly.len();
    (0..n)
        .map(|k| {
            let (a, b) = (poly[k], poly[(k + 1) % n]);
            let (d, w) = (sub(b, a), sub(x, a));
            let t = ((w[0] * d[0] + w[1] * d[1]) / (d[0] * d[0] + d[1] * d[1]).max(1e-300))
                .clamp(0.0, 1.0);
            let e = sub(w, [d[0] * t, d[1] * t]);
            (e[0] * e[0] + e[1] * e[1]).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

/// Negative inside the polygon.
fn signed(x: P2, poly: &[P2]) -> f64 {
    let n = poly.len();
    let mut inside = false;
    for k in 0..n {
        let (a, b) = (poly[k], poly[(k + 1) % n]);
        if (a[1] > x[1]) != (b[1] > x[1])
            && x[0] < a[0] + (x[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            inside = !inside;
        }
    }
    let d = distance(x, poly);
    if inside {
        -d
    } else {
        d
    }
}

/// Closest approach of two boundaries.
fn clearance(a: &[P2], b: &[P2]) -> f64 {
    let d = |x: &[P2], y: &[P2]| {
        x.iter()
            .map(|&p| distance(p, y))
            .fold(f64::INFINITY, f64::min)
    };
    d(a, b).min(d(b, a))
}

/// Gap between two regions; 0 where they overlap.
fn between(a: &[P2], b: &[P2]) -> f64 {
    if a.iter().any(|&p| signed(p, b) < 0.0) || b.iter().any(|&p| signed(p, a) < 0.0) {
        return 0.0;
    }
    clearance(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(extra: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, {extra} }}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_sound_part_passes() {
        let spec = part(
            r#""base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                         "holes": [{ "center": [50, 30], "diameter": 8 }] },
               "flanges": [{ "id": "f", "on": { "base": 2 }, "angle": 90, "length": 30,
                             "datum": "outer", "position": "outside", "side": "up" }]"#,
        );
        assert_eq!(check(&spec, &Rules::default()).unwrap(), vec![]);
    }

    #[test]
    fn each_rule_names_the_feature_value_and_limit() {
        // T 2. The flange: R 1 (< 2), 8 long, so 8 - (R + T) = 5 straight
        // (< 6); its bend starts at y = 59. Hole 1: 1.5 across (< 2), its
        // edge 1 from the outline (< 2). Hole 2: its edge at y = 58, 1 from
        // the bend (< 4). The cut-out crosses the bend.
        let spec = part(
            r#""base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                         "holes": [{ "center": [10, 1.75], "diameter": 1.5 },
                                   { "center": [50, 55], "diameter": 6 }],
                         "cutouts": [{ "outline": [[70,58],[80,58],[80,62],[70,62]] }] },
               "flanges": [{ "id": "f", "on": { "base": 2 }, "angle": 90, "length": 8,
                             "datum": "outer", "position": "outside", "side": "up",
                             "bend_radius": 1 }]"#,
        );
        let findings = check(&spec, &Rules::default()).unwrap();
        let got: Vec<(Rule, &str, f64)> = findings
            .iter()
            .map(|f| (f.rule, f.feature.as_str(), (f.value * 1e6).round() / 1e6))
            .collect();
        assert_eq!(
            got,
            vec![
                (Rule::BendRadius, "f", 1.0),
                (Rule::FlangeLength, "f", 5.0),
                (Rule::HoleSize, "hole 1", 1.5),
                (Rule::HoleToEdge, "hole 1", 1.0),
                (Rule::HoleToBend, "hole 2", 1.0),
                (Rule::HoleToBend, "cut-out 1", 0.0),
            ]
        );
        assert_eq!(
            findings[3].to_string(),
            "hole 1: distance to the outline 1.00 mm, at least 2.00 mm"
        );
    }
}
