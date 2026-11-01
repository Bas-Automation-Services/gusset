//! Flanges, contour flanges and hems as one list of bend steps: a bend of some
//! angle and inner radius, starting `setback` from where it sits, followed by a
//! straight piece. Contours and hems are short chains of steps; the rules are
//! these.

use crate::section::section;
use crate::spec::{BaseFace, Contour, EdgeRef, HemKind, PartSpec, Side, Width};
use crate::SheetError;

#[derive(Clone, Debug)]
pub(crate) struct Step {
    pub id: String,
    pub on: EdgeRef,
    pub width: Option<Width>,
    pub side: Side,
    /// Bend angle in degrees.
    pub angle: f64,
    pub radius: f64,
    pub k_factor: f64,
    pub setback: f64,
    pub bend_allowance: f64,
    pub straight: f64,
}

impl Step {
    pub fn reach(&self) -> f64 {
        self.setback + self.bend_allowance + self.straight
    }
}

pub(crate) fn steps(spec: &PartSpec) -> Result<Vec<Step>, SheetError> {
    let t = spec.thickness;
    let mut steps = Vec::new();
    // A profile part's own contours are its first feature.
    for c in spec
        .contours
        .iter()
        .filter(|c| c.id.starts_with("profile."))
    {
        steps.extend(contour(spec, c)?);
    }
    for f in &spec.flanges {
        let r = f.bend_radius.unwrap_or(spec.bend_radius);
        let (bend_allowance, k) = spec.bend_allowance(f.angle, r, f.k_factor)?;
        let s = section(t, r, k, f)?;
        steps.push(Step {
            id: f.id.clone(),
            on: f.on.clone(),
            width: f.width,
            side: f.side,
            angle: f.angle,
            radius: r,
            k_factor: k,
            setback: s.setback,
            bend_allowance,
            straight: s.straight,
        });
    }
    for c in spec
        .contours
        .iter()
        .filter(|c| !c.id.starts_with("profile."))
    {
        steps.extend(contour(spec, c)?);
    }
    for h in &spec.hems {
        // Bend allowance and K-factor of a fold of the hem.
        let fold = |angle: f64, r: f64| spec.bend_allowance(angle, r, h.k_factor);
        let invalid = |what: &str| Err(SheetError::Invalid(format!("hem {}: {what}", h.id)));
        // Every hem starts with its outer bend face flush with the edge.
        let first = |angle: f64, r: f64, straight: f64| -> Result<Step, SheetError> {
            let (bend_allowance, k_factor) = fold(angle, r)?;
            Ok(Step {
                id: h.id.clone(),
                on: h.on.clone(),
                width: h.width,
                side: h.side,
                angle,
                radius: r,
                k_factor,
                setback: -(r + t),
                bend_allowance,
                straight,
            })
        };
        // Further folds continue straight on, turning the same way.
        let next = |n: usize,
                    parent: String,
                    angle: f64,
                    r: f64,
                    straight: f64|
         -> Result<Step, SheetError> {
            let (bend_allowance, k_factor) = fold(angle, r)?;
            Ok(Step {
                id: format!("{}.{n}", h.id),
                on: EdgeRef::Tip(parent),
                width: None,
                side: Side::Up,
                angle,
                radius: r,
                k_factor,
                setback: 0.0,
                bend_allowance,
                straight,
            })
        };
        match h.kind {
            HemKind::Single { gap, length } => {
                let r = gap / 2.0;
                let straight = length - (r + t);
                if !(gap > 0.0 && straight > 0.0) {
                    return invalid("needs a gap and a length beyond the fold");
                }
                steps.push(first(180.0, r, straight)?);
            }
            HemKind::Double { gap, length } => {
                // Up 90, across T + gap, over 90, back along the top, down 180, back.
                let r = gap / 2.0;
                let across = length - 2.0 * (r + t);
                if !(gap > 0.0 && across > 0.0) {
                    return invalid("needs a gap and a length beyond both folds");
                }
                steps.push(first(90.0, r, t + gap)?);
                steps.push(next(2, h.id.clone(), 90.0, r, across)?);
                steps.push(next(3, format!("{}.2", h.id), 180.0, r, across)?);
            }
            HemKind::Teardrop { radius: r, angle } => {
                if !(r > 0.0 && angle > 180.0 && angle < 360.0) {
                    return invalid("a teardrop needs a radius and an angle between 180 and 360");
                }
                // Section with the outer plane at z = 0: the bend's centre is at
                // z = T + R. The tail runs on until its lower end corner is T/4
                // above the sheet (z = 1.25 T).
                let (sin, cos) = angle.to_radians().sin_cos();
                let inner = t + r - r * cos;
                let outer = t + r - (r + t) * cos;
                let tail = if sin < -1e-12 {
                    ((inner.min(outer) - 1.25 * t) / -sin).max(0.0)
                } else {
                    0.0
                };
                steps.push(first(angle, r, tail)?);
            }
            HemKind::Rolled { radius: r, angle } => {
                if !(r > 0.0 && angle > 0.0 && angle < 360.0) {
                    return invalid("a rolled hem needs a radius and an angle below 360");
                }
                steps.push(first(angle, r, 0.0)?);
            }
        }
    }
    // A tip reference to a contour means its free end: its last step.
    for c in &spec.contours {
        let own = |id: &str| id == c.id || id.starts_with(&format!("{}.", c.id));
        let Some(last) = steps
            .iter()
            .rev()
            .find(|s| own(&s.id))
            .map(|s| s.id.clone())
        else {
            continue;
        };
        for s in steps.iter_mut().filter(|s| !own(&s.id)) {
            if s.on == EdgeRef::Tip(c.id.clone()) {
                s.on = EdgeRef::Tip(last.clone());
            }
        }
    }
    // A step on a tip goes after its parent (a flange may sit on a contour's
    // tip); otherwise the order stays.
    let mut ordered: Vec<Step> = Vec::new();
    while !steps.is_empty() {
        let ready = steps
            .iter()
            .position(|s| match &s.on {
                EdgeRef::Tip(p) => {
                    ordered.iter().any(|o| &o.id == p) || !steps.iter().any(|q| &q.id == p)
                }
                EdgeRef::Base(_) => true,
            })
            .unwrap_or(0);
        ordered.push(steps.remove(ready));
    }
    Ok(ordered)
}

/// A contour flange as a chain of steps: the profile is one face of the sheet and the base
/// continues it. A sharp corner is a bend of the part's radius whose inner
/// (material outside the turn) or outer (inside) face has its corner there; a
/// tangent arc is a bend whose face on the profile is that arc.
pub(crate) fn contour(spec: &PartSpec, c: &Contour) -> Result<Vec<Step>, SheetError> {
    let t = spec.thickness;
    let r = c.bend_radius.unwrap_or(spec.bend_radius);
    let fail = |what: &str| SheetError::Invalid(format!("contour {}: {what}", c.id));
    let unsupported = |what: &str| SheetError::Unsupported(format!("contour {}: {what}", c.id));
    if c.width.is_some() && !matches!(c.on, EdgeRef::Base(_)) {
        return Err(fail("a contour on a tip spans its parent"));
    }
    let p = &c.profile;
    if p.len() < 2 {
        return Err(fail("a profile needs two points"));
    }
    let bulge = |i: usize| c.bulges.get(i).copied().unwrap_or(0.0);
    // The first line meets the base face the profile continues at the first
    // corner; the base runs toward +s up to it. A first arc starts on that
    // face, tangent to it.
    let (z_face, material_left) = match c.face {
        BaseFace::Bottom => (0.0, true),
        BaseFace::Top => (t, false),
    };
    let pts = if bulge(0) != 0.0 {
        if (p[0][1] - z_face).abs() > 1e-9 {
            return Err(fail(
                "a profile starting with an arc must start on the base face",
            ));
        }
        p.clone()
    } else {
        let d = [p[1][0] - p[0][0], p[1][1] - p[0][1]];
        if d[1].abs() < 1e-9 {
            return Err(unsupported("a first line along the base"));
        }
        let w = (z_face - p[0][1]) / d[1];
        let mut pts = vec![[p[0][0] + w * d[0], z_face]];
        pts.extend_from_slice(&p[1..]);
        pts
    };

    // The section as lines and bends. A bend: signed turn (left, seen with s
    // right and z up, positive), inner radius, and how far its tangent points
    // sit from the corner (a sharp corner trims its neighbouring lines).
    enum Item {
        Line(f64),
        Bend { turn: f64, radius: f64, trim: f64 },
    }
    let rotate = |v: [f64; 2], a: f64| {
        let (sin, cos) = a.sin_cos();
        [v[0] * cos - v[1] * sin, v[0] * sin + v[1] * cos]
    };
    let angle_between =
        |a: [f64; 2], b: [f64; 2]| (a[0] * b[1] - a[1] * b[0]).atan2(a[0] * b[0] + a[1] * b[1]);
    // The profile is a bend's outer face where the material is inside the turn.
    let outer = |turn: f64| (turn > 0.0) == material_left;
    let mut items = Vec::new();
    let mut previous = ([1.0, 0.0], false);
    for (i, q) in pts.windows(2).enumerate() {
        let chord = [q[1][0] - q[0][0], q[1][1] - q[0][1]];
        let length = (chord[0] * chord[0] + chord[1] * chord[1]).sqrt();
        let dir = [chord[0] / length, chord[1] / length];
        let sweep = 4.0 * bulge(i).atan();
        let (start, end) = (rotate(dir, -sweep / 2.0), rotate(dir, sweep / 2.0));
        let arc = sweep != 0.0;
        let (before, before_arc) = previous;
        let turn = angle_between(before, start);
        if before_arc || arc {
            if turn.abs() > 1e-9 {
                return Err(unsupported("an arc that isn't tangent to its neighbours"));
            }
        } else {
            if !(turn.abs() > 1e-9 && turn.abs() < std::f64::consts::PI - 1e-9) {
                return Err(fail(
                    "every corner must turn by more than 0 and less than 180 degrees",
                ));
            }
            let radius = if outer(turn) { r + t } else { r };
            items.push(Item::Bend {
                turn,
                radius: r,
                trim: radius * (turn.abs() / 2.0).tan(),
            });
        }
        if arc {
            let on_profile = length / (2.0 * (sweep.abs() / 2.0).sin());
            let radius = if outer(sweep) {
                on_profile - t
            } else {
                on_profile
            };
            if radius <= 0.0 {
                return Err(fail("an arc smaller than the thickness on its inside"));
            }
            items.push(Item::Bend {
                turn: sweep,
                radius,
                trim: 0.0,
            });
        } else {
            items.push(Item::Line(length));
        }
        previous = (end, arc);
    }

    // One step per bend: the straight piece runs to the next bend.
    let bends: Vec<usize> = (0..items.len())
        .filter(|&i| matches!(items[i], Item::Bend { .. }))
        .collect();
    let mut steps = Vec::new();
    for (n, &i) in bends.iter().enumerate() {
        let Item::Bend { turn, radius, trim } = items[i] else {
            unreachable!()
        };
        let next = bends.get(n + 1).copied().unwrap_or(items.len());
        let next_trim = match items.get(next) {
            Some(&Item::Bend { trim, .. }) => trim,
            _ => 0.0,
        };
        let lines: f64 = items[i + 1..next]
            .iter()
            .map(|item| match item {
                Item::Line(l) => *l,
                _ => 0.0,
            })
            .sum();
        let straight = lines - trim - next_trim;
        if straight < -1e-9 {
            return Err(fail("a profile line is too short for its bends"));
        }
        let (bend_allowance, k) =
            spec.bend_allowance(turn.abs().to_degrees(), radius, c.k_factor)?;
        let id = |n: usize| match n {
            0 => c.id.clone(),
            _ => format!("{}.{}", c.id, n + 1),
        };
        let parent_turn = n.checked_sub(1).map(|m| match items[bends[m]] {
            Item::Bend { turn, .. } => turn,
            _ => unreachable!(),
        });
        steps.push(Step {
            id: id(n),
            on: if n == 0 {
                c.on.clone()
            } else {
                EdgeRef::Tip(id(n - 1))
            },
            width: if n == 0 { c.width } else { None },
            // On the edge: up for a left turn. On a tip: up turns as the parent.
            side: if parent_turn.map_or(turn > 0.0, |pt| (pt > 0.0) == (turn > 0.0)) {
                Side::Up
            } else {
                Side::Down
            },
            angle: turn.abs().to_degrees(),
            radius,
            k_factor: k,
            setback: if n == 0 { pts[0][0] - trim } else { 0.0 },
            bend_allowance,
            straight: straight.max(0.0),
        });
    }
    Ok(steps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::Hem;

    /// Flat width of a 60 mm deep plate (T 2 unless given, K 0.4) with this hem
    /// on its far edge.
    fn width(kind: HemKind, t: f64) -> f64 {
        let spec: PartSpec = serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": {t}, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }} }}"#
        ))
        .unwrap();
        let spec = PartSpec {
            hems: vec![Hem {
                id: "h".into(),
                on: EdgeRef::Base(2),
                width: None,
                kind,
                side: Side::Up,
                k_factor: None,
            }],
            ..spec
        };
        // A chain: every step starts where the previous one ends.
        60.0 + steps(&spec).unwrap().iter().map(|s| s.reach()).sum::<f64>()
    }

    #[test]
    fn every_hem_kind_gives_its_flat_width() {
        use HemKind::*;
        for (kind, t, expected) in [
            (
                Single {
                    gap: 1.0,
                    length: 10.0,
                },
                2.0,
                69.084070,
            ),
            (
                Single {
                    gap: 2.0,
                    length: 15.0,
                },
                1.0,
                75.398230,
            ),
            (
                Double {
                    gap: 1.0,
                    length: 10.0,
                },
                2.0,
                78.668141,
            ),
            (
                Double {
                    gap: 2.0,
                    length: 15.0,
                },
                2.0,
                90.309734,
            ),
            (
                Double {
                    gap: 0.5,
                    length: 8.0,
                },
                1.0,
                75.334070,
            ),
            (
                Teardrop {
                    radius: 3.0,
                    angle: 270.0,
                },
                2.0,
                75.407078,
            ),
            (
                Teardrop {
                    radius: 3.0,
                    angle: 240.0,
                },
                2.0,
                75.536205,
            ),
            (
                Teardrop {
                    radius: 4.0,
                    angle: 270.0,
                },
                2.0,
                80.119467,
            ),
            (
                Teardrop {
                    radius: 3.0,
                    angle: 270.0,
                },
                1.0,
                74.772123,
            ),
            (
                Teardrop {
                    radius: 3.0,
                    angle: 300.0,
                },
                2.0,
                74.896753,
            ),
            (
                Rolled {
                    radius: 5.0,
                    angle: 240.0,
                },
                2.0,
                77.294983,
            ),
            (
                Rolled {
                    radius: 4.0,
                    angle: 200.0,
                },
                1.0,
                70.358897,
            ),
        ] {
            let w = width(kind, t);
            assert!(
                (w - expected).abs() < 1e-5,
                "{kind:?} T {t}: {w} vs {expected}"
            );
        }
    }
}
