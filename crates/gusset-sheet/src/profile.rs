//! A part made from a profile (a contour flange as the first
//! feature) as a base plate with a contour flange on each end: the same
//! sheet, bent the same way.

use std::borrow::Cow;

use crate::hem::contour;
use crate::spec::{BaseFace, Contour, EdgeRef, Hand, PartSpec, Side};
use crate::SheetError;

impl PartSpec {
    /// This spec with each jog as its two folds (after the part's own).
    fn without_jogs(&self) -> Result<PartSpec, SheetError> {
        let mut spec = self.clone();
        let (t, r) = (self.thickness, self.bend_radius);
        for j in std::mem::take(&mut spec.jogs) {
            let angle = j.angle.unwrap_or(90.0);
            let theta = angle.to_radians();
            let run = (j.offset - (2.0 * r + t) * (1.0 - theta.cos())) / theta.sin();
            if !(angle > 0.0 && angle < 180.0 && run >= 0.0) {
                return Err(SheetError::Invalid(format!(
                    "jog {}: an offset of {} is too small for two bends of {angle} degrees",
                    j.id, j.offset
                )));
            }
            let (allowance, _) = self.bend_allowance(angle, r, None)?;
            let [a, b] = j.line;
            let d = [b[0] - a[0], b[1] - a[1]];
            let len = (d[0] * d[0] + d[1] * d[1]).sqrt();
            // The turning side: right of the line.
            let n = [d[1] / len, -d[0] / len];
            let shift = allowance + run;
            let fold = |id: String, line: [[f64; 2]; 2], side: Side| crate::spec::FoldLine {
                id,
                line,
                angle,
                location: crate::spec::BendLocation::Start,
                side,
                bend_radius: None,
                k_factor: None,
                bounded: false,
            };
            let back = if j.side == Side::Up {
                Side::Down
            } else {
                Side::Up
            };
            spec.folds.push(fold(format!("{}.1", j.id), j.line, j.side));
            spec.folds.push(fold(
                format!("{}.2", j.id),
                [
                    [a[0] + shift * n[0], a[1] + shift * n[1]],
                    [b[0] + shift * n[0], b[1] + shift * n[1]],
                ],
                back,
            ));
        }
        Ok(spec)
    }

    /// This spec, with a profile base replaced by its base plate and contours.
    /// A contour along several edges as one contour per edge, `id` on the
    /// first and `id@e` on edge e: one feature, mitered where they meet.
    fn without_along(&self) -> Result<PartSpec, SheetError> {
        let n = self.base.outline.len();
        let mut part = self.clone();
        let mut contours = Vec::new();
        for c in &self.contours {
            if c.id.contains('@') {
                return Err(SheetError::Invalid(format!(
                    "contour {}: '@' is reserved in ids",
                    c.id
                )));
            }
            if c.along.is_empty() {
                contours.push(c.clone());
                continue;
            }
            let unsupported =
                |what: &str| Err(SheetError::Unsupported(format!("contour {}: {what}", c.id)));
            let EdgeRef::Base(first) = c.on else {
                return unsupported("along several edges from a flange's tip");
            };
            if c.width.is_some() {
                return unsupported("a partial width along several edges");
            }
            let mut previous = first;
            for &e in &c.along {
                if n == 0 || e != (previous + 1) % n || e == first {
                    return Err(SheetError::Invalid(format!(
                        "contour {}: `along` must name the following base edges in turn",
                        c.id
                    )));
                }
                previous = e;
            }
            let mut head = c.clone();
            head.along = Vec::new();
            contours.push(head.clone());
            for &e in &c.along {
                contours.push(Contour {
                    id: format!("{}@{e}", c.id),
                    on: EdgeRef::Base(e),
                    ..head.clone()
                });
            }
        }
        part.contours = contours;
        Ok(part)
    }

    pub fn resolved(&self) -> Result<Cow<'_, PartSpec>, SheetError> {
        if !self.faces.is_empty() || !self.bends.is_empty() || self.has_face_edges() {
            return Ok(Cow::Owned(
                crate::faces::without_faces(self)?.resolved()?.into_owned(),
            ));
        }
        if !self.jogs.is_empty() {
            return Ok(Cow::Owned(self.without_jogs()?.resolved()?.into_owned()));
        }
        if self.contours.iter().any(|c| !c.along.is_empty()) {
            return Ok(Cow::Owned(self.without_along()?.resolved()?.into_owned()));
        }
        let Some(p) = &self.base.profile else {
            return Ok(Cow::Borrowed(self));
        };
        let invalid = |what: &str| Err(SheetError::Invalid(format!("profile: {what}")));
        if !self.base.outline.is_empty() {
            return invalid("a base has an outline or a profile, not both");
        }
        if !(p.width > 0.0) {
            return invalid("needs a width");
        }
        let (pts, t) = (&p.points, self.thickness);
        let bulge = |i: usize| p.bulges.get(i).copied().unwrap_or(0.0);
        // The base: the first straight line on z = 0 with the material above
        // it, or on z = T with the material below it.
        let base = (0..pts.len().saturating_sub(1)).find(|&i| {
            let (a, b) = (pts[i], pts[i + 1]);
            let above = (b[0] > a[0]) == (p.material == Hand::Left);
            bulge(i) == 0.0
                && (b[0] - a[0]).abs() > 1e-9
                && (b[1] - a[1]).abs() < 1e-9
                && ((a[1].abs() < 1e-9 && above) || ((a[1] - t).abs() < 1e-9 && !above))
        });
        let Some(k) = base else {
            return invalid(
                "needs a straight line on z = 0 with the material above it, \
                 or on z = T with the material below it (the base)",
            );
        };
        let (lo, hi) = (pts[k][0].min(pts[k + 1][0]), pts[k][0].max(pts[k + 1][0]));
        let face = if pts[k][1].abs() < 1e-9 {
            BaseFace::Bottom
        } else {
            BaseFace::Top
        };
        // A chain of the profile walked from the base: on edge 2 (at y = hi,
        // outward +y) or edge 0 (at y = lo, outward -y: a mirror, which turns
        // the arcs the other way).
        let chain = |id: String, points: Vec<[f64; 2]>, bulges: Vec<f64>| {
            let high = points[0][0] == hi;
            Contour {
                id,
                on: EdgeRef::Base(if high { 2 } else { 0 }),
                width: None,
                profile: points
                    .iter()
                    .map(|q| [if high { q[0] - hi } else { lo - q[0] }, q[1]])
                    .collect(),
                bulges: bulges.iter().map(|b| if high { *b } else { -b }).collect(),
                face,
                bend_radius: None,
                k_factor: None,
                along: Vec::new(),
            }
        };
        let mut contours = Vec::new();
        if k > 0 {
            // Walked backwards, which turns the arcs the other way too.
            contours.push(chain(
                "profile.start".into(),
                pts[..=k].iter().rev().copied().collect(),
                (0..k).rev().map(|i| -bulge(i)).collect(),
            ));
        }
        if k + 2 < pts.len() {
            contours.push(chain(
                "profile.end".into(),
                pts[k + 1..].to_vec(),
                (k + 1..pts.len() - 1).map(bulge).collect(),
            ));
        }
        // Flanges and hems sit on the profile's free ends, {"tip": "start"} and
        // {"tip": "end"}: the last bend of that end's contour.
        let mut ends = Vec::new();
        for c in &contours {
            let last = contour(self, c)?.pop().map(|s| s.id);
            ends.push((c.id.trim_start_matches("profile."), last));
        }
        let end = |on: &EdgeRef, _id: &str| match on {
            EdgeRef::Tip(t) if t == "start" || t == "end" => {
                match ends.iter().find(|e| e.0 == t).and_then(|e| e.1.clone()) {
                    Some(step) => Ok::<_, SheetError>(EdgeRef::Tip(step)),
                    // An end where the base line is: that base edge (up is +z).
                    None => {
                        let y = if t == "start" {
                            pts[k][0]
                        } else {
                            pts[k + 1][0]
                        };
                        Ok(EdgeRef::Base(if y == hi { 2 } else { 0 }))
                    }
                }
            }
            // The base plate's edges: 0 at y = lo, 1 at x = width, 2 at
            // y = hi, 3 at x = 0 (0 and 2 only where the profile ends).
            on => Ok(on.clone()),
        };
        let mut spec = self.clone();
        for f in &mut spec.flanges {
            f.on = end(&f.on, &f.id)?;
        }
        for h in &mut spec.hems {
            h.on = end(&h.on, &h.id)?;
        }
        for c in &mut spec.contours {
            c.on = end(&c.on, &c.id)?;
        }
        spec.base.profile = None;
        spec.base.outline = vec![[0.0, lo], [p.width, lo], [p.width, hi], [0.0, hi]];
        // The profile's own contours first, the part's after them.
        contours.append(&mut spec.contours);
        spec.contours = contours;
        Ok(Cow::Owned(spec))
    }
}

#[cfg(test)]
mod tests {
    use crate::fold::{exact_mass, fold};
    use crate::spec::{EdgeRef, PartSpec};

    fn part(profile: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "profile": {profile} }} }}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_profile_folds_with_their_volume_and_faces() {
        // The base is the first line on z = 0 (material above) or z = T
        // (material below); the rest are contours on its ends.
        // Expected: volume, faces.
        for (name, profile, volume, faces) in [
            (
                "L, outer face",
                r#"{ "points": [[40,0],[0,0],[0,20]], "material": "right", "width": 100 }"#,
                11342.4778,
                14,
            ),
            (
                "L, inner face",
                r#"{ "points": [[0,22],[0,2],[40,2]], "material": "right", "width": 100 }"#,
                12142.4778,
                14,
            ),
            (
                "channel",
                r#"{ "points": [[0,40],[0,0],[60,0],[60,40]], "material": "left", "width": 100 }"#,
                26684.9556,
                22,
            ),
            (
                "hat",
                r#"{ "points": [[0,0],[20,0],[20,30],[60,30],[60,0],[80,0]],
                     "material": "left", "width": 100 }"#,
                26969.9112,
                38,
            ),
            (
                "Z, base on top",
                r#"{ "points": [[0,30],[0,2],[40,2],[40,-26]], "material": "right", "width": 100 }"#,
                18684.9556,
                22,
            ),
            (
                "slanted",
                r#"{ "points": [[-10,30],[0,0],[50,0],[60,30]], "material": "left", "width": 100 }"#,
                21841.5361,
                22,
            ),
            (
                "rounded channel",
                &format!(
                    r#"{{ "points": [[0,40],[0,10],[10,0],[50,0],[60,10],[60,40]],
                         "bulges": [0,{q},0,{q},0], "material": "left", "width": 100 }}"#,
                    q = (std::f64::consts::PI / 8.0).tan()
                ),
                25654.8668,
                22,
            ),
        ] {
            let spec = part(&profile);
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn hems_and_flanges_sit_on_a_profiles_free_ends() {
        let channel = r#""base": { "profile": { "points": [[0,40],[0,0],[60,0],[60,40]],
                                                 "material": "left", "width": 100 } }"#;
        let part = |features: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     {channel}, {features} }}"#
            ))
            .unwrap()
        };
        // Up turns as the profile's last bend: the end's goes inward.
        let hems = part(
            r#""hems": [
                { "id": "h1", "on": { "tip": "end" }, "kind": { "single": { "gap": 1, "length": 8 } }, "side": "down" },
                { "id": "h2", "on": { "tip": "start" }, "kind": { "single": { "gap": 1, "length": 8 } }, "side": "up" }]"#,
        );
        let lips = part(
            r#""flanges": [
                { "id": "l1", "on": { "tip": "end" }, "angle": 90, "length": 15,
                  "datum": "outer", "position": "outside", "side": "up" },
                { "id": "l2", "on": { "tip": "start" }, "angle": 90, "length": 15,
                  "datum": "outer", "position": "outside", "side": "down" }]"#,
        );
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("hems", hems, 29769.9112, 38),
            ("lips", lips, 32169.9112, 38),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn flanges_sit_on_a_profile_parts_base_edges() {
        let part = |profile: &str, flanges: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "profile": {profile} }}, "flanges": [{flanges}] }}"#
            ))
            .unwrap()
        };
        let flange = |id: &str, edge: usize, length: f64, side: &str| {
            format!(
                r#"{{ "id": "{id}", "on": {{ "base": {edge} }}, "angle": 90, "length": {length},
                     "datum": "outer", "position": "outside", "side": "{side}" }}"#
            )
        };
        let channel =
            r#"{ "points": [[0,40],[0,0],[60,0],[60,40]], "material": "left", "width": 100 }"#;
        let l = r#"{ "points": [[40,0],[0,0],[0,20]], "material": "right", "width": 100 }"#;
        // Expected: volume, faces. The profile is built first, so the end
        // flanges stop at the legs' bend starts.
        for (name, spec, volume, faces) in [
            (
                "channel ends",
                part(
                    channel,
                    &format!(
                        "{},{}",
                        flange("e1", 1, 25.0, "up"),
                        flange("e3", 3, 25.0, "up")
                    ),
                ),
                31541.7343,
                42,
            ),
            (
                "L's free base edge",
                part(l, &flange("f", 2, 15.0, "down")),
                14084.9556,
                22,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn contours_and_folds_sit_on_a_profile_part() {
        let part = |profile: &str, features: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "profile": {profile} }}, {features} }}"#
            ))
            .unwrap()
        };
        let channel =
            r#"{ "points": [[0,40],[0,0],[60,0],[60,40]], "material": "left", "width": 100 }"#;
        let l = r#"{ "points": [[40,0],[0,0],[0,20]], "material": "right", "width": 100 }"#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "a contour on the channel's end",
                part(
                    channel,
                    r#""contours": [{ "id": "c", "on": { "base": 1 },
                                      "profile": [[0,2],[0,22],[10,22]], "face": "bottom" }]"#,
                ),
                29669.7343,
                40,
            ),
            (
                "a fold across the L's base",
                part(
                    l,
                    r#""folds": [{ "id": "fold", "line": [[100,25],[0,25]], "angle": 90,
                                   "location": "centerline", "side": "up" }]"#,
                ),
                11405.3096,
                22,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn an_end_without_a_bend_is_the_base_edge_there() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "profile": { "points": [[40,0],[0,0],[0,20]], "material": "right",
                                        "width": 100 } },
                 "hems": [{ "id": "h", "on": { "tip": "start" },
                            "kind": { "single": { "gap": 1, "length": 8 } }, "side": "up" }] }"#,
        )
        .unwrap();
        assert_eq!(spec.resolved().unwrap().hems[0].on, EdgeRef::Base(2));
        // Expected: volume, faces.
        let m = exact_mass(&spec).unwrap();
        assert!((m.volume - 12884.9556).abs() < 1e-3, "{}", m.volume);
        assert_eq!(fold(&spec).unwrap().topology.faces.len(), 22);
    }

    #[test]
    fn a_profile_without_a_base_line_is_refused() {
        // The material below z = 0: no base plate at z = 0..T.
        let spec =
            part(r#"{ "points": [[0,20],[0,0],[40,0]], "material": "right", "width": 100 }"#);
        assert!(matches!(
            exact_mass(&spec),
            Err(crate::SheetError::Invalid(_))
        ));
    }
}
