//! Formed embosses (`base.embosses`): the sheet pressed out of the base plate
//! at its thickness. Seen from the side (the meridian half-plane, r out from
//! the emboss's axis, z up) each face of the sheet runs from the plate over a
//! rounded corner, along the wall, over a second rounded corner into the
//! raised face. Turned around the axis that gives tori, a cone (or a
//! cylinder) and a disk, all exact; the plate keeps a round hole, the
//! footprint, where the sheet leaves it.

use crate::flat::{Arc, Curve, Ring, P2};
use crate::fold::{Mass, GAUSS, P3};
use crate::spec::{Emboss, EmbossShape, Side};
use crate::SheetError;
use opengeometry::brep::{
    boundary, plane_boundary, uv_line, Builder, CurveGeometry, Frame3, GeometryError, Orientation,
    Surface, SurfaceGeometry,
};
use opengeometry::math::Interval;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const EPS: f64 = 1e-9;
const X: P3 = [1.0, 0.0, 0.0];
const Y: P3 = [0.0, 1.0, 0.0];
const Z: P3 = [0.0, 0.0, 1.0];

/// A piece of a face's profile in the meridian half-plane, (r, z).
#[derive(Clone, Copy, Debug)]
enum Meridian {
    /// Around `c` with radius `rho`, from angle `a0` to `a1` (from +r towards
    /// +z).
    Arc {
        c: P2,
        rho: f64,
        a0: f64,
        a1: f64,
    },
    Line {
        a: P2,
        b: P2,
    },
}

impl Meridian {
    /// The point a fraction s along, and its derivative.
    fn at(&self, s: f64) -> (P2, P2) {
        match *self {
            Meridian::Arc { c, rho, a0, a1 } => {
                let a = a0 + s * (a1 - a0);
                let (sn, cs) = a.sin_cos();
                (
                    [c[0] + rho * cs, c[1] + rho * sn],
                    [-rho * sn * (a1 - a0), rho * cs * (a1 - a0)],
                )
            }
            Meridian::Line { a, b } => (
                [a[0] + s * (b[0] - a[0]), a[1] + s * (b[1] - a[1])],
                [b[0] - a[0], b[1] - a[1]],
            ),
        }
    }

    /// Mirrored about z = t / 2.
    fn mirrored(self, t: f64) -> Meridian {
        match self {
            Meridian::Arc { c, rho, a0, a1 } => Meridian::Arc {
                c: [c[0], t - c[1]],
                rho,
                a0: -a0,
                a1: -a1,
            },
            Meridian::Line { a, b } => Meridian::Line {
                a: [a[0], t - a[1]],
                b: [b[0], t - b[1]],
            },
        }
    }
}

/// An emboss, worked out.
#[derive(Clone, Debug)]
pub(crate) struct Formed {
    pub id: String,
    pub center: P2,
    /// Where the sheet leaves the plate.
    pub footprint: f64,
    /// Per face of the sheet, by the plate face it leaves (0: z = 0, 1:
    /// z = T): its profile from the footprint inward, the last piece the
    /// raised face out to the axis.
    chains: [[Meridian; 4]; 2],
}

impl Formed {
    pub fn new(e: &Emboss, t: f64, bend_radius: f64) -> Result<Formed, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("emboss {}: {what}", e.id)));
        let EmbossShape::Round { diameter } = e.shape;
        let r = e.radius.unwrap_or(bend_radius);
        let (h, d) = (e.height, e.draft.to_radians());
        if !(diameter > 0.0 && h > 0.0 && r > 0.0) {
            return invalid("needs a positive diameter, height and radius");
        }
        if !(d >= 0.0 && d < FRAC_PI_2) {
            return invalid("the draft must lie between 0 and 90 degrees");
        }
        // Tangent length of each corner: they turn by 90 degrees less the
        // draft.
        let k = ((FRAC_PI_2 - d) / 2.0).tan();
        let half = diameter / 2.0;
        // Plate, wall and raised face of the face on the raised side; its
        // first corner has the inside radius r, the second r + t.
        let foot = half + h * d.tan() + r * k;
        let top = half - (r + t) * k;
        let c1 = [foot, t + r];
        let c2 = [top, t + h - (r + t)];
        let nw = [d.cos(), d.sin()];
        let off = |c: P2, s: f64| [c[0] + s * nw[0], c[1] + s * nw[1]];
        if top <= r + t + EPS {
            return invalid("too small across for its radius and thickness");
        }
        let (a, b) = (off(c1, -r), off(c2, r + t));
        if (b[1] - a[1]) * d.cos() - (b[0] - a[0]) * d.sin() <= EPS {
            return invalid("too low for its radius and thickness");
        }
        let corner = |c: P2, rho: f64, a0: f64, a1: f64| Meridian::Arc { c, rho, a0, a1 };
        let up = [
            [
                corner(c1, r + t, 1.5 * PI, PI + d),
                Meridian::Line {
                    a: off(c1, -(r + t)),
                    b: off(c2, r),
                },
                corner(c2, r, d, FRAC_PI_2),
                Meridian::Line {
                    a: [top, h],
                    b: [0.0, h],
                },
            ],
            [
                corner(c1, r, 1.5 * PI, PI + d),
                Meridian::Line { a, b },
                corner(c2, r + t, d, FRAC_PI_2),
                Meridian::Line {
                    a: [top, t + h],
                    b: [0.0, t + h],
                },
            ],
        ];
        let chains = match e.side {
            Side::Up => up,
            // Mirrored, the face that left z = T leaves z = 0.
            Side::Down => [up[1].map(|m| m.mirrored(t)), up[0].map(|m| m.mirrored(t))],
        };
        Ok(Formed {
            id: e.id.clone(),
            center: e.at,
            footprint: foot,
            chains,
        })
    }

    /// The footprint as a ring, for the flat pattern.
    pub fn outline(&self) -> Ring {
        let (c, r) = (self.center, self.footprint);
        let arc = Some(Curve::Arc(Arc {
            center: c,
            ccw: true,
        }));
        Ring {
            points: vec![[c[0] + r, c[1]], [c[0] - r, c[1]]],
            arcs: vec![arc, arc],
        }
    }

    /// Where the meridian point q lies at angle u around the axis (the base
    /// frame).
    fn point(&self, q: P2, u: f64) -> P3 {
        [
            self.center[0] + q[0] * u.cos(),
            self.center[1] + q[0] * u.sin(),
            q[1],
        ]
    }

    /// The sheet's material between the two faces, piece by piece, turned
    /// around the axis; and both faces' areas. Base frame.
    pub fn mass(&self) -> Mass {
        let (mut r1, mut rz, mut area) = (0.0, 0.0, 0.0);
        for i in 0..4 {
            let [p, q] = [self.chains[0][i], self.chains[1][i]];
            // The region between the two pieces: along p, across to q's
            // end, back along q, across to p's start. Green's theorem with
            // F = r^2 / 2 and r^2 z / 2 gives the integrals of r and r z.
            let across = |a: P2, b: P2| Meridian::Line { a, b };
            let (mut i1, mut i2) = (0.0, 0.0);
            for (m, sign) in [
                (p, 1.0),
                (across(p.at(1.0).0, q.at(1.0).0), 1.0),
                (q, -1.0),
                (across(p.at(0.0).0, q.at(0.0).0), -1.0),
            ] {
                let (a, b) = green(m);
                i1 += sign * a;
                i2 += sign * b;
            }
            // The loop runs either way round.
            let s = i1.signum();
            r1 += s * i1;
            rz += s * i2;
            for m in [p, q] {
                area += TAU * length_moment(m);
            }
        }
        let volume = TAU * r1;
        Mass {
            volume,
            area,
            center_of_mass: [self.center[0], self.center[1], rz / r1],
        }
    }

    /// Builds both faces of the sheet in `b`, on the footprint's circles
    /// (`circles[level]`, each closed at vertex `v[level]`, starting along +x).
    pub fn build(
        &self,
        b: &mut Builder,
        circles: [u32; 2],
        v: [u32; 2],
    ) -> Result<(), GeometryError> {
        for level in 0..2 {
            let chain = &self.chains[level];
            let other = &self.chains[1 - level];
            // Circles and their vertices where the pieces meet, from the
            // footprint inward.
            let mut rings = vec![(circles[level], v[level])];
            for m in &chain[..3] {
                let q = m.at(1.0).0;
                let vertex = b.vertex(self.point(q, 0.0));
                let circle = b.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: [self.center[0], self.center[1], q[1]],
                            x: X,
                            y: Y,
                            z: Z,
                        },
                        radius: q[0],
                    },
                    Interval::new(0.0, TAU)?,
                    false,
                );
                rings.push((circle, vertex));
            }
            for (i, m) in chain.iter().enumerate() {
                // Outward: away from the other face.
                let (here, there) = (m.at(0.5).0, other[i].at(0.5).0);
                let out = [here[0] - there[0], here[1] - there[1]];
                let (circle_a, va) = rings[i];
                if i == 3 {
                    let frame = Frame3 {
                        origin: [self.center[0], self.center[1], here[1]],
                        x: X,
                        y: cross([0.0, 0.0, out[1].signum()], X),
                        z: [0.0, 0.0, out[1].signum()],
                    };
                    let use_ = plane_boundary(b, frame, circle_a, va, va, Orientation::Forward)?;
                    b.face(
                        &format!("{}.{level}.{i}", self.id),
                        SurfaceGeometry::Plane { frame },
                        [[-chain[3].at(0.0).0[0] - 1.0, chain[3].at(0.0).0[0] + 1.0]; 2],
                        vec![use_],
                    )?;
                    continue;
                }
                let (circle_b, vb) = rings[i + 1];
                self.band(
                    b,
                    &format!("{}.{level}.{i}", self.id),
                    *m,
                    out,
                    (circle_a, va),
                    (circle_b, vb),
                )?;
            }
        }
        Ok(())
    }

    /// A full turn of a corner or the wall between the circles at its two
    /// ends, with a seam along the piece at angle 0.
    fn band(
        &self,
        b: &mut Builder,
        key: &str,
        m: Meridian,
        out: P2,
        (circle_a, va): (u32, u32),
        (circle_b, vb): (u32, u32),
    ) -> Result<(), GeometryError> {
        let [cx, cy] = self.center;
        // The surface, its parameter direction around the axis (u0, du: u
        // = u0 + du * circle parameter), each end's height parameter, and the
        // seam: an edge from end a to end b (or b to a) with its pcurve
        // origin and direction in v.
        let (surface, (u0, du), [ha, hb], seam, seam_v, forward) = match m {
            Meridian::Arc { c, rho, a0, a1 } => {
                let surface = SurfaceGeometry::Torus {
                    frame: Frame3 {
                        origin: [cx, cy, c[1]],
                        x: X,
                        y: Y,
                        z: Z,
                    },
                    major_radius: c[0],
                    minor_radius: rho,
                };
                // The seam arc runs with its angle, from the lower end.
                let seam = b.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: [cx + c[0], cy, c[1]],
                            x: X,
                            y: Z,
                            z: [0.0, -1.0, 0.0],
                        },
                        radius: rho,
                    },
                    Interval::new(a0.min(a1), a0.max(a1))?,
                    false,
                );
                (
                    surface,
                    (0.0, 1.0),
                    [a0, a1],
                    seam,
                    ([0.0, 0.0], [0.0, 1.0]),
                    a0 < a1,
                )
            }
            Meridian::Line { a, b: e } => {
                let len = (e[0] - a[0]).hypot(e[1] - a[1]);
                let dir = [(e[0] - a[0]) / len, (e[1] - a[1]) / len];
                let seam = b.edge(
                    CurveGeometry::Line {
                        origin: self.point(a, 0.0),
                        direction: [dir[0], 0.0, dir[1]],
                    },
                    Interval::new(0.0, len)?,
                    false,
                );
                if dir[0].abs() < EPS {
                    let surface = SurfaceGeometry::Cylinder {
                        frame: Frame3 {
                            origin: [cx, cy, 0.0],
                            x: X,
                            y: Y,
                            z: Z,
                        },
                        radius: a[0],
                    };
                    (
                        surface,
                        (0.0, 1.0),
                        [a[1], e[1]],
                        seam,
                        ([0.0, a[1]], [0.0, dir[1]]),
                        true,
                    )
                } else {
                    // The apex on the axis; the cone opens the way r grows.
                    let apex = a[1] - a[0] / dir[0] * dir[1];
                    let up = (dir[0] * dir[1]).signum();
                    let z = [0.0, 0.0, up];
                    let surface = SurfaceGeometry::Cone {
                        frame: Frame3 {
                            origin: [cx, cy, apex],
                            x: X,
                            y: cross(z, X),
                            z,
                        },
                        semi_angle: (dir[0].abs() / dir[1].abs()).atan(),
                    };
                    let (va_, vb_) = (up * (a[1] - apex), up * (e[1] - apex));
                    let turn = if up > 0.0 { (0.0, 1.0) } else { (TAU, -1.0) };
                    (
                        surface,
                        turn,
                        [va_, vb_],
                        seam,
                        ([0.0, va_], [0.0, up * dir[1]]),
                        true,
                    )
                }
            }
        };
        // Arc seams run from the lower angle: name the ends in edge order.
        let ((ca, pa, ha), (cb, pb, hb)) = if forward {
            ((circle_a, va, ha), (circle_b, vb, hb))
        } else {
            ((circle_b, vb, hb), (circle_a, va, ha))
        };
        let seam_at = |u: f64| uv_line([u + seam_v.0[0], seam_v.0[1]], seam_v.1);
        let u1 = u0 + du * TAU;
        let uses = vec![
            boundary(
                ca,
                pa,
                pa,
                Orientation::Forward,
                uv_line([u0, ha], [du, 0.0]),
            ),
            boundary(seam, pa, pb, Orientation::Forward, seam_at(u1)),
            boundary(
                cb,
                pb,
                pb,
                Orientation::Reverse,
                uv_line([u0, hb], [du, 0.0]),
            ),
            boundary(seam, pb, pa, Orientation::Reverse, seam_at(u0)),
        ];
        let mid = [PI, (ha + hb) / 2.0];
        let normal = surface.normal_at(mid)?;
        let p = surface.point_at(mid)?;
        let angle = (p[1] - cy).atan2(p[0] - cx);
        let outward = [out[0] * angle.cos(), out[0] * angle.sin(), out[1]];
        let face = b.face(key, surface, [[0.0, TAU], [ha.min(hb), ha.max(hb)]], uses)?;
        if dot(normal, outward) < 0.0 {
            b.set_face_sense(face, Orientation::Reverse)?;
        }
        Ok(())
    }
}

/// Over the piece, the line integrals of r^2 / 2 dz and r^2 z / 2 dz.
fn green(m: Meridian) -> (f64, f64) {
    let (mut a, mut b) = (0.0, 0.0);
    let parts = 8;
    for k in 0..parts {
        for (node, w) in GAUSS {
            let s = (k as f64 + node) / parts as f64;
            let (p, d) = m.at(s);
            let w = w / parts as f64;
            a += p[0] * p[0] / 2.0 * d[1] * w;
            b += p[0] * p[0] * p[1] / 2.0 * d[1] * w;
        }
    }
    (a, b)
}

/// Over the piece, the integral of r ds (its area turned once, over 2 pi).
fn length_moment(m: Meridian) -> f64 {
    let mut a = 0.0;
    let parts = 8;
    for k in 0..parts {
        for (node, w) in GAUSS {
            let (p, d) = m.at((k as f64 + node) / parts as f64);
            a += p[0] * d[0].hypot(d[1]) * w / parts as f64;
        }
    }
    a
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

#[cfg(test)]
mod tests {
    use crate::spec::PartSpec;
    use crate::SheetError;

    fn part(emboss: &str, flanges: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{"id": "e", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
            "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                "embosses": [{emboss}]}},
            "flanges": [{flanges}]}}"#
        ))
        .unwrap()
    }

    fn round(at: [f64; 2], diameter: f64, height: f64, draft: f64, side: &str) -> String {
        format!(
            r#"{{"id": "d", "at": {at:?}, "shape": {{"round": {{"diameter": {diameter}}}}},
            "height": {height}, "draft": {draft}, "side": "{side}"}}"#
        )
    }

    /// The solid's mesh agrees with the exact mass, up and down.
    #[test]
    fn a_round_emboss_is_an_exact_solid() {
        let mut volumes = Vec::new();
        for (side, draft) in [("up", 30.0), ("down", 30.0), ("up", 0.0), ("down", 10.0)] {
            let spec = part(&round([50.0, 30.0], 30.0, 8.0, draft, side), "");
            let what = format!("{side} {draft}");
            let solid = crate::fold::fold(&spec).unwrap_or_else(|err| panic!("{what}: {err}"));
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            assert!(
                (mesh.volume - exact.volume).abs() / exact.volume < 2e-5,
                "{what}: {} vs {}",
                mesh.volume,
                exact.volume
            );
            assert!(
                (mesh.area - exact.area).abs() / exact.area < 2e-5,
                "{what}: area {} vs {}",
                mesh.area,
                exact.area
            );
            for j in 0..3 {
                assert!(
                    (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs() < 1e-3,
                    "{what}: {:?} vs {:?}",
                    mesh.center_of_mass,
                    exact.center_of_mass
                );
            }
            volumes.push((exact.volume, exact.center_of_mass[2]));
        }
        // Down is up mirrored about the plate's middle.
        assert!((volumes[0].0 - volumes[1].0).abs() < 1e-9);
        assert!((volumes[0].1 + volumes[1].1 - 2.0).abs() < 1e-9);
        // The plate stays whole in the flat; the footprint is marked.
        let f = crate::flat::flat_pattern(&part(&round([50.0, 30.0], 30.0, 6.0, 30.0, "up"), ""))
            .unwrap();
        assert!((f.area() - 6000.0).abs() < 1e-9 && f.holes.is_empty());
        assert_eq!(f.forms.len(), 1);
    }

    #[test]
    fn what_an_emboss_refuses() {
        let invalid = |e: &str, flanges: &str| {
            matches!(
                crate::flat::flat_pattern(&part(e, flanges)),
                Err(SheetError::Invalid(_))
            )
        };
        // Too small across for its corners, too low, across the outline.
        assert!(invalid(&round([50.0, 30.0], 8.0, 6.0, 30.0, "up"), ""));
        assert!(invalid(&round([50.0, 30.0], 30.0, 3.0, 30.0, "up"), ""));
        assert!(invalid(&round([50.0, 10.0], 30.0, 6.0, 30.0, "up"), ""));
        // Over a bend, or on a flange.
        let flange = r#"{"id": "f", "on": {"base": 2}, "angle": 90, "length": 50,
            "datum": "outer", "position": "outside", "side": "up"}"#;
        assert!(crate::flat::flat_pattern(&part(
            &round([50.0, 30.0], 30.0, 6.0, 30.0, "up"),
            flange
        ))
        .is_ok());
        for y in [55.0, 85.0] {
            assert!(matches!(
                crate::flat::flat_pattern(&part(&round([50.0, y], 30.0, 6.0, 30.0, "up"), flange)),
                Err(SheetError::Unsupported(_))
            ));
        }
    }

    /// The profile meets its definition: the raised face's wall, extended,
    /// meets its top at the diameter's mold line, the height above the
    /// plate's face, at the draft; the sheet keeps its thickness.
    #[test]
    fn the_profile_keeps_its_measures() {
        use super::{Formed, Meridian};
        use crate::spec::Emboss;
        for (side, draft) in [("up", 30.0), ("down", 12.5), ("up", 0.0)] {
            let e: Emboss =
                serde_json::from_str(&round([0.0, 0.0], 40.0, 9.0, draft, side)).unwrap();
            let f = Formed::new(&e, 2.0, 1.5).unwrap();
            // The raised face: the one leaving the plate's face on that side.
            let (raised, other) = if side == "up" { (1, 0) } else { (0, 1) };
            let up = if side == "up" { 1.0 } else { -1.0 };
            let Meridian::Line { a, b } = f.chains[raised][1] else {
                panic!()
            };
            let Meridian::Line { a: top, .. } = f.chains[raised][3] else {
                panic!()
            };
            let plate = if side == "up" { 2.0 } else { 0.0 };
            assert!(
                (top[1] - (plate + up * 9.0)).abs() < 1e-12,
                "{side}: height"
            );
            let s = (top[1] - a[1]) / (b[1] - a[1]);
            let mold = a[0] + s * (b[0] - a[0]);
            assert!((mold - 20.0).abs() < 1e-12, "{side}: mold line at {mold}");
            let slope = ((a[0] - b[0]) / (b[1] - a[1]) * up).atan().to_degrees();
            assert!((slope - draft).abs() < 1e-9, "{side}: draft {slope}");
            for i in 0..4 {
                for k in 1..8 {
                    let p = f.chains[raised][i].at(k as f64 / 8.0).0;
                    let q = f.chains[other][i].at(k as f64 / 8.0).0;
                    let gap = (p[0] - q[0]).hypot(p[1] - q[1]);
                    assert!(
                        (gap - 2.0).abs() < 1e-12,
                        "{side}: piece {i} thickness {gap}"
                    );
                }
            }
        }
    }
}
