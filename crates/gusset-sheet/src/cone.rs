//! A part rolled to a cone or cylinder segment (`base.cone`): its flat
//! pattern, its exact solid and its mass, all in closed form.
//!
//! In the meridian half-plane (r, z) the sheet is a rectangle: the slant line
//! through the two arcs, `L` long, and the thickness square to it. Turned
//! around the z axis by the sweep, its long sides give the sheet's faces
//! (cones, or cylinders), its short sides the ends (cones square to the
//! sheet, or planes), and the two meridian planes close it.

use crate::flat::{Arc, Curve, FlatPattern, Layout, P2};
use crate::fold::{accuracy, Mass, P3};
use crate::spec::{Cone, PartSpec, Radial};
use crate::SheetError;
use opengeometry::brep::{
    plane_boundary, projected_boundary, BrepEnvelope, Builder, CurveGeometry, Frame3,
    GeometryError, Orientation, Surface, SurfaceGeometry,
};
use opengeometry::math::Interval;
use std::f64::consts::PI;

const EPS: f64 = 1e-9;

/// The sheet's rectangle in the meridian half-plane.
struct Section {
    sweep: f64,
    /// Where the arcs' face meets z = 0, (r, z).
    p0: P2,
    /// Unit along the slant, from z = 0 up.
    d: P2,
    /// Unit square to it, away from the axis.
    m: P2,
    l: f64,
    /// The two faces' offsets along m: the arcs' face and the other one.
    o: [f64; 2],
    t: f64,
    /// The neutral surface's offset along m.
    neutral: f64,
}

impl Section {
    fn new(spec: &PartSpec, cone: &Cone) -> Result<Section, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("cone: {what}")));
        let b = &spec.base;
        if !b.outline.is_empty() || b.profile.is_some() {
            return invalid("a base has an outline, a profile or a cone, one of them");
        }
        let [r0, r1] = cone.radii;
        if !(r0 > 0.0 && r1 > 0.0 && cone.height > 0.0) {
            return invalid("needs positive radii and height");
        }
        if !(cone.sweep > 0.0 && cone.sweep < 360.0) {
            return invalid("the sweep must lie between 0 and 360 degrees");
        }
        let features = [
            (!spec.flanges.is_empty(), "flanges"),
            (!spec.contours.is_empty(), "contour flanges"),
            (!spec.hems.is_empty(), "hems"),
            (!spec.folds.is_empty(), "folds"),
            (!spec.jogs.is_empty(), "jogs"),
            (!b.holes.is_empty() || !spec.face_cuts.is_empty(), "holes"),
            (!b.cutouts.is_empty(), "cut-outs"),
            (spec.bend_table.is_some(), "a bend table"),
        ];
        if let Some((_, what)) = features.iter().find(|f| f.0) {
            return Err(SheetError::Unsupported(format!("{what} on a cone")));
        }
        let t = spec.thickness;
        let (dr, h) = (r1 - r0, cone.height);
        let l = dr.hypot(h);
        let side = match cone.material {
            Radial::Outside => t,
            Radial::Inside => -t,
        };
        let s = Section {
            sweep: cone.sweep.to_radians(),
            p0: [r0, 0.0],
            d: [dr / l, h / l],
            m: [h / l, -dr / l],
            l,
            o: [0.0, side],
            t,
            neutral: side.min(0.0) + spec.k_factor * t,
        };
        let corners = [0.0, l].map(|v| s.o.map(|o| s.at(v, o)[0]));
        if corners.iter().flatten().any(|&r| r < EPS) {
            return invalid("the sheet reaches the axis");
        }
        Ok(s)
    }

    /// The meridian point at slant position v and offset o.
    fn at(&self, v: f64, o: f64) -> P2 {
        [
            self.p0[0] + v * self.d[0] + o * self.m[0],
            self.p0[1] + v * self.d[1] + o * self.m[1],
        ]
    }

    /// The developed neutral surface: the flat point of slant position v at
    /// turn angle u.
    fn flat(&self, u: f64, v: f64) -> P2 {
        let (rn0, dr) = (self.at(0.0, self.neutral)[0], self.d[0]);
        if dr.abs() < EPS {
            return [u * rn0, v];
        }
        // The apex lies on the first slant edge, at y = a.
        let a = -rn0 / dr;
        let phi = u * rn0 / a;
        [(a - v) * phi.sin(), a - (a - v) * phi.cos()]
    }

    fn pattern(&self) -> FlatPattern {
        let (b, l) = (self.sweep, self.l);
        let corners = [(0.0, 0.0), (b, 0.0), (b, l), (0.0, l)].map(|(u, v)| self.flat(u, v));
        let arcs = if self.d[0].abs() < EPS {
            vec![None; 4]
        } else {
            let a = -self.at(0.0, self.neutral)[0] / self.d[0];
            let arc = |ccw| {
                Some(Curve::Arc(Arc {
                    center: [0.0, a],
                    ccw,
                }))
            };
            vec![arc(a > 0.0), None, arc(a < 0.0), None]
        };
        FlatPattern {
            outline: corners.to_vec(),
            outline_arcs: arcs,
            holes: Vec::new(),
            cutouts: Vec::new(),
            bends: Vec::new(),
        }
    }

    /// The surface the meridian line through q along dir turns into. Its
    /// angle parameter has the seam opposite the sweep's middle.
    fn surface(&self, q: P2, dir: P2) -> SurfaceGeometry {
        let c = self.sweep / 2.0 + PI;
        let x = [c.cos(), c.sin(), 0.0];
        let frame = |z0: f64, up: f64| Frame3 {
            origin: [0.0, 0.0, z0],
            x,
            y: cross([0.0, 0.0, up], x),
            z: [0.0, 0.0, up],
        };
        if dir[0].abs() < EPS {
            return SurfaceGeometry::Cylinder {
                frame: frame(0.0, 1.0),
                radius: q[0],
            };
        }
        if dir[1].abs() < EPS {
            return SurfaceGeometry::Plane {
                frame: frame(q[1], 1.0),
            };
        }
        // The apex on the axis; the cone opens the way r grows.
        let s = -q[0] / dir[0];
        let up = (dir[1] * dir[0]).signum();
        SurfaceGeometry::Cone {
            frame: frame(q[1] + s * dir[1], up),
            semi_angle: (dir[0].abs() / dir[1].abs()).atan(),
        }
    }

    fn point(&self, q: P2, theta: f64) -> P3 {
        [q[0] * theta.cos(), q[0] * theta.sin(), q[1]]
    }

    fn solid(&self) -> Result<BrepEnvelope, GeometryError> {
        let (b, l) = (self.sweep, self.l);
        let mut builder = Builder::new("gusset".into(), accuracy())?;
        let corner = |th: usize, vi: usize, oi: usize| {
            self.point(self.at(vi as f64 * l, self.o[oi]), th as f64 * b)
        };
        // Vertex [theta end][slant end][face].
        let mut v = [[[0u32; 2]; 2]; 2];
        for (th, vi, oi) in triples() {
            v[th][vi][oi] = builder.vertex(corner(th, vi, oi));
        }
        let line = |builder: &mut Builder, p: P3, q: P3| {
            let length = norm(sub(q, p));
            builder.edge(
                CurveGeometry::Line {
                    origin: p,
                    direction: sub(q, p).map(|x| x / length),
                },
                Interval::new(0.0, length).expect("a positive length"),
                false,
            )
        };
        // Arcs [slant end][face], slant lines [theta end][face], lines
        // through the sheet [theta end][slant end]; each from its lower index.
        let (mut arc, mut slant, mut through) = ([[0u32; 2]; 2], [[0u32; 2]; 2], [[0u32; 2]; 2]);
        for vi in 0..2 {
            for oi in 0..2 {
                let q = self.at(vi as f64 * l, self.o[oi]);
                arc[vi][oi] = builder.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: [0.0, 0.0, q[1]],
                            x: [1.0, 0.0, 0.0],
                            y: [0.0, 1.0, 0.0],
                            z: [0.0, 0.0, 1.0],
                        },
                        radius: q[0],
                    },
                    Interval::new(0.0, b)?,
                    false,
                );
            }
        }
        for th in 0..2 {
            for oi in 0..2 {
                slant[th][oi] = line(&mut builder, corner(th, 0, oi), corner(th, 1, oi));
            }
            for vi in 0..2 {
                through[th][vi] = line(&mut builder, corner(th, vi, 0), corner(th, vi, 1));
            }
        }

        let mid = |vf: f64, of: f64| self.at(vf * l, self.o[0] + of * (self.o[1] - self.o[0]));
        let theta_mid = b / 2.0;
        let toward = |n: P2, theta: f64| [n[0] * theta.cos(), n[0] * theta.sin(), n[1]];
        let neg = |p: P2| [-p[0], -p[1]];
        let outward_o = if self.o[1] > self.o[0] {
            [self.m, neg(self.m)]
        } else {
            [neg(self.m), self.m]
        };
        // Points over a face turned from the meridian segment a..b.
        let turned = |a: P2, c: P2| -> Vec<P3> {
            (0..=32)
                .flat_map(|i| {
                    (0..=2).map(move |j| {
                        let f = j as f64 / 2.0;
                        (i, [a[0] + f * (c[0] - a[0]), a[1] + f * (c[1] - a[1])])
                    })
                })
                .map(|(i, q)| self.point(q, b * i as f64 / 32.0))
                .collect()
        };
        // (key, surface, a point inside, its outward normal, points over it,
        // loop as (edge, from, to, forward)).
        type Loop = Vec<(u32, u32, u32, bool)>;
        let mut faces: Vec<(String, SurfaceGeometry, P3, P3, Vec<P3>, Loop)> = Vec::new();
        for oi in 0..2 {
            let q = self.at(0.0, self.o[oi]);
            let inside = mid(0.5, oi as f64);
            faces.push((
                format!("sheet.{oi}"),
                self.surface(q, self.d),
                self.point(inside, theta_mid),
                toward(outward_o[1 - oi], theta_mid),
                turned(q, self.at(l, self.o[oi])),
                vec![
                    (arc[0][oi], v[0][0][oi], v[1][0][oi], true),
                    (slant[1][oi], v[1][0][oi], v[1][1][oi], true),
                    (arc[1][oi], v[1][1][oi], v[0][1][oi], false),
                    (slant[0][oi], v[0][1][oi], v[0][0][oi], false),
                ],
            ));
        }
        for vi in 0..2 {
            let q = self.at(vi as f64 * l, 0.0);
            let out = if vi == 0 { neg(self.d) } else { self.d };
            faces.push((
                format!("end.{vi}"),
                self.surface(q, self.m),
                self.point(mid(vi as f64, 0.5), theta_mid),
                toward(out, theta_mid),
                turned(q, self.at(vi as f64 * l, self.o[1])),
                vec![
                    (arc[vi][0], v[0][vi][0], v[1][vi][0], true),
                    (through[1][vi], v[1][vi][0], v[1][vi][1], true),
                    (arc[vi][1], v[1][vi][1], v[0][vi][1], false),
                    (through[0][vi], v[0][vi][1], v[0][vi][0], false),
                ],
            ));
        }
        for th in 0..2 {
            let theta = th as f64 * b;
            let normal = if th == 0 {
                [0.0, -1.0, 0.0]
            } else {
                [-theta.sin(), theta.cos(), 0.0]
            };
            let radial = [theta.cos(), theta.sin(), 0.0];
            let frame = Frame3 {
                origin: [0.0; 3],
                x: radial,
                y: cross(normal, radial),
                z: normal,
            };
            faces.push((
                format!("side.{th}"),
                SurfaceGeometry::Plane { frame },
                self.point(mid(0.5, 0.5), theta),
                normal,
                triples()
                    .filter(|c| c.0 == th)
                    .map(|(th, vi, oi)| corner(th, vi, oi))
                    .collect(),
                vec![
                    (slant[th][0], v[th][0][0], v[th][1][0], true),
                    (through[th][1], v[th][1][0], v[th][1][1], true),
                    (slant[th][1], v[th][1][1], v[th][0][1], false),
                    (through[th][0], v[th][0][1], v[th][0][0], false),
                ],
            ));
        }

        for (key, surface, inside, outward, samples, edges) in faces {
            let hint = surface.project(inside, None)?;
            let sense = |forward: bool| {
                if forward {
                    Orientation::Forward
                } else {
                    Orientation::Reverse
                }
            };
            let mut uses = Vec::new();
            for &(e, from, to, forward) in &edges {
                uses.push(match &surface {
                    SurfaceGeometry::Plane { frame } => {
                        plane_boundary(&builder, *frame, e, from, to, sense(forward))?
                    }
                    _ => projected_boundary(
                        &builder,
                        &surface,
                        e,
                        from,
                        to,
                        sense(forward),
                        Some(hint),
                    )?,
                });
            }
            // The face's parameter box.
            let mut bounds = [[f64::MAX, f64::MIN]; 2];
            for &point in &samples {
                let uv = surface.project(point, Some(hint))?;
                for j in 0..2 {
                    bounds[j][0] = bounds[j][0].min(uv[j]);
                    bounds[j][1] = bounds[j][1].max(uv[j]);
                }
            }
            let pad = 1e-6;
            let bounds = bounds.map(|[lo, hi]| [lo - pad, hi + pad]);
            let normal = surface.normal_at(hint)?;
            let face = builder.face(&key, surface, bounds, uses)?;
            if dot(normal, outward) < 0.0 {
                builder.set_face_sense(face, Orientation::Reverse)?;
            }
        }
        builder.finish_solid()
    }

    fn mass(&self) -> Mass {
        let (b, l, t) = (self.sweep, self.l, self.t);
        let c = self.at(l / 2.0, (self.o[0] + self.o[1]) / 2.0);
        let area = l * t;
        // Integrals over the rectangle: of r, r^2 and r z.
        let r = area * c[0];
        let rr =
            area * (c[0] * c[0] + (self.d[0] * l).powi(2) / 12.0 + (self.m[0] * t).powi(2) / 12.0);
        let rz = area
            * (c[0] * c[1]
                + self.d[0] * self.d[1] * l * l / 12.0
                + self.m[0] * self.m[1] * t * t / 12.0);
        let volume = b * r;
        // Each face turned: its length times the sweep times its mid radius.
        let faces = self
            .o
            .iter()
            .map(|&o| l * self.at(l / 2.0, o)[0])
            .sum::<f64>()
            + [0.0, l]
                .iter()
                .map(|&v| t * self.at(v, (self.o[0] + self.o[1]) / 2.0)[0])
                .sum::<f64>();
        Mass {
            volume,
            area: b * faces + 2.0 * area,
            center_of_mass: [
                rr * b.sin() / volume,
                rr * (1.0 - b.cos()) / volume,
                b * rz / volume,
            ],
        }
    }
}

fn triples() -> impl Iterator<Item = (usize, usize, usize)> {
    (0..8).map(|i| (i >> 2, (i >> 1) & 1, i & 1))
}

fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}
fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn section(spec: &PartSpec) -> Option<Result<Section, SheetError>> {
    spec.base.cone.as_ref().map(|c| Section::new(spec, c))
}

/// The layout of a cone part: its flat pattern, nothing to fold.
pub(crate) fn layout(spec: &PartSpec) -> Option<Result<Layout, SheetError>> {
    section(spec).map(|s| {
        Ok(Layout {
            flat: s?.pattern(),
            cells: Vec::new(),
            folds: Vec::new(),
            holes: Vec::new(),
        })
    })
}

pub(crate) fn fold(spec: &PartSpec) -> Option<Result<BrepEnvelope, SheetError>> {
    section(spec).map(|s| Ok(s?.solid()?))
}

pub(crate) fn mass(spec: &PartSpec) -> Option<Result<Mass, SheetError>> {
    section(spec).map(|s| Ok(s?.mass()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(radii: [f64; 2], height: f64, sweep: f64, material: Radial) -> PartSpec {
        let mut spec: PartSpec = serde_json::from_str(
            r#"{"id": "cone", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, "base": {}}"#,
        )
        .unwrap();
        spec.base.cone = Some(Cone {
            radii,
            height,
            sweep,
            material,
        });
        spec
    }

    /// Flat area, volume and surface area as a lofted flange between the
    /// same arcs measures them.
    #[test]
    fn cones_and_cylinders_unroll_on_their_neutral_surface() {
        use Radial::*;
        let cases = [
            // r 50 to 30 over 60, half a turn.
            (
                part([50.0, 30.0], 60.0, 180.0, Outside),
                8098.467060,
                16272.332344,
                17039.890887,
            ),
            (
                part([50.0, 30.0], 60.0, 180.0, Inside),
                7721.475942,
                15518.350107,
                16262.065638,
            ),
            (
                part([30.0, 50.0], 60.0, 120.0, Outside),
                5398.978040,
                10848.221562,
                11444.254662,
            ),
            (
                part([30.0, 50.0], 60.0, 120.0, Inside),
                5147.650628,
                10345.566738,
                10925.704496,
            ),
            (
                part([50.0, 30.0], 60.0, 270.0, Outside),
                12147.700590,
                24408.498516,
                25433.345224,
            ),
            (
                part([50.0, 10.0], 20.0, 120.0, Outside),
                2843.436214,
                5703.627589,
                6137.587007,
            ),
            (
                part([50.0, 20.0], 1.0, 90.0, Outside),
                1651.508824,
                3303.645966,
                3643.833423,
            ),
            (
                part([50.0, 50.0], 60.0, 180.0, Outside),
                9575.574408,
                19226.547040,
                20107.431941,
            ),
        ];
        for (spec, flat, volume, area) in cases {
            let what = format!("{:?}", spec.base.cone);
            let f = crate::flat::flat_pattern(&spec).unwrap();
            assert!((f.area() - flat).abs() < 1e-5, "{what}: flat {}", f.area());
            let m = crate::fold::exact_mass(&spec).unwrap();
            assert!(
                (m.volume - volume).abs() < 1e-5,
                "{what}: volume {}",
                m.volume
            );
            assert!((m.area - area).abs() < 1e-5, "{what}: area {}", m.area);
            let solid = crate::fold::fold(&spec).unwrap();
            let meshed = crate::fold::mass(&solid, 1e-3).unwrap();
            assert!(
                (meshed.volume - volume).abs() / volume < 1e-4,
                "{what}: meshed {}",
                meshed.volume
            );
            for j in 0..3 {
                assert!(
                    (meshed.center_of_mass[j] - m.center_of_mass[j]).abs() < 1e-2,
                    "{what}: centre {:?} vs {:?}",
                    meshed.center_of_mass,
                    m.center_of_mass
                );
            }
        }
    }

    #[test]
    fn the_first_slant_edge_runs_up_the_y_axis() {
        let spec = part([50.0, 30.0], 60.0, 180.0, Radial::Outside);
        let f = crate::flat::flat_pattern(&spec).unwrap();
        let l = 20f64.hypot(60.0);
        assert!(f.outline[0][0].abs() < 1e-12 && f.outline[0][1].abs() < 1e-12);
        assert!(f.outline[3][0].abs() < 1e-12 && (f.outline[3][1] - l).abs() < 1e-12);
        // The bottom edge is the neutral arc, rn0 = 50 + 0.8 cos(alpha), long.
        let rn0 = 50.0 + 0.8 * 60.0 / l;
        let Some(Curve::Arc(arc)) = f.outline_arcs[0] else {
            panic!()
        };
        let r = (arc.center[1]).abs();
        assert!((r * arc.sweep(f.outline[0], f.outline[1]).abs() - PI * rn0).abs() < 1e-9);
    }

    #[test]
    fn a_cone_passes_the_shop_checks_and_has_no_clashes() {
        let spec = part([50.0, 30.0], 60.0, 180.0, Radial::Outside);
        assert!(crate::check::check(&spec, &Default::default())
            .unwrap()
            .is_empty());
        assert!(crate::clash::clashes(&spec).unwrap().is_empty());
        assert!(crate::clash::preview(&spec).is_ok());
    }

    #[test]
    fn what_a_cone_refuses() {
        let mut spec = part([50.0, 30.0], 60.0, 180.0, Radial::Inside);
        spec.base.cone.as_mut().unwrap().radii = [1.0, 30.0];
        assert!(matches!(
            crate::flat::flat_pattern(&spec),
            Err(SheetError::Invalid(_))
        ));
        let mut spec = part([50.0, 30.0], 60.0, 360.0, Radial::Outside);
        assert!(matches!(
            crate::flat::flat_pattern(&spec),
            Err(SheetError::Invalid(_))
        ));
        spec.base.cone.as_mut().unwrap().sweep = 90.0;
        spec.base.holes.push(crate::spec::Hole {
            center: [10.0, 10.0],
            diameter: 5.0,
        });
        assert!(matches!(
            crate::fold::fold(&spec),
            Err(SheetError::Unsupported(_))
        ));
    }
}
