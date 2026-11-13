//! A part rolled to a cone or cylinder segment (`base.cone`): its flat
//! pattern, its exact solid and its mass, all in closed form.
//!
//! In the meridian half-plane (r, z) the sheet is a rectangle: the slant line
//! through the two arcs, `L` long, and the thickness square to it. Turned
//! around the z axis by the sweep, its long sides give the sheet's faces
//! (cones, or cylinders), its short sides the ends (cones square to the
//! sheet, or planes), and the two meridian planes close it.

use crate::flat::{Arc, Cell, Curve, FlatPattern, FlatPunch, Layout, Ring, Zone, P2};
use crate::fold::{accuracy, Mass, GAUSS, P3};
use crate::spec::{Cone, PartSpec, Radial};
use crate::SheetError;
use opengeometry::brep::{
    boundary, plane_boundary, projected_boundary, uv_line, BrepEnvelope, Builder, CurveGeometry,
    Frame3, GeometryError, Orientation, Surface, SurfaceGeometry, Use, Wrap,
};
use opengeometry::math::Interval;
use std::f64::consts::{PI, TAU};

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
            (!spec.face_cuts.is_empty(), "face cuts"),
            (
                b.punches.iter().any(|p| p.on.is_some()),
                "a punch on a face",
            ),
            (!b.embosses.is_empty(), "embosses"),
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

    /// The neutral surface's development as an OpenGeometry wrap, its z the
    /// neutral radius at the first arc plus the offset along m. Its frame is
    /// turned to the middle of the sweep (its inverse takes angles within
    /// half a turn of that); `to_wrap` takes flat points there.
    fn wrap(&self) -> Wrap {
        let n0 = self.at(0.0, self.neutral);
        let (s, c) = (self.sweep / 2.0).sin_cos();
        Wrap {
            frame: Frame3 {
                origin: [0.0, 0.0, n0[1]],
                x: [c, s, 0.0],
                y: [-s, c, 0.0],
                z: [0.0, 0.0, 1.0],
            },
            radius: n0[0],
            semi_angle: if self.d[0].abs() < EPS {
                0.0
            } else {
                self.d[0].atan2(self.d[1])
            },
        }
    }

    /// A flat point in the wrap's flat: turned about the apex by half the
    /// sweep times sin(semi-angle), or for a cylinder moved back by half the
    /// sweep's arc.
    fn to_wrap(&self, p: P2) -> P2 {
        let r = self.at(0.0, self.neutral)[0];
        let s = self.d[0];
        if s.abs() < EPS {
            return [p[0] - self.sweep / 2.0 * r, p[1]];
        }
        let apex = [0.0, -r / s];
        let v = self.to_wrap_dir([p[0] - apex[0], p[1] - apex[1]]);
        [apex[0] + v[0], apex[1] + v[1]]
    }

    /// A flat direction in the wrap's flat.
    fn to_wrap_dir(&self, v: P2) -> P2 {
        if self.d[0].abs() < EPS {
            return v;
        }
        let (sn, cs) = (self.sweep / 2.0 * self.d[0]).sin_cos();
        [cs * v[0] - sn * v[1], sn * v[0] + cs * v[1]]
    }

    /// The wrap's z of face oi.
    fn level(&self, oi: usize) -> f64 {
        self.at(0.0, self.neutral)[0] + self.o[oi] - self.neutral
    }

    /// At a flat point: the turn angle, the neutral radius, and the unit
    /// direction away from the apex (along x, along y).
    fn developed(&self, p: P2) -> (f64, f64, f64, f64) {
        let r = self.at(0.0, self.neutral)[0];
        let s = self.d[0];
        if s.abs() < EPS {
            return (p[0] / r, r, 0.0, 1.0);
        }
        let k = s / r;
        let (a, b) = (1.0 + k * p[1], k * p[0]);
        let e = a.hypot(b);
        (b.atan2(a) / s, e * r, b / e, a / e)
    }

    /// The flat pattern with its holes, cut-outs and punch stamps, which
    /// must lie inside the sheet, clear of its edges and each other.
    fn flat_pattern(&self, spec: &PartSpec) -> Result<FlatPattern, SheetError> {
        let mut f = self.pattern();
        let mut holes = spec.base.holes.clone();
        let mut cutouts: Vec<Ring> = spec.base.cutouts.iter().map(Ring::from_cutout).collect();
        for p in &spec.base.punches {
            if p.at.is_empty() || p.holes.is_empty() && p.cutouts.is_empty() {
                return Err(SheetError::Invalid(format!(
                    "punch {}: needs a shape and a place",
                    p.id
                )));
            }
            for &at in &p.at {
                let (h, k) = p.stamped(at);
                holes.extend(h);
                cutouts.extend(k.iter().map(Ring::from_cutout));
                f.punches.push(FlatPunch {
                    id: p.id.clone(),
                    center: at,
                    angle: p.angle,
                });
            }
        }
        if holes.iter().any(|h| crate::cut::crosses_outline(&f, h)) {
            return Err(SheetError::Unsupported(
                "a hole across a cone's edge".into(),
            ));
        }
        let n = cutouts.len();
        f.holes = holes.clone();
        f.cutouts = cutouts;
        let mut layout = Layout {
            cells: vec![Cell {
                name: "sheet".into(),
                polygon: f.outline.clone(),
                arcs: f.outline_arcs.clone(),
                holes: Vec::new(),
                inner: Vec::new(),
                zone: Zone::Base,
            }],
            flat: f,
            folds: Vec::new(),
            holes: holes.clone(),
            embosses: Vec::new(),
        };
        crate::cut::cut_holes(&mut layout, &holes)?;
        crate::cut::cut_cutouts(&mut layout)?;
        if layout.flat.cutouts.len() != n {
            return Err(SheetError::Unsupported(
                "a cut-out across a cone's edge".into(),
            ));
        }
        // Their edges go onto the sheet's faces by projection, which takes
        // them only within a quarter turn around the axis.
        let span = |points: Vec<P2>| {
            let turns = points.iter().map(|&p| self.developed(p).0);
            let (lo, hi) = turns.fold((f64::MAX, f64::MIN), |(lo, hi), t| (lo.min(t), hi.max(t)));
            hi - lo
        };
        let circles = layout.flat.holes.iter().map(|h| {
            (0..64)
                .map(|i| {
                    let a = TAU * i as f64 / 64.0;
                    let r = h.diameter / 2.0;
                    [h.center[0] + r * a.cos(), h.center[1] + r * a.sin()]
                })
                .collect::<Vec<_>>()
        });
        let rings = layout.flat.cutouts.iter().map(|ring| {
            let n = ring.points.len();
            (0..n)
                .flat_map(|k| {
                    let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
                    let arc = ring.arcs[k].and_then(Curve::arc);
                    (0..16).map(move |i| {
                        let f = i as f64 / 16.0;
                        match arc {
                            None => [p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1])],
                            Some(arc) => {
                                let c = arc.center;
                                let r = (p[0] - c[0]).hypot(p[1] - c[1]);
                                let a = (p[1] - c[1]).atan2(p[0] - c[0]) + f * arc.sweep(p, q);
                                [c[0] + r * a.cos(), c[1] + r * a.sin()]
                            }
                        }
                    })
                })
                .collect::<Vec<_>>()
        });
        if circles.chain(rings).any(|points| span(points) > 0.45 * PI) {
            return Err(SheetError::Unsupported(
                "a hole or cut-out reaching a quarter turn around a cone's axis".into(),
            ));
        }
        Ok(layout.flat)
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
            punches: Vec::new(),
            forms: Vec::new(),
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

    fn solid(&self, flat: &FlatPattern) -> Result<BrepEnvelope, GeometryError> {
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

        let (inner, walls) = self.cuts(&mut builder, flat)?;

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
            // The holes and cut-outs on a face of the sheet, clockwise in the flat.
            let mut holes = Vec::new();
            if let Some(oi) = key
                .strip_prefix("sheet.")
                .map(|k| k.parse::<usize>().unwrap())
            {
                for ring in &inner[oi] {
                    let mut uses = Vec::new();
                    for &(e, from, to, sense, at) in ring {
                        let near = surface.project(at, Some(hint))?;
                        uses.push(projected_boundary(
                            &builder,
                            &surface,
                            e,
                            from,
                            to,
                            sense,
                            Some(near),
                        )?);
                    }
                    holes.push(uses);
                }
            }
            let face = builder.face_with_holes(&key, surface, bounds, uses, holes)?;
            if dot(normal, outward) < 0.0 {
                builder.set_face_sense(face, Orientation::Reverse)?;
            }
        }
        for (key, surface, bounds, uses, mid, into) in walls {
            let normal = surface.normal_at(mid)?;
            let face = builder.face(&key, surface, bounds, uses)?;
            // The wall faces out of the material, into the cut.
            if dot(normal, into) < 0.0 {
                builder.set_face_sense(face, Orientation::Reverse)?;
            }
        }
        builder.finish_solid()
    }

    /// The edges of the holes and cut-outs on both faces of the sheet, and
    /// their walls: per face, the inner loops as (edge, from, to, sense, a
    /// point at from); the walls as (key, surface, bounds, uses, a UV in
    /// the middle, the way into the cut there).
    #[allow(clippy::type_complexity)]
    fn cuts(
        &self,
        builder: &mut Builder,
        flat: &FlatPattern,
    ) -> Result<
        (
            [Vec<Vec<(u32, u32, u32, Orientation, P3)>>; 2],
            Vec<(
                String,
                SurfaceGeometry,
                [[f64; 2]; 2],
                Vec<Use>,
                [f64; 2],
                P3,
            )>,
        ),
        GeometryError,
    > {
        let wrap = self.wrap();
        let t = self.t;
        let z = [self.level(0), self.level(1)];
        let (lo, hi) = if z[0] < z[1] { (0, 1) } else { (1, 0) };
        let to = |p: P2| self.to_wrap(p);
        let dir = |v: P2| self.to_wrap_dir(v);
        let at = |p: P2, oi: usize| {
            let q = to(p);
            wrap.point([q[0], q[1], z[oi]])
        };
        // The way into the cut from flat point p along flat direction d, mid-sheet.
        let into = |p: P2, d: P2| {
            let mid = (z[0] + z[1]) / 2.0;
            let h = 1e-6;
            let (q, e) = (to(p), dir(d));
            sub(
                wrap.point([q[0] + h * e[0], q[1] + h * e[1], mid]),
                wrap.point([q[0], q[1], mid]),
            )
        };
        let flat_frame = |c: P2, zz: f64| {
            let (o, x, y) = (to(c), dir([1.0, 0.0]), dir([0.0, 1.0]));
            Frame3 {
                origin: [o[0], o[1], zz],
                x: [x[0], x[1], 0.0],
                y: [y[0], y[1], 0.0],
                z: [0.0, 0.0, 1.0],
            }
        };
        let wrapped = |curve: CurveGeometry| CurveGeometry::Wrapped {
            wrap,
            curve: Box::new(curve),
        };
        let vertical = |builder: &mut Builder, p: P2, seam: bool| -> Result<u32, GeometryError> {
            let (a, c) = (at(p, lo), at(p, hi));
            Ok(builder.edge(
                CurveGeometry::Line {
                    origin: a,
                    direction: sub(c, a).map(|x| x / norm(sub(c, a))),
                },
                Interval::new(0.0, t)?,
                seam,
            ))
        };
        let mut inner: [Vec<Vec<_>>; 2] = [Vec::new(), Vec::new()];
        let mut walls = Vec::new();

        for (k, h) in flat.holes.iter().enumerate() {
            let r = h.diameter / 2.0;
            let sp = [h.center[0] + r, h.center[1]];
            let v = [0, 1].map(|oi| builder.vertex(at(sp, oi)));
            let mut circle = [0; 2];
            for oi in 0..2 {
                circle[oi] = builder.edge(
                    wrapped(CurveGeometry::Circle {
                        frame: flat_frame(h.center, z[oi]),
                        radius: r,
                    }),
                    Interval::new(0.0, TAU)?,
                    false,
                );
                inner[oi].push(vec![(
                    circle[oi],
                    v[oi],
                    v[oi],
                    Orientation::Reverse,
                    at(sp, oi),
                )]);
            }
            let seam = vertical(builder, sp, true)?;
            walls.push((
                format!("hole.{k}"),
                SurfaceGeometry::Wrapped {
                    wrap,
                    surface: Box::new(SurfaceGeometry::Cylinder {
                        frame: flat_frame(h.center, z[lo]),
                        radius: r,
                    }),
                },
                [[0.0, TAU], [0.0, t]],
                vec![
                    boundary(
                        circle[lo],
                        v[lo],
                        v[lo],
                        Orientation::Forward,
                        uv_line([0.0, 0.0], [1.0, 0.0]),
                    ),
                    boundary(
                        seam,
                        v[lo],
                        v[hi],
                        Orientation::Forward,
                        uv_line([TAU, 0.0], [0.0, 1.0]),
                    ),
                    boundary(
                        circle[hi],
                        v[hi],
                        v[hi],
                        Orientation::Reverse,
                        uv_line([0.0, t], [1.0, 0.0]),
                    ),
                    boundary(
                        seam,
                        v[hi],
                        v[lo],
                        Orientation::Reverse,
                        uv_line([0.0, 0.0], [0.0, 1.0]),
                    ),
                ],
                [PI, t / 2.0],
                into([h.center[0] - r, h.center[1]], [1.0, 0.0]),
            ));
        }

        for (k, ring) in flat.cutouts.iter().enumerate() {
            let n = ring.points.len();
            let v: Vec<[u32; 2]> = ring
                .points
                .iter()
                .map(|&p| [0, 1].map(|oi| builder.vertex(at(p, oi))))
                .collect();
            let mut verticals = Vec::with_capacity(n);
            for &p in &ring.points {
                verticals.push(vertical(builder, p, false)?);
            }
            // Per segment: its edge on each face, and whether it runs p -> q.
            let mut segments = Vec::with_capacity(n);
            for i in 0..n {
                let j = (i + 1) % n;
                let (p, q) = (ring.points[i], ring.points[j]);
                let (edges, forward) = match ring.arcs[i].and_then(Curve::arc) {
                    None => {
                        let d = [q[0] - p[0], q[1] - p[1]];
                        let len = d[0].hypot(d[1]);
                        let u = [d[0] / len, d[1] / len];
                        let (o, uw) = (to(p), dir(u));
                        let mut edges = [0; 2];
                        for oi in 0..2 {
                            edges[oi] = builder.edge(
                                wrapped(CurveGeometry::Line {
                                    origin: [o[0], o[1], z[oi]],
                                    direction: [uw[0], uw[1], 0.0],
                                }),
                                Interval::new(0.0, len)?,
                                false,
                            );
                        }
                        walls.push((
                            format!("cut.{k}.{i}"),
                            SurfaceGeometry::Wrapped {
                                wrap,
                                surface: Box::new(SurfaceGeometry::Plane {
                                    frame: Frame3 {
                                        origin: [o[0], o[1], z[lo]],
                                        x: [uw[0], uw[1], 0.0],
                                        y: [0.0, 0.0, 1.0],
                                        z: [uw[1], -uw[0], 0.0],
                                    },
                                }),
                            },
                            [[0.0, len], [0.0, t]],
                            vec![
                                boundary(
                                    edges[lo],
                                    v[i][lo],
                                    v[j][lo],
                                    Orientation::Forward,
                                    uv_line([0.0, 0.0], [1.0, 0.0]),
                                ),
                                boundary(
                                    verticals[j],
                                    v[j][lo],
                                    v[j][hi],
                                    Orientation::Forward,
                                    uv_line([len, 0.0], [0.0, 1.0]),
                                ),
                                boundary(
                                    edges[hi],
                                    v[j][hi],
                                    v[i][hi],
                                    Orientation::Reverse,
                                    uv_line([0.0, t], [1.0, 0.0]),
                                ),
                                boundary(
                                    verticals[i],
                                    v[i][hi],
                                    v[i][lo],
                                    Orientation::Reverse,
                                    uv_line([0.0, 0.0], [0.0, 1.0]),
                                ),
                            ],
                            [len / 2.0, t / 2.0],
                            into([p[0] + d[0] / 2.0, p[1] + d[1] / 2.0], [-u[1], u[0]]),
                        ));
                        (edges, true)
                    }
                    Some(arc) => {
                        let c = arc.center;
                        let r = (p[0] - c[0]).hypot(p[1] - c[1]);
                        let a_p = (p[1] - c[1]).atan2(p[0] - c[0]);
                        let sweep = arc.sweep(p, q);
                        // The edge runs with the angle, from s to e.
                        let (a0, a1, s, e) = if sweep > 0.0 {
                            (a_p, a_p + sweep, i, j)
                        } else {
                            (a_p + sweep, a_p, j, i)
                        };
                        let mut edges = [0; 2];
                        for oi in 0..2 {
                            edges[oi] = builder.edge(
                                wrapped(CurveGeometry::Circle {
                                    frame: flat_frame(c, z[oi]),
                                    radius: r,
                                }),
                                Interval::new(a0, a1)?,
                                false,
                            );
                        }
                        let am = (a0 + a1) / 2.0;
                        let (sm, cm) = am.sin_cos();
                        // Left of p -> q there.
                        let left = if sweep > 0.0 { [-cm, -sm] } else { [cm, sm] };
                        walls.push((
                            format!("cut.{k}.{i}"),
                            SurfaceGeometry::Wrapped {
                                wrap,
                                surface: Box::new(SurfaceGeometry::Cylinder {
                                    frame: flat_frame(c, z[lo]),
                                    radius: r,
                                }),
                            },
                            [[a0, a1], [0.0, t]],
                            vec![
                                boundary(
                                    edges[lo],
                                    v[s][lo],
                                    v[e][lo],
                                    Orientation::Forward,
                                    uv_line([0.0, 0.0], [1.0, 0.0]),
                                ),
                                boundary(
                                    verticals[e],
                                    v[e][lo],
                                    v[e][hi],
                                    Orientation::Forward,
                                    uv_line([a1, 0.0], [0.0, 1.0]),
                                ),
                                boundary(
                                    edges[hi],
                                    v[e][hi],
                                    v[s][hi],
                                    Orientation::Reverse,
                                    uv_line([0.0, t], [1.0, 0.0]),
                                ),
                                boundary(
                                    verticals[s],
                                    v[s][hi],
                                    v[s][lo],
                                    Orientation::Reverse,
                                    uv_line([a0, 0.0], [0.0, 1.0]),
                                ),
                            ],
                            [am, t / 2.0],
                            into([c[0] + r * cm, c[1] + r * sm], left),
                        ));
                        (edges, sweep > 0.0)
                    }
                };
                segments.push((edges, forward));
            }
            for oi in 0..2 {
                let ring_uses = (0..n)
                    .rev()
                    .map(|i| {
                        let j = (i + 1) % n;
                        let (edges, forward) = segments[i];
                        let sense = if forward {
                            Orientation::Reverse
                        } else {
                            Orientation::Forward
                        };
                        (edges[oi], v[j][oi], v[i][oi], sense, at(ring.points[j], oi))
                    })
                    .collect();
                inner[oi].push(ring_uses);
            }
        }
        Ok((inner, walls))
    }

    /// The mass of the sheet less its holes and cut-outs, with their walls.
    fn mass(&self, flat: &FlatPattern) -> Mass {
        let whole = self.plain_mass();
        let mut volume = whole.volume;
        let mut area = whole.area;
        let mut moment = whole.center_of_mass.map(|c| c * whole.volume);
        let holes = flat
            .holes
            .iter()
            .map(|h| vec![Piece::Arc(h.center, h.diameter / 2.0, 0.0, TAU)]);
        let rings = flat.cutouts.iter().map(|ring| {
            let n = ring.points.len();
            (0..n)
                .map(|k| {
                    let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
                    match ring.arcs[k].and_then(Curve::arc) {
                        None => Piece::Line(p, q),
                        Some(arc) => {
                            let c = arc.center;
                            Piece::Arc(
                                c,
                                (p[0] - c[0]).hypot(p[1] - c[1]),
                                (p[1] - c[1]).atan2(p[0] - c[0]),
                                arc.sweep(p, q),
                            )
                        }
                    }
                })
                .collect::<Vec<_>>()
        });
        for pieces in holes.chain(rings) {
            let (cut, wall) = self.removed(&pieces);
            volume -= cut[0];
            area += wall - cut[1] - cut[2];
            for j in 0..3 {
                moment[j] -= cut[3 + j];
            }
        }
        Mass {
            volume,
            area,
            center_of_mass: moment.map(|m| m / volume),
        }
    }

    /// Per unit of flat area at p: the sheet's volume, the areas of its two
    /// faces, and its first moment.
    fn density(&self, wrap: &Wrap, p: P2) -> [f64; 6] {
        let (s, c) = (self.d[0], self.d[1]);
        let (theta, rn, _, _) = self.developed(p);
        let w = self.o.map(|o| o - self.neutral);
        let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
        let w2 = (hi * hi - lo * lo) / 2.0;
        let w3 = (hi.powi(3) - lo.powi(3)) / 3.0;
        let volume = self.t + c * w2 / rn;
        let q = self.to_wrap(p);
        let neutral = wrap.point([q[0], q[1], wrap.radius]);
        let normal = [c * theta.cos(), c * theta.sin(), -s];
        let off = w2 + c * w3 / rn;
        [
            volume,
            1.0 + c * w[0] / rn,
            1.0 + c * w[1] / rn,
            neutral[0] * volume + normal[0] * off,
            neutral[1] * volume + normal[1] * off,
            neutral[2] * volume + normal[2] * off,
        ]
    }

    /// What a counter-clockwise flat region takes from the sheet (its
    /// integrated density, by Green's theorem with Gauss-Legendre along
    /// its boundary and across it), and the area of its wall.
    fn removed(&self, pieces: &[Piece]) -> ([f64; 6], f64) {
        let wrap = self.wrap();
        let start = pieces[0].at(0.0).0;
        let h = 0.05 * self.developed(start).1;
        let x0 = start[0];
        let w = self.o.map(|o| o - self.neutral);
        let (lo, hi) = (w[0].min(w[1]), w[0].max(w[1]));
        let c = self.d[1];
        let (mut total, mut wall) = ([0.0; 6], 0.0);
        for piece in pieces {
            let parts = piece.parts(h);
            for part in 0..parts {
                for (node, weight) in GAUSS {
                    let (p, d) = piece.at((part as f64 + node) / parts as f64);
                    let d = d.map(|x| x / parts as f64);
                    // The density integrated along x from x0 to p.
                    let span = p[0] - x0;
                    let chunks = (span.abs() / h).ceil().max(1.0) as usize;
                    for chunk in 0..chunks {
                        for (inner, iw) in GAUSS {
                            let x = x0 + span * (chunk as f64 + inner) / chunks as f64;
                            let f = self.density(&wrap, [x, p[1]]);
                            let scale = span / chunks as f64 * iw * d[1] * weight;
                            for j in 0..6 {
                                total[j] += f[j] * scale;
                            }
                        }
                    }
                    let (_, rn, sigma, gamma) = self.developed(p);
                    let along = sigma * d[0] + gamma * d[1];
                    let across = gamma * d[0] - sigma * d[1];
                    for (wn, ww) in GAUSS {
                        let lambda = 1.0 + c * (lo + (hi - lo) * wn) / rn;
                        wall += (along * along + (lambda * across).powi(2)).sqrt()
                            * (hi - lo)
                            * ww
                            * weight;
                    }
                }
            }
        }
        (total, wall)
    }

    fn plain_mass(&self) -> Mass {
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

/// A piece of a flat boundary: a line, or an arc (centre, radius, start
/// angle, signed sweep).
#[derive(Clone, Copy)]
enum Piece {
    Line(P2, P2),
    Arc(P2, f64, f64, f64),
}

impl Piece {
    /// Point and derivative at s in [0, 1].
    fn at(&self, s: f64) -> (P2, P2) {
        match *self {
            Piece::Line(p, q) => (
                [p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])],
                [q[0] - p[0], q[1] - p[1]],
            ),
            Piece::Arc(c, r, a0, sweep) => {
                let (sn, cs) = (a0 + s * sweep).sin_cos();
                (
                    [c[0] + r * cs, c[1] + r * sn],
                    [-r * sweep * sn, r * sweep * cs],
                )
            }
        }
    }

    /// How many parts to integrate it in: none longer than h, no arc part
    /// over 22.5 degrees.
    fn parts(&self, h: f64) -> usize {
        let (length, turn) = match *self {
            Piece::Line(p, q) => ((q[0] - p[0]).hypot(q[1] - p[1]), 0.0),
            Piece::Arc(_, r, _, sweep) => (r * sweep.abs(), sweep.abs() / (PI / 8.0)),
        };
        (length / h).ceil().max(turn.ceil()).max(1.0) as usize
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
            flat: s?.flat_pattern(spec)?,
            cells: Vec::new(),
            folds: Vec::new(),
            holes: Vec::new(),
            embosses: Vec::new(),
        })
    })
}

pub(crate) fn fold(spec: &PartSpec) -> Option<Result<BrepEnvelope, SheetError>> {
    section(spec).map(|s| {
        let s = s?;
        Ok(s.solid(&s.flat_pattern(spec)?)?)
    })
}

pub(crate) fn mass(spec: &PartSpec) -> Option<Result<Mass, SheetError>> {
    section(spec).map(|s| {
        let s = s?;
        Ok(s.mass(&s.flat_pattern(spec)?))
    })
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

    /// A hole's measures depend only on its distance from the apex: the
    /// reference CAD's hole of 10 at 133.07 from it, cut in the unfolded
    /// part and refolded.
    #[test]
    fn a_hole_in_a_cone_takes_what_the_reference_cad_takes() {
        let mut spec = part([50.0, 30.0], 60.0, 180.0, Radial::Outside);
        let rho = (25.298221281f64 + 97.268329805).hypot(51.822165126);
        let l = 20f64.hypot(60.0);
        let a = (50.0 + 0.8 * 60.0 / l) * l / 20.0;
        let phi = 0.4f64;
        spec.base.holes.push(crate::spec::Hole {
            center: [rho * phi.sin(), a - rho * phi.cos()],
            diameter: 10.0,
        });
        let f = crate::flat::flat_pattern(&spec).unwrap();
        assert!((f.area() - 8019.927244).abs() < 1e-5, "flat {}", f.area());
        let m = crate::fold::exact_mass(&spec).unwrap();
        assert!(
            (m.volume - 16114.544338).abs() < 1e-5,
            "volume {}",
            m.volume
        );
        assert!((m.area - 16945.077151).abs() < 1e-5, "area {}", m.area);
    }

    /// Holes, punches and cut-outs build and mesh (the reference CAD checks
    /// them exactly, in the harness); a cut-out drawn as a circle takes what
    /// the hole takes; on a cylinder a hole takes pi r^2 T, stretched by
    /// the mid-sheet radius over the neutral one.
    #[test]
    fn cut_cones_and_cylinders_take_what_they_cut() {
        use crate::spec::{Cutout, Hole, Punch};
        let punch = Punch {
            id: "P".into(),
            on: None,
            at: vec![[60.0, 35.0], [72.0, 35.0]],
            angle: 30.0,
            holes: vec![Hole {
                center: [0.0, 0.0],
                diameter: 4.0,
            }],
            cutouts: vec![Cutout {
                outline: vec![[-1.0, 3.0], [1.0, 3.0], [1.0, 6.0], [-1.0, 6.0]],
                bulges: Vec::new(),
            }],
        };
        for (radii, material) in [
            ([50.0, 30.0], Radial::Outside),
            ([30.0, 50.0], Radial::Inside),
            ([50.0, 50.0], Radial::Outside),
        ] {
            let what = format!("{radii:?} {material:?}");
            let mut holed = part(radii, 60.0, 180.0, material);
            holed.base.holes.push(Hole {
                center: [30.0, 30.0],
                diameter: 8.0,
            });
            let mut ringed = holed.clone();
            ringed.base.holes.clear();
            ringed.base.cutouts.push(Cutout {
                outline: vec![[34.0, 30.0], [26.0, 30.0]],
                bulges: vec![1.0, 1.0],
            });
            let (a, b) = (
                crate::fold::exact_mass(&holed).unwrap(),
                crate::fold::exact_mass(&ringed).unwrap(),
            );
            assert!((a.volume - b.volume).abs() < 1e-9, "{what}: {a:?} {b:?}");
            assert!((a.area - b.area).abs() < 1e-9, "{what}: {a:?} {b:?}");
            for j in 0..3 {
                assert!(
                    (a.center_of_mass[j] - b.center_of_mass[j]).abs() < 1e-9,
                    "{what}"
                );
            }
            let mut spec = holed.clone();
            spec.base.cutouts.push(Cutout {
                outline: vec![[38.0, 12.0], [50.0, 13.0], [48.0, 22.0], [39.0, 21.0]],
                bulges: vec![0.0, 0.3, 0.0, -0.2],
            });
            spec.base.punches.push(punch.clone());
            let f = crate::flat::flat_pattern(&spec).unwrap();
            assert_eq!((f.holes.len(), f.cutouts.len(), f.punches.len()), (3, 3, 2));
            let solid = crate::fold::fold(&spec).unwrap();
            assert!(crate::fold::mass(&solid, 1e-1).is_ok(), "{what}");
        }
        let mut cylinder = part([50.0, 50.0], 60.0, 180.0, Radial::Outside);
        let plain = crate::fold::exact_mass(&cylinder).unwrap();
        cylinder.base.holes.push(Hole {
            center: [30.0, 30.0],
            diameter: 8.0,
        });
        let m = crate::fold::exact_mass(&cylinder).unwrap();
        let taken = 16.0 * PI * 2.0 * 51.0 / 50.8;
        assert!((plain.volume - m.volume - taken).abs() < 1e-9);
    }

    /// The wrap takes a flat point, through to_wrap, where the sheet has it.
    #[test]
    fn the_wrap_rolls_the_flat_onto_the_sheet() {
        for (radii, sweep) in [
            ([50.0, 30.0], 270.0),
            ([30.0, 50.0], 120.0),
            ([50.0, 50.0], 200.0),
        ] {
            let spec = part(radii, 60.0, sweep, Radial::Outside);
            let s = Section::new(&spec, spec.base.cone.as_ref().unwrap()).unwrap();
            let wrap = s.wrap();
            for (u, v, o) in [(0.1, 3.0, 0.0), (3.5, 50.0, 2.0), (4.6, 20.0, 1.0)] {
                let u = u * sweep / 270.0_f64;
                let q = s.to_wrap(s.flat(u, v));
                let at = s.point(s.at(v, o), u);
                let got = wrap.point([q[0], q[1], s.level(0) + o]);
                let d = (0..3).map(|j| (got[j] - at[j]).powi(2)).sum::<f64>().sqrt();
                assert!(d < 1e-9, "{radii:?} {sweep}: {got:?} vs {at:?}");
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
        // A hole across the sheet's edge (a notch) isn't built yet.
        spec.base.holes.push(crate::spec::Hole {
            center: [10.0, 1.0],
            diameter: 5.0,
        });
        assert!(matches!(
            crate::fold::fold(&spec),
            Err(SheetError::Unsupported(_))
        ));
        spec.base.holes.clear();
        spec.base.cutouts.push(crate::spec::Cutout {
            outline: vec![[10.0, -2.0], [20.0, -2.0], [20.0, 5.0], [10.0, 5.0]],
            bulges: Vec::new(),
        });
        assert!(matches!(
            crate::flat::flat_pattern(&spec),
            Err(SheetError::Unsupported(_))
        ));
    }
}
