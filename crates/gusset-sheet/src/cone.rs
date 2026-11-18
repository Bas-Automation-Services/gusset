//! A part rolled to a cone or cylinder segment (`base.cone`): its flat
//! pattern, its exact solid and its mass, all in closed form.
//!
//! In the meridian half-plane (r, z) the sheet is a rectangle: the slant line
//! through the two arcs, `L` long, and the thickness square to it. Turned
//! around the z axis by the sweep, its long sides give the sheet's faces
//! (cones, or cylinders), its short sides the ends (cones square to the
//! sheet, or planes), and the two meridian planes close it.

use crate::flat::{Arc, Cell, Curve, FlatPattern, Layout, Zone, P2};
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

/// A piece of a cut's edge: its edge on each face, its vertices at either
/// end, and its flat points there.
type Piece3 = ([u32; 2], [u32; 2], [u32; 2], P2, P2);

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
        if !b.outline.is_empty() || b.profile.is_some() || b.loft.is_some() {
            return invalid("a base has an outline, a profile, a cone or a loft, one of them");
        }
        let [r0, r1] = cone.radii;
        if !(r0 > 0.0 && r1 > 0.0 && cone.height > 0.0) {
            return invalid("needs positive radii and height");
        }
        if !(cone.sweep > 0.0 && cone.sweep < 360.0) {
            return invalid("the sweep must lie between 0 and 360 degrees");
        }
        bare(spec, "cone")?;
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
        let f = self.pattern();
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
            holes: Vec::new(),
            embosses: Vec::new(),
        };
        crate::cut::cut_sheet(spec, &mut layout)?;
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
        let b = self.sweep;
        let mut builder = Builder::new("gusset".into(), accuracy())?;
        let r = Rolled::new(self);
        let (outer, mut walls) = self.rim(&r, &mut builder, flat)?;
        let (inner, cut_walls) = self.cuts(&r, &mut builder, flat)?;
        walls.extend(cut_walls);
        let theta_mid = b / 2.0;
        let toward = |n: P2, theta: f64| [n[0] * theta.cos(), n[0] * theta.sin(), n[1]];
        let outward = if self.o[1] > self.o[0] {
            [[-self.m[0], -self.m[1]], self.m]
        } else {
            [self.m, [-self.m[0], -self.m[1]]]
        };
        for oi in 0..2 {
            let q = self.at(0.0, self.o[oi]);
            let surface = self.surface(q, self.d);
            let inside = self.point(self.at(self.l / 2.0, self.o[oi]), theta_mid);
            let hint = surface.project(inside, None)?;
            // The face's parameter box: points over the whole turned sheet.
            let mut bounds = [[f64::MAX, f64::MIN]; 2];
            for i in 0..=32 {
                for j in 0..=2 {
                    let v = self.l * j as f64 / 2.0;
                    let point = self.point(self.at(v, self.o[oi]), b * i as f64 / 32.0);
                    let uv = surface.project(point, Some(hint))?;
                    for k in 0..2 {
                        bounds[k][0] = bounds[k][0].min(uv[k]);
                        bounds[k][1] = bounds[k][1].max(uv[k]);
                    }
                }
            }
            let bounds = bounds.map(|[lo, hi]| [lo - 1e-6, hi + 1e-6]);
            let along = |loop_: &[(u32, u32, u32, Orientation, P3)]| {
                loop_
                    .iter()
                    .map(|&(e, from, to, sense, at)| {
                        let near = surface.project(at, Some(hint))?;
                        projected_boundary(&builder, &surface, e, from, to, sense, Some(near))
                    })
                    .collect::<Result<Vec<Use>, GeometryError>>()
            };
            let uses = along(&outer[oi])?;
            // The holes and cut-outs, clockwise in the flat.
            let holes = inner[oi]
                .iter()
                .map(|ring| along(ring))
                .collect::<Result<Vec<_>, _>>()?;
            let normal = surface.normal_at(hint)?;
            let key = format!("sheet.{oi}");
            let face = builder.face_with_holes(&key, surface, bounds, uses, holes)?;
            if dot(normal, toward(outward[oi], theta_mid)) < 0.0 {
                builder.set_face_sense(face, Orientation::Reverse)?;
            }
        }
        for wall in walls {
            let normal = wall.surface.normal_at(wall.mid)?;
            let face = builder.face(&wall.key, wall.surface, wall.bounds, wall.uses)?;
            if dot(normal, wall.out) < 0.0 {
                builder.set_face_sense(face, Orientation::Reverse)?;
            }
        }
        builder.finish_solid()
    }

    /// The flat outline's edges on both faces of the sheet, in its order
    /// (per face: edge, from, to, sense, the point where the edge starts), and its walls:
    /// the ends (cones square to the sheet), the sides (meridian planes) and
    /// the walls of notches cut into them.
    #[allow(clippy::type_complexity)]
    fn rim(
        &self,
        r: &Rolled,
        builder: &mut Builder,
        flat: &FlatPattern,
    ) -> Result<([Vec<(u32, u32, u32, Orientation, P3)>; 2], Vec<Wall>), GeometryError> {
        let (b, l) = (self.sweep, self.l);
        let points = &flat.outline;
        let n = points.len();
        let corners: Vec<[u32; 2]> = points
            .iter()
            .map(|&p| [0, 1].map(|oi| builder.vertex(r.at(p, oi))))
            .collect();
        let mut verticals = Vec::with_capacity(n);
        for &p in points {
            verticals.push(r.vertical(builder, p, false)?);
        }
        // Up through the sheet from face 0 to face 1.
        let up = if r.lo == 0 {
            Orientation::Forward
        } else {
            Orientation::Reverse
        };
        let down = if r.lo == 0 {
            Orientation::Reverse
        } else {
            Orientation::Forward
        };
        let same = |a: f64, c: f64| (a - c).abs() < 1e-7;
        let mut walls = Vec::new();
        let mut keys = std::collections::HashMap::<String, usize>::new();
        let mut key = |base: String| {
            let count = keys.entry(base.clone()).or_insert(0);
            *count += 1;
            if *count == 1 {
                base
            } else {
                format!("{base}#{}", *count - 1)
            }
        };
        let mut segments: Vec<(Vec<Piece3>, bool)> = Vec::with_capacity(n);
        for k in 0..n {
            let j = (k + 1) % n;
            let (p, q) = (points[k], points[j]);
            let ((tp, vp), (tq, vq)) = (r.place(p), r.place(q));
            // Its middle too: a notch's arc may end on an end or a side.
            let middle = match flat.outline_arcs[k].and_then(Curve::arc) {
                None => [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0],
                Some(a) => {
                    let c = a.center;
                    let turn = (p[1] - c[1]).atan2(p[0] - c[0]) + a.sweep(p, q) / 2.0;
                    let radius = (p[0] - c[0]).hypot(p[1] - c[1]);
                    [c[0] + radius * turn.cos(), c[1] + radius * turn.sin()]
                }
            };
            let (tm, vm) = r.place(middle);
            let end = [0.0, l]
                .iter()
                .position(|&v| same(vp, v) && same(vq, v) && same(vm, v));
            let side = [0.0, b]
                .iter()
                .position(|&t| same(tp, t) && same(tq, t) && same(tm, t));
            if let Some(vi) = end {
                // An end: arcs round the axis, run with the angle.
                let forward = tp < tq;
                let (t0, t1) = (tp.min(tq), tp.max(tq));
                let (s, e) = if forward { (k, j) } else { (j, k) };
                let span = Interval::new(t0, t1)?;
                let edges = self.o.map(|o| {
                    let c = self.at(vi as f64 * l, o);
                    builder.edge(
                        CurveGeometry::Circle {
                            frame: Frame3 {
                                origin: [0.0, 0.0, c[1]],
                                x: [1.0, 0.0, 0.0],
                                y: [0.0, 1.0, 0.0],
                                z: [0.0, 0.0, 1.0],
                            },
                            radius: c[0],
                        },
                        span,
                        false,
                    )
                });
                let q0 = self.at(vi as f64 * l, 0.0);
                let surface = self.surface(q0, self.m);
                let through = |f: f64| self.o[0] + f * (self.o[1] - self.o[0]);
                let mid = (t0 + t1) / 2.0;
                let samples: Vec<P3> = (0..=32)
                    .flat_map(|i| {
                        (0..=2).map(move |jj| (t0 + (t1 - t0) * i as f64 / 32.0, jj as f64 / 2.0))
                    })
                    .map(|(t, f)| self.point(self.at(vi as f64 * l, through(f)), t))
                    .collect();
                let inside = self.point(self.at(vi as f64 * l, through(0.5)), mid);
                let out = if vi == 0 {
                    [-self.d[0], -self.d[1]]
                } else {
                    self.d
                };
                let loop_ = [
                    (edges[0], corners[s][0], corners[e][0], Orientation::Forward),
                    (verticals[e], corners[e][0], corners[e][1], up),
                    (edges[1], corners[e][1], corners[s][1], Orientation::Reverse),
                    (verticals[s], corners[s][1], corners[s][0], down),
                ];
                walls.push(self.wall(
                    builder,
                    key(format!("end.{vi}")),
                    surface,
                    &loop_,
                    inside,
                    [out[0] * mid.cos(), out[0] * mid.sin(), out[1]],
                    &samples,
                )?);
                let piece = (edges, corners[s], corners[e], points[s], points[e]);
                segments.push((vec![piece], forward));
            } else if let Some(th) = side {
                // A side: straight up the slant, from its lower end.
                let forward = vp < vq;
                let (s, e) = if forward { (k, j) } else { (j, k) };
                let edges = [
                    line(builder, r.at(points[s], 0), r.at(points[e], 0))?,
                    line(builder, r.at(points[s], 1), r.at(points[e], 1))?,
                ];
                let theta = th as f64 * b;
                let normal = if th == 0 {
                    [0.0, -1.0, 0.0]
                } else {
                    [-theta.sin(), theta.cos(), 0.0]
                };
                let radial = [theta.cos(), theta.sin(), 0.0];
                let surface = SurfaceGeometry::Plane {
                    frame: Frame3 {
                        origin: [0.0; 3],
                        x: radial,
                        y: cross(normal, radial),
                        z: normal,
                    },
                };
                let samples: Vec<P3> = [s, e]
                    .iter()
                    .flat_map(|&c| [0, 1].map(|oi| r.at(points[c], oi)))
                    .collect();
                let (va, vb) = (vp.min(vq), vp.max(vq));
                let inside = self.point(
                    self.at((va + vb) / 2.0, (self.o[0] + self.o[1]) / 2.0),
                    theta,
                );
                let loop_ = [
                    (edges[0], corners[s][0], corners[e][0], Orientation::Forward),
                    (verticals[e], corners[e][0], corners[e][1], up),
                    (edges[1], corners[e][1], corners[s][1], Orientation::Reverse),
                    (verticals[s], corners[s][1], corners[s][0], down),
                ];
                walls.push(self.wall(
                    builder,
                    key(format!("side.{th}")),
                    surface,
                    &loop_,
                    inside,
                    normal,
                    &samples,
                )?);
                let piece = (edges, corners[s], corners[e], points[s], points[e]);
                segments.push((vec![piece], forward));
            } else {
                // A notch's edge: the material lies on its left.
                let arc = flat.outline_arcs[k].and_then(Curve::arc);
                let (pieces, forward, wall) = r.segment(
                    builder,
                    [p, q],
                    arc,
                    [corners[k], corners[j]],
                    [verticals[k], verticals[j]],
                    key(format!("notch.{k}")),
                    -1.0,
                )?;
                walls.push(wall);
                segments.push((pieces, forward));
            }
        }
        let mut outer: [Vec<_>; 2] = [Vec::new(), Vec::new()];
        for (oi, uses) in outer.iter_mut().enumerate() {
            for (list, forward) in &segments {
                if *forward {
                    uses.extend(list.iter().map(|&(e, a, c, x, _)| {
                        (e[oi], a[oi], c[oi], Orientation::Forward, r.at(x, oi))
                    }));
                } else {
                    uses.extend(list.iter().rev().map(|&(e, a, c, x, _)| {
                        (e[oi], c[oi], a[oi], Orientation::Reverse, r.at(x, oi))
                    }));
                }
            }
        }
        Ok((outer, walls))
    }

    /// An end or side wall: its loop on an analytic surface, its parameter
    /// box from points over it, facing `out` at `inside`.
    #[allow(clippy::too_many_arguments)]
    fn wall(
        &self,
        builder: &Builder,
        key: String,
        surface: SurfaceGeometry,
        loop_: &[(u32, u32, u32, Orientation)],
        inside: P3,
        out: P3,
        samples: &[P3],
    ) -> Result<Wall, GeometryError> {
        let hint = surface.project(inside, None)?;
        let mut uses = Vec::with_capacity(loop_.len());
        for &(e, from, to, sense) in loop_ {
            uses.push(match &surface {
                SurfaceGeometry::Plane { frame } => {
                    plane_boundary(builder, *frame, e, from, to, sense)?
                }
                _ => projected_boundary(builder, &surface, e, from, to, sense, Some(hint))?,
            });
        }
        let mut bounds = [[f64::MAX, f64::MIN]; 2];
        for &point in samples {
            let uv = surface.project(point, Some(hint))?;
            for k in 0..2 {
                bounds[k][0] = bounds[k][0].min(uv[k]);
                bounds[k][1] = bounds[k][1].max(uv[k]);
            }
        }
        Ok(Wall {
            key,
            surface,
            bounds: bounds.map(|[lo, hi]| [lo - 1e-6, hi + 1e-6]),
            uses,
            mid: hint,
            out,
        })
    }

    /// The edges of the holes and cut-outs on both faces of the sheet (per
    /// face, the inner loops as (edge, from, to, sense, the edge's start)),
    /// and their walls.
    #[allow(clippy::type_complexity)]
    fn cuts(
        &self,
        r: &Rolled,
        builder: &mut Builder,
        flat: &FlatPattern,
    ) -> Result<([Vec<Vec<(u32, u32, u32, Orientation, P3)>>; 2], Vec<Wall>), GeometryError> {
        let (t, lo, hi, z) = (self.t, r.lo, r.hi, r.z);
        let mut inner: [Vec<Vec<_>>; 2] = [Vec::new(), Vec::new()];
        let mut walls = Vec::new();

        for (k, h) in flat.holes.iter().enumerate() {
            let rad = h.diameter / 2.0;
            let on = |a: f64| [h.center[0] + rad * a.cos(), h.center[1] + rad * a.sin()];
            let m = r.pieces(&|f| on(TAU * f));
            let angle = |j: usize| TAU * j as f64 / m as f64;
            let v: Vec<[u32; 2]> = (0..m)
                .map(|j| [0, 1].map(|oi| builder.vertex(r.at(on(angle(j)), oi))))
                .collect();
            let mut arcs = Vec::with_capacity(m);
            for j in 0..m {
                let mut e = [0; 2];
                for oi in 0..2 {
                    e[oi] = builder.edge(
                        r.wrapped(CurveGeometry::Circle {
                            frame: r.flat_frame(h.center, z[oi]),
                            radius: rad,
                        }),
                        Interval::new(angle(j), angle(j + 1))?,
                        false,
                    );
                }
                arcs.push(e);
            }
            let next = |j: usize| (j + 1) % m;
            for oi in 0..2 {
                inner[oi].push(
                    (0..m)
                        .rev()
                        .map(|j| {
                            let start = r.at(on(angle(j)), oi);
                            (
                                arcs[j][oi],
                                v[next(j)][oi],
                                v[j][oi],
                                Orientation::Reverse,
                                start,
                            )
                        })
                        .collect(),
                );
            }
            let seam = r.vertical(builder, on(0.0), true)?;
            let mut uses: Vec<Use> = (0..m)
                .map(|j| {
                    let line = uv_line([0.0, 0.0], [1.0, 0.0]);
                    boundary(
                        arcs[j][lo],
                        v[j][lo],
                        v[next(j)][lo],
                        Orientation::Forward,
                        line,
                    )
                })
                .collect();
            uses.push(boundary(
                seam,
                v[0][lo],
                v[0][hi],
                Orientation::Forward,
                uv_line([TAU, 0.0], [0.0, 1.0]),
            ));
            uses.extend((0..m).rev().map(|j| {
                let line = uv_line([0.0, t], [1.0, 0.0]);
                boundary(
                    arcs[j][hi],
                    v[next(j)][hi],
                    v[j][hi],
                    Orientation::Reverse,
                    line,
                )
            }));
            uses.push(boundary(
                seam,
                v[0][hi],
                v[0][lo],
                Orientation::Reverse,
                uv_line([0.0, 0.0], [0.0, 1.0]),
            ));
            walls.push(Wall {
                key: format!("hole.{k}"),
                surface: SurfaceGeometry::Wrapped {
                    wrap: r.wrap,
                    surface: Box::new(SurfaceGeometry::Cylinder {
                        frame: r.flat_frame(h.center, z[lo]),
                        radius: rad,
                    }),
                },
                bounds: [[0.0, TAU], [0.0, t]],
                uses,
                mid: [PI, t / 2.0],
                out: r.into([h.center[0] - rad, h.center[1]], [1.0, 0.0]),
            });
        }

        for (k, ring) in flat.cutouts.iter().enumerate() {
            let n = ring.points.len();
            let corners: Vec<[u32; 2]> = ring
                .points
                .iter()
                .map(|&p| [0, 1].map(|oi| builder.vertex(r.at(p, oi))))
                .collect();
            let mut verticals = Vec::with_capacity(n);
            for &p in &ring.points {
                verticals.push(r.vertical(builder, p, false)?);
            }
            let mut segments = Vec::with_capacity(n);
            for i in 0..n {
                let j = (i + 1) % n;
                let (pieces, forward, wall) = r.segment(
                    builder,
                    [ring.points[i], ring.points[j]],
                    ring.arcs[i].and_then(Curve::arc),
                    [corners[i], corners[j]],
                    [verticals[i], verticals[j]],
                    format!("cut.{k}.{i}"),
                    1.0,
                )?;
                walls.push(wall);
                segments.push((pieces, forward));
            }
            // Each face's loop runs the ring backwards, q -> p.
            for oi in 0..2 {
                let mut ring_uses = Vec::new();
                for (list, forward) in segments.iter().rev() {
                    if *forward {
                        ring_uses.extend(list.iter().rev().map(|&(e, a, b, x, _)| {
                            (e[oi], b[oi], a[oi], Orientation::Reverse, r.at(x, oi))
                        }));
                    } else {
                        ring_uses.extend(list.iter().map(|&(e, a, b, x, _)| {
                            (e[oi], a[oi], b[oi], Orientation::Forward, r.at(x, oi))
                        }));
                    }
                }
                inner[oi].push(ring_uses);
            }
        }
        Ok((inner, walls))
    }

    /// The mass of the sheet less its holes and cut-outs, with their walls.
    fn mass(&self, flat: &FlatPattern) -> Mass {
        // Closed form, or with notches in the edge, over the outline.
        let whole = if flat.outline == self.pattern().outline {
            self.plain_mass()
        } else {
            let (sheet, wall) = self.removed(&pieces(&flat.outline, &flat.outline_arcs));
            Mass {
                volume: sheet[0],
                area: sheet[1] + sheet[2] + wall,
                center_of_mass: [3, 4, 5].map(|j| sheet[j] / sheet[0]),
            }
        };
        let mut volume = whole.volume;
        let mut area = whole.area;
        let mut moment = whole.center_of_mass.map(|c| c * whole.volume);
        let holes = flat
            .holes
            .iter()
            .map(|h| vec![Piece::Arc(h.center, h.diameter / 2.0, 0.0, TAU)]);
        let rings = flat
            .cutouts
            .iter()
            .map(|ring| pieces(&ring.points, &ring.arcs));
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

/// A face made later: its loop is built, its sense set from `out`, the way
/// out of the material at `mid`.
struct Wall {
    key: String,
    surface: SurfaceGeometry,
    bounds: [[f64; 2]; 2],
    uses: Vec<Use>,
    mid: [f64; 2],
    out: P3,
}

/// The sheet as the wrap rolls it: flat points onto either face, and the
/// edges and walls of cuts through it.
struct Rolled<'a> {
    s: &'a Section,
    wrap: Wrap,
    z: [f64; 2],
    lo: usize,
    hi: usize,
}

impl Rolled<'_> {
    fn new(s: &Section) -> Rolled<'_> {
        let z = [s.level(0), s.level(1)];
        let (lo, hi) = if z[0] < z[1] { (0, 1) } else { (1, 0) };
        Rolled {
            s,
            wrap: s.wrap(),
            z,
            lo,
            hi,
        }
    }

    /// Flat point p on face oi.
    fn at(&self, p: P2, oi: usize) -> P3 {
        let q = self.s.to_wrap(p);
        self.wrap.point([q[0], q[1], self.z[oi]])
    }

    /// A flat point's turn angle (from 0 to the sweep) and slant position.
    fn place(&self, p: P2) -> (f64, f64) {
        let q = self.s.to_wrap(p);
        let n = self.wrap.point([q[0], q[1], self.wrap.radius]);
        let b = self.s.sweep;
        let mut theta = n[1].atan2(n[0]).rem_euclid(TAU);
        if theta > (b + TAU) / 2.0 {
            theta -= TAU;
        }
        let n0 = self.s.at(0.0, self.s.neutral);
        let v = (n[0].hypot(n[1]) - n0[0]) * self.s.d[0] + (n[2] - n0[1]) * self.s.d[1];
        (theta, v)
    }

    /// The way into the cut from flat point p along flat direction d, mid-sheet.
    fn into(&self, p: P2, d: P2) -> P3 {
        let mid = (self.z[0] + self.z[1]) / 2.0;
        let h = 1e-6;
        let (q, e) = (self.s.to_wrap(p), self.s.to_wrap_dir(d));
        sub(
            self.wrap.point([q[0] + h * e[0], q[1] + h * e[1], mid]),
            self.wrap.point([q[0], q[1], mid]),
        )
    }

    fn flat_frame(&self, c: P2, z: f64) -> Frame3 {
        let (o, x, y) = (
            self.s.to_wrap(c),
            self.s.to_wrap_dir([1.0, 0.0]),
            self.s.to_wrap_dir([0.0, 1.0]),
        );
        Frame3 {
            origin: [o[0], o[1], z],
            x: [x[0], x[1], 0.0],
            y: [y[0], y[1], 0.0],
            z: [0.0, 0.0, 1.0],
        }
    }

    fn wrapped(&self, curve: CurveGeometry) -> CurveGeometry {
        CurveGeometry::Wrapped {
            wrap: self.wrap,
            curve: Box::new(curve),
        }
    }

    /// The straight edge through the sheet at flat point p.
    fn vertical(&self, builder: &mut Builder, p: P2, seam: bool) -> Result<u32, GeometryError> {
        let (a, c) = (self.at(p, self.lo), self.at(p, self.hi));
        Ok(builder.edge(
            CurveGeometry::Line {
                origin: a,
                direction: sub(c, a).map(|x| x / norm(sub(c, a))),
            },
            Interval::new(0.0, self.s.t)?,
            seam,
        ))
    }

    /// Edges on the sheet's faces take projected pcurves, which follow an
    /// edge only while it strays less than a quarter turn from a straight
    /// guess between its ends: how many pieces of at most a fifth of a half
    /// turn around the axis the curve through `point` (on 0..1) needs.
    fn pieces(&self, point: &dyn Fn(f64) -> P2) -> usize {
        let turns: Vec<f64> = (0..=16)
            .map(|i| self.s.developed(point(i as f64 / 16.0)).0)
            .collect();
        let (lo, hi) = turns
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), &t| (lo.min(t), hi.max(t)));
        (((hi - lo) / (0.2 * PI)).ceil() as usize).max(1)
    }

    /// A cut's straight or arc segment from flat p to q, between vertices
    /// `ends` and through-edges `verticals` there: its pieces in the edge's
    /// direction, whether that runs p -> q, and its wall, the cut on the
    /// `side` of p -> q (1 left, -1 right).
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn segment(
        &self,
        builder: &mut Builder,
        [p, q]: [P2; 2],
        arc: Option<Arc>,
        ends: [[u32; 2]; 2],
        verticals: [u32; 2],
        key: String,
        side: f64,
    ) -> Result<(Vec<Piece3>, bool, Wall), GeometryError> {
        let (t, z, lo, hi) = (self.s.t, self.z, self.lo, self.hi);
        // The edge's curve on each face, its parameter range, its flat point
        // at a parameter, which of p and q it starts at, its wall's surface
        // and the way into the cut at its middle.
        let (curves, range, point, starts_at_p, surface, into, forward): (
            [CurveGeometry; 2],
            [f64; 2],
            Box<dyn Fn(f64) -> P2>,
            bool,
            SurfaceGeometry,
            P3,
            bool,
        ) = match arc {
            None => {
                let d = [q[0] - p[0], q[1] - p[1]];
                let len = d[0].hypot(d[1]);
                let u = [d[0] / len, d[1] / len];
                let (o, uw) = (self.s.to_wrap(p), self.s.to_wrap_dir(u));
                (
                    [0, 1].map(|oi| {
                        self.wrapped(CurveGeometry::Line {
                            origin: [o[0], o[1], z[oi]],
                            direction: [uw[0], uw[1], 0.0],
                        })
                    }),
                    [0.0, len],
                    Box::new(move |s| [p[0] + s * u[0], p[1] + s * u[1]]),
                    true,
                    SurfaceGeometry::Wrapped {
                        wrap: self.wrap,
                        surface: Box::new(SurfaceGeometry::Plane {
                            frame: Frame3 {
                                origin: [o[0], o[1], z[lo]],
                                x: [uw[0], uw[1], 0.0],
                                y: [0.0, 0.0, 1.0],
                                z: [uw[1], -uw[0], 0.0],
                            },
                        }),
                    },
                    self.into(
                        [p[0] + d[0] / 2.0, p[1] + d[1] / 2.0],
                        [-u[1] * side, u[0] * side],
                    ),
                    true,
                )
            }
            Some(arc) => {
                let c = arc.center;
                let r = (p[0] - c[0]).hypot(p[1] - c[1]);
                let a_p = (p[1] - c[1]).atan2(p[0] - c[0]);
                let sweep = arc.sweep(p, q);
                // The edge runs with the angle.
                let (a0, a1) = if sweep > 0.0 {
                    (a_p, a_p + sweep)
                } else {
                    (a_p + sweep, a_p)
                };
                let (sm, cm) = ((a0 + a1) / 2.0).sin_cos();
                // Left of p -> q there.
                let left = if sweep > 0.0 { [-cm, -sm] } else { [cm, sm] };
                (
                    [0, 1].map(|oi| {
                        self.wrapped(CurveGeometry::Circle {
                            frame: self.flat_frame(c, z[oi]),
                            radius: r,
                        })
                    }),
                    [a0, a1],
                    Box::new(move |a: f64| [c[0] + r * a.cos(), c[1] + r * a.sin()]),
                    sweep > 0.0,
                    SurfaceGeometry::Wrapped {
                        wrap: self.wrap,
                        surface: Box::new(SurfaceGeometry::Cylinder {
                            frame: self.flat_frame(c, z[lo]),
                            radius: r,
                        }),
                    },
                    self.into(
                        [c[0] + r * cm, c[1] + r * sm],
                        [left[0] * side, left[1] * side],
                    ),
                    sweep > 0.0,
                )
            }
        };
        let (start, end, start_vertical, end_vertical) = if starts_at_p {
            (ends[0], ends[1], verticals[0], verticals[1])
        } else {
            (ends[1], ends[0], verticals[1], verticals[0])
        };
        let [s0, s1] = range;
        let param = |f: f64| s0 + f * (s1 - s0);
        let m = self.pieces(&|f| point(param(f)));
        let mut v = vec![start];
        for piece in 1..m {
            let x = point(param(piece as f64 / m as f64));
            v.push([0, 1].map(|oi| builder.vertex(self.at(x, oi))));
        }
        v.push(end);
        let mut list = Vec::with_capacity(m);
        for piece in 0..m {
            let span = Interval::new(
                param(piece as f64 / m as f64),
                param((piece + 1) as f64 / m as f64),
            )?;
            let edges = [0, 1].map(|oi| builder.edge(curves[oi].clone(), span, false));
            let flat = [span.lo, span.hi].map(&point);
            list.push((edges, v[piece], v[piece + 1], flat[0], flat[1]));
        }
        let mut uses: Vec<Use> = list
            .iter()
            .map(|&(e, a, b, _, _)| {
                let line = uv_line([0.0, 0.0], [1.0, 0.0]);
                boundary(e[lo], a[lo], b[lo], Orientation::Forward, line)
            })
            .collect();
        uses.push(boundary(
            end_vertical,
            end[lo],
            end[hi],
            Orientation::Forward,
            uv_line([s1, 0.0], [0.0, 1.0]),
        ));
        uses.extend(list.iter().rev().map(|&(e, a, b, _, _)| {
            let line = uv_line([0.0, t], [1.0, 0.0]);
            boundary(e[hi], b[hi], a[hi], Orientation::Reverse, line)
        }));
        uses.push(boundary(
            start_vertical,
            start[hi],
            start[lo],
            Orientation::Reverse,
            uv_line([s0, 0.0], [0.0, 1.0]),
        ));
        let wall = Wall {
            key,
            surface,
            bounds: [[s0, s1], [0.0, t]],
            uses,
            mid: [(s0 + s1) / 2.0, t / 2.0],
            out: into,
        };
        Ok((list, forward, wall))
    }
}

/// A closed path's segments for integration.
fn pieces(points: &[P2], arcs: &[Option<Curve>]) -> Vec<Piece> {
    let n = points.len();
    (0..n)
        .map(|k| {
            let (p, q) = (points[k], points[(k + 1) % n]);
            match arcs[k].and_then(Curve::arc) {
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
        .collect()
}

fn line(builder: &mut Builder, p: P3, q: P3) -> Result<u32, GeometryError> {
    let length = norm(sub(q, p));
    Ok(builder.edge(
        CurveGeometry::Line {
            origin: p,
            direction: sub(q, p).map(|x| x / length),
        },
        Interval::new(0.0, length)?,
        false,
    ))
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
    spec.base
        .cone
        .as_ref()
        .filter(|c| c.facets.is_none())
        .map(|c| Section::new(spec, c))
}

/// Why a cone can't be built, if it can't.
/// A cone or loft carries holes, cut-outs and punches in its flat pattern,
/// no other features.
pub(crate) fn bare(spec: &PartSpec, what: &str) -> Result<(), SheetError> {
    let b = &spec.base;
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
    match features.iter().find(|f| f.0) {
        Some((_, feature)) => Err(SheetError::Unsupported(format!("{feature} on a {what}"))),
        None => Ok(()),
    }
}

pub(crate) fn check(spec: &PartSpec, cone: &Cone) -> Result<(), SheetError> {
    Section::new(spec, cone).map(|_| ())
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
            facets: None,
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
    fn notches_in_a_cones_edges_take_what_they_cut() {
        use crate::spec::{Cutout, Hole};
        for (radii, material) in [
            ([50.0, 30.0], Radial::Outside),
            ([30.0, 50.0], Radial::Inside),
            ([50.0, 50.0], Radial::Outside),
        ] {
            let what = format!("{radii:?} {material:?}");
            let mut spec = part(radii, 60.0, 180.0, material);
            let plain = crate::fold::exact_mass(&spec).unwrap();
            // A hole across the first side, a slot across the first arc.
            spec.base.holes.push(Hole {
                center: [3.0, 40.0],
                diameter: 10.0,
            });
            spec.base.cutouts.push(Cutout {
                outline: vec![[15.0, -6.0], [25.0, -6.0], [25.0, 4.0], [15.0, 4.0]],
                bulges: vec![0.0, 0.0, 0.5, 0.0],
            });
            let f = crate::flat::flat_pattern(&spec).unwrap();
            assert!(f.holes.is_empty() && f.cutouts.is_empty(), "{what}");
            assert!(f.outline.len() > 4, "{what}");
            let exact = crate::fold::exact_mass(&spec).unwrap();
            assert!(exact.volume < plain.volume - 100.0, "{what}");
            let mesh = crate::fold::mass(&crate::fold::fold(&spec).unwrap(), 1e-3).unwrap();
            assert!(
                ((mesh.volume - exact.volume) / exact.volume).abs() < 2e-5,
                "{what}: mesh {} exact {}",
                mesh.volume,
                exact.volume
            );
            assert!(
                ((mesh.area - exact.area) / exact.area).abs() < 2e-5,
                "{what}: mesh {} exact {}",
                mesh.area,
                exact.area
            );
        }
    }

    /// The flat outline is the closed form's sheet: integrated over it, the
    /// sheet measures as the closed form does.
    #[test]
    fn the_sheet_integrated_over_its_outline_is_the_closed_form() {
        for (radii, material) in [
            ([50.0, 30.0], Radial::Outside),
            ([30.0, 50.0], Radial::Inside),
            ([50.0, 50.0], Radial::Outside),
        ] {
            let spec = part(radii, 60.0, 270.0, material);
            let Some(Ok(s)) = section(&spec) else {
                panic!()
            };
            let f = s.pattern();
            let (sheet, wall) = s.removed(&pieces(&f.outline, &f.outline_arcs));
            let plain = s.plain_mass();
            assert!((sheet[0] - plain.volume).abs() < 1e-9 * plain.volume);
            assert!((sheet[1] + sheet[2] + wall - plain.area).abs() < 1e-9 * plain.area);
            for j in 0..3 {
                let c = sheet[3 + j] / sheet[0];
                assert!((c - plain.center_of_mass[j]).abs() < 1e-8, "{radii:?} {j}");
            }
        }
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
        // A cut across the whole sheet would split it.
        spec.base.cutouts.push(crate::spec::Cutout {
            outline: vec![[-5.0, 20.0], [200.0, 20.0], [200.0, 25.0], [-5.0, 25.0]],
            bulges: Vec::new(),
        });
        assert!(matches!(
            crate::flat::flat_pattern(&spec),
            Err(SheetError::Unsupported(_))
        ));
    }
}
