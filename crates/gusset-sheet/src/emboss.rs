//! Formed embosses (`base.embosses`): the sheet pressed out of the base plate
//! at its thickness. Seen in a section square to its edge, each face of the
//! sheet runs from the plate over a rounded corner, along the wall, over a
//! second rounded corner into the raised face. That section, (s, z) with s
//! out from the raised face's mold line, is the same all round: a round
//! emboss turns it around its axis (tori, a cone or a cylinder, a disk); a
//! rectangular one runs it along each side (cylinders and a plane) and turns
//! it around each corner's axis (tori and a cone). All exact. The plate keeps
//! a hole the shape of the footprint, where the sheet leaves it.

use crate::flat::{Arc, Curve, Ring, P2};
use crate::fold::{Mass, GAUSS, P3};
use crate::spec::{Emboss, EmbossShape, Hole, Side};
use crate::SheetError;
use opengeometry::brep::{
    boundary, plane_boundary, projected_boundary, uv_line, Builder, CurveGeometry, Frame3,
    GeometryError, Orientation, Surface, SurfaceGeometry,
};
use opengeometry::math::Interval;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

const EPS: f64 = 1e-9;
const Z: P3 = [0.0, 0.0, 1.0];

/// A piece of a face's section, (s, z).
#[derive(Clone, Copy, Debug)]
enum Meridian {
    /// Around `c` with radius `rho`, from angle `a0` to `a1` (from +s towards
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
    /// The point a fraction f along, and its derivative.
    fn at(&self, f: f64) -> (P2, P2) {
        match *self {
            Meridian::Arc { c, rho, a0, a1 } => {
                let a = a0 + f * (a1 - a0);
                let (sn, cs) = a.sin_cos();
                (
                    [c[0] + rho * cs, c[1] + rho * sn],
                    [-rho * sn * (a1 - a0), rho * cs * (a1 - a0)],
                )
            }
            Meridian::Line { a, b } => (
                [a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])],
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

    /// Moved by ds along s.
    fn shifted(self, ds: f64) -> Meridian {
        match self {
            Meridian::Arc { c, rho, a0, a1 } => Meridian::Arc {
                c: [c[0] + ds, c[1]],
                rho,
                a0,
                a1,
            },
            Meridian::Line { a, b } => Meridian::Line {
                a: [a[0] + ds, a[1]],
                b: [b[0] + ds, b[1]],
            },
        }
    }
}

/// The raised face's outline at its mold line.
#[derive(Clone, Copy, Debug)]
enum Plan {
    /// A circle of this radius.
    Round(f64),
    /// A rectangle 2 hw by 2 hl with corners of radius rc, turned by `angle`
    /// (radians).
    Rect {
        hw: f64,
        hl: f64,
        rc: f64,
        angle: f64,
    },
}

/// One of a rectangle's eight pieces of outline, counter-clockwise from the
/// right side's lower end: a side (even) or a corner (odd).
#[derive(Clone, Copy)]
struct Segment {
    /// Where it starts at the mold line, and the outward direction there.
    start: P2,
    n: P2,
    /// A side's direction; a corner's centre and start angle.
    d: P2,
    corner: Option<(P2, f64)>,
}

/// An emboss, worked out.
#[derive(Clone, Debug)]
pub(crate) struct Formed {
    pub id: String,
    pub center: P2,
    plan: Plan,
    /// Per face of the sheet, by the plate face it leaves (0: z = 0, 1:
    /// z = T): its section from the footprint inward, and the height of its
    /// raised face.
    chains: [[Meridian; 3]; 2],
    tops: [f64; 2],
    /// s where the sheet leaves the plate, and where the raised face starts.
    foot: f64,
    top: f64,
}

impl Formed {
    pub fn new(e: &Emboss, t: f64, bend_radius: f64) -> Result<Formed, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("emboss {}: {what}", e.id)));
        let (plan, rc) = match e.shape {
            EmbossShape::Round { diameter } => (Plan::Round(diameter / 2.0), diameter / 2.0),
            EmbossShape::Rectangle {
                width,
                length,
                corner_radius,
            } => {
                if !(corner_radius > 0.0
                    && width > 2.0 * corner_radius
                    && length > 2.0 * corner_radius)
                {
                    return invalid("needs corners of a positive radius, less than half its sides");
                }
                let plan = Plan::Rect {
                    hw: width / 2.0,
                    hl: length / 2.0,
                    rc: corner_radius,
                    angle: e.angle.to_radians(),
                };
                (plan, corner_radius)
            }
        };
        let r = e.radius.unwrap_or(bend_radius);
        let (h, d) = (e.height, e.draft.to_radians());
        if !(rc > 0.0 && h > 0.0 && r > 0.0) {
            return invalid("needs a positive size, height and radius");
        }
        if !(d >= 0.0 && d < FRAC_PI_2) {
            return invalid("the draft must lie between 0 and 90 degrees");
        }
        // Tangent length of each corner: they turn by 90 degrees less the
        // draft. The face on the raised side: its first corner has the inside
        // radius r, the second r + t.
        let k = ((FRAC_PI_2 - d) / 2.0).tan();
        let foot = h * d.tan() + r * k;
        let top = -(r + t) * k;
        if rc + top <= r + t + EPS {
            return invalid("too small across (or at its corners) for its radius and thickness");
        }
        let c1 = [foot, t + r];
        let c2 = [top, t + h - (r + t)];
        let nw = [d.cos(), d.sin()];
        let off = |c: P2, s: f64| [c[0] + s * nw[0], c[1] + s * nw[1]];
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
            ],
            [
                corner(c1, r, 1.5 * PI, PI + d),
                Meridian::Line { a, b },
                corner(c2, r + t, d, FRAC_PI_2),
            ],
        ];
        let (chains, tops) = match e.side {
            Side::Up => (up, [h, t + h]),
            // Mirrored, the face that left z = T leaves z = 0.
            Side::Down => (
                [up[1].map(|m| m.mirrored(t)), up[0].map(|m| m.mirrored(t))],
                [-h, t - h],
            ),
        };
        Ok(Formed {
            id: e.id.clone(),
            center: e.at,
            plan,
            chains,
            tops,
            foot,
            top,
        })
    }

    /// A round emboss's footprint, as the hole the plate keeps.
    pub fn hole(&self) -> Option<Hole> {
        match self.plan {
            Plan::Round(half) => Some(Hole {
                center: self.center,
                diameter: 2.0 * (half + self.foot),
            }),
            Plan::Rect { .. } => None,
        }
    }

    /// The footprint in the flat pattern, counter-clockwise; for a
    /// rectangle its points are the eight segments' starts.
    pub fn ring(&self) -> Ring {
        match self.plan {
            Plan::Round(half) => {
                let (c, r) = (self.center, half + self.foot);
                let arc = Some(Curve::Arc(Arc {
                    center: c,
                    ccw: true,
                }));
                Ring {
                    points: vec![[c[0] + r, c[1]], [c[0] - r, c[1]]],
                    arcs: vec![arc, arc],
                }
            }
            Plan::Rect { .. } => Ring {
                points: (0..8).map(|j| self.flat_at(j, 0.0, self.foot)).collect(),
                arcs: (0..8)
                    .map(|j| {
                        self.segment(j).corner.map(|(c, _)| {
                            Curve::Arc(Arc {
                                center: self.place(c),
                                ccw: true,
                            })
                        })
                    })
                    .collect(),
            },
        }
    }

    /// The rectangle's corner radius, or the circle's radius.
    fn corner_radius(&self) -> f64 {
        match self.plan {
            Plan::Round(half) => half,
            Plan::Rect { rc, .. } => rc,
        }
    }

    /// The outline's length, offset by s.
    fn perimeter(&self, s: f64) -> f64 {
        match self.plan {
            Plan::Round(half) => TAU * (half + s),
            Plan::Rect { hw, hl, rc, .. } => 4.0 * (hw + hl) - 8.0 * rc + TAU * (rc + s),
        }
    }

    /// The area inside the outline offset by s.
    fn plan_area(&self, s: f64) -> f64 {
        match self.plan {
            Plan::Round(half) => PI * (half + s).powi(2),
            Plan::Rect { hw, hl, rc, .. } => {
                4.0 * (hw + s) * (hl + s) - (4.0 - PI) * (rc + s).powi(2)
            }
        }
    }

    fn turn(&self) -> f64 {
        match self.plan {
            Plan::Rect { angle, .. } => angle,
            Plan::Round(_) => 0.0,
        }
    }

    /// A point in the emboss's own frame placed in the flat.
    fn place(&self, p: P2) -> P2 {
        let (s, c) = self.turn().sin_cos();
        [
            self.center[0] + c * p[0] - s * p[1],
            self.center[1] + s * p[0] + c * p[1],
        ]
    }

    fn place_dir(&self, d: P2) -> P3 {
        let (s, c) = self.turn().sin_cos();
        [c * d[0] - s * d[1], s * d[0] + c * d[1], 0.0]
    }

    fn segment(&self, j: usize) -> Segment {
        let Plan::Rect { hw, hl, rc, .. } = self.plan else {
            unreachable!()
        };
        let (a, b) = (hw - rc, hl - rc);
        // Sides: start, outward; corners: centre, start angle.
        let sides = [
            ([hw, -b], [1.0, 0.0], [0.0, 1.0]),
            ([a, hl], [0.0, 1.0], [-1.0, 0.0]),
            ([-hw, b], [-1.0, 0.0], [0.0, -1.0]),
            ([-a, -hl], [0.0, -1.0], [1.0, 0.0]),
        ];
        let corners = [
            ([a, b], 0.0),
            ([-a, b], FRAC_PI_2),
            ([-a, -b], PI),
            ([a, -b], 1.5 * PI),
        ];
        if j % 2 == 0 {
            let (start, n, d) = sides[j / 2];
            Segment {
                start,
                n,
                d,
                corner: None,
            }
        } else {
            let (c, a0) = corners[j / 2];
            Segment {
                start: [c[0] + rc * a0.cos(), c[1] + rc * a0.sin()],
                n: [a0.cos(), a0.sin()],
                d: [0.0, 0.0],
                corner: Some((c, a0)),
            }
        }
    }

    /// A fraction f along segment j, offset by s: in the flat, and the
    /// outward direction there.
    fn along(&self, j: usize, f: f64, s: f64) -> (P2, P2) {
        let g = self.segment(j);
        match g.corner {
            None => {
                let next = self.segment((j + 1) % 8).start;
                let p = [
                    g.start[0] + f * (next[0] - g.start[0]) + s * g.n[0],
                    g.start[1] + f * (next[1] - g.start[1]) + s * g.n[1],
                ];
                (p, g.n)
            }
            Some((c, a0)) => {
                let a = a0 + f * FRAC_PI_2;
                let n = [a.cos(), a.sin()];
                let rho = self.corner_radius() + s;
                ([c[0] + rho * n[0], c[1] + rho * n[1]], n)
            }
        }
    }

    fn flat_at(&self, j: usize, f: f64, s: f64) -> P2 {
        self.place(self.along(j, f, s).0)
    }

    /// The section point q on segment j a fraction f along, in 3D.
    fn point3(&self, j: usize, f: f64, q: P2) -> P3 {
        let p = self.flat_at(j, f, q[0]);
        [p[0], p[1], q[1]]
    }

    /// The material between the two faces and both faces' areas: the
    /// section swept along the outline (offset outline length P(s) = P0 +
    /// 2 pi s), and the raised face's slab. Base frame.
    pub fn mass(&self) -> Mass {
        let p0 = self.perimeter(0.0);
        let (mut volume, mut mz, mut area) = (0.0, 0.0, 0.0);
        for i in 0..3 {
            let [p, q] = [self.chains[0][i], self.chains[1][i]];
            // The region between the two pieces: along p, across to q's end,
            // back along q, across to p's start. Green's theorem gives the
            // integrals of 1, s, z and s z over it.
            let across = |a: P2, b: P2| Meridian::Line { a, b };
            let mut g = [0.0; 4];
            for (m, sign) in [
                (p, 1.0),
                (across(p.at(1.0).0, q.at(1.0).0), 1.0),
                (q, -1.0),
                (across(p.at(0.0).0, q.at(0.0).0), -1.0),
            ] {
                for (k, v) in green(m).into_iter().enumerate() {
                    g[k] += sign * v;
                }
            }
            // The loop runs either way round.
            let g = g.map(|v| v * g[0].signum());
            volume += p0 * g[0] + TAU * g[1];
            mz += p0 * g[2] + TAU * g[3];
            for m in [p, q] {
                area += length_moment(m, p0);
            }
        }
        let slab = self.plan_area(self.top);
        let t = (self.tops[1] - self.tops[0]).abs();
        volume += slab * t;
        mz += slab * t * (self.tops[0] + self.tops[1]) / 2.0;
        area += 2.0 * slab;
        Mass {
            volume,
            area,
            center_of_mass: [self.center[0], self.center[1], mz / volume],
        }
    }

    /// Builds a round emboss's two faces in `b`, on the footprint's circles
    /// (`circles[level]`, each closed at vertex `v[level]`, starting along
    /// +x).
    pub fn build_round(
        &self,
        b: &mut Builder,
        circles: [u32; 2],
        v: [u32; 2],
    ) -> Result<(), GeometryError> {
        let Plan::Round(half) = self.plan else {
            unreachable!()
        };
        let [cx, cy] = self.center;
        for level in 0..2 {
            let chain = self.chains[level].map(|m| m.shifted(half));
            let other = self.chains[1 - level].map(|m| m.shifted(half));
            // Circles and their vertices where the pieces meet, from the
            // footprint inward.
            let mut rings = vec![(circles[level], v[level])];
            for m in &chain {
                let q = m.at(1.0).0;
                let vertex = b.vertex([cx + q[0], cy, q[1]]);
                let circle = b.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: [cx, cy, q[1]],
                            x: [1.0, 0.0, 0.0],
                            y: [0.0, 1.0, 0.0],
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
                self.band(
                    b,
                    &format!("{}.{level}.{i}", self.id),
                    *m,
                    out,
                    rings[i],
                    rings[i + 1],
                )?;
            }
            // The raised face.
            let up = (self.tops[level] - self.tops[1 - level]).signum();
            let frame = Frame3 {
                origin: [cx, cy, self.tops[level]],
                x: [1.0, 0.0, 0.0],
                y: [0.0, up, 0.0],
                z: [0.0, 0.0, up],
            };
            let (circle, vertex) = rings[3];
            let use_ = plane_boundary(b, frame, circle, vertex, vertex, Orientation::Forward)?;
            let r = half + self.top + 1.0;
            b.face(
                &format!("{}.{level}.top", self.id),
                SurfaceGeometry::Plane { frame },
                [[-r, r]; 2],
                vec![use_],
            )?;
        }
        Ok(())
    }

    /// A full turn of a corner or the wall (in r, z) between the circles at
    /// its two ends, with a seam along it at angle 0.
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
        let x = [1.0, 0.0, 0.0];
        // The surface, its parameter around the axis (u = u0 + du * circle
        // parameter), each end's other parameter, the seam edge and its
        // pcurve (origin, direction) at u = 0, and whether the seam runs from
        // end a.
        let (surface, (u0, du), [ha, hb], seam, seam_v, forward) = match m {
            Meridian::Arc { c, rho, a0, a1 } => {
                let surface = SurfaceGeometry::Torus {
                    frame: Frame3 {
                        origin: [cx, cy, c[1]],
                        x,
                        y: [0.0, 1.0, 0.0],
                        z: Z,
                    },
                    major_radius: c[0],
                    minor_radius: rho,
                };
                let seam = b.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: [cx + c[0], cy, c[1]],
                            x,
                            y: Z,
                            z: [0.0, -1.0, 0.0],
                        },
                        radius: rho,
                    },
                    Interval::new(a0.min(a1), a0.max(a1))?,
                    false,
                );
                let seam_v = ([0.0, 0.0], [0.0, 1.0]);
                (surface, (0.0, 1.0), [a0, a1], seam, seam_v, a0 < a1)
            }
            Meridian::Line { a, b: e } => {
                let len = (e[0] - a[0]).hypot(e[1] - a[1]);
                let dir = [(e[0] - a[0]) / len, (e[1] - a[1]) / len];
                let seam = b.edge(
                    CurveGeometry::Line {
                        origin: [cx + a[0], cy, a[1]],
                        direction: [dir[0], 0.0, dir[1]],
                    },
                    Interval::new(0.0, len)?,
                    false,
                );
                let (surface, up, apex) = revolved([cx, cy, 0.0], x, a, dir);
                let turn = if up > 0.0 { (0.0, 1.0) } else { (TAU, -1.0) };
                let (va_, vb_) = (up * (a[1] - apex), up * (e[1] - apex));
                let seam_v = ([0.0, va_], [0.0, up * dir[1]]);
                (surface, turn, [va_, vb_], seam, seam_v, true)
            }
        };
        // Name the ends in the seam's order.
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

    /// Builds a rectangular emboss's two faces in `b`, on the footprint's
    /// edges: per level, the vertex at each of the eight segments' starts and
    /// each segment's edge (id, start vertex).
    pub fn build_rect(
        &self,
        b: &mut Builder,
        vertices: [[u32; 8]; 2],
        edges: [[(u32, u32); 8]; 2],
    ) -> Result<(), GeometryError> {
        let Plan::Rect { rc, .. } = self.plan else {
            unreachable!()
        };
        for level in 0..2 {
            let chain = self.chains[level];
            let other = self.chains[1 - level];
            // Section points where the pieces meet, from the footprint in.
            let qs: Vec<P2> = std::iter::once(chain[0].at(0.0).0)
                .chain(chain.iter().map(|m| m.at(1.0).0))
                .collect();
            // Vertices [breakpoint][segment start] and the outline's edges
            // there [breakpoint][segment], as (id, start vertex).
            let mut vs = vec![vertices[level]];
            let mut rings = vec![edges[level]];
            for &q in &qs[1..] {
                let v: [u32; 8] = std::array::from_fn(|j| b.vertex(self.point3(j, 0.0, q)));
                let e: [(u32, u32); 8] = std::array::from_fn(|j| {
                    let curve_range = match self.segment(j).corner {
                        None => {
                            let (p, n) = (self.point3(j, 0.0, q), self.point3((j + 1) % 8, 0.0, q));
                            let len = norm(sub(n, p));
                            (
                                CurveGeometry::Line {
                                    origin: p,
                                    direction: sub(n, p).map(|x| x / len),
                                },
                                (0.0, len),
                            )
                        }
                        Some((c, a0)) => {
                            let c = self.place(c);
                            let x = self.place_dir([1.0, 0.0]);
                            (
                                CurveGeometry::Circle {
                                    frame: Frame3 {
                                        origin: [c[0], c[1], q[1]],
                                        x,
                                        y: cross(Z, x),
                                        z: Z,
                                    },
                                    radius: rc + q[0],
                                },
                                (a0, a0 + FRAC_PI_2),
                            )
                        }
                    };
                    let id = b.edge(
                        curve_range.0,
                        Interval::new(curve_range.1 .0, curve_range.1 .1).unwrap(),
                        false,
                    );
                    (id, v[j])
                });
                vs.push(v);
                rings.push(e);
            }
            // The pieces' sections at each segment start, (id, start vertex).
            let mut across = Vec::new();
            for (i, m) in chain.iter().enumerate() {
                let e: [(u32, u32); 8] = std::array::from_fn(|j| {
                    let (_, n) = self.along(j, 0.0, 0.0);
                    let n3 = self.place_dir(n);
                    match *m {
                        Meridian::Arc { c, rho, a0, a1 } => {
                            let id = b.edge(
                                CurveGeometry::Circle {
                                    frame: Frame3 {
                                        origin: self.point3(j, 0.0, c),
                                        x: n3,
                                        y: Z,
                                        z: cross(n3, Z),
                                    },
                                    radius: rho,
                                },
                                Interval::new(a0.min(a1), a0.max(a1)).unwrap(),
                                false,
                            );
                            (id, if a0 < a1 { vs[i][j] } else { vs[i + 1][j] })
                        }
                        Meridian::Line { a, b: e } => {
                            let (p, q) = (self.point3(j, 0.0, a), self.point3(j, 0.0, e));
                            let len = norm(sub(q, p));
                            let id = b.edge(
                                CurveGeometry::Line {
                                    origin: p,
                                    direction: sub(q, p).map(|x| x / len),
                                },
                                Interval::new(0.0, len).unwrap(),
                                false,
                            );
                            (id, vs[i][j])
                        }
                    }
                });
                across.push(e);
            }
            for (i, m) in chain.iter().enumerate() {
                let (here, there) = (m.at(0.5).0, other[i].at(0.5).0);
                let out = [here[0] - there[0], here[1] - there[1]];
                for j in 0..8 {
                    let next = (j + 1) % 8;
                    let surface = self.piece_surface(j, *m);
                    let mid = self.point3(j, 0.5, m.at(0.5).0);
                    let hint = surface.project(mid, None)?;
                    let mut uses = Vec::new();
                    for ((id, start), from, to) in [
                        (rings[i][j], vs[i][j], vs[i][next]),
                        (across[i][next], vs[i][next], vs[i + 1][next]),
                        (rings[i + 1][j], vs[i + 1][next], vs[i + 1][j]),
                        (across[i][j], vs[i + 1][j], vs[i][j]),
                    ] {
                        let sense = if from == start {
                            Orientation::Forward
                        } else {
                            Orientation::Reverse
                        };
                        uses.push(projected_boundary(
                            b,
                            &surface,
                            id,
                            from,
                            to,
                            sense,
                            Some(hint),
                        )?);
                    }
                    let mut bounds = [[f64::MAX, f64::MIN]; 2];
                    for fu in 0..=8 {
                        for fv in 0..=8 {
                            let p = self.point3(j, fu as f64 / 8.0, m.at(fv as f64 / 8.0).0);
                            let uv = surface.project(p, Some(hint))?;
                            for k in 0..2 {
                                bounds[k][0] = bounds[k][0].min(uv[k]);
                                bounds[k][1] = bounds[k][1].max(uv[k]);
                            }
                        }
                    }
                    let bounds = bounds.map(|[lo, hi]| [lo - 1e-6, hi + 1e-6]);
                    let (_, n) = self.along(j, 0.5, 0.0);
                    let n3 = self.place_dir(n);
                    let outward = [n3[0] * out[0], n3[1] * out[0], out[1]];
                    let normal = surface.normal_at(hint)?;
                    let face = b.face(
                        &format!("{}.{level}.{i}.{j}", self.id),
                        surface,
                        bounds,
                        uses,
                    )?;
                    if dot(normal, outward) < 0.0 {
                        b.set_face_sense(face, Orientation::Reverse)?;
                    }
                }
            }
            // The raised face, inside the last ring.
            let up = (self.tops[level] - self.tops[1 - level]).signum();
            let x = self.place_dir([1.0, 0.0]);
            let frame = Frame3 {
                origin: [self.center[0], self.center[1], self.tops[level]],
                x,
                y: cross([0.0, 0.0, up], x),
                z: [0.0, 0.0, up],
            };
            let mut uses = Vec::new();
            for j in 0..8 {
                let (id, start) = rings[3][j];
                let (from, to) = (vs[3][j], vs[3][(j + 1) % 8]);
                let sense = if from == start {
                    Orientation::Forward
                } else {
                    Orientation::Reverse
                };
                uses.push(plane_boundary(b, frame, id, from, to, sense)?);
            }
            let Plan::Rect { hw, hl, .. } = self.plan else {
                unreachable!()
            };
            let r = hw.max(hl) + 1.0;
            b.face(
                &format!("{}.{level}.top", self.id),
                SurfaceGeometry::Plane { frame },
                [[-r, r]; 2],
                uses,
            )?;
        }
        Ok(())
    }

    /// The surface a section piece sweeps along segment j: along a side a
    /// cylinder or a plane, around a corner a torus or a cone (cylinder).
    /// Angle parameters are centred on the face, away from the seam.
    fn piece_surface(&self, j: usize, m: Meridian) -> SurfaceGeometry {
        let seg = self.segment(j);
        let (_, n_mid) = self.along(j, 0.5, 0.0);
        let away = self.place_dir([-n_mid[0], -n_mid[1]]);
        match (seg.corner, m) {
            (None, Meridian::Arc { c, rho, a0, a1 }) => {
                let n3 = self.place_dir(seg.n);
                let mid = (a0 + a1) / 2.0;
                let x = [-(mid.cos() * n3[0]), -(mid.cos() * n3[1]), -mid.sin()];
                let z = self.place_dir(seg.d);
                SurfaceGeometry::Cylinder {
                    frame: Frame3 {
                        origin: self.point3(j, 0.0, c),
                        x,
                        y: cross(z, x),
                        z,
                    },
                    radius: rho,
                }
            }
            (None, Meridian::Line { a, b }) => {
                let (p, q) = (self.point3(j, 0.0, a), self.point3(j, 0.0, b));
                let x = self.place_dir(seg.d);
                let y = sub(q, p).map(|v| v / norm(sub(q, p)));
                SurfaceGeometry::Plane {
                    frame: Frame3 {
                        origin: p,
                        x,
                        y,
                        z: cross(x, y),
                    },
                }
            }
            (Some((c, _)), Meridian::Arc { c: k, rho, .. }) => {
                let c = self.place(c);
                SurfaceGeometry::Torus {
                    frame: Frame3 {
                        origin: [c[0], c[1], k[1]],
                        x: away,
                        y: cross(Z, away),
                        z: Z,
                    },
                    major_radius: self.corner_radius() + k[0],
                    minor_radius: rho,
                }
            }
            (Some((c, _)), Meridian::Line { a, b }) => {
                let c = self.place(c);
                let rc = self.corner_radius();
                let len = (b[0] - a[0]).hypot(b[1] - a[1]);
                let dir = [(b[0] - a[0]) / len, (b[1] - a[1]) / len];
                revolved([c[0], c[1], 0.0], away, [rc + a[0], a[1]], dir).0
            }
        }
    }
}

/// The surface the line through q (r, z) along dir turns into around the
/// vertical axis through `axis`, its angle measured from `x`: a cylinder, or
/// a cone opening the way r grows. With the cone's axis direction (+1 up)
/// and its apex height.
fn revolved(axis: P3, x: P3, q: P2, dir: P2) -> (SurfaceGeometry, f64, f64) {
    if dir[0].abs() < EPS {
        let surface = SurfaceGeometry::Cylinder {
            frame: Frame3 {
                origin: axis,
                x,
                y: cross(Z, x),
                z: Z,
            },
            radius: q[0],
        };
        return (surface, 1.0, 0.0);
    }
    let apex = q[1] - q[0] / dir[0] * dir[1];
    let up = (dir[0] * dir[1]).signum();
    let z = [0.0, 0.0, up];
    let surface = SurfaceGeometry::Cone {
        frame: Frame3 {
            origin: [axis[0], axis[1], apex],
            x,
            y: cross(z, x),
            z,
        },
        semi_angle: (dir[0].abs() / dir[1].abs()).atan(),
    };
    (surface, up, apex)
}

/// Over the piece, by Green's theorem, the integrals over the region it
/// bounds of 1, s, z and s z: of s dz, s^2 / 2 dz, s z dz, s^2 z / 2 dz.
fn green(m: Meridian) -> [f64; 4] {
    let mut g = [0.0; 4];
    let parts = 8;
    for k in 0..parts {
        for (node, w) in GAUSS {
            let (p, d) = m.at((k as f64 + node) / parts as f64);
            let w = w / parts as f64 * d[1];
            let (s, z) = (p[0], p[1]);
            g[0] += s * w;
            g[1] += s * s / 2.0 * w;
            g[2] += s * z * w;
            g[3] += s * s * z / 2.0 * w;
        }
    }
    g
}

/// The area the piece sweeps along the outline: the integral of (p0 + 2 pi
/// s) along it.
fn length_moment(m: Meridian, p0: f64) -> f64 {
    let mut a = 0.0;
    let parts = 8;
    for k in 0..parts {
        for (node, w) in GAUSS {
            let (p, d) = m.at((k as f64 + node) / parts as f64);
            a += (p0 + TAU * p[0]) * d[0].hypot(d[1]) * w / parts as f64;
        }
    }
    a
}

fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
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
    use super::{Formed, Meridian};
    use crate::spec::{Emboss, PartSpec};
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

    fn rect(size: [f64; 3], height: f64, draft: f64, side: &str, angle: f64) -> String {
        format!(
            r#"{{"id": "r", "at": [50, 30], "shape": {{"rectangle": {{"width": {}, "length": {},
            "corner_radius": {}}}}}, "height": {height}, "draft": {draft}, "side": "{side}",
            "angle": {angle}}}"#,
            size[0], size[1], size[2]
        )
    }

    /// The solid's mesh agrees with the exact mass; down is up mirrored about
    /// the plate's middle.
    fn exact_solid(embosses: &[(String, &str)]) {
        let mut seen = Vec::new();
        for (e, what) in embosses {
            let spec = part(e, "");
            let solid = crate::fold::fold(&spec).unwrap_or_else(|err| panic!("{what}: {err}"));
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(
                rel(mesh.volume, exact.volume) < 2e-5,
                "{what}: {} vs {}",
                mesh.volume,
                exact.volume
            );
            assert!(
                rel(mesh.area, exact.area) < 2e-5,
                "{what}: area {} vs {}",
                mesh.area,
                exact.area
            );
            for j in 0..3 {
                let d = (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs();
                assert!(
                    d < 1e-3,
                    "{what}: {:?} vs {:?}",
                    mesh.center_of_mass,
                    exact.center_of_mass
                );
            }
            seen.push((exact.volume, exact.center_of_mass[2]));
        }
        assert!((seen[0].0 - seen[1].0).abs() < 1e-9);
        assert!((seen[0].1 + seen[1].1 - 2.0).abs() < 1e-9);
    }

    #[test]
    fn a_round_emboss_is_an_exact_solid() {
        exact_solid(&[
            (round([50.0, 30.0], 30.0, 8.0, 30.0, "up"), "up 30"),
            (round([50.0, 30.0], 30.0, 8.0, 30.0, "down"), "down 30"),
            (round([50.0, 30.0], 30.0, 8.0, 0.0, "up"), "up 0"),
            (round([50.0, 30.0], 30.0, 8.0, 10.0, "down"), "down 10"),
        ]);
        // The plate stays whole in the flat; the footprint is marked.
        let f = crate::flat::flat_pattern(&part(&round([50.0, 30.0], 30.0, 6.0, 30.0, "up"), ""))
            .unwrap();
        assert!((f.area() - 6000.0).abs() < 1e-9 && f.holes.is_empty() && f.cutouts.is_empty());
        assert_eq!(f.forms.len(), 1);
    }

    #[test]
    fn a_rectangular_emboss_is_an_exact_solid() {
        exact_solid(&[
            (rect([40.0, 24.0, 8.0], 5.0, 20.0, "up", 0.0), "up 20"),
            (rect([40.0, 24.0, 8.0], 5.0, 20.0, "down", 0.0), "down 20"),
            (rect([40.0, 24.0, 9.0], 7.0, 0.0, "up", 0.0), "up 0"),
            (rect([36.0, 20.0, 9.0], 5.0, 35.0, "up", 30.0), "turned"),
        ]);
        let f =
            crate::flat::flat_pattern(&part(&rect([40.0, 24.0, 8.0], 5.0, 20.0, "up", 0.0), ""))
                .unwrap();
        assert!((f.area() - 6000.0).abs() < 1e-9 && f.cutouts.is_empty());
        assert_eq!(f.forms[0].outline.points.len(), 8);
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
        // Plan corners too tight for the section, or past half a side; and
        // a rectangle across the outline.
        assert!(invalid(&rect([40.0, 24.0, 3.0], 5.0, 20.0, "up", 0.0), ""));
        assert!(invalid(&rect([40.0, 24.0, 13.0], 5.0, 20.0, "up", 0.0), ""));
        let across =
            crate::flat::flat_pattern(&part(&rect([98.0, 24.0, 8.0], 5.0, 20.0, "up", 0.0), ""));
        assert!(matches!(across, Err(SheetError::Invalid(_))), "{across:?}");
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

    /// The section meets its definition: the raised face's wall, extended,
    /// meets its top at the mold line, the height above the plate's face, at
    /// the draft; the sheet keeps its thickness.
    #[test]
    fn the_profile_keeps_its_measures() {
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
            let plate = if side == "up" { 2.0 } else { 0.0 };
            let top = f.tops[raised];
            assert!((top - (plate + up * 9.0)).abs() < 1e-12, "{side}: height");
            let s = (top - a[1]) / (b[1] - a[1]);
            let mold = a[0] + s * (b[0] - a[0]);
            assert!(mold.abs() < 1e-12, "{side}: mold line {mold} off");
            let slope = ((a[0] - b[0]) / (b[1] - a[1]) * up).atan().to_degrees();
            assert!((slope - draft).abs() < 1e-9, "{side}: draft {slope}");
            assert!(((f.tops[raised] - f.tops[other]).abs() - 2.0).abs() < 1e-12);
            for i in 0..3 {
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
