//! The folded part as an exact B-rep, built face by face from the layout
//! (no booleans). Every flat cell maps to 3D by its own fold map: the base
//! stays, a bend strip wraps onto a cylinder around the bend axis, a flange
//! strip turns rigidly about that axis (lips chain on their parent's map).
//!
//! Flat frame: base cell in XY, sheet from z = 0 to z = T.

use crate::flat::{layout, Curve, Cut, Fold, Layout, Zone, P2};
use crate::spec::PartSpec;
use crate::SheetError;
use opengeometry::brep::{
    boundary, conoid_ellipse, plane_boundary, projected_boundary, uv_line, Accuracy, BrepEnvelope,
    Builder, CurveGeometry, Frame3, GeometryError, Orientation, PcurveGeometry, SurfaceGeometry,
    Use, Wrap,
};
use opengeometry::math::Interval;
use opengeometry::tessellation::tessellate;
use std::collections::HashMap;
use std::f64::consts::TAU;

pub type P3 = [f64; 3];

/// Model accuracy in millimetres.
pub fn accuracy() -> Accuracy {
    Accuracy {
        geometric: 1e-6,
        intersection: 1e-7,
        tessellation: 1e-3,
        exchange: 1e-5,
    }
}

impl From<GeometryError> for SheetError {
    fn from(e: GeometryError) -> Self {
        SheetError::Kernel(e.to_string())
    }
}

pub fn fold(spec: &PartSpec) -> Result<BrepEnvelope, SheetError> {
    let spec = &*spec.resolved()?;
    let layout = layout(spec)?;
    let model = Model::new(spec, &layout);
    model.refuse_clashes()?;
    model.build().map_err(|e| match e {
        // What the kernel sees when features fold into each other.
        SheetError::Kernel(m) if m.contains("pointing into the solid") || m.contains("nonmanifold") => {
            SheetError::Kernel(format!(
                "the folded part isn't a valid solid; most likely features fold into each other ({m})"
            ))
        }
        e => e,
    })
}

/// Volume, surface area and centre of mass, from a tessellation with the
/// given chord deviation (mm).
#[derive(Clone, Copy, Debug)]
pub struct Mass {
    pub volume: f64,
    pub area: f64,
    pub center_of_mass: P3,
}

pub fn mass(brep: &BrepEnvelope, deflection: f64) -> Result<Mass, SheetError> {
    let mesh = tessellate(brep, deflection, 10_000_000)?;
    let p = |i: u32| -> P3 { std::array::from_fn(|j| mesh.positions[i as usize * 3 + j]) };
    let (mut volume, mut area, mut moment) = (0.0, 0.0, [0.0; 3]);
    for t in mesh.indices.chunks_exact(3) {
        let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
        let v = dot(a, cross(b, c)) / 6.0;
        volume += v;
        area += norm(cross(sub(b, a), sub(c, a))) / 2.0;
        for j in 0..3 {
            moment[j] += v * (a[j] + b[j] + c[j]) / 4.0;
        }
    }
    Ok(Mass {
        volume,
        area,
        center_of_mass: moment.map(|m| m / volume),
    })
}

/// Rigid map p -> r p + t (r row-major).
#[derive(Clone, Copy)]
struct Rigid {
    r: [[f64; 3]; 3],
    t: P3,
}

impl Rigid {
    const IDENTITY: Rigid = Rigid {
        r: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        t: [0.0; 3],
    };
    fn rot(&self, v: P3) -> P3 {
        std::array::from_fn(|i| dot(self.r[i], v))
    }
    fn apply(&self, p: P3) -> P3 {
        add(self.rot(p), self.t)
    }
    /// self after other: p -> self(other(p)).
    fn after(&self, other: &Rigid) -> Rigid {
        let r = std::array::from_fn(|i| {
            std::array::from_fn(|j| (0..3).map(|k| self.r[i][k] * other.r[k][j]).sum())
        });
        Rigid {
            r,
            t: self.apply(other.t),
        }
    }
    /// Rotation by `angle` about the line through `c` with unit direction `a`.
    fn about(c: P3, a: P3, angle: f64) -> Rigid {
        let (s, k) = (angle.sin(), 1.0 - angle.cos());
        let r = [
            [
                1.0 - k * (a[1] * a[1] + a[2] * a[2]),
                k * a[0] * a[1] - s * a[2],
                k * a[0] * a[2] + s * a[1],
            ],
            [
                k * a[0] * a[1] + s * a[2],
                1.0 - k * (a[0] * a[0] + a[2] * a[2]),
                k * a[1] * a[2] - s * a[0],
            ],
            [
                k * a[0] * a[2] - s * a[1],
                k * a[1] * a[2] + s * a[0],
                1.0 - k * (a[0] * a[0] + a[1] * a[1]),
            ],
        ];
        let rigid = Rigid { r, t: [0.0; 3] };
        Rigid {
            r,
            t: sub(c, rigid.rot(c)),
        }
    }
    fn shift(d: P3) -> Rigid {
        Rigid {
            t: d,
            ..Rigid::IDENTITY
        }
    }
}

/// One fold in 3D: the map of the flange after it and its bend cylinder.
pub(crate) struct FoldMap {
    /// Map of the zone this fold sits in.
    parent: Rigid,
    /// Map of the flange strip after the bend.
    child: Rigid,
    pub(crate) fold: Fold,
    /// +1 bends toward +z of the parent zone, -1 toward -z.
    sigma: f64,
    /// Height of the bend axis in the parent's flat frame.
    zc: f64,
}

pub(crate) struct Model<'a> {
    pub(crate) t: f64,
    r: f64,
    pub(crate) layout: &'a Layout,
    pub(crate) folds: Vec<FoldMap>,
}

impl<'a> Model<'a> {
    pub(crate) fn new(spec: &PartSpec, layout: &'a Layout) -> Self {
        let t = spec.thickness;
        // Largest bend radius, for padding UV bounds.
        let r = layout
            .folds
            .iter()
            .map(|f| f.radius)
            .fold(spec.bend_radius, f64::max);
        let mut folds: Vec<FoldMap> = Vec::new();
        for f in &layout.folds {
            let parent = f.parent.map_or(Rigid::IDENTITY, |p| folds[p].child);
            let sigma = if f.up { 1.0 } else { -1.0 };
            let zc = if f.up { t + f.radius } else { -f.radius };
            let axis = [
                f.edge_point[0] + f.s0 * f.n[0],
                f.edge_point[1] + f.s0 * f.n[1],
                zc,
            ];
            let u = [f.u[0], f.u[1], 0.0];
            // Turning n toward +z is a negative rotation about u (u x n = -z).
            let child = parent
                .after(&Rigid::about(axis, u, -sigma * f.angle))
                .after(&Rigid::shift([
                    -f.bend_allowance * f.n[0],
                    -f.bend_allowance * f.n[1],
                    0.0,
                ]));
            folds.push(FoldMap {
                parent,
                child,
                fold: f.clone(),
                sigma,
                zc,
            });
        }
        Model {
            t,
            r,
            layout,
            folds,
        }
    }

    /// Radius of the sheet layer at height z in a bend.
    fn radius(&self, f: &FoldMap, z: f64) -> f64 {
        if f.sigma > 0.0 {
            self.t + f.fold.radius - z
        } else {
            f.fold.radius + z
        }
    }

    /// Bend angle at flat point p.
    fn phi(&self, f: &FoldMap, p: P2) -> f64 {
        let s = dot2(sub2(p, f.fold.edge_point), f.fold.n);
        (s - f.fold.s0) / (f.fold.radius + f.fold.k_factor * self.t)
    }

    /// The bend cylinder's frame: point = origin + r (cos phi x + sin phi y) + v z.
    fn cylinder_frame(&self, f: &FoldMap) -> Frame3 {
        let n = [f.fold.n[0], f.fold.n[1], 0.0];
        let x = [0.0, 0.0, -f.sigma];
        let origin = [
            f.fold.edge_point[0] + f.fold.s0 * f.fold.n[0],
            f.fold.edge_point[1] + f.fold.s0 * f.fold.n[1],
            f.zc,
        ];
        Frame3 {
            origin: f.parent.apply(origin),
            x: f.parent.rot(x),
            y: f.parent.rot(n),
            z: f.parent.rot(cross(x, n)),
        }
    }

    pub(crate) fn map(&self, zone: Zone, p: P2, z: f64) -> P3 {
        match zone {
            Zone::Base => [p[0], p[1], z],
            Zone::Flange(i) => self.folds[i].child.apply([p[0], p[1], z]),
            Zone::Bend(i) => {
                let f = &self.folds[i];
                let frame = self.cylinder_frame(f);
                let phi = self.phi(f, p);
                let v = dot2(sub2(p, f.fold.edge_point), f.fold.u);
                let (s, c) = phi.sin_cos();
                // v runs along frame.z, which is -u for up bends and +u for down.
                let along = v * -f.sigma;
                frame.point([self.radius(f, z) * c, self.radius(f, z) * s, along])
            }
        }
    }

    /// The flat point a point of the zone's sheet in 3D comes from.
    pub(crate) fn unmap(&self, zone: Zone, x: P3) -> P2 {
        match zone {
            Zone::Bend(i) => {
                let f = &self.folds[i];
                let frame = self.cylinder_frame(f);
                let d = sub(x, frame.origin);
                let (lx, ly) = (dot(d, frame.x), dot(d, frame.y));
                // Angles from the middle of the bend, so a hem's 180 degrees don't wrap.
                let (s, c) = (f.fold.angle / 2.0).sin_cos();
                let phi = f.fold.angle / 2.0 + (ly * c - lx * s).atan2(lx * c + ly * s);
                let along = f.fold.s0 + phi * (f.fold.radius + f.fold.k_factor * self.t);
                let v = -f.sigma * dot(d, frame.z);
                let (e, n, u) = (f.fold.edge_point, f.fold.n, f.fold.u);
                [
                    e[0] + along * n[0] + v * u[0],
                    e[1] + along * n[1] + v * u[1],
                ]
            }
            _ => {
                let m = self.rigid(zone);
                let d = sub(x, m.t);
                [
                    (0..3).map(|k| m.r[k][0] * d[k]).sum(),
                    (0..3).map(|k| m.r[k][1] * d[k]).sum(),
                ]
            }
        }
    }

    /// The rotation part of a flat zone's map.
    fn rigid(&self, zone: Zone) -> Rigid {
        match zone {
            Zone::Base => Rigid::IDENTITY,
            Zone::Flange(i) => self.folds[i].child,
            Zone::Bend(i) => self.folds[i].parent,
        }
    }

    /// A bend's wrap: its cylinder frame turned to the middle of the bend (the
    /// wrap's inverse takes angles within half a turn of that), with the
    /// neutral radius, where flat arc length is kept.
    fn wrap(&self, f: &FoldMap) -> Wrap {
        let c = self.cylinder_frame(f);
        let (s, co) = (f.fold.angle / 2.0).sin_cos();
        Wrap {
            frame: Frame3 {
                origin: c.origin,
                x: add(scale(c.x, co), scale(c.y, s)),
                y: sub(scale(c.y, co), scale(c.x, s)),
                z: c.z,
            },
            radius: f.fold.radius + f.fold.k_factor * self.t,
        }
    }

    /// A flat point of a bend in its wrap's flat space: arc length from the
    /// middle of the bend on the neutral layer, along the axis, radius.
    fn to_wrap(&self, f: &FoldMap, p: P2, z: f64) -> P3 {
        let d = sub2(p, f.fold.edge_point);
        [
            dot2(d, f.fold.n) - f.fold.s0 - f.fold.bend_allowance / 2.0,
            -f.sigma * dot2(d, f.fold.u),
            self.radius(f, z),
        ]
    }

    /// A flat direction in a bend's wrap flat space.
    fn to_wrap_dir(&self, f: &FoldMap, v: P2) -> P3 {
        [dot2(v, f.fold.n), -f.sigma * dot2(v, f.fold.u), 0.0]
    }

    /// The frame of a flat circle around c at height z (x along x_dir, z up
    /// the sheet), in 3D or, in a bend, in the bend's wrap flat space.
    fn circle_frame(&self, zone: Zone, c: P2, x_dir: P2, z: f64) -> (Option<Wrap>, Frame3) {
        let y_dir = [-x_dir[1], x_dir[0]];
        match zone {
            Zone::Bend(i) => {
                let f = &self.folds[i];
                (
                    Some(self.wrap(f)),
                    Frame3 {
                        origin: self.to_wrap(f, c, z),
                        x: self.to_wrap_dir(f, x_dir),
                        y: self.to_wrap_dir(f, y_dir),
                        z: [0.0, 0.0, -f.sigma],
                    },
                )
            }
            _ => {
                let m = self.rigid(zone);
                (
                    None,
                    Frame3 {
                        origin: self.map(zone, c, z),
                        x: m.rot([x_dir[0], x_dir[1], 0.0]),
                        y: m.rot([y_dir[0], y_dir[1], 0.0]),
                        z: m.rot([0.0, 0.0, 1.0]),
                    },
                )
            }
        }
    }

    /// A flat circle at height z: parameter 0 at x_dir, counter-clockwise.
    fn circle(&self, zone: Zone, c: P2, r: f64, x_dir: P2, z: f64) -> CurveGeometry {
        let (wrap, frame) = self.circle_frame(zone, c, x_dir, z);
        let circle = CurveGeometry::Circle { frame, radius: r };
        match wrap {
            Some(wrap) => CurveGeometry::Wrapped {
                wrap,
                curve: Box::new(circle),
            },
            None => circle,
        }
    }

    /// The wall through the sheet along a flat circle: UV = (the circle's
    /// parameter, height in the sheet).
    fn wall(&self, zone: Zone, c: P2, r: f64, x_dir: P2) -> SurfaceGeometry {
        let (wrap, frame) = self.circle_frame(zone, c, x_dir, 0.0);
        let cylinder = SurfaceGeometry::Cylinder { frame, radius: r };
        match wrap {
            Some(wrap) => SurfaceGeometry::Wrapped {
                wrap,
                surface: Box::new(cylinder),
            },
            None => cylinder,
        }
    }

    /// Where a flat circle's parameter starts. In a bend: along the fold
    /// normal, so its pcurve on the bend faces is an axis-aligned ellipse.
    fn circle_axis(&self, zone: Zone) -> P2 {
        match zone {
            Zone::Bend(i) => self.folds[i].fold.n,
            _ => [1.0, 0.0],
        }
    }

    /// How an edge along a flat circle sits on the faces of its zone.
    fn circle_kind(&self, zone: Zone, c: P2, r: f64, x_dir: P2) -> Kind {
        match zone {
            Zone::Bend(i) => {
                // Bend face UV: (angle from the bend start, along the axis);
                // the angle is wrap x / radius + half the bend angle.
                let f = &self.folds[i];
                let rn = f.fold.radius + f.fold.k_factor * self.t;
                let o = self.to_wrap(f, c, 0.0);
                let (x, y) = (
                    self.to_wrap_dir(f, x_dir),
                    self.to_wrap_dir(f, [-x_dir[1], x_dir[0]]),
                );
                Kind::Conic {
                    origin: [o[0] / rn + f.fold.angle / 2.0, o[1]],
                    a: [r * x[0] / rn, r * x[1]],
                    b: [r * y[0] / rn, r * y[1]],
                }
            }
            _ => Kind::Circle,
        }
    }

    fn build(&self) -> Result<BrepEnvelope, SheetError> {
        let (t, layout) = (self.t, self.layout);
        let mut b = Builder::new("gusset".into(), accuracy())?;
        let cells = &layout.cells;

        // Shared 2D vertices and edges.
        let key = |p: P2| ((p[0] * 1e6).round() as i64, (p[1] * 1e6).round() as i64);
        let mut ids: HashMap<(i64, i64), usize> = HashMap::new();
        let mut points: Vec<(P2, Zone)> = Vec::new();
        let mut loops: Vec<(Vec<usize>, Vec<Option<Curve>>)> = Vec::new();
        // Cut-outs a cell holds whole: (cell, its loop), clockwise.
        let mut inner: Vec<(usize, Vec<usize>, Vec<Option<Curve>>)> = Vec::new();
        let rings = cells.iter().enumerate().flat_map(|(c, cell)| {
            std::iter::once((c, &cell.polygon, &cell.arcs, false)).chain(
                cell.inner
                    .iter()
                    .map(move |r| (c, &r.points, &r.arcs, true)),
            )
        });
        for (c, polygon, arcs, is_inner) in rings {
            let cell = &cells[c];
            let (mut lp, mut la): (Vec<usize>, Vec<Option<Curve>>) = (Vec::new(), Vec::new());
            for (k, &p) in polygon.iter().enumerate() {
                let id = *ids.entry(key(p)).or_insert_with(|| {
                    points.push((p, cell.zone));
                    points.len() - 1
                });
                if lp.last() == Some(&id) {
                    // A repeated point: the segment from it is the one that counts.
                    *la.last_mut().unwrap() = arcs[k];
                } else {
                    lp.push(id);
                    la.push(arcs[k]);
                }
            }
            if lp.first() == lp.last() && lp.len() > 1 {
                lp.pop();
                la.pop();
            }
            if is_inner {
                inner.push((c, lp, la));
            } else {
                loops.push((lp, la));
            }
        }
        // An undirected flat edge: its end points, and for an arc a point on
        // it (a hole crossing one line leaves two arcs between the same ends).
        let ek = |p: usize, q: usize, arc: Option<Curve>| -> EdgeKey {
            let mid = arc.map(|c| {
                let a = match c {
                    Curve::Arc(a) => a,
                    Curve::Cut(cut) => {
                        let (pp, qq) = (points[p].0, points[q].0);
                        return key(cut.point((cut.s_of(pp) + cut.s_of(qq)) / 2.0));
                    }
                };
                let (pp, qq) = (points[p].0, points[q].0);
                let r = dot2(sub2(pp, a.center), sub2(pp, a.center)).sqrt();
                let angle =
                    (pp[1] - a.center[1]).atan2(pp[0] - a.center[0]) + a.sweep(pp, qq) / 2.0;
                key([a.center[0] + r * angle.cos(), a.center[1] + r * angle.sin()])
            });
            (p.min(q), p.max(q), mid)
        };
        // Undirected edge -> (cells using it, owner cell, owner direction, arc).
        type Use2 = (usize, usize, (usize, usize), Option<Curve>);
        let mut uses: HashMap<EdgeKey, Use2> = HashMap::new();
        let all = loops
            .iter()
            .enumerate()
            .map(|(c, (lp, la))| (c, lp, la))
            .chain(inner.iter().map(|(c, lp, la)| (*c, lp, la)));
        for (c, lp, la) in all {
            for i in 0..lp.len() {
                let (p, q) = (lp[i], lp[(i + 1) % lp.len()]);
                uses.entry(ek(p, q, la[i]))
                    .and_modify(|e| e.0 += 1)
                    .or_insert((1, c, (p, q), la[i]));
            }
        }
        if let Some(bad) = uses.values().find(|e| e.0 > 2) {
            return Err(SheetError::Kernel(format!(
                "flat edge used by {} cells",
                bad.0
            )));
        }

        // 3D vertices, two per flat vertex.
        let vertex: Vec<[u32; 2]> = points
            .iter()
            .map(|(p, zone)| {
                [
                    b.vertex(self.map(*zone, *p, 0.0)),
                    b.vertex(self.map(*zone, *p, t)),
                ]
            })
            .collect();

        // 3D edges, two per flat edge; curves follow the owner cell's zone.
        let mut edges: HashMap<(EdgeKey, usize), Edge> = HashMap::new();
        // Walls along arcs: surface, parameter at the start vertex, sweep, start vertex.
        let mut arc_walls: HashMap<EdgeKey, (SurfaceGeometry, f64, f64, usize)> = HashMap::new();
        // Walls along mitered bend ends: fold, conoid coefficients, angle
        // range, the vertex at its start.
        let mut cut_walls: HashMap<EdgeKey, (usize, [f64; 3], f64, f64, usize)> = HashMap::new();
        for (&k, &(_, owner, (p, q), arc)) in &uses {
            let (lo, hi) = (k.0, k.1);
            let zone = cells[owner].zone;
            if let Some(Curve::Cut(cut)) = arc {
                // A mitered bend's end: the cut face is a conoid around the
                // bend's axis (along = -sigma t(phi)); its edges on the
                // bend's layers are ellipses.
                let Zone::Bend(fi) = zone else {
                    return Err(SheetError::Kernel("a mitered end outside a bend".into()));
                };
                let f = &self.folds[fi];
                let frame = self.cylinder_frame(f);
                let c3 = cut.c.map(|x| -f.sigma * x);
                let phi = |v: usize| cut.phi(cut.s_of(points[v].0));
                let (start, lo_phi, hi_phi) = if phi(p) <= phi(q) {
                    (p, phi(p), phi(q))
                } else {
                    (q, phi(q), phi(p))
                };
                for level in 0..2 {
                    let rho = self.radius(f, level as f64 * t);
                    let (ellipse, phase) = conoid_ellipse(frame, c3, rho)?;
                    let id = b.edge(
                        ellipse,
                        Interval::new(lo_phi - phase, hi_phi - phase).map_err(GeometryError::from)?,
                        false,
                    );
                    edges.insert(
                        (k, level),
                        Edge {
                            id,
                            start: vertex[start][level],
                            kind: Kind::Cut {
                                phi: [lo_phi, hi_phi],
                                phase,
                            },
                        },
                    );
                }
                cut_walls.insert(k, (fi, c3, lo_phi, hi_phi, start));
                continue;
            }
            if let Some(Curve::Arc(arc)) = arc {
                // Counter-clockwise from s, in the zone's own frame.
                let s = if arc.ccw { p } else { q };
                let sweep = arc.sweep(points[p].0, points[q].0).abs();
                let (c, sp) = (arc.center, points[s].0);
                let r = dot2(sub2(sp, c), sub2(sp, c)).sqrt();
                // Planar zones start the parameter at s; bends at their axis.
                let (x_dir, t0) = match zone {
                    Zone::Bend(_) => {
                        let x = self.circle_axis(zone);
                        let d = sub2(sp, c);
                        (x, dot2(d, [-x[1], x[0]]).atan2(dot2(d, x)))
                    }
                    _ => (unit2(sub2(sp, c)), 0.0),
                };
                for level in 0..2 {
                    let id = b.edge(
                        self.circle(zone, c, r, x_dir, level as f64 * t),
                        Interval::new(t0, t0 + sweep).map_err(GeometryError::from)?,
                        false,
                    );
                    edges.insert(
                        (k, level),
                        Edge {
                            id,
                            start: vertex[s][level],
                            kind: self.circle_kind(zone, c, r, x_dir),
                        },
                    );
                }
                arc_walls.insert(k, (self.wall(zone, c, r, x_dir), t0, sweep, s));
                continue;
            }
            for level in 0..2 {
                let z = level as f64 * t;
                let e = self.edge(
                    &mut b,
                    zone,
                    (points[lo].0, vertex[lo][level]),
                    (points[hi].0, vertex[hi][level]),
                    z,
                )?;
                edges.insert((k, level), e);
            }
        }
        // Vertical edges at outline vertices (straight, radial in a bend).
        let mut vertical: HashMap<usize, Edge> = HashMap::new();
        for (&(lo, hi, _), &(count, _, _, _)) in &uses {
            if count == 1 {
                for v in [lo, hi] {
                    vertical.entry(v).or_insert_with(|| {
                        let (p, zone) = points[v];
                        let (a, c) = (self.map(zone, p, 0.0), self.map(zone, p, t));
                        let id = b.edge(
                            CurveGeometry::Line {
                                origin: a,
                                direction: unit(sub(c, a)),
                            },
                            Interval::new(0.0, t).unwrap(),
                            false,
                        );
                        Edge {
                            id,
                            start: vertex[v][0],
                            kind: Kind::Line {
                                direction: unit(sub(c, a)),
                            },
                        }
                    });
                }
            }
        }

        // Holes a cell holds whole: circles top and bottom, a seam, the wall.
        let mut holes: Vec<HoleEdges> = Vec::new();
        for (ci, cell) in cells.iter().enumerate() {
            for h in &cell.holes {
                let zone = cell.zone;
                let r = h.diameter / 2.0;
                let x_dir = self.circle_axis(zone);
                let seam_point = [h.center[0] + r * x_dir[0], h.center[1] + r * x_dir[1]];
                let (a0, a1) = (
                    self.map(zone, seam_point, 0.0),
                    self.map(zone, seam_point, t),
                );
                let v = [b.vertex(a0), b.vertex(a1)];
                let mut circles = [0; 2];
                for (level, id) in circles.iter_mut().enumerate() {
                    *id = b.edge(
                        self.circle(zone, h.center, r, x_dir, level as f64 * t),
                        Interval::new(0.0, TAU).unwrap(),
                        false,
                    );
                }
                let seam = b.edge(
                    CurveGeometry::Line {
                        origin: a0,
                        direction: unit(sub(a1, a0)),
                    },
                    Interval::new(0.0, t).unwrap(),
                    true,
                );
                holes.push(HoleEdges {
                    cell: ci,
                    seam,
                    circles,
                    v,
                    wall: self.wall(zone, h.center, r, x_dir),
                    kind: self.circle_kind(zone, h.center, r, x_dir),
                });
            }
        }

        // Top and bottom faces of every cell.
        for (c, cell) in cells.iter().enumerate() {
            let zone = cell.zone;
            // Arcs can bulge past the loop's vertices by up to their radius.
            let pad = loops[c]
                .0
                .iter()
                .zip(&loops[c].1)
                .filter_map(|(&v, a)| {
                    a.map(|a| match a {
                        Curve::Arc(a) => {
                            dot2(sub2(points[v].0, a.center), sub2(points[v].0, a.center)).sqrt()
                        }
                        Curve::Cut(c) => c.c[1].abs() + c.c[2].abs(),
                    })
                })
                .fold(0.0, f64::max);
            for level in 0..2 {
                let z = level as f64 * t;
                let surface = self.surface(zone, z, level == 1);
                let lp = &loops[c].0;
                let mut outer = Vec::new();
                for i in 0..lp.len() {
                    let (p, q) = (lp[i], lp[(i + 1) % lp.len()]);
                    let e = edges[&(ek(p, q, loops[c].1[i]), level)];
                    outer.push(self.use_on(&b, &surface, e, vertex[p][level], vertex[q][level])?);
                }
                let mut hole_loops = Vec::new();
                for hole in holes.iter().filter(|h| h.cell == c) {
                    let e = Edge {
                        id: hole.circles[level],
                        start: hole.v[level],
                        kind: hole.kind,
                    };
                    hole_loops.push(vec![self.use_on(
                        &b,
                        &surface,
                        e,
                        hole.v[level],
                        hole.v[level],
                    )?]);
                }
                for (_, lp, la) in inner.iter().filter(|r| r.0 == c) {
                    let mut ring = Vec::new();
                    for i in 0..lp.len() {
                        let (p, q) = (lp[i], lp[(i + 1) % lp.len()]);
                        let e = edges[&(ek(p, q, la[i]), level)];
                        ring.push(self.use_on(
                            &b,
                            &surface,
                            e,
                            vertex[p][level],
                            vertex[q][level],
                        )?);
                    }
                    hole_loops.push(ring);
                }
                let key = format!(
                    "{}.{}",
                    cell.name,
                    if level == 1 { "top" } else { "bottom" }
                );
                // Holes held whole lie within the loop's own box.
                let bounds = surface.bounds(&b, lp.iter().map(|&v| vertex[v][level]));
                let bounds = bounds.map(|[lo, hi]| [lo - pad, hi + pad]);
                // A bend face spans exactly its bend angle in u (atan2 would wrap
                // past 180 degrees).
                let bounds = match zone {
                    Zone::Bend(i) => [[-1e-3, self.folds[i].fold.angle + 1e-3], bounds[1]],
                    _ => bounds,
                };
                let face =
                    b.face_with_holes(&key, surface.geometry(), bounds, outer, hole_loops)?;
                if surface.is_reversed() {
                    b.set_face_sense(face, Orientation::Reverse)?;
                }
            }
        }

        // Side walls along the outline.
        let mut side = 0;
        let mut outline: Vec<_> = uses
            .iter()
            .filter(|(_, e)| e.0 == 1)
            .map(|(k, e)| (e.1, e.2, *k))
            .collect();
        outline.sort();
        for (owner, (p, q), k2) in outline {
            let zone = cells[owner].zone;
            if let Some(&(fi, c3, lo_phi, hi_phi, start)) = cut_walls.get(&k2) {
                let f = &self.folds[fi];
                let frame = self.cylinder_frame(f);
                let (r0, r1) = (self.radius(f, 0.0), self.radius(f, t));
                let phi = |v: usize| if v == start { lo_phi } else { hi_phi };
                let (bottom, top) = (edges[&(k2, 0)], edges[&(k2, 1)]);
                let phase = |e: &Edge| match e.kind {
                    Kind::Cut { phase, .. } => phase,
                    _ => 0.0,
                };
                let along = |e: &Edge, from: u32| {
                    if from == e.start {
                        Orientation::Forward
                    } else {
                        Orientation::Reverse
                    }
                };
                let rate = (r1 - r0) / t;
                let wall = vec![
                    boundary(
                        bottom.id,
                        vertex[p][0],
                        vertex[q][0],
                        along(&bottom, vertex[p][0]),
                        uv_line([phase(&bottom), r0], [1.0, 0.0]),
                    ),
                    boundary(
                        vertical[&q].id,
                        vertex[q][0],
                        vertex[q][1],
                        Orientation::Forward,
                        uv_line([phi(q), r0], [0.0, rate]),
                    ),
                    boundary(
                        top.id,
                        vertex[q][1],
                        vertex[p][1],
                        along(&top, vertex[q][1]),
                        uv_line([phase(&top), r1], [1.0, 0.0]),
                    ),
                    boundary(
                        vertical[&p].id,
                        vertex[p][1],
                        vertex[p][0],
                        Orientation::Reverse,
                        uv_line([phi(p), r0], [0.0, rate]),
                    ),
                ];
                let face = b.face(
                    &format!("side.{side}"),
                    SurfaceGeometry::Conoid { frame, c: c3 },
                    [[lo_phi, hi_phi], [r0.min(r1), r0.max(r1)]],
                    wall,
                )?;
                // Its normal S_u x S_v at the middle against the way into the
                // cell (left of p -> q in the flat).
                let (um, vm) = ((lo_phi + hi_phi) / 2.0, (r0 + r1) / 2.0);
                let (sn, cs) = um.sin_cos();
                let su = add(
                    add(scale(frame.x, -vm * sn), scale(frame.y, vm * cs)),
                    scale(frame.z, -c3[1] * sn + c3[2] * cs),
                );
                let sv = add(scale(frame.x, cs), scale(frame.y, sn));
                let normal = cross(su, sv);
                let Some(Curve::Cut(cut)) = uses[&k2].3 else {
                    unreachable!()
                };
                let s_mid = cut.s0 + um * cut.radius;
                let m = cut.point(s_mid);
                let d = sub2(points[q].0, points[p].0);
                let left = unit2([-d[1], d[0]]).map(|x| x * 1e-4);
                let z_mid = t / 2.0;
                let inward = sub(
                    self.map(zone, [m[0] + left[0], m[1] + left[1]], z_mid),
                    self.map(zone, m, z_mid),
                );
                if dot(normal, inward) > 0.0 {
                    b.set_face_sense(face, Orientation::Reverse)?;
                }
                side += 1;
                continue;
            }
            if let Some((wall, t0, sweep, start)) = arc_walls.get(&k2) {
                // The cell lies left of p -> q: inside the circle when the arc
                // runs counter-clockwise.
                let ccw = uses[&k2].3.and_then(Curve::arc).unwrap().ccw;
                let phi = |v: usize| if v == *start { *t0 } else { t0 + sweep };
                let (bottom, top) = (edges[&(k2, 0)], edges[&(k2, 1)]);
                let along = |e: &Edge, from: u32| {
                    if from == e.start {
                        Orientation::Forward
                    } else {
                        Orientation::Reverse
                    }
                };
                let uses = vec![
                    boundary(
                        bottom.id,
                        vertex[p][0],
                        vertex[q][0],
                        along(&bottom, vertex[p][0]),
                        uv_line([0.0, 0.0], [1.0, 0.0]),
                    ),
                    boundary(
                        vertical[&q].id,
                        vertex[q][0],
                        vertex[q][1],
                        Orientation::Forward,
                        uv_line([phi(q), 0.0], [0.0, 1.0]),
                    ),
                    boundary(
                        top.id,
                        vertex[q][1],
                        vertex[p][1],
                        along(&top, vertex[q][1]),
                        uv_line([0.0, t], [1.0, 0.0]),
                    ),
                    boundary(
                        vertical[&p].id,
                        vertex[p][1],
                        vertex[p][0],
                        Orientation::Reverse,
                        uv_line([phi(p), 0.0], [0.0, 1.0]),
                    ),
                ];
                let face = b.face(
                    &format!("side.{side}"),
                    wall.clone(),
                    [[*t0, t0 + sweep], [0.0, t]],
                    uses,
                )?;
                if !ccw {
                    b.set_face_sense(face, Orientation::Reverse)?;
                }
                side += 1;
                continue;
            }
            if let Some(face) = self.slanted_wall(
                &mut b,
                zone,
                (points[p].0, vertex[p]),
                (points[q].0, vertex[q]),
                [edges[&(k2, 0)], edges[&(k2, 1)]],
                [vertical[&p], vertical[&q]],
                side,
            )? {
                let _ = face;
                side += 1;
                continue;
            }
            let (p0, q0) = (
                self.map(zone, points[p].0, 0.0),
                self.map(zone, points[q].0, 0.0),
            );
            // The cell lies left of p -> q, so outward is to the right.
            let d = sub2(points[q].0, points[p].0);
            let out = self
                .rigid_at(zone, points[p].0)
                .rot(unit([d[1], -d[0], 0.0]));
            let x = unit(sub(sub(q0, p0), scale(out, dot(sub(q0, p0), out))));
            let frame = Frame3 {
                origin: p0,
                x,
                y: cross(out, x),
                z: out,
            };
            let surface = Surface::Plane(frame);
            let uses = vec![
                self.use_on(&b, &surface, edges[&(k2, 0)], vertex[p][0], vertex[q][0])?,
                self.use_on(&b, &surface, vertical[&q], vertex[q][0], vertex[q][1])?,
                self.use_on(&b, &surface, edges[&(k2, 1)], vertex[q][1], vertex[p][1])?,
                self.use_on(&b, &surface, vertical[&p], vertex[p][1], vertex[p][0])?,
            ];
            let bounds = surface.bounds(
                &b,
                [vertex[p][0], vertex[q][0], vertex[q][1], vertex[p][1]].into_iter(),
            );
            // Arcs bulge past their end points: widen by the sagitta bound.
            // A bend end sector of up to 360 degrees reaches 2 (R + T) past its corners.
            let pad = 2.0 * (self.t + self.r);
            let bounds = [
                [bounds[0][0] - pad, bounds[0][1] + pad],
                [bounds[1][0] - pad, bounds[1][1] + pad],
            ];
            b.face(&format!("side.{side}"), surface.geometry(), bounds, uses)?;
            side += 1;
        }

        // Hole walls: UV = (circle parameter, height), seam at parameter 0.
        for (k, h) in holes.iter().enumerate() {
            let face = b.face(
                &format!("hole.{k}"),
                h.wall.clone(),
                [[0.0, TAU], [0.0, t]],
                vec![
                    boundary(
                        h.circles[0],
                        h.v[0],
                        h.v[0],
                        Orientation::Forward,
                        uv_line([0.0, 0.0], [1.0, 0.0]),
                    ),
                    boundary(
                        h.seam,
                        h.v[0],
                        h.v[1],
                        Orientation::Forward,
                        uv_line([TAU, 0.0], [0.0, 1.0]),
                    ),
                    boundary(
                        h.circles[1],
                        h.v[1],
                        h.v[1],
                        Orientation::Reverse,
                        uv_line([0.0, t], [1.0, 0.0]),
                    ),
                    boundary(
                        h.seam,
                        h.v[1],
                        h.v[0],
                        Orientation::Reverse,
                        uv_line([0.0, 0.0], [0.0, 1.0]),
                    ),
                ],
            )?;
            b.set_face_sense(face, Orientation::Reverse)?;
        }

        Ok(b.finish_solid()?)
    }

    /// The wall along a flat edge that crosses a bend at an angle: the flat
    /// vertical plane through it, wrapped. UV: (along p -> q, across the
    /// sheet). None for any other edge.
    #[allow(clippy::too_many_arguments)]
    fn slanted_wall(
        &self,
        b: &mut Builder,
        zone: Zone,
        (p, vp): (P2, [u32; 2]),
        (q, vq): (P2, [u32; 2]),
        levels: [Edge; 2],
        verticals: [Edge; 2],
        side: usize,
    ) -> Result<Option<u32>, SheetError> {
        let Zone::Bend(i) = zone else { return Ok(None) };
        let f = &self.folds[i];
        let d = unit2(sub2(q, p));
        if dot2(d, f.fold.u).abs() >= 1.0 - 1e-9 || dot2(d, f.fold.n).abs() >= 1.0 - 1e-9 {
            return Ok(None);
        }
        // The cell lies left of p -> q: outward is to the right.
        let x = self.to_wrap_dir(f, d);
        let z = self.to_wrap_dir(f, [d[1], -d[0]]);
        let frame = Frame3 {
            origin: self.to_wrap(f, p, 0.0),
            x,
            y: cross(z, x),
            z,
        };
        let uv = |w: P3| {
            let l = frame.to_local(w);
            [l[0], l[1]]
        };
        let pc = |e: &Edge, flat: P2, level: f64, dir: P3| {
            let o = uv(self.to_wrap(f, flat, level));
            let _ = e;
            uv_line(o, [dot(dir, frame.x), dot(dir, frame.y)])
        };
        let sense = |e: &Edge, from: u32| {
            if from == e.start {
                Orientation::Forward
            } else {
                Orientation::Reverse
            }
        };
        let t = self.t;
        let up = [0.0, 0.0, -f.sigma]; // wrap flat direction of a radial edge, per mm
                                       // Each level edge starts at its own start vertex: p or q.
        let level_use = |e: &Edge, level: usize, from: u32, to: u32| {
            let at_p = e.start == vp[level];
            let (start, dir) = if at_p { (p, x) } else { (q, scale(x, -1.0)) };
            boundary(
                e.id,
                from,
                to,
                sense(e, from),
                pc(e, start, level as f64 * t, dir),
            )
        };
        let uses = vec![
            level_use(&levels[0], 0, vp[0], vq[0]),
            boundary(
                verticals[1].id,
                vq[0],
                vq[1],
                Orientation::Forward,
                pc(&verticals[1], q, 0.0, up),
            ),
            level_use(&levels[1], 1, vq[1], vp[1]),
            boundary(
                verticals[0].id,
                vp[1],
                vp[0],
                Orientation::Reverse,
                pc(&verticals[0], p, 0.0, up),
            ),
        ];
        let corners = [
            uv(self.to_wrap(f, p, 0.0)),
            uv(self.to_wrap(f, q, 0.0)),
            uv(self.to_wrap(f, q, t)),
            uv(self.to_wrap(f, p, t)),
        ];
        let span = |j: usize| {
            let (lo, hi) = corners.iter().fold((f64::MAX, f64::MIN), |(lo, hi), c| {
                (lo.min(c[j]), hi.max(c[j]))
            });
            [lo, hi]
        };
        let face = b.face(
            &format!("side.{side}"),
            SurfaceGeometry::Wrapped {
                wrap: self.wrap(f),
                surface: Box::new(SurfaceGeometry::Plane { frame }),
            },
            [span(0), span(1)],
            uses,
        )?;
        Ok(Some(face))
    }

    /// A flat edge mapped to 3D at height z.
    fn edge(
        &self,
        b: &mut Builder,
        zone: Zone,
        (p, vp): (P2, u32),
        (q, vq): (P2, u32),
        z: f64,
    ) -> Result<Edge, SheetError> {
        let (a, c) = (self.map(zone, p, z), self.map(zone, q, z));
        if let Zone::Bend(i) = zone {
            let f = &self.folds[i];
            let d = unit2(sub2(q, p));
            if dot2(d, f.fold.u).abs() < 1.0 - 1e-9 {
                if dot2(d, f.fold.n).abs() < 1.0 - 1e-9 {
                    // Slanted: a wrapped line, i.e. a helix around the axis.
                    let (o, dir) = (self.to_wrap(f, p, z), self.to_wrap_dir(f, d));
                    let rn = f.fold.radius + f.fold.k_factor * self.t;
                    let id = b.edge(
                        CurveGeometry::Wrapped {
                            wrap: self.wrap(f),
                            curve: Box::new(CurveGeometry::Line {
                                origin: o,
                                direction: dir,
                            }),
                        },
                        Interval::new(0.0, norm2(sub2(q, p))).map_err(GeometryError::from)?,
                        false,
                    );
                    return Ok(Edge {
                        id,
                        start: vp,
                        kind: Kind::Helix {
                            origin: [o[0] / rn + f.fold.angle / 2.0, o[1]],
                            direction: [dir[0] / rn, dir[1]],
                        },
                    });
                }
                // A strip end: an arc around the bend axis.
                let frame = self.cylinder_frame(f);
                let v = dot(sub(a, frame.origin), frame.z);
                let (fp, fq) = (self.phi(f, p), self.phi(f, q));
                let id = b.edge(
                    CurveGeometry::Circle {
                        frame: Frame3 {
                            origin: add(frame.origin, scale(frame.z, v)),
                            ..frame
                        },
                        radius: self.radius(f, z),
                    },
                    Interval::new(fp.min(fq), fp.max(fq)).map_err(GeometryError::from)?,
                    false,
                );
                let start = if fp <= fq { vp } else { vq };
                return Ok(Edge {
                    id,
                    start,
                    kind: Kind::Arc { v },
                });
            }
        }
        let length = norm(sub(c, a));
        let direction = scale(sub(c, a), 1.0 / length);
        let id = b.edge(
            CurveGeometry::Line {
                origin: a,
                direction,
            },
            Interval::new(0.0, length).map_err(GeometryError::from)?,
            false,
        );
        // In a bend, a line runs along the axis at a known bend angle.
        let kind = match zone {
            Zone::Bend(i) => Kind::Axial {
                phi: self.phi(&self.folds[i], p),
                direction,
            },
            _ => Kind::Line { direction },
        };
        Ok(Edge {
            id,
            start: vp,
            kind,
        })
    }

    /// The rigid map that is valid at flat point p: for a bend, the parent's map
    /// turned to p's bend angle (its tangent plane there).
    fn rigid_at(&self, zone: Zone, p: P2) -> Rigid {
        match zone {
            Zone::Bend(i) => {
                let f = &self.folds[i];
                let axis = [
                    f.fold.edge_point[0] + f.fold.s0 * f.fold.n[0],
                    f.fold.edge_point[1] + f.fold.s0 * f.fold.n[1],
                    f.zc,
                ];
                let u = [f.fold.u[0], f.fold.u[1], 0.0];
                f.parent
                    .after(&Rigid::about(axis, u, -f.sigma * self.phi(f, p)))
            }
            _ => self.rigid(zone),
        }
    }

    fn surface(&self, zone: Zone, z: f64, top: bool) -> Surface {
        match zone {
            Zone::Bend(i) => {
                let f = &self.folds[i];
                // Outward is away from the axis on the outer layer only.
                let outer = (f.sigma > 0.0) != top;
                Surface::Cylinder {
                    frame: self.cylinder_frame(f),
                    radius: self.radius(f, z),
                    reversed: !outer,
                }
            }
            _ => {
                let m = self.rigid(zone);
                let s = if top { 1.0 } else { -1.0 };
                Surface::Plane(Frame3 {
                    origin: m.apply([0.0, 0.0, z]),
                    x: m.rot([1.0, 0.0, 0.0]),
                    y: m.rot([0.0, s, 0.0]),
                    z: m.rot([0.0, 0.0, s]),
                })
            }
        }
    }

    /// One use of `e` from vertex `from` to `to`, with its pcurve on `surface`.
    fn use_on(
        &self,
        b: &Builder,
        surface: &Surface,
        e: Edge,
        from: u32,
        to: u32,
    ) -> Result<Use, SheetError> {
        let sense = if from == e.start {
            Orientation::Forward
        } else {
            Orientation::Reverse
        };
        Ok(match surface {
            Surface::Plane(frame) => plane_boundary(b, *frame, e.id, from, to, sense)?,
            Surface::Cylinder { frame, .. } => {
                if let Kind::Cut { phi, .. } = e.kind {
                    let from_phi = if from == e.start { phi[0] } else { phi[1] };
                    let local = frame.to_local(b.brep().topology.vertices[from as usize].position);
                    return Ok(projected_boundary(
                        b,
                        &surface.geometry(),
                        e.id,
                        from,
                        to,
                        sense,
                        Some([from_phi, local[2]]),
                    )?);
                }
                let pcurve = match e.kind {
                    Kind::Arc { v } => uv_line([0.0, v], [1.0, 0.0]),
                    Kind::Helix { origin, direction } => uv_line(origin, direction),
                    Kind::Conic { origin, a, b } => PcurveGeometry::Conic2 {
                        origin,
                        axis_a: a,
                        axis_b: b,
                    },
                    Kind::Axial { phi, direction } => {
                        // A line along the axis at bend angle phi (taken from the
                        // flat, so angles past 180 degrees don't wrap).
                        let start = b.brep().topology.vertices[e.start as usize].position;
                        let local = frame.to_local(start);
                        uv_line([phi, local[2]], [0.0, dot(direction, frame.z)])
                    }
                    Kind::Line { direction } => {
                        // A line along the axis at a fixed angle.
                        let start = b.brep().topology.vertices[e.start as usize].position;
                        let local = frame.to_local(start);
                        uv_line(
                            [local[1].atan2(local[0]), local[2]],
                            [0.0, dot(direction, frame.z)],
                        )
                    }
                    Kind::Circle | Kind::Cut { .. } => unreachable!(),
                };
                boundary(e.id, from, to, sense, pcurve)
            }
        })
    }
}

#[derive(Clone, Copy)]
struct Edge {
    id: u32,
    /// Vertex at the start of the curve's parameter range.
    start: u32,
    kind: Kind,
}

#[derive(Clone, Copy)]
enum Kind {
    Line {
        direction: [f64; 3],
    },
    Circle,
    /// Line along a bend axis at bend angle phi.
    Axial {
        phi: f64,
        direction: [f64; 3],
    },
    /// Arc around a bend axis, at height v along the cylinder frame's z.
    Arc {
        v: f64,
    },
    /// A slanted flat line in a bend: a line in the bend faces' UV.
    Helix {
        origin: [f64; 2],
        direction: [f64; 2],
    },
    /// A flat circle in a bend: an ellipse in the bend faces' UV.
    Conic {
        origin: [f64; 2],
        a: [f64; 2],
        b: [f64; 2],
    },
    /// A mitered bend's end on one of its layers: an ellipse (parameter =
    /// bend angle - phase), its pcurve on the bend face projected.
    Cut {
        phi: [f64; 2],
        phase: f64,
    },
}

/// A flat edge: its end vertices (low, high) and, for an arc, a point on it.
type EdgeKey = (usize, usize, Option<(i64, i64)>);

/// A hole held whole by a cell.
struct HoleEdges {
    cell: usize,
    seam: u32,
    circles: [u32; 2],
    v: [u32; 2],
    wall: SurfaceGeometry,
    kind: Kind,
}

enum Surface {
    Plane(Frame3),
    Cylinder {
        frame: Frame3,
        radius: f64,
        reversed: bool,
    },
}

impl Surface {
    fn geometry(&self) -> SurfaceGeometry {
        match self {
            Surface::Plane(frame) => SurfaceGeometry::Plane { frame: *frame },
            Surface::Cylinder { frame, radius, .. } => SurfaceGeometry::Cylinder {
                frame: *frame,
                radius: *radius,
            },
        }
    }

    /// UV box around the loop vertices (and hole circles), padded.
    fn bounds(&self, b: &Builder, vertices: impl Iterator<Item = u32>) -> [[f64; 2]; 2] {
        let pos = |v: u32| b.brep().topology.vertices[v as usize].position;
        let mut uv: Vec<[f64; 2]> = Vec::new();
        match self {
            Surface::Plane(frame) => {
                for v in vertices {
                    let l = frame.to_local(pos(v));
                    uv.push([l[0], l[1]]);
                }
            }
            Surface::Cylinder { frame, .. } => {
                for v in vertices {
                    let l = frame.to_local(pos(v));
                    uv.push([l[1].atan2(l[0]), l[2]]);
                }
            }
        }
        let pad = 1e-3;
        let lo = |j: usize| uv.iter().map(|p| p[j]).fold(f64::INFINITY, f64::min) - pad;
        let hi = |j: usize| uv.iter().map(|p| p[j]).fold(f64::NEG_INFINITY, f64::max) + pad;
        [[lo(0), hi(0)], [lo(1), hi(1)]]
    }

    /// Material on the side opposite the surface normal (the inner layer of a bend).
    fn is_reversed(&self) -> bool {
        matches!(self, Surface::Cylinder { reversed: true, .. })
    }
}

fn add(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: P3, s: f64) -> P3 {
    a.map(|x| x * s)
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
fn norm(a: P3) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: P3) -> P3 {
    scale(a, 1.0 / norm(a))
}
fn sub2(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn dot2(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn norm2(a: P2) -> f64 {
    dot2(a, a).sqrt()
}
fn unit2(a: P2) -> P2 {
    let l = dot2(a, a).sqrt();
    [a[0] / l, a[1] / l]
}

/// Volume, surface area and centre of mass of the folded part, integrated
/// exactly over the flat layout (Green's theorem along the cell boundaries,
/// Gauss-Legendre quadrature on pieces short enough to be exact to rounding).
/// No mesh: use this where `mass` would need a very fine tessellation.
pub fn exact_mass(spec: &PartSpec) -> Result<Mass, SheetError> {
    let spec = &*spec.resolved()?;
    let layout = layout(spec)?;
    let model = Model::new(spec, &layout);
    model.refuse_clashes()?;
    Ok(model.exact_mass())
}

/// 8-point Gauss-Legendre on [0, 1]: (node, weight).
pub(crate) const GAUSS: [(f64, f64); 8] = {
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let mut g = [(0.0, 0.0); 8];
    let mut i = 0;
    while i < 4 {
        g[2 * i] = ((1.0 - X[i]) / 2.0, W[i] / 2.0);
        g[2 * i + 1] = ((1.0 + X[i]) / 2.0, W[i] / 2.0);
        i += 1;
    }
    g
};

/// A boundary piece of a flat cell: a line, or an arc (centre, radius, start
/// angle, signed sweep).
#[derive(Clone, Copy)]
enum Piece {
    Line(P2, P2),
    Arc(P2, f64, f64, f64),
    /// A mitered bend's end from flat s to s.
    Cut(Cut, f64, f64),
}

impl Piece {
    /// Point and derivative at s in [0, 1].
    fn at(&self, s: f64) -> (P2, P2) {
        match *self {
            Piece::Line(p, q) => (
                [p[0] + s * (q[0] - p[0]), p[1] + s * (q[1] - p[1])],
                sub2(q, p),
            ),
            Piece::Arc(c, r, a0, sweep) => {
                let (sn, cs) = (a0 + s * sweep).sin_cos();
                (
                    [c[0] + r * cs, c[1] + r * sn],
                    [-r * sweep * sn, r * sweep * cs],
                )
            }
            Piece::Cut(c, a, b) => {
                let x = a + s * (b - a);
                let (sn, cs) = c.phi(x).sin_cos();
                let dt = (-c.c[1] * sn + c.c[2] * cs) / c.radius;
                (
                    c.point(x),
                    [(c.u[0] * dt + c.n[0]) * (b - a), (c.u[1] * dt + c.n[1]) * (b - a)],
                )
            }
        }
    }

    /// Split so no piece turns more than pi/8 or moves more than `x_step`
    /// along `x_dir`.
    fn split(&self, x_dir: P2, x_step: f64) -> Vec<Piece> {
        let (turn, along) = match *self {
            Piece::Line(p, q) => (0.0, dot2(sub2(q, p), x_dir).abs()),
            Piece::Arc(_, r, _, sweep) => (sweep.abs(), r * sweep.abs()),
            Piece::Cut(c, a, b) => ((b - a).abs() / c.radius, (b - a).abs()),
        };
        let n = ((turn / (std::f64::consts::PI / 8.0))
            .max(along / x_step)
            .ceil() as usize)
            .max(1);
        (0..n)
            .map(|i| {
                let (s0, s1) = (i as f64 / n as f64, (i + 1) as f64 / n as f64);
                match *self {
                    Piece::Line(_, _) => Piece::Line(self.at(s0).0, self.at(s1).0),
                    Piece::Arc(c, r, a0, sweep) => {
                        Piece::Arc(c, r, a0 + s0 * sweep, sweep / n as f64)
                    }
                    Piece::Cut(c, a, b) => Piece::Cut(c, a + s0 * (b - a), a + s1 * (b - a)),
                }
            })
            .collect()
    }
}

/// A boundary piece's identity: its end points and, for an arc, its middle.
type PieceKey = ((i64, i64), (i64, i64), Option<(i64, i64)>);

fn piece_key(piece: &Piece) -> PieceKey {
    let key = |p: P2| ((p[0] * 1e6).round() as i64, (p[1] * 1e6).round() as i64);
    let (p, q) = (key(piece.at(0.0).0), key(piece.at(1.0).0));
    let mid = matches!(piece, Piece::Arc(..) | Piece::Cut(..)).then(|| key(piece.at(0.5).0));
    (p.min(q), p.max(q), mid)
}

impl Model<'_> {
    fn exact_mass(&self) -> Mass {
        let t = self.t;
        // Every boundary piece with its cell and whether it is a hole held
        // whole; a piece two cells share is no wall.
        let mut pieces: Vec<(usize, Piece, bool)> = Vec::new();
        let mut count: HashMap<PieceKey, usize> = HashMap::new();
        for (c, cell) in self.layout.cells.iter().enumerate() {
            let rings = std::iter::once((&cell.polygon, &cell.arcs))
                .chain(cell.inner.iter().map(|r| (&r.points, &r.arcs)));
            for (polygon, arcs) in rings {
                let n = polygon.len();
                for k in 0..n {
                    let (p, q) = (polygon[k], polygon[(k + 1) % n]);
                    let piece = match arcs[k] {
                        None => Piece::Line(p, q),
                        Some(Curve::Arc(arc)) => {
                            let r = dot2(sub2(p, arc.center), sub2(p, arc.center)).sqrt();
                            let a0 = (p[1] - arc.center[1]).atan2(p[0] - arc.center[0]);
                            Piece::Arc(arc.center, r, a0, arc.sweep(p, q))
                        }
                        Some(Curve::Cut(c)) => Piece::Cut(c, c.s_of(p), c.s_of(q)),
                    };
                    *count.entry(piece_key(&piece)).or_default() += 1;
                    pieces.push((c, piece, false));
                }
            }
            for h in &cell.holes {
                // Clockwise: the hole is not material.
                pieces.push((c, Piece::Arc(h.center, h.diameter / 2.0, 0.0, -TAU), true));
            }
        }
        let is_wall = |piece: &Piece, hole: bool| hole || count[&piece_key(piece)] == 1;

        let (mut volume, mut area, mut moment) = (0.0, 0.0, [0.0; 3]);
        for (c, cell) in self.layout.cells.iter().enumerate() {
            let mine = pieces.iter().filter(|(owner, _, _)| *owner == c);
            match cell.zone {
                Zone::Bend(i) => {
                    let f = &self.folds[i];
                    let rn = f.fold.radius + f.fold.k_factor * t;
                    let (a, b) = (self.radius(f, 0.0), self.radius(f, t));
                    let (r0, r1) = (a.min(b), a.max(b));
                    // Wrap flat coordinates: x around the axis, y along it.
                    let to = |p: P2| {
                        let w = self.to_wrap(f, p, 0.0);
                        [w[0], w[1]]
                    };
                    let to_dir = |v: P2| {
                        let w = self.to_wrap_dir(f, v);
                        [w[0], w[1]]
                    };
                    // Integrals over the cell of 1, y, cos(x / rn), sin(x / rn).
                    let mut s = [0.0; 4];
                    for (_, piece, hole) in mine {
                        let wall = is_wall(piece, *hole);
                        for part in piece.split(f.fold.n, rn * std::f64::consts::PI / 8.0) {
                            for (node, w) in GAUSS {
                                let (p, d) = part.at(node);
                                let (x, dx) = (to(p), to_dir(d));
                                s[0] += x[0] * dx[1] * w;
                                s[1] -= x[1] * x[1] / 2.0 * dx[0] * w;
                                s[2] += rn * (x[0] / rn).sin() * dx[1] * w;
                                s[3] -= rn * (x[0] / rn).cos() * dx[1] * w;
                                if wall {
                                    // |dW/ds| over the radius.
                                    for (rho_node, rw) in GAUSS {
                                        let rho = r0 + (r1 - r0) * rho_node;
                                        let speed =
                                            ((rho / rn * dx[0]).powi(2) + dx[1] * dx[1]).sqrt();
                                        area += speed * (r1 - r0) * w * rw;
                                    }
                                }
                            }
                        }
                    }
                    // The layout -> wrap map turns areas over for up bends.
                    let s = s.map(|v| -f.sigma * v);
                    let m1 = (r1 * r1 - r0 * r0) / (2.0 * rn);
                    let m2 = (r1.powi(3) - r0.powi(3)) / (3.0 * rn);
                    let v = s[0] * m1;
                    volume += v;
                    area += s[0] * (r0 + r1) / rn;
                    let frame = self.wrap(f).frame;
                    let local = frame.to_world([m2 * s[2], m2 * s[3], m1 * s[1]]);
                    for j in 0..3 {
                        moment[j] += frame.origin[j] * v + local[j];
                    }
                }
                zone => {
                    // Integrals over the cell of 1, x, y.
                    let mut s = [0.0; 3];
                    for (_, piece, hole) in mine {
                        let wall = is_wall(piece, *hole);
                        for part in piece.split([1.0, 0.0], f64::INFINITY) {
                            for (node, w) in GAUSS {
                                let (p, d) = part.at(node);
                                s[0] += p[0] * d[1] * w;
                                s[1] += p[0] * p[0] / 2.0 * d[1] * w;
                                s[2] -= p[1] * p[1] / 2.0 * d[0] * w;
                                if wall {
                                    area += (d[0] * d[0] + d[1] * d[1]).sqrt() * t * w;
                                }
                            }
                        }
                    }
                    let v = s[0] * t;
                    volume += v;
                    area += 2.0 * s[0];
                    let m = self.rigid(zone);
                    let local = m.rot([s[1] * t, s[2] * t, s[0] * t * t / 2.0]);
                    for j in 0..3 {
                        moment[j] += m.t[j] * v + local[j];
                    }
                }
            }
        }
        Mass {
            volume,
            area,
            center_of_mass: moment.map(|m| m / volume),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::flat::Zone;

    fn part(json: &str) -> PartSpec {
        serde_json::from_str(json).unwrap()
    }

    fn base(flanges: &str, holes: &str) -> String {
        format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "holes": [{holes}] }},
                "flanges": [{flanges}] }}"#
        )
    }

    const F1: &str = r#"{ "id": "f1", "on": { "base": 2 }, "angle": 90, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }"#;

    /// Volume from the layout: flat cells times T; a bend cell's flat area
    /// times the radial integral of r / (neutral radius) over the sheet.
    fn analytic_volume(spec: &PartSpec) -> f64 {
        let l = layout(spec).unwrap();
        let t = spec.thickness;
        l.cells
            .iter()
            .map(|c| {
                let area = crate::flat::path_area(&c.polygon, &c.arcs)
                    - c.holes
                        .iter()
                        .map(|h| std::f64::consts::PI * h.diameter * h.diameter / 4.0)
                        .sum::<f64>();
                match c.zone {
                    Zone::Bend(i) => {
                        let (r, k) = (l.folds[i].radius, l.folds[i].k_factor);
                        area * ((r + t).powi(2) - r * r) / (2.0 * (r + k * t))
                    }
                    _ => area * t,
                }
            })
            .sum()
    }

    #[test]
    fn every_part_type_folds_into_a_valid_solid_with_the_right_volume() {
        let lip = r#"{ "id": "lip", "on": { "tip": "f1" }, "angle": 90, "length": 15,
                       "datum": "outer", "position": "outside", "side": "down" }"#;
        let partial = r#"{ "id": "f1", "on": { "base": 2 }, "width": { "offset_start": 30, "offset_end": 20 },
                           "angle": 135, "length": 30, "datum": "outer", "position": "inside", "side": "up" }"#;
        let down = r#"{ "id": "f2", "on": { "base": 0 }, "angle": 90, "length": 30,
                        "datum": "outer", "position": "outside", "side": "down" }"#;
        for (name, json) in [
            ("plate", base("", "")),
            (
                "holes",
                base(
                    "",
                    r#"{ "center": [20, 30], "diameter": 8 }, { "center": [75, 20], "diameter": 12.5 }"#,
                ),
            ),
            ("flange", base(F1, "")),
            ("z-lip", base(&format!("{F1}, {lip}"), "")),
            (
                "partial 135 inside",
                base(partial, r#"{ "center": [20, 20], "diameter": 8 }"#),
            ),
            ("channel up/down", base(&format!("{F1}, {down}"), "")),
            (
                "holes across, inside and past the bend",
                base(
                    F1,
                    r#"{ "center": [50, 59], "diameter": 6 }, { "center": [25, 60.2], "diameter": 3 },
                       { "center": [80, 60.2], "diameter": 8 }, { "center": [70, 75], "diameter": 8 }"#,
                ),
            ),
            (
                "holes in a down bend and a lip",
                base(
                    &format!("{F1}, {lip}, {down}"),
                    r#"{ "center": [50, 1], "diameter": 5 }, { "center": [50, 89], "diameter": 4 }"#,
                ),
            ),
        ] {
            let spec = part(&json);
            let brep = fold(&spec).unwrap_or_else(|e| panic!("{name}: {e}"));
            let m = mass(&brep, 1e-4).unwrap();
            let expected = analytic_volume(&spec);
            assert!(
                ((m.volume - expected) / expected).abs() < 1e-5,
                "{name}: mesh {} vs analytic {expected}",
                m.volume
            );
        }
    }

    #[test]
    fn a_90_degree_flange_stands_on_its_edge() {
        // Expected: y 0..62, z 0..30 for this flange.
        let brep = fold(&part(&base(F1, ""))).unwrap();
        // Extents from the mesh: the B-rep's own bounds are padded UV boxes.
        let mesh = tessellate(&brep, 1e-4, 10_000_000).unwrap();
        for (axis, lo, hi) in [(0, 0.0, 100.0), (1, 0.0, 62.0), (2, 0.0, 30.0)] {
            // Only vertices that triangles use: the mesh also carries unused positions.
            let values = mesh
                .indices
                .iter()
                .map(|&i| &mesh.positions[i as usize * 3 + axis]);
            let (min, max) = values.fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
            assert!(
                (min - lo).abs() < 1e-6 && (max - hi).abs() < 1e-6,
                "axis {axis}: {min} .. {max}"
            );
        }
    }
}

#[cfg(test)]
mod corner_tests {
    use super::*;

    #[test]
    fn a_corner_folds_into_a_valid_solid() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
                 "flanges": [
                   { "id": "f1", "on": { "base": 1 }, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" },
                   { "id": "f2", "on": { "base": 2 }, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" } ] }"#,
        )
        .unwrap();
        let brep = fold(&spec).unwrap();
        let m = mass(&brep, 1e-5).unwrap();
        // Expected: 21011.1872 mm3, 24 faces.
        assert!((m.volume - 21011.1872).abs() < 1e-2, "{}", m.volume);
        assert!((exact_mass(&spec).unwrap().volume - 21011.1872).abs() < 1e-4);
        assert_eq!(brep.topology.faces.len(), 24);
    }
}

#[cfg(test)]
mod relief_tests {
    use super::*;

    #[test]
    fn a_round_relief_folds_with_their_volume_and_faces() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "relief": { "shape": "round", "width": 2, "depth": 1, "minimum_remnant": 4 },
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
                 "flanges": [{ "id": "f1", "on": { "base": 2 }, "width": { "offset_start": 30, "offset_end": 20 },
                               "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }] }"#,
        )
        .unwrap();
        let brep = fold(&spec).unwrap();
        let m = mass(&brep, 1e-5).unwrap();
        // Expected: 14848.9557 mm3, 20 faces (4 cylinders: two bend, two relief).
        assert!((m.volume - 14848.9557).abs() < 1e-2, "{}", m.volume);
        assert!((exact_mass(&spec).unwrap().volume - 14848.9557).abs() < 1e-4);
        assert_eq!(brep.topology.faces.len(), 20);
    }
}

#[cfg(test)]
mod hem_tests {
    use super::*;

    #[test]
    fn every_hem_kind_folds_with_their_volume_and_faces() {
        for (kind, volume, faces) in [
            (
                r#"{ "single": { "gap": 1, "length": 10 } }"#,
                13942.4778,
                14,
            ),
            (
                r#"{ "double": { "gap": 1, "length": 10 } }"#,
                15984.9556,
                30,
            ),
            (
                r#"{ "teardrop": { "radius": 3, "angle": 270 } }"#,
                15269.9112,
                14,
            ),
            (
                r#"{ "rolled": { "radius": 5, "angle": 240 } }"#,
                15626.5482,
                10,
            ),
        ] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "hems": [{{ "id": "h", "on": {{ "base": 2 }}, "kind": {kind}, "side": "up" }}] }}"#
            ))
            .unwrap();
            let brep = fold(&spec).unwrap_or_else(|e| panic!("{kind}: {e}"));
            let m = mass(&brep, 1e-5).unwrap();
            // Expected values.
            assert!(
                (m.volume - volume).abs() < 1e-2,
                "{kind}: {} vs {volume}",
                m.volume
            );
            assert!(
                (exact_mass(&spec).unwrap().volume - volume).abs() < 1e-4,
                "{kind}"
            );
            assert_eq!(brep.topology.faces.len(), faces, "{kind}");
        }
    }
}

#[cfg(test)]
mod hole_tests {
    use super::*;

    /// A corner with holes across both bends: the mesh agrees with the exact
    /// mass properties to its chord error (this mesh used to fail on a
    /// zero-area sliver along the bend end).
    #[test]
    fn holes_across_bends_mesh_close_to_the_exact_mass() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                           "holes": [{ "center": [60, 59], "diameter": 5 }, { "center": [99, 30], "diameter": 5 }] },
                 "flanges": [
                   { "id": "f1", "on": { "base": 1 }, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" },
                   { "id": "f2", "on": { "base": 2 }, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" } ] }"#,
        )
        .unwrap();
        let m = mass(&fold(&spec).unwrap(), 1e-3).unwrap();
        let e = exact_mass(&spec).unwrap();
        assert!(((m.volume - e.volume) / e.volume).abs() < 2e-5);
        assert!(((m.area - e.area) / e.area).abs() < 2e-4);
        assert!(norm(sub(m.center_of_mass, e.center_of_mass)) < 1e-3);
    }

    #[test]
    fn a_wrap_carries_the_flat_exactly_like_the_bend_map() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
                 "flanges": [
                   { "id": "f1", "on": { "base": 2 }, "angle": 120, "length": 30, "datum": "outer", "position": "outside", "side": "up" },
                   { "id": "f2", "on": { "base": 0 }, "angle": 70, "length": 30, "datum": "outer", "position": "outside", "side": "down" } ] }"#,
        )
        .unwrap();
        let l = layout(&spec).unwrap();
        let model = Model::new(&spec, &l);
        for (i, f) in model.folds.iter().enumerate() {
            let w = model.wrap(f);
            for p in [[10.0, 59.0], [50.0, 61.0], [90.0, 1.0], [30.0, -2.0]] {
                for z in [0.0, 0.7, 2.0] {
                    let a = model.map(Zone::Bend(i), p, z);
                    let b = w.point(model.to_wrap(f, p, z));
                    assert!(norm(sub(a, b)) < 1e-12, "{a:?} vs {b:?}");
                }
            }
        }
    }
}

#[cfg(test)]
mod corner_relief_tests {
    use super::*;

    fn corner(corner: &str, relief_shape: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "relief": {{ "shape": "{relief_shape}", "width": 2, "depth": 1, "minimum_remnant": 4
                              {corner} }},
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [
                   {{ "id": "f1", "on": {{ "base": 1 }}, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }},
                   {{ "id": "f2", "on": {{ "base": 2 }}, "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }} ] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn corner_reliefs_fold_with_their_volume_flat_and_faces() {
        // Expected: volume, flat area, faces.
        for (shape, volume, flat, faces) in [
            ("round", 20973.7123, 10438.8458, 25),
            ("square", 20953.5624, 10428.9039, 28),
            ("linear_weld", 21016.7081, 10459.1612, 22),
            ("tear", 21016.7081, 10461.3944, 22),
        ] {
            let spec = corner(
                &format!(r#", "corner": {{ "shape": "{shape}" }}"#),
                "straight",
            );
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{shape}: {}", m.volume);
            let area = crate::flat::flat_pattern(&spec).unwrap().area();
            assert!((area - flat).abs() < 1e-3, "{shape}: flat {area}");
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{shape}");
        }
    }

    #[test]
    fn a_small_corner_relief_sits_on_the_bisector_through_q_capped_at_the_flanges_first_corners() {
        // Expected: q = (98, 58); the flanges'
        // first corners (102.398, 58.586) and (98.586, 62.398).
        let spec = |shape: &str, size: f64| {
            corner(
                &format!(r#", "corner": {{ "shape": "{shape}", "size": {size} }}"#),
                "straight",
            )
        };
        let centre = |size: f64| {
            let flat = crate::flat::flat_pattern(&spec("round", size)).unwrap();
            flat.outline_arcs
                .iter()
                .flatten()
                .find_map(|c| c.arc())
                .unwrap()
                .center
        };
        // Not capped: the circle runs through q.
        let d = 1.5 / 2f64.sqrt();
        let c = centre(3.0);
        assert!((c[0] - 98.0 - d).abs() < 1e-12 && (c[1] - 58.0 - d).abs() < 1e-12);
        // Capped: the centre of the circle through q and both first corners.
        let c = centre(5.6);
        assert!((c[0] - 99.975071).abs() < 1e-6 && (c[1] - 59.975071).abs() < 1e-6);
        // A relief that isn't capped folds as trim to bend (its flat
        // 10457.7190 mm2 against 10456.5847); a capped one is cut.
        for (shape, size, volume, flat) in [
            ("round", 3.0, 21011.1872, 10457.7190),
            ("round", 5.6, 21002.0051, 10452.2998),
            ("square", 5.5, 20993.6001, 10448.2781),
        ] {
            let m = exact_mass(&spec(shape, size)).unwrap();
            assert!(
                (m.volume - volume).abs() < 1e-3,
                "{shape} {size}: {}",
                m.volume
            );
            let area = crate::flat::flat_pattern(&spec(shape, size))
                .unwrap()
                .area();
            assert!((area - flat).abs() < 1e-3, "{shape} {size}: {area}");
        }
        // Too small to reach the flanges' end lines: centred where they cross,
        // (e, e) past q with e = R - gap / sqrt 2.
        let e = 2.0 - 2.0 / 2f64.sqrt();
        let c = centre(0.6);
        assert!((c[0] - 98.0 - e).abs() < 1e-9 && (c[1] - 58.0 - e).abs() < 1e-9);
        // Refused: a capped square short of the bend's far side.
        assert!(matches!(
            fold(&spec("square", 4.5)),
            Err(SheetError::Unsupported(_))
        ));
    }

    #[test]
    fn small_corner_reliefs_at_sloped_flanges_follow_the_flanges_end_lines() {
        let spec = |angle: f64, size: f64| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "relief": {{ "width": 2, "depth": 1, "minimum_remnant": 4,
                                  "corner": {{ "shape": "round", "size": {size} }} }},
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "flanges": [
                       {{ "id": "f1", "on": {{ "base": 1 }}, "angle": {angle}, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }},
                       {{ "id": "f2", "on": {{ "base": 2 }}, "angle": {angle}, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }}] }}"#
            ))
            .unwrap()
        };
        let centre = |angle: f64, size: f64| {
            let flat = crate::flat::flat_pattern(&spec(angle, size)).unwrap();
            flat.outline_arcs
                .iter()
                .flatten()
                .find_map(|c| c.arc())
                .unwrap()
                .center
        };
        // 135 degrees: on the bisector, its circle through q (95.171573,
        // 55.171573), reaching the flanges' end lines.
        let c = centre(135.0, 4.0);
        assert!(
            (c[0] - 96.585786).abs() < 1e-6 && (c[1] - 56.585786).abs() < 1e-6,
            "{c:?}"
        );
        // 45 degrees: the end lines cross in the base, behind q; a relief
        // that doesn't reach them is centred there, whatever its size, and cut
        // from the flat only.
        for size in [1.5, 2.8] {
            let c = centre(45.0, size);
            assert!(
                (c[0] - 93.86244).abs() < 1e-5 && (c[1] - 53.86244).abs() < 1e-5,
                "{c:?}"
            );
            let m = exact_mass(&spec(45.0, size)).unwrap();
            assert!((m.volume - 22594.42182753766).abs() < 1e-3, "{}", m.volume);
        }
    }

    #[test]
    fn a_tear_relief_is_a_hairline_notch_in_3d_and_none_in_the_flat() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "relief": { "shape": "tear", "width": 2, "depth": 1, "minimum_remnant": 4 },
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
                 "flanges": [{ "id": "f1", "on": { "base": 2 }, "width": { "offset_start": 30, "offset_end": 20 },
                               "angle": 90, "length": 30, "datum": "outer", "position": "outside", "side": "up" }] }"#,
        )
        .unwrap();
        // Expected: 14871.1585 mm3, flat 7419.9115 mm2 (plate plus flange), 22 faces.
        let m = exact_mass(&spec).unwrap();
        assert!((m.volume - 14871.1585).abs() < 1e-3, "{}", m.volume);
        let area = crate::flat::flat_pattern(&spec).unwrap().area();
        assert!((area - 7419.9115).abs() < 1e-3, "{area}");
        assert_eq!(fold(&spec).unwrap().topology.faces.len(), 22);
    }
}

#[cfg(test)]
mod cutout_tests {
    use super::*;

    fn part(cutouts: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "cutouts": [{cutouts}] }},
                 "flanges": [{{ "id": "f1", "on": {{ "base": 2 }}, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }}] }}"#
        ))
        .unwrap()
    }

    /// The flange alone.
    const FLANGE: f64 = 17742.4778;

    #[test]
    fn cut_outs_take_their_flat_area_times_the_thickness_or_the_bend_factor() {
        for (cutout, removed) in [
            // A 20 x 15 rectangle on the base: 300 mm2 x T.
            (r#"{ "outline": [[30,20],[50,20],[50,35],[30,35]] }"#, 600.0),
            // A triangle on the flange: 75 mm2 x T.
            (r#"{ "outline": [[40,70],[55,70],[45,80]] }"#, 150.0),
            // 7 x 2.8 inside the bend: 19.6 mm2 x ((R+T)^2 - R^2) / (2 (R + kT)).
            (
                r#"{ "outline": [[85,58.8],[92,58.8],[92,61.6],[85,61.6]] }"#,
                42.0,
            ),
        ] {
            let spec = part(cutout);
            let m = exact_mass(&spec).unwrap();
            assert!(
                (FLANGE - m.volume - removed).abs() < 1e-3,
                "{cutout}: {}",
                m.volume
            );
            let mesh = mass(&fold(&spec).unwrap(), 1e-3).unwrap();
            assert!(
                ((mesh.volume - m.volume) / m.volume).abs() < 2e-5,
                "{cutout}"
            );
        }
    }

    #[test]
    fn slots_rounded_and_non_convex_cut_outs_cross_bends() {
        let k = (std::f64::consts::PI / 8.0).tan();
        for cutout in [
            // A slot across the whole bend.
            r#"{ "outline": [[67.5,54],[72.5,54],[72.5,66],[67.5,66]], "bulges": [1,0,1,0] }"#
                .to_string(),
            // An L across the bend.
            r#"{ "outline": [[20,50],[35,50],[35,56],[26,56],[26,64],[20,64]] }"#.to_string(),
            // Rounded corners whose circles touch the bend start off the arcs.
            format!(
                r#"{{ "outline": [[82,54],[90,54],[92,56],[92,60],[90,62],[82,62],[80,60],[80,56]],
                     "bulges": [0,{k},0,{k},0,{k},0,{k}] }}"#
            ),
        ] {
            let spec = part(&cutout);
            let m = exact_mass(&spec).unwrap();
            let mesh = mass(&fold(&spec).unwrap(), 1e-3).unwrap();
            assert!(
                ((mesh.volume - m.volume) / m.volume).abs() < 2e-5,
                "{cutout}"
            );
        }
    }

    #[test]
    fn a_cut_out_across_the_outline_is_a_notch() {
        // 20 x 10 of it lies on the part: 200 mm2 x T.
        let spec = part(r#"{ "outline": [[30,-5],[50,-5],[50,10],[30,10]] }"#);
        let m = exact_mass(&spec).unwrap();
        assert!((FLANGE - m.volume - 400.0).abs() < 1e-3, "{}", m.volume);
        let flat = crate::flat::flat_pattern(&spec).unwrap();
        assert!(flat.cutouts.is_empty(), "a notch is part of the outline");
        fold(&spec).unwrap();
        // Wholly outside: wrong.
        let spec = part(r#"{ "outline": [[30,-25],[50,-25],[50,-10],[30,-10]] }"#);
        assert!(matches!(fold(&spec), Err(SheetError::Invalid(_))));
    }
}

#[cfg(test)]
mod corner_angle_tests {
    use super::*;
    use crate::spec::{CornerRelief, CornerShape, Relief, ReliefShape};

    fn part(outline: &str, edges: &[usize], angle: f64) -> PartSpec {
        let flanges: Vec<String> = edges
            .iter()
            .map(|e| {
                format!(
                    r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, "angle": {angle}, "length": 30,
                         "datum": "outer", "position": "outside", "side": "up" }}"#
                )
            })
            .collect();
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": {outline} }}, "flanges": [{}] }}"#,
            flanges.join(",")
        ))
        .unwrap()
    }

    #[test]
    fn corner_reliefs_at_any_angle_fold_with_their_volume_and_faces() {
        let rect = "[[0,0],[100,0],[100,60],[0,60]]";
        let trapezoid = "[[0,0],[100,0],[100,60],[65.358983848622,80],[0,80]]";
        let relieved = |outline: &str, angle: f64, shape: CornerShape, size: Option<f64>| {
            let mut spec = part(outline, &[1, 2], angle);
            spec.relief = Some(Relief {
                shape: ReliefShape::Straight,
                width: 2.0,
                depth: 1.0,
                minimum_remnant: 4.0,
                miter_gap: None,
                corner: Some(CornerRelief { shape, size }),
                seam: None,
                transition: None,
            });
            spec
        };
        use CornerShape::*;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "round, 45",
                relieved(rect, 45.0, Round, None),
                22525.4758,
                25,
            ),
            (
                "round, 135",
                relieved(rect, 135.0, Round, Some(10.0)),
                18121.0757,
                25,
            ),
            (
                "round, 120 corner",
                relieved(trapezoid, 90.0, Round, Some(10.0)),
                20817.1412,
                26,
            ),
            (
                "square, 45",
                relieved(rect, 45.0, Square, None),
                22504.7132,
                28,
            ),
            // Turned with the corner: its diagonal on the normals' bisector.
            (
                "square, 120 corner",
                relieved(trapezoid, 90.0, Square, Some(10.0)),
                20784.5204,
                29,
            ),
            (
                "weld, 135",
                relieved(rect, 135.0, LinearWeld, None),
                18158.5987,
                22,
            ),
            (
                "tear, 120 corner",
                relieved(trapezoid, 90.0, Tear, None),
                20921.8178,
                23,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn corners_at_any_angle_fold_with_their_volume_and_faces() {
        let rect = "[[0,0],[100,0],[100,60],[0,60]]";
        let trapezoid = "[[0,0],[100,0],[100,60],[65.358983848622,80],[0,80]]";
        let hexagon = format!(
            "[{}]",
            (0..6)
                .map(|i| {
                    let a = std::f64::consts::PI / 3.0 * i as f64;
                    format!("[{},{}]", 50.0 + 40.0 * a.cos(), 50.0 + 40.0 * a.sin())
                })
                .collect::<Vec<_>>()
                .join(",")
        );
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("45", part(rect, &[1, 2], 45.0), 22594.4218, 22),
            ("135", part(rect, &[1, 2], 135.0), 18158.5987, 22),
            ("120", part(trapezoid, &[1, 2], 90.0), 20921.8178, 23),
            // The flange starts short of the bend: the bends stop with it.
            ("120 / 45", part(trapezoid, &[1, 2], 45.0), 21840.0027, 25),
            (
                "60",
                part(
                    "[[65.358983848622,0],[100,0],[100,60],[48.038475772934,30]]",
                    &[1, 2],
                    90.0,
                ),
                10788.4549,
                24,
            ),
            (
                "hexagonal tray",
                part(&hexagon, &[0, 1, 2, 3, 4, 5], 90.0),
                21272.3767,
                56,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod mixed_corner_tests {
    use super::*;

    fn flange(id: &str, edge: usize, length: f64, side: &str) -> String {
        format!(
            r#"{{ "id": "{id}", "on": {{ "base": {edge} }}, "angle": 90, "length": {length},
                 "datum": "outer", "position": "outside", "side": "{side}" }}"#
        )
    }

    fn part(flanges: &[String]) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }}, "flanges": [{}] }}"#,
            flanges.join(",")
        ))
        .unwrap()
    }

    #[test]
    fn flanges_with_different_settings_meet_as_separate_features() {
        // Expected: the flange built first keeps its full edge.
        for (name, spec, volume, faces) in [
            (
                "lengths",
                part(&[flange("f1", 1, 30.0, "up"), flange("f2", 2, 20.0, "up")]),
                19094.2654,
                23,
            ),
            (
                "the other order",
                part(&[flange("f2", 2, 20.0, "up"), flange("f1", 1, 30.0, "up")]),
                19054.2654,
                23,
            ),
            (
                "up and down",
                part(&[flange("f1", 1, 30.0, "up"), flange("f2", 2, 30.0, "down")]),
                21054.2654,
                23,
            ),
            (
                "a tray",
                part(&[
                    flange("f0", 0, 30.0, "up"),
                    flange("f1", 1, 20.0, "up"),
                    flange("f2", 2, 30.0, "up"),
                    flange("f3", 3, 20.0, "up"),
                ]),
                27561.1325,
                42,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod corner_lip_tests {
    use super::*;

    fn flange(id: &str, on: &str, length: f64) -> String {
        format!(
            r#"{{ "id": "{id}", "on": {on}, "angle": 90, "length": {length},
                 "datum": "outer", "position": "outside", "side": "up" }}"#
        )
    }

    fn part(flanges: &[String]) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }}, "flanges": [{}] }}"#,
            flanges.join(",")
        ))
        .unwrap()
    }

    #[test]
    fn a_lip_continues_its_flanges_corner_edge() {
        let base = |e: usize| flange(&format!("f{e}"), &format!(r#"{{ "base": {e} }}"#), 30.0);
        let lip = |id: &str, on: &str, length: f64| {
            flange(id, &format!(r#"{{ "tip": "{on}" }}"#), length)
        };
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "corner",
                part(&[base(1), base(2), lip("l2", "f2", 15.0)]),
                23714.8805,
                32,
            ),
            (
                "tray",
                part(&[
                    base(0),
                    base(1),
                    base(2),
                    base(3),
                    lip("l0", "f0", 12.0),
                    lip("l2", "f2", 12.0),
                ]),
                33832.5784,
                62,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod corner_hem_tests {
    use super::*;

    fn part(flanges: &[usize], hems: &[(&str, &str, &str)]) -> PartSpec {
        let flanges: Vec<String> = flanges
            .iter()
            .map(|e| {
                format!(
                    r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, "angle": 90, "length": 30,
                         "datum": "outer", "position": "outside", "side": "up" }}"#
                )
            })
            .collect();
        let hems: Vec<String> = hems
            .iter()
            .map(|(id, on, side)| {
                format!(
                    r#"{{ "id": "{id}", "on": {on}, "side": "{side}",
                         "kind": {{ "single": {{ "gap": 1, "length": 8 }} }} }}"#
                )
            })
            .collect();
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [{}], "hems": [{}] }}"#,
            flanges.join(","),
            hems.join(",")
        ))
        .unwrap()
    }

    #[test]
    fn hems_at_corners_fold() {
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            // On a flange's tip, like a lip.
            (
                "on a tip",
                part(&[1, 2], &[("h2", r#"{ "tip": "f2" }"#, "up")]),
                22531.851,
                32,
            ),
            // Next to a flange: built after it, so the flange keeps its edge.
            (
                "next to a flange",
                part(&[1], &[("h2", r#"{ "base": 2 }"#, "up")]),
                16933.553,
                23,
            ),
            (
                "down, before the flange's edge",
                part(&[2], &[("h1", r#"{ "base": 1 }"#, "down")]),
                18613.553,
                23,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn hems_that_meet_are_refused() {
        // These fold, but don't unfold.
        for spec in [
            part(
                &[],
                &[
                    ("h1", r#"{ "base": 1 }"#, "up"),
                    ("h2", r#"{ "base": 2 }"#, "up"),
                ],
            ),
            part(
                &[1, 2],
                &[
                    ("h1", r#"{ "tip": "f1" }"#, "up"),
                    ("h2", r#"{ "tip": "f2" }"#, "up"),
                ],
            ),
        ] {
            assert!(matches!(fold(&spec), Err(SheetError::Unsupported(_))));
        }
    }
}

#[cfg(test)]
mod contour_tests {
    use super::*;
    use crate::spec::{Relief, ReliefShape, Width};

    #[test]
    fn a_contour_along_two_edges_is_mitered_at_their_corner() {
        let along = |profile: &str, gap: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4 {gap},
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "contours": [{{ "id": "c", "on": {{ "base": 1 }}, "along": [2],
                                    "profile": {profile}, "face": "bottom" }}] }}"#
            ))
            .unwrap()
        };
        // An L (up 30, back 15 over the base): every leg ends where the first
        // of its faces is half the miter gap from the corner's symmetry plane;
        // across its second bend that is a curve. Expected: flat area, the
        // return's mitered end at the corner, volume.
        let l = along("[[0,2],[0,32],[-15,32]]", "");
        let flat = crate::flat::flat_pattern(&l).unwrap();
        assert!(
            (flat.area() - 11907.826811774863).abs() < 1e-6,
            "{}",
            flat.area()
        );
        assert!(
            flat.outline
                .iter()
                .any(|p| dist2(*p, [94.585786, 88.796459]) < 1e-5),
            "{:?}",
            flat.outline
        );
        assert!(flat
            .outline_arcs
            .iter()
            .flatten()
            .any(|c| matches!(c, crate::flat::Curve::Cut(_))));
        let v = exact_mass(&l).unwrap().volume;
        assert!((v - 24006.485303712703).abs() < 1e-6, "{v}");
        // One bend: two inside flanges with the same settings.
        let wall = along("[[0,2],[0,32]]", "");
        let v = exact_mass(&wall).unwrap().volume;
        assert!((v - 20762.174331051156).abs() < 1e-6, "{v}");
        // Along edges that don't follow each other: refused.
        let mut skip = l.clone();
        skip.contours[0].along = vec![3];
        assert!(matches!(fold(&skip), Err(SheetError::Invalid(_))));
    }

    fn part(profile: &str, face: &str) -> PartSpec {
        arcs(profile, "[]", face)
    }

    fn arcs(profile: &str, bulges: &str, face: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "contours": [{{ "id": "c", "on": {{ "base": 2 }}, "profile": {profile},
                                "bulges": {bulges}, "face": "{face}" }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn contours_at_corners_follow_build_order() {
        let spec = |features: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }}, {features} }}"#
            ))
            .unwrap()
        };
        let flange = |e: usize| {
            format!(
                r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, "angle": 90, "length": 30,
                     "datum": "outer", "position": "outside", "side": "up" }}"#
            )
        };
        let contour = |id: &str, e: usize, profile: &str| {
            format!(
                r#"{{ "id": "{id}", "on": {{ "base": {e} }}, "profile": {profile}, "face": "bottom" }}"#
            )
        };
        let (l, z) = ("[[0,2],[0,32],[-15,32]]", "[[0,2],[0,22],[20,22],[20,42]]");
        let out = "[[0,2],[0,32],[15,32]]";
        let tray: Vec<String> = (0..4).map(|e| contour(&format!("c{e}"), e, out)).collect();
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "a flange, then a one-bend contour",
                spec(&format!(
                    r#""flanges": [{}], "contours": [{}]"#,
                    flange(1),
                    contour("c", 2, "[[0,2],[0,32]]")
                )),
                21035.4158,
                23,
            ),
            (
                "a flange, then an L",
                spec(&format!(
                    r#""flanges": [{}], "contours": [{}]"#,
                    flange(1),
                    contour("c", 2, l)
                )),
                23331.044,
                31,
            ),
            (
                "flanges on both sides of a Z",
                spec(&format!(
                    r#""flanges": [{}, {}], "contours": [{}]"#,
                    flange(1),
                    flange(3),
                    contour("c", 2, z)
                )),
                29593.9112,
                48,
            ),
            (
                "an L, then a hem",
                spec(&format!(
                    r#""contours": [{}], "hems": [{{ "id": "h", "on": {{ "base": 1 }},
                        "kind": {{ "single": {{ "gap": 1, "length": 8 }} }}, "side": "up" }}]"#,
                    contour("c", 2, l)
                )),
                20925.1812,
                31,
            ),
            (
                "an L, then a Z",
                spec(&format!(
                    r#""contours": [{}, {}]"#,
                    contour("c", 2, l),
                    contour("c1", 1, z)
                )),
                26334.6192,
                47,
            ),
            (
                "a tray of four contours",
                spec(&format!(r#""contours": [{}]"#, tray.join(","))),
                38859.4686,
                74,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn contours_sit_on_flange_tips_in_the_flanges_frame() {
        let spec = |angle: f64, side: &str, width: &str, profile: &str, face: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "flanges": [{{ "id": "f", "on": {{ "base": 2 }}, {width} "angle": {angle},
                                    "length": 30, "datum": "outer", "position": "outside",
                                    "side": "{side}" }}],
                     "contours": [{{ "id": "c", "on": {{ "tip": "f" }}, "profile": {profile},
                                     "face": "{face}" }}] }}"#
            ))
            .unwrap()
        };
        let step = "[[0,0],[0,15],[10,15]]";
        let partial = r#""width": { "offset_start": 20, "offset_end": 30 },"#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("up", spec(90.0, "up", "", step, "bottom"), 22227.4334, 30),
            (
                "down",
                spec(90.0, "down", "", step, "bottom"),
                22227.4334,
                30,
            ),
            (
                "45 degrees, top face",
                spec(45.0, "up", "", "[[0,2],[0,12],[15,12],[15,24]]", "top"),
                25001.616,
                38,
            ),
            (
                "on a partial flange",
                spec(90.0, "up", partial, "[[0,0],[0,-12]]", "top"),
                15918.4778,
                30,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn partial_contours_get_bend_reliefs_like_flanges() {
        let width = |mut spec: PartSpec, start: f64, end: f64| {
            spec.contours[0].width = Some(Width {
                offset_start: start,
                offset_end: end,
            });
            spec
        };
        let l = "[[0,2],[0,32],[-15,32]]";
        let mut round = width(part(l, "bottom"), 25.0, 25.0);
        round.relief = Some(Relief {
            shape: ReliefShape::Round,
            width: 4.0,
            depth: 3.0,
            minimum_remnant: 4.0,
            miter_gap: None,
            corner: None,
            seam: None,
            transition: None,
        });
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("L", width(part(l, "bottom"), 20.0, 30.0), 16002.4778, 30),
            (
                "Z",
                width(part("[[0,2],[0,22],[20,22],[20,42]]", "bottom"), 10.0, 10.0),
                20941.9467,
                38,
            ),
            (
                "arc, flush with a corner",
                width(
                    arcs(
                        "[[0,2],[0,27],[-5,32],[-15,32]]",
                        &format!("[0,{},0]", (std::f64::consts::PI / 8.0).tan()),
                        "top",
                    ),
                    0.0,
                    40.0,
                ),
                17404.46,
                26,
            ),
            ("round reliefs", round, 15937.3451, 30),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn arcs_in_a_contour_are_bends_with_their_radius() {
        let q = (std::f64::consts::PI / 8.0).tan();
        let l = "[[0,2],[0,27],[-5,32],[-15,32]]";
        let quarter = format!("[0,{q},0]");
        let wide = format!("[0,{},0]", (135f64.to_radians() / 4.0).tan());
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            // The arc is the outer face: inner radius 3.
            ("R5 corner", arcs(l, &quarter, "bottom"), 19999.1149, 22),
            (
                "R5 corner, top face",
                arcs(l, &quarter, "top"),
                21027.4334,
                22,
            ),
            (
                "135 degrees",
                arcs(
                    "[[0,2],[0,20],[-10.242640687119286,24.242640687119286],
                      [-17.31370849898476,17.17157287525381]]",
                    &wide,
                    "bottom",
                ),
                19698.6723,
                22,
            ),
            (
                "an arc and a sharp corner",
                arcs(
                    "[[0,2],[0,17],[5,22],[20,22],[20,42]]",
                    &format!("[0,{},0,0]", -q),
                    "bottom",
                ),
                22969.9112,
                30,
            ),
            (
                "starting with an arc",
                arcs("[[0,0],[10,10],[10,30]]", &format!("[{q},0]"), "bottom"),
                18827.4334,
                14,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn contour_flanges_fold_with_their_volume_and_faces() {
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "L",
                part("[[0,2],[0,32],[-15,32]]", "bottom"),
                20084.9556,
                22,
            ),
            (
                "L, top face",
                part("[[0,2],[0,32],[-15,32]]", "top"),
                21284.9556,
                22,
            ),
            (
                "L outward",
                part("[[0,2],[0,32],[15,32]]", "bottom"),
                20884.9556,
                22,
            ),
            (
                "L from z 0",
                part("[[0,0],[0,30],[-15,30]]", "bottom"),
                19684.9556,
                22,
            ),
            (
                "Z",
                part("[[0,2],[0,22],[20,22],[20,42]]", "bottom"),
                23227.4334,
                30,
            ),
            // A slanted line: its corner is where it meets the face.
            ("slanted", part("[[0,2],[10,22]]", "bottom"), 16394.7844, 14),
            (
                "slanted from z 0",
                part("[[0,0],[10,20]]", "bottom"),
                16147.5708,
                14,
            ),
            (
                "slanted, top face",
                part("[[0,2],[10,22]]", "top"),
                16641.9980,
                14,
            ),
            (
                "slanted from z 0, top",
                part("[[0,0],[10,20]]", "top"),
                16394.7844,
                14,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod rounded_tests {
    use super::*;

    const ROUNDED: &str = "[[10,0],[90,0],[100,10],[100,50],[90,60],[10,60],[0,50],[0,10]]";

    fn part(outline: &str, bulges: &str, features: &str, relief: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, {relief}
                 "base": {{ "outline": {outline}, "bulges": {bulges} }}, {features} }}"#
        ))
        .unwrap()
    }

    fn flanges(edges: &[usize], angle: f64, position: &str, side: &str) -> String {
        let f: Vec<String> = edges
            .iter()
            .map(|e| {
                format!(
                    r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, "angle": {angle}, "length": 30,
                         "datum": "outer", "position": "{position}", "side": "{side}" }}"#
                )
            })
            .collect();
        format!(r#""flanges": [{}]"#, f.join(","))
    }

    #[test]
    fn flanges_next_to_rounded_corners_get_reliefs_running_to_the_arc() {
        let q = (std::f64::consts::PI / 8.0).tan();
        let corners = format!("[0,{q},0,{q},0,{q},0,{q}]");
        let narrow = r#""relief": { "width": 1, "depth": 1, "minimum_remnant": 2 },"#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "remnant at the minimum: the notch runs on",
                part(ROUNDED, &corners, &flanges(&[4], 90.0, "outside", "up"), ""),
                16363.201,
                22,
            ),
            (
                "remnant above it: a notch of the relief's width",
                part(
                    ROUNDED,
                    &corners,
                    &flanges(&[4], 90.0, "inside", "down"),
                    "",
                ),
                16062.8373,
                24,
            ),
            (
                "narrow relief",
                part(
                    ROUNDED,
                    &corners,
                    &flanges(&[4], 90.0, "outside", "up"),
                    narrow,
                ),
                16410.3675,
                24,
            ),
            (
                "45 degrees",
                part(ROUNDED, &corners, &flanges(&[4], 45.0, "outside", "up"), ""),
                16579.0021,
                22,
            ),
            (
                "tray: notches crossing before the arc",
                part(
                    ROUNDED,
                    &corners,
                    &flanges(&[0, 2, 4, 6], 90.0, "outside", "up"),
                    "",
                ),
                25373.9467,
                54,
            ),
            (
                "tray, R20: both notches trim the arc",
                part(
                    "[[20,0],[80,0],[100,20],[100,40],[80,60],[20,60],[0,40],[0,20]]",
                    &corners,
                    &flanges(&[0, 2, 4, 6], 90.0, "outside", "up"),
                    "",
                ),
                20406.3069,
                66,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }

    #[test]
    fn a_flange_on_an_arc_edge_is_refused() {
        let q = (std::f64::consts::PI / 8.0).tan();
        let spec = part(
            ROUNDED,
            &format!("[0,{q},0,{q},0,{q},0,{q}]"),
            &flanges(&[1], 90.0, "outside", "up"),
            "",
        );
        assert!(matches!(exact_mass(&spec), Err(SheetError::Unsupported(_))));
    }
}

#[cfg(test)]
mod seam_tests {
    use super::*;

    fn part(edges: &[usize], position: &str, side: &str, relief: &str) -> PartSpec {
        let flanges: Vec<String> = edges
            .iter()
            .map(|e| {
                format!(
                    r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, "angle": 90, "length": 30,
                         "datum": "outer", "position": "{position}", "side": "{side}" }}"#
                )
            })
            .collect();
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "relief": {{ "width": 2, "depth": 1, "minimum_remnant": 4, {relief} }},
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [{}] }}"#,
            flanges.join(",")
        ))
        .unwrap()
    }

    #[test]
    fn overlapping_seams_fold_with_their_volume_and_faces() {
        let overlap = r#""seam": "overlap""#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "overlap",
                part(&[1, 2], "outside", "up", overlap),
                21158.2654,
                23,
            ),
            (
                "reverse",
                part(&[1, 2], "outside", "up", r#""seam": "reverse_overlap""#),
                21158.2654,
                23,
            ),
            (
                "gap 1",
                part(
                    &[1, 2],
                    "outside",
                    "up",
                    r#""seam": "overlap", "miter_gap": 1"#,
                ),
                21210.2654,
                24,
            ),
            (
                "inside",
                part(&[1, 2], "inside", "up", overlap),
                20296.5663,
                23,
            ),
            (
                "down",
                part(&[1, 2], "outside", "down", overlap),
                21158.2654,
                23,
            ),
            (
                "tray",
                part(&[0, 1, 2, 3], "outside", "up", overlap),
                30257.1325,
                42,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod tip_corner_tests {
    use super::*;

    fn part(angle: f64, position: &str, side: &str, width: &str, cut: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [{{ "id": "f", "on": {{ "base": 2 }}, {width} "angle": {angle},
                                "length": 30, "datum": "outer", "position": "{position}",
                                "side": "{side}", "tip_corners": {cut} }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn rounded_and_chamfered_tips_fold_with_their_volume_and_faces() {
        let round = r#"{ "round": { "radius": 5 } }"#;
        let partial = r#""width": { "offset_start": 20, "offset_end": 30 },"#;
        // Expected: volume, faces. Chamfer and round differ by
        // 2 x (12.5 - (25 - 25 pi / 4)) x T = 28.54.
        for (name, spec, volume, faces) in [
            (
                "round",
                part(90.0, "outside", "up", "", round),
                17721.0176,
                16,
            ),
            (
                "chamfer",
                part(
                    90.0,
                    "outside",
                    "up",
                    "",
                    r#"{ "chamfer": { "distance": 5 } }"#,
                ),
                17692.4778,
                16,
            ),
            (
                "round, partial",
                part(
                    90.0,
                    "outside",
                    "up",
                    partial,
                    r#"{ "round": { "radius": 8 } }"#,
                ),
                14792.3008,
                24,
            ),
            (
                "round, 45",
                part(45.0, "outside", "up", "", round),
                17952.7224,
                16,
            ),
            (
                "round, down",
                part(90.0, "inside", "down", "", round),
                17321.0176,
                16,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod partial_hem_tests {
    use super::*;

    fn part(kind: &str, side: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "hems": [{{ "id": "h", "on": {{ "base": 2 }},
                            "width": {{ "offset_start": 20, "offset_end": 30 }},
                            "kind": {kind}, "side": "{side}" }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn partial_hems_get_bend_reliefs_like_flanges() {
        let single = r#"{ "single": { "gap": 1, "length": 8 } }"#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("single", part(single, "up"), 12743.2389, 22),
            ("single, down", part(single, "down"), 12743.2389, 22),
            (
                "double",
                part(r#"{ "double": { "gap": 1, "length": 10 } }"#, "up"),
                13964.4778,
                38,
            ),
            (
                "teardrop",
                part(r#"{ "teardrop": { "radius": 2, "angle": 270 } }"#, "up"),
                13123.7167,
                22,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod own_settings_tests {
    use super::*;

    fn part(features: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }}, {features} }}"#
        ))
        .unwrap()
    }

    fn flange(edge: usize, own: &str) -> String {
        format!(
            r#"{{ "id": "f{edge}", "on": {{ "base": {edge} }}, "angle": 90, "length": 30,
                 "datum": "outer", "position": "outside", "side": "up" {own} }}"#
        )
    }

    #[test]
    fn features_bend_with_their_own_radius_and_k_factor() {
        let own = r#", "bend_radius": 4, "k_factor": 0.45"#;
        let tray: Vec<String> = (0..4).map(|e| flange(e, own)).collect();
        // Expected: volume, faces, flat extent across the bend.
        for (name, spec, volume, faces, flat) in [
            (
                "flange, K 0.3: a shorter flat",
                part(&format!(
                    r#""flanges": [{}]"#,
                    flange(2, r#", "k_factor": 0.3"#)
                )),
                17742.4778,
                14,
                88.0841,
            ),
            (
                "flange, R5",
                part(&format!(
                    r#""flanges": [{}]"#,
                    flange(2, r#", "bend_radius": 5"#)
                )),
                17484.9556,
                14,
                f64::NAN,
            ),
            (
                "R5 next to the part's R2: different settings",
                part(&format!(
                    r#""flanges": [{}, {}]"#,
                    flange(1, ""),
                    flange(2, r#", "bend_radius": 5"#)
                )),
                20773.6193,
                23,
                f64::NAN,
            ),
            (
                "tray, R4 K 0.45",
                part(&format!(r#""flanges": [{}]"#, tray.join(","))),
                28908.8354,
                46,
                f64::NAN,
            ),
            (
                "contour, R4 K 0.45",
                part(
                    r#""contours": [{ "id": "c", "on": { "base": 2 }, "profile": [[0,2],[0,32],[-15,32]],
                                      "face": "bottom", "bend_radius": 4, "k_factor": 0.45 }]"#,
                ),
                19741.5927,
                22,
                98.3938,
            ),
            (
                "hem, K 0.3",
                part(
                    r#""hems": [{ "id": "h", "on": { "base": 2 },
                                  "kind": { "single": { "gap": 1, "length": 8 } }, "side": "up",
                                  "k_factor": 0.3 }]"#,
                ),
                13542.4778,
                14,
                66.4558,
            ),
            (
                "fold, R6 K 0.33",
                part(
                    r#""folds": [{ "id": "f", "line": [[100,40],[0,40]], "angle": 90,
                                   "location": "centerline", "side": "up",
                                   "bend_radius": 6, "k_factor": 0.33 }]"#,
                ),
                12106.8142,
                14,
                60.0,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
            if flat.is_finite() {
                let ys: Vec<f64> = crate::flat::flat_pattern(&spec)
                    .unwrap()
                    .outline
                    .iter()
                    .map(|p| p[1])
                    .collect();
                let extent = ys.iter().cloned().fold(f64::MIN, f64::max)
                    - ys.iter().cloned().fold(f64::MAX, f64::min);
                assert!((extent - flat).abs() < 1e-3, "{name}: flat {extent}");
            }
        }
    }
}

#[cfg(test)]
mod remnant_tests {
    use super::*;

    fn part(start: f64, end: f64) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [{{ "id": "f", "on": {{ "base": 2 }},
                                "width": {{ "offset_start": {start}, "offset_end": {end} }},
                                "angle": 90, "length": 30, "datum": "outer",
                                "position": "outside", "side": "up" }}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn a_remnant_no_longer_than_the_minimum_runs_the_notch_to_the_corner() {
        // Expected: volume, faces. Relief 2 wide, minimum remnant 4.
        for (name, spec, volume, faces) in [
            ("remnant 1", part(3.0, 30.0), 15817.4601, 20),
            ("remnants 3 and 1", part(5.0, 3.0), 17235.0796, 18),
            ("remnant exactly 4", part(6.0, 30.0), 15627.1858, 20),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod bend_table_tests {
    use super::*;

    /// A DIN 6935 bend table for 2 mm sheet.
    const TABLE: &str = include_str!("testdata/din-bend-table-2mm.json");

    fn part(features: &str, radius: f64) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": {radius}, "k_factor": 0.4,
                 "bend_table": {TABLE},
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }}, {features} }}"#
        ))
        .unwrap()
    }

    fn flange(angle: f64) -> String {
        format!(
            r#""flanges": [{{ "id": "f", "on": {{ "base": 2 }}, "angle": {angle}, "length": 30,
                             "datum": "outer", "position": "outside", "side": "up" }}]"#
        )
    }

    #[test]
    fn bend_tables_give_their_flat_lengths() {
        let hem = r#""hems": [{ "id": "h", "on": { "base": 2 },
                                 "kind": { "single": { "gap": 1, "length": 8 } }, "side": "up" }]"#;
        // Expected: the flat blank across the bend; K 0.4 gives 88.3982 at 90.
        for (name, spec, flat) in [
            (
                "90: legs to the virtual sharp",
                part(&flange(90.0), 2.0),
                88.1626,
            ),
            (
                "135: legs to the outer tangent",
                part(&flange(135.0), 2.0),
                81.7586,
            ),
            ("72: between two rows", part(&flange(72.0), 2.0), 88.9667),
            ("hem: opening 0", part(hem, 2.0), 65.6671),
        ] {
            let ys: Vec<f64> = crate::flat::flat_pattern(&spec)
                .unwrap()
                .outline
                .iter()
                .map(|p| p[1])
                .collect();
            let extent = ys.iter().cloned().fold(f64::MIN, f64::max)
                - ys.iter().cloned().fold(f64::MAX, f64::min);
            assert!((extent - flat).abs() < 1e-3, "{name}: {extent}");
            exact_mass(&spec).unwrap();
        }
        // Off the table's radii: refused.
        assert!(matches!(
            exact_mass(&part(&flange(90.0), 8.0)),
            Err(SheetError::Invalid(_))
        ));
    }
}

#[cfg(test)]
mod notch_hole_tests {
    use super::*;

    #[test]
    fn holes_across_the_outline_are_round_notches() {
        let part = |holes: &str, flanges: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]], "holes": {holes} }},
                     "flanges": [{flanges}] }}"#
            ))
            .unwrap()
        };
        let flange = r#"{ "id": "f", "on": { "base": 2 }, "angle": 90, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }"#;
        // Expected: volume, faces. Half a hole on an edge, a quarter at a
        // corner; on a flanged plate the base's edge, the flange's side and tip.
        for (name, spec, volume, faces) in [
            (
                "plate",
                part(
                    r#"[{ "center": [50, 0], "diameter": 10 }, { "center": [100, 60], "diameter": 8 }]"#,
                    "",
                ),
                11896.3274,
                9,
            ),
            (
                "flanged",
                part(
                    r#"[{ "center": [0, 30], "diameter": 10 }, { "center": [100, 75], "diameter": 8 },
                        { "center": [50, 88.3982297150257], "diameter": 8 }]"#,
                    flange,
                ),
                17563.407,
                20,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
            assert!(crate::flat::flat_pattern(&spec).unwrap().holes.is_empty());
        }
    }
}

#[cfg(test)]
mod mixed_relief_tests {
    use super::*;

    #[test]
    fn two_features_corner_takes_no_corner_relief() {
        let part = |corner: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "relief": {{ "width": 2, "depth": 1, "minimum_remnant": 4 {corner} }},
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "flanges": [
                       {{ "id": "f1", "on": {{ "base": 1 }}, "angle": 90, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }},
                       {{ "id": "f2", "on": {{ "base": 2 }}, "angle": 90, "length": 20,
                          "datum": "outer", "position": "outside", "side": "up" }}] }}"#
            ))
            .unwrap()
        };
        // Expected: the same part whatever the relief says.
        for corner in [
            "",
            r#", "corner": { "shape": "round" }"#,
            r#", "corner": { "shape": "linear_weld" }"#,
        ] {
            let m = exact_mass(&part(corner)).unwrap();
            assert!(
                (m.volume - 19094.2654).abs() < 1e-3,
                "{corner}: {}",
                m.volume
            );
        }
    }
}

#[cfg(test)]
mod relieved_lip_tests {
    use super::*;

    fn part(relief: &str, lips: &str, hems: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, {relief}
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "flanges": [
                   {{ "id": "f1", "on": {{ "base": 1 }}, "angle": 90, "length": 30,
                      "datum": "outer", "position": "outside", "side": "up" }} {lips}],
                 "hems": [{hems}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn lips_at_relieved_corners_and_double_hems_next_to_a_flange() {
        let f2_lip = r#", { "id": "f2", "on": { "base": 2 }, "angle": 90, "length": 30,
                            "datum": "outer", "position": "outside", "side": "up" },
                          { "id": "lip", "on": { "tip": "f2" }, "angle": 90, "length": 15,
                            "datum": "outer", "position": "outside", "side": "up" }"#;
        let relief = |shape: &str| {
            format!(
                r#""relief": {{ "width": 2, "depth": 1, "minimum_remnant": 4,
                               "corner": {{ "shape": "{shape}", "size": 10 }} }},"#
            )
        };
        let double = r#"{ "id": "h", "on": { "base": 2 },
                          "kind": { "double": { "gap": 1, "length": 10 } }, "side": "up" }"#;
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            (
                "lip, welded corner",
                part(&relief("linear_weld"), f2_lip, ""),
                23720.4014,
                30,
            ),
            (
                "lip, round relief",
                part(&relief("round"), f2_lip, ""),
                23640.2826,
                33,
            ),
            (
                "double hem next to a flange",
                part("", "", double),
                19327.1812,
                39,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod mixed_lip_tests {
    use super::*;

    #[test]
    fn lips_at_corners_between_separate_flanges() {
        let part = |lip_on: &str, length: f64| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "flanges": [
                       {{ "id": "f1", "on": {{ "base": 1 }}, "angle": 90, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }},
                       {{ "id": "f2", "on": {{ "base": 2 }}, "angle": 90, "length": 20,
                          "datum": "outer", "position": "outside", "side": "up" }},
                       {{ "id": "lip", "on": {{ "tip": "{lip_on}" }}, "angle": 90,
                          "length": {length}, "datum": "outer", "position": "outside",
                          "side": "up" }}] }}"#
            ))
            .unwrap()
        };
        // Expected: volume, faces.
        for (name, spec, volume, faces) in [
            ("on the first-built", part("f1", 12.0), 20379.752, 31),
            ("on the second-built", part("f2", 10.0), 20801.8936, 31),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod contour_end_tests {
    use super::*;

    #[test]
    fn features_sit_on_a_contours_free_end() {
        let part = |profile: &str, features: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                     "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                     "contours": [{{ "id": "c", "on": {{ "base": 2 }}, "profile": {profile},
                                     "face": "bottom" }}], {features} }}"#
            ))
            .unwrap()
        };
        // Expected: volume, faces. {"tip": "c"} is the profile's end.
        for (name, spec, volume, faces) in [
            (
                "a hem on an L",
                part(
                    "[[0,2],[0,32],[-15,32]]",
                    r#""hems": [{ "id": "h", "on": { "tip": "c" },
                                  "kind": { "single": { "gap": 1, "length": 8 } }, "side": "down" }]"#,
                ),
                21627.4334,
                30,
            ),
            (
                "a lip on a Z",
                part(
                    "[[0,2],[0,22],[20,22],[20,42]]",
                    r#""flanges": [{ "id": "lip", "on": { "tip": "c" }, "angle": 90, "length": 10,
                                     "datum": "outer", "position": "outside", "side": "up" }]"#,
                ),
                24969.9112,
                38,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            assert_eq!(fold(&spec).unwrap().topology.faces.len(), faces, "{name}");
        }
    }
}

#[cfg(test)]
mod jog_tests {
    use super::*;

    fn part(jog: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0,0],[100,0],[100,60],[0,60]] }},
                 "jogs": [{jog}] }}"#
        ))
        .unwrap()
    }

    #[test]
    fn jogs_step_by_their_offset() {
        // Expected (two folds): volume, faces, top of the part.
        for (name, spec, volume, faces, top) in [
            (
                "up 15",
                part(r#"{ "id": "j", "line": [[100,25],[0,25]], "offset": 15, "side": "up" }"#),
                12125.6637,
                22,
                17.0,
            ),
            (
                "down 12 at 60 degrees",
                part(
                    r#"{ "id": "j", "line": [[100,25],[0,25]], "offset": 12, "angle": 60,
                         "side": "down" }"#,
                ),
                12083.7758,
                22,
                2.0,
            ),
        ] {
            let m = exact_mass(&spec).unwrap();
            assert!((m.volume - volume).abs() < 1e-3, "{name}: {}", m.volume);
            let brep = fold(&spec).unwrap();
            assert_eq!(brep.topology.faces.len(), faces, "{name}");
            let zmax = brep
                .topology
                .vertices
                .iter()
                .map(|v| v.position[2])
                .fold(f64::MIN, f64::max);
            assert!((zmax - top).abs() < 1e-9, "{name}: top {zmax}");
        }
        // Too small for its two bends.
        let low = part(r#"{ "id": "j", "line": [[100,25],[0,25]], "offset": 3, "side": "up" }"#);
        assert!(matches!(exact_mass(&low), Err(SheetError::Invalid(_))));
    }
}

/// Frame conversions the kernel keeps to itself.
trait FrameExt {
    fn to_local(&self, p: [f64; 3]) -> [f64; 3];
    fn to_world(&self, v: [f64; 3]) -> [f64; 3];
}

impl FrameExt for Frame3 {
    fn to_local(&self, p: [f64; 3]) -> [f64; 3] {
        let d = sub(p, self.origin);
        [dot(d, self.x), dot(d, self.y), dot(d, self.z)]
    }

    fn to_world(&self, v: [f64; 3]) -> [f64; 3] {
        add(
            add(scale(self.x, v[0]), scale(self.y, v[1])),
            scale(self.z, v[2]),
        )
    }
}

#[cfg(test)]
fn dist2(a: P2, b: P2) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2)).sqrt()
}
