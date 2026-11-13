//! The flat pattern, built directly in flat space (flat-first): the base stays
//! where the spec puts it, and each flange adds a bend strip and a flange strip
//! beyond its base edge. Exact by construction; nothing is unfolded.
//!
//! Each base edge `a -> b` is laid out in its own coordinates: `t` along the
//! edge from `a`, `s` along the outward normal from the edge.

use crate::cut::Clip;
use crate::spec::{
    CornerCut, CornerShape, EdgeRef, Frame2, Hole, PartSpec, Punch, ReliefShape, Seam, Side,
    Transition,
};
use crate::SheetError;

pub type P2 = [f64; 2];

#[derive(Clone, Debug, PartialEq)]
pub struct FlatPattern {
    /// Counter-clockwise outer boundary, base frame (z = 0 face seen from +z).
    pub outline: Vec<P2>,
    /// Per outline segment (`outline[k]` to `outline[k + 1]`): an arc, or a line.
    pub outline_arcs: Vec<Option<Curve>>,
    pub holes: Vec<Hole>,
    /// Cut-outs, counter-clockwise.
    pub cutouts: Vec<Ring>,
    pub bends: Vec<FlatBend>,
    /// Every punch stamp, in spec order.
    pub punches: Vec<FlatPunch>,
    /// Every emboss's footprint (where the sheet leaves the plate). Not cut:
    /// the plate stays whole in the flat.
    pub forms: Vec<FlatForm>,
}

/// Where a formed feature leaves the flat plate.
#[derive(Clone, Debug, PartialEq)]
pub struct FlatForm {
    pub id: String,
    pub center: P2,
    pub outline: Ring,
}

/// Where a punch stamps the flat pattern.
#[derive(Clone, Debug, PartialEq)]
pub struct FlatPunch {
    pub id: String,
    pub center: P2,
    /// Degrees, counter-clockwise.
    pub angle: f64,
}

/// A closed boundary of lines and arcs: `arcs[k]` is the segment starting at
/// `points[k]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Ring {
    pub points: Vec<P2>,
    pub arcs: Vec<Option<Curve>>,
}

impl Ring {
    /// From a spec cut-out: bulges to arc centres.
    pub fn from_cutout(c: &crate::spec::Cutout) -> Ring {
        let n = c.outline.len();
        let arcs = (0..n)
            .map(|k| {
                let b = c.bulges.get(k).copied().unwrap_or(0.0);
                if b.abs() < 1e-12 {
                    return None;
                }
                let (p, q) = (c.outline[k], c.outline[(k + 1) % n]);
                let sweep = 4.0 * b.atan();
                let d = sub(q, p);
                let chord = dot(d, d).sqrt();
                let left = [-d[1] / chord, d[0] / chord];
                let mid = [(p[0] + q[0]) / 2.0, (p[1] + q[1]) / 2.0];
                let h = chord / 2.0 / (sweep / 2.0).tan();
                Some(Curve::Arc(Arc {
                    center: add(mid, scale(left, h)),
                    ccw: b > 0.0,
                }))
            })
            .collect();
        Ring {
            points: c.outline.clone(),
            arcs,
        }
    }

    /// A hole's circle, as two half circles.
    pub fn from_hole(h: &Hole) -> Ring {
        let (c, r) = (h.center, h.diameter / 2.0);
        let arc = Some(Curve::Arc(Arc {
            center: c,
            ccw: true,
        }));
        Ring {
            points: vec![[c[0] + r, c[1]], [c[0] - r, c[1]]],
            arcs: vec![arc, arc],
        }
    }

    /// The same boundary run the other way.
    pub fn reversed(&self) -> Ring {
        let n = self.points.len();
        let points: Vec<P2> = (0..n).map(|k| self.points[(n - k) % n]).collect();
        // Segment k of the reversed ring is segment (n - 1 - k) of this one.
        let arcs = (0..n)
            .map(|k| {
                self.arcs[n - 1 - k].map(|c| match c {
                    Curve::Arc(a) => Curve::Arc(Arc {
                        center: a.center,
                        ccw: !a.ccw,
                    }),
                    // The same curve; its ends swap with the ring's points.
                    Curve::Cut(_) => c,
                })
            })
            .collect();
        Ring { points, arcs }
    }

    pub fn area(&self) -> f64 {
        path_area(&self.points, &self.arcs)
    }
}

/// Rounds the corner at `v`, coming in along unit `d_in` and leaving along
/// unit `d_out`, with radius `r`: a line to where the arc starts, then the
/// arc.
fn fillet(path: &mut Path, v: P2, d_in: P2, d_out: P2, r: f64) {
    let (p, q, center) = fillet_points(v, d_in, d_out, r);
    let mid = sub(scale(add(p, q), 0.5), center);
    let via = add(center, scale(mid, r / dot(mid, mid).sqrt()));
    path.line(p);
    path.arc(q, center, via);
}

/// A fillet's ends (where it leaves the incoming line, where it meets the
/// outgoing one) and its centre.
fn fillet_points(v: P2, d_in: P2, d_out: P2, r: f64) -> (P2, P2, P2) {
    let turn = cross2(d_in, d_out);
    let along = r * (turn.abs().atan2(dot(d_in, d_out)) / 2.0).tan();
    let (p, q) = (sub(v, scale(d_in, along)), add(v, scale(d_out, along)));
    let left = [-d_in[1], d_in[0]];
    (p, q, add(p, scale(left, r * turn.signum())))
}

/// Lowest and highest distance of a ring from `origin` along unit `d`,
/// arcs exact.
fn extent_along(ring: &Ring, d: P2, origin: P2) -> (f64, f64) {
    let n = ring.points.len();
    let at = |p: P2| dot(sub(p, origin), d);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for k in 0..n {
        let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
        let mut seen = vec![at(p)];
        if let Some(Curve::Arc(a)) = ring.arcs[k] {
            // The arc's points furthest along +d and -d, where it passes them.
            let r = dot(sub(p, a.center), sub(p, a.center)).sqrt();
            let sweep = a.sweep(p, q);
            for sign in [1.0, -1.0] {
                let x = add(a.center, scale(d, sign * r));
                let to = a.sweep(p, x);
                if to.signum() == sweep.signum() && to.abs() <= sweep.abs() {
                    seen.push(at(x));
                }
            }
        }
        for v in seen {
            (lo, hi) = (lo.min(v), hi.max(v));
        }
    }
    (lo, hi)
}

/// A boundary segment that runs around `center` instead of straight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arc {
    pub center: P2,
    /// Counter-clockwise from the segment's start to its end (seen from +z).
    pub ccw: bool,
}

impl Arc {
    /// Signed sweep from `from` to `to` (positive counter-clockwise).
    pub fn sweep(&self, from: P2, to: P2) -> f64 {
        let angle = |p: P2| (p[1] - self.center[1]).atan2(p[0] - self.center[0]);
        let d = (angle(to) - angle(from)).rem_euclid(std::f64::consts::TAU);
        let d = if d < EPS { std::f64::consts::TAU } else { d };
        if self.ccw {
            d
        } else {
            d - std::f64::consts::TAU
        }
    }
}

/// A boundary segment that isn't straight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Curve {
    Arc(Arc),
    /// The end of a bend mitered at a corner.
    Cut(Cut),
}

impl Curve {
    /// The arc, if it is one.
    pub fn arc(self) -> Option<Arc> {
        match self {
            Curve::Arc(a) => Some(a),
            Curve::Cut(_) => None,
        }
    }
}

/// The end of a bend mitered at a corner, in the flat: the point at flat
/// distance s along `n` (from the bend's edge through `origin`) lies at
/// t = c0 + c1 cos(phi) + c2 sin(phi) along `u`, phi = (s - s0) / radius
/// (the bend angle; `radius` its neutral radius).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cut {
    pub origin: P2,
    pub u: P2,
    pub n: P2,
    pub s0: f64,
    pub radius: f64,
    pub c: [f64; 3],
}

impl Cut {
    /// Bend angle at flat s.
    pub fn phi(&self, s: f64) -> f64 {
        (s - self.s0) / self.radius
    }
    pub fn t(&self, s: f64) -> f64 {
        let (sin, cos) = self.phi(s).sin_cos();
        self.c[0] + self.c[1] * cos + self.c[2] * sin
    }
    pub fn point(&self, s: f64) -> P2 {
        add(self.origin, add(scale(self.u, self.t(s)), scale(self.n, s)))
    }
    /// Flat s of a point on it.
    pub fn s_of(&self, p: P2) -> f64 {
        dot(sub(p, self.origin), self.n)
    }
    /// Points from `from` to `to` (both on it), `pieces` + 1 of them.
    pub fn points(&self, from: P2, to: P2, pieces: usize) -> Vec<P2> {
        let (a, b) = (self.s_of(from), self.s_of(to));
        (0..=pieces)
            .map(|j| self.point(a + (b - a) * j as f64 / pieces as f64))
            .collect()
    }
}

/// Boundary points with an optional arc per segment; `arcs[k]` is the segment
/// starting at `pts[k]`.
#[derive(Clone, Debug, Default)]
struct Path {
    pts: Vec<P2>,
    arcs: Vec<Option<Curve>>,
}

impl Path {
    /// A straight segment from the last point to `p`.
    fn line(&mut self, p: P2) {
        if self
            .pts
            .last()
            .is_some_and(|q| dot(sub(*q, p), sub(*q, p)) < EPS * EPS)
        {
            return;
        }
        self.pts.push(p);
        self.arcs.push(None);
    }
    /// An arc around `center` from the last point to `p`, passing `via`.
    fn arc(&mut self, p: P2, center: P2, via: P2) {
        let from = *self.pts.last().unwrap();
        let (a, b) = (sub(from, center), sub(via, center));
        *self.arcs.last_mut().unwrap() = Some(Curve::Arc(Arc {
            center,
            ccw: a[0] * b[1] - a[1] * b[0] > 0.0,
        }));
        self.line(p);
    }
    /// A curve (a mitered bend's end) from the last point to `p`.
    fn curve(&mut self, p: P2, c: Curve) {
        *self.arcs.last_mut().unwrap() = Some(c);
        self.line(p);
    }
    fn close(&mut self) {
        while self.pts.len() > 1
            && dot(
                sub(self.pts[0], self.pts[self.pts.len() - 1]),
                sub(self.pts[0], self.pts[self.pts.len() - 1]),
            ) < EPS * EPS
        {
            self.pts.pop();
            self.arcs.pop();
        }
    }
}

/// Exact signed area of a closed boundary with arc segments (Green's theorem).
pub fn path_area(pts: &[P2], arcs: &[Option<Curve>]) -> f64 {
    (0..pts.len())
        .map(|k| {
            let (p, q) = (pts[k], pts[(k + 1) % pts.len()]);
            match arcs.get(k).copied().flatten() {
                None => (p[0] * q[1] - q[0] * p[1]) / 2.0,
                Some(Curve::Cut(c)) => {
                    // (x y' - y x') / 2 over s: smooth, so Gauss-Legendre in
                    // a few pieces is exact to rounding.
                    let (a, b) = (c.s_of(p), c.s_of(q));
                    let pieces = 8;
                    let mut sum = 0.0;
                    for j in 0..pieces {
                        let (s0, s1) = (
                            a + (b - a) * j as f64 / pieces as f64,
                            a + (b - a) * (j + 1) as f64 / pieces as f64,
                        );
                        for (x, w) in crate::fold::GAUSS {
                            let s = s0 + (s1 - s0) * x;
                            let (pt, h) = (c.point(s), 1e-6);
                            let d = scale(sub(c.point(s + h), c.point(s - h)), 0.5 / h);
                            sum += w * (s1 - s0) * (pt[0] * d[1] - pt[1] * d[0]) / 2.0;
                        }
                    }
                    sum
                }
                Some(Curve::Arc(arc)) => {
                    let c = arc.center;
                    let r = dot(sub(p, c), sub(p, c)).sqrt();
                    let a0 = (p[1] - c[1]).atan2(p[0] - c[0]);
                    let a1 = a0 + arc.sweep(p, q);
                    (r * r * (a1 - a0) + r * c[0] * (a1.sin() - a0.sin())
                        - r * c[1] * (a1.cos() - a0.cos()))
                        / 2.0
                }
            }
        })
        .sum()
}

#[derive(Clone, Debug, PartialEq)]
pub struct FlatBend {
    pub flange: String,
    /// When this bend is made, from 1. The pieces of one fold line share it.
    pub order: usize,
    /// Bend centre line across the flange width.
    pub line: [P2; 2],
    pub angle: f64,
    pub inner_radius: f64,
    pub k_factor: f64,
    /// Bends toward +z (the base's top face).
    pub up: bool,
}

impl FlatBend {
    /// The bend's note: order, direction, angle and inner radius, as
    /// "2 UP 90° R2".
    pub fn note(&self) -> String {
        let short = |v: f64| {
            let s = format!("{v:.2}");
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        };
        format!(
            "{} {} {}° R{}",
            self.order,
            if self.up { "UP" } else { "DOWN" },
            short(self.angle),
            short(self.inner_radius)
        )
    }

    /// Where a note reads along the bend line: its midpoint, a distance d to
    /// the left of the line as it reads, and the reading direction in degrees
    /// (-90, 90].
    pub fn note_place(&self, d: f64) -> (P2, f64) {
        let [a, b] = self.line;
        let mut angle = (b[1] - a[1]).atan2(b[0] - a[0]).to_degrees();
        if angle > 90.0 + 1e-9 {
            angle -= 180.0;
        } else if angle <= -90.0 + 1e-9 {
            angle += 180.0;
        }
        let (s, c) = angle.to_radians().sin_cos();
        let mid = [(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0];
        ([mid[0] - s * d, mid[1] + c * d], angle)
    }
}

/// What a base edge carries: one flange strip, possibly with lips on it.
#[derive(Clone, Copy)]
struct Strip {
    /// Flange span along the edge.
    t0: f64,
    t1: f64,
    /// Where the outermost flange of the chain ends, along the normal.
    end: f64,
    /// Where the first bend starts (negative: inside the base outline).
    setback: f64,
}

/// How one end of a strip meets its surroundings.
#[derive(Clone, Copy, Debug)]
enum End {
    /// Flush with a square base edge that carries no flange there.
    Free,
    /// A partial flange's end: a relief notch.
    Relief,
    /// A strip's end next to a base edge that isn't square to it and carries
    /// no flange there, or an arc at an angle to it: the flange ends square
    /// at the vertex; the bend ends in a slanted cut from `bend` (where its
    /// start line meets the neighbouring edge) to the vertex at its far edge.
    /// `flat`: how the flat pattern shows it (the bend transition).
    Open { bend: f64, flat: Flat },
    /// Such an end trimmed to the bend: the bend ends square at `bend` (where
    /// it first meets the neighbouring edge or the vertex), the flange at the
    /// vertex; the base runs along the bend start line to `q` on the
    /// neighbouring edge.
    Square { bend: f64, q: P2 },
    /// A relief notch that runs on to the neighbouring edge: next to a
    /// tangent arc (a rounded corner), a notch to `corner` (its far bottom
    /// corner) whose far side runs on to the arc at `p`; without `corner` (the
    /// remnant would be no longer than the minimum, also next to a straight
    /// edge) the notch's bottom runs on to the edge at `p`. `p` trims an arc.
    /// `round`: a round relief's run-on, a quarter circle from the strip's
    /// side down to the notch bottom.
    RunOn {
        corner: Option<P2>,
        p: P2,
        round: bool,
    },
    /// Two flanges meet: the bend runs to `bend`, the flange to `flange` at
    /// its start (positions along this edge), and on at `slope` per unit
    /// along the edge normal. When the flange starts short of where the bend
    /// start lines cross (`short`, that point), the bend stops there too and
    /// the base keeps its corner at that point.
    Corner {
        bend: f64,
        flange: f64,
        slope: f64,
        short: Option<P2>,
        /// The corner relief: trim to bend between separate features.
        shape: CornerShape,
        /// The bend runs on past q to its flange's end (a small round or
        /// square relief); the base still ends at q.
        run_on: RunOn,
        /// A run-on bend's end turns onto the relief: its centre, where the
        /// end meets it, and the end's own corner before that (the meeting
        /// point when it meets the end line itself). The bends split along q
        /// to the centre; the relief cuts the rest.
        turn: (P2, P2, P2),
        /// One contour along both edges: every strip past the first bend ends
        /// on the miter too.
        miter: Option<Miter>,
    },
}

/// How the flat pattern shows an open end (the folded part is the same).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Flat {
    /// As folded: the bend's slanted cut.
    Slanted,
    /// The strip ending square through the vertex.
    Square,
    /// The same, its corner at the vertex rounded with the bend radius.
    Fillet,
}

/// A mitered corner of a contour along several edges: each strip ends where
/// the first of its faces lies `half_gap` from the corner's symmetry plane
/// (through `vertex`, unit normal `plane` toward this edge).
#[derive(Clone, Copy, Debug)]
struct Miter {
    vertex: P2,
    plane: P2,
    half_gap: f64,
}

/// Where a bend at a corner with a small round or square relief runs past q
///.
#[derive(Clone, Copy, PartialEq, Debug)]
enum RunOn {
    No,
    /// In the flat pattern only: a relief that doesn't reach
    /// the circle through q and both flanges' first corners folds as trim to bend.
    Flat,
    /// In the flat pattern and folded, with the relief cut from both.
    Both,
}

const EPS: f64 = 1e-9;
/// Lines along a mitered bend's curved end in the cells, until the 3D build
/// has a surface for it.
const MITER_SAMPLES: usize = 96;

/// A boundary for formats without the mitered bend end's curve: those ends
/// as fine lines, arcs kept.
pub fn lines_and_arcs(pts: &[P2], curves: &[Option<Curve>]) -> (Vec<P2>, Vec<Option<Arc>>) {
    let (mut p, mut c) = (pts.to_vec(), curves.to_vec());
    lines_for_cuts(&mut p, &mut c);
    (p, c.into_iter().map(|c| c.and_then(Curve::arc)).collect())
}

/// A boundary's mitered bend ends as fine lines.
fn lines_for_cuts(pts: &mut Vec<P2>, arcs: &mut Vec<Option<Curve>>) {
    if !arcs.iter().any(|c| matches!(c, Some(Curve::Cut(_)))) {
        return;
    }
    let n = pts.len();
    let (mut p2, mut a2) = (Vec::new(), Vec::new());
    for k in 0..n {
        let (p, q) = (pts[k], pts[(k + 1) % n]);
        match arcs[k] {
            Some(Curve::Cut(c)) => {
                let line = c.points(p, q, MITER_SAMPLES);
                p2.extend(&line[..line.len() - 1]);
                a2.extend(std::iter::repeat_n(None, line.len() - 1));
            }
            other => {
                p2.push(p);
                a2.push(other);
            }
        }
    }
    *pts = p2;
    *arcs = a2;
}

/// x in rows of [a, b, c, d] with a x0 + b x1 + c x2 = d (Cramer).
fn solve3(r: [[f64; 4]; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let a = r.map(|row| [row[0], row[1], row[2]]);
    let d = det(a);
    std::array::from_fn(|j| {
        let mut m = a;
        for (i, row) in m.iter_mut().enumerate() {
            row[j] = r[i][3];
        }
        det(m) / d
    })
}
/// Width and depth of a tear relief's notch in 3D (fixed, whatever the
/// thickness or the relief settings).
const TEAR: f64 = 0.01;

/// The flat pattern split into the cells that fold differently, plus how each
/// flange folds. Input for the 3D model.
#[derive(Clone, Debug)]
pub struct Layout {
    pub flat: FlatPattern,
    /// The base (what stays in place), the bend strips and the flange strips,
    /// with the holes cut out. Together they tile the flat pattern.
    pub cells: Vec<Cell>,
    pub folds: Vec<Fold>,
    /// Every hole in the flat, in spec order: the base's, then those on
    /// flange faces.
    pub holes: Vec<Hole>,
    /// The embosses, each held whole by a base cell as a hole the size of
    /// its footprint.
    pub(crate) embosses: Vec<crate::emboss::Formed>,
}

#[derive(Clone, Debug)]
pub struct Cell {
    pub name: String,
    /// Counter-clockwise, with every junction with a neighbouring cell as a
    /// vertex.
    pub polygon: Vec<P2>,
    /// Per segment: an arc (relief ends, holes crossing the cell), or a line.
    pub arcs: Vec<Option<Curve>>,
    /// Holes that lie wholly inside the cell.
    pub holes: Vec<Hole>,
    /// Cut-outs that lie wholly inside the cell, clockwise.
    pub inner: Vec<Ring>,
    pub zone: Zone,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Zone {
    /// The part of the base that stays in place.
    Base,
    /// Bend strip of `folds[i]`.
    Bend(usize),
    /// Flat flange strip after `folds[i]`.
    Flange(usize),
}

/// One bend, in flat coordinates: the strip `s0 .. s0 + bend_allowance` along
/// `n` from the base edge through `edge_point` (direction `u`).
#[derive(Clone, Debug)]
pub struct Fold {
    pub flange: String,
    /// The fold this one sits on (its tip), `None` on the base.
    pub parent: Option<usize>,
    pub edge_point: P2,
    pub u: P2,
    pub n: P2,
    pub s0: f64,
    pub bend_allowance: f64,
    /// Bend angle in radians.
    pub angle: f64,
    /// Inner bend radius.
    pub radius: f64,
    pub k_factor: f64,
    /// Bends toward +z of the flat frame.
    pub up: bool,
}

pub fn flat_pattern(spec: &PartSpec) -> Result<FlatPattern, SheetError> {
    layout(spec).map(|l| l.flat)
}

pub fn layout(spec: &PartSpec) -> Result<Layout, SheetError> {
    if let Some(l) = crate::cone::layout(spec) {
        return l;
    }
    let spec = &*spec.resolved()?;
    let base = &spec.base.outline;
    let n = base.len();
    if n < 3 || signed_area(base) <= 0.0 {
        return Err(SheetError::Invalid(
            "base outline must be a counter-clockwise polygon".into(),
        ));
    }
    let edge = |i: usize| {
        let (a, b) = (base[i], base[(i + 1) % n]);
        let length = dot(sub(b, a), sub(b, a)).sqrt();
        let u = scale(sub(b, a), 1.0 / length);
        (a, u, [u[1], -u[0]], length)
    };
    let at = |i: usize, t: f64, s: f64| {
        let (a, u, normal, _) = edge(i);
        add(a, add(scale(u, t), scale(normal, s)))
    };
    // Per outline edge, its arc (rounded corners and the like).
    let base_arcs = Ring::from_cutout(&crate::spec::Cutout {
        outline: base.clone(),
        bulges: spec.base.bulges.clone(),
    })
    .arcs
    .into_iter()
    .map(|c| c.and_then(Curve::arc))
    .collect::<Vec<Option<Arc>>>();

    let mut strips: Vec<Option<Strip>> = vec![None; n];
    // Per placed flange: id, base edge, its flat end, bends toward +z.
    let mut placed: Vec<(&str, usize, f64, bool)> = Vec::new();
    let mut occupied = std::collections::HashSet::new();
    let mut bends = Vec::new();
    let mut folds: Vec<Fold> = Vec::new();
    // Per fold: base edge, its flat end.
    let mut fold_ends: Vec<(usize, f64)> = Vec::new();
    let steps = crate::hem::steps(spec)?;
    for step in &steps {
        if !occupied.insert(step.on.clone()) {
            return Err(SheetError::Invalid(format!(
                "{}: its edge already has a flange",
                step.id
            )));
        }
        let (i, start, up) = match &step.on {
            EdgeRef::Base(i) => {
                let i = *i;
                if i >= n {
                    return Err(SheetError::Invalid(format!(
                        "{}: no base edge {i}",
                        step.id
                    )));
                }
                if base_arcs[i].is_some() {
                    return Err(SheetError::Unsupported(format!(
                        "{}: on an arc edge",
                        step.id
                    )));
                }
                let length = edge(i).3;
                let (t0, t1) = step
                    .width
                    .map_or((0.0, length), |w| (w.offset_start, length - w.offset_end));
                if !(t0 >= 0.0 && t1 <= length && t1 - t0 > EPS) {
                    return Err(SheetError::Invalid(format!(
                        "{}: width doesn't fit on edge {i}",
                        step.id
                    )));
                }
                strips[i] = Some(Strip {
                    t0,
                    t1,
                    end: 0.0,
                    setback: step.setback,
                });
                (i, 0.0, step.side == Side::Up)
            }
            EdgeRef::Tip(parent) => {
                if step.width.is_some() {
                    return Err(SheetError::Invalid(format!(
                        "{}: a flange on a tip spans its parent",
                        step.id
                    )));
                }
                let &(_, i, end, parent_up) = placed
                    .iter()
                    .find(|p| p.0 == parent.as_str())
                    .ok_or_else(|| {
                        SheetError::Invalid(format!("{}: no earlier flange {parent}", step.id))
                    })?;
                // Same turning sense as the parent means the same flat bend direction.
                (i, end, parent_up == (step.side == Side::Up))
            }
            // Resolved to base edges before.
            EdgeRef::Face(e) => {
                return Err(SheetError::Invalid(format!(
                    "{}: face {} has no bend",
                    step.id, e.face
                )))
            }
        };
        let strip = strips[i].as_mut().unwrap();
        strip.end = start + step.reach();
        let (a, u, normal, _) = edge(i);
        folds.push(Fold {
            flange: step.id.clone(),
            parent: match &step.on {
                EdgeRef::Base(_) | EdgeRef::Face(_) => None,
                EdgeRef::Tip(p) => folds.iter().position(|f| &f.flange == p),
            },
            edge_point: a,
            u,
            n: normal,
            s0: start + step.setback,
            bend_allowance: step.bend_allowance,
            angle: step.angle.to_radians(),
            radius: step.radius,
            k_factor: step.k_factor,
            up,
        });
        fold_ends.push((i, strip.end));
        let centre = start + step.setback + step.bend_allowance / 2.0;
        bends.push(FlatBend {
            flange: step.id.clone(),
            order: 0,
            line: [at(i, strip.t0, centre), at(i, strip.t1, centre)],
            angle: step.angle,
            inner_radius: step.radius,
            k_factor: step.k_factor,
            up,
        });
        placed.push((&step.id, i, strip.end, up));
    }

    let relief = spec.relief();
    // A bend relief's notch depth past the bend start (a tear's is fixed).
    let notch_depth = if relief.shape == ReliefShape::Tear {
        TEAR
    } else {
        relief.depth
    };
    // The fold that starts on each base edge, and whether a lip sits on it.
    let root: Vec<Option<usize>> = (0..n)
        .map(|i| (0..folds.len()).find(|&k| folds[k].parent.is_none() && fold_ends[k].0 == i))
        .collect();
    let has_lip = |k: usize| folds.iter().any(|f| f.parent == Some(k));

    // How each end of each strip meets its surroundings.
    let mut sides: Vec<[End; 2]> = vec![[End::Free, End::Free]; n];
    // Each with whether it cuts the cells too, or the flat pattern only.
    let mut corner_clips: Vec<(Clip, bool)> = Vec::new();
    // An open end by the bend transition: its bend's start line meets the
    // neighbouring edge at q (t = bend), the vertex is at t = vertex.
    let transition = match relief.transition.unwrap_or_default() {
        Transition::StraightLine => Transition::None,
        t => t,
    };
    let open_end = |bend: f64, vertex: f64, side: usize, q: P2| match transition {
        Transition::TrimToBend => End::Square {
            bend: if side == 0 {
                bend.max(vertex)
            } else {
                bend.min(vertex)
            },
            q,
        },
        Transition::Intersection => End::Open {
            bend,
            flat: Flat::Square,
        },
        Transition::Arc => End::Open {
            bend,
            flat: Flat::Fillet,
        },
        _ => End::Open {
            bend,
            flat: Flat::Slanted,
        },
    };
    // Strip ends cut back by an earlier flange's run-on notch: free.
    let mut cut = vec![[false; 2]; n];
    for i in 0..n {
        let Some(st) = strips[i] else { continue };
        let length = edge(i).3;
        let unsupported = |what: &str| {
            Err(SheetError::Unsupported(format!(
                "flange on edge {i}: {what}"
            )))
        };
        for (side, neighbour) in [(0, (i + n - 1) % n), (1, (i + 1) % n)] {
            if cut[i][side] {
                continue;
            }
            let at_corner = if side == 0 {
                st.t0 < EPS
            } else {
                st.t1 > length - EPS
            };
            let width = if relief.shape == ReliefShape::Tear {
                TEAR
            } else {
                relief.width
            };
            if !at_corner {
                let remnant = if side == 0 { st.t0 } else { length - st.t1 };
                // Next to an arc the remnant stays, however short.
                if base_arcs[neighbour].is_some() && remnant - width < EPS {
                    return unsupported("a relief notch reaching past its edge into an arc");
                }
                // An earlier flange on the neighbouring edge left this edge
                // shorter: to its bend start, or to its own run-on notch's
                // bottom.
                let shortened = match strips[neighbour] {
                    Some(o) if root[neighbour] < root[i] && base_arcs[neighbour].is_none() => {
                        let lj = edge(neighbour).3;
                        let (full, rest) = if side == 0 {
                            (o.t1 > lj - EPS, lj - o.t1)
                        } else {
                            (o.t0 < EPS, o.t0)
                        };
                        let s = if full {
                            Some(o.setback)
                        } else if rest - width <= relief.minimum_remnant + 1e-6 {
                            Some(o.setback - notch_depth)
                        } else {
                            None
                        };
                        let out = scale(edge(i).1, if side == 0 { -1.0 } else { 1.0 });
                        s.map_or(0.0, |s| (-s / dot(out, edge(neighbour).2)).max(0.0))
                    }
                    _ => 0.0,
                };
                if remnant - shortened - width <= relief.minimum_remnant + 1e-6
                    && base_arcs[neighbour].is_none()
                {
                    // The notch's bottom runs on to the neighbouring edge
                    // instead (also at exactly the minimum).
                    let vertex = if side == 0 {
                        base[i]
                    } else {
                        base[(i + 1) % n]
                    };
                    // Where the notch bottom (relief depth past the bend
                    // start) meets the neighbouring edge's line.
                    let (u, uj) = (edge(i).1, edge(neighbour).1);
                    let from = at(
                        i,
                        if side == 0 { st.t0 } else { st.t1 },
                        st.setback - notch_depth,
                    );
                    let w = cross2(sub(vertex, from), uj) / cross2(u, uj);
                    let mut p = add(from, scale(u, w));
                    let round = relief.shape == ReliefShape::Round;
                    // A flange on the neighbouring edge reaching past the
                    // notch's bottom line. Flanges are built in spec order:
                    // the later one meets the corner the earlier one left.
                    let t_p = dot(sub(p, edge(neighbour).0), uj);
                    let reached = strips[neighbour].filter(|o| {
                        if side == 0 {
                            o.t1 > t_p + EPS
                        } else {
                            o.t0 < t_p - EPS
                        }
                    });
                    if let Some(o) = reached {
                        let (aj, _, nj, lj) = edge(neighbour);
                        let full = if side == 0 {
                            o.t1 > lj - EPS
                        } else {
                            o.t0 < EPS
                        };
                        let (k, m) = (root[i].unwrap(), root[neighbour].unwrap());
                        let plain = |k: usize| {
                            spec.flanges.iter().any(|f| f.id == folds[k].flange) && !has_lip(k)
                        };
                        if !full || dot(u, uj).abs() > EPS || !plain(k) || !plain(m) {
                            return unsupported(
                                "a short remnant next to a flange on the neighbouring edge",
                            );
                        }
                        // Where a line along this edge through `x` meets the
                        // neighbour's bend start line.
                        let onto = |x: P2| {
                            add(x, scale(u, (o.setback - dot(sub(x, aj), nj)) / dot(u, nj)))
                        };
                        let end_t = if side == 0 { st.t0 } else { st.t1 };
                        // How far past `x` toward the corner the neighbour's
                        // bend starts: this strip's end must stay clear of it.
                        let past =
                            |x: P2| dot(sub(onto(x), x), u) * if side == 0 { -1.0 } else { 1.0 };
                        if m < k {
                            // Built first, the neighbour keeps its corner as
                            // next to a flange with other settings: its bend
                            // ends where both bend start lines cross, its
                            // flange at the vertex. This strip ends square,
                            // the base along its bend start line, no notch.
                            let q = onto(at(i, end_t, st.setback));
                            if past(at(i, end_t, st.setback)) < -EPS {
                                return unsupported(
                                    "its end reaching into the neighbouring flange's bend",
                                );
                            }
                            sides[neighbour][1 - side] = End::Corner {
                                bend: dot(sub(q, aj), uj),
                                flange: if side == 0 { lj } else { 0.0 },
                                slope: 0.0,
                                short: None,
                                shape: CornerShape::TrimToBend,
                                run_on: RunOn::No,
                                turn: (q, q, q),
                                miter: None,
                            };
                            if transition != Transition::None {
                                return unsupported(
                                    "a bend transition at a short remnant next to a flange",
                                );
                            }
                            sides[i][side] = End::Open {
                                bend: end_t,
                                flat: Flat::Slanted,
                            };
                            continue;
                        }
                        // Built later, it is cut back to the notch's bottom
                        // line, which runs on to its bend start.
                        p = onto(from);
                        // A round relief's quarter circle too.
                        let clear = if round { relief.width / 2.0 } else { 0.0 };
                        if past(from) < clear - EPS {
                            return unsupported(
                                "its notch reaching into the neighbouring flange's bend",
                            );
                        }
                        let t = dot(sub(p, aj), uj);
                        let other = strips[neighbour].as_mut().unwrap();
                        if side == 0 {
                            other.t1 = t;
                        } else {
                            other.t0 = t;
                        }
                        cut[neighbour][1 - side] = true;
                    }
                    sides[i][side] = End::RunOn {
                        corner: None,
                        p,
                        round,
                    };
                    continue;
                }
                // The notch must stay clear of a bend along the neighbouring
                // edge (where its start line crosses this edge).
                if let Some(other) = strips[neighbour] {
                    let (a, u, _, _) = edge(i);
                    let (aj, _, nj, _) = edge(neighbour);
                    let cross = (other.setback - dot(sub(a, aj), nj)) / dot(u, nj);
                    let reach = if side == 0 {
                        st.t0 - width
                    } else {
                        st.t1 + width
                    };
                    let into = if side == 0 {
                        reach < cross + EPS
                    } else {
                        reach > cross - EPS
                    };
                    // ...and that bend reaches the notch's depth from the corner.
                    let deep = notch_depth - st.setback;
                    let near = if side == 0 {
                        other.t1 > edge(neighbour).3 - deep
                    } else {
                        other.t0 < deep
                    };
                    if dot(u, nj).abs() > EPS && into && near {
                        return unsupported(
                            "a relief notch running into the neighbouring flange's bend",
                        );
                    }
                }
                sides[i][side] = End::Relief;
                continue;
            }
            if let Some(arc) = base_arcs[neighbour] {
                // A rounded corner: the relief's notch line (its depth past
                // the bend start) runs on until it meets the arc.
                let (a, u, ni, _) = edge(i);
                let vertex = if side == 0 { a } else { base[(i + 1) % n] };
                let r = sub(vertex, arc.center);
                let tangent = if arc.ccw {
                    [-r[1], r[0]]
                } else {
                    [r[1], -r[0]]
                };
                let radius = dot(r, r).sqrt();
                // The arc runs from `from` to `to`.
                let (from, to) = if side == 0 {
                    (base[neighbour], vertex)
                } else {
                    (vertex, base[(neighbour + 1) % n])
                };
                let along = |q: P2| {
                    let (x, y) = (sub(from, arc.center), sub(q, arc.center));
                    let t = cross2(x, y).atan2(dot(x, y)) * if arc.ccw { 1.0 } else { -1.0 };
                    if t < -EPS {
                        t + std::f64::consts::TAU
                    } else {
                        t
                    }
                };
                // Where a line through `o` along `d` meets the circle, on the
                // far side (`sign` 1) or the near one (-1).
                let meet = |o: P2, d: P2, sign: f64| {
                    let w = sub(o, arc.center);
                    let b = dot(w, d);
                    let disc = b * b - dot(w, w) + radius * radius;
                    (disc > 0.0).then(|| add(o, scale(d, -b + sign * disc.sqrt())))
                };
                // Along the edge, away from the strip.
                let (out, end_t) = if side == 0 {
                    (scale(u, -1.0), st.t0)
                } else {
                    (u, st.t1)
                };
                if cross2(tangent, u).abs() > 1e-6 * radius || dot(tangent, u) < 0.0 {
                    // An arc at an angle to the edge: an open end, the bend
                    // start line cut off where it meets the arc (the crossing
                    // nearest the vertex).
                    let o = at(i, end_t, st.setback);
                    let p = [1.0, -1.0]
                        .into_iter()
                        .filter_map(|sign| meet(o, out, sign))
                        .filter(|&p| along(p) <= along(to) + EPS)
                        .min_by(|p, q| dist(*p, o).total_cmp(&dist(*q, o)));
                    let Some(p) = p else {
                        return unsupported("its bend start runs past the neighbouring arc");
                    };
                    if transition == Transition::Arc {
                        return unsupported("an arc bend transition next to an arc");
                    }
                    sides[i][side] = open_end(dot(sub(p, a), u), end_t, side, p);
                    continue;
                }
                let depth = st.setback - notch_depth;
                let far = end_t + if side == 0 { -width } else { width };
                // The remnant: from the notch's far side to the arc, along the
                // bend start line. The notch runs on to the arc when it is
                // no longer than the minimum.
                let remnant = meet(at(i, far, st.setback), out, 1.0)
                    .map_or(0.0, |q| dot(sub(q, at(i, far, st.setback)), out));
                let round = relief.shape == ReliefShape::Round;
                let (corner, p) = if remnant <= relief.minimum_remnant + 1e-6 {
                    // A round notch keeps its quarter circle on the strip's
                    // side.
                    (None, meet(at(i, end_t, depth), out, 1.0))
                } else if round {
                    // Its whole half circle, then on from its far side at the
                    // centre line.
                    let c = at(i, far, depth + relief.width / 2.0);
                    (Some(c), meet(c, ni, 1.0))
                } else {
                    (Some(at(i, far, depth)), meet(at(i, far, depth), ni, 1.0))
                };
                let Some(p) = p.filter(|&p| along(p) <= along(to) + EPS) else {
                    return unsupported("running into an arc too small for its relief");
                };
                sides[i][side] = End::RunOn { corner, p, round };
                continue;
            }
            // A strip's free end next to an edge that isn't square: the bend
            // start line meets that edge at P. Where two flanges meet, any
            // angle.
            let square_edges = dot(edge(neighbour).1, edge(i).1).abs() <= EPS;
            let open = || -> Result<End, SheetError> {
                let (a, u, ni, _) = edge(i);
                let (aj, uj, _, lj) = edge(neighbour);
                // a + t u + s0 ni = aj + w uj.
                let det = -cross2(u, uj);
                let rhs = sub(aj, add(a, scale(ni, st.setback)));
                let (t, w) = (-cross2(rhs, uj) / det, cross2(u, rhs) / det);
                if !(w > EPS && w < lj - EPS) {
                    return Err(SheetError::Unsupported(format!(
                        "flange on edge {i}: its bend start runs past the neighbouring edge"
                    )));
                }
                Ok(open_end(
                    t,
                    if side == 0 { 0.0 } else { length },
                    side,
                    add(aj, scale(uj, w)),
                ))
            };
            let Some(other) = strips[neighbour] else {
                if !square_edges {
                    sides[i][side] = open()?;
                }
                continue;
            };
            let other_length = edge(neighbour).3;
            let meets = if side == 0 {
                other.t1 > other_length - EPS
            } else {
                other.t0 < EPS
            };
            if !meets {
                if !square_edges {
                    sides[i][side] = open()?;
                }
                continue;
            }
            // Two flanges meet here: a trim-to-bend corner with a gap.
            let (k, m) = (root[i].unwrap(), root[neighbour].unwrap());
            let contour = |k: usize| spec.contours.iter().any(|c| c.id == folds[k].flange);
            // One flange feature with two edges: same settings.
            // A hem has no spec flange, so it never has the same settings.
            let spec_of = |k: usize| spec.flanges.iter().find(|f| f.id == folds[k].flange);
            // Copies of one contour along several edges (`id`, `id@e`).
            let feature = |k: usize| {
                let id = &folds[k].flange;
                let copies = id.contains('@')
                    || spec
                        .contours
                        .iter()
                        .any(|c| c.id.starts_with(&format!("{id}@")));
                copies.then(|| id.split('@').next().unwrap().to_string())
            };
            let one_contour =
                contour(k) && contour(m) && feature(k).is_some() && feature(k) == feature(m);
            let same = one_contour
                || match (spec_of(k), spec_of(m)) {
                    (Some(a), Some(b)) => {
                        a.angle == b.angle
                            && a.bend_radius == b.bend_radius
                            && a.k_factor == b.k_factor
                            && a.length == b.length
                            && a.datum == b.datum
                            && a.position == b.position
                            && a.side == b.side
                    }
                    _ => false,
                };
            // A corner seam between flanges of separate features: they meet
            // as alike ones do (measured on 90-degree flanges that differ in
            // length, datum or side).
            let seamed = match (spec_of(k), spec_of(m)) {
                (Some(a), Some(b)) => spec
                    .corner_seams
                    .iter()
                    .find(|c| c.flanges.contains(&a.id) && c.flanges.contains(&b.id))
                    .map(|c| (c, a, b)),
                _ => None,
            };
            if let Some((_, a, b)) = seamed {
                if same {
                    return Err(SheetError::Invalid(format!(
                        "corner seam {} / {}: one feature, already seamed",
                        a.id, b.id
                    )));
                }
                if a.angle != 90.0
                    || b.angle != 90.0
                    || a.bend_radius != b.bend_radius
                    || a.k_factor != b.k_factor
                    || a.position != b.position
                    || has_lip(k)
                    || has_lip(m)
                {
                    return unsupported(
                        "a corner seam but between 90-degree flanges of the same radius, K-factor and position, without lips",
                    );
                }
            }
            let seamed = seamed.map(|(c, _, _)| c);
            let same = same || seamed.is_some();
            let gap = seamed
                .and_then(|c| c.gap)
                .or(relief.miter_gap)
                .unwrap_or(spec.thickness);
            // A lip on one of the two flanges continues its flange's corner
            // edge. Flanges alike with lips alike on both (one feature on
            // two edges, as for one flange) are mitered at the corner like a
            // contour along both edges; other lips on both don't unfold.
            let lips = has_lip(k) || has_lip(m);
            let alike = same && !contour(k) && !contour(m) && has_lip(k) && has_lip(m) && {
                let child = |k: usize| folds.iter().position(|f| f.parent == Some(k));
                let (mut a, mut b) = (Some(k), Some(m));
                loop {
                    match (a, b) {
                        (None, None) => break st.end == other.end,
                        (Some(x), Some(y)) => {
                            let (f, g) = (&folds[x], &folds[y]);
                            if f.angle != g.angle
                                || f.radius != g.radius
                                || f.k_factor != g.k_factor
                                || f.s0 != g.s0
                                || f.bend_allowance != g.bend_allowance
                                || f.up != g.up
                            {
                                break false;
                            }
                            (a, b) = (child(x), child(y));
                        }
                        _ => break false,
                    }
                }
            };
            if same && has_lip(k) && has_lip(m) && !(contour(k) && contour(m)) && !alike {
                return unsupported("lips or hems on both flanges at a corner, unless both flanges and their lips are alike");
            }
            // Two hems meeting fold, but don't unfold.
            if spec_of(k).is_none() && spec_of(m).is_none() && !contour(k) && !contour(m) {
                return unsupported(if contour(k) || contour(m) {
                    "a contour flange meeting a contour flange or a hem at a corner"
                } else {
                    "two hems meeting at a corner"
                });
            }
            let (a, u, ni, _) = edge(i);
            let (aj, uj, nj, _) = edge(neighbour);
            // The bend ends on the line square to its edge through the point
            // where both bend start lines cross.
            let det = ni[0] * nj[1] - ni[1] * nj[0];
            let (ci, cj) = (st.setback + dot(ni, a), other.setback + dot(nj, aj));
            let q = [
                (ci * nj[1] - cj * ni[1]) / det,
                (ni[0] * cj - nj[0] * ci) / det,
            ];
            let bend = dot(sub(q, a), u);
            let shape = relief.corner.map_or(CornerShape::TrimToBend, |c| c.shape);
            if !same {
                // Different settings: two features, no seam. The
                // one built first (flanges in spec order, then contours, then
                // hems) keeps its flange to the vertex; the other's edge ends
                // where the first one's bend starts. Both bends stop at that
                // line. Folds
                // are in build order. A contour built second runs its whole
                // profile over its shortened edge.
                // Two features' corner gets no relief (round and weld
                // settings change nothing).
                // A contour's own bends are fine either way; lips on a flange
                // aren't built here.
                // Lips on one of the two follow their flange (on the
                // first-built and on the second-built one).
                let vertex_t = if side == 0 { 0.0 } else { length };
                sides[i][side] = End::Corner {
                    bend,
                    flange: if k < m { vertex_t } else { bend },
                    slope: 0.0,
                    short: None,
                    shape: CornerShape::TrimToBend,
                    run_on: RunOn::No,
                    turn: (q, q, q),
                    miter: None,
                };
                continue;
            }
            // The flange's inner face stops half the miter gap from the
            // corner's vertical symmetry plane, cut square to the sheet. In
            // the flat that is a line.
            let vertex = if side == 0 { a } else { base[(i + 1) % n] };
            let (into_i, into_j) = if side == 0 {
                (u, scale(uj, -1.0))
            } else {
                (scale(u, -1.0), uj)
            };
            let b = add(into_i, into_j);
            let b = scale(b, 1.0 / dot(b, b).sqrt());
            let plane = if dot([b[1], -b[0]], into_i) > 0.0 {
                [b[1], -b[0]]
            } else {
                [-b[1], b[0]]
            };
            let f = &folds[k];
            let s_end = f.s0 + f.bend_allowance;
            // The inner face, seen from above, at flat s past the bend.
            let inner = |s: f64| f.s0 + f.radius * f.angle.sin() + (s - s_end) * f.angle.cos();
            let half_gap = gap / 2.0;
            let t_at = |s: f64| {
                (half_gap - dot(sub(a, vertex), plane) - inner(s) * dot(ni, plane)) / dot(u, plane)
            };
            let (mut flange, mut slope) = (t_at(s_end), t_at(s_end + 1.0) - t_at(s_end));
            // An overlapping seam (square corners of 90-degree flanges):
            // the bends end as for a miter;
            // the overlapping flange runs on to the other's outer face, the
            // other stops the gap short of the overlapping one's inner face.
            let seam = seamed.map_or(relief.seam.unwrap_or_default(), |c| c.seam);
            if seam != Seam::Miter {
                let right = square_edges
                    && (f.angle - std::f64::consts::FRAC_PI_2).abs() < EPS
                    && shape == CornerShape::TrimToBend
                    && !lips;
                if !right {
                    return unsupported(
                        "an overlapping seam needs a square corner of 90-degree flanges without lips or corner reliefs",
                    );
                }
                // The neighbour's inner face, past its edge (a 90-degree bend).
                let theirs = other.setback + folds[m].radius;
                let past = if (side == 1) == (seam == Seam::Overlap) {
                    theirs + spec.thickness
                } else {
                    theirs - gap
                };
                flange = if side == 0 { -past } else { length + past };
                slope = 0.0;
                // A gap past the bend radius stops a flange short of its
                // bend's end (not measured).
                if if side == 0 {
                    flange > bend + EPS
                } else {
                    flange < bend - EPS
                } {
                    return unsupported(
                        "an overlapping seam whose gap stops a flange short of its bend",
                    );
                }
            }
            // A flange that starts short of the bend's end takes the bend
            // with it; the base keeps its corner at q.
            let short = if side == 0 {
                flange > bend + EPS
            } else {
                flange < bend - EPS
            };
            let (bend, short) = if short {
                (flange, Some(q))
            } else {
                (bend, None)
            };
            // A round or square relief. Larger than the diagonal of where the
            // bends overlap, it
            // is centred where the bend centre lines cross (a square of 6.25
            // is, one of 6.2 isn't, the diagonal being 6.22). Smaller, each
            // bend runs on to its flange's end and the relief sits on the
            // bisector through q: its inscribed circle through q, until its
            // corner (a circle's radius) reaches the circle through q and both
            // flanges' first corners; then at that circle's centre. Clipped
            // once per corner.
            let mut run_on = RunOn::No;
            let mut turn = (q, q, q);
            let mut shape = shape;
            if matches!(shape, CornerShape::Round | CornerShape::Square) {
                let size = relief
                    .corner
                    .and_then(|c| c.size)
                    .unwrap_or(4.0 * spec.thickness);
                let nn = edge(neighbour).2;
                let (ba, bb) = (folds[k].bend_allowance, folds[m].bend_allowance);
                // The point x past edge i and y past its neighbour (outward).
                let det = ni[0] * nn[1] - ni[1] * nn[0];
                let point = |x: f64, y: f64| {
                    add(
                        vertex,
                        [(x * nn[1] - y * ni[1]) / det, (ni[0] * y - nn[0] * x) / det],
                    )
                };
                let crossing = point(st.setback + ba / 2.0, other.setback + bb / 2.0);
                let h = size / 2.0;
                let clip_at = |centre: P2| {
                    if shape == CornerShape::Round {
                        Clip::Disc { c: centre, r: h }
                    } else {
                        // Its diagonal on the bisector of the two edge normals
                        // (along the edges at a square corner).
                        let d = add(ni, nn);
                        let d = scale(d, 1.0 / dot(d, d).sqrt());
                        let r = std::f64::consts::FRAC_1_SQRT_2;
                        let (x, y) = (
                            [r * (d[0] + d[1]), r * (d[1] - d[0])],
                            [r * (d[0] - d[1]), r * (d[0] + d[1])],
                        );
                        let mut q: Vec<P2> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
                            .iter()
                            .map(|&(a, b)| add(centre, add(scale(x, a * h), scale(y, b * h))))
                            .collect();
                        if signed_area(&q) < 0.0 {
                            q.reverse();
                        }
                        Clip::Profile(Ring {
                            arcs: vec![None; q.len()],
                            points: q,
                        })
                    }
                };
                let overlap = [
                    point(st.setback, other.setback),
                    point(st.setback + ba, other.setback),
                    point(st.setback, other.setback + bb),
                    point(st.setback + ba, other.setback + bb),
                ];
                let reach = overlap
                    .iter()
                    .map(|&p| dist(p, crossing))
                    .fold(0.0, f64::max);
                let clip = if let Some(kept) = short {
                    // The flanges end short of q, where the base keeps its
                    // corner: the relief is centred where the bend centre
                    // lines cross once it covers that corner, else dropped.
                    if shape == CornerShape::Round && h > dist(crossing, kept) + EPS {
                        Some(clip_at(crossing))
                    } else if h > reach + EPS {
                        return unsupported(
                            "a square corner relief where the flange starts short of its bend",
                        );
                    } else {
                        shape = CornerShape::TrimToBend;
                        None
                    }
                } else if h > reach + EPS {
                    Some(clip_at(crossing))
                } else {
                    let mirror = |p: P2| sub(p, scale(plane, 2.0 * dot(sub(p, vertex), plane)));
                    let first = at(i, flange, s_end);
                    // Capped at the centre of the circle through q and both
                    // flanges' first corners; at a corner that isn't square,
                    // of flanges that aren't at 90 degrees, at the middle of
                    // those two corners.
                    let cap = if square_edges || slope.abs() < EPS {
                        circumcentre(q, first, mirror(first))
                    } else {
                        scale(add(first, mirror(first)), 0.5)
                    };
                    let corner = if shape == CornerShape::Round {
                        h
                    } else {
                        h * std::f64::consts::SQRT_2
                    };
                    let capped = corner > dist(cap, q);
                    let mut centre = if capped { cap } else { sub(q, scale(b, h)) };
                    let mut clip = clip_at(centre);
                    let covers = |clip: &Clip, p: P2| clip.contains(p) || clip.on_boundary(p);
                    // Each bend runs on to its flange's end line (sloped at
                    // flanges other than 90 degrees); those lines cross at X
                    // on the symmetry plane (behind q, inside the base, under
                    // 90 degrees).
                    let run = |s: f64| at(i, flange + slope * (s - s_end), s);
                    let side_of = |s: f64| dot(sub(run(s), vertex), plane);
                    let (g0, g1) = (side_of(st.setback), side_of(s_end));
                    let x = run(st.setback + (s_end - st.setback) * g0 / (g0 - g1));
                    // Where the end line, from the flange's first corner on,
                    // first meets the relief: the outline turns onto the
                    // relief there.
                    let meet = |clip: &Clip| {
                        if covers(clip, first) {
                            // Capped: through the flange's first corner.
                            return Some(first);
                        }
                        clip.crossings(first, x)
                            .ok()
                            .and_then(|ts| ts.first().copied())
                            .map(|t| add(first, scale(sub(x, first), t)))
                    };
                    // Whether it does, a square by its inscribed circle.
                    if !capped && meet(&Clip::Disc { c: centre, r: h }).is_none() {
                        // It doesn't reach the end line: it is centred on X,
                        // in the flat only (X is in the base under 90
                        // degrees).
                        centre = x;
                        clip = clip_at(x);
                    }
                    // A capped one the end line doesn't reach before the bend
                    // start line (a sharp corner of sloped flanges): its end
                    // takes another line there (not built).
                    let Some(met) = meet(&clip) else {
                        return unsupported(if capped {
                            "a capped corner relief its flanges' end lines don't reach"
                        } else {
                            "a corner relief over its flange's first corner"
                        });
                    };
                    turn = (centre, met, met);
                    // A capped square short of the flange's first corner:
                    // the bend's end would follow the miter there (not built).
                    if capped && shape == CornerShape::Square && !covers(&clip, first) {
                        return unsupported(
                            "a square corner relief too large to touch q and too small to reach past the bend",
                        );
                    }
                    run_on = if capped { RunOn::Both } else { RunOn::Flat };
                    Some(clip)
                };
                if let (Some(clip), 1) = (clip, side) {
                    corner_clips.push((clip, run_on != RunOn::Flat));
                }
            }
            // A lip at a relieved corner continues its flange's edge as at a
            // trimmed one; at a sloped
            // flange edge it does something else.
            if lips && (slope.abs() > EPS || short.is_some()) {
                return unsupported("a lip at a corner whose flange edge isn't square to its edge");
            }
            sides[i][side] = End::Corner {
                bend,
                flange,
                slope,
                short,
                shape,
                run_on,
                turn,
                miter: (one_contour || alike).then_some(Miter {
                    vertex,
                    plane,
                    half_gap,
                }),
            };
        }
    }
    // An arc between two strips whose notches cross before reaching it is
    // gone: both notches end where they cross.
    let mut crossed = vec![false; n];
    for j in 0..n {
        let (i, k) = ((j + n - 1) % n, (j + 1) % n);
        let (Some(arc), Some(si), Some(sk)) = (base_arcs[j], strips[i], strips[k]) else {
            continue;
        };
        let (
            End::RunOn {
                corner: ci, p: pi, ..
            },
            End::RunOn {
                corner: ck, p: pk, ..
            },
        ) = (sides[i][1], sides[k][0])
        else {
            continue;
        };
        let (a1, a2) = (
            ci.unwrap_or(at(i, si.t1, si.setback - notch_depth)),
            ck.unwrap_or(at(k, sk.t0, sk.setback - notch_depth)),
        );
        // Past each other along the arc: pk comes before pi.
        let (x, y) = (sub(pi, arc.center), sub(pk, arc.center));
        if (cross2(x, y) > 0.0) == arc.ccw {
            continue;
        }
        let (d1, d2) = (sub(pi, a1), sub(pk, a2));
        let x = add(a1, scale(d1, cross2(sub(a2, a1), d2) / cross2(d1, d2)));
        sides[i][1] = End::RunOn {
            corner: ci,
            p: x,
            round: false,
        };
        sides[k][0] = End::RunOn {
            corner: ck,
            p: x,
            round: false,
        };
        crossed[j] = true;
    }
    // Two notches running on into the same corner: both end where their
    // bottoms cross, whichever flange comes first (straight and tear
    // reliefs at a square corner; a round one's run-on meets the other's
    // quarter circle).
    for i in 0..n {
        let k = (i + 1) % n;
        let (End::RunOn { corner: None, .. }, End::RunOn { corner: None, .. }) =
            (sides[i][1], sides[k][0])
        else {
            continue;
        };
        if relief.shape == ReliefShape::Round || dot(edge(i).1, edge(k).1).abs() > EPS {
            return Err(SheetError::Unsupported(format!(
                "flanges on edges {i} and {k}: short remnants on both sides of the corner"
            )));
        }
        let (si, sk) = (strips[i].unwrap(), strips[k].unwrap());
        let (a1, d1) = (at(i, si.t1, si.setback - notch_depth), edge(i).1);
        let (a2, d2) = (at(k, sk.t0, sk.setback - notch_depth), edge(k).1);
        let x = add(a1, scale(d1, cross2(sub(a2, a1), d2) / cross2(d1, d2)));
        let p = End::RunOn {
            corner: None,
            p: x,
            round: false,
        };
        (sides[i][1], sides[k][0]) = (p, p);
    }

    // Span of a strip's bend and of its flange along the edge.
    let span = |i: usize, st: &Strip| {
        let pick = |side: usize, t: f64| match sides[i][side] {
            End::Corner { bend, flange, .. } => (bend, flange),
            End::Open { bend, .. } | End::Square { bend, .. } => (bend, t),
            _ => (t, t),
        };
        let ((b0, f0), (b1, f1)) = (pick(0, st.t0), pick(1, st.t1));
        ([b0, b1], [f0, f1])
    };
    let bend_end = |i: usize| {
        let f = &folds[root[i].unwrap()];
        f.s0 + f.bend_allowance
    };
    // The fold maps alone (no cells yet), for where a mitered chain ends.
    let shape = Layout {
        flat: FlatPattern {
            outline: Vec::new(),
            outline_arcs: Vec::new(),
            holes: Vec::new(),
            cutouts: Vec::new(),
            bends: Vec::new(),
            punches: Vec::new(),
            forms: Vec::new(),
        },
        cells: Vec::new(),
        folds: folds.clone(),
        holes: Vec::new(),
        embosses: Vec::new(),
    };
    let maps = crate::fold::Model::new(spec, &shape);
    // A mitered chain's end along edge i at flat s, past its first bend:
    // where the miter's face of the strip there lies half the gap from the
    // symmetry plane. Every fold turns about an axis along the edge, so
    // moving along the edge moves the point along it: linear in t.
    let mitered = |i: usize, side: usize, s: f64| -> Option<f64> {
        let End::Corner { miter: Some(m), .. } = sides[i][side] else {
            return None;
        };
        // The first bend keeps its own corner (trim, relief, weld, tear).
        if s < bend_end(i) - EPS {
            return None;
        }
        let k = (0..folds.len())
            .rev()
            .find(|&k| fold_ends[k].0 == i && folds[k].s0 <= s + EPS)?;
        let f = &folds[k];
        let zone = if s < f.s0 + f.bend_allowance - EPS {
            Zone::Bend(k)
        } else {
            Zone::Flange(k)
        };
        // Whichever face reaches it first: the shorter strip.
        let t = |z: f64| {
            let p = maps.map(zone, at(i, 0.0, s), z);
            (m.half_gap - dot(sub([p[0], p[1]], m.vertex), m.plane)) / dot(edge(i).1, m.plane)
        };
        let (a, b) = (t(0.0), t(spec.thickness));
        Some(if side == 0 { a.max(b) } else { a.min(b) })
    };
    // Where the flange ends along the edge at flat s (past its bend).
    let flange_at = |i: usize, side: usize, t: f64, s: f64| {
        if let Some(t) = mitered(i, side, s) {
            return t;
        }
        match sides[i][side] {
            End::Corner { flange, slope, .. } => flange + slope * (s - bend_end(i)),
            _ => t,
        }
    };
    // A bend past the first in a mitered chain: its end, one curve from the
    // bend's start to its end (on the face the miter takes mid-bend).
    let mitered_cut = |i: usize, side: usize, k: usize| -> Option<(P2, P2, Cut)> {
        let End::Corner { miter: Some(m), .. } = sides[i][side] else {
            return None;
        };
        let f = &folds[k];
        let (origin, u, normal, _) = edge(i);
        let radius = f.radius + f.k_factor * spec.thickness;
        let face = |z: f64, s: f64| {
            let p = maps.map(Zone::Bend(k), at(i, 0.0, s), z);
            (m.half_gap - dot(sub([p[0], p[1]], m.vertex), m.plane)) / dot(u, m.plane)
        };
        let mid = f.s0 + f.bend_allowance / 2.0;
        let z = if (face(0.0, mid) - mitered(i, side, mid)?).abs() < 1e-9 {
            0.0
        } else {
            spec.thickness
        };
        let ss = [f.s0, mid, f.s0 + f.bend_allowance];
        let rows = ss.map(|s| {
            let (sin, cos) = ((s - f.s0) / radius).sin_cos();
            [1.0, cos, sin, face(z, s)]
        });
        let cut = Cut {
            origin,
            u,
            n: normal,
            s0: f.s0,
            radius,
            c: solve3(rows),
        };
        Some((cut.point(ss[0]), cut.point(ss[2]), cut))
    };
    // The later folds of edge i's chain, outward.
    let later = |i: usize| -> Vec<usize> {
        (0..folds.len())
            .filter(|&k| fold_ends[k].0 == i && folds[k].parent.is_some())
            .collect()
    };
    // One face takes the miter through each bend (where they'd swap mid-bend
    // the end isn't one curve: not built).
    for i in 0..n {
        for side in 0..2 {
            for k in later(i) {
                let Some((_, _, cut)) = mitered_cut(i, side, k) else {
                    continue;
                };
                let f = &folds[k];
                for j in 1..8 {
                    let s = f.s0 + f.bend_allowance * j as f64 / 8.0;
                    if (cut.t(s) - mitered(i, side, s).unwrap()).abs() > 1e-7 {
                        return Err(SheetError::Unsupported(format!(
                            "{}: a mitered bend whose faces swap at the miter",
                            f.flange
                        )));
                    }
                }
            }
        }
    }
    // A flange's free tip, (f0, end) to (f1, end), with its corners cut.
    let tip =
        |path: &mut Path, i: usize, f0: f64, f1: f64, end: f64, cut: Option<CornerCut>| match cut {
            None => {
                path.line(at(i, flange_at(i, 0, f0, end), end));
                path.line(at(i, flange_at(i, 1, f1, end), end));
            }
            Some(CornerCut::Chamfer { distance: d }) => {
                for p in [
                    at(i, f0, end - d),
                    at(i, f0 + d, end),
                    at(i, f1 - d, end),
                    at(i, f1, end - d),
                ] {
                    path.line(p);
                }
            }
            Some(CornerCut::Round { radius: r }) => {
                for (from, to, center) in [
                    (
                        at(i, f0, end - r),
                        at(i, f0 + r, end),
                        at(i, f0 + r, end - r),
                    ),
                    (
                        at(i, f1 - r, end),
                        at(i, f1, end - r),
                        at(i, f1 - r, end - r),
                    ),
                ] {
                    path.line(from);
                    *path.arcs.last_mut().unwrap() = Some(Curve::Arc(Arc { center, ccw: true }));
                    path.line(to);
                }
            }
        };
    // The corner cut on a flange's tip, where nothing sits on it.
    let tip_cut = |k: usize| {
        if has_lip(k) {
            return None;
        }
        spec.flanges
            .iter()
            .find(|f| f.id == folds[k].flange)
            .and_then(|f| f.tip_corners)
    };
    let last_fold = |i: usize| (0..folds.len()).find(|&k| fold_ends[k].0 == i && !has_lip(k));
    for (k, fold) in folds.iter().enumerate() {
        let Some(cut) = tip_cut(k) else { continue };
        let i = fold_ends[k].0;
        let ([_, _], [f0, f1]) = span(i, &strips[i].unwrap());
        let size = match cut {
            CornerCut::Round { radius } => radius,
            CornerCut::Chamfer { distance } => distance,
        };
        let square = sides[i]
            .iter()
            .all(|e| matches!(e, End::Free | End::Relief));
        if !square {
            return Err(SheetError::Unsupported(format!(
                "{}: tip corners where its ends aren't square to it",
                fold.flange
            )));
        }
        let straight = fold_ends[k].1 - (fold.s0 + fold.bend_allowance);
        if !(size > 0.0 && 2.0 * size < f1 - f0 - EPS && size < straight - EPS) {
            return Err(SheetError::Invalid(format!(
                "{}: tip corners must be smaller than half its width and its straight part",
                fold.flange
            )));
        }
    }

    // A relief notch next to a strip end, walked in the edge's direction:
    // `into` from the base edge down to the notch end and back up at `t`
    // (start side), or the reverse (end side). Returns the path pieces.
    if relief.shape == ReliefShape::Round && relief.depth < relief.width / 2.0 - EPS {
        return Err(SheetError::Invalid(
            "a round relief's depth must be at least half its width".into(),
        ));
    }
    let notch = |path: &mut Path, i: usize, st: &Strip, start: bool, outline: bool| {
        if relief.shape == ReliefShape::Tear && outline {
            // Torn, not cut: the flat shows no notch.
            path.line(at(i, if start { st.t0 } else { st.t1 }, 0.0));
            return;
        }
        let (w, depth) = if relief.shape == ReliefShape::Tear {
            (TEAR, TEAR)
        } else {
            (relief.width, relief.depth)
        };
        let deepest = st.setback - depth;
        let (near, far) = if start {
            (st.t0 - w, st.t0)
        } else {
            (st.t1, st.t1 + w)
        };
        match relief.shape {
            ReliefShape::Straight | ReliefShape::Tear => {
                if start {
                    path.line(at(i, near, 0.0));
                }
                path.line(at(i, near, deepest));
                path.line(at(i, far, deepest));
                if !start {
                    path.line(at(i, far, 0.0));
                }
            }
            ReliefShape::Round => {
                // Straight sides down to the centre line, a half circle beyond.
                let centre = deepest + w / 2.0;
                if start {
                    path.line(at(i, near, 0.0));
                }
                path.line(at(i, near, centre));
                path.arc(
                    at(i, far, centre),
                    at(i, near + w / 2.0, centre),
                    at(i, near + w / 2.0, deepest),
                );
                if !start {
                    path.line(at(i, far, 0.0));
                }
            }
        }
    };

    // A round relief's notch running on to the corner: a quarter circle of
    // half its width from the strip's side down to the notch bottom
    //. From the bottom up at a start, the
    // other way at an end.
    // With `half` the notch stops short of the corner: its whole half circle,
    // from (or to) its far side at the centre line.
    let quarter = |path: &mut Path, i: usize, st: &Strip, start: bool, half: bool| {
        let (r, depth) = (relief.width / 2.0, st.setback - relief.depth);
        let h = r * std::f64::consts::FRAC_1_SQRT_2;
        if half {
            let (t, c) = if start {
                (st.t0, st.t0 - r)
            } else {
                path.line(at(i, st.t1, depth + r));
                (st.t1 + 2.0 * r, st.t1 + r)
            };
            path.arc(at(i, t, depth + r), at(i, c, depth + r), at(i, c, depth));
        } else if start {
            path.line(at(i, st.t0 - r, depth));
            path.arc(
                at(i, st.t0, depth + r),
                at(i, st.t0 - r, depth + r),
                at(i, st.t0 - r + h, depth + r - h),
            );
        } else {
            path.line(at(i, st.t1, depth + r));
            path.arc(
                at(i, st.t1 + r, depth),
                at(i, st.t1 + r, depth + r),
                at(i, st.t1 + r - h, depth + r - h),
            );
        }
    };
    // A tear's flat closes the corner after edge i: the two flange edges run
    // on until they meet.
    let meet = |i: usize| {
        let line = |e: usize, side: usize| {
            let f = |s: f64| at(e, flange_at(e, side, 0.0, s), s);
            (f(bend_end(e)), f(bend_end(e) + 1.0))
        };
        let ((p, p2), (q, q2)) = (line(i, 1), line((i + 1) % n, 0));
        let (d, e) = (sub(p2, p), sub(q2, q));
        add(p, scale(d, cross2(sub(q, p), e) / cross2(d, e)))
    };
    // Welded and torn corners end each bend in a straight cut from where both
    // bends start to its flange's corner.
    let welded = |e: End| {
        matches!(
            e,
            End::Corner {
                shape: CornerShape::LinearWeld | CornerShape::Tear,
                ..
            }
        )
    };

    // The outline, edge by edge: base corner, then the strip's profile.
    let torn = relief.shape == ReliefShape::Tear;
    let mut outline = Path::default();
    for i in 0..n {
        // A corner, or an open end on either side, cuts the base vertex off.
        let prev = (i + n - 1) % n;
        let corner_start = strips[i].is_some()
            && matches!(
                sides[i][0],
                End::Corner { .. } | End::Open { .. } | End::Square { .. } | End::RunOn { .. }
            )
            || strips[prev].is_some()
                && matches!(
                    sides[prev][1],
                    End::Open { .. } | End::Square { .. } | End::RunOn { .. }
                );
        if !corner_start {
            outline.line(base[i]);
            if strips[i].is_none() && base_arcs[i].is_some() {
                *outline.arcs.last_mut().unwrap() = base_arcs[i].map(Curve::Arc);
            }
        }
        let Some(st) = strips[i] else { continue };
        let ([b0, b1], [f0, f1]) = span(i, &st);
        let be = bend_end(i);
        // An open end's corner at the vertex, rounded (an arc transition).
        let radius = folds[root[i].unwrap()].radius;
        match sides[i][0] {
            End::Free => {}
            End::Open { flat, .. } => {
                match flat {
                    Flat::Slanted => outline.line(at(i, b0, st.setback)),
                    Flat::Square => outline.line(base[i]),
                    Flat::Fillet => fillet(&mut outline, base[i], edge(prev).1, edge(i).2, radius),
                }
                // Drawn square, the end runs straight on into the flange's
                // (a fillet may reach past the bend there).
                if flat == Flat::Slanted {
                    outline.line(at(i, st.t0, be));
                }
            }
            End::Square { q, .. } => {
                for p in [q, at(i, b0, st.setback), at(i, b0, be), at(i, st.t0, be)] {
                    outline.line(p);
                }
            }
            End::Relief => notch(&mut outline, i, &st, true, true),
            // A tear's slit next to an arc: torn, not cut, as `notch` does.
            End::RunOn {
                corner: Some(_), ..
            } if torn => outline.line(at(i, st.t0, 0.0)),
            End::RunOn { corner, p, round } => {
                outline.line(p);
                if let Some(c) = corner {
                    outline.line(c);
                }
                if round {
                    quarter(&mut outline, i, &st, true, corner.is_some());
                } else {
                    outline.line(at(i, st.t0, st.setback - notch_depth));
                }
            }
            End::Corner {
                short,
                shape,
                run_on,
                turn,
                ..
            } => {
                if let Some(q) = short {
                    outline.line(q);
                }
                let pieces = match shape {
                    // The flange edges run on until they meet.
                    CornerShape::Tear => vec![meet((i + n - 1) % n)],
                    CornerShape::LinearWeld => vec![at(i, b0, st.setback), at(i, f0, be)],
                    // Bends run on: their ends cross on the symmetry plane.
                    _ if run_on != RunOn::No => vec![turn.0, turn.1, turn.2, at(i, f0, be)],
                    _ => vec![at(i, b0, st.setback), at(i, b0, be), at(i, f0, be)],
                };
                for p in pieces {
                    outline.line(p);
                }
            }
        }
        for k in later(i) {
            if let Some((p, q, cut)) = mitered_cut(i, 0, k) {
                outline.line(p);
                outline.curve(q, Curve::Cut(cut));
            }
        }
        tip(
            &mut outline,
            i,
            f0,
            f1,
            st.end,
            last_fold(i).and_then(tip_cut),
        );
        for k in later(i).into_iter().rev() {
            if let Some((p, q, cut)) = mitered_cut(i, 1, k) {
                outline.line(q);
                outline.curve(p, Curve::Cut(cut));
            }
        }
        match sides[i][1] {
            End::Free => {}
            End::Open { flat, .. } => {
                if flat == Flat::Slanted {
                    outline.line(at(i, st.t1, be));
                }
                let next = (i + 1) % n;
                match flat {
                    Flat::Slanted => outline.line(at(i, b1, st.setback)),
                    Flat::Square => outline.line(base[next]),
                    Flat::Fillet => fillet(
                        &mut outline,
                        base[next],
                        scale(edge(i).2, -1.0),
                        edge(next).1,
                        radius,
                    ),
                }
                if flat != Flat::Fillet {
                    *outline.arcs.last_mut().unwrap() = base_arcs[next].map(Curve::Arc);
                }
            }
            End::Square { q, .. } => {
                for p in [at(i, st.t1, be), at(i, b1, be), at(i, b1, st.setback), q] {
                    outline.line(p);
                }
                *outline.arcs.last_mut().unwrap() = base_arcs[(i + 1) % n].map(Curve::Arc);
            }
            End::Relief => notch(&mut outline, i, &st, false, true),
            End::RunOn {
                corner: Some(_), ..
            } if torn => {
                outline.line(at(i, st.t1, 0.0));
                *outline.arcs.last_mut().unwrap() = base_arcs[(i + 1) % n].map(Curve::Arc);
            }
            End::RunOn { corner, p, round } => {
                if round {
                    quarter(&mut outline, i, &st, false, corner.is_some());
                } else {
                    outline.line(at(i, st.t1, st.setback - notch_depth));
                }
                if let Some(c) = corner {
                    outline.line(c);
                }
                outline.line(p);
                if !crossed[(i + 1) % n] {
                    *outline.arcs.last_mut().unwrap() = base_arcs[(i + 1) % n].map(Curve::Arc);
                }
            }
            End::Corner {
                shape,
                run_on,
                turn,
                ..
            } => {
                let pieces = match shape {
                    CornerShape::Tear => vec![meet(i)],
                    CornerShape::LinearWeld => vec![at(i, f1, be), at(i, b1, st.setback)],
                    _ if run_on != RunOn::No => vec![at(i, f1, be), turn.2, turn.1, turn.0],
                    _ => vec![at(i, f1, be), at(i, b1, be), at(i, b1, st.setback)],
                };
                for p in pieces {
                    outline.line(p);
                }
            }
        }
    }
    outline.close();
    simplify(&mut outline.pts, &mut outline.arcs, true);

    // The base cell: the outline with every flanged span cut back to its bend
    // start, keeping the junctions with the strips as vertices.
    let mut base_cell = Path::default();
    for i in 0..n {
        let prev = (i + n - 1) % n;
        if let (Some(p), End::Corner { short: Some(q), .. }) = (strips[prev], sides[i][0]) {
            // The neighbour's bend start, the corner, then this one's.
            base_cell.line(at(prev, span(prev, &p).0[1], p.setback));
            base_cell.line(q);
        }
        let depth = |st: &Strip| st.setback - notch_depth;
        // From the previous strip's notch on past this vertex.
        let prev_runs = if let (
            Some(p),
            End::RunOn {
                corner,
                p: q,
                round,
            },
        ) = (strips[prev], sides[prev][1])
        {
            base_cell.line(at(prev, span(prev, &p).0[1], p.setback));
            if round {
                quarter(&mut base_cell, prev, &p, false, corner.is_some());
            } else {
                base_cell.line(at(prev, p.t1, depth(&p)));
            }
            if let Some(c) = corner {
                base_cell.line(c);
            }
            base_cell.line(q);
            true
        } else {
            false
        };
        // From the previous strip's end trimmed to its bend, along its bend
        // start to the neighbouring edge.
        let prev_square = if let (Some(p), End::Square { q, .. }) = (strips[prev], sides[prev][1]) {
            base_cell.line(at(prev, span(prev, &p).0[1], p.setback));
            base_cell.line(q);
            true
        } else {
            false
        };
        match (strips[i], strips[prev], sides[i][0], sides[prev][1]) {
            // Along the neighbouring edge to a notch that runs on.
            (Some(st), _, End::RunOn { corner, p, round }, _) => {
                base_cell.line(p);
                if let Some(c) = corner {
                    base_cell.line(c);
                }
                if round {
                    quarter(&mut base_cell, i, &st, true, corner.is_some());
                } else {
                    base_cell.line(at(i, st.t0, depth(&st)));
                }
                base_cell.line(at(i, span(i, &st).0[0], st.setback));
            }
            (Some(st), prev_strip, end, prev_end) if !matches!(end, End::Relief) => {
                // A square end next to a trimmed corner: both bend starts.
                if let (
                    Some(p),
                    (End::Open { .. }, End::Corner { .. }) | (End::Corner { .. }, End::Open { .. }),
                ) = (prev_strip, (prev_end, end))
                {
                    base_cell.line(at(prev, span(prev, &p).0[1], p.setback));
                }
                if let End::Square { q, .. } = end {
                    base_cell.line(q);
                }
                base_cell.line(at(i, span(i, &st).0[0], st.setback));
            }
            _ if prev_runs || prev_square => {}
            (_, Some(p), _, end) if !matches!(end, End::Relief) => {
                base_cell.line(at(prev, span(prev, &p).0[1], p.setback));
            }
            _ => base_cell.line(base[i]),
        }
        if strips[i].is_none() && base_arcs[i].is_some() && !crossed[i] {
            *base_cell.arcs.last_mut().unwrap() = base_arcs[i].map(Curve::Arc);
        }
        let Some(st) = strips[i] else { continue };
        let [b0, b1] = span(i, &st).0;
        if matches!(sides[i][0], End::Relief) {
            notch(&mut base_cell, i, &st, true, false);
            base_cell.line(at(i, b0, st.setback));
        }
        if matches!(sides[i][1], End::Relief) {
            base_cell.line(at(i, b1, st.setback));
            base_cell.line(at(i, st.t1, st.setback));
            notch(&mut base_cell, i, &st, false, false);
        }
    }
    base_cell.close();

    // Strips: per fold its bend, then its flat flange up to the next bend
    // (a lip on its tip) or its own end. A flange wider than its bend (at a
    // corner) keeps the bend's ends as vertices on its bottom edge.
    let mut cells = vec![Cell {
        name: "base".into(),
        polygon: base_cell.pts,
        arcs: base_cell.arcs,
        holes: Vec::new(),
        inner: Vec::new(),
        zone: Zone::Base,
    }];
    for (k, fold) in folds.iter().enumerate() {
        let (i, end) = fold_ends[k];
        let st = strips[i].unwrap();
        let ([b0, b1], [f0, f1]) = span(i, &st);
        // A lip spans its parent flange, not the base bend.
        let (b0, b1) = if fold.parent.is_some() {
            (f0, f1)
        } else {
            (b0, b1)
        };
        // Where the bend ends along the edge at its far side: at a welded
        // corner, where the flange ends.
        let corner = |side: usize| welded(sides[i][side]);
        let open = |side: usize| matches!(sides[i][side], End::Open { .. });
        // A bend that runs on to its flange's end, start and far side: the
        // cell where it folds so, its bend line wherever it does.
        let runs = |side: usize, how: &[RunOn]| {
            fold.parent.is_none()
                && matches!(sides[i][side], End::Corner { run_on, .. } if how.contains(&run_on))
        };
        let run_on = |side: usize| runs(side, &[RunOn::Both]);
        let line_on = |side: usize| runs(side, &[RunOn::Both, RunOn::Flat]);
        let (q0, q1) = (b0, b1);
        let (b0, b1) = (
            if run_on(0) { f0 } else { b0 },
            if run_on(1) { f1 } else { b1 },
        );
        // A run-on bend's cell: from q along the corner's diagonal to where
        // both bends' ends cross, then along its own end.
        let start = |side: usize, t: f64| match sides[i][side] {
            End::Corner { turn, .. } if run_on(side) => {
                if side == 0 {
                    vec![turn.0, turn.1, turn.2]
                } else {
                    vec![turn.2, turn.1, turn.0]
                }
            }
            _ => vec![at(i, t, fold.s0)],
        };
        let e0 = if corner(0) || run_on(0) {
            f0
        } else if open(0) {
            st.t0
        } else {
            b0
        };
        let e1 = if corner(1) || run_on(1) {
            f1
        } else if open(1) {
            st.t1
        } else {
            b1
        };
        let fold_bend_end = fold.s0 + fold.bend_allowance;
        let flange_end = folds
            .iter()
            .find(|f| f.parent == Some(k))
            .map_or(end, |child| child.s0);
        // (t, s) is a mirrored frame, so this order is counter-clockwise in x, y.
        // A flange narrower than its bend at a corner leaves its ends as
        // vertices on the bend's far edge.
        let mut bend = start(0, b0);
        bend.push(at(i, e0, fold_bend_end));
        for t in [f0, f1] {
            if t > e0 + EPS && t < e1 - EPS {
                bend.push(at(i, t, fold_bend_end));
            }
        }
        bend.push(at(i, e1, fold_bend_end));
        bend.extend(start(1, b1));
        let mut bend_arcs = vec![None; bend.len()];
        // A later bend of a mitered chain: its ends are curves.
        if fold.parent.is_some() {
            let (a, c) = (mitered_cut(i, 0, k), mitered_cut(i, 1, k));
            if a.is_some() || c.is_some() {
                let (p0, q0, c0) = match a {
                    Some((p, q, cut)) => (p, q, Some(Curve::Cut(cut))),
                    None => (at(i, b0, fold.s0), at(i, e0, fold_bend_end), None),
                };
                let (p1, q1, c1) = match c {
                    Some((p, q, cut)) => (p, q, Some(Curve::Cut(cut))),
                    None => (at(i, b1, fold.s0), at(i, e1, fold_bend_end), None),
                };
                bend = vec![p0, q0, q1, p1];
                bend_arcs = vec![c0, None, c1, None];
            }
        }
        // A bend run on past q keeps q, where the base ends, as a vertex.
        if run_on(1) {
            bend.push(at(i, q1, fold.s0));
        }
        if run_on(0) {
            bend.push(at(i, q0, fold.s0));
        }
        bend_arcs.resize(bend.len(), None);
        cells.push(Cell {
            name: format!("{}.bend", fold.flange),
            arcs: bend_arcs,
            polygon: bend,
            holes: Vec::new(),
            inner: Vec::new(),
            zone: Zone::Bend(k),
        });
        let mut flange = Path::default();
        let first = fold.parent.is_none();
        let from = |side: usize, t: f64| {
            mitered(i, side, fold_bend_end)
                .filter(|_| !first)
                .unwrap_or(t)
        };
        flange.line(at(i, from(0, f0), fold_bend_end));
        tip(&mut flange, i, f0, f1, flange_end, tip_cut(k));
        flange.line(at(i, from(1, f1), fold_bend_end));
        for t in [e1, e0] {
            if (t - f0).abs() > EPS && (t - f1).abs() > EPS {
                flange.line(at(i, t, fold_bend_end));
            }
        }
        // A bend with nothing straight after it (a rolled hem) has no flange cell.
        if flange_end - fold_bend_end > EPS {
            cells.push(Cell {
                name: fold.flange.clone(),
                arcs: flange.arcs,
                polygon: flange.pts,
                holes: Vec::new(),
                inner: Vec::new(),
                zone: Zone::Flange(k),
            });
        }
        let centre = fold.s0 + fold.bend_allowance / 2.0;
        // At a welded corner the bend line ends halfway along the slanted
        // end; at a torn one, on the flange edge run back over the bend (the
        // flat shows the corner closed).
        // Where the bend line ends at an open end the flat shows square
        // through the vertex: where its centre line meets that outline (the
        // square line, its rounded corner, or the neighbouring edge).
        let drawn = |side: usize| -> Option<f64> {
            let End::Open { flat, .. } = sides[i][side] else {
                return None;
            };
            if flat == Flat::Slanted {
                return None;
            }
            let (a, u, ni, _) = edge(i);
            let (v, d_in, d_out, away) = if side == 0 {
                let along = edge((i + n - 1) % n).1;
                (base[i], along, ni, scale(along, -1.0))
            } else {
                let along = edge((i + 1) % n).1;
                (base[(i + 1) % n], scale(ni, -1.0), along, along)
            };
            let r = folds[root[i].unwrap()].radius;
            let (p, q, c) = if flat == Flat::Fillet {
                fillet_points(v, d_in, d_out, r)
            } else {
                (v, v, v)
            };
            // Where the outline leaves the square line, and joins the
            // neighbouring edge.
            let (on_square, on_edge) = if side == 0 { (q, p) } else { (p, q) };
            if centre >= dot(sub(on_square, a), ni) - EPS {
                return Some(dot(sub(v, a), u));
            }
            let from = add(a, scale(ni, centre));
            let t = cross2(sub(on_edge, from), away) / cross2(u, away);
            if dot(sub(add(from, scale(u, t)), on_edge), away) >= -EPS {
                return Some(t);
            }
            // On the rounded corner: of the circle's two crossings, the one
            // nearer the vertex.
            let h = centre - dot(sub(c, a), ni);
            let w = (r * r - h * h).sqrt();
            let (tc, tv) = (dot(sub(c, a), u), dot(sub(v, a), u));
            Some(if (tc + w - tv).abs() < (tc - w - tv).abs() {
                tc + w
            } else {
                tc - w
            })
        };
        let mid = |side: usize, b: f64, e: f64| {
            if let Some(t) = drawn(side) {
                t
            } else if matches!(
                sides[i][side],
                End::Corner {
                    shape: CornerShape::Tear,
                    ..
                }
            ) {
                flange_at(i, side, e, centre)
            } else {
                (b + e) / 2.0
            }
        };
        // A run-on bend's line runs to its end: the end line, then the turn
        // onto the relief (which then cuts it).
        let ends = |side: usize, f: f64| -> f64 {
            let End::Corner {
                turn: (c, h, v), ..
            } = sides[i][side]
            else {
                return f;
            };
            let p0 = at(i, 0.0, centre);
            let u = sub(at(i, 1.0, centre), p0);
            for (a, b) in [(at(i, f, fold_bend_end), v), (v, h), (h, c)] {
                let d = sub(b, a);
                let den = cross2(u, d);
                if den.abs() < EPS {
                    continue;
                }
                let w = sub(a, p0);
                if (-EPS..=1.0 + EPS).contains(&(cross2(w, u) / den)) {
                    return cross2(w, d) / den;
                }
            }
            flange_at(i, side, f, centre)
        };
        let (l0, l1) = (
            if line_on(0) {
                ends(0, f0)
            } else if let Some(t) = mitered(i, 0, centre).filter(|_| fold.parent.is_some()) {
                t
            } else {
                mid(0, b0, e0)
            },
            if line_on(1) {
                ends(1, f1)
            } else if let Some(t) = mitered(i, 1, centre).filter(|_| fold.parent.is_some()) {
                t
            } else {
                mid(1, b1, e1)
            },
        );
        bends[k].line = [at(i, l0, centre), at(i, l1, centre)];
    }

    // Holes and cut-outs on a flange's face, into the flat: along the edge as
    // they are, up from the bend's end less its outer setback (the flange's
    // face is flat from there to the next bend).
    let mut holes = spec.base.holes.clone();
    let mut cutouts: Vec<Ring> = spec.base.cutouts.iter().map(Ring::from_cutout).collect();
    // A flange's face: [along, up] into the flat, along the edge as they are,
    // up from the bend's end less its outer setback (the face is flat from
    // there to the next bend); and a check that a flat extent across the
    // edge (lowest, highest distance from it) stays on that flat.
    let face = |on: &str, what: &str| {
        let mut found = folds.iter().enumerate().filter(|(_, f)| f.flange == on);
        let (Some((k, f)), None, true) = (
            found.next(),
            found.next(),
            spec.flanges.iter().any(|g| g.id == on),
        ) else {
            return Err(SheetError::Invalid(format!("{what}: no flange {on}")));
        };
        let start = f.s0 + f.bend_allowance;
        let setback = (f.radius + spec.thickness) * (f.angle / 2.0).tan();
        let next = folds
            .iter()
            .filter(|g| g.parent == Some(k))
            .map(|g| g.s0)
            .fold(f64::INFINITY, f64::min);
        let frame = Frame2 {
            o: add(f.edge_point, scale(f.n, start - setback)),
            x: f.u,
            y: f.n,
        };
        let (edge, n) = (f.edge_point, f.n);
        let what = format!("{what} on {on}");
        let on_face = move |ring: &Ring| {
            let (lo, hi) = extent_along(ring, n, edge);
            if lo < start - EPS || hi > next + EPS {
                return Err(SheetError::Unsupported(format!(
                    "{what}: reaching into a bend"
                )));
            }
            Ok(())
        };
        Ok((frame, on_face))
    };
    let embosses = spec
        .base
        .embosses
        .iter()
        .map(|e| {
            let placed = match &e.on {
                Some(on) => e.placed(&face(on, &format!("emboss {}", e.id))?.0),
                None => e.clone(),
            };
            crate::emboss::Formed::new(&placed, spec.thickness, spec.bend_radius)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut punches = Vec::new();
    for p in &spec.base.punches {
        if p.at.is_empty() || p.holes.is_empty() && p.cutouts.is_empty() {
            return Err(SheetError::Invalid(format!(
                "punch {}: needs a shape and a place",
                p.id
            )));
        }
        let (p, on_face) = match &p.on {
            Some(on) => {
                let (frame, on_face) = face(on, &format!("punch {}", p.id))?;
                (p.placed(&frame), Some(on_face))
            }
            None => (p.clone(), None),
        };
        for &at in &p.at {
            let (h, k) = p.stamped(at);
            let k: Vec<Ring> = k.iter().map(Ring::from_cutout).collect();
            if let Some(on_face) = &on_face {
                for ring in h.iter().map(Ring::from_hole).chain(k.iter().cloned()) {
                    on_face(&ring)?;
                }
            }
            holes.extend(h);
            cutouts.extend(k);
            punches.push(FlatPunch {
                id: p.id.clone(),
                center: at,
                angle: p.angle,
            });
        }
    }
    for fc in &spec.face_cuts {
        let (frame, on_face) = face(&fc.on, "face cuts")?;
        // The cuts as one stamp at the frame's origin.
        let placed = Punch {
            id: String::new(),
            on: None,
            at: vec![[0.0, 0.0]],
            angle: 0.0,
            holes: fc.holes.clone(),
            cutouts: fc.cutouts.clone(),
        }
        .placed(&frame);
        let (h, k) = placed.stamped(placed.at[0]);
        let k: Vec<Ring> = k.iter().map(Ring::from_cutout).collect();
        for ring in h.iter().map(Ring::from_hole).chain(k.iter().cloned()) {
            on_face(&ring)?;
        }
        holes.extend(h);
        cutouts.extend(k);
    }

    let mut layout = Layout {
        flat: FlatPattern {
            outline: outline.pts,
            outline_arcs: outline.arcs,
            holes: holes.clone(),
            cutouts,
            bends,
            punches,
            forms: embosses
                .iter()
                .map(|e| FlatForm {
                    id: e.id.clone(),
                    center: e.center,
                    outline: e.ring(),
                })
                .collect(),
        },
        cells,
        folds,
        holes: holes.clone(),
        embosses,
    };
    for (clip, cells) in &corner_clips {
        crate::cut::cut_corner(&mut layout, clip, *cells)?;
    }
    if !spec.folds.is_empty() {
        crate::foldline::apply_folds(&mut layout, spec)?;
    }
    // A hole across the outline is a round notch: a circular cut-out.
    let (across, inside): (Vec<Hole>, Vec<Hole>) = holes
        .into_iter()
        .partition(|h| crate::cut::crosses_outline(&layout.flat, h));
    layout.flat.holes = inside.clone();
    let notches = across.iter().map(|h| {
        let (c, r) = (h.center, h.diameter / 2.0);
        let arc = Some(Curve::Arc(Arc {
            center: c,
            ccw: true,
        }));
        // Two half circles, from an odd angle: a hole centred on an edge
        // would otherwise have its vertices on that edge.
        let (sin, cos) = 1f64.sin_cos();
        Ring {
            points: vec![
                [c[0] + r * cos, c[1] + r * sin],
                [c[0] - r * cos, c[1] - r * sin],
            ],
            arcs: vec![arc, arc],
        }
    });
    layout.flat.cutouts.extend(notches);
    // An emboss's footprint goes into the cells (round ones as holes): the
    // formed sheet starts at its edge.
    let discs: Vec<Hole> = layout.embosses.iter().filter_map(|e| e.hole()).collect();
    for (e, h) in layout
        .embosses
        .iter()
        .filter_map(|e| e.hole().map(|h| (e, h)))
    {
        if crate::cut::crosses_outline(&layout.flat, &h) {
            return Err(SheetError::Invalid(format!(
                "emboss {}: its footprint crosses the outline",
                e.id
            )));
        }
    }
    crate::cut::cut_holes(&mut layout, &[inside.clone(), discs.clone()].concat())?;
    crate::cut::cut_cutouts(&mut layout)?;
    let rings: Vec<(String, Ring)> = layout
        .embosses
        .iter()
        .filter(|e| e.hole().is_none())
        .map(|e| (e.id.clone(), e.ring()))
        .collect();
    for (k, (id, ring)) in rings.iter().enumerate() {
        let others: Vec<Ring> = rings[..k].iter().map(|r| r.1.clone()).collect();
        crate::cut::cut_inner(&mut layout, ring, &discs, &others, &format!("emboss {id}"))?;
    }
    for e in &layout.embosses {
        let held = layout.cells.iter().any(|c| {
            !matches!(c.zone, Zone::Bend(_))
                && match e.hole() {
                    Some(h) => c.holes.contains(&h),
                    None => c
                        .inner
                        .iter()
                        .any(|r| r.points.contains(&e.ring().points[0])),
                }
        });
        if !held {
            return Err(SheetError::Unsupported(format!(
                "emboss {}: only on one flat plate, clear of bends and cut-outs",
                e.id
            )));
        }
    }
    // No cell runs over itself: two features' pieces lie over each other,
    // such as a relief notch running into a neighbouring bend.
    let crossed = layout.cells.iter().find_map(|c| {
        crate::faces::overlap(&Ring {
            points: c.polygon.clone(),
            arcs: c.arcs.clone(),
        })
    });
    if let Some((a, b)) = crossed {
        return Err(SheetError::Unsupported(format!(
            "features overlap in the flat pattern near [{:.3}, {:.3}] and [{:.3}, {:.3}] (a relief or cut running into a bend)",
            a[0], a[1], b[0], b[1]
        )));
    }
    simplify(
        &mut layout.flat.outline,
        &mut layout.flat.outline_arcs,
        true,
    );
    // In the cells only arcs: a straight run's points may be junctions with
    // a neighbouring cell, a hole's edge never is.
    for c in &mut layout.cells {
        simplify(&mut c.polygon, &mut c.arcs, false);
    }
    number_bends(&spec.bend_order, &mut layout)?;
    Ok(layout)
}

/// Numbers the bends in the order they're made: `first` in its order, then
/// the bends furthest from the base, then as the flat pattern lists them.
fn number_bends(first: &[String], layout: &mut Layout) -> Result<(), SheetError> {
    let depth = |name: &str| {
        let mut fold = layout.folds.iter().position(|f| f.flange == name);
        let mut d = 0;
        while let Some(i) = fold {
            d += 1;
            fold = layout.folds[i].parent;
        }
        d
    };
    let mut names: Vec<&str> = Vec::new();
    for b in &layout.flat.bends {
        if !names.contains(&b.flange.as_str()) {
            names.push(&b.flange);
        }
    }
    names.sort_by_key(|n| std::cmp::Reverse(depth(n)));
    for (k, name) in first.iter().enumerate() {
        if first[..k].contains(name) {
            return Err(SheetError::Invalid(format!(
                "bend_order lists {name} twice"
            )));
        }
        let Some(i) = names.iter().position(|n| n == name) else {
            return Err(SheetError::Invalid(format!("bend_order: no bend {name}")));
        };
        let n = names.remove(i);
        names.insert(k, n);
    }
    let names: Vec<String> = names.into_iter().map(String::from).collect();
    for b in &mut layout.flat.bends {
        b.order = 1 + names.iter().position(|n| *n == b.flange).unwrap();
    }
    Ok(())
}

impl FlatPattern {
    /// The blank's bounding box, [lowest, highest] x and y, arcs exact.
    pub fn bounds(&self) -> [P2; 2] {
        let ring = Ring {
            points: self.outline.clone(),
            arcs: self.outline_arcs.clone(),
        };
        let (x, y) = (
            extent_along(&ring, [1.0, 0.0], [0.0, 0.0]),
            extent_along(&ring, [0.0, 1.0], [0.0, 0.0]),
        );
        [[x.0, y.0], [x.1, y.1]]
    }

    /// Length of every cut, mm: the outline, the cut-outs (punch stamps
    /// included) and the holes, measured on the lines and arcs the DXF writes.
    pub fn cut_length(&self) -> f64 {
        let ring_length = |pts: &[P2], curves: &[Option<Curve>]| -> f64 {
            let (pts, arcs) = lines_and_arcs(pts, curves);
            (0..pts.len())
                .map(|k| {
                    let (a, b) = (pts[k], pts[(k + 1) % pts.len()]);
                    match arcs[k] {
                        None => dot(sub(b, a), sub(b, a)).sqrt(),
                        Some(arc) => {
                            dot(sub(a, arc.center), sub(a, arc.center)).sqrt()
                                * arc.sweep(a, b).abs()
                        }
                    }
                })
                .sum()
        };
        ring_length(&self.outline, &self.outline_arcs)
            + self
                .cutouts
                .iter()
                .map(|c| ring_length(&c.points, &c.arcs))
                .sum::<f64>()
            + self
                .holes
                .iter()
                .map(|h| std::f64::consts::PI * h.diameter)
                .sum::<f64>()
    }

    /// Closed cuts the laser starts: the outline, each cut-out and each hole.
    pub fn pierces(&self) -> usize {
        1 + self.cutouts.len() + self.holes.len()
    }

    /// Area of the flat blank: outline minus holes and cut-outs.
    pub fn area(&self) -> f64 {
        path_area(&self.outline, &self.outline_arcs)
            - self.cutouts.iter().map(Ring::area).sum::<f64>()
            - self
                .holes
                .iter()
                .map(|h| std::f64::consts::PI * h.diameter * h.diameter / 4.0)
                .sum::<f64>()
    }
}

fn signed_area(p: &[P2]) -> f64 {
    (0..p.len())
        .map(|i| {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        / 2.0
}
fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
}
fn cross2(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
/// Drops points in the middle of a straight run (e.g. a base corner that a
/// full-width flange strip continues past) or of one circle's arc (a round
/// notch's two halves).
fn simplify(pts: &mut Vec<P2>, arcs: &mut Vec<Option<Curve>>, lines: bool) {
    loop {
        let len = pts.len();
        let drop = (0..len).find(|&k| {
            let prev = (k + len - 1) % len;
            let (p, q, r) = (pts[prev], pts[k], pts[(k + 1) % len]);
            let (d1, d2) = (sub(q, p), sub(r, q));
            let straight = lines
                && arcs[prev].is_none()
                && arcs[k].is_none()
                && (d1[0] * d2[1] - d1[1] * d2[0]).abs() < EPS
                && dot(d1, d2) > 0.0;
            let one_circle = match (arcs[prev], arcs[k]) {
                (Some(Curve::Arc(a)), Some(Curve::Arc(b))) => {
                    dist(a.center, b.center) < EPS && a.ccw == b.ccw
                }
                _ => false,
            };
            straight || one_circle
        });
        let Some(k) = drop else { break };
        pts.remove(k);
        arcs.remove(k);
    }
}

fn dist(a: P2, b: P2) -> f64 {
    dot(sub(a, b), sub(a, b)).sqrt()
}

/// The centre of the circle through three points.
fn circumcentre(a: P2, b: P2, c: P2) -> P2 {
    let (b, c) = (sub(b, a), sub(c, a));
    let d = 2.0 * cross2(b, c);
    let (bb, cc) = (dot(b, b), dot(c, c));
    add(
        a,
        [(c[1] * bb - b[1] * cc) / d, (b[0] * cc - c[0] * bb) / d],
    )
}
fn sub(a: P2, b: P2) -> P2 {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale(a: P2, s: f64) -> P2 {
    [a[0] * s, a[1] * s]
}
fn dot(a: P2, b: P2) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::ReliefShape;
    use crate::spec::{Base, BendPosition, Flange, HeightDatum, Width};

    fn plate(flanges: Vec<Flange>) -> PartSpec {
        PartSpec {
            id: "t".into(),
            thickness: 2.0,
            bend_radius: 2.0,
            k_factor: 0.4,
            bend_table: None,
            relief: None,
            base: Base {
                outline: vec![[0.0, 0.0], [100.0, 0.0], [100.0, 60.0], [0.0, 60.0]],
                holes: vec![Hole {
                    center: [20.0, 30.0],
                    diameter: 8.0,
                }],
                ..Default::default()
            },
            flanges,
            contours: Vec::new(),
            hems: Vec::new(),
            folds: Vec::new(),

            jogs: Vec::new(),
            face_cuts: Vec::new(),
            corner_seams: Vec::new(),
            faces: Vec::new(),
            bends: Vec::new(),
            bend_order: Vec::new(),
        }
    }

    fn flange(edge: usize) -> Flange {
        Flange {
            id: format!("f{edge}"),
            on: EdgeRef::Base(edge),
            width: None,
            angle: 90.0,
            length: 30.0,
            datum: HeightDatum::Outer,
            position: BendPosition::Outside,
            side: Side::Up,
            tip_corners: None,
            bend_radius: None,
            k_factor: None,
        }
    }

    const REACH: f64 = 28.3982297150257; // Expected: 88.398... - 60

    fn has(outline: &[P2], p: P2) -> bool {
        outline
            .iter()
            .any(|q| (q[0] - p[0]).abs() < 1e-9 && (q[1] - p[1]).abs() < 1e-9)
    }

    #[test]
    fn a_plate_is_its_own_flat_pattern() {
        let flat = flat_pattern(&plate(vec![])).unwrap();
        assert_eq!(flat.outline, plate(vec![]).base.outline);
        assert!((flat.area() - (6000.0 - std::f64::consts::PI * 16.0)).abs() < 1e-9);
    }

    #[test]
    fn a_flange_moves_its_edge_out_and_puts_the_bend_line_mid_strip() {
        let flat = flat_pattern(&plate(vec![flange(2)])).unwrap();
        assert!(has(&flat.outline, [100.0, 60.0 + REACH]));
        assert!(has(&flat.outline, [0.0, 60.0 + REACH]));
        // Bend starts 2 mm inside the edge, strip is 4.398 mm wide.
        let centre = 60.0 - 2.0 + 0.5 * std::f64::consts::FRAC_PI_2 * 2.8;
        assert!((flat.bends[0].line[0][1] - centre).abs() < 1e-9);
        assert!(flat.bends[0].up);
    }

    #[test]
    fn a_partial_flange_gets_relief_notches() {
        // Edge 2 runs (100,60) -> (0,60); offsets 30 from (100,60), 20 from (0,60).
        let partial = Flange {
            width: Some(Width {
                offset_start: 30.0,
                offset_end: 20.0,
            }),
            ..flange(2)
        };
        let flat = flat_pattern(&plate(vec![partial])).unwrap();
        // Expected: notches 2 wide reaching to y = 57 (bend start 58, depth 1).
        for p in [
            [72.0, 60.0],
            [72.0, 57.0],
            [70.0, 57.0],
            [70.0, 60.0 + REACH],
            [20.0, 60.0 + REACH],
            [20.0, 57.0],
            [18.0, 57.0],
            [18.0, 60.0],
        ] {
            assert!(
                has(&flat.outline, p),
                "{p:?} missing from {:?}",
                flat.outline
            );
        }
        let bend = &flat.bends[0];
        assert!((bend.line[0][0] - 70.0).abs() < 1e-9 && (bend.line[1][0] - 20.0).abs() < 1e-9);
    }

    #[test]
    fn a_lip_on_a_flange_adds_its_reach_and_keeps_or_flips_the_bend_direction() {
        for (side, lip_up) in [(Side::Up, true), (Side::Down, false)] {
            let lip = Flange {
                id: "lip".into(),
                on: EdgeRef::Tip("f2".into()),
                length: 15.0,
                side,
                ..flange(2)
            };
            let flat = flat_pattern(&plate(vec![flange(2), lip])).unwrap();
            // Expected: 101.7965 mm for either lip direction.
            let top = flat.outline.iter().map(|p| p[1]).fold(f64::MIN, f64::max);
            assert!((top - 101.79645943005142).abs() < 1e-9, "{top}");
            assert_eq!(flat.bends[1].up, lip_up);
            let centre = 60.0 + REACH - 2.0 + 0.5 * std::f64::consts::FRAC_PI_2 * 2.8;
            assert!((flat.bends[1].line[0][1] - centre).abs() < 1e-9);
        }
        let orphan = Flange {
            on: EdgeRef::Tip("nope".into()),
            ..flange(2)
        };
        assert!(matches!(
            flat_pattern(&plate(vec![orphan])),
            Err(SheetError::Invalid(_))
        ));
    }

    #[test]
    fn two_flanges_meeting_at_a_corner_trim_to_their_bends() {
        let flat = flat_pattern(&plate(vec![flange(1), flange(2)])).unwrap();
        // Expected (trim to bend, miter gap T): each bend runs to the other's
        // bend start (x 98 / y 58), each flange stops T/sqrt(2) short of the
        // other's inner face (x 98.586 / y 58.586).
        let g = 2.0 - std::f64::consts::SQRT_2;
        let bend_end = 98.0 + 4.3982297150257;
        for p in [
            [0.0, 0.0],
            [100.0 + REACH, 0.0],
            [100.0 + REACH, 58.0 + g],
            [bend_end, 58.0 + g],
            [bend_end, 58.0],
            [98.0, 58.0],
            [98.0, 58.0 + 4.3982297150257],
            [98.0 + g, 58.0 + 4.3982297150257],
            [98.0 + g, 60.0 + REACH],
            [0.0, 60.0 + REACH],
        ] {
            assert!(
                has(&flat.outline, p),
                "{p:?} missing from {:?}",
                flat.outline
            );
        }
        assert_eq!(flat.outline.len(), 10);
    }
    fn round_relief(depth: f64) -> PartSpec {
        let mut part = plate(vec![Flange {
            width: Some(Width {
                offset_start: 30.0,
                offset_end: 20.0,
            }),
            ..flange(2)
        }]);
        part.relief = Some(crate::spec::Relief {
            shape: ReliefShape::Round,
            width: 2.0,
            depth,
            minimum_remnant: 4.0,
            miter_gap: None,
            corner: None,
            seam: None,
            transition: None,
        });
        part
    }

    #[test]
    fn a_round_relief_ends_in_a_half_circle() {
        // Expected: half circles of radius w/2 whose bottom is `depth` past
        // the bend start (y 58): centres at y 58 for depth 1, y 56 for depth 3.
        for (depth, cy) in [(1.0, 58.0), (3.0, 56.0)] {
            let flat = flat_pattern(&round_relief(depth)).unwrap();
            let arcs: Vec<_> = flat
                .outline_arcs
                .iter()
                .flatten()
                .filter_map(|c| c.arc())
                .collect();
            assert_eq!(arcs.len(), 2);
            for (cx, p) in [(71.0, [72.0, cy]), (19.0, [20.0, cy])] {
                assert!(arcs
                    .iter()
                    .any(|a| (a.center[0] - cx).abs() < 1e-9 && (a.center[1] - cy).abs() < 1e-9));
                assert!(
                    has(&flat.outline, p),
                    "{p:?} missing from {:?}",
                    flat.outline
                );
            }
            // Area: the straight notch of the same depth, minus the two corners
            // the half circle leaves: 2 * (w^2 - pi w^2 / 4) / 2.
            let mut straight = round_relief(depth);
            straight.relief.as_mut().unwrap().shape = ReliefShape::Straight;
            let corners = 2.0 * (4.0 - std::f64::consts::PI) / 2.0;
            let diff = flat.area() - flat_pattern(&straight).unwrap().area();
            assert!((diff - corners).abs() < 1e-9, "{diff} vs {corners}");
        }
        assert!(matches!(
            flat_pattern(&round_relief(0.5)),
            Err(SheetError::Invalid(_))
        ));
    }

    #[test]
    fn a_round_relief_at_a_short_remnant_runs_on_with_a_quarter_circle() {
        // Expected (remnant-round-start / -end): the notch runs on to the
        // corner and keeps a quarter circle of radius w/2 on the strip's side,
        // so it leaves one corner of the straight run-on, r^2 - pi r^2 / 4,
        // plus the two of the half circle at the flange's other end.
        for (start, end) in [(3.0, 30.0), (30.0, 4.0)] {
            let mut round = round_relief(1.5);
            round.flanges[0].width = Some(Width {
                offset_start: start,
                offset_end: end,
            });
            let mut straight = round.clone();
            straight.relief.as_mut().unwrap().shape = ReliefShape::Straight;
            let flat = flat_pattern(&round).unwrap();
            assert_eq!(flat.outline_arcs.iter().flatten().count(), 2);
            let diff = flat.area() - flat_pattern(&straight).unwrap().area();
            let corner = 3.0 * (1.0 - std::f64::consts::PI / 4.0);
            assert!((diff - corner).abs() < 1e-9, "{diff} vs {corner}");
        }
    }

    #[test]
    fn short_remnants_on_both_sides_of_a_corner_cross() {
        // Both notches run on into the corner at (0, 0): their bottoms (3
        // inside both edges) end where they cross, in either order (a fuzz
        // find, once refused). A round relief's run-on meets the other's
        // quarter circle instead (not built).
        let part = |relief: &str, flipped: bool| -> PartSpec {
            let (f0, f3) = (
                r#"{ "id": "f0", "on": { "base": 0 }, "width": { "offset_start": 5, "offset_end": 10 },
                     "angle": 90, "length": 20, "datum": "outer", "position": "outside", "side": "up" }"#,
                r#"{ "id": "f3", "on": { "base": 3 }, "width": { "offset_start": 10, "offset_end": 5 },
                     "angle": 90, "length": 20, "datum": "outer", "position": "outside", "side": "up" }"#,
            );
            let flanges = if flipped { [f3, f0] } else { [f0, f3] };
            serde_json::from_str(&format!(
                r#"{{ "id": "a", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "relief": {{ "shape": "{relief}", "width": 2, "depth": 1, "minimum_remnant": 4 }},
                "base": {{ "outline": [[0,0],[70,0],[70,40],[0,40]] }},
                "flanges": [{}, {}] }}"#,
                flanges[0], flanges[1]
            ))
            .unwrap()
        };
        let flat = flat_pattern(&part("straight", false)).unwrap();
        assert!(flat.outline.iter().any(|p| dist(*p, [3.0, 3.0]) < 1e-9));
        let flipped = flat_pattern(&part("straight", true)).unwrap();
        assert!((flat.area() - flipped.area()).abs() < 1e-9);
        assert_eq!(
            flat_pattern(&part("round", false)).unwrap_err(),
            SheetError::Unsupported(
                "flanges on edges 3 and 0: short remnants on both sides of the corner".into()
            )
        );
    }

    #[test]
    fn a_short_remnant_next_to_a_flange_follows_the_build_order() {
        // f on edge 2 ends 3 from the corner (100, 60), g runs along edge 1
        // (remnant-next-flange and -rev). g built first keeps its
        // corner: its bend stops at f's bend start line (y 58), f ends
        // square with no notch. Built later, g is cut back to f's notch
        // bottom (y 57).
        let f = Flange {
            width: Some(Width {
                offset_start: 3.0,
                offset_end: 30.0,
            }),
            ..flange(2)
        };
        let g = flange(1);
        let later = flat_pattern(&plate(vec![f.clone(), g.clone()])).unwrap();
        let first = flat_pattern(&plate(vec![g, f])).unwrap();
        assert!(has(&later.outline, [97.0, 57.0]) && !has(&later.outline, [97.0, 58.0]));
        assert!(has(&first.outline, [97.0, 58.0]) && !has(&first.outline, [97.0, 57.0]));
        assert!(first.area() > later.area());
    }

    #[test]
    fn a_round_relief_runs_on_to_an_arc_with_a_quarter_circle() {
        // Expected (flange-into-arcs-round / -relief): flat areas
        // 8126.315972323419 and 8124.599157630599, two corners r^2 - pi r^2 / 4
        // apart (r 2). With a minimum remnant of 1 the notch stops short of
        // the arc (flange-into-arcs-round-short): 8151.63898178028.
        let part = |shape: &str, minimum: f64| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "a", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "relief": {{ "shape": "{shape}", "width": 4, "depth": 3, "minimum_remnant": {minimum} }},
                "base": {{ "outline": [[10,0],[90,0],[100,10],[100,50],[90,60],[10,60],[0,50],[0,10]],
                           "bulges": [0,0.41421356237309503,0,0.41421356237309503,0,0.41421356237309503,0,0.41421356237309503] }},
                "flanges": [{{ "id": "f", "on": {{ "base": 4 }}, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }}] }}"#
            ))
            .unwrap()
        };
        let area = |shape: &str, minimum: f64| flat_pattern(&part(shape, minimum)).unwrap().area();
        let round = area("round", 4.0);
        assert!((round - 8126.315972323419).abs() < 1e-9, "{round}");
        let corners = round - area("straight", 4.0);
        assert!((corners - 2.0 * (4.0 - std::f64::consts::PI)).abs() < 1e-9);
        let short = area("round", 1.0);
        assert!((short - 8151.63898178028).abs() < 1e-9, "{short}");
    }

    #[test]
    fn flanges_next_to_arcs() {
        let part = |outline: &str, bulges: &str, relief: &str, flange: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{ "id": "a", "thickness": 2, "bend_radius": 2, "k_factor": 0.4, "relief": {{ {relief} }},
                "base": {{ "outline": {outline}, "bulges": {bulges} }},
                "flanges": [{{ "id": "f", {flange}, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }}] }}"#
            ))
            .unwrap()
        };
        let q = 0.41421356237309503;
        let rounded = |relief: &str, flange: &str| {
            part(
                "[[10,0],[90,0],[100,10],[100,50],[90,60],[10,60],[0,50],[0,10]]",
                &format!("[0,{q},0,{q},0,{q},0,{q}]"),
                relief,
                flange,
            )
        };
        let area = |p: &PartSpec| flat_pattern(p).unwrap().area();
        let volume = |p: &PartSpec| crate::fold::exact_mass(p).unwrap().volume;
        // Each flat: the rounded plate plus the strip past its edge.
        let plate = 6000.0 - (4.0 - std::f64::consts::PI) * 100.0;
        // (the bend starts 2 inside the edge; 26 of flange past it).
        let strip = 24.0 + std::f64::consts::FRAC_PI_2 * 2.8;

        // A tear: its 0.01 slit up to the arc, torn in the flat; in 3D 2 x
        // 0.01 by ~2.01 through the sheet.
        let tear = rounded(
            r#""shape": "tear", "width": 2, "depth": 1, "minimum_remnant": 4"#,
            r#""on": { "base": 4 }"#,
        );
        assert!((area(&tear) - (plate + 80.0 * strip)).abs() < 1e-9);
        let v = volume(&tear);
        assert!((v - 16422.220368103895).abs() < 1e-6, "{v}");
        // ...whose remnant to the arc is short: it runs on at its depth.
        let tear = rounded(
            r#""shape": "tear", "width": 2, "depth": 1, "minimum_remnant": 8"#,
            r#""on": { "base": 4 }"#,
        );
        assert!((area(&tear) - 8169.547399).abs() < 1e-6);

        // A partial flange keeps its remnant to an arc, however short: two
        // plain 4 x 5 notches.
        let partial = rounded(
            r#""width": 4, "depth": 3, "minimum_remnant": 4"#,
            r#""on": { "base": 4 }, "width": { "offset_start": 6, "offset_end": 6 }"#,
        );
        assert!((area(&partial) - (plate + 68.0 * strip - 40.0)).abs() < 1e-9);

        // An arc at an angle to the edge (convex or concave): the bend
        // start line is cut off where it meets the arc.
        for (bulge, expected) in [(0.2, 9082.682068), (-0.2, 8596.963875)] {
            let angled = part(
                "[[0,0],[100,0],[100,60],[0,60]]",
                &format!("[0,{bulge},0,0]"),
                r#""width": 2, "depth": 1, "minimum_remnant": 4"#,
                r#""on": { "base": 0 }"#,
            );
            let f = flat_pattern(&angled).unwrap();
            assert!((f.area() - expected).abs() < 1e-6, "{}", f.area());
            // On the arc (centre (28, 30) or (172, 30), radius 78) at y = 2.
            let meet = 28.0 + 5300f64.sqrt();
            let meet = if bulge > 0.0 { meet } else { 200.0 - meet };
            assert!(f.outline.iter().any(|p| dist(*p, [meet, 2.0]) < 1e-9));
        }
    }

    #[test]
    fn a_tear_relief_at_a_short_remnant_runs_on_at_the_tear_depth() {
        // Expected: the notch runs on to the corner 0.01 past the bend
        // start, so the flat loses the remnant strip. The flat areas of
        // remnant-tear-start (offsets 3 / 30) and relief-tear (30 / 20) differ
        // by 476.739905155436.
        let torn = |start: f64, end: f64| {
            let mut part = round_relief(1.0);
            part.relief.as_mut().unwrap().shape = ReliefShape::Tear;
            part.flanges[0].width = Some(Width {
                offset_start: start,
                offset_end: end,
            });
            flat_pattern(&part).unwrap().area()
        };
        let diff = torn(3.0, 30.0) - torn(30.0, 20.0);
        assert!((diff - 476.739905155436).abs() < 1e-9, "{diff}");
    }

    #[test]
    fn the_cells_tile_the_flat_outline_exactly() {
        let partial = Flange {
            width: Some(Width {
                offset_start: 30.0,
                offset_end: 20.0,
            }),
            ..flange(2)
        };
        let lip = |parent: &str| Flange {
            id: format!("lip-{parent}"),
            on: EdgeRef::Tip(parent.into()),
            length: 15.0,
            ..flange(2)
        };
        for flanges in [
            vec![flange(2)],
            vec![partial.clone()],
            vec![partial, lip("f2")],
            vec![flange(0), flange(2), lip("f0"), lip("f2")],
            vec![flange(1), flange(2), lip("f1"), lip("f2")],
        ] {
            let l = layout(&plate(flanges)).unwrap();
            let cells: f64 = l.cells.iter().map(|c| signed_area(&c.polygon)).sum();
            let outline = signed_area(&l.flat.outline);
            assert!(l.cells.iter().all(|c| signed_area(&c.polygon) > 0.0));
            assert!((cells - outline).abs() < 1e-9, "{cells} vs {outline}");
        }
    }

    fn with_holes(holes: &[([f64; 2], f64)]) -> PartSpec {
        let mut part = plate(vec![flange(2)]);
        part.base.holes = holes
            .iter()
            .map(|&(center, diameter)| Hole { center, diameter })
            .collect();
        part
    }

    #[test]
    fn holes_anywhere_in_the_flat_cut_the_cells_exactly() {
        // Edge 2 at y = 60 carries the flange: its bend runs y 58 .. 62.398.
        for holes in [
            vec![([50.0, 59.0], 6.0)], // across the bend start
            vec![([50.0, 60.2], 8.0)], // across the whole bend
            vec![([30.0, 60.2], 3.0)], // inside the bend
            vec![([70.0, 75.0], 8.0)], // on the flange
            vec![
                ([20.0, 30.0], 8.0),
                ([50.0, 60.2], 8.0),
                ([80.0, 58.0], 3.0),
            ],
        ] {
            let l = layout(&with_holes(&holes)).unwrap();
            let cells: f64 = l
                .cells
                .iter()
                .map(|c| {
                    path_area(&c.polygon, &c.arcs)
                        - c.holes
                            .iter()
                            .map(|h| std::f64::consts::PI * h.diameter * h.diameter / 4.0)
                            .sum::<f64>()
                })
                .sum();
            assert!(l.cells.iter().all(|c| path_area(&c.polygon, &c.arcs) > 0.0));
            assert!(
                (cells - l.flat.area()).abs() < 1e-9,
                "{holes:?}: {cells} vs {}",
                l.flat.area()
            );
        }
        // Across the whole bend, the hole splits the strip into two pieces.
        let l = layout(&with_holes(&[([50.0, 60.2], 8.0)])).unwrap();
        assert_eq!(
            l.cells.iter().filter(|c| c.zone == Zone::Bend(0)).count(),
            2
        );
    }

    #[test]
    fn a_corner_seam_joins_flanges_of_separate_features() {
        // Flanges 30 (edge 1) and 20 (edge 2) long. Expected flat areas and
        // volumes: no seam 9498.123836 / 19094.265362 (by build order),
        // miter 9470.726866 / 19039.471423, overlap 9550.123836 /
        // 19198.265362, reverse overlap 9510.123836 / 19118.265362.
        let part = |seam: Option<Seam>| {
            let mut p = plate(vec![
                flange(1),
                Flange {
                    length: 20.0,
                    ..flange(2)
                },
            ]);
            p.base.holes.clear();
            p.corner_seams = seam
                .map(|seam| crate::spec::CornerSeam {
                    flanges: ["f1".into(), "f2".into()],
                    seam,
                    gap: None,
                })
                .into_iter()
                .collect();
            p
        };
        for (seam, area, volume) in [
            (None, 9498.123836, 19094.265362),
            (Some(Seam::Miter), 9470.726866, 19039.471423),
            (Some(Seam::Overlap), 9550.123836, 19198.265362),
            (Some(Seam::ReverseOverlap), 9510.123836, 19118.265362),
        ] {
            let a = flat_pattern(&part(seam)).unwrap().area();
            assert!((a - area).abs() < 1e-6, "{seam:?}: {a}");
            let v = crate::fold::exact_mass(&part(seam)).unwrap().volume;
            assert!((v - volume).abs() < 1e-6, "{seam:?}: {v}");
        }
    }

    #[test]
    fn bend_transitions_at_an_open_end() {
        // A flange on edge 2 ending at the obtuse vertex (30, 60). Expected
        // flat areas: none 7089.075195, intersection 7087.876080 (the strip
        // square through the vertex), arc 7087.893057 (that corner rounded
        // R 2), trim to bend 7086.876080; volume 14222.446846, trimmed
        // 14217.734457 (the bend's slanted wedge gone).
        let part = |transition: Option<Transition>| {
            let mut p = plate(vec![flange(2)]);
            p.base.holes.clear();
            p.base.outline = vec![[0.0, 0.0], [100.0, 0.0], [100.0, 60.0], [30.0, 60.0]];
            p.relief = Some(crate::spec::Relief {
                transition,
                ..p.relief()
            });
            p
        };
        for (t, area, volume) in [
            (None, 7089.075195, 14222.446846),
            (Some(Transition::Intersection), 7087.876080, 14222.446846),
            (Some(Transition::Arc), 7087.893057, 14222.446846),
            (Some(Transition::TrimToBend), 7086.876080, 14217.734457),
        ] {
            let a = flat_pattern(&part(t)).unwrap().area();
            assert!((a - area).abs() < 1e-6, "{t:?}: {a}");
            let v = crate::fold::exact_mass(&part(t)).unwrap().volume;
            assert!((v - volume).abs() < 1e-6, "{t:?}: {v}");
        }
    }

    #[test]
    fn holes_and_cutouts_on_a_flange_face_land_in_the_flat() {
        // Flange on edge 2 (along -x from (100, 60)): its bend ends 2.398
        // past the edge, its outer setback is 4. A hole 30 along, 15 up lands
        // at (70, 73.398); a slot [50..70] x [10..14] with round ends the same
        // way, still counter-clockwise.
        let part = |up: f64| -> PartSpec {
            let mut p = plate(vec![flange(2)]);
            p.face_cuts = vec![crate::spec::FaceCuts {
                on: "f2".into(),
                holes: vec![Hole {
                    center: [30.0, up],
                    diameter: 6.0,
                }],
                cutouts: vec![crate::spec::Cutout {
                    outline: vec![[50.0, 10.0], [70.0, 10.0], [70.0, 14.0], [50.0, 14.0]],
                    bulges: vec![0.0, 1.0, 0.0, 1.0],
                }],
            }];
            p
        };
        let l = layout(&part(15.0)).unwrap();
        let s = 60.0 + std::f64::consts::FRAC_PI_2 * 2.8 - 2.0 - 4.0;
        assert!(l
            .flat
            .holes
            .iter()
            .any(|h| dist(h.center, [70.0, s + 15.0]) < 1e-9));
        let slot = l.flat.cutouts.last().unwrap();
        assert!(slot.area() > 0.0);
        assert!((slot.area() - (20.0 * 4.0 + std::f64::consts::PI * 4.0)).abs() < 1e-9);
        assert!(has(&slot.points, [50.0, s + 10.0]) && has(&slot.points, [30.0, s + 14.0]));
        // 3 up, the hole reaches into the bend (the face starts 4 up).
        assert_eq!(
            layout(&part(6.5)).unwrap_err(),
            SheetError::Unsupported("face cuts on f2: reaching into a bend".into())
        );
        layout(&part(7.5)).unwrap();
    }

    #[test]
    fn a_punch_on_a_flange_stamps_as_its_face_cuts() {
        // An L-shaped cut-out and a hole, turned 90 degrees in the face's
        // [along, up], at two points: the same as face cuts stamped by hand.
        let punch = crate::spec::Punch {
            id: "L".into(),
            on: Some("f2".into()),
            at: vec![[30.0, 15.0], [60.0, 20.0]],
            angle: 90.0,
            holes: vec![Hole {
                center: [4.0, 4.0],
                diameter: 2.0,
            }],
            cutouts: vec![crate::spec::Cutout {
                outline: vec![
                    [-3.0, -3.0],
                    [5.0, -3.0],
                    [5.0, 0.0],
                    [0.0, 0.0],
                    [0.0, 3.0],
                    [-3.0, 3.0],
                ],
                bulges: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            }],
        };
        let mut punched = plate(vec![flange(2)]);
        punched.base.holes.clear();
        punched.base.punches = vec![punch.clone()];
        let mut cut = punched.clone();
        cut.base.punches.clear();
        for &at in &punch.at {
            let (holes, cutouts) = punch.stamped(at);
            cut.face_cuts.push(crate::spec::FaceCuts {
                on: "f2".into(),
                holes,
                cutouts,
            });
        }
        let (a, b) = (layout(&punched).unwrap().flat, layout(&cut).unwrap().flat);
        assert_eq!(a.holes.len(), 2);
        for (h, k) in a.holes.iter().zip(&b.holes) {
            assert!(dist(h.center, k.center) < 1e-9 && h.diameter == k.diameter);
        }
        assert_eq!(a.cutouts.len(), 2);
        for (p, q) in a.cutouts.iter().zip(&b.cutouts) {
            assert!(p.area() > 0.0 && (p.area() - q.area()).abs() < 1e-9);
            assert!(p.points.iter().all(|&x| has(&q.points, x)));
        }
        // The stamps: at their flat centres, and the flange's mirrored frame
        // (along -x, up +y) turns them the other way from its along.
        assert_eq!(a.punches.len(), 2);
        let s = 60.0 + std::f64::consts::FRAC_PI_2 * 2.8 - 2.0 - 4.0;
        assert!(dist(a.punches[0].center, [70.0, s + 15.0]) < 1e-9);
        assert!((a.punches[0].angle - 90.0).abs() < 1e-9);
        // Past the face into the bend: refused.
        punched.base.punches[0].at[0] = [30.0, 6.0];
        assert_eq!(
            layout(&punched).unwrap_err(),
            SheetError::Unsupported("punch L on f2: reaching into a bend".into())
        );
    }

    #[test]
    fn holes_must_cross_or_clear_the_outline_and_clear_each_other() {
        let invalid = |holes: &[([f64; 2], f64)]| {
            matches!(layout(&with_holes(holes)), Err(SheetError::Invalid(_)))
        };
        assert!(!invalid(&[([50.0, 3.0], 8.0)])); // through the outline: a round notch
        assert!(invalid(&[([50.0, 4.0], 8.0)])); // touching it
        assert!(invalid(&[([50.0, 120.0], 8.0)])); // off the part
        assert!(invalid(&[([50.0, 30.0], 8.0), ([55.0, 30.0], 8.0)])); // overlapping
        assert!(!invalid(&[([50.0, 30.0], 8.0), ([60.0, 30.0], 8.0)]));
    }

    #[test]
    fn unsupported_layouts_are_reported() {
        let unsupported =
            |part: &PartSpec| matches!(flat_pattern(part), Err(SheetError::Unsupported(_)));
        // Flanges alike meeting at a corner: a lip on one of them, or lips
        // alike on both (mitered), not different lips on both.
        let lip = |id: &str, on: &str, length: f64| Flange {
            id: id.into(),
            on: EdgeRef::Tip(on.into()),
            length,
            ..flange(2)
        };
        assert!(!unsupported(&plate(vec![
            flange(1),
            flange(2),
            lip("l2", "f2", 15.0)
        ])));
        assert!(!unsupported(&plate(vec![
            flange(1),
            flange(2),
            lip("l1", "f1", 15.0),
            lip("l2", "f2", 15.0)
        ])));
        assert!(unsupported(&plate(vec![
            flange(1),
            flange(2),
            lip("l1", "f1", 15.0),
            lip("l2", "f2", 10.0)
        ])));
        // A remnant shorter than 2T next to a relief runs the notch on to
        // the corner, also next to a flange on the other edge; not where
        // that flange ends past the notch's bottom short of the corner.
        let tight = Flange {
            width: Some(Width {
                offset_start: 5.0,
                offset_end: 20.0,
            }),
            ..flange(2)
        };
        assert!(!unsupported(&plate(vec![tight.clone()])));
        assert!(!unsupported(&plate(vec![flange(1), tight.clone()])));
        let past = Flange {
            width: Some(Width {
                offset_start: 30.0,
                offset_end: 1.0,
            }),
            ..flange(1)
        };
        assert!(unsupported(&plate(vec![past, tight])));
        // A relief notch reaching into the neighbouring flange's bend (that
        // bend starts 9.7 from the corner: inside, 135 degrees).
        let leaning = Flange {
            angle: 135.0,
            position: BendPosition::Inside,
            ..flange(2)
        };
        let near_corner = |start: f64, end: f64, edge: usize| Flange {
            width: Some(Width {
                offset_start: start,
                offset_end: end,
            }),
            ..flange(edge)
        };
        assert!(unsupported(&plate(vec![
            near_corner(10.0, 7.0, 1),
            leaning.clone()
        ])));
        assert!(!unsupported(&plate(vec![
            near_corner(10.0, 15.0, 1),
            leaning
        ])));
        // A notch running on to the corner across a flange on the other edge.
        assert!(unsupported(&plate(vec![
            near_corner(10.0, 1.0, 0),
            near_corner(5.0, 10.0, 1)
        ])));
    }

    fn parsed(json: &str) -> PartSpec {
        serde_json::from_str(json).unwrap()
    }

    #[test]
    fn bends_are_numbered_from_the_outside_in_unless_listed() {
        // A hem on a flange's tip, and a fold across the base.
        let json = |order: &str| {
            format!(
                r#"{{"id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]]}},
                "flanges": [{{"id": "f1", "on": {{"base": 2}}, "angle": 90, "length": 30,
                    "datum": "outer", "position": "outside", "side": "up"}}],
                "hems": [{{"id": "h", "on": {{"tip": "f1"}},
                    "kind": {{"single": {{"gap": 1, "length": 8}}}}, "side": "up"}}],
                "folds": [{{"id": "x", "line": [[0, 20], [100, 20]], "angle": 45,
                    "location": "centerline", "side": "down"}}],
                "bend_order": {order}}}"#
            )
        };
        let orders = |order: &str| {
            let f = flat_pattern(&parsed(&json(order))).unwrap();
            let mut o: Vec<(String, usize)> = f
                .bends
                .iter()
                .map(|b| (b.flange.clone(), b.order))
                .collect();
            o.sort();
            o
        };
        let o = |v: &[(&str, usize)]| {
            v.iter()
                .map(|(n, k)| (n.to_string(), *k))
                .collect::<Vec<_>>()
        };
        // The hem (two folds out) first, then the flange and the fold in
        // flat order.
        let default = orders("[]");
        assert_eq!(default.iter().find(|b| b.0 == "h").unwrap().1, 1);
        assert_eq!(default.iter().map(|b| b.1).max(), Some(3));
        assert_eq!(
            orders(r#"["x", "f1"]"#),
            o(&[("f1", 2), ("h", 3), ("x", 1)])
        );
        assert_eq!(orders(r#"["f1"]"#), o(&[("f1", 1), ("h", 2), ("x", 3)]));
        for bad in [r#"["nope"]"#, r#"["h", "h"]"#] {
            assert!(matches!(
                flat_pattern(&parsed(&json(bad))),
                Err(SheetError::Invalid(_))
            ));
        }
        let f = flat_pattern(&parsed(&json("[]"))).unwrap();
        let hem = f.bends.iter().find(|b| b.flange == "h").unwrap();
        assert_eq!(hem.note(), "1 UP 180° R0.5");
        let fold = f.bends.iter().find(|b| b.flange == "x").unwrap();
        assert!(fold.note().ends_with("DOWN 45° R2"));
    }

    #[test]
    fn the_pieces_of_one_fold_line_share_their_number() {
        // An E: the line crosses its three prongs.
        let f = flat_pattern(&parsed(
            r#"{"id": "e", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
            "base": {"outline": [[0, 0], [100, 0], [100, 60], [80, 60], [80, 20], [60, 20],
                [60, 60], [40, 60], [40, 20], [20, 20], [20, 60], [0, 60]]},
            "folds": [{"id": "f", "line": [[100, 40], [0, 40]], "angle": 90,
                "location": "centerline", "side": "up"}]}"#,
        ))
        .unwrap();
        assert_eq!(f.bends.len(), 3);
        assert!(f.bends.iter().all(|b| b.order == 1));
    }

    #[test]
    fn a_note_reads_along_its_bend_left_to_right() {
        let mut b = FlatBend {
            flange: "f".into(),
            order: 2,
            line: [[100.0, 10.0], [0.0, 10.0]],
            angle: 90.0,
            inner_radius: 2.0,
            k_factor: 0.4,
            up: true,
        };
        assert_eq!(b.note_place(1.0), ([50.0, 11.0], 0.0));
        b.line = [[5.0, 0.0], [5.0, 40.0]];
        assert_eq!(b.note_place(1.0), ([4.0, 20.0], 90.0));
    }

    #[test]
    fn a_punch_stamps_its_shape_at_every_point() {
        // A 4 mm hole and a 10 x 4 obround beside it, turned a quarter.
        let part = |at: &str| {
            parsed(&format!(
                r#"{{"id": "p", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                    "punches": [{{"id": "P1", "at": {at}, "angle": 90,
                        "holes": [{{"center": [0, 0], "diameter": 4}}],
                        "cutouts": [{{"outline": [[-3, 6], [3, 6], [3, 10], [-3, 10]],
                            "bulges": [0, 1, 0, 1]}}]}}]}}}}"#
            ))
        };
        let f = flat_pattern(&part("[[20, 30], [70, 20]]")).unwrap();
        assert_eq!(
            f.punches
                .iter()
                .map(|p| (p.id.as_str(), p.center))
                .collect::<Vec<_>>(),
            [("P1", [20.0, 30.0]), ("P1", [70.0, 20.0])]
        );
        assert_eq!(
            f.holes.iter().map(|h| h.center).collect::<Vec<_>>(),
            [[20.0, 30.0], [70.0, 20.0]]
        );
        // Turned a quarter: the obround's straight side from (3 - 6, ...)
        // lies left of the centre.
        let first = &f.cutouts[0].points;
        assert!(
            (first[0][0] - 14.0).abs() < 1e-12 && (first[0][1] - 27.0).abs() < 1e-12,
            "{first:?}"
        );
        let obround = 6.0 * 4.0 + std::f64::consts::PI * 4.0;
        let hole = std::f64::consts::PI * 4.0;
        assert!((f.area() - (6000.0 - 2.0 * (obround + hole))).abs() < 1e-9);
        // Mass agrees with the cut flat.
        let m = crate::fold::exact_mass(&part("[[20, 30], [70, 20]]")).unwrap();
        assert!((m.volume - 2.0 * f.area()).abs() < 1e-9);
        assert!(matches!(
            flat_pattern(&part("[]")),
            Err(SheetError::Invalid(_))
        ));
    }

    #[test]
    fn cut_length_of_a_plate_is_its_outline_and_hole() {
        let f = flat_pattern(&plate(Vec::new())).unwrap();
        let expected = 2.0 * (100.0 + 60.0) + std::f64::consts::PI * 8.0;
        assert!(
            (f.cut_length() - expected).abs() < 1e-9,
            "{}",
            f.cut_length()
        );
        assert_eq!(f.pierces(), 2);
    }

    #[test]
    fn cut_length_and_pierces_count_every_punch_stamp() {
        // Two stamps of a 20 x 6 slot with round ends: 2 * (2 * 14 + 6 pi) each.
        let mut spec = plate(Vec::new());
        spec.base.punches = vec![crate::spec::Punch {
            id: "slot".into(),
            on: None,
            at: vec![[50.0, 20.0], [50.0, 40.0]],
            angle: 0.0,
            holes: Vec::new(),
            cutouts: vec![crate::spec::Cutout {
                outline: vec![[-7.0, -3.0], [7.0, -3.0], [7.0, 3.0], [-7.0, 3.0]],
                bulges: vec![0.0, 1.0, 0.0, 1.0],
            }],
        }];
        let f = flat_pattern(&spec).unwrap();
        let slot = 2.0 * 14.0 + std::f64::consts::PI * 6.0;
        let expected = 2.0 * (100.0 + 60.0) + std::f64::consts::PI * 8.0 + 2.0 * slot;
        assert!(
            (f.cut_length() - expected).abs() < 1e-9,
            "{}",
            f.cut_length()
        );
        assert_eq!(f.pierces(), 4);
    }
}
