//! Faces on planes of their own, joined to the part by bends. The part
//! becomes its base plate with every face unrolled onto it, and a fold line
//! per bend. A bend runs where its two edges overlap; each face is extended
//! (by a strip as wide as the bend) or trimmed, square to the bend, to where
//! the bend meets it. Trimmed beside the bend, a face gets the part's bend
//! reliefs where it carries on past the bend, and that part stays.

use crate::cut::{cut, Clip};
use crate::flat::{Arc, Cell, Curve, Ring, Zone, P2};
use crate::spec::{
    Bend, BendLocation, Cutout, DoubleBend, EdgeRef, Face, FaceEdge, FixedEdge, FoldLine, Hole,
    PartSpec, Relief, ReliefShape, Side,
};
use crate::SheetError;

type V3 = [f64; 3];
const EPS: f64 = 1e-6;

/// The spec with its faces unrolled onto the base and its bends as fold
/// lines (before the part's own folds), every face edge as a base edge.
pub(crate) fn without_faces(spec: &PartSpec) -> Result<PartSpec, SheetError> {
    if spec.base.profile.is_some() || spec.base.cone.is_some() {
        return Err(SheetError::Unsupported(
            "faces and bends on a profile or cone part".into(),
        ));
    }
    let base_ring = Ring::from_cutout(&Cutout {
        outline: spec.base.outline.clone(),
        bulges: spec.base.bulges.clone(),
    });
    let mut placed = vec![Placed {
        id: "base".into(),
        o: [0.0; 3],
        x: [1.0, 0.0, 0.0],
        y: [0.0, 1.0, 0.0],
        n: [0.0, 0.0, 1.0],
        material: 1.0,
        top: [0.0, 0.0, 1.0],
        original: base_ring,
        map: Affine::IDENTITY,
        edits: Vec::new(),
        holes: Vec::new(),
        cutouts: Vec::new(),
    }];
    for (k, f) in spec.faces.iter().enumerate() {
        if f.id == "base" || spec.faces[..k].iter().any(|g| g.id == f.id) {
            return Err(SheetError::Invalid(format!(
                "face {}: ids must be unique and not \"base\"",
                f.id
            )));
        }
    }

    let mut folds = Vec::new();
    let mut joins = Vec::new();
    for bend in &spec.bends {
        let invalid = |what: String| SheetError::Invalid(format!("bend {}: {what}", bend.id));
        let p = placed
            .iter()
            .position(|f| f.id == bend.from.face)
            .ok_or_else(|| {
                invalid(format!(
                    "{} is neither the base nor a face joined by an earlier bend",
                    bend.from.face
                ))
            })?;
        if placed.iter().any(|f| f.id == bend.to.face) {
            return Err(SheetError::Unsupported(format!(
                "bend {}: face {} is joined already (a closed ring of faces)",
                bend.id, bend.to.face
            )));
        }
        let face = spec
            .faces
            .iter()
            .find(|f| f.id == bend.to.face)
            .ok_or_else(|| invalid(format!("no face {}", bend.to.face)))?;
        let mut c = Placed::face(face)?;
        let plan = plan(spec, bend, &placed[p], &c)?;

        // The child lands in the flat beyond the bends, along the parent's
        // edge as it runs in 3D.
        let par = &placed[p];
        let (v0, a2, _) = par.edge(bend.from.edge, &bend.id)?;
        let w2 = [a2[1], -a2[0]];
        let (a, d2) = (par.dir(a2), plan.d2);
        let (a0, b0) = (par.at(v0), plan.b0);
        let local = |t: f64, s: f64| add(v0, add(scale(a2, t), scale(w2, s)));
        let to_flat = |q2: P2| {
            let q = c.at(q2);
            let (t, h) = (dot3(sub3(q, a0), a), dot3(sub3(q, b0), d2));
            par.map.apply(local(t, plan.offset + h))
        };
        c.map = Affine::from_fn(to_flat);
        c.top = plan.top;

        for (k, s) in plan.segments.iter().enumerate() {
            let (mut l0, mut l1) = (
                par.map.apply(local(plan.lo, s.start)),
                par.map.apply(local(plan.hi, s.start)),
            );
            // The turning side right of the line.
            let out = par.map.linear(w2);
            let d = sub(l1, l0);
            if dot([d[1], -d[0]], out) < 0.0 {
                std::mem::swap(&mut l0, &mut l1);
            }
            folds.push(FoldLine {
                id: if k == 0 {
                    bend.id.clone()
                } else {
                    format!("{}.{}", bend.id, k + 1)
                },
                line: [l0, l1],
                angle: s.angle,
                location: BendLocation::Start,
                side: if s.up { Side::Up } else { Side::Down },
                bend_radius: Some(s.radius),
                k_factor: bend.k_factor,
                bounded: true,
            });
        }
        let end = plan.offset - plan.e_c;
        joins.push(Join {
            child: placed.len(),
            x: [
                par.map.apply(local(plan.lo, plan.e_p)),
                par.map.apply(local(plan.hi, plan.e_p)),
            ],
            y: [
                par.map.apply(local(plan.lo, end)),
                par.map.apply(local(plan.hi, end)),
            ],
        });
        if let Some(other) = placed[p]
            .edits
            .iter()
            .find(|e| e.edge == bend.from.edge && e.t0 < plan.hi - EPS && plan.lo < e.t1 - EPS)
        {
            return Err(SheetError::Unsupported(format!(
                "bend {}: bend {} runs on the same edge there",
                bend.id, other.bend
            )));
        }
        placed[p].edits.push(Edit {
            bend: bend.id.clone(),
            edge: bend.from.edge,
            t0: plan.lo,
            t1: plan.hi,
            e: plan.e_p,
        });
        c.edits.push(Edit {
            bend: bend.id.clone(),
            edge: bend.to.edge,
            t0: plan.c_lo,
            t1: plan.c_hi,
            e: plan.e_c,
        });
        placed.push(c);
    }
    if let Some(f) = spec
        .faces
        .iter()
        .find(|f| placed.iter().all(|p| p.id != f.id))
    {
        return Err(SheetError::Invalid(format!(
            "face {}: no bend joins it to the part",
            f.id
        )));
    }

    // Every plate as it lies in the flat, then all of them as one outline.
    let relief = spec.relief();
    let mut flat = Vec::new();
    for f in &placed {
        flat.push(f.map.ring(&f.edited(&relief)?));
    }
    let mut outline = flat[0].clone();
    for j in &joins {
        outline = splice(&outline, &flat[j.child], j)?;
    }
    // A strip's side runs on from its plates' sides: one edge.
    let outline = straight(outline);
    if let Some((a, b)) = overlap(&outline) {
        return Err(SheetError::Invalid(format!(
            "the faces overlap in the flat pattern, near [{:.3}, {:.3}] and [{:.3}, {:.3}]",
            a[0], a[1], b[0], b[1]
        )));
    }

    let mut part = spec.clone();
    let (points, bulges) = to_bulges(&outline);
    part.base.outline = points;
    part.base.bulges = bulges;
    for f in &placed[1..] {
        for h in &f.holes {
            part.base.holes.push(Hole {
                center: f.map.apply(h.center),
                diameter: h.diameter,
            });
        }
        for c in &f.cutouts {
            let (outline, bulges) = to_bulges(&f.map.ring(&Ring::from_cutout(c)));
            part.base.cutouts.push(Cutout { outline, bulges });
        }
    }
    part.faces = Vec::new();
    part.bends = Vec::new();
    folds.append(&mut part.folds);
    part.folds = folds;

    // Edges by face: where they are now in the outline.
    let edge_index = |e: &FaceEdge, id: &str| -> Result<usize, SheetError> {
        let f = placed
            .iter()
            .find(|f| f.id == e.face)
            .ok_or_else(|| SheetError::Invalid(format!("{id}: no face {} with a bend", e.face)))?;
        let n = f.original.points.len();
        if e.edge >= n {
            return Err(SheetError::Invalid(format!(
                "{id}: no edge {} on {}",
                e.edge, e.face
            )));
        }
        let (p, q) = (
            f.map.apply(f.original.points[e.edge]),
            f.map.apply(f.original.points[(e.edge + 1) % n]),
        );
        let m = outline.points.len();
        (0..m)
            .find(|&k| {
                let (a, b) = (outline.points[k], outline.points[(k + 1) % m]);
                (near(a, p) && near(b, q)) || (near(a, q) && near(b, p))
            })
            .ok_or_else(|| {
                SheetError::Unsupported(format!(
                    "{id}: edge {} of {} carries a bend or relief",
                    e.edge, e.face
                ))
            })
    };
    let remap = |on: &mut EdgeRef, id: &str| -> Result<(), SheetError> {
        let face = match on {
            EdgeRef::Base(i) => FaceEdge {
                face: "base".into(),
                edge: *i,
            },
            EdgeRef::Face(e) => e.clone(),
            EdgeRef::Tip(_) => return Ok(()),
        };
        *on = EdgeRef::Base(edge_index(&face, id)?);
        Ok(())
    };
    for f in &mut part.flanges {
        remap(&mut f.on, &f.id)?;
    }
    for h in &mut part.hems {
        remap(&mut h.on, &h.id)?;
    }
    for c in &mut part.contours {
        remap(&mut c.on, &c.id)?;
        for e in &mut c.along {
            let mut on = EdgeRef::Base(*e);
            remap(&mut on, &c.id)?;
            if let EdgeRef::Base(i) = on {
                *e = i;
            }
        }
    }
    Ok(part)
}

impl PartSpec {
    /// Whether a flange, hem or contour sits on a face's edge.
    pub(crate) fn has_face_edges(&self) -> bool {
        let on = |e: &EdgeRef| matches!(e, EdgeRef::Face(_));
        self.flanges.iter().any(|f| on(&f.on))
            || self.hems.iter().any(|h| on(&h.on))
            || self.contours.iter().any(|c| on(&c.on))
    }
}

/// A plate: the base, or a face once a bend has joined it.
struct Placed {
    id: String,
    o: V3,
    x: V3,
    y: V3,
    n: V3,
    /// The sheet on the normal's side (1) or the other (-1).
    material: f64,
    /// Where the flat pattern's +z points on this plate.
    top: V3,
    original: Ring,
    /// Its own x, y to the flat.
    map: Affine,
    edits: Vec<Edit>,
    holes: Vec<Hole>,
    cutouts: Vec<Cutout>,
}

/// A bend on a plate's edge: its span along the edge (from the edge's
/// start) and how far the plate moves out there (negative: trimmed back).
struct Edit {
    bend: String,
    edge: usize,
    t0: f64,
    t1: f64,
    e: f64,
}

/// Where a child joins: the parent's edge at the bend's start and the
/// child's at its end, at both ends of the bend, in the flat.
struct Join {
    child: usize,
    x: [P2; 2],
    y: [P2; 2],
}

impl Placed {
    fn face(f: &Face) -> Result<Placed, SheetError> {
        let invalid = |what: &str| SheetError::Invalid(format!("face {}: {what}", f.id));
        let (x, y) = (f.x_axis, f.y_axis);
        if (norm3(x) - 1.0).abs() > 1e-6 || (norm3(y) - 1.0).abs() > 1e-6 || dot3(x, y).abs() > 1e-6
        {
            return Err(invalid("x and y axes must be unit length and square"));
        }
        let ring = Ring::from_cutout(&Cutout {
            outline: f.outline.clone(),
            bulges: f.bulges.clone(),
        });
        if ring.points.len() < 3 || ring.area() <= 0.0 {
            return Err(invalid(
                "outline must be a counter-clockwise polygon, seen from the normal",
            ));
        }
        Ok(Placed {
            id: f.id.clone(),
            o: f.origin,
            x,
            y,
            n: cross3(x, y),
            material: if f.material == Side::Up { 1.0 } else { -1.0 },
            top: [0.0; 3],
            original: ring,
            map: Affine::IDENTITY,
            edits: Vec::new(),
            holes: f.holes.clone(),
            cutouts: f.cutouts.clone(),
        })
    }

    fn at(&self, p: P2) -> V3 {
        add3(self.o, self.dir(p))
    }

    fn dir(&self, v: P2) -> V3 {
        add3(scale3(self.x, v[0]), scale3(self.y, v[1]))
    }

    /// Edge `i`: its start, unit direction and length (lines only).
    fn edge(&self, i: usize, bend: &str) -> Result<(P2, P2, f64), SheetError> {
        let n = self.original.points.len();
        if i >= n {
            return Err(SheetError::Invalid(format!(
                "bend {bend}: no edge {i} on {}",
                self.id
            )));
        }
        if self.original.arcs[i].is_some() {
            return Err(SheetError::Unsupported(format!(
                "bend {bend}: on an arc edge of {}",
                self.id
            )));
        }
        let (a, b) = (self.original.points[i], self.original.points[(i + 1) % n]);
        let len = dist(a, b);
        Ok((a, scale(sub(b, a), 1.0 / len), len))
    }

    /// The outline with its bends' extensions, trims and reliefs.
    fn edited(&self, relief: &Relief) -> Result<Ring, SheetError> {
        let mut ring = self.original.clone();
        for ed in self.edits.iter().filter(|e| e.e <= EPS) {
            let (v, a, len) = self.edge(ed.edge, &ed.bend)?;
            let what = format!("bend {}", ed.bend);
            let w = [a[1], -a[0]];
            let at = |t: f64, s: f64| add(v, add(scale(a, t), scale(w, s)));
            if ring.points.iter().all(|&p| dot(sub(p, v), w) >= ed.e - EPS) {
                return Err(SheetError::Unsupported(format!(
                    "{what}: trimmed to the bend, nothing is left of {}",
                    self.id
                )));
            }
            let big = 10.0 + ring.points.iter().map(|&p| dist(p, v)).fold(0.0, f64::max);
            for clip in clips(ed, len, big, relief, &what)? {
                let pts: Vec<(P2, Option<P2>)> = clip
                    .into_iter()
                    .map(|(p, via)| (at(p[0], p[1]), via.map(|q| at(q[0], q[1]))))
                    .collect();
                let cell = Cell {
                    name: self.id.clone(),
                    polygon: ring.points.clone(),
                    arcs: ring.arcs.clone(),
                    holes: Vec::new(),
                    inner: Vec::new(),
                    zone: Zone::Base,
                };
                let mut pieces = cut(cell, &Clip::Profile(clip_ring(&pts)), false, &what)?;
                if pieces.len() != 1 {
                    return Err(SheetError::Unsupported(format!(
                        "{what}: trimming {} to the bend cuts it apart",
                        self.id
                    )));
                }
                let piece = pieces.remove(0);
                ring = Ring {
                    points: piece.polygon,
                    arcs: piece.arcs,
                };
            }
        }
        // Each trimmed edge's ends, all trims made.
        for ed in self.edits.iter().filter(|e| e.e <= EPS) {
            let (v, a, _) = self.edge(ed.edge, &ed.bend)?;
            let w = [a[1], -a[0]];
            for t in [ed.t0, ed.t1] {
                if !insert_vertex(&mut ring, add(v, add(scale(a, t), scale(w, ed.e)))) {
                    return Err(SheetError::Unsupported(format!(
                        "bend {}: {} doesn't reach across the bend where it is trimmed",
                        ed.bend, self.id
                    )));
                }
            }
        }
        // Then the extensions: a strip on the edge where the bend runs, the
        // edge as the trims left it.
        for ed in self.edits.iter().filter(|e| e.e > EPS) {
            let (v, a, _) = self.edge(ed.edge, &ed.bend)?;
            let w = [a[1], -a[0]];
            let at = |t: f64, s: f64| add(v, add(scale(a, t), scale(w, s)));
            let n = ring.points.len();
            let (t, s) = (|p: P2| dot(sub(p, v), a), |p: P2| dot(sub(p, v), w));
            let k = (0..n).find(|&k| {
                let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
                ring.arcs[k].is_none()
                    && s(p).abs() < EPS
                    && s(q).abs() < EPS
                    && t(p) <= ed.t0 + EPS
                    && t(q) >= ed.t1 - EPS
            });
            let Some(k) = k else {
                return Err(SheetError::Unsupported(format!(
                    "bend {}: another bend's trim cuts into the edge of {} it runs on",
                    ed.bend, self.id
                )));
            };
            let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
            let mut new = Vec::new();
            if t(p) < ed.t0 - EPS {
                new.push(at(ed.t0, 0.0));
            }
            new.extend([at(ed.t0, ed.e), at(ed.t1, ed.e)]);
            if t(q) > ed.t1 + EPS {
                new.push(at(ed.t1, 0.0));
            }
            for (m, p) in new.into_iter().enumerate() {
                ring.points.insert(k + 1 + m, p);
                ring.arcs.insert(k + 1 + m, None);
            }
        }
        Ok(ring)
    }
}

/// The regions a trim at `ed` takes off, in the edge's (t, s): s outward.
/// Each a closed path of points, each with the point an arc from it passes
/// through (a round relief's end).
#[allow(clippy::type_complexity)]
fn clips(
    ed: &Edit,
    len: f64,
    big: f64,
    relief: &Relief,
    what: &str,
) -> Result<Vec<Vec<(P2, Option<P2>)>>, SheetError> {
    let (w, d, e) = (relief.width, relief.depth, ed.e);
    let band = e < -EPS;
    // Per side: none (the bend reaches the edge's end), a notch, or a notch
    // running on past a short remnant. `side` is -1 before t0, +1 after t1.
    let side = |t: f64, remnant: f64, dir: f64| -> Result<Vec<(P2, Option<P2>)>, SheetError> {
        if remnant < EPS {
            return Ok(if band {
                vec![([t + dir * big, big], None), ([t + dir * big, e], None)]
            } else {
                Vec::new()
            });
        }
        if relief.shape == ReliefShape::Tear {
            return Err(SheetError::Unsupported(format!(
                "{what}: tear reliefs beside a bend between faces"
            )));
        }
        let mut path;
        if remnant - w <= relief.minimum_remnant + 1e-6 {
            if relief.shape == ReliefShape::Round {
                return Err(SheetError::Unsupported(format!(
                    "{what}: a round relief next to a short remnant"
                )));
            }
            path = vec![
                ([t + dir * big, big], None),
                ([t + dir * big, e - d], None),
                ([t, e - d], None),
            ];
        } else if relief.shape == ReliefShape::Round {
            let mid = e - d + w / 2.0;
            path = vec![
                ([t + dir * w, big], None),
                ([t + dir * w, mid], Some([t + dir * w / 2.0, e - d])),
                ([t, mid], None),
            ];
            if !band {
                path.push(([t, big], None));
            }
            return Ok(path);
        } else {
            path = vec![
                ([t + dir * w, big], None),
                ([t + dir * w, e - d], None),
                ([t, e - d], None),
            ];
        }
        if !band {
            path.push(([t, big], None));
        }
        Ok(path)
    };
    let left = side(ed.t0, ed.t0, -1.0)?;
    let mut right = side(ed.t1, len - ed.t1, 1.0)?;
    if !band {
        return Ok([left, right]
            .into_iter()
            .filter(|c| !c.is_empty())
            .collect());
    }
    // One region: down the left side, along the trim line, up the right.
    // The right side's path runs down too: turn it round (an arc's
    // through-point moves to the segment before it).
    let mut path = left;
    if path
        .last()
        .is_some_and(|p| p.0[0] == ed.t0 && p.0[1] < e - EPS)
    {
        path.push(([ed.t0, e], None));
    }
    let mut up: Vec<(P2, Option<P2>)> = Vec::new();
    for k in (0..right.len()).rev() {
        let via = if k > 0 { right[k - 1].1 } else { None };
        up.push((right[k].0, via));
    }
    right.clear();
    if up
        .first()
        .is_some_and(|p| p.0[0] == ed.t1 && p.0[1] < e - EPS)
    {
        path.push(([ed.t1, e], None));
    }
    path.extend(up);
    path.dedup_by(|a, b| near(a.0, b.0));
    Ok(vec![path])
}

/// A clip region from its path: counter-clockwise, arcs through the given
/// points.
fn clip_ring(path: &[(P2, Option<P2>)]) -> Ring {
    let n = path.len();
    let points: Vec<P2> = path.iter().map(|p| p.0).collect();
    let arcs = (0..n)
        .map(|k| {
            path[k].1.map(|via| {
                let (p, q) = (points[k], points[(k + 1) % n]);
                Curve::Arc(Arc {
                    center: circle_center(p, via, q),
                    ccw: cross(sub(via, p), sub(q, via)) > 0.0,
                })
            })
        })
        .collect();
    let ring = Ring { points, arcs };
    if ring.area() < 0.0 {
        ring.reversed()
    } else {
        ring
    }
}

/// A bend's layout: the plates' moves, the bends, and where the child lands.
struct Plan {
    /// The span along the parent's edge, from its start.
    lo: f64,
    hi: f64,
    /// The span along the child's edge.
    c_lo: f64,
    c_hi: f64,
    /// How far the parent and the child move out (+) or back (-).
    e_p: f64,
    e_c: f64,
    segments: Vec<Segment>,
    /// The child's edge, out from the parent's edge in the flat.
    offset: f64,
    /// The child's edge start and its unit direction into the child (3D).
    b0: V3,
    d2: V3,
    top: V3,
}

struct Segment {
    angle: f64,
    radius: f64,
    up: bool,
    /// Out from the parent's edge in the flat.
    start: f64,
}

fn plan(spec: &PartSpec, bend: &Bend, p: &Placed, c: &Placed) -> Result<Plan, SheetError> {
    let invalid = |what: &str| SheetError::Invalid(format!("bend {}: {what}", bend.id));
    let unsupported = |what: &str| SheetError::Unsupported(format!("bend {}: {what}", bend.id));
    let t = spec.thickness;
    let r = bend.bend_radius.unwrap_or(spec.bend_radius);
    let allowance = |angle: f64, r: f64| {
        spec.bend_allowance(angle, r, bend.k_factor)
            .map(|(ba, _)| ba)
    };

    let (v0, a2, lp) = p.edge(bend.from.edge, &bend.id)?;
    let (u0, b2, lc) = c.edge(bend.to.edge, &bend.id)?;
    let (a0, a, w, np) = (p.at(v0), p.dir(a2), p.dir([a2[1], -a2[0]]), p.n);
    let (b0, b) = (c.at(u0), c.dir(b2));
    let d2 = scale3(c.dir([b2[1], -b2[0]]), -1.0);
    if norm3(cross3(a, b)) > 1e-6 {
        return Err(invalid("its edges aren't parallel"));
    }
    // The span: where both edges run, along the parent's.
    let along = dot3(a, b);
    let (tb0, tb1) = (dot3(sub3(b0, a0), a), dot3(sub3(b0, a0), a) + lc * along);
    let (lo, hi) = (tb0.min(tb1).max(0.0), tb0.max(tb1).min(lp));
    if hi - lo < EPS {
        return Err(invalid("its edges don't overlap along the bend"));
    }
    let on_child = |t: f64| (t - tb0) * along;
    let (c_lo, c_hi) = (
        on_child(lo).min(on_child(hi)),
        on_child(lo).max(on_child(hi)),
    );

    // The section square to the bend: xi out from the parent's edge, eta
    // along its normal; the parent's sheet at eta 0 .. material * T.
    let sec = |v: V3| [dot3(v, w), dot3(v, np)];
    let bs = sec(sub3(b0, a0));
    let ds = sec(d2);
    let (p_lo, p_hi) = (0.0f64.min(p.material * t), 0.0f64.max(p.material * t));
    let cn = sec(c.n);

    let mut segments = Vec::new();
    let (e_p, e_c, top);
    if norm3(cross3(np, c.n)) > 1e-7 {
        if bend.double.is_some() || bend.fixed.is_some() {
            return Err(invalid("`double` and `fixed` are for parallel faces"));
        }
        let (cx, cy) = (ds[0], ds[1]);
        let theta = cx.clamp(-1.0, 1.0).acos();
        let turn = cy.signum();
        // Each sheet's face on the inside of the bend, the bend's centre R in
        // from both.
        let p_in = if turn > 0.0 { p_hi } else { p_lo };
        let eta = p_in + turn * r;
        let m = unit2([cx * cx - 1.0, cx * cy]);
        let c_in = (dot(cn, m) * c.material * t).max(0.0);
        let xi = bs[0] + (c_in + r - (eta - bs[1]) * m[1]) / m[0];
        e_p = xi;
        e_c = -dot(sub([xi, eta], bs), ds);
        let up = dot3(scale3(np, turn), p.top) > 0.0;
        let m3 = add3(scale3(w, m[0]), scale3(np, m[1]));
        top = if up { m3 } else { scale3(m3, -1.0) };
        segments.push(Segment {
            angle: theta.to_degrees(),
            radius: r,
            up,
            start: e_p,
        });
    } else {
        // Parallel: the child's sheet above or below the parent's.
        let s = dot3(c.n, np).signum();
        let (c_a, c_b) = (bs[1], bs[1] + c.material * t * s);
        let (cs_lo, cs_hi) = (c_a.min(c_b), c_a.max(c_b));
        let dir = if cs_lo >= p_hi - EPS {
            1.0
        } else if cs_hi <= p_lo + EPS {
            -1.0
        } else {
            return Err(invalid("its faces overlap through their thickness"));
        };
        // Facing sheet faces, seen with the child above.
        let (p_in, c_in) = if dir > 0.0 {
            (p_hi, cs_lo)
        } else {
            (-p_lo, -cs_hi)
        };
        let side_by_side = ds[0] > 0.0;
        let kind = bend.double.unwrap_or(if side_by_side {
            DoubleBend::Deg45
        } else {
            DoubleBend::Deg90
        });
        let keep_from = bend.fixed.unwrap_or(FixedEdge::From) == FixedEdge::From;
        let up = dot3(scale3(np, dir), p.top) > 0.0;
        let bx = bs[0];
        let rho = r + t / 2.0;
        match (side_by_side, kind) {
            (true, DoubleBend::FixedEdges | DoubleBend::Deg45) => {
                // Two bends the opposite way, joined by a straight run
                // tangent to both (their centre lines at mid thickness).
                let (phi, run, x1, x2);
                if kind == DoubleBend::FixedEdges {
                    let delta = [bx, (c_in - r) - (p_in + r)];
                    let dd = norm2(delta);
                    if dd < 2.0 * rho {
                        return Err(invalid("its faces are too close for two bends"));
                    }
                    phi = delta[1].atan2(delta[0]) + (2.0 * rho / dd).asin();
                    run = 2.0 * ((dd / 2.0).powi(2) - rho * rho).sqrt();
                    (x1, x2) = (0.0, bx);
                } else {
                    phi = std::f64::consts::FRAC_PI_4;
                    let dh = (c_in - r) - (p_in + r);
                    let dx = (dh * phi.cos() + 2.0 * rho) / phi.sin();
                    run = dx * phi.cos() + dh * phi.sin();
                    (x1, x2) = if keep_from { (0.0, dx) } else { (bx - dx, bx) };
                }
                if !(phi > 0.0 && phi < std::f64::consts::PI) || run < -EPS {
                    return Err(invalid("its faces are too close for two bends"));
                }
                let deg = phi.to_degrees();
                let ba = allowance(deg, r)?;
                e_p = x1;
                e_c = bx - x2;
                segments.push(Segment {
                    angle: deg,
                    radius: r,
                    up,
                    start: e_p,
                });
                segments.push(Segment {
                    angle: deg,
                    radius: r,
                    up: !up,
                    start: e_p + ba + run.max(0.0),
                });
                top = p.top;
            }
            (false, DoubleBend::Deg90) => {
                let wall = c_in - p_in - 2.0 * r;
                if wall < -EPS {
                    return Err(invalid("its faces are too close for two bends"));
                }
                let x1 = if keep_from { 0.0 } else { bx };
                let ba = allowance(90.0, r)?;
                e_p = x1;
                e_c = x1 - bx;
                segments.push(Segment {
                    angle: 90.0,
                    radius: r,
                    up,
                    start: e_p,
                });
                segments.push(Segment {
                    angle: 90.0,
                    radius: r,
                    up,
                    start: e_p + ba + wall.max(0.0),
                });
                top = scale3(p.top, -1.0);
            }
            (false, DoubleBend::FullRadius) => {
                let rf = (c_in - p_in) / 2.0;
                if rf < EPS {
                    return Err(invalid("its faces are too close for a full radius"));
                }
                let x1 = if keep_from { 0.0 } else { bx };
                e_p = x1;
                e_c = x1 - bx;
                segments.push(Segment {
                    angle: 180.0,
                    radius: rf,
                    up,
                    start: e_p,
                });
                top = scale3(p.top, -1.0);
            }
            (true, _) => {
                return Err(unsupported(
                    "faces side by side take fixed edges or 45 degrees",
                ))
            }
            (false, _) => return Err(unsupported("facing faces take 90 degrees or a full radius")),
        }
    }
    let last = segments.last().unwrap();
    let offset = last.start + allowance(last.angle, last.radius)? + e_c;
    Ok(Plan {
        lo,
        hi,
        c_lo,
        c_hi,
        e_p,
        e_c,
        segments,
        offset,
        b0,
        d2,
        top,
    })
}

/// The outline with the child's spliced in at its join: from the parent's
/// edge along the bend's side to the child, round the child, and back.
fn splice(outline: &Ring, child: &Ring, j: &Join) -> Result<Ring, SheetError> {
    let mismatch = || SheetError::Kernel("a bend's plates don't meet in the flat".into());
    let n = outline.points.len();
    let k = (0..n)
        .find(|&k| {
            let (a, b) = (outline.points[k], outline.points[(k + 1) % n]);
            (near(a, j.x[0]) && near(b, j.x[1])) || (near(a, j.x[1]) && near(b, j.x[0]))
        })
        .ok_or_else(mismatch)?;
    let first = if near(outline.points[k], j.x[0]) {
        0
    } else {
        1
    };
    let (y_first, y_second) = (j.y[first], j.y[1 - first]);
    let m = child.points.len();
    let find = |y: P2| (0..m).find(|&i| near(child.points[i], y));
    let (i_f, i_s) = (
        find(y_first).ok_or_else(mismatch)?,
        find(y_second).ok_or_else(mismatch)?,
    );
    if (i_s + 1) % m != i_f {
        return Err(mismatch());
    }
    // The outline from the edge's second point round to its first.
    let mut points: Vec<P2> = (0..n).map(|i| outline.points[(k + 1 + i) % n]).collect();
    let mut arcs: Vec<Option<Curve>> = (0..n - 1).map(|i| outline.arcs[(k + 1 + i) % n]).collect();
    arcs.push(None);
    let mut i = i_f;
    loop {
        points.push(child.points[i]);
        if i == i_s {
            break;
        }
        arcs.push(child.arcs[i]);
        i = (i + 1) % m;
    }
    arcs.push(None);
    Ok(Ring { points, arcs })
}

/// The ring without corners between two lines running straight on.
fn straight(ring: Ring) -> Ring {
    let (mut points, mut arcs) = (ring.points, ring.arcs);
    let mut k = 0;
    while k < points.len() && points.len() > 3 {
        let n = points.len();
        let prev = (k + n - 1) % n;
        let (a, p, b) = (points[prev], points[k], points[(k + 1) % n]);
        let (d1, d2) = (sub(p, a), sub(b, p));
        let on = arcs[prev].is_none()
            && arcs[k].is_none()
            && cross(d1, d2).abs() <= 1e-9 * norm2(d1) * norm2(d2)
            && dot(d1, d2) > 0.0;
        if on {
            points.remove(k);
            arcs.remove(k);
        } else {
            k += 1;
        }
    }
    Ring { points, arcs }
}

/// Two segments of the outline that cross or touch, if any (arcs as chords
/// of eight pieces).
pub(crate) fn overlap(ring: &Ring) -> Option<(P2, P2)> {
    let n = ring.points.len();
    let mut segs: Vec<(usize, P2, P2)> = Vec::new();
    for k in 0..n {
        let (p, q) = (ring.points[k], ring.points[(k + 1) % n]);
        match ring.arcs[k].and_then(Curve::arc) {
            Some(a) => {
                let (c, sweep) = (a.center, a.sweep(p, q));
                let (r, a0) = (dist(p, c), (p[1] - c[1]).atan2(p[0] - c[0]));
                let mut prev = p;
                for i in 1..=8 {
                    let x = if i == 8 {
                        q
                    } else {
                        let phi = a0 + sweep * i as f64 / 8.0;
                        [c[0] + r * phi.cos(), c[1] + r * phi.sin()]
                    };
                    segs.push((k, prev, x));
                    prev = x;
                }
            }
            None => segs.push((k, p, q)),
        }
    }
    for i in 0..segs.len() {
        for j in i + 1..segs.len() {
            let (ki, a, b) = segs[i];
            let (kj, c, d) = segs[j];
            // Neighbours share a point.
            let adjacent = ki == kj || (ki + 1) % n == kj || (kj + 1) % n == ki;
            if adjacent {
                continue;
            }
            if segments_meet(a, b, c, d) {
                return Some((a, c));
            }
        }
    }
    None
}

fn segments_meet(a: P2, b: P2, c: P2, d: P2) -> bool {
    let (d1, d2) = (sub(b, a), sub(d, c));
    let den = cross(d1, d2);
    if den.abs() < 1e-12 {
        // Parallel: touching only if collinear and overlapping.
        if cross(d1, sub(c, a)).abs() > 1e-9 * norm2(d1).max(1.0) {
            return false;
        }
        let l = dot(d1, d1);
        let (s0, s1) = (dot(sub(c, a), d1) / l, dot(sub(d, a), d1) / l);
        return s0.max(s1) > 1e-9 && s0.min(s1) < 1.0 - 1e-9;
    }
    let s = cross(sub(c, a), d2) / den;
    let u = cross(sub(c, a), d1) / den;
    (-1e-9..=1.0 + 1e-9).contains(&s) && (-1e-9..=1.0 + 1e-9).contains(&u)
}

/// Splits a line segment of the ring at p; false when p isn't on its
/// boundary's lines.
fn insert_vertex(ring: &mut Ring, p: P2) -> bool {
    if ring.points.iter().any(|&q| near(q, p)) {
        return true;
    }
    let n = ring.points.len();
    for k in 0..n {
        let (a, b) = (ring.points[k], ring.points[(k + 1) % n]);
        if ring.arcs[k].is_none() && segment_distance(p, a, b) < EPS {
            ring.points.insert(k + 1, p);
            ring.arcs.insert(k + 1, None);
            return true;
        }
    }
    false
}

fn to_bulges(ring: &Ring) -> (Vec<[f64; 2]>, Vec<f64>) {
    let n = ring.points.len();
    let bulges: Vec<f64> = (0..n)
        .map(|k| {
            ring.arcs[k].and_then(Curve::arc).map_or(0.0, |a| {
                (a.sweep(ring.points[k], ring.points[(k + 1) % n]) / 4.0).tan()
            })
        })
        .collect();
    let any = bulges.iter().any(|b| *b != 0.0);
    (ring.points.clone(), if any { bulges } else { Vec::new() })
}

/// A 2D affine map: m p + b.
#[derive(Clone, Copy)]
struct Affine {
    m: [[f64; 2]; 2],
    b: P2,
}

impl Affine {
    const IDENTITY: Affine = Affine {
        m: [[1.0, 0.0], [0.0, 1.0]],
        b: [0.0, 0.0],
    };

    fn from_fn(f: impl Fn(P2) -> P2) -> Affine {
        let o = f([0.0, 0.0]);
        let (ex, ey) = (sub(f([1.0, 0.0]), o), sub(f([0.0, 1.0]), o));
        Affine {
            m: [[ex[0], ey[0]], [ex[1], ey[1]]],
            b: o,
        }
    }

    fn linear(&self, v: P2) -> P2 {
        [
            self.m[0][0] * v[0] + self.m[0][1] * v[1],
            self.m[1][0] * v[0] + self.m[1][1] * v[1],
        ]
    }

    fn apply(&self, p: P2) -> P2 {
        add(self.linear(p), self.b)
    }

    fn det(&self) -> f64 {
        self.m[0][0] * self.m[1][1] - self.m[0][1] * self.m[1][0]
    }

    /// A ring mapped, counter-clockwise again if the map mirrors it.
    fn ring(&self, r: &Ring) -> Ring {
        let mirror = self.det() < 0.0;
        let mapped = Ring {
            points: r.points.iter().map(|&p| self.apply(p)).collect(),
            arcs: r
                .arcs
                .iter()
                .map(|a| {
                    a.and_then(Curve::arc).map(|a| {
                        Curve::Arc(Arc {
                            center: self.apply(a.center),
                            ccw: a.ccw != mirror,
                        })
                    })
                })
                .collect(),
        };
        if mirror {
            mapped.reversed()
        } else {
            mapped
        }
    }
}

fn circle_center(a: P2, b: P2, c: P2) -> P2 {
    let d = 2.0 * (a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1]));
    let (a2, b2, c2) = (dot(a, a), dot(b, b), dot(c, c));
    [
        (a2 * (b[1] - c[1]) + b2 * (c[1] - a[1]) + c2 * (a[1] - b[1])) / d,
        (a2 * (c[0] - b[0]) + b2 * (a[0] - c[0]) + c2 * (b[0] - a[0])) / d,
    ]
}

fn segment_distance(x: P2, a: P2, b: P2) -> f64 {
    let d = sub(b, a);
    let s = (dot(sub(x, a), d) / dot(d, d)).clamp(0.0, 1.0);
    dist(x, add(a, scale(d, s)))
}

fn near(a: P2, b: P2) -> bool {
    dist(a, b) < 1e-6
}
fn add(a: P2, b: P2) -> P2 {
    [a[0] + b[0], a[1] + b[1]]
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
fn cross(a: P2, b: P2) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn norm2(a: P2) -> f64 {
    dot(a, a).sqrt()
}
fn unit2(a: P2) -> P2 {
    scale(a, 1.0 / norm2(a))
}
fn dist(a: P2, b: P2) -> f64 {
    norm2(sub(a, b))
}
fn add3(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn scale3(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot3(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross3(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm3(a: V3) -> f64 {
    dot3(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use crate::spec::PartSpec;

    /// A base 100 wide and `depth` deep (edge 2 at y = depth) and one face
    /// joined to that edge.
    fn part(depth: f64, face: &str, to_edge: usize, extra: &str) -> PartSpec {
        let json = format!(
            r#"{{ "id": "p", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0, 0], [100, 0], [100, {depth}], [0, {depth}]] }},
                 "faces": [{{ "id": "b", {face} }}],
                 "bends": [{{ "id": "bend", "from": {{ "face": "base", "edge": 2 }},
                              "to": {{ "face": "b", "edge": {to_edge} }} {extra} }}] }}"#
        );
        serde_json::from_str(&json).unwrap()
    }

    fn wall(y: f64, outline: &str, material: &str) -> String {
        format!(
            r#""origin": [0, {y}, 0], "x_axis": [1, 0, 0], "y_axis": [0, 0, 1],
               "outline": {outline}, "material": "{material}""#
        )
    }

    fn above(z: f64, outline: &str) -> String {
        format!(
            r#""origin": [0, 0, {z}], "x_axis": [1, 0, 0], "y_axis": [0, 1, 0],
               "outline": {outline}, "material": "up""#
        )
    }

    const WALL: &str = "[[0, 10], [100, 10], [100, 50], [0, 50]]";

    #[test]
    fn faces_joined_by_bends_measure_as_the_reference_parts() {
        // Expected: the reference CAD's bend feature between the same faces
        // (volume, area, flat area, centre of mass). The solid meshes close
        // to it.
        let s60 = r#""origin": [0, 70, 0], "x_axis": [1, 0, 0], "y_axis": [0, 0.5, 0.8660254037844386],
                     "outline": [[0, 10], [100, 10], [100, 50], [0, 50]], "material": "up""#;
        let s120 = s60.replace("[0, 0.5,", "[0, -0.5,");
        let fixed = |kind: &str, keep: &str| format!(r#", "double": "{kind}", "fixed": "{keep}""#);
        let z = above(20.0, "[[0, 70], [100, 70], [100, 120], [0, 120]]");
        let u = above(20.0, WALL);
        let u2 = above(20.0, "[[0, 10], [100, 10], [100, 70], [0, 70]]");
        let slant = PartSpec {
            base: crate::spec::Base {
                outline: vec![[0.0, 0.0], [100.0, 0.0], [90.0, 60.0], [10.0, 60.0]],
                ..Default::default()
            },
            ..part(
                60.0,
                &wall(70.0, "[[0, 10], [100, 10], [90, 50], [10, 50]]", "up"),
                0,
                "",
            )
        };
        let cases: Vec<(&str, PartSpec, f64, f64, f64, [f64; 3])> = vec![
            (
                "l90",
                part(60.0, &wall(70.0, WALL, "up"), 0, ""),
                23342.477796,
                24209.327352,
                11639.822972,
                [50.0, 48.6011, 11.2886],
            ),
            (
                "short",
                part(
                    60.0,
                    &wall(70.0, "[[20, 10], [80, 10], [80, 50], [20, 50]]", "up"),
                    0,
                    "",
                ),
                18805.486678,
                19672.336234,
                9383.893783,
                [50.0, 43.8533, 8.6625],
            ),
            (
                "shift",
                part(
                    60.0,
                    &wall(70.0, "[[50, 10], [150, 10], [150, 50], [50, 50]]", "up"),
                    0,
                    "",
                ),
                21671.238898,
                22738.088454,
                10819.911486,
                [70.3856, 47.2163, 11.8937],
            ),
            (
                "slant",
                slant,
                20673.982237,
                21549.066519,
                10311.858377,
                [50.0, 47.7303, 11.1575],
            ),
            (
                "trim",
                part(80.0, &wall(70.0, WALL, "up"), 0, ""),
                23342.477796,
                24209.327352,
                11639.822972,
                [50.0, 48.6011, 11.2886],
            ),
            (
                "neg",
                part(60.0, &wall(70.0, WALL, "down"), 0, ""),
                23742.477796,
                24617.327352,
                11839.822972,
                [50.0, 49.7655, 11.1152],
            ),
            (
                "a60",
                part(60.0, s60, 0, ""),
                23704.558100,
                24578.649262,
                11831.335099,
                [50.0, 54.2431, 9.9242],
            ),
            (
                "a120",
                part(60.0, &s120, 0, ""),
                22485.355769,
                23335.062885,
                11200.789983,
                [50.0, 42.3746, 9.9650],
            ),
            (
                "shorttrim",
                part(
                    80.0,
                    &wall(70.0, "[[20, 10], [80, 10], [80, 50], [20, 50]]", "up"),
                    0,
                    "",
                ),
                20285.486678,
                21272.336234,
                10123.893783,
                [50.0, 45.7461, 8.1034],
            ),
            (
                "tight",
                part(
                    60.0,
                    &wall(62.0, "[[0, 0.5], [100, 0.5], [100, 50], [0, 50]]", "up"),
                    0,
                    "",
                ),
                21742.477796,
                22577.327352,
                10839.822972,
                [50.0, 43.8832, 12.0457],
            ),
            (
                "zero",
                part(
                    60.0,
                    &wall(64.0, "[[20, 10], [80, 10], [80, 50], [20, 50]]", "up"),
                    0,
                    "",
                ),
                18077.486678,
                18928.336234,
                9019.893783,
                [50.0, 41.0640, 8.9710],
            ),
            (
                "small",
                part(
                    60.0,
                    &wall(64.5, "[[20, 10], [80, 10], [80, 50], [20, 50]]", "up"),
                    0,
                    "",
                ),
                18145.486678,
                18990.336234,
                9053.893783,
                [50.0, 41.3032, 8.9412],
            ),
            (
                "narrow",
                part(
                    80.0,
                    &wall(70.0, "[[3, 10], [80, 10], [80, 50], [3, 50]]", "up"),
                    0,
                    "",
                ),
                21503.707903,
                22430.557459,
                10727.663688,
                [47.8817, 46.9625, 9.5996],
            ),
            (
                "oneear",
                part(
                    80.0,
                    &wall(70.0, "[[0, 10], [80, 10], [80, 50], [0, 50]]", "up"),
                    0,
                    "",
                ),
                21813.982237,
                22740.831793,
                10881.858377,
                [47.2220, 47.2736, 9.8076],
            ),
            (
                "z-fixed",
                part(60.0, &z, 0, &fixed("fixed_edges", "from")),
                26793.034680,
                27728.895373,
                13344.246287,
                [50.0, 60.5212, 10.2535],
            ),
            (
                "z-45",
                part(60.0, &z, 0, ""),
                25605.219496,
                26517.323886,
                12771.193821,
                [50.0, 60.7048, 9.2437],
            ),
            (
                "z-45-to",
                part(60.0, &z, 0, &fixed("deg45", "to")),
                25605.219496,
                26517.323886,
                12771.193821,
                [50.0, 59.9221, 11.1941],
            ),
            (
                "u-90",
                part(60.0, &u, 2, ""),
                26684.955592,
                27618.654704,
                13279.645943,
                [50.0, 37.5954, 10.2505],
            ),
            (
                "u-90-to",
                part(60.0, &u, 2, &fixed("deg90", "to")),
                22684.955592,
                23538.654704,
                11279.645943,
                [50.0, 32.4612, 10.1184],
            ),
            (
                "u-full",
                part(60.0, &u, 2, &fixed("full_radius", "from")),
                28283.185307,
                29248.849013,
                14078.760801,
                [50.0, 39.8514, 10.2929],
            ),
            (
                "u-full-to",
                part(60.0, &u, 2, &fixed("full_radius", "to")),
                24283.185307,
                25168.849013,
                12078.760801,
                [50.0, 34.7686, 10.1764],
            ),
            (
                "u2-90",
                part(60.0, &u2, 2, ""),
                26684.955592,
                27618.654704,
                13279.645943,
                [50.0, 37.5954, 10.2505],
            ),
            (
                "u2-90-to",
                part(60.0, &u2, 2, &fixed("deg90", "to")),
                30684.955592,
                31698.654704,
                15279.645943,
                [50.0, 42.6945, 10.3482],
            ),
            (
                "u2-full-to",
                part(60.0, &u2, 2, &fixed("full_radius", "to")),
                32283.185307,
                33328.849013,
                16078.760801,
                [50.0, 44.9137, 10.3805],
            ),
        ];
        let mut bad = Vec::new();
        for (name, spec, volume, area, flat_area, centre) in cases {
            let result = (|| {
                let flat = crate::flat::flat_pattern(&spec)?;
                let mass = crate::fold::exact_mass(&spec)?;
                let mesh = crate::fold::mass(&crate::fold::fold(&spec)?, 1e-3)?;
                Ok::<_, crate::SheetError>((flat.area(), mass, mesh))
            })();
            match result {
                Ok((flat, mass, mesh)) => {
                    let c = mass.center_of_mass;
                    let off = (0..3).map(|i| (c[i] - centre[i]).abs()).fold(0.0, f64::max);
                    if (mass.volume - volume).abs() > 1e-5
                        || (mass.area - area).abs() > 1e-5
                        || (flat - flat_area).abs() > 1e-5
                        || off > 1e-4
                        || ((mesh.volume - volume) / volume).abs() > 5e-5
                    {
                        bad.push(format!(
                            "{name}: volume {:.6} ({volume}), area {:.6} ({area}), flat {flat:.6} ({flat_area}), centre {c:.4?} ({centre:?}), mesh {:.3}",
                            mass.volume, mass.area, mesh.volume
                        ));
                    }
                }
                Err(e) => bad.push(format!("{name}: {e}")),
            }
        }
        assert!(bad.is_empty(), "{}", bad.join("\n"));
    }

    #[test]
    fn a_face_carries_holes_cut_outs_and_flanges_like_a_flange_would() {
        // The wall as a face joined by a bend, and as a flange on a base
        // reaching to its inner face: a hole and a slot on it, a lip on its
        // top edge.
        // Its plane seen from either side (the flat mirrors it back; the
        // normal turns round, and the material with it).
        let cut = |x0: f64, x1: f64| {
            format!(
                r#""holes": [{{ "center": [50, 30], "diameter": 10 }}],
                   "cutouts": [{{ "outline": [[{x0}, 20], [{x1}, 20], [{x1}, 26], [{x0}, 26]], "bulges": [0, 1, 0, 1] }}]"#
            )
        };
        let lip = r#"{ "id": "lip", "on": ON, "angle": 90, "length": 20, "datum": "outer",
                       "position": "outside", "side": "up" }"#;
        let frames = [
            format!("{}, {}", wall(70.0, WALL, "up"), cut(10.0, 30.0)),
            format!(
                r#""origin": [100, 70, 0], "x_axis": [-1, 0, 0], "y_axis": [0, 0, 1], "outline": {WALL}, "material": "down", {}"#,
                cut(70.0, 90.0)
            ),
        ];
        let flanged: PartSpec = serde_json::from_str(&format!(
            r#"{{ "id": "f", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": {{ "outline": [[0, 0], [100, 0], [100, 68], [0, 68]] }},
                 "flanges": [{{ "id": "w", "on": {{ "base": 2 }}, "angle": 90, "length": 50,
                                "datum": "outer", "position": "outside", "side": "up" }}, {}],
                 "face_cuts": [{{ "on": "w", "holes": [{{ "center": [50, 30], "diameter": 10 }}],
                                  "cutouts": [{{ "outline": [[70, 20], [90, 20], [90, 26], [70, 26]], "bulges": [0, 1, 0, 1] }}] }}] }}"#,
            lip.replace("ON", r#"{ "tip": "w" }"#)
        ))
        .unwrap();
        let b = crate::fold::exact_mass(&flanged).unwrap();
        let other = crate::flat::flat_pattern(&flanged).unwrap();
        for face in frames {
            let mut faces = part(60.0, &face, 0, "");
            faces.flanges = vec![serde_json::from_str(
                &lip.replace("ON", r#"{ "face": { "face": "b", "edge": 2 } }"#),
            )
            .unwrap()];
            let a = crate::fold::exact_mass(&faces).unwrap();
            assert!(
                (a.volume - b.volume).abs() < 1e-8,
                "{} {}",
                a.volume,
                b.volume
            );
            assert!((a.area - b.area).abs() < 1e-8, "{} {}", a.area, b.area);
            for i in 0..3 {
                assert!((a.center_of_mass[i] - b.center_of_mass[i]).abs() < 1e-9);
            }
            let flat = crate::flat::flat_pattern(&faces).unwrap();
            assert!((flat.area() - other.area()).abs() < 1e-8);
            crate::fold::fold(&faces).unwrap();
        }
    }

    #[test]
    fn what_bends_between_faces_refuse() {
        let refused = |spec: PartSpec, what: &str| {
            let e = crate::flat::flat_pattern(&spec).unwrap_err().to_string();
            assert!(e.contains(what), "{e}");
        };
        let wall = wall(70.0, WALL, "up");
        refused(
            part(60.0, &wall, 0, r#", "double": "deg45""#),
            "for parallel faces",
        );
        refused(
            part(60.0, &above(20.0, WALL), 2, r#", "double": "deg45""#),
            "facing faces take 90 degrees or a full radius",
        );
        refused(
            part(
                60.0,
                &above(20.0, "[[0, 70], [100, 70], [100, 120], [0, 120]]"),
                0,
                r#", "double": "full_radius""#,
            ),
            "faces side by side take fixed edges or 45 degrees",
        );
        refused(
            part(60.0, &above(1.0, WALL), 2, ""),
            "overlap through their thickness",
        );
        refused(part(60.0, &wall, 1, ""), "aren't parallel");
        refused(
            part(
                60.0,
                &wall
                    .replace("[[0, 10], [100, 10]", "[[200, 10], [300, 10]")
                    .replace("[100, 50], [0, 50]", "[300, 50], [200, 50]"),
                0,
                "",
            ),
            "don't overlap",
        );
        // A face no bend joins, and one joined twice.
        let mut lone = part(60.0, &wall, 0, "");
        lone.faces
            .push(PartSpec { ..lone.clone() }.faces[0].clone());
        lone.faces[1].id = "c".into();
        refused(lone, "no bend joins it");
        let mut twice = part(60.0, &wall, 0, "");
        let again = twice.bends[0].clone();
        twice.bends.push(again);
        refused(twice, "joined already");
        // Bent back flat onto the base: they overlap in the flat.
        refused(
            part(
                60.0,
                r#""origin": [0, 70, 0], "x_axis": [1, 0, 0], "y_axis": [0, -0.9998476951563913, 0.01745240643728351],
                   "outline": [[0, 10], [100, 10], [100, 50], [0, 50]], "material": "up""#,
                0,
                "",
            ),
            "",
        );
    }
}
