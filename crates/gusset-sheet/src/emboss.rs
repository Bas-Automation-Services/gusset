//! Formed embosses (`base.embosses`): the sheet pressed out of the base plate
//! at its thickness. Seen in a section square to its edge, each face of the
//! sheet runs from the plate over a rounded corner, along the wall, over a
//! second rounded corner into the raised face. That section, (s, z) with s
//! out from the raised face's mold line, is the same all round: a round
//! emboss turns it around its axis (tori, a cone or a cylinder, a disk); a
//! rectangular one runs it along each side (cylinders and a plane) and turns
//! it around each corner's axis (tori and a cone). All exact. The plate keeps
//! a hole the shape of the footprint, where the sheet leaves it.
//!
//! Louvers (`base.louvers`) are worked out the same way: their section runs
//! from the footprint over a corner and round the hood down to the lance;
//! it runs along the hood's straight (cylinders) and turns a quarter turn
//! round each end's axis (tori, spheres). Two plane faces close it at the
//! lance.

use crate::flat::{Arc, Curve, Ring, P2};
use crate::fold::{Mass, Rigid, GAUSS, P3};
use crate::spec::{Emboss, EmbossShape, FlangedHole, Hole, HoleShape, Louver, Side};
use crate::SheetError;
use opengeometry::brep::{
    boundary, plane_boundary, projected_boundary, uv_line, Builder, CurveGeometry, Frame3,
    GeometryError, Orientation, Surface, SurfaceGeometry, Use,
};
use opengeometry::math::Interval;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_3, FRAC_PI_4, PI, TAU};

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
    /// A louver's hood, its straight `length` along x between its ends' axes
    /// at (-length / 2, 0) and (length / 2, 0), x along the lance and y the
    /// way it opens (the hood at y < 0); turned by `angle` (radians).
    Louver { length: f64, angle: f64 },
    /// A bead's centre line, `length` along x between its ends' centres,
    /// turned by `angle` (radians); its section reaches the line.
    Slot { length: f64, angle: f64 },
}

/// Where a louver's section meets its lance.
#[derive(Clone, Copy, Debug)]
struct Lance {
    /// The level of the face that leaves the plate.
    out: usize,
    /// Where that face's section crosses the plate's other face: the piece,
    /// how far along it, and s there.
    piece: usize,
    f: f64,
    cross: f64,
    /// The two plane faces' area.
    area: f64,
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
    /// z = T): its section from the footprint inward (corner, wall, corner;
    /// no wall where the corners meet), and the height of its raised face.
    chains: [Vec<Meridian>; 2],
    tops: [f64; 2],
    /// s where the sheet leaves the plate, and where the raised face starts.
    foot: f64,
    top: f64,
    /// A louver's lance.
    lance: Option<Lance>,
    /// A dome's or a bead's: the section ends in the crown at the centre
    /// (line), no raised plane face.
    crowned: bool,
    /// A flanged hole's: the section ends at the collar's end, square to its
    /// wall (no raised face); its pre-hole's radius in the flat.
    open: Option<f64>,
}

impl Formed {
    pub fn new(e: &Emboss, t: f64, bend_radius: f64) -> Result<Formed, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("emboss {}: {what}", e.id)));
        let (plan, rc) = match e.shape {
            EmbossShape::Dome { radius } => {
                return Formed::crowned(e, t, bend_radius, radius, None)
            }
            EmbossShape::Bead { length, radius } => {
                return Formed::crowned(e, t, bend_radius, radius, Some(length))
            }
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
        // The wall's length; none where the corners meet (as tools with
        // tangent corners have them).
        let wall = (b[1] - a[1]) * d.cos() - (b[0] - a[0]) * d.sin();
        if wall < -EPS {
            return invalid("too low for its radius and thickness");
        }
        let corner = |c: P2, rho: f64, a0: f64, a1: f64| Meridian::Arc { c, rho, a0, a1 };
        let mut up = [
            vec![
                corner(c1, r + t, 1.5 * PI, PI + d),
                Meridian::Line {
                    a: off(c1, -(r + t)),
                    b: off(c2, r),
                },
                corner(c2, r, d, FRAC_PI_2),
            ],
            vec![
                corner(c1, r, 1.5 * PI, PI + d),
                Meridian::Line { a, b },
                corner(c2, r + t, d, FRAC_PI_2),
            ],
        ];
        if wall <= EPS {
            for chain in &mut up {
                chain.remove(1);
            }
        }
        let (chains, tops) = match e.side {
            Side::Up => (up, [h, t + h]),
            // Mirrored, the face that left z = T leaves z = 0.
            Side::Down => {
                let flip = |c: &Vec<Meridian>| c.iter().map(|m| m.mirrored(t)).collect();
                ([flip(&up[1]), flip(&up[0])], [-h, t - h])
            }
        };
        Ok(Formed {
            id: e.id.clone(),
            center: e.at,
            plan,
            chains,
            tops,
            foot,
            top,
            lance: None,
            crowned: false,
            open: None,
        })
    }

    /// A dome (no length) or a bead: from the footprint over the corner
    /// straight into the crown, an arc of inside radius `crown` reaching the
    /// centre (line). Worked out rising, s from the centre (line).
    fn crowned(
        e: &Emboss,
        t: f64,
        bend_radius: f64,
        crown: f64,
        length: Option<f64>,
    ) -> Result<Formed, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("emboss {}: {what}", e.id)));
        let (r, h) = (e.radius.unwrap_or(bend_radius), e.height);
        if !(h > 0.0 && r > 0.0 && crown > 0.0 && length.is_none_or(|l| l > 0.0)) {
            return invalid("needs a positive height, radius, crown radius and length");
        }
        if e.draft != 0.0 {
            return invalid("a dome or a bead has no wall to draft");
        }
        // The corner around c1 over the plate's face at T (inside r on that
        // face), the crown around its centre on the axis at zc (inside on
        // the face at z = 0), outermost at T + h; the two touch.
        let (c1z, zc) = (t + r, h - crown);
        let d = crown + r + t;
        if !(zc < c1z - EPS) {
            return invalid("too high for its crown: it would turn back under itself");
        }
        let foot = (d * d - (c1z - zc).powi(2)).sqrt();
        let (c1, cc) = ([foot, c1z], [0.0, zc]);
        let bend = (zc - c1z).atan2(-foot) + TAU;
        let cap = (c1z - zc).atan2(foot);
        let arc = |c: P2, rho: f64, a0: f64, a1: f64| Meridian::Arc { c, rho, a0, a1 };
        let up = [
            vec![
                arc(c1, r + t, 1.5 * PI, bend),
                arc(cc, crown, cap, FRAC_PI_2),
            ],
            vec![
                arc(c1, r, 1.5 * PI, bend),
                arc(cc, crown + t, cap, FRAC_PI_2),
            ],
        ];
        let (chains, tops) = match e.side {
            Side::Up => (up, [h, t + h]),
            Side::Down => {
                let flip = |c: &Vec<Meridian>| c.iter().map(|m| m.mirrored(t)).collect();
                ([flip(&up[1]), flip(&up[0])], [-h, t - h])
            }
        };
        let plan = match length {
            None => Plan::Round(0.0),
            Some(length) => Plan::Slot {
                length,
                angle: e.angle.to_radians(),
            },
        };
        Ok(Formed {
            id: e.id.clone(),
            center: e.at,
            plan,
            chains,
            tops,
            foot,
            top: 0.0,
            lance: None,
            crowned: true,
            open: None,
        })
    }

    /// A louver: its section from the footprint to the lance, worked out
    /// pressed down from z = T and mirrored for one that rises.
    pub fn louver(l: &Louver, t: f64, bend_radius: f64) -> Result<Formed, SheetError> {
        let invalid = |what: &str| Err(SheetError::Invalid(format!("louver {}: {what}", l.id)));
        let (big, r, len) = (l.radius, l.bend_radius.unwrap_or(bend_radius), l.length);
        if !(len > 0.0 && r >= 0.0) {
            return invalid("needs a positive length and a bend radius of 0 or more");
        }
        if !(big > t + EPS) {
            return invalid("its radius must be more than the thickness, or it doesn't open");
        }
        // (s, z), s from the lance towards the hood: the corner where it
        // leaves the plate around fc, the hood around its axis on the lance
        // in the face at z = T; the two touch.
        let d = big + r + t;
        let a = (d * d - (r + t).powi(2)).sqrt();
        // Turned round an end's axis the corner's outside must stay clear of
        // it (a ring torus): R > (sqrt 2 - 1)(r + T).
        if a <= r + t + EPS {
            return invalid("its bend radius is too large for its radius: the corner closes on itself at the ends");
        }
        let (fc, hc) = ([a, -r], [0.0, t]);
        let (bend, hood) = ((t + r).atan2(-a), (-(t + r)).atan2(a));
        let arc = |c: P2, rho: f64, a0: f64, a1: f64| Meridian::Arc { c, rho, a0, a1 };
        let down = [
            vec![
                arc(fc, r, FRAC_PI_2, bend),
                arc(hc, big + t, hood, -FRAC_PI_2),
            ],
            vec![
                arc(fc, r + t, FRAC_PI_2, bend),
                arc(hc, big, hood, -FRAC_PI_2),
            ],
        ];
        // The face that leaves z = T crosses z = 0 on its corner, where the
        // sine comes to r / (r + T), or on the hood.
        let (piece, angle) = if r / (r + t) >= (t + r) / d {
            (0, PI - (r / (r + t)).asin())
        } else {
            (1, -(t / big).asin())
        };
        let Meridian::Arc { a0, a1, .. } = down[1][piece] else {
            unreachable!()
        };
        let f = (angle - a0) / (a1 - a0);
        let cross = down[1][piece].at(f).0[0];
        // The lance faces: the plate's side at the lance (len + 2a by T) and
        // the hood's section there (len by T, and the section at each end),
        // less twice where they overlap: under that face, above z = 0.
        let section: f64 = (0..2).map(|i| between(down[0][i], down[1][i])[0]).sum();
        let mut under = a * t;
        for (i, m) in down[1].iter().enumerate().take(piece + 1) {
            let m = match *m {
                Meridian::Arc { c, rho, a0, .. } if i == piece => arc(c, rho, a0, angle),
                m => m,
            };
            under += green(m)[0];
        }
        let area = (len + 2.0 * a) * t + len * t + 2.0 * section - 4.0 * under.abs();
        let (chains, out) = match l.side {
            Side::Down => (down, 1),
            // Mirrored, the face that left z = T leaves z = 0.
            Side::Up => {
                let flip = |c: &Vec<Meridian>| c.iter().map(|m| m.mirrored(t)).collect();
                ([flip(&down[1]), flip(&down[0])], 0)
            }
        };
        Ok(Formed {
            id: l.id.clone(),
            center: l.at,
            plan: Plan::Louver {
                length: len,
                angle: l.angle.to_radians() - FRAC_PI_2,
            },
            chains,
            tops: [0.0; 2],
            foot: a,
            top: 0.0,
            lance: Some(Lance {
                out,
                piece,
                f,
                cross,
                area,
            }),
            crowned: false,
            open: None,
        })
    }

    /// A flanged hole: from the footprint over the corner along the collar's
    /// wall to its end, cut square to the wall. Worked out rising, s from the
    /// axis: the inside wall (the face leaving z = 0) meets z = 0 at d / 2.
    pub fn flanged(
        h: &FlangedHole,
        t: f64,
        bend_radius: f64,
        k: f64,
    ) -> Result<Formed, SheetError> {
        let invalid = |what: &str| {
            Err(SheetError::Invalid(format!(
                "flanged hole {}: {what}",
                h.id
            )))
        };
        let (r, d) = (h.radius.unwrap_or(bend_radius), h.draft.to_radians());
        // A round hole's s is the distance from its axis; a rectangle's from
        // its inside wall's outline at the face (plan corners rc).
        let (plan, base, rc) = match h.shape {
            HoleShape::Round { diameter } => (Plan::Round(0.0), diameter / 2.0, f64::INFINITY),
            HoleShape::Rectangle {
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
                    angle: h.angle.to_radians(),
                };
                (plan, 0.0, corner_radius)
            }
        };
        if !(base + rc.min(1.0) > 0.0 && h.height > 0.0 && r >= 0.0) {
            return invalid("needs a positive size and height, and a radius of 0 or more");
        }
        if !(d >= 0.0 && d < FRAC_PI_2) {
            return invalid("the draft must lie between 0 and 90 degrees");
        }
        // The inside wall: (d / 2 - u sin, u cos) for u along it; its normal
        // out from the axis n. The corner's centre lies r + T off it, at
        // z = T + r; the wall's end at z = height, the end face T along n.
        let (sn, cs) = d.sin_cos();
        let n = [cs, sn];
        let wall = |u: f64| [base - u * sn, u * cs];
        let u0 = (t + r) * (1.0 - sn) / cs;
        let u1 = h.height / cs;
        if !(u1 > u0 + EPS) {
            return invalid("too low for its radius and thickness: no wall");
        }
        let end = wall(u1);
        // What s may come down to: a round hole's axis, a rectangle's corner
        // centres.
        let floor = if rc.is_finite() { -rc } else { 0.0 };
        if !(end[0] - floor > EPS) {
            return invalid("its wall closes the hole (or its corners) before it ends");
        }
        let at = wall(u0);
        let c1 = [at[0] + (r + t) * n[0], t + r];
        let off = |p: P2, by: f64| [p[0] + by * n[0], p[1] + by * n[1]];
        let corner = |rho: f64| Meridian::Arc {
            c: c1,
            rho,
            a0: 1.5 * PI,
            a1: PI + d,
        };
        let up = [
            vec![corner(r + t), Meridian::Line { a: at, b: end }],
            vec![
                corner(r),
                Meridian::Line {
                    a: off(at, t),
                    b: off(end, t),
                },
            ],
        ];
        // The pre-hole: the collar along its neutral fibre, unrolled from the
        // footprint inwards.
        let developed = (r + k * t) * (FRAC_PI_2 - d) + (u1 - u0);
        let pre = c1[0] - developed;
        // A rectangle's pre-hole may go in past its plan corners (sharp then).
        let half = match h.shape {
            HoleShape::Rectangle { width, length, .. } => width.min(length) / 2.0,
            HoleShape::Round { .. } => 0.0,
        };
        if !(pre - floor.min(-half) > EPS) {
            return invalid("too tall to develop from the plate: no pre-hole is left");
        }
        let (chains, tops) = match h.side {
            Side::Up => (up, [0.0; 2]),
            Side::Down => {
                let flip = |c: &Vec<Meridian>| c.iter().map(|m| m.mirrored(t)).collect();
                ([flip(&up[1]), flip(&up[0])], [0.0; 2])
            }
        };
        Ok(Formed {
            id: h.id.clone(),
            center: h.at,
            plan,
            chains,
            tops,
            foot: c1[0],
            top: 0.0,
            lance: None,
            crowned: false,
            open: Some(pre),
        })
    }

    /// A round flanged hole's pre-hole in the flat pattern.
    pub fn pre_hole(&self) -> Option<Hole> {
        match self.plan {
            Plan::Round(_) => self.open.map(|r| Hole {
                center: self.center,
                diameter: 2.0 * r,
            }),
            _ => None,
        }
    }

    /// A rectangular flanged hole's pre-hole in the flat pattern: the
    /// outline offset to it.
    pub fn pre_cutout(&self) -> Option<Ring> {
        match (self.plan, self.open) {
            (Plan::Rect { rc, .. }, Some(s)) if rc + s > EPS => Some(self.ring_at(s)),
            // In past the plan corners: a sharp rectangle.
            (Plan::Rect { hw, hl, .. }, Some(s)) => Some(Ring {
                points: [[1.0, -1.0], [1.0, 1.0], [-1.0, 1.0], [-1.0, -1.0]]
                    .map(|[x, y]| self.place([x * (hw + s), y * (hl + s)]))
                    .to_vec(),
                arcs: vec![None; 4],
            }),
            _ => None,
        }
    }

    /// What it is, for messages: "emboss e1", "louver l1".
    pub fn what(&self) -> String {
        match self.plan {
            Plan::Louver { .. } => format!("louver {}", self.id),
            _ if self.open.is_some() => format!("flanged hole {}", self.id),
            _ => format!("emboss {}", self.id),
        }
    }

    /// A round emboss's footprint, as the hole the plate keeps.
    pub fn hole(&self) -> Option<Hole> {
        match self.plan {
            Plan::Round(half) => Some(Hole {
                center: self.center,
                diameter: 2.0 * (half + self.foot),
            }),
            Plan::Rect { .. } | Plan::Louver { .. } | Plan::Slot { .. } => None,
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
            Plan::Rect { .. } => self.ring_at(self.foot),
            Plan::Louver { .. } | Plan::Slot { .. } => self.ring_other(),
        }
    }

    /// A rectangle's outline offset by s, counter-clockwise; its points the
    /// eight segments' starts.
    fn ring_at(&self, s: f64) -> Ring {
        Ring {
            points: (0..8).map(|j| self.flat_at(j, 0.0, s)).collect(),
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
        }
    }

    /// A louver's or a bead's footprint.
    fn ring_other(&self) -> Ring {
        match self.plan {
            // The far side, the far end, the lance (its end beyond the
            // hood, the part across from it, its near end), the near end.
            Plan::Louver { length, .. } => {
                let (h, a) = (length / 2.0, self.foot);
                let c = self.lance.expect("a louver's lance").cross;
                let arc = |x: f64| {
                    Some(Curve::Arc(Arc {
                        center: self.place([x, 0.0]),
                        ccw: true,
                    }))
                };
                Ring {
                    points: [
                        [-h, -a],
                        [h, -a],
                        [h + a, 0.0],
                        [h + c, 0.0],
                        [-h - c, 0.0],
                        [-h - a, 0.0],
                    ]
                    .map(|p| self.place(p))
                    .to_vec(),
                    arcs: vec![None, arc(h), None, None, None, arc(-h)],
                }
            }
            // A stadium round the centre line.
            Plan::Slot { length, .. } => {
                let (h, a) = (length / 2.0, self.foot);
                let arc = |x: f64| {
                    Some(Curve::Arc(Arc {
                        center: self.place([x, 0.0]),
                        ccw: true,
                    }))
                };
                Ring {
                    points: [[-h, -a], [h, -a], [h, a], [-h, a]]
                        .map(|p| self.place(p))
                        .to_vec(),
                    arcs: vec![None, arc(h), None, arc(-h)],
                }
            }
            _ => unreachable!("a round or rectangular footprint"),
        }
    }

    /// The rectangle's corner radius, or the circle's radius.
    fn corner_radius(&self) -> f64 {
        match self.plan {
            Plan::Round(half) => half,
            Plan::Rect { rc, .. } => rc,
            Plan::Louver { .. } | Plan::Slot { .. } => 0.0,
        }
    }

    /// The outline's length, offset by s.
    fn perimeter(&self, s: f64) -> f64 {
        match self.plan {
            Plan::Round(half) => TAU * (half + s),
            Plan::Rect { hw, hl, rc, .. } => 4.0 * (hw + hl) - 8.0 * rc + TAU * (rc + s),
            Plan::Louver { length, .. } => length + PI * s,
            Plan::Slot { length, .. } => 2.0 * length + TAU * s,
        }
    }

    /// The area inside the outline offset by s.
    fn plan_area(&self, s: f64) -> f64 {
        match self.plan {
            Plan::Round(half) => PI * (half + s).powi(2),
            Plan::Rect { hw, hl, rc, .. } => {
                4.0 * (hw + s) * (hl + s) - (4.0 - PI) * (rc + s).powi(2)
            }
            Plan::Louver { .. } | Plan::Slot { .. } => 0.0,
        }
    }

    fn turn(&self) -> f64 {
        match self.plan {
            Plan::Rect { angle, .. } | Plan::Louver { angle, .. } | Plan::Slot { angle, .. } => {
                angle
            }
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
    /// 2 pi s; a louver's L + pi s), and the raised face's slab or a
    /// louver's lance faces. Base frame.
    pub fn mass(&self) -> Mass {
        let p0 = self.perimeter(0.0);
        let turn = self.perimeter(1.0) - p0;
        let (mut volume, mut mz, mut area) = (0.0, 0.0, 0.0);
        let (mut ms, mut ms2) = (0.0, 0.0);
        for i in 0..self.chains[0].len() {
            let [p, q] = [self.chains[0][i], self.chains[1][i]];
            let g = between(p, q);
            volume += p0 * g[0] + turn * g[1];
            mz += p0 * g[2] + turn * g[3];
            ms += g[1];
            ms2 += g[4];
            for m in [p, q] {
                area += length_moment(m, p0, turn);
            }
        }
        let mut center = self.center;
        match (self.plan, self.lance) {
            (Plan::Louver { length, .. }, Some(lance)) => {
                area += lance.area;
                // Towards the hood: the straight's moment length * s, each
                // end's s^2 (cos phi over a quarter turn).
                let back = (length * ms + 2.0 * ms2) / volume;
                let w = self.place_dir([0.0, 1.0]);
                center = [center[0] - back * w[0], center[1] - back * w[1]];
            }
            _ if self.crowned => {}
            _ if self.open.is_some() => {
                // The end face, square to the wall.
                let (p, q) = (self.chains[0].last(), self.chains[1].last());
                let (p, q) = (
                    p.expect("a section").at(1.0).0,
                    q.expect("a section").at(1.0).0,
                );
                area += length_moment(Meridian::Line { a: p, b: q }, p0, turn);
            }
            _ => {
                let slab = self.plan_area(self.top);
                let t = (self.tops[1] - self.tops[0]).abs();
                volume += slab * t;
                mz += slab * t * (self.tops[0] + self.tops[1]) / 2.0;
                area += 2.0 * slab;
            }
        }
        Mass {
            volume,
            area,
            center_of_mass: [center[0], center[1], mz / volume],
        }
    }

    /// Builds a round emboss's two faces in `b`, on the footprint's circles
    /// (`circles[level]`, each closed at vertex `v[level]`, starting along
    /// +x).
    pub fn build_round(
        &self,
        b: &mut Builder,
        m: Rigid,
        circles: [u32; 2],
        v: [u32; 2],
    ) -> Result<(), GeometryError> {
        let b = &mut Moved { b, m };
        let Plan::Round(half) = self.plan else {
            unreachable!()
        };
        let [cx, cy] = self.center;
        // A collar's last circle per level, and its section point.
        let mut ends: Vec<((u32, u32), P2)> = Vec::new();
        for level in 0..2 {
            let shift = |c: &Vec<Meridian>| {
                self.split(c)
                    .iter()
                    .map(|m| m.shifted(half))
                    .collect::<Vec<_>>()
            };
            let chain = shift(&self.chains[level]);
            let other = shift(&self.chains[1 - level]);
            // Circles and their vertices where the pieces meet, from the
            // footprint inward.
            let mut rings = vec![(circles[level], v[level])];
            let banded = chain.len() - usize::from(self.crowned);
            // A sharp corner (no radius) is no face: its two circles are one.
            let sharp = |m: &Meridian| matches!(m, Meridian::Arc { rho, .. } if *rho < EPS);
            for m in &chain[..banded] {
                if sharp(m) {
                    let last = *rings.last().expect("the footprint's circle");
                    rings.push(last);
                    continue;
                }
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
                if sharp(m) {
                    continue;
                }
                // Outward: away from the other face.
                let (here, there) = (m.at(0.5).0, other[i].at(0.5).0);
                let out = [here[0] - there[0], here[1] - there[1]];
                if i == banded {
                    // The crown: a cap inside the last circle.
                    let Meridian::Arc { c, rho, a1, .. } = *m else {
                        unreachable!("a crown is an arc")
                    };
                    let up = a1.sin().signum();
                    let (surface, samples) =
                        crown_sphere([cx, cy, c[1]], [0.0, 1.0, 0.0], up, rho, |f, g| {
                            let q = m.at(g).0;
                            let a = f * TAU;
                            [cx + q[0] * a.cos(), cy + q[0] * a.sin(), q[1]]
                        });
                    let (circle, vertex) = rings[i];
                    let q = m.at(0.5).0;
                    louver_face(
                        b,
                        &format!("{}.{level}.cap", self.id),
                        surface,
                        &[(circle, vertex, vertex, vertex)],
                        &samples,
                        [cx + q[0], cy, q[1]],
                        [out[0], 0.0, out[1]],
                    )?;
                    continue;
                }
                self.band(
                    b,
                    &format!("{}.{level}.{i}", self.id),
                    *m,
                    out,
                    rings[i],
                    rings[i + 1],
                )?;
            }
            if self.crowned {
                continue;
            }
            if self.open.is_some() {
                ends.push((rings[chain.len()], chain[chain.len() - 1].at(1.0).0));
                continue;
            }
            // The raised face.
            let up = (self.tops[level] - self.tops[1 - level]).signum();
            let frame = Frame3 {
                origin: [cx, cy, self.tops[level]],
                x: [1.0, 0.0, 0.0],
                y: [0.0, up, 0.0],
                z: [0.0, 0.0, up],
            };
            let (circle, vertex) = rings[chain.len()];
            let use_ = b.plane_boundary(frame, circle, vertex, vertex, Orientation::Forward)?;
            let r = half + self.top + 1.0;
            b.face(
                &format!("{}.{level}.top", self.id),
                SurfaceGeometry::Plane { frame },
                [[-r, r]; 2],
                vec![use_],
            )?;
        }

        // A collar's end: a ring (a plane, or a cone round a drafted wall)
        // between the two levels' last circles.
        if let [(ring_a, a), (ring_b, e)] = ends[..] {
            let up = (e[1] - a[1]).abs() > EPS;
            if up {
                // Outward: away from the material, along the wall.
                let (ma, mb) = (self.chains[0].last(), self.chains[1].last());
                let along = ma.expect("a section").at(1.0).1;
                let len = along[0].hypot(along[1]);
                let out = [along[0] / len, along[1] / len];
                let _ = mb;
                self.band(
                    b,
                    &format!("{}.end", self.id),
                    Meridian::Line { a, b: e },
                    out,
                    ring_a,
                    ring_b,
                )?;
            } else {
                let z = a[1];
                let along = self.chains[0].last().expect("a section").at(1.0).1;
                let sense = along[1].signum();
                let frame = Frame3 {
                    origin: [cx, cy, z],
                    x: [1.0, 0.0, 0.0],
                    y: [0.0, sense, 0.0],
                    z: [0.0, 0.0, sense],
                };
                let (outer, inner) = if a[0] > e[0] {
                    (ring_a, ring_b)
                } else {
                    (ring_b, ring_a)
                };
                let outer_use =
                    b.plane_boundary(frame, outer.0, outer.1, outer.1, Orientation::Forward)?;
                let inner_use =
                    b.plane_boundary(frame, inner.0, inner.1, inner.1, Orientation::Reverse)?;
                let r = a[0].max(e[0]) + 1.0;
                b.face_with_holes(
                    &format!("{}.end", self.id),
                    SurfaceGeometry::Plane { frame },
                    [[-r, r]; 2],
                    vec![outer_use],
                    vec![vec![inner_use]],
                )?;
            }
        }
        Ok(())
    }

    /// A full turn of a corner or the wall (in r, z) between the circles at
    /// its two ends, with a seam along it at angle 0.
    fn band(
        &self,
        b: &mut Moved,
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
                let frame = Frame3 {
                    origin: [cx, cy, c[1]],
                    x,
                    y: [0.0, 1.0, 0.0],
                    z: Z,
                };
                // About the axis itself (a dome's crown): a sphere, its
                // latitude the arc's angle.
                let surface = if c[0].abs() < EPS {
                    SurfaceGeometry::Sphere { frame, radius: rho }
                } else {
                    SurfaceGeometry::Torus {
                        frame,
                        major_radius: c[0],
                        minor_radius: rho,
                    }
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

    /// Builds a rectangular emboss's or a louver's faces in `b`, on the
    /// footprint's edges: per level, the vertex at each ring point and each
    /// segment's edge (id, start vertex).
    pub fn build_ring(
        &self,
        b: &mut Builder,
        m: Rigid,
        vertices: [Vec<u32>; 2],
        edges: [Vec<(u32, u32)>; 2],
    ) -> Result<(), GeometryError> {
        match self.plan {
            Plan::Louver { .. } => return self.build_louver(b, m, vertices, edges),
            Plan::Slot { .. } => return self.build_slot(b, m, vertices, edges),
            _ => {}
        }
        let vertices = std::array::from_fn(|l| std::array::from_fn(|j| vertices[l][j]));
        let edges = std::array::from_fn(|l| std::array::from_fn(|j| edges[l][j]));
        self.build_rect(b, m, vertices, edges)
    }

    /// A section to build: a crown whose edge lies low on its sphere in two,
    /// a band up to 60 degrees and the cap, so the cap's edge stays clear of
    /// its sphere's poles.
    fn split(&self, c: &[Meridian]) -> Vec<Meridian> {
        let mut c = c.to_vec();
        if let Some(&Meridian::Arc { c: k, rho, a0, a1 }) = c.last() {
            if self.crowned && a0.abs() < FRAC_PI_4 {
                let mid = a1.signum() * FRAC_PI_3;
                c.pop();
                c.push(Meridian::Arc {
                    c: k,
                    rho,
                    a0,
                    a1: mid,
                });
                c.push(Meridian::Arc {
                    c: k,
                    rho,
                    a0: mid,
                    a1,
                });
            }
        }
        c
    }

    /// Builds a bead's faces: per face of the sheet each piece of its section
    /// along both sides of the centre line (cylinders) and round each end (a
    /// half turn: a torus, or a sphere where it turns about the end's axis).
    /// The crown's sides meet at the centre line.
    fn build_slot(
        &self,
        b: &mut Builder,
        m: Rigid,
        vertices: [Vec<u32>; 2],
        edges: [Vec<(u32, u32)>; 2],
    ) -> Result<(), GeometryError> {
        let b = &mut Moved { b, m };
        let Plan::Slot { length, .. } = self.plan else {
            unreachable!()
        };
        let h = length / 2.0;
        let (u, w) = (self.place_dir([1.0, 0.0]), self.place_dir([0.0, 1.0]));
        let o = self.center;
        // A section point q (s, z) at p along the centre line on side tau
        // (+1: along w), or on the end at sigma turned phi from side -1.
        let side = |p: f64, tau: f64, q: P2| -> P3 {
            [
                o[0] + p * u[0] + tau * q[0] * w[0],
                o[1] + p * u[1] + tau * q[0] * w[1],
                q[1],
            ]
        };
        let end = |sigma: f64, phi: f64, q: P2| -> P3 {
            let (along, across) = (sigma * (h + q[0] * phi.sin()), -q[0] * phi.cos());
            [
                o[0] + along * u[0] + across * w[0],
                o[1] + along * u[1] + across * w[1],
                q[1],
            ]
        };
        type Step = (u32, u32, u32, u32);
        let rev = |path: &[Step]| -> Vec<Step> {
            path.iter()
                .rev()
                .map(|&(e, s, f, t)| (e, s, t, f))
                .collect()
        };
        let signs = [-1.0, 1.0];
        for level in 0..2 {
            let chain = &self.split(&self.chains[level]);
            let other = &self.split(&self.chains[1 - level]);
            let n = chain.len();
            let mut qs = vec![chain[0].at(0.0).0];
            qs.extend(chain.iter().map(|m| m.at(1.0).0));
            let ring = &vertices[level];
            // Vertices [sigma][tau][station]; at the last (s = 0) both sides
            // are one point.
            let mut at = [
                [vec![ring[0]], vec![ring[3]]],
                [vec![ring[1]], vec![ring[2]]],
            ];
            for (k, q) in qs.iter().enumerate().skip(1) {
                for (si, &sigma) in signs.iter().enumerate() {
                    let a = b.vertex(side(sigma * h, -1.0, *q));
                    let c = if k == n {
                        a
                    } else {
                        b.vertex(side(sigma * h, 1.0, *q))
                    };
                    at[si][0].push(a);
                    at[si][1].push(c);
                }
            }
            // Lines per side and station, half circles per end and station.
            let ring_step = |j: usize, from: u32, to: u32| -> Step {
                (edges[level][j].0, edges[level][j].1, from, to)
            };
            let mut lines: [Vec<Step>; 2] = [
                vec![ring_step(0, at[0][0][0], at[1][0][0])],
                vec![ring_step(2, at[0][1][0], at[1][1][0])],
            ];
            let mut halves: [Vec<Option<Step>>; 2] = [
                vec![Some(ring_step(3, at[0][0][0], at[0][1][0]))],
                vec![Some(ring_step(1, at[1][0][0], at[1][1][0]))],
            ];
            for (k, q) in qs.iter().enumerate().skip(1) {
                for (ti, &tau) in signs.iter().enumerate() {
                    // The crown runs across the centre line: no edge there.
                    if k == n {
                        continue;
                    }
                    let id = b.edge(
                        CurveGeometry::Line {
                            origin: side(-h, tau, *q),
                            direction: u,
                        },
                        Interval::new(0.0, length)?,
                        false,
                    );
                    lines[ti].push((id, at[0][ti][k], at[0][ti][k], at[1][ti][k]));
                }
                for (si, &sigma) in signs.iter().enumerate() {
                    if k == n {
                        halves[si].push(None);
                        continue;
                    }
                    let x = [-w[0], -w[1], 0.0];
                    let y = u.map(|v| sigma * v);
                    let id = b.edge(
                        CurveGeometry::Circle {
                            frame: Frame3 {
                                origin: end(sigma, 0.0, [0.0, q[1]]),
                                x,
                                y,
                                z: cross(x, y),
                            },
                            radius: q[0],
                        },
                        Interval::new(0.0, PI)?,
                        false,
                    );
                    halves[si].push(Some((id, at[si][0][k], at[si][0][k], at[si][1][k])));
                }
            }
            for (k, piece) in chain.iter().enumerate() {
                let Meridian::Arc { c, rho, a0, a1 } = *piece else {
                    unreachable!("a bead's section is arcs")
                };
                let (here, there) = (piece.at(0.5).0, other[k].at(0.5).0);
                let out = [here[0] - there[0], here[1] - there[1]];
                // Where the sides join the ends: per end and side, an arc in
                // the plane square to the centre line.
                let mut join = [[(0, 0, 0, 0); 2]; 2];
                for (si, &sigma) in signs.iter().enumerate() {
                    for (ti, &tau) in signs.iter().enumerate() {
                        let x = w.map(|v| tau * v);
                        let (from, to) = ((a0, at[si][ti][k]), (a1, at[si][ti][k + 1]));
                        let (lo, hi) = if from.0 < to.0 {
                            (from, to)
                        } else {
                            (to, from)
                        };
                        let id = b.edge(
                            CurveGeometry::Circle {
                                frame: Frame3 {
                                    origin: side(sigma * h, tau, [c[0], c[1]]),
                                    x,
                                    y: Z,
                                    z: cross(x, Z),
                                },
                                radius: rho,
                            },
                            Interval::new(lo.0, hi.0)?,
                            false,
                        );
                        join[si][ti] = (id, lo.1, from.1, to.1);
                    }
                }
                // The crown: one face across the centre line, its seam on
                // the far side.
                if k + 1 == n {
                    let x = [0.0, 0.0, -a1.sin().signum()];
                    let surface = SurfaceGeometry::Cylinder {
                        frame: Frame3 {
                            origin: side(0.0, 1.0, c),
                            x,
                            y: w,
                            z: cross(x, w),
                        },
                        radius: rho,
                    };
                    let steps = [
                        vec![lines[0][k], join[1][0]],
                        rev(&[lines[1][k], join[1][1]]),
                        vec![join[0][1]],
                        rev(&[join[0][0]]),
                    ]
                    .concat();
                    let samples: Vec<P3> = (0..=8)
                        .flat_map(|i| (0..=16).map(move |j| (i as f64 / 8.0, j as f64 / 16.0)))
                        .map(|(fp, fq)| {
                            let (tau, g) = if fq < 0.5 {
                                (-1.0, 2.0 * fq)
                            } else {
                                (1.0, 2.0 - 2.0 * fq)
                            };
                            side(-h + fp * length, tau, piece.at(g).0)
                        })
                        .collect();
                    let q = piece.at(1.0).0;
                    louver_face(
                        b,
                        &format!("{}.{level}.crown", self.id),
                        surface,
                        &steps,
                        &samples,
                        side(0.0, 1.0, q),
                        [0.0, 0.0, out[1]],
                    )?;
                }
                for (ti, &tau) in signs.iter().enumerate() {
                    if k + 1 == n {
                        break;
                    }
                    let x = w.map(|v| tau * v);
                    let surface = SurfaceGeometry::Cylinder {
                        frame: Frame3 {
                            origin: side(0.0, tau, c),
                            x,
                            y: Z,
                            z: cross(x, Z),
                        },
                        radius: rho,
                    };
                    let steps = [
                        vec![lines[ti][k], join[1][ti]],
                        rev(&[join[0][ti], lines[ti][k + 1]]),
                    ]
                    .concat();
                    let samples: Vec<P3> = (0..=8)
                        .flat_map(|i| (0..=8).map(move |j| (i as f64 / 8.0, j as f64 / 8.0)))
                        .map(|(fp, fq)| side(-h + fp * length, tau, piece.at(fq).0))
                        .collect();
                    louver_face(
                        b,
                        &format!("{}.{level}.{k}.side{ti}", self.id),
                        surface,
                        &steps,
                        &samples,
                        side(0.0, tau, piece.at(0.5).0),
                        [tau * out[0] * w[0], tau * out[0] * w[1], out[1]],
                    )?;
                }
                for (si, &sigma) in signs.iter().enumerate() {
                    let centre = end(sigma, 0.0, [0.0, c[1]]);
                    let (surface, samples) = if c[0].abs() > EPS || k + 1 < n {
                        let x = u.map(|v| sigma * v);
                        let frame = Frame3 {
                            origin: centre,
                            x,
                            y: cross(Z, x),
                            z: Z,
                        };
                        // A band of the crown's sphere about the end's axis.
                        let surface = if c[0].abs() > EPS {
                            SurfaceGeometry::Torus {
                                frame,
                                major_radius: c[0],
                                minor_radius: rho,
                            }
                        } else {
                            SurfaceGeometry::Sphere { frame, radius: rho }
                        };
                        let samples: Vec<P3> = (0..=8)
                            .flat_map(|i| (0..=8).map(move |j| (i as f64 / 8.0, j as f64 / 8.0)))
                            .map(|(f, g)| end(sigma, f * PI, piece.at(g).0))
                            .collect();
                        (surface, samples)
                    } else {
                        crown_sphere(centre, w, a1.sin().signum(), rho, |f, g| {
                            end(sigma, f * PI, piece.at(g).0)
                        })
                    };
                    let mut steps = Vec::new();
                    steps.extend(halves[si][k]);
                    steps.push(join[si][1]);
                    if let Some(q) = halves[si][k + 1] {
                        steps.extend(rev(&[q]));
                    }
                    steps.extend(rev(&[join[si][0]]));
                    let radial = u.map(|v| sigma * v);
                    louver_face(
                        b,
                        &format!("{}.{level}.{k}.end{si}", self.id),
                        surface,
                        &steps,
                        &samples,
                        end(sigma, FRAC_PI_2, piece.at(0.5).0),
                        [out[0] * radial[0], out[0] * radial[1], out[1]],
                    )?;
                }
            }
        }
        Ok(())
    }

    /// Builds a louver's faces: per face of the sheet each piece of its
    /// section along the straight (a cylinder) and round each end (a torus,
    /// or a sphere where it turns about the end's axis), and the two plane
    /// faces at the lance: the plate's side above the other face's level,
    /// and the hood's section below it.
    fn build_louver(
        &self,
        b: &mut Builder,
        m: Rigid,
        vertices: [Vec<u32>; 2],
        edges: [Vec<(u32, u32)>; 2],
    ) -> Result<(), GeometryError> {
        let b = &mut Moved { b, m };
        let Plan::Louver { length, .. } = self.plan else {
            unreachable!()
        };
        let lance = self.lance.expect("a louver's lance");
        let h = length / 2.0;
        let (u, w) = (self.place_dir([1.0, 0.0]), self.place_dir([0.0, 1.0]));
        let o = self.center;
        // A section point q (s, z) at p along the straight, or on the end
        // at sigma (+1: along u), turned phi from the straight towards the
        // lance.
        let straight = |p: f64, q: P2| -> P3 {
            [
                o[0] + p * u[0] - q[0] * w[0],
                o[1] + p * u[1] - q[0] * w[1],
                q[1],
            ]
        };
        let end = |sigma: f64, phi: f64, q: P2| -> P3 {
            let along = sigma * (h + q[0] * phi.sin());
            let back = q[0] * phi.cos();
            [
                o[0] + along * u[0] - back * w[0],
                o[1] + along * u[1] - back * w[1],
                q[1],
            ]
        };
        let minus_w = [-w[0], -w[1], 0.0];
        // A step round a loop: (edge, its start vertex, from, to).
        type Step = (u32, u32, u32, u32);
        let rev = |path: &[Step]| -> Vec<Step> {
            path.iter()
                .rev()
                .map(|&(e, s, f, t)| (e, s, t, f))
                .collect()
        };
        let sigmas = [-1.0, 1.0];
        // Per level: the faces' steps, kept for the lance faces.
        struct Level {
            /// Down each end's lance plane, from the footprint (sigma index).
            lance_path: [Vec<(usize, Vec<Step>)>; 2],
            /// The straight's last line (at the lance), and its ends.
            last: Step,
        }
        let mut levels = Vec::new();
        for level in 0..2 {
            let all = &self.chains[level];
            let other = &self.chains[1 - level];
            // Pieces of no size (a sharp inside corner) aren't faces.
            let pieces: Vec<(usize, Meridian)> = all
                .iter()
                .copied()
                .enumerate()
                .filter(|(_, m)| !matches!(m, Meridian::Arc { rho, .. } if *rho < EPS))
                .collect();
            let n = pieces.len();
            let mut qs = vec![pieces[0].1.at(0.0).0];
            qs.extend(pieces.iter().map(|(_, m)| m.at(1.0).0));
            let ring = &vertices[level];
            // Vertices per sigma and station: where it joins the straight,
            // and at the lance; at the last station (s = 0) one point.
            let mut join = [vec![ring[0]], vec![ring[1]]];
            let mut at_lance = [vec![ring[5]], vec![ring[2]]];
            for (k, q) in qs.iter().enumerate().skip(1) {
                for (si, &sigma) in sigmas.iter().enumerate() {
                    let j = b.vertex(end(sigma, 0.0, *q));
                    let l = if k == n {
                        j
                    } else {
                        b.vertex(end(sigma, FRAC_PI_2, *q))
                    };
                    join[si].push(j);
                    at_lance[si].push(l);
                }
            }
            // The straight's lines, and each end's quarter circles, per
            // station.
            let mut lines: Vec<Step> =
                vec![(edges[level][0].0, edges[level][0].1, join[0][0], join[1][0])];
            let mut quarters: [Vec<Option<Step>>; 2] = [
                vec![Some((
                    edges[level][5].0,
                    edges[level][5].1,
                    join[0][0],
                    at_lance[0][0],
                ))],
                vec![Some((
                    edges[level][1].0,
                    edges[level][1].1,
                    join[1][0],
                    at_lance[1][0],
                ))],
            ];
            for (k, q) in qs.iter().enumerate().skip(1) {
                let id = b.edge(
                    CurveGeometry::Line {
                        origin: straight(-h, *q),
                        direction: u,
                    },
                    Interval::new(0.0, length)?,
                    false,
                );
                lines.push((id, join[0][k], join[0][k], join[1][k]));
                for (si, &sigma) in sigmas.iter().enumerate() {
                    if k == n {
                        quarters[si].push(None);
                        continue;
                    }
                    let x = minus_w;
                    let y = u.map(|v| sigma * v);
                    let id = b.edge(
                        CurveGeometry::Circle {
                            frame: Frame3 {
                                origin: end(sigma, 0.0, [0.0, q[1]]),
                                x,
                                y,
                                z: cross(x, y),
                            },
                            radius: q[0],
                        },
                        Interval::new(0.0, FRAC_PI_2)?,
                        false,
                    );
                    quarters[si].push(Some((id, join[si][k], join[si][k], at_lance[si][k])));
                }
            }
            let mut lance_path: [Vec<(usize, Vec<Step>)>; 2] = [Vec::new(), Vec::new()];
            for (k, &(index, piece)) in pieces.iter().enumerate() {
                let Meridian::Arc { c, rho, a0, a1 } = piece else {
                    unreachable!("a louver's section is arcs")
                };
                let (here, there) = (piece.at(0.5).0, other[index].at(0.5).0);
                let out = [here[0] - there[0], here[1] - there[1]];
                // An arc of the section in a vertical plane, centred at the
                // section's centre placed by `at`, from `from` to `to` (angles
                // between a0 and a1, from x towards z) as one edge.
                let arc_edge =
                    |b: &mut Moved, centre: P3, x: P3, from: (f64, u32), to: (f64, u32)| {
                        let (lo, hi) = if from.0 < to.0 {
                            (from, to)
                        } else {
                            (to, from)
                        };
                        let id = b.edge(
                            CurveGeometry::Circle {
                                frame: Frame3 {
                                    origin: centre,
                                    x,
                                    y: Z,
                                    z: cross(x, Z),
                                },
                                radius: rho,
                            },
                            Interval::new(lo.0, hi.0)?,
                            false,
                        );
                        Ok::<Step, GeometryError>((id, lo.1, from.1, to.1))
                    };
                let mut join_m = Vec::new();
                for (si, &sigma) in sigmas.iter().enumerate() {
                    // Where it joins the straight.
                    let step = arc_edge(
                        b,
                        end(sigma, 0.0, c),
                        minus_w,
                        (a0, join[si][k]),
                        (a1, join[si][k + 1]),
                    )?;
                    join_m.push(step);
                    // In the lance's plane: split where the face that leaves
                    // the plate crosses its other face (a ring point there).
                    let x = u.map(|v| sigma * v);
                    let centre = end(sigma, FRAC_PI_2, c);
                    let (from, to) = ((a0, at_lance[si][k]), (a1, at_lance[si][k + 1]));
                    let path = if level == lance.out && index == lance.piece {
                        let crossing = (
                            a0 + lance.f * (a1 - a0),
                            vertices[1 - level][if sigma > 0.0 { 3 } else { 4 }],
                        );
                        vec![
                            arc_edge(b, centre, x, from, crossing)?,
                            arc_edge(b, centre, x, crossing, to)?,
                        ]
                    } else {
                        vec![arc_edge(b, centre, x, from, to)?]
                    };
                    lance_path[si].push((index, path));
                }
                // Along the straight: a cylinder about the section's centre.
                let x = minus_w;
                let surface = SurfaceGeometry::Cylinder {
                    frame: Frame3 {
                        origin: straight(0.0, c),
                        x,
                        y: Z,
                        z: cross(x, Z),
                    },
                    radius: rho,
                };
                let steps = [
                    vec![lines[k]],
                    vec![join_m[1]],
                    rev(&[lines[k + 1]]),
                    rev(&[join_m[0]]),
                ]
                .concat();
                let samples: Vec<P3> = (0..=8)
                    .flat_map(|i| (0..=8).map(move |j| (i as f64 / 8.0, j as f64 / 8.0)))
                    .map(|(fp, fq)| straight(-h + fp * length, piece.at(fq).0))
                    .collect();
                let outward = [out[0] * minus_w[0], out[0] * minus_w[1], out[1]];
                louver_face(
                    b,
                    &format!("{}.{level}.{k}.straight", self.id),
                    surface,
                    &steps,
                    &samples,
                    straight(0.0, piece.at(0.5).0),
                    outward,
                )?;
                // Round each end: a torus, or a sphere about the end's axis.
                for (si, &sigma) in sigmas.iter().enumerate() {
                    let x = unit3([sigma * u[0] - w[0], sigma * u[1] - w[1], 0.0]);
                    let surface = if c[0].abs() > EPS {
                        SurfaceGeometry::Torus {
                            frame: Frame3 {
                                origin: end(sigma, 0.0, [0.0, c[1]]),
                                x,
                                y: cross(Z, x),
                                z: Z,
                            },
                            major_radius: c[0],
                            minor_radius: rho,
                        }
                    } else {
                        // Its poles off the face, level with the centre.
                        let z = unit3([sigma * u[0] + w[0], sigma * u[1] + w[1], 0.0]);
                        SurfaceGeometry::Sphere {
                            frame: Frame3 {
                                origin: end(sigma, 0.0, c),
                                x,
                                y: cross(z, x),
                                z,
                            },
                            radius: rho,
                        }
                    };
                    let path = &lance_path[si].last().expect("just pushed").1;
                    let mut steps = Vec::new();
                    steps.extend(quarters[si][k]);
                    steps.extend(path.iter().copied());
                    if let Some(q) = quarters[si][k + 1] {
                        steps.extend(rev(&[q]));
                    }
                    steps.extend(rev(&[join_m[si]]));
                    let samples: Vec<P3> = (0..=8)
                        .flat_map(|i| (0..=8).map(move |j| (i as f64 / 8.0, j as f64 / 8.0)))
                        .map(|(fp, fq)| end(sigma, fp * FRAC_PI_2, piece.at(fq).0))
                        .collect();
                    let phi = FRAC_PI_2 / 2.0;
                    let radial = [
                        sigma * phi.sin() * u[0] - phi.cos() * w[0],
                        sigma * phi.sin() * u[1] - phi.cos() * w[1],
                    ];
                    let outward = [out[0] * radial[0], out[0] * radial[1], out[1]];
                    louver_face(
                        b,
                        &format!("{}.{level}.{k}.end{si}", self.id),
                        surface,
                        &steps,
                        &samples,
                        end(sigma, phi, piece.at(0.5).0),
                        outward,
                    )?;
                }
            }
            levels.push(Level {
                lance_path,
                last: lines[n],
            });
        }
        // The lance faces, in the plane through the lance square to w.
        let (out, inner) = (lance.out, 1 - lance.out);
        let ring = |level: usize, j: usize| -> Step {
            let (id, start) = edges[level][j];
            (id, start, vertices[level][j], vertices[level][(j + 1) % 6])
        };
        // The face leaving the plate, down each end's lance plane: above the
        // crossing and below it.
        let split = |si: usize| -> (Vec<Step>, Vec<Step>) {
            let (mut above, mut below) = (Vec::new(), Vec::new());
            for (index, path) in &levels[out].lance_path[si] {
                if *index < lance.piece {
                    above.extend(path.iter().copied());
                } else if *index == lance.piece {
                    above.push(path[0]);
                    below.push(path[1]);
                } else {
                    below.extend(path.iter().copied());
                }
            }
            (above, below)
        };
        let ((above_near, below_near), (above_far, below_far)) = (split(0), split(1));
        let side = [
            vec![ring(out, 2), ring(out, 3), ring(out, 4)],
            above_near,
            rev(&[ring(inner, 3)]),
            rev(&above_far),
        ]
        .concat();
        let whole = |level: usize, si: usize| -> Vec<Step> {
            levels[level].lance_path[si]
                .iter()
                .flat_map(|(_, p)| p.iter().copied())
                .collect()
        };
        let section = [
            vec![ring(inner, 2)],
            below_far,
            rev(&[levels[out].last]),
            rev(&below_near),
            vec![ring(inner, 4)],
            whole(inner, 0),
            vec![levels[inner].last],
            rev(&whole(inner, 1)),
        ]
        .concat();
        let frame = Frame3 {
            origin: [o[0], o[1], 0.0],
            x: u,
            y: Z,
            z: cross(u, Z),
        };
        let r = h
            + self.foot
            + self
                .chains
                .iter()
                .flatten()
                .map(|m| match m {
                    Meridian::Arc { c, rho, .. } => c[1].abs() + rho,
                    Meridian::Line { .. } => 0.0,
                })
                .fold(0.0, f64::max)
            + 1.0;
        // The plate's side faces -w (the frame's normal), the hood's +w.
        for (name, steps, sense) in [
            ("side", side, Orientation::Forward),
            ("section", section, Orientation::Reverse),
        ] {
            let mut uses = Vec::new();
            for (id, start, from, to) in steps {
                let along = if from == start {
                    Orientation::Forward
                } else {
                    Orientation::Reverse
                };
                uses.push(b.plane_boundary(frame, id, from, to, along)?);
            }
            let face = b.face(
                &format!("{}.lance.{name}", self.id),
                SurfaceGeometry::Plane { frame },
                [[-r, r]; 2],
                uses,
            )?;
            b.set_face_sense(face, sense)?;
        }
        Ok(())
    }

    /// Builds a rectangular emboss's two faces in `b`, on the footprint's
    /// edges: per level, the vertex at each of the eight segments' starts and
    /// each segment's edge (id, start vertex).
    fn build_rect(
        &self,
        b: &mut Builder,
        m: Rigid,
        vertices: [[u32; 8]; 2],
        edges: [[(u32, u32); 8]; 2],
    ) -> Result<(), GeometryError> {
        let b = &mut Moved { b, m };
        let Plan::Rect { rc, .. } = self.plan else {
            unreachable!()
        };
        // A collar's last vertices, ring edges and section end (point, along)
        // per level.
        type End = ([u32; 8], [(u32, u32); 8], (P2, P2));
        let mut ends: Vec<End> = Vec::new();
        for level in 0..2 {
            let chain = &self.chains[level];
            let other = &self.chains[1 - level];
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
                        uses.push(b.projected_boundary(
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
            if self.open.is_some() {
                let last = chain.len();
                ends.push((vs[last], rings[last], chain[last - 1].at(1.0)));
                continue;
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
                let last = chain.len();
                let (id, start) = rings[last][j];
                let (from, to) = (vs[last][j], vs[last][(j + 1) % 8]);
                let sense = if from == start {
                    Orientation::Forward
                } else {
                    Orientation::Reverse
                };
                uses.push(b.plane_boundary(frame, id, from, to, sense)?);
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
        if let [(v0, e0, (a, along)), (v1, e1, (q, _))] = ends[..] {
            let len = along[0].hypot(along[1]);
            let out = [along[0] / len, along[1] / len];
            let loop_uses = |b: &Moved, frame: Frame3, vs: [u32; 8], es: [(u32, u32); 8]| {
                (0..8)
                    .map(|j| {
                        let (id, start) = es[j];
                        let (from, to) = (vs[j], vs[(j + 1) % 8]);
                        let sense = if from == start {
                            Orientation::Forward
                        } else {
                            Orientation::Reverse
                        };
                        b.plane_boundary(frame, id, from, to, sense)
                    })
                    .collect::<Result<Vec<_>, _>>()
            };
            if (q[1] - a[1]).abs() < EPS {
                // Square to a straight wall: one plane ring.
                let up = out[1].signum();
                let x = self.place_dir([1.0, 0.0]);
                let frame = Frame3 {
                    origin: [self.center[0], self.center[1], a[1]],
                    x,
                    y: cross([0.0, 0.0, up], x),
                    z: [0.0, 0.0, up],
                };
                let (outer, inner) = if a[0] > q[0] {
                    ((v0, e0), (v1, e1))
                } else {
                    ((v1, e1), (v0, e0))
                };
                let outer = loop_uses(b, frame, outer.0, outer.1)?;
                let inner = loop_uses(b, frame, inner.0, inner.1)?;
                let Plan::Rect { hw, hl, .. } = self.plan else {
                    unreachable!()
                };
                let r = hw.max(hl) + q[0].abs().max(a[0].abs()) + 1.0;
                b.face_with_holes(
                    &format!("{}.end", self.id),
                    SurfaceGeometry::Plane { frame },
                    [[-r, r]; 2],
                    outer,
                    vec![inner],
                )?;
            } else {
                // Round a drafted wall: per segment a band, across at each
                // segment's start.
                let across: [(u32, u32); 8] = std::array::from_fn(|j| {
                    let (p0, p1) = (self.point3(j, 0.0, a), self.point3(j, 0.0, q));
                    let len = norm(sub(p1, p0));
                    let id = b.edge(
                        CurveGeometry::Line {
                            origin: p0,
                            direction: sub(p1, p0).map(|x| x / len),
                        },
                        Interval::new(0.0, len).unwrap(),
                        false,
                    );
                    (id, v0[j])
                });
                let m = Meridian::Line { a, b: q };
                for j in 0..8 {
                    let next = (j + 1) % 8;
                    let surface = self.piece_surface(j, m);
                    let mid = self.point3(j, 0.5, m.at(0.5).0);
                    let hint = surface.project(mid, None)?;
                    let mut uses = Vec::new();
                    for ((id, start), from, to) in [
                        (e0[j], v0[j], v0[next]),
                        (across[next], v0[next], v1[next]),
                        (e1[j], v1[next], v1[j]),
                        (across[j], v1[j], v0[j]),
                    ] {
                        let sense = if from == start {
                            Orientation::Forward
                        } else {
                            Orientation::Reverse
                        };
                        uses.push(b.projected_boundary(
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
                    let face = b.face(&format!("{}.end.{j}", self.id), surface, bounds, uses)?;
                    if dot(normal, outward) < 0.0 {
                        b.set_face_sense(face, Orientation::Reverse)?;
                    }
                }
            }
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

/// A builder that moves everything built into it by `m`: an emboss is
/// worked out in the flat pattern's frame and moved with its plate. Surface
/// parameters don't change, so projections onto the unmoved surface hold.
struct Moved<'a> {
    b: &'a mut Builder,
    m: Rigid,
}

impl Moved<'_> {
    fn frame(&self, f: Frame3) -> Frame3 {
        Frame3 {
            origin: self.m.apply(f.origin),
            x: self.m.rot(f.x),
            y: self.m.rot(f.y),
            z: self.m.rot(f.z),
        }
    }

    fn surface(&self, s: &SurfaceGeometry) -> SurfaceGeometry {
        let mut s = s.clone();
        match &mut s {
            SurfaceGeometry::Plane { frame }
            | SurfaceGeometry::Cylinder { frame, .. }
            | SurfaceGeometry::Cone { frame, .. }
            | SurfaceGeometry::Torus { frame, .. }
            | SurfaceGeometry::Sphere { frame, .. } => *frame = self.frame(*frame),
            _ => unreachable!("an emboss has no other surfaces"),
        }
        s
    }

    fn vertex(&mut self, p: P3) -> u32 {
        self.b.vertex(self.m.apply(p))
    }

    fn edge(&mut self, curve: CurveGeometry, range: Interval, seam: bool) -> u32 {
        let curve = match curve {
            CurveGeometry::Line { origin, direction } => CurveGeometry::Line {
                origin: self.m.apply(origin),
                direction: self.m.rot(direction),
            },
            CurveGeometry::Circle { frame, radius } => CurveGeometry::Circle {
                frame: self.frame(frame),
                radius,
            },
            _ => unreachable!("an emboss has no other curves"),
        };
        self.b.edge(curve, range, seam)
    }

    fn face(
        &mut self,
        key: &str,
        surface: SurfaceGeometry,
        bounds: [[f64; 2]; 2],
        uses: Vec<Use>,
    ) -> Result<u32, GeometryError> {
        let surface = self.surface(&surface);
        self.b.face(key, surface, bounds, uses)
    }

    fn set_face_sense(&mut self, face: u32, sense: Orientation) -> Result<(), GeometryError> {
        self.b.set_face_sense(face, sense)
    }

    fn face_with_holes(
        &mut self,
        key: &str,
        surface: SurfaceGeometry,
        bounds: [[f64; 2]; 2],
        outer: Vec<Use>,
        holes: Vec<Vec<Use>>,
    ) -> Result<u32, GeometryError> {
        let surface = self.surface(&surface);
        self.b.face_with_holes(key, surface, bounds, outer, holes)
    }

    fn plane_boundary(
        &self,
        frame: Frame3,
        edge: u32,
        from: u32,
        to: u32,
        sense: Orientation,
    ) -> Result<Use, GeometryError> {
        plane_boundary(self.b, self.frame(frame), edge, from, to, sense)
    }

    fn projected_boundary(
        &self,
        surface: &SurfaceGeometry,
        edge: u32,
        from: u32,
        to: u32,
        sense: Orientation,
        hint: Option<[f64; 2]>,
    ) -> Result<Use, GeometryError> {
        projected_boundary(self.b, &self.surface(surface), edge, from, to, sense, hint)
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

/// Between two pieces of the section (the region along p, across to q's
/// end, back along q, across to p's start): the integrals of 1, s, z, s z
/// and s^2 over it, by Green's theorem.
fn between(p: Meridian, q: Meridian) -> [f64; 5] {
    let across = |a: P2, b: P2| Meridian::Line { a, b };
    let mut g = [0.0; 5];
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
    g.map(|v| v * g[0].signum())
}

/// Over the piece, by Green's theorem, the integrals over the region it
/// bounds of 1, s, z, s z and s^2: of s dz, s^2 / 2 dz, s z dz, s^2 z / 2
/// dz, s^3 / 3 dz.
fn green(m: Meridian) -> [f64; 5] {
    let mut g = [0.0; 5];
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
            g[4] += s * s * s / 3.0 * w;
        }
    }
    g
}

/// The area the piece sweeps along the outline: the integral of (p0 + turn
/// s) along it (turn 2 pi round a closed outline).
fn length_moment(m: Meridian, p0: f64, turn: f64) -> f64 {
    let mut a = 0.0;
    let parts = 8;
    for k in 0..parts {
        for (node, w) in GAUSS {
            let (p, d) = m.at((k as f64 + node) / parts as f64);
            a += (p0 + turn * p[0]) * d[0].hypot(d[1]) * w / parts as f64;
        }
    }
    a
}

/// A louver's curved face: its loop's steps (edge, start, from, to) on the
/// surface, its bounds from the samples, turned to face `outward` at `mid`.
fn louver_face(
    b: &mut Moved,
    key: &str,
    surface: SurfaceGeometry,
    steps: &[(u32, u32, u32, u32)],
    samples: &[P3],
    mid: P3,
    outward: P3,
) -> Result<(), GeometryError> {
    let hint = surface.project(mid, None)?;
    let mut uses = Vec::new();
    for &(id, start, from, to) in steps {
        let sense = if from == start {
            Orientation::Forward
        } else {
            Orientation::Reverse
        };
        uses.push(b.projected_boundary(&surface, id, from, to, sense, Some(hint))?);
    }
    let mut bounds = [[f64::MAX, f64::MIN]; 2];
    for &p in samples {
        let uv = surface.project(p, Some(hint))?;
        for k in 0..2 {
            bounds[k][0] = bounds[k][0].min(uv[k]);
            bounds[k][1] = bounds[k][1].max(uv[k]);
        }
    }
    let bounds = bounds.map(|[lo, hi]| [lo - 1e-6, hi + 1e-6]);
    let normal = surface.normal_at(hint)?;
    let face = b.face(key, surface, bounds, uses)?;
    if dot(normal, outward) < 0.0 {
        b.set_face_sense(face, Orientation::Reverse)?;
    }
    Ok(())
}

/// A crown's sphere about `centre`, its poles along `axis` (level with the
/// centre, off the crown) and its seam on the far side from the crown
/// (`up`: +1 above the centre); with samples of the crown from `at`.
fn crown_sphere(
    centre: P3,
    axis: P3,
    up: f64,
    radius: f64,
    at: impl Fn(f64, f64) -> P3,
) -> (SurfaceGeometry, Vec<P3>) {
    let x = [0.0, 0.0, -up];
    let surface = SurfaceGeometry::Sphere {
        frame: Frame3 {
            origin: centre,
            x,
            y: cross(axis, x),
            z: axis,
        },
        radius,
    };
    let samples = (0..=16)
        .flat_map(|i| (0..=8).map(move |j| (i as f64 / 16.0, j as f64 / 8.0)))
        .map(|(f, g)| at(f, g))
        .collect();
    (surface, samples)
}

fn unit3(a: P3) -> P3 {
    let n = norm(a);
    a.map(|v| v / n)
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

    /// A louver on the 100 x 60 plate, the hood's middle on the lance at
    /// (53, 30), opening along +x.
    fn louver_part(t: f64, extra: &str) -> PartSpec {
        serde_json::from_str(&format!(
            r#"{{"id": "l", "thickness": {t}, "bend_radius": 2, "k_factor": 0.4,
            "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                "louvers": [{{"id": "v", "at": [53, 30], "length": 20, "radius": 6,
                    "bend_radius": {}{extra}}}]}}}}"#,
            2.0 - t
        ))
        .unwrap()
    }

    /// The reference CAD's louver tool: a hood of inside radius 6, 20 long
    /// between its ends' axes, leaving the plate over a corner of outside
    /// radius 2 (inside 0 on T 2, 1 on T 1).
    #[test]
    fn a_louver_is_an_exact_solid() {
        for (t, volume, area, center, faces) in [
            (
                2.0,
                12213.809432,
                12994.346914,
                [49.967314, 30.0, 0.808252],
                17,
            ),
            (
                1.0,
                6104.093121,
                12599.854554,
                [49.972392, 30.0, 0.316575],
                20,
            ),
        ] {
            let spec = louver_part(t, r#", "side": "down""#);
            let exact = crate::fold::exact_mass(&spec).unwrap();
            assert!(
                (exact.volume - volume).abs() < 1e-6,
                "T {t}: {}",
                exact.volume
            );
            assert!(
                (exact.area - area).abs() < 1e-6,
                "T {t}: area {}",
                exact.area
            );
            for j in 0..3 {
                let d = (exact.center_of_mass[j] - center[j]).abs();
                assert!(d < 1e-6, "T {t}: {:?}", exact.center_of_mass);
            }
            let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("T {t}: {e}"));
            assert_eq!(solid.topology.faces.len(), faces, "T {t}");
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(
                rel(mesh.volume, exact.volume) < 2e-5,
                "T {t}: {}",
                mesh.volume
            );
            assert!(
                rel(mesh.area, exact.area) < 2e-5,
                "T {t}: area {}",
                mesh.area
            );
        }
        // Up, opening along +y: the same, mirrored and turned.
        let spec = louver_part(2.0, r#", "side": "up", "angle": 90"#);
        let m = crate::fold::exact_mass(&spec).unwrap();
        assert!((m.volume - 12213.809432).abs() < 1e-6, "{}", m.volume);
        assert!((m.area - 12994.346914).abs() < 1e-6, "{}", m.area);
        // The louver's moment about the lance's middle (the down case's, less
        // the whole plate's) turns a quarter turn; the plate's stays.
        let (v, plate) = (12213.809432, 12000.0);
        let louver = v * (49.967314 - 53.0) - plate * (50.0 - 53.0);
        let want = [
            53.0 + plate * (50.0 - 53.0) / v,
            30.0 + louver / v,
            2.0 - 0.808252,
        ];
        for j in 0..3 {
            assert!(
                (m.center_of_mass[j] - want[j]).abs() < 1e-6,
                "{:?}",
                m.center_of_mass
            );
        }
        let solid = crate::fold::fold(&spec).unwrap();
        let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
        assert!(
            (mesh.volume - m.volume).abs() / m.volume < 2e-5,
            "{}",
            mesh.volume
        );
        // The flat keeps the plate whole and marks the footprint.
        let f = crate::flat::flat_pattern(&spec).unwrap();
        assert!((f.area() - 6000.0).abs() < 1e-9 && f.cutouts.is_empty());
        assert_eq!(f.forms[0].outline.points.len(), 6);
    }

    /// A large bend radius: the face leaving the plate crosses the other
    /// face's level on its corner, not on the hood. Also rising, turned.
    #[test]
    fn a_louver_crossing_on_its_corner() {
        for (side, angle) in [("down", 0.0), ("up", 30.0)] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{"id": "l", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                    "louvers": [{{"id": "v", "at": [50, 30], "length": 12, "radius": 6,
                        "bend_radius": 10, "side": "{side}", "angle": {angle}}}]}}}}"#
            ))
            .unwrap();
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{side}: {e}"));
            assert_eq!(solid.topology.faces.len(), 20, "{side}");
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(
                rel(mesh.volume, exact.volume) < 2e-5,
                "{side}: {} vs {}",
                mesh.volume,
                exact.volume
            );
            assert!(
                rel(mesh.area, exact.area) < 2e-5,
                "{side}: area {} vs {}",
                mesh.area,
                exact.area
            );
            for j in 0..3 {
                let d = (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs();
                assert!(
                    d < 1e-3,
                    "{side}: {:?} vs {:?}",
                    mesh.center_of_mass,
                    exact.center_of_mass
                );
            }
        }
    }

    /// A louver on a flange moves with it, up or down, turned: the solid
    /// meshes to its exact mass.
    #[test]
    fn a_louver_on_a_flange_moves_with_it() {
        for (angle, fside, side, turn) in [
            (90.0, "up", "up", 0.0),
            (90.0, "down", "down", 90.0),
            (120.0, "up", "down", 200.0),
        ] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{"id": "l", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                    "louvers": [{{"id": "v", "on": "f", "at": [50, 25], "length": 12,
                        "radius": 6, "bend_radius": 0, "side": "{side}", "angle": {turn}}}]}},
                "flanges": [{{"id": "f", "on": {{"base": 0}}, "angle": {angle}, "length": 50,
                    "datum": "outer", "position": "outside", "side": "{fside}"}}]}}"#
            ))
            .unwrap();
            let what = format!("{angle} {fside} {side} {turn}");
            let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{what}: {e}"));
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(rel(mesh.volume, exact.volume) < 2e-5, "{what}: volume");
            assert!(rel(mesh.area, exact.area) < 2e-5, "{what}: area");
            for j in 0..3 {
                let d = (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs();
                assert!(
                    d < 1e-3,
                    "{what}: {:?} vs {:?}",
                    mesh.center_of_mass,
                    exact.center_of_mass
                );
            }
        }
    }

    /// The reference CAD's dimple and single rib tools on a 200 x 200 plate,
    /// T 2, pressed down at (100, 100): a dome (crown inside radius 30, 5
    /// deep, corner inside 8) and a bead along y (crown 14.5, 80 between its
    /// ends, 4 deep, corner 2). Up mirrors them.
    #[test]
    fn domes_and_beads_are_exact_solids() {
        let part = |shape: &str, height: f64, radius: f64, side: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "d", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "embosses": [{{"id": "e", "at": [100, 100], "shape": {shape},
                        "height": {height}, "radius": {radius}, "side": "{side}", "angle": 90}}]}}}}"#
            ))
            .unwrap()
        };
        let dome = r#"{"dome": {"radius": 30}}"#;
        let bead = r#"{"bead": {"length": 80, "radius": 14.5}}"#;
        for (what, shape, height, radius, volume, area, z, faces) in [
            (
                "dome",
                dome,
                5.0,
                8.0,
                80132.906683,
                81732.906683,
                0.937962,
                10,
            ),
            (
                "bead",
                bead,
                4.0,
                2.0,
                80379.402501,
                81979.402501,
                0.860153,
                20,
            ),
        ] {
            for side in ["down", "up"] {
                let spec = part(shape, height, radius, side);
                let exact = crate::fold::exact_mass(&spec).unwrap();
                assert!(
                    (exact.volume - volume).abs() < 1e-6,
                    "{what}: {}",
                    exact.volume
                );
                assert!(
                    (exact.area - area).abs() < 1e-6,
                    "{what}: area {}",
                    exact.area
                );
                let want = if side == "down" { z } else { 2.0 - z };
                let c = exact.center_of_mass;
                assert!(
                    (c[0] - 100.0).abs() < 1e-9 && (c[1] - 100.0).abs() < 1e-9,
                    "{what}: {c:?}"
                );
                assert!((c[2] - want).abs() < 1e-6, "{what} {side}: {c:?}");
                let solid =
                    crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{what} {side}: {e}"));
                assert_eq!(solid.topology.faces.len(), faces, "{what} {side}");
                let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
                let rel = |a: f64, b: f64| (a - b).abs() / b;
                assert!(
                    rel(mesh.volume, exact.volume) < 2e-5,
                    "{what} {side}: {}",
                    mesh.volume
                );
                assert!(
                    rel(mesh.area, exact.area) < 2e-5,
                    "{what} {side}: area {}",
                    mesh.area
                );
            }
        }
    }

    #[test]
    fn what_a_dome_or_bead_refuses() {
        let part = |shape: &str, height: f64, draft: f64| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "d", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "embosses": [{{"id": "e", "at": [100, 100], "shape": {shape},
                        "height": {height}, "draft": {draft}, "side": "up"}}]}}}}"#
            ))
            .unwrap()
        };
        let invalid = |spec: &PartSpec| {
            matches!(crate::flat::flat_pattern(spec), Err(SheetError::Invalid(_)))
        };
        let dome = r#"{"dome": {"radius": 10}}"#;
        // No wall to draft; so high it would turn back under itself.
        assert!(invalid(&part(dome, 5.0, 10.0)));
        assert!(invalid(&part(dome, 15.0, 0.0)));
        assert!(!invalid(&part(dome, 12.0, 0.0)));
        // A bead needs a length.
        assert!(invalid(&part(
            r#"{"bead": {"length": 0, "radius": 10}}"#,
            5.0,
            0.0
        )));
    }

    /// A deep dome or bead: the crown's edge low on its sphere, built as a
    /// band and a cap; still its exact mass.
    #[test]
    fn deep_domes_and_beads_split_their_crowns() {
        for (shape, faces) in [
            (r#"{"dome": {"radius": 6}}"#, 12),
            (r#"{"bead": {"length": 30, "radius": 6}}"#, 28),
        ] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{"id": "d", "thickness": 0.8, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "embosses": [{{"id": "e", "at": [100, 100], "shape": {shape},
                        "height": 8, "radius": 4, "side": "down", "angle": 35}}]}}}}"#
            ))
            .unwrap();
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{shape}: {e}"));
            assert_eq!(solid.topology.faces.len(), faces, "{shape}");
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(
                rel(mesh.volume, exact.volume) < 2e-5,
                "{shape}: {}",
                mesh.volume
            );
            assert!(
                rel(mesh.area, exact.area) < 2e-5,
                "{shape}: area {}",
                mesh.area
            );
        }
    }

    /// The reference CAD's round flange and counter sink emboss tools on a
    /// 200 x 200 plate, T 2, pressed down at (100, 100): a straight collar
    /// (inside diameter 12.7, 4 high, corner 0.5) and one drafted 20 degrees
    /// (12 at the top face, 5 high, corner 2). Up mirrors them.
    #[test]
    fn flanged_holes_are_exact_solids() {
        let part = |hole: &str, side: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "h", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "flanged_holes": [{{"id": "c", "at": [100, 100], {hole}, "side": "{side}"}}]}}}}"#
            ))
            .unwrap()
        };
        let straight = r#""shape": {"round": {"diameter": 12.7}}, "height": 4, "radius": 0.5"#;
        let drafted =
            r#""shape": {"round": {"diameter": 12}}, "height": 5, "radius": 2, "draft": 20"#;
        for (what, hole, volume, area, z) in [
            ("straight", straight, 79876.004328, 81559.989571, 0.995),
            ("drafted", drafted, 80017.298994, 81676.12454, 0.990825),
        ] {
            for side in ["down", "up"] {
                let spec = part(hole, side);
                let exact = crate::fold::exact_mass(&spec).unwrap();
                assert!(
                    (exact.volume - volume).abs() < 1e-6,
                    "{what}: {}",
                    exact.volume
                );
                assert!(
                    (exact.area - area).abs() < 1e-6,
                    "{what}: area {}",
                    exact.area
                );
                let want = if side == "down" { z } else { 2.0 - z };
                assert!(
                    (exact.center_of_mass[2] - want).abs() < 1e-6,
                    "{what} {side}: {:?}",
                    exact.center_of_mass
                );
                let solid =
                    crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{what} {side}: {e}"));
                assert_eq!(solid.topology.faces.len(), 11, "{what} {side}");
                let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
                let rel = |a: f64, b: f64| (a - b).abs() / b;
                assert!(
                    rel(mesh.volume, exact.volume) < 2e-5,
                    "{what} {side}: {}",
                    mesh.volume
                );
                assert!(
                    rel(mesh.area, exact.area) < 2e-5,
                    "{what} {side}: area {}",
                    mesh.area
                );
            }
        }
        // The flat: the pre-hole (the straight collar developed at K 0.4:
        // 8.85 less 1.3 pi / 2 and the 1.5 wall), the footprint marked.
        let f = crate::flat::flat_pattern(&part(straight, "down")).unwrap();
        assert_eq!(f.holes.len(), 1);
        let pre = 2.0 * (8.85 - 1.3 * std::f64::consts::FRAC_PI_2 - 1.5);
        assert!(
            (f.holes[0].diameter - pre).abs() < 1e-9,
            "{}",
            f.holes[0].diameter
        );
        assert_eq!(f.forms.len(), 1);
    }

    /// The reference CAD's rectangular flange on the 200 x 200 plate, T 2:
    /// a straight collar round a 20 x 30 hole with plan corners 2.5, 4 high,
    /// corner 1; and a drafted one, turned, meshed to its exact mass.
    #[test]
    fn rectangular_flanged_holes_are_exact_solids() {
        let part = |hole: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "h", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "flanged_holes": [{{"id": "c", "at": [100, 100], {hole}}}]}}}}"#
            ))
            .unwrap()
        };
        let rect = r#""shape": {"rectangle": {"width": 20, "length": 30, "corner_radius": 2.5}}"#;
        let spec = part(&format!(
            r#"{rect}, "height": 4, "radius": 1, "side": "down""#
        ));
        let exact = crate::fold::exact_mass(&spec).unwrap();
        assert!(
            (exact.volume - 79049.24779).abs() < 1e-6,
            "{}",
            exact.volume
        );
        assert!((exact.area - 80844.852507).abs() < 1e-6, "{}", exact.area);
        assert!(
            (exact.center_of_mass[2] - 0.98848).abs() < 1e-6,
            "{:?}",
            exact.center_of_mass
        );
        let solid = crate::fold::fold(&spec).unwrap();
        assert_eq!(solid.topology.faces.len(), 39);
        let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
        assert!(
            (mesh.volume - exact.volume).abs() / exact.volume < 2e-5,
            "{}",
            mesh.volume
        );
        // The flat: a rounded rectangle pre-hole, the footprint marked.
        let f = crate::flat::flat_pattern(&spec).unwrap();
        assert_eq!((f.holes.len(), f.cutouts.len(), f.forms.len()), (0, 1, 1));
        let drafted = part(&format!(
            r#"{rect}, "height": 6, "radius": 1.5, "draft": 15, "side": "up", "angle": 30"#
        ));
        let exact = crate::fold::exact_mass(&drafted).unwrap();
        let solid = crate::fold::fold(&drafted).unwrap();
        assert_eq!(solid.topology.faces.len(), 46);
        let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
        let rel = |a: f64, b: f64| (a - b).abs() / b;
        assert!(rel(mesh.volume, exact.volume) < 2e-5, "{}", mesh.volume);
        assert!(rel(mesh.area, exact.area) < 2e-5, "{}", mesh.area);
        for j in 0..3 {
            let d = (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs();
            assert!(
                d < 1e-3,
                "{:?} vs {:?}",
                mesh.center_of_mass,
                exact.center_of_mass
            );
        }
    }

    #[test]
    fn what_a_flanged_hole_refuses() {
        let part = |hole: &str| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "h", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "flanged_holes": [{{"id": "c", "at": [100, 100], {hole}, "side": "up"}}]}}}}"#
            ))
            .unwrap()
        };
        let invalid = |spec: &PartSpec| {
            matches!(crate::flat::flat_pattern(spec), Err(SheetError::Invalid(_)))
        };
        // No wall above the corner; a wall closing the hole; too tall to
        // develop from the plate.
        assert!(invalid(&part(
            r#""shape": {"round": {"diameter": 12}}, "height": 3, "radius": 2"#
        )));
        assert!(invalid(&part(
            r#""shape": {"round": {"diameter": 6}}, "height": 20, "radius": 1, "draft": 20"#
        )));
        assert!(invalid(&part(
            r#""shape": {"round": {"diameter": 6}}, "height": 20, "radius": 1"#
        )));
        assert!(!invalid(&part(
            r#""shape": {"round": {"diameter": 6}}, "height": 6, "radius": 1"#
        )));
    }

    /// A collar with a sharp inside corner (radius 0): that corner is an
    /// edge, not a face.
    #[test]
    fn a_sharp_flanged_hole() {
        for (draft, side, faces) in [(0.0, "down", 10), (15.0, "up", 10)] {
            let spec: PartSpec = serde_json::from_str(&format!(
                r#"{{"id": "h", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [200, 0], [200, 200], [0, 200]],
                    "flanged_holes": [{{"id": "c", "at": [100, 100], "shape": {{"round": {{"diameter": 20}}}}, "height": 6,
                        "radius": 0, "draft": {draft}, "side": "{side}"}}]}}}}"#
            ))
            .unwrap();
            let exact = crate::fold::exact_mass(&spec).unwrap();
            let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{side}: {e}"));
            assert_eq!(solid.topology.faces.len(), faces, "{side}");
            let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
            let rel = |a: f64, b: f64| (a - b).abs() / b;
            assert!(
                rel(mesh.volume, exact.volume) < 2e-5,
                "{side}: {}",
                mesh.volume
            );
            assert!(
                rel(mesh.area, exact.area) < 2e-5,
                "{side}: area {}",
                mesh.area
            );
        }
    }

    #[test]
    fn what_a_louver_refuses() {
        let louver = |t: f64, at: [f64; 2], length: f64, radius: f64| -> PartSpec {
            serde_json::from_str(&format!(
                r#"{{"id": "l", "thickness": {t}, "bend_radius": 2, "k_factor": 0.4,
                "base": {{"outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
                    "louvers": [{{"id": "v", "at": {at:?}, "length": {length},
                        "radius": {radius}, "side": "down"}}]}}}}"#
            ))
            .unwrap()
        };
        let invalid = |spec: &PartSpec| {
            matches!(crate::flat::flat_pattern(spec), Err(SheetError::Invalid(_)))
        };
        // A hood no deeper than the sheet doesn't open; no length.
        assert!(invalid(&louver(2.0, [50.0, 30.0], 20.0, 2.0)));
        // A corner closing on itself at the ends (R 3, r 6 on T 2).
        let mut spec = louver(2.0, [50.0, 30.0], 20.0, 3.0);
        spec.base.louvers[0].bend_radius = Some(6.0);
        assert!(invalid(&spec));
        assert!(invalid(&louver(2.0, [50.0, 30.0], 0.0, 6.0)));
        // Its footprint across the outline.
        assert!(crate::flat::flat_pattern(&louver(2.0, [50.0, 5.0], 60.0, 6.0)).is_err());
        assert!(crate::flat::flat_pattern(&louver(2.0, [8.0, 30.0], 20.0, 6.0)).is_err());
        // Its default bend radius is the part's.
        assert!(crate::flat::flat_pattern(&louver(2.0, [50.0, 30.0], 20.0, 6.0)).is_ok());
    }

    /// A drafted wall of no length: the corners meet. The reference CAD's
    /// circular emboss tool on T 2: inside radii 3, 5 deep, the corners
    /// tangent (draft asin 3/8), the raised face r 22 at its inside.
    #[test]
    fn corners_that_meet_leave_no_wall() {
        let draft = (3.0f64 / 8.0).asin().to_degrees();
        let diameter = 44.0 + 50.0 / 55f64.sqrt();
        let tool = |side: &str| {
            format!(
                r#"{{"id": "c", "at": [50, 30], "shape": {{"round": {{"diameter": {diameter}}}}},
                "height": 5, "draft": {draft}, "radius": 3, "side": "{side}"}}"#
            )
        };
        exact_solid(&[(tool("down"), "down"), (tool("up"), "up")]);
        let spec = part(&tool("down"), "");
        let m = crate::fold::exact_mass(&spec).unwrap();
        assert!((m.volume - 12670.343650).abs() < 1e-6, "{}", m.volume);
        assert!((m.area - 13310.343650).abs() < 1e-6, "{}", m.area);
        assert!(
            (m.center_of_mass[2] + 0.773639).abs() < 1e-6,
            "{:?}",
            m.center_of_mass
        );
        let solid = crate::fold::fold(&spec).unwrap();
        assert_eq!(solid.topology.faces.len(), 12);
        exact_solid(&[
            (rect_tangent(5.0, 3.0, "up"), "rectangle up"),
            (rect_tangent(5.0, 3.0, "down"), "rectangle down"),
        ]);
    }

    /// A rectangle whose corners meet, for height h and inside radius r on
    /// T 2: the draft where the wall's length comes to nothing.
    fn rect_tangent(h: f64, r: f64, side: &str) -> String {
        // The corners' centres are 2r + T apart, (2r + T - h) apart in z.
        let d = ((2.0 * r + 2.0 - h) / (2.0 * r + 2.0)).asin().to_degrees();
        format!(
            r#"{{"id": "r", "at": [50, 30], "shape": {{"rectangle": {{"width": 50, "length": 36,
            "corner_radius": 16}}}}, "height": {h}, "draft": {d}, "radius": {r}, "side": "{side}"}}"#
        )
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
        // Too small across for its corners, too low (at 3 its corners just
        // meet), across the outline.
        assert!(invalid(&round([50.0, 30.0], 8.0, 6.0, 30.0, "up"), ""));
        assert!(invalid(&round([50.0, 30.0], 30.0, 2.5, 30.0, "up"), ""));
        assert!(invalid(&round([50.0, 10.0], 30.0, 6.0, 30.0, "up"), ""));
        // Plan corners too tight for the section, or past half a side; and
        // a rectangle across the outline.
        assert!(invalid(&rect([40.0, 24.0, 3.0], 5.0, 20.0, "up", 0.0), ""));
        assert!(invalid(&rect([40.0, 24.0, 13.0], 5.0, 20.0, "up", 0.0), ""));
        let across =
            crate::flat::flat_pattern(&part(&rect([98.0, 24.0, 8.0], 5.0, 20.0, "up", 0.0), ""));
        assert!(matches!(across, Err(SheetError::Invalid(_))), "{across:?}");
        // Over a bend; on the base or the flange (in the flat) it builds.
        let flange = r#"{"id": "f", "on": {"base": 2}, "angle": 90, "length": 50,
            "datum": "outer", "position": "outside", "side": "up"}"#;
        for y in [30.0, 85.0] {
            let spec = part(&round([50.0, y], 30.0, 6.0, 30.0, "up"), flange);
            crate::fold::fold(&spec).unwrap();
        }
        assert!(matches!(
            crate::flat::flat_pattern(&part(&round([50.0, 55.0], 30.0, 6.0, 30.0, "up"), flange)),
            Err(SheetError::Unsupported(_))
        ));
        // On a flange, reaching into its bend.
        let on = round([50.0, 12.0], 30.0, 6.0, 30.0, "up")
            .replace(r#""id": "d","#, r#""id": "d", "on": "f","#);
        assert!(matches!(
            crate::flat::flat_pattern(&part(&on, flange)),
            Err(SheetError::Unsupported(_))
        ));
    }

    /// An emboss on a flange (edge 0, along +x; up from the outer mold line)
    /// moves with it: the solid meshes to its exact mass, and the raised
    /// face lies `height` off the flange's face on its side, centred where
    /// it was put.
    #[test]
    fn an_emboss_on_a_flange_moves_with_it() {
        let flange = |angle: f64, side: &str| {
            format!(
                r#"{{"id": "f", "on": {{"base": 0}}, "angle": {angle}, "length": 50,
                "datum": "outer", "position": "outside", "side": "{side}"}}"#
            )
        };
        let round = |side: &str| {
            format!(
                r#"{{"id": "d", "on": "f", "at": [50, 25], "shape": {{"round": {{"diameter": 24}}}},
                "height": 8, "draft": 20, "side": "{side}"}}"#
            )
        };
        let rect = |side: &str| {
            format!(
                r#"{{"id": "r", "on": "f", "at": [50, 25], "shape": {{"rectangle": {{"width": 30,
                "length": 20, "corner_radius": 8}}}}, "height": 8, "draft": 20, "side": "{side}",
                "angle": 0}}"#
            )
        };
        // Flange angle and side, emboss side; for the square flanges the
        // raised face's plane (y) and centre (x, z).
        let cases = [
            // Up: the face's normal points out (-y); the outer face at y = -2.
            (90.0, "up", "up", Some((-10.0, 25.0))),
            (90.0, "up", "down", Some((8.0, 25.0))),
            // Down: normal in (+y), the inner face at y = 0; up from the top.
            (90.0, "down", "up", Some((8.0, -23.0))),
            (90.0, "down", "down", Some((-10.0, -23.0))),
            (60.0, "up", "up", None),
            (120.0, "down", "down", None),
        ];
        for (angle, fside, eside, expect) in cases {
            for emboss in [round(eside), rect(eside)] {
                let spec = part(&emboss, &flange(angle, fside));
                let what = format!("{angle} {fside} {eside} {}", &emboss[8..9]);
                let solid = crate::fold::fold(&spec).unwrap_or_else(|e| panic!("{what}: {e}"));
                let exact = crate::fold::exact_mass(&spec).unwrap();
                let mesh = crate::fold::mass(&solid, 2e-4).unwrap();
                let rel = |a: f64, b: f64| (a - b).abs() / b;
                assert!(rel(mesh.volume, exact.volume) < 2e-5, "{what}: volume");
                assert!(rel(mesh.area, exact.area) < 2e-5, "{what}: area");
                for j in 0..3 {
                    let d = (mesh.center_of_mass[j] - exact.center_of_mass[j]).abs();
                    assert!(d < 1e-3, "{what}: centre");
                }
                let Some((y, z)): Option<(f64, f64)> = expect else {
                    continue;
                };
                // The raised face: the vertices furthest along y that way,
                // clear of the base plate.
                let m = opengeometry::tessellation::tessellate(&solid, 1e-3, 10_000_000).unwrap();
                let pts: Vec<[f64; 3]> = m
                    .positions
                    .chunks_exact(3)
                    .map(|p| [p[0], p[1], p[2]])
                    .filter(|p| (p[2] - z).abs() < 20.0)
                    .collect();
                let far = pts
                    .iter()
                    .map(|p| p[1] * y.signum())
                    .fold(f64::MIN, f64::max)
                    * y.signum();
                assert!((far - y).abs() < 1e-9, "{what}: raised face at y = {far}");
                let top: Vec<_> = pts.iter().filter(|p| (p[1] - y).abs() < 1e-9).collect();
                for (k, want) in [(0, 50.0), (2, z)] {
                    let lo = top.iter().map(|p| p[k]).fold(f64::MAX, f64::min);
                    let hi = top.iter().map(|p| p[k]).fold(f64::MIN, f64::max);
                    assert!(((lo + hi) / 2.0 - want).abs() < 1e-6, "{what}: centre {k}");
                }
            }
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
