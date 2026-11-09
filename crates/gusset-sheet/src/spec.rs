//! The part description: Gusset's input format. Lengths in millimetres, angles in degrees.
//!
//! Frame: the base outline lies in the XY plane and the base plate occupies
//! z = 0 .. thickness. `Up` flanges bend toward +z (from the top face edge),
//! `Down` flanges toward -z (from the bottom face edge).

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartSpec {
    pub id: String,
    pub thickness: f64,
    /// Inner bend radius, used by every bend.
    pub bend_radius: f64,
    /// Linear unfold: neutral layer at `k_factor * thickness` from the inside.
    pub k_factor: f64,
    /// Bend compensation from a table instead, for every bend without a
    /// `k_factor` of its own.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_table: Option<BendTable>,
    /// Bend reliefs for flanges narrower than their edge. Default: width T,
    /// depth T/2, minimum remnant 2T.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relief: Option<Relief>,
    pub base: Base,
    #[serde(default)]
    pub flanges: Vec<Flange>,
    /// Contour flanges, built after the flanges.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contours: Vec<Contour>,
    /// Hems, built after the flanges (so a hem can sit on a flange's tip).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub hems: Vec<Hem>,
    /// Bends along a line across the base.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub folds: Vec<FoldLine>,
    /// Z-steps across the base: two folds each.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub jogs: Vec<Jog>,
    /// Holes and cut-outs placed on a flange's face instead of in the flat.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub face_cuts: Vec<FaceCuts>,
    /// Corner seams between two flanges built as separate features (with
    /// different settings), which otherwise meet by build order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub corner_seams: Vec<CornerSeam>,
    /// Faces on planes of their own, each joined to the base or an earlier
    /// face by one of `bends`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub faces: Vec<Face>,
    /// Bends joining an edge of the base or a face to an edge of another
    /// face, in build order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bends: Vec<Bend>,
    /// Bends to make first, in this order, by the names the flat pattern
    /// gives them (a feature's id; its further bends `id.2`, `id.3`, ...). The
    /// rest follow in the default order: the bends furthest from the base
    /// first, then as listed in the flat pattern.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bend_order: Vec<String>,
}

/// A corner seam between two flanges on neighbouring edges that differ in
/// length, height datum or side: they meet as flanges with the same settings
/// do, each mitered by its own bend, or overlapping.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CornerSeam {
    pub flanges: [String; 2],
    #[serde(default, skip_serializing_if = "Seam::is_miter")]
    pub seam: Seam,
    /// Default: the relief's miter gap.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gap: Option<f64>,
}

/// A face: a plate on a plane of its own (in the base frame), joined to the
/// part by a bend.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Face {
    pub id: String,
    /// A point on the plane.
    pub origin: [f64; 3],
    /// The plane's x and y directions: unit length, square to each other.
    pub x_axis: [f64; 3],
    pub y_axis: [f64; 3],
    /// Counter-clockwise polygon in the plane's x, y, seen from its normal
    /// (x cross y). Edge `i` runs from `outline[i]` to `outline[i + 1]`.
    pub outline: Vec<[f64; 2]>,
    /// Per outline edge: 0 for a line, otherwise an arc's bulge, as on the
    /// base. Bends sit on lines only. All lines when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bulges: Vec<f64>,
    /// The sheet lies on the normal's side of the plane (`up`) or the other.
    pub material: Side,
    /// Holes and cut-outs in the plane's x, y.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Hole>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cutouts: Vec<Cutout>,
}

/// A bend from an edge of the base or of a face joined earlier (`from`) to
/// an edge of a face (`to`), which joins that face to the part. The edges
/// must be parallel; the bend runs where they overlap along it. Each face
/// is extended or trimmed, square to the bend, to where the bend meets it;
/// a face trimmed beside the bend gets the part's bend reliefs there, where
/// it carries on past the bend.
///
/// Faces at an angle get one bend, on the inside of the corner their sheets
/// make. Parallel faces get two (`double`).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Bend {
    pub id: String,
    pub from: FaceEdge,
    pub to: FaceEdge,
    /// Between parallel faces. Default: `deg45` for faces side by side,
    /// `deg90` for faces that face each other.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub double: Option<DoubleBend>,
    /// The edge that stays where it is in a 45 or 90 degree or full radius
    /// double bend; the other face is extended or trimmed. Default: `from`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fixed: Option<FixedEdge>,
    /// Inner bend radius for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_radius: Option<f64>,
    /// K-factor for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_factor: Option<f64>,
}

/// Edge `edge` of the base's outline (face "base") or of a face's.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FaceEdge {
    pub face: String,
    pub edge: usize,
}

/// Two bends between parallel faces.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum DoubleBend {
    /// Faces side by side: both edges stay, two equal bends of whatever
    /// angle joins them.
    FixedEdges,
    /// Faces side by side: two bends of 45 degrees.
    Deg45,
    /// Faces facing each other: two bends of 90 degrees and a straight wall.
    Deg90,
    /// Faces facing each other: one bend of 180 degrees, its inner radius
    /// half the gap between them.
    FullRadius,
}

/// Which of a bend's edges stays.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum FixedEdge {
    From,
    To,
}

/// Holes and cut-outs on a flange's face, cut square through it. Positions
/// are [along, up]: along the edge the flange stands on (from its start, as
/// in the flat), and along the face from the outer mold line, where the
/// flange's outer face meets the plane of the outer face it is bent from (for
/// a 90-degree flange: the height over that face). Cut-outs run
/// counter-clockwise in [along, up]. They stay on the flat of the face, clear
/// of its bends; across its tip or ends they are notches.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FaceCuts {
    /// The flange (or a flange on its tip) whose face they are on.
    pub on: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Hole>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cutouts: Vec<Cutout>,
}

/// A Z-step: the part to the right of `line[0] -> line[1]` (seen from +z)
/// steps up or down by `offset` and runs on parallel. Two folds of `angle`:
/// the first starts at the line, the second turns back after a straight
/// run, so that the faces move by `offset` (bottom faces for up, top faces
/// for down: (2R + T)(1 - cos angle) + run sin angle).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Jog {
    pub id: String,
    pub line: [[f64; 2]; 2],
    pub offset: f64,
    /// Default 90.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub angle: Option<f64>,
    pub side: Side,
}

/// Bend compensation by table (DIN 6935): a bend
/// unfolds to A + B - x, where A and B are the two legs measured to the bend's
/// outside (to the outer faces' virtual sharp for an opening angle of 90
/// degrees or more, else to the outer face's tangent square to the other leg).
/// One table, for the part's thickness; x is interpolated linearly between
/// rows and columns, and a value within 0.004 of the table's range counts as
/// its end.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BendTable {
    /// Inner bend radii, ascending.
    pub radii: Vec<f64>,
    /// Opening angles (180 minus the bend angle), ascending, degrees.
    pub angles: Vec<f64>,
    /// x per opening angle (a row each), per radius.
    pub values: Vec<Vec<f64>>,
}

impl BendTable {
    fn correction(&self, angle: f64, radius: f64) -> Option<f64> {
        // Index of the interval holding v, and how far along it.
        let find = |xs: &[f64], v: f64| -> Option<(usize, usize, f64)> {
            const TOL: f64 = 0.004;
            let n = xs.len();
            if n == 0 || v < xs[0] - TOL || v > xs[n - 1] + TOL {
                return None;
            }
            let v = v.clamp(xs[0], xs[n - 1]);
            let i = (0..n.saturating_sub(1))
                .find(|&i| v <= xs[i + 1])
                .unwrap_or(0);
            let j = (i + 1).min(n - 1);
            let f = if j == i {
                0.0
            } else {
                (v - xs[i]) / (xs[j] - xs[i])
            };
            Some((i, j, f))
        };
        let (i0, i1, fi) = find(&self.angles, angle)?;
        let (j0, j1, fj) = find(&self.radii, radius)?;
        let v = |i: usize, j: usize| self.values.get(i).and_then(|row| row.get(j)).copied();
        let row = |i: usize| Some((1.0 - fj) * v(i, j0)? + fj * v(i, j1)?);
        Some((1.0 - fi) * row(i0)? + fi * row(i1)?)
    }
}

/// A bend along a straight line across the base, or across a flange, hem or
/// contour strip in the flat (it folds that strip): the part to the right of
/// `line[0] -> line[1]` (seen from +z) turns, the rest stays. On a strip the
/// side away from the strip's own bend turns, whichever way the line runs.
/// Uses the part's bend radius.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FoldLine {
    pub id: String,
    pub line: [[f64; 2]; 2],
    pub angle: f64,
    pub location: BendLocation,
    pub side: Side,
    /// Inner bend radius for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_radius: Option<f64>,
    /// K-factor for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_factor: Option<f64>,
    /// Bend only the sheet the line crosses between its two points (up to
    /// 180 degrees); sheet beside them stays, or turns with what it sits on.
    /// Otherwise the line runs on across the part.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub bounded: bool,
}

/// Where the bend sits relative to a fold line.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum BendLocation {
    /// The line is the middle of the bend.
    Centerline,
    /// The bend starts at the line (on the turning side).
    Start,
    /// The bend ends at the line.
    End,
}

/// A profile run along a base edge (a contour flange):
/// every corner of the profile becomes a bend of the part's radius.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Contour {
    pub id: String,
    /// A base edge, or a flange's tip. On a tip the section frame is the
    /// parent's: s continues along the parent, z runs from the parent's outer
    /// face (outside its last bend: "bottom") toward its inner face ("top").
    pub on: EdgeRef,
    /// Partial width on a base edge, with bend reliefs as for a flange; the
    /// full edge when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub width: Option<Width>,
    /// An open polyline in the edge's section: s outward from the edge (along
    /// its outward normal), z up from the base's bottom face. Its first line
    /// leaves the base face it continues; where that line meets the face, the
    /// first bend's corner is. The start point may be anywhere on the line.
    pub profile: Vec<[f64; 2]>,
    /// Per segment: 0 for a line, otherwise tan(sweep / 4) of an arc (as a
    /// DXF bulge; positive turns counter-clockwise, seen with s right and z
    /// up). An arc must be tangent to its neighbours; it is a bend whose face
    /// on the profile is the arc. All lines when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bulges: Vec<f64>,
    /// The base face the profile continues. The profile is that face of the
    /// contour flange; the material lies on the base's side of it.
    pub face: BaseFace,
    /// Inner bend radius for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_radius: Option<f64>,
    /// K-factor for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_factor: Option<f64>,
    /// Further base edges the same contour runs along, each next to the one
    /// before, mitered at the corners between them. Full width only.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub along: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum BaseFace {
    Top,
    Bottom,
}

/// An edge folded back over itself.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Hem {
    pub id: String,
    pub on: EdgeRef,
    /// Partial width on a base edge, with bend reliefs as for a flange; the
    /// full edge when absent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub width: Option<Width>,
    pub kind: HemKind,
    /// The side it folds over: as for flanges.
    pub side: Side,
    /// K-factor for its bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_factor: Option<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum HemKind {
    /// One 180-degree fold: inner radius gap/2, `length` from the edge to the end.
    Single { gap: f64, length: f64 },
    /// Folded twice into three layers with `gap` between them; `length` from
    /// the edge to the far fold.
    Double { gap: f64, length: f64 },
    /// One bend of `angle` (180..360) around `radius`, closed by a straight
    /// tail that ends T/4 above the sheet.
    Teardrop { radius: f64, angle: f64 },
    /// One bend of `angle` around `radius`, no tail.
    Rolled { radius: f64, angle: f64 },
}

/// A bend relief: a notch next to each end of a partial flange, outside its
/// width.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Relief {
    /// Straight (rectangular) by default.
    #[serde(default, skip_serializing_if = "ReliefShape::is_straight")]
    pub shape: ReliefShape,
    /// Notch width along the edge.
    pub width: f64,
    /// How far the notch reaches past the start of the bend, into the base.
    pub depth: f64,
    /// A notch leaving no more than this of the base edge to its corner runs
    /// on to the neighbouring edge instead (next to an arc the remnant stays).
    pub minimum_remnant: f64,
    /// Where two flanges meet at a corner: the gap between their inner
    /// corners. Default T.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub miter_gap: Option<f64>,
    /// How a corner where two flanges meet is relieved. Default: trim to bend.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub corner: Option<CornerRelief>,
    /// How two flanges meeting at a corner close it.
    /// Default: mitered.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seam: Option<Seam>,
    /// How a strip's free end next to an edge that isn't square to it (or an
    /// arc at an angle) shows. Default: none.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transition: Option<Transition>,
}

/// A bend transition: the end of a strip next to an edge that isn't square
/// to it. With none (or a straight line) the flange ends square at the
/// vertex and the bend in a straight cut from where its start line meets the
/// neighbouring edge to the vertex at its far edge.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Transition {
    #[default]
    None,
    /// The same as none.
    StraightLine,
    /// The flat shows the strip ending square through the vertex (the
    /// folded part as with none).
    Intersection,
    /// As intersection, the corner at the vertex rounded with the bend
    /// radius.
    Arc,
    /// The bend ends square where it first meets the neighbouring edge or
    /// the vertex, the flange at the vertex; the base runs along the bend
    /// start line to the neighbouring edge. Flat and folded.
    TrimToBend,
}

impl Seam {
    fn is_miter(&self) -> bool {
        *self == Seam::Miter
    }
}

/// How two flanges meeting at a corner close it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Seam {
    /// Each flange stops half the miter gap from the corner's symmetry plane.
    #[default]
    Miter,
    /// The flange on the edge that ends at the corner (in outline order) runs
    /// on over the other's end; the other stops the miter gap short of it.
    Overlap,
    /// The flange on the edge that starts at the corner runs on over the
    /// other's end.
    ReverseOverlap,
}

/// The relief at a corner where two flanges meet.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CornerRelief {
    pub shape: CornerShape,
    /// Diameter (round) or side (square). Default 4T. Larger than the
    /// diagonal of where the bends overlap, centred where the bend centre
    /// lines cross. Smaller, on the bisector through q (where the bend start
    /// lines cross): its inscribed circle through q, until its corner reaches
    /// the circle through q and both flanges' first corners, then at that
    /// circle's centre; until then it is cut from the flat
    /// pattern only (folded, the corner is trimmed to bend).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<f64>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum CornerShape {
    /// Each bend runs to the other's bend start; the flanges stop the miter
    /// gap apart.
    #[default]
    TrimToBend,
    /// Trim to bend, minus a disc.
    Round,
    /// Trim to bend, minus a square: its diagonal on the bisector of the two
    /// edge normals (along the edges at a square corner).
    Square,
    /// Each bend ends in a straight cut from where both bends start to the
    /// corner of its flange; the flat shows that cut (a weld preparation).
    LinearWeld,
    /// As the linear weld in 3D; the flat shows the corner closed, its flange
    /// edges meeting (torn when folded).
    Tear,
}

/// The notch's end, past the start of the bend.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ReliefShape {
    /// Square end at `depth` past the bend start.
    #[default]
    Straight,
    /// Semicircular end of diameter `width`, reaching `depth` past the bend
    /// start (so `depth` must be at least `width / 2`).
    Round,
    /// Torn, not cut: a 0.01 x 0.01 mm notch in 3D (whatever `width` and
    /// `depth` say), none in the flat.
    Tear,
}

impl ReliefShape {
    fn is_straight(&self) -> bool {
        *self == ReliefShape::Straight
    }
}

impl PartSpec {
    /// A bend's length in the flat and the K-factor that gives it: the
    /// feature's own K-factor, else the bend table, else the part's.
    pub(crate) fn bend_allowance(
        &self,
        angle: f64,
        radius: f64,
        own_k: Option<f64>,
    ) -> Result<(f64, f64), crate::SheetError> {
        let (t, theta) = (self.thickness, angle.to_radians());
        let Some(table) = self.bend_table.as_ref().filter(|_| own_k.is_none()) else {
            let k = own_k.unwrap_or(self.k_factor);
            return Ok((theta * (radius + k * t), k));
        };
        let x = table.correction(180.0 - angle, radius).ok_or_else(|| {
            crate::SheetError::Invalid(format!(
                "a bend of {angle} degrees with radius {radius} is outside the bend table"
            ))
        })?;
        // The legs run to the bend's outside: its virtual sharp, or for a
        // bend past 90 degrees the outer face's tangent.
        let leg = if angle <= 90.0 + 1e-9 {
            (radius + t) * (theta / 2.0).tan()
        } else {
            radius + t
        };
        let allowance = 2.0 * leg - x;
        Ok((allowance, (allowance / theta - radius) / t))
    }

    pub fn relief(&self) -> Relief {
        self.relief.unwrap_or(Relief {
            shape: ReliefShape::Straight,
            width: self.thickness,
            depth: self.thickness / 2.0,
            minimum_remnant: 2.0 * self.thickness,
            miter_gap: None,
            corner: None,
            seam: None,
            transition: None,
        })
    }
}

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Base {
    /// Counter-clockwise polygon seen from +z. Edge `i` runs from
    /// `outline[i]` to `outline[i + 1]` (wrapping). Empty for a profile.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub outline: Vec<[f64; 2]>,
    /// Per outline edge: 0 for a line, otherwise tan(sweep / 4) of an arc (as
    /// a DXF bulge; positive turns counter-clockwise). No feature sits on an
    /// arc edge. A flange running into a tangent arc (a rounded corner) gets
    /// a relief whose notch runs on to the arc; one meeting an arc at an
    /// angle ends its bend where the bend's start line meets the arc. All
    /// lines when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bulges: Vec<f64>,
    /// In place of an outline: the part is a profile run along a width.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub profile: Option<Profile>,
    /// In place of an outline: the part is a rolled cone or cylinder segment.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cone: Option<Cone>,
    /// Holes in the flat pattern (the base frame, which the flat pattern
    /// shares): on the base, across or inside a bend, or on a flange. Each
    /// must lie inside the flat outline or cross it (a round notch), without
    /// touching it, and clear of the others.
    #[serde(default)]
    pub holes: Vec<Hole>,
    /// Cut-outs in the flat pattern, like holes but any shape.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cutouts: Vec<Cutout>,
    /// Punched shapes, stamped in the flat pattern like holes and cut-outs.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub punches: Vec<Punch>,
    /// Formed embosses on the base plate.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub embosses: Vec<Emboss>,
}

/// An emboss: the sheet pressed out of the base plate at its thickness,
/// around a centre in the flat pattern. Seen from the side it rises out of
/// the plate through a rounded corner, along a straight wall, and over a
/// second rounded corner into its raised face, `height` above the plate's
/// face on that side. Both corners have `radius` on their inside. The flat
/// pattern keeps the plate whole (the emboss stretches it) and marks the
/// footprint, where the sheet leaves the plate.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Emboss {
    pub id: String,
    /// The centre, in the flat pattern.
    pub at: [f64; 2],
    pub shape: EmbossShape,
    pub height: f64,
    /// The wall's slope from square to the plate, degrees: 0 for a straight
    /// wall, below 90.
    #[serde(default)]
    pub draft: f64,
    /// Inside radius of both corners. Default: the part's bend radius.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub radius: Option<f64>,
    /// `up`: out of the plate's top face (+z); `down`: out of its bottom.
    pub side: Side,
    /// A rectangle's turn, counter-clockwise, degrees.
    #[serde(default)]
    pub angle: f64,
}

/// An emboss's raised face, measured where the wall meets it (its mold
/// line).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum EmbossShape {
    Round {
        diameter: f64,
    },
    /// `width` along x, `length` along y (before the turn), corners rounded
    /// in plan with `corner_radius`, which must be large enough for the
    /// section to turn round them.
    Rectangle {
        width: f64,
        length: f64,
        corner_radius: f64,
    },
}

/// A punch: holes and cut-outs drawn around (0, 0), stamped in the flat
/// pattern at each point and turned counter-clockwise by `angle` degrees. The
/// flat pattern lists every stamp (id, centre, angle) for the punch press.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Punch {
    /// The punch's ID, as the flat pattern marks it.
    pub id: String,
    /// The stamps' centres in the flat pattern.
    pub at: Vec<[f64; 2]>,
    #[serde(default)]
    pub angle: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub holes: Vec<Hole>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cutouts: Vec<Cutout>,
}

impl Punch {
    /// Its holes and cut-outs at one stamp's centre.
    pub fn stamped(&self, at: [f64; 2]) -> (Vec<Hole>, Vec<Cutout>) {
        let (s, c) = self.angle.to_radians().sin_cos();
        let place = |p: [f64; 2]| [at[0] + c * p[0] - s * p[1], at[1] + s * p[0] + c * p[1]];
        let holes = self
            .holes
            .iter()
            .map(|h| Hole {
                center: place(h.center),
                diameter: h.diameter,
            })
            .collect();
        let cutouts = self
            .cutouts
            .iter()
            .map(|k| Cutout {
                outline: k.outline.iter().map(|&p| place(p)).collect(),
                bulges: k.bulges.clone(),
            })
            .collect();
        (holes, cutouts)
    }
}

/// A part made from a profile run along a width (a contour flange as the
/// first feature). The profile is one face of the sheet; every
/// sharp corner is a bend of the part's radius and a tangent arc a bend with
/// the arc's radius, as for a contour flange on an edge. Its first straight
/// line on z = 0 with the material above it, or on z = T with the material
/// below it, is the base plate (z = 0..T): the flat pattern, holes and
/// cut-outs use the base's frame, as for an outline. Flanges and hems may sit
/// on its free ends, `{"tip": "start"}` and `{"tip": "end"}` (`up` turns as
/// the profile's last bend before that end; an end without a bend is the base
/// plate's edge there), on each other's tips, and on the
/// base plate's edges: 0 at y = lo, 1 at x = width, 2 at y = hi, 3 at x = 0
/// (the profile is built first: those meet its bends by build order).
/// Contour flanges and folds go on it as on an outline.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Profile {
    /// An open polyline in the YZ plane: [y, z].
    pub points: Vec<[f64; 2]>,
    /// Per segment: 0 for a line, otherwise tan(sweep / 4) of an arc (as a
    /// DXF bulge; positive turns counter-clockwise, seen with y right and z
    /// up). All lines when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bulges: Vec<f64>,
    /// The side of the profile the material lies on, going along the points
    /// seen with y right and z up.
    pub material: Hand,
    /// Run along +x from x = 0.
    pub width: f64,
}

/// A part rolled to a cone or cylinder segment around the z axis: the sheet
/// between two coaxial arcs, one on z = 0 and one on z = height, both from +x
/// counter-clockwise by `sweep` (a lofted flange between two such arcs, or a
/// straight line rolled around an axis). The arcs lie on one face of the
/// sheet; the sheet's ends are cut square to it, so an end of a cone is a
/// cone too. The flat pattern is the developed neutral surface (at the
/// K-factor): the edge from the arcs' start points runs from (0, 0) to (0, L),
/// L the slant length, and the sheet unrolls towards +x, seen from outside.
/// Holes and features on it aren't supported yet.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Cone {
    /// The arcs' radii: on z = 0, then on z = height. Equal for a cylinder.
    pub radii: [f64; 2],
    pub height: f64,
    /// Degrees, below 360.
    pub sweep: f64,
    /// The side of the arcs' face the material lies on: away from the axis
    /// or towards it.
    pub material: Radial,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Radial {
    Outside,
    Inside,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Hand {
    Left,
    Right,
}

/// A cut-out: a closed profile of lines and arcs, counter-clockwise seen
/// from +z, cut through the flat sheet. Segment `k` runs from `outline[k]`
/// to `outline[k + 1]` (wrapping).
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Cutout {
    pub outline: Vec<[f64; 2]>,
    /// Per segment: 0 for a line, otherwise tan(sweep / 4) of an arc (as a
    /// DXF bulge; positive turns counter-clockwise). All lines when empty.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub bulges: Vec<f64>,
}

/// A circular hole, cut through the flat sheet.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Hole {
    pub center: [f64; 2],
    pub diameter: f64,
}

/// A straight flange along the edge it is built on: the full length, or a
/// part of it on a base edge. For a flange on a flange, "base" in the datum
/// and position docs means the parent flange.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Flange {
    pub id: String,
    pub on: EdgeRef,
    /// Partial width on a base edge; the full edge when absent. A flange on
    /// a tip always spans its parent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub width: Option<Width>,
    /// Bend angle: how far the flange turns away from the base plane.
    /// 0 continues the base, 90 is perpendicular, 135 leans back over the base.
    pub angle: f64,
    /// Flange height, measured as `datum` says.
    pub length: f64,
    pub datum: HeightDatum,
    pub position: BendPosition,
    pub side: Side,
    /// Both corners of its free tip rounded or chamfered, where the flange's ends are square to it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tip_corners: Option<CornerCut>,
    /// Inner bend radius for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bend_radius: Option<f64>,
    /// K-factor for this feature's bends, in place of the part's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub k_factor: Option<f64>,
}

/// A cut across a sheet corner.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum CornerCut {
    Round {
        radius: f64,
    },
    /// The same distance back along both edges.
    Chamfer {
        distance: f64,
    },
}

/// Where the flange height is measured from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum HeightDatum {
    /// From the outside of the base, along the flange.
    Outer,
    /// From the inside of the base, along the flange.
    Inner,
    /// From the tangent at the end of the bend, along the flange.
    Tangent,
    /// From the outside of the base, perpendicular to the base.
    OuterOrtho,
    /// From the inside of the base, perpendicular to the base.
    InnerOrtho,
}

/// Where the bend sits relative to the base edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum BendPosition {
    /// The flange's outer face stays within the base outline.
    Inside,
    /// The flange's inner face is flush with the base edge.
    Outside,
}

/// Where a partial flange sits on its base edge `a -> b`.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct Width {
    /// Gap between `a` (outline[i]) and the flange.
    pub offset_start: f64,
    /// Gap between the flange and `b` (outline[i + 1]).
    pub offset_end: f64,
}

/// The edge a flange is built on.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum EdgeRef {
    /// Base outline edge `i`.
    Base(usize),
    /// The free end of the flange with this id (defined earlier in the list),
    /// of a contour flange (its profile's end), or of a profile part: "start"
    /// or "end".
    Tip(String),
    /// Edge `edge` of the face with this id (or of "base").
    Face(FaceEdge),
}

/// Which way the flange bends. On a base edge: `Up` toward +z. On a tip:
/// `Up` turns the same way as the parent flange (a return lip, C-profile),
/// `Down` the opposite way (a Z-step).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Side {
    Up,
    Down,
}

#[cfg(test)]
mod tests {
    use super::*;

    const SCHEMA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../schema/part.schema.json");

    /// The committed JSON Schema matches the types. Regenerate with
    /// `UPDATE_SCHEMA=1 cargo test -p gusset-sheet`.
    #[test]
    fn json_schema_is_current() {
        let schema = serde_json::to_string_pretty(&schemars::schema_for!(PartSpec)).unwrap() + "\n";
        if std::env::var("UPDATE_SCHEMA").is_ok() {
            std::fs::write(SCHEMA, &schema).unwrap();
        }
        let committed = std::fs::read_to_string(SCHEMA)
            .unwrap_or_default()
            .replace("\r\n", "\n");
        assert_eq!(
            committed, schema,
            "schema out of date: rerun with UPDATE_SCHEMA=1"
        );
    }

    #[test]
    fn a_part_round_trips_through_json() {
        let json = r#"{
            "id": "l-bracket",
            "thickness": 2.0, "bend_radius": 2.0, "k_factor": 0.4,
            "base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                      "holes": [{ "center": [20, 30], "diameter": 8 }] },
            "flanges": [{ "id": "f1", "on": { "base": 2 }, "angle": 90, "length": 30,
                          "datum": "outer", "position": "outside", "side": "up" }]
        }"#;
        let part: PartSpec = serde_json::from_str(json).unwrap();
        assert_eq!(part.flanges[0].datum, HeightDatum::Outer);
        let again: PartSpec = serde_json::from_str(&serde_json::to_string(&part).unwrap()).unwrap();
        assert_eq!(part, again);
    }
}
