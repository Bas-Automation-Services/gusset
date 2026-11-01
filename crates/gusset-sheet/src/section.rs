//! The cross-section of one straight flange, perpendicular to its edge.
//!
//! Work in (s, z): s along the base's outward edge normal (s = 0 at the base
//! edge), z across the sheet with the outer plane at z = 0 and the inner plane at
//! z = T. The bend's centre is at (setback, T + R) and the flange leaves in the
//! direction (cos θ, sin θ). A down flange is the same section mirrored, so the
//! numbers don't depend on the side.

use crate::spec::{BendPosition, Flange, HeightDatum};
use crate::SheetError;

/// Where the flange sits in the flat layout, measured along the edge normal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Section {
    /// Start of the bend relative to the base edge (negative: inside the base).
    pub setback: f64,
    /// Length of the bend strip in the flat: θ (R + k T).
    pub bend_allowance: f64,
    /// Length of the flat flange beyond the bend.
    pub straight: f64,
}

impl Section {
    /// Where the flat flange ends, relative to the base edge.
    pub fn reach(&self) -> f64 {
        self.setback + self.bend_allowance + self.straight
    }
}

pub fn section(
    thickness: f64,
    radius: f64,
    k: f64,
    flange: &Flange,
) -> Result<Section, SheetError> {
    let (t, r) = (thickness, radius);
    if !(flange.angle > 0.0 && flange.angle < 180.0) {
        return Err(SheetError::Invalid(format!(
            "flange {}: bend angle must be between 0 and 180 degrees",
            flange.id
        )));
    }
    let theta = flange.angle.to_radians();
    let (sin, cos) = theta.sin_cos();
    let half_tan = (theta / 2.0).tan();

    // The flange face that ends flush with the base edge: inner for `outside`
    // (its virtual sharp with the inner plane at s = 0), outer for `inside`
    // (it crosses the inner plane at s = 0: its sharp with the outer plane is
    // T cot θ in; the same at 90 degrees).
    let setback = match flange.position {
        BendPosition::Outside => -r * half_tan,
        BendPosition::Inside => -(r + t) * half_tan - t * cos / sin,
    };
    let l = flange.length;
    let straight = match flange.datum {
        // Along the flange from the virtual sharp of its outer / inner face.
        HeightDatum::Outer => l - (r + t) * half_tan,
        HeightDatum::Inner => l - r * half_tan,
        // Along the flange from the plane normal to it that touches the outer bend.
        HeightDatum::Tangent => l - (r + t),
        // Height of the outer / inner face end above the matching base plane.
        HeightDatum::OuterOrtho => (l - t - r + (r + t) * cos) / sin,
        HeightDatum::InnerOrtho => (l - r + r * cos) / sin,
    };
    if straight <= 0.0 {
        return Err(SheetError::Invalid(format!(
            "flange {}: length {l} leaves no straight flange after the bend",
            flange.id
        )));
    }
    Ok(Section {
        setback,
        bend_allowance: theta * (r + k * t),
        straight,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::Side;

    fn flange(angle: f64, datum: HeightDatum, position: BendPosition) -> Flange {
        Flange {
            id: "f".into(),
            on: crate::spec::EdgeRef::Base(0),
            width: None,
            angle,
            length: 30.0,
            datum,
            position,
            side: Side::Up,
            tip_corners: None,
            bend_radius: None,
            k_factor: None,
        }
    }

    /// Flat width of a 60 mm deep plate with this flange on its far edge.
    /// T 2, R 2, K 0.4.
    fn width(angle: f64, datum: HeightDatum, position: BendPosition) -> f64 {
        60.0 + section(2.0, 2.0, 0.4, &flange(angle, datum, position))
            .unwrap()
            .reach()
    }

    #[test]
    fn every_datum_and_position_gives_its_flat_width() {
        use BendPosition::*;
        use HeightDatum::*;
        for (angle, datum, position, expected) in [
            (90.0, Outer, Outside, 88.3982297150257),
            (90.0, Outer, Inside, 86.3982297150257),
            (90.0, Inner, Outside, 90.3982297150257),
            (90.0, Inner, Inside, 88.3982297150257),
            (90.0, Tangent, Outside, 88.3982297150257),
            (90.0, Tangent, Inside, 86.39822971502569),
            (45.0, Outer, Outside, 89.71383348327429),
            (135.0, Outer, Outside, 82.1120631983),
            (135.0, Tangent, Outside, 87.76891744779238),
            (135.0, OuterOrtho, Outside, 94.53847006949287),
            (135.0, InnerOrtho, Outside, 99.36689719423903),
        ] {
            let w = width(angle, datum, position);
            assert!(
                (w - expected).abs() < 1e-9,
                "{angle} {datum:?} {position:?}: {w} vs {expected}"
            );
        }
    }

    #[test]
    fn a_flange_too_short_for_its_bend_is_rejected() {
        let mut f = flange(90.0, HeightDatum::Outer, BendPosition::Outside);
        f.length = 3.0;
        assert!(section(2.0, 2.0, 0.4, &f).is_err());
    }
}
