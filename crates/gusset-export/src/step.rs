//! The folded part as STEP (AP214, millimetres) through the kernel's exporter.

use gusset_sheet::fold::fold;
use gusset_sheet::spec::PartSpec;
use gusset_sheet::SheetError;
use opengeometry::exchange::export_step_in;

pub fn part_step(spec: &PartSpec) -> Result<String, SheetError> {
    let brep = fold(spec)?;
    // Gusset models in millimetres.
    let (text, _report) =
        export_step_in(&brep, "millimetre", "millimetre").map_err(SheetError::from)?;
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flanged_plate_with_a_hole_exports_one_solid_with_exact_faces() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "outline": [[0,0],[100,0],[100,60],[0,60]],
                           "holes": [{ "center": [20, 30], "diameter": 8 }] },
                 "flanges": [{ "id": "f1", "on": { "base": 2 }, "angle": 90, "length": 30,
                               "datum": "outer", "position": "outside", "side": "up" }] }"#,
        )
        .unwrap();
        let step = part_step(&spec).unwrap();
        assert_eq!(step.matches("MANIFOLD_SOLID_BREP").count(), 1);
        // Two bend faces and the hole wall are exact cylinders.
        assert_eq!(step.matches("CYLINDRICAL_SURFACE").count(), 3);
        assert!(step.contains("'f1.bend.bottom'"));
        // Millimetres, unscaled: the plate is 100 wide, not 100000.
        assert!(step.contains("SI_UNIT(.MILLI.,.METRE.)"));
        assert!(step.contains("1.00000000000000000E2") && !step.contains("1.00000000000000000E5"));
    }

    #[test]
    fn a_cone_exports_with_conical_faces() {
        let spec: PartSpec = serde_json::from_str(
            r#"{ "id": "c", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
                 "base": { "cone": { "radii": [50, 30], "height": 60, "sweep": 180,
                                     "material": "outside" } } }"#,
        )
        .unwrap();
        let step = part_step(&spec).unwrap();
        assert_eq!(step.matches("MANIFOLD_SOLID_BREP").count(), 1);
        // Both faces of the sheet and both ends.
        assert_eq!(step.matches("CONICAL_SURFACE").count(), 4);
        assert_eq!(step.matches("PLANE(").count(), 2);
    }
}
