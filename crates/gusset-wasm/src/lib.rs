//! Gusset for the browser. Every function takes a part description as JSON
//! (schema: `schema()`) and fails with the kernel's message.

use gusset_export::{dxf::flat_dxf, step::part_step, svg::flat_svg_marked};
use gusset_sheet::flat::flat_pattern;
use gusset_sheet::fold::{exact_mass, fold};
use gusset_sheet::spec::PartSpec;
use opengeometry::tessellation::tessellate;
use wasm_bindgen::prelude::*;

fn part(json: &str) -> Result<PartSpec, JsError> {
    serde_json::from_str(json).map_err(|e| JsError::new(&format!("spec: {e}")))
}

fn kernel<T, E: std::fmt::Display>(r: Result<T, E>) -> Result<T, JsError> {
    r.map_err(|e| JsError::new(&e.to_string()))
}

/// The flat pattern as DXF (R12, millimetres).
#[wasm_bindgen(js_name = flatDxf)]
pub fn flat_dxf_js(json: &str) -> Result<String, JsError> {
    Ok(flat_dxf(&kernel(flat_pattern(&part(json)?))?))
}

/// The flat pattern as SVG, with any areas that fold into each other marked
/// (class `clash`).
#[wasm_bindgen(js_name = flatSvg)]
pub fn flat_svg_js(json: &str) -> Result<String, JsError> {
    let spec = part(json)?;
    let flat = kernel(flat_pattern(&spec))?;
    let marks: Vec<Vec<[f64; 2]>> = gusset_sheet::clash::clashes(&spec)
        .unwrap_or_default()
        .into_iter()
        .flat_map(|c| c.areas.into_iter().flatten())
        .collect();
    Ok(flat_svg_marked(&flat, &marks))
}

/// Where the part folds into itself, as JSON: `[{features, depth, areas,
/// message}]`, deepest first (areas: convex pieces on the flat pattern, per
/// feature); empty when it doesn't.
#[wasm_bindgen]
pub fn clashes(json: &str) -> Result<String, JsError> {
    let found = kernel(gusset_sheet::clash::clashes(&part(json)?))?;
    let list: Vec<serde_json::Value> = found
        .iter()
        .map(|c| {
            let mut v = serde_json::to_value(c).unwrap();
            v["message"] = c.to_string().into();
            v
        })
        .collect();
    Ok(serde_json::Value::from(list).to_string())
}

/// A part that folds into itself, shown anyway: the sheet's two faces and,
/// apart, the areas that clash, as loose triangles.
#[wasm_bindgen]
pub struct Preview {
    sheet: Mesh,
    clash: Mesh,
}

#[wasm_bindgen]
impl Preview {
    #[wasm_bindgen(getter)]
    pub fn sheet(&self) -> Mesh {
        self.sheet.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn clash(&self) -> Mesh {
        self.clash.clone()
    }
}

#[wasm_bindgen]
pub fn preview(json: &str) -> Result<Preview, JsError> {
    let p = kernel(gusset_sheet::clash::preview(&part(json)?))?;
    Ok(Preview {
        sheet: Mesh::loose(&p.sheet),
        clash: Mesh::loose(&p.clashes),
    })
}

/// The flat blank's size (mm, along x and y) and its bends (feature, order,
/// angle, inner radius, K-factor, up, length), as JSON.
#[wasm_bindgen]
pub fn flat(json: &str) -> Result<String, JsError> {
    let f = kernel(flat_pattern(&part(json)?))?;
    let [lo, hi] = f.bounds();
    let bends: Vec<String> = f
        .bends
        .iter()
        .map(|b| {
            let d = [b.line[1][0] - b.line[0][0], b.line[1][1] - b.line[0][1]];
            format!(
                r#"{{"feature":{},"order":{},"angle":{},"inner_radius":{},"k_factor":{},"up":{},"length":{}}}"#,
                serde_json::to_string(&b.flange).unwrap(),
                b.order,
                b.angle,
                b.inner_radius,
                b.k_factor,
                b.up,
                (d[0] * d[0] + d[1] * d[1]).sqrt()
            )
        })
        .collect();
    Ok(format!(
        r#"{{"width":{},"height":{},"bends":[{}]}}"#,
        hi[0] - lo[0],
        hi[1] - lo[1],
        bends.join(",")
    ))
}

/// Manufacturability findings as JSON: `[{rule, feature, value, limit,
/// message}]`. `rules` is JSON for `gusset_sheet::check::Rules` (multiples of
/// the thickness); empty for the defaults.
#[wasm_bindgen]
pub fn check(json: &str, rules: &str) -> Result<String, JsError> {
    let rules: gusset_sheet::check::Rules = if rules.trim().is_empty() {
        Default::default()
    } else {
        serde_json::from_str(rules).map_err(|e| JsError::new(&format!("rules: {e}")))?
    };
    let findings = kernel(gusset_sheet::check::check(&part(json)?, &rules))?;
    let list: Vec<serde_json::Value> = findings
        .iter()
        .map(|f| {
            let mut v = serde_json::to_value(f).unwrap();
            v["message"] = f.to_string().into();
            v
        })
        .collect();
    Ok(serde_json::Value::from(list).to_string())
}

/// The folded part as STEP (AP214, millimetres).
#[wasm_bindgen]
pub fn step(json: &str) -> Result<String, JsError> {
    kernel(part_step(&part(json)?))
}

/// Volume (mm3), surface area (mm2) and centre of mass (mm) as JSON.
#[wasm_bindgen]
pub fn mass(json: &str) -> Result<String, JsError> {
    let m = kernel(exact_mass(&part(json)?))?;
    Ok(format!(
        r#"{{"volume":{},"area":{},"center_of_mass":[{},{},{}]}}"#,
        m.volume, m.area, m.center_of_mass[0], m.center_of_mass[1], m.center_of_mass[2]
    ))
}

/// The folded part as triangles, within `deflection` mm of its surfaces.
#[wasm_bindgen]
#[derive(Clone)]
pub struct Mesh {
    positions: Vec<f32>,
    normals: Vec<f32>,
    indices: Vec<u32>,
}

#[wasm_bindgen]
impl Mesh {
    #[wasm_bindgen(getter)]
    pub fn positions(&self) -> Vec<f32> {
        self.positions.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn normals(&self) -> Vec<f32> {
        self.normals.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn indices(&self) -> Vec<u32> {
        self.indices.clone()
    }
}

impl Mesh {
    /// Unshared triangles, each with its own normal.
    fn loose(triangles: &[[[f64; 3]; 3]]) -> Mesh {
        let mut m = Mesh {
            positions: Vec::new(),
            normals: Vec::new(),
            indices: Vec::new(),
        };
        for t in triangles {
            let (u, v) = (
                std::array::from_fn::<f64, 3, _>(|j| t[1][j] - t[0][j]),
                std::array::from_fn::<f64, 3, _>(|j| t[2][j] - t[0][j]),
            );
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let l = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-300);
            for p in t {
                m.indices.push(m.indices.len() as u32);
                m.positions.extend(p.iter().map(|&x| x as f32));
                m.normals.extend(n.iter().map(|&x| (x / l) as f32));
            }
        }
        m
    }
}

#[wasm_bindgen]
pub fn mesh(json: &str, deflection: f64) -> Result<Mesh, JsError> {
    let brep = kernel(fold(&part(json)?))?;
    let t = kernel(tessellate(&brep, deflection, 2_000_000))?;
    Ok(Mesh {
        positions: t.positions.iter().map(|&x| x as f32).collect(),
        normals: t.normals,
        indices: t.indices,
    })
}

/// The part description's JSON Schema.
#[wasm_bindgen]
pub fn schema() -> String {
    include_str!("../../../schema/part.schema.json").into()
}

#[cfg(test)]
mod tests {
    use super::*;

    const PLATE: &str = r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
        "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
        "flanges": [{ "id": "f", "on": { "base": 2 }, "angle": 90, "length": 30,
                      "datum": "outer", "position": "outside", "side": "up" }] }"#;

    #[test]
    fn a_flanged_plate_comes_out_in_every_format() {
        assert!(flat_dxf_js(PLATE).unwrap().contains("BendUp"));
        assert!(flat_svg_js(PLATE).unwrap().contains("bend up"));
        let f: serde_json::Value = serde_json::from_str(&flat(PLATE).unwrap()).unwrap();
        assert_eq!(f["bends"][0]["angle"].as_f64(), Some(90.0));
        assert!(step(PLATE).unwrap().contains("ISO-10303-21"));
        let m: serde_json::Value = serde_json::from_str(&mass(PLATE).unwrap()).unwrap();
        assert!((m["volume"].as_f64().unwrap() - 17742.4778).abs() < 1e-3);
        let m = mesh(PLATE, 0.05).unwrap();
        assert_eq!(m.positions.len(), m.normals.len());
        assert_eq!(m.indices.len() % 3, 0);
        assert!(schema().contains("PartSpec"));
        assert_eq!(check(PLATE, "").unwrap(), "[]");
        assert_eq!(clashes(PLATE).unwrap(), "[]");
    }

    #[test]
    fn a_cone_shows_its_unrolled_blank() {
        let part = r#"{ "id": "c", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
            "base": { "cone": { "radii": [50, 30], "height": 60, "sweep": 180, "material": "outside" } } }"#;
        let f: serde_json::Value = serde_json::from_str(&flat(part).unwrap()).unwrap();
        // The neutral sector about the apex, 134.4976 x 107.4253.
        let l = 20f64.hypot(60.0);
        let outer = (50.0 + 0.8 * 60.0 / l) * l / 20.0;
        let (inner, angle) = (outer - l, std::f64::consts::PI * 20.0 / l);
        assert!((f["width"].as_f64().unwrap() - outer * angle.sin()).abs() < 1e-9);
        assert!((f["height"].as_f64().unwrap() - (outer - inner * angle.cos())).abs() < 1e-9);
        assert_eq!(f["bends"].as_array().unwrap().len(), 0);
        assert!(!mesh(part, 0.05).unwrap().indices.is_empty());
    }

    #[test]
    fn a_part_that_folds_into_itself_shows_where() {
        // Two flaps folded 120 degrees over a 10 wide strip cross.
        let part = r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
            "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
            "folds": [{ "id": "a", "line": [[0,25],[100,25]], "angle": 120, "location": "centerline", "side": "up" },
                      { "id": "b", "line": [[100,35],[0,35]], "angle": 120, "location": "centerline", "side": "up" }] }"#;
        let c: serde_json::Value = serde_json::from_str(&clashes(part).unwrap()).unwrap();
        assert_eq!(c[0]["features"], serde_json::json!(["a", "b"]));
        assert!(c[0]["message"]
            .as_str()
            .unwrap()
            .starts_with("a and b fold into each other"));
        assert!(flat_svg_js(part).unwrap().contains(r#"class="clash""#));
        let p = preview(part).unwrap();
        assert!(!p.sheet.indices.is_empty() && !p.clash.indices.is_empty());
        assert_eq!(p.clash.positions.len(), p.clash.normals.len());
    }
}
