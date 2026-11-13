//! What the tools do, apart from the protocol: a part description in, text an
//! agent can act on out. Every failure is the kernel's message, verbatim.

use gusset_export::{
    dxf::flat_dxf,
    step::part_step,
    svg::{flat_svg, flat_svg_marked},
};
use gusset_sheet::flat::{flat_pattern, FlatPattern};
use gusset_sheet::fold::exact_mass;
use gusset_sheet::spec::PartSpec;
use gusset_sheet::{check, clash};
use serde_json::{json, Value};

pub const SCHEMA: &str = include_str!("../../../schema/part.schema.json");

/// Where `build` writes its files, if anywhere.
#[derive(Default)]
pub struct Outputs {
    pub step: Option<String>,
    pub dxf: Option<String>,
    pub svg: Option<String>,
}

/// A part as a JSON object, or as JSON text.
pub fn parse(part: &Value) -> Result<PartSpec, String> {
    match part {
        Value::String(text) => serde_json::from_str(text),
        v => serde_json::from_value(v.clone()),
    }
    .map_err(|e| format!("spec: {e}"))
}

/// Runs the kernel, turning a panic into an error instead of a request that
/// never answers.
pub fn guarded<T>(f: impl FnOnce() -> Result<T, String>) -> Result<T, String> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f))
        .unwrap_or_else(|_| Err("kernel: Gusset panicked on this part (a bug in Gusset)".into()))
}

fn round(v: f64) -> f64 {
    (v * 1e4).round() / 1e4
}

/// Builds the part and sums it up: flat pattern, folded mass, bend table and
/// shop checks; writes the files `out` asks for.
pub fn build(spec: &PartSpec, out: &Outputs) -> Result<Value, String> {
    let flat = flat_pattern(spec).map_err(|e| e.to_string())?;
    let mass = exact_mass(spec).map_err(|e| with_clashes(spec, e.to_string()))?;
    let findings = check::check(spec, &Default::default()).map_err(|e| e.to_string())?;

    let mut written = Vec::new();
    let mut write = |path: &String, text: String| {
        std::fs::write(path, text).map_err(|e| format!("writing {path}: {e}"))?;
        written.push(path.clone());
        Ok::<_, String>(())
    };
    if let Some(p) = &out.step {
        write(p, part_step(spec).map_err(|e| e.to_string())?)?;
    }
    if let Some(p) = &out.dxf {
        write(p, flat_dxf(&flat))?;
    }
    if let Some(p) = &out.svg {
        write(p, flat_svg(&flat))?;
    }

    let [lo, hi] = flat.bounds();
    Ok(json!({
        "id": spec.id,
        "units": "mm, mm2, mm3, degrees",
        "flat": {
            "width": round(hi[0] - lo[0]),
            "height": round(hi[1] - lo[1]),
            "min": lo.map(round),
            "max": hi.map(round),
            "area": round(flat.area()),
            "holes": flat.holes.len(),
            "cutouts": flat.cutouts.len(),
            "punches": flat.punches.len(),
        },
        "folded": {
            "volume": round(mass.volume),
            "area": round(mass.area),
            "center_of_mass": mass.center_of_mass.map(round),
        },
        "bends": bends(&flat),
        "checks": findings.iter().map(|f| {
            let mut v = serde_json::to_value(f).unwrap();
            v["message"] = f.to_string().into();
            v
        }).collect::<Vec<_>>(),
        "written": written,
    }))
}

/// The bend table, in the flat pattern's order.
fn bends(flat: &FlatPattern) -> Vec<Value> {
    flat.bends
        .iter()
        .map(|b| {
            let [p, q] = b.line;
            json!({
                "feature": b.flange,
                "order": b.order,
                "direction": if b.up { "up" } else { "down" },
                "angle": round(b.angle),
                "inner_radius": round(b.inner_radius),
                "k_factor": b.k_factor,
                "length": round((q[0] - p[0]).hypot(q[1] - p[1])),
                "note": b.note(),
            })
        })
        .collect()
}

/// A failed fold's message, with every place the part folds into itself that
/// it doesn't already name.
fn with_clashes(spec: &PartSpec, message: String) -> String {
    let more: Vec<String> = clash::clashes(spec)
        .unwrap_or_default()
        .iter()
        .map(|c| c.to_string())
        .filter(|c| !message.contains(c.as_str()))
        .collect();
    if more.is_empty() {
        return message;
    }
    format!(
        "{message}\nAlso folds into itself:\n- {}",
        more.join("\n- ")
    )
}

/// The flat pattern as SVG, with any areas that fold into each other marked
/// (class `clash`).
pub fn svg(spec: &PartSpec) -> Result<String, String> {
    let flat = flat_pattern(spec).map_err(|e| e.to_string())?;
    let marks: Vec<Vec<[f64; 2]>> = clash::clashes(spec)
        .unwrap_or_default()
        .into_iter()
        .flat_map(|c| c.areas.into_iter().flatten())
        .collect();
    Ok(flat_svg_marked(&flat, &marks))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bracket() -> PartSpec {
        parse(&Value::String(bracket!().into())).unwrap()
    }

    #[test]
    fn the_bracket_sums_up() {
        let s = build(&bracket(), &Outputs::default()).unwrap();
        // The flanged plate's 17742.4778 less the 8 mm hole through 2 mm.
        let hole = std::f64::consts::PI * 16.0 * 2.0;
        assert!((s["folded"]["volume"].as_f64().unwrap() - (17742.4778 - hole)).abs() < 1e-3);
        assert_eq!(s["flat"]["width"].as_f64(), Some(100.0));
        assert_eq!(s["flat"]["holes"], 1);
        assert_eq!(s["bends"].as_array().unwrap().len(), 1);
        let b = &s["bends"][0];
        assert_eq!(
            (b["feature"].as_str(), b["order"].as_u64()),
            (Some("f"), Some(1))
        );
        assert_eq!(
            (b["angle"].as_f64(), b["direction"].as_str()),
            (Some(90.0), Some("up"))
        );
        assert_eq!(b["inner_radius"].as_f64(), Some(2.0));
        assert_eq!(s["checks"], json!([]));
        assert_eq!(s["written"], json!([]));
    }

    #[test]
    fn an_object_and_its_text_are_the_same_part() {
        let v: Value = serde_json::from_str(bracket!()).unwrap();
        assert_eq!(parse(&v).unwrap(), bracket());
    }

    #[test]
    fn the_bracket_writes_its_files() {
        let dir = std::env::temp_dir().join(format!("gusset-mcp-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = |name: &str| Some(dir.join(name).to_string_lossy().into_owned());
        let out = Outputs {
            step: path("b.step"),
            dxf: path("b.dxf"),
            svg: path("b.svg"),
        };
        let s = build(&bracket(), &out).unwrap();
        assert_eq!(s["written"].as_array().unwrap().len(), 3);
        let read = |name: &str| std::fs::read_to_string(dir.join(name)).unwrap();
        assert!(read("b.step").contains("ISO-10303-21"));
        assert!(read("b.dxf").contains("BendUp"));
        assert!(read("b.svg").starts_with("<svg"));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_bad_part_fails_with_the_kernels_message() {
        assert!(parse(&json!({ "id": "x" }))
            .unwrap_err()
            .starts_with("spec:"));
        assert!(parse(&json!("{")).unwrap_err().starts_with("spec:"));
        let bad = parse(&json!(bracket!().replace(r#""base": 2"#, r#""base": 9"#))).unwrap();
        assert!(build(&bad, &Outputs::default())
            .unwrap_err()
            .contains("no base edge 9"));
    }

    #[test]
    fn a_part_that_folds_into_itself_says_where() {
        // Two flaps folded 120 degrees over a 10 wide strip cross.
        let part = parse(&json!({ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
            "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
            "folds": [{ "id": "a", "line": [[0,25],[100,25]], "angle": 120, "location": "centerline", "side": "up" },
                      { "id": "b", "line": [[100,35],[0,35]], "angle": 120, "location": "centerline", "side": "up" }] }))
        .unwrap();
        let e = build(&part, &Outputs::default()).unwrap_err();
        assert!(e.contains("a and b fold into each other"), "{e}");
        assert!(svg(&part).unwrap().contains(r#"class="clash""#));
    }

    #[test]
    fn a_panic_comes_back_as_an_error() {
        let r: Result<(), String> = guarded(|| panic!("boom"));
        assert!(r.unwrap_err().starts_with("kernel:"));
    }
}
