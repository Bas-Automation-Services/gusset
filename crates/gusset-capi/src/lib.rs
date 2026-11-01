//! Gusset as a C library. Every function takes a part description as UTF-8
//! JSON (schema/part.schema.json) and returns a new UTF-8 JSON string,
//! `{"ok": ...}` or `{"error": "..."}`, that the caller frees with
//! `gusset_free`. Header: include/gusset.h.

// The pointers come from C callers, who own the contract the header states.
#![allow(clippy::not_unsafe_ptr_arg_deref)]

use std::ffi::{c_char, CStr, CString};

use gusset_export::{dxf::flat_dxf, step::part_step, svg::flat_svg};
use gusset_sheet::flat::flat_pattern;
use gusset_sheet::fold::exact_mass;
use gusset_sheet::spec::PartSpec;
use serde_json::{json, Value};

fn call(spec: *const c_char, f: impl FnOnce(&PartSpec) -> Result<Value, String>) -> *mut c_char {
    let result = (|| {
        if spec.is_null() {
            return Err("no part description".to_string());
        }
        let text = unsafe { CStr::from_ptr(spec) }
            .to_str()
            .map_err(|e| format!("not UTF-8: {e}"))?;
        let part: PartSpec = serde_json::from_str(text).map_err(|e| format!("spec: {e}"))?;
        f(&part)
    })();
    let out = match result {
        Ok(v) => json!({ "ok": v }),
        Err(e) => json!({ "error": e }),
    };
    // serde_json escapes control characters, so the text has no NUL.
    CString::new(out.to_string()).unwrap().into_raw()
}

/// The flat pattern as DXF (R12, millimetres): `{"ok": "<dxf>"}`.
#[no_mangle]
pub extern "C" fn gusset_flat_dxf(spec: *const c_char) -> *mut c_char {
    call(spec, |p| {
        flat_pattern(p)
            .map(|f| flat_dxf(&f).into())
            .map_err(|e| e.to_string())
    })
}

/// The flat pattern as SVG: `{"ok": "<svg>"}`.
#[no_mangle]
pub extern "C" fn gusset_flat_svg(spec: *const c_char) -> *mut c_char {
    call(spec, |p| {
        flat_pattern(p)
            .map(|f| flat_svg(&f).into())
            .map_err(|e| e.to_string())
    })
}

/// The folded part as STEP (AP214, millimetres): `{"ok": "<step>"}`.
#[no_mangle]
pub extern "C" fn gusset_step(spec: *const c_char) -> *mut c_char {
    call(spec, |p| {
        part_step(p).map(Value::from).map_err(|e| e.to_string())
    })
}

/// Volume (mm3), area (mm2) and centre of mass (mm) of the folded part:
/// `{"ok": {"volume": .., "area": .., "center_of_mass": [x, y, z]}}`.
#[no_mangle]
pub extern "C" fn gusset_mass(spec: *const c_char) -> *mut c_char {
    call(spec, |p| {
        exact_mass(p)
            .map(|m| json!({ "volume": m.volume, "area": m.area, "center_of_mass": m.center_of_mass }))
            .map_err(|e| e.to_string())
    })
}

/// Manufacturability findings with the default rules:
/// `{"ok": [{"rule", "feature", "value", "limit", "message"}]}`.
#[no_mangle]
pub extern "C" fn gusset_check(spec: *const c_char) -> *mut c_char {
    call(spec, |p| {
        gusset_sheet::check::check(p, &Default::default())
            .map(|findings| {
                findings
                    .iter()
                    .map(|f| {
                        let mut v = serde_json::to_value(f).unwrap();
                        v["message"] = f.to_string().into();
                        v
                    })
                    .collect()
            })
            .map_err(|e| e.to_string())
    })
}

/// Frees a string returned by any gusset_ function. Null is ignored.
#[no_mangle]
pub extern "C" fn gusset_free(text: *mut c_char) {
    if !text.is_null() {
        drop(unsafe { CString::from_raw(text) });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(f: extern "C" fn(*const c_char) -> *mut c_char, spec: &str) -> Value {
        let input = CString::new(spec).unwrap();
        let out = f(input.as_ptr());
        let value = serde_json::from_str(unsafe { CStr::from_ptr(out) }.to_str().unwrap()).unwrap();
        gusset_free(out);
        value
    }

    const PLATE: &str = r#"{ "id": "t", "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
        "base": { "outline": [[0,0],[100,0],[100,60],[0,60]] },
        "flanges": [{ "id": "f", "on": { "base": 2 }, "angle": 90, "length": 30,
                      "datum": "outer", "position": "outside", "side": "up" }] }"#;

    #[test]
    fn a_flanged_plate_comes_back_wrapped_in_ok() {
        let v = (run(gusset_mass, PLATE)["ok"]["volume"].as_f64().unwrap() - 17742.4778).abs();
        assert!(v < 1e-3);
        assert!(run(gusset_flat_dxf, PLATE)["ok"]
            .as_str()
            .unwrap()
            .contains("BendUp"));
        assert!(run(gusset_flat_svg, PLATE)["ok"]
            .as_str()
            .unwrap()
            .starts_with("<svg"));
        assert!(run(gusset_step, PLATE)["ok"]
            .as_str()
            .unwrap()
            .contains("ISO-10303-21"));
        assert_eq!(run(gusset_check, PLATE)["ok"], json!([]));
    }

    #[test]
    fn errors_come_back_as_messages() {
        assert!(run(gusset_mass, "{")["error"]
            .as_str()
            .unwrap()
            .starts_with("spec:"));
        let bad = PLATE.replace(r#""base": 2"#, r#""base": 9"#);
        assert!(run(gusset_mass, &bad)["error"]
            .as_str()
            .unwrap()
            .contains("no base edge 9"));
        let out = gusset_mass(std::ptr::null());
        assert!(unsafe { CStr::from_ptr(out) }
            .to_str()
            .unwrap()
            .contains("no part description"));
        gusset_free(out);
    }
}
