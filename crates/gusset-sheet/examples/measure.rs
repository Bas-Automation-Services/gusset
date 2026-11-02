//! Flat size, flat area and folded volume of parts, for comparing with another CAD system.
//!
//!     cargo run --release --example measure -- parts.json
//!     OUTLINE=1 ...   (also print the flat outline's corners)
//!
//! `parts.json` holds an array of parts. The flat size is the extent of the outline's corner
//! points (arcs bulging past them aren't counted).
use gusset_sheet::{flat, fold, spec::PartSpec};

fn main() {
    let path = std::env::args().nth(1).expect("parts.json");
    let parts: Vec<PartSpec> =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    for part in &parts {
        match (flat::flat_pattern(part), fold::exact_mass(part)) {
            (Ok(f), Ok(m)) => {
                let extent = |i: usize| {
                    let v = f.outline.iter().map(|p| p[i]);
                    v.clone().fold(f64::MIN, f64::max) - v.fold(f64::MAX, f64::min)
                };
                let holes: f64 = f
                    .holes
                    .iter()
                    .map(|h| std::f64::consts::PI * h.diameter * h.diameter / 4.0)
                    .sum();
                let area = flat::path_area(&f.outline, &f.outline_arcs).abs() - holes;
                println!(
                    "{}\tflat {:.6} x {:.6}\tarea {:.6}\tvolume {:.6}",
                    part.id,
                    extent(0),
                    extent(1),
                    area,
                    m.volume
                );
                if std::env::var("OUTLINE").is_ok() {
                    for (p, c) in f.outline.iter().zip(&f.outline_arcs) {
                        let curve = if c.is_some() { " curve" } else { "" };
                        println!("    {:.4}, {:.4}{curve}", p[0], p[1]);
                    }
                }
            }
            (Err(e), _) | (_, Err(e)) => println!("{}\t{e}", part.id),
        }
    }
}
