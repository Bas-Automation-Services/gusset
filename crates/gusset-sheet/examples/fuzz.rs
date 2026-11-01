//! Random parts through the flat pattern, the fold, the exact mass and a mesh.
//! Invalid and Unsupported are fine answers; a panic, a kernel error or a mesh
//! volume off the exact one is reported with the part's JSON.
//!
//!     cargo run --release --example fuzz -- <count> <seed>
//!     SAFE=1 cargo run --release --example fuzz -- 1000 7   (fewer self-intersecting parts)
//!     WHY=1 ...                                             (print each refusal and its part)
//!
//! Without SAFE many parts fold into themselves (flanges past 90 degrees into a
//! neighbour); Gusset refuses those (`clash`).
use std::panic;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next()
    }
    fn pick<'a>(&mut self, xs: &[&'a str]) -> &'a str {
        xs[(self.next() * xs.len() as f64) as usize % xs.len()]
    }
    fn chance(&mut self, p: f64) -> bool {
        self.next() < p
    }
}

fn safe() -> bool {
    std::env::var("SAFE").is_ok()
}

fn spec(r: &mut Rng) -> String {
    let t = [1.0, 1.5, 2.0, 3.0][(r.next() * 4.0) as usize % 4];
    let rad = r.range(0.5, 5.0);
    let (w, h) = (r.range(40.0, 150.0), r.range(30.0, 100.0));
    let n = 4;
    let outline = format!("[[0,0],[{w},0],[{w},{h}],[0,{h}]]");
    let mut flanges = Vec::new();
    let mut hems = Vec::new();
    let mut contours = Vec::new();
    for e in 0..n {
        if r.chance(0.45) {
            let width = if r.chance(0.3) {
                let len = if e % 2 == 0 { w } else { h };
                let a = r.range(0.0, len * 0.3);
                let b = r.range(0.0, len * 0.3);
                format!(r#""width": {{ "offset_start": {a}, "offset_end": {b} }},"#)
            } else {
                String::new()
            };
            let tips = if r.chance(0.15) {
                format!(
                    r#", "tip_corners": {{ "round": {{ "radius": {} }} }}"#,
                    r.range(1.0, 6.0)
                )
            } else {
                String::new()
            };
            flanges.push(format!(
                r#"{{ "id": "f{e}", "on": {{ "base": {e} }}, {width} "angle": {}, "length": {},
                     "datum": "{}", "position": "{}", "side": "{}" {tips} }}"#,
                if safe() {
                    r.range(30.0, 120.0)
                } else {
                    r.range(20.0, 160.0)
                },
                r.range(8.0, 40.0),
                r.pick(&["outer", "inner", "tangent", "outer_ortho", "inner_ortho"]),
                r.pick(&["inside", "outside"]),
                r.pick(&["up", "down"]),
            ));
            if r.chance(0.25) {
                flanges.push(format!(
                    r#"{{ "id": "l{e}", "on": {{ "tip": "f{e}" }}, "angle": {}, "length": {},
                         "datum": "outer", "position": "outside", "side": "{}" }}"#,
                    if safe() {
                        r.range(30.0, 60.0)
                    } else {
                        r.range(30.0, 150.0)
                    },
                    r.range(6.0, 20.0),
                    if safe() {
                        "down"
                    } else {
                        r.pick(&["up", "down"])
                    },
                ));
            }
        } else if r.chance(0.2) {
            hems.push(format!(
                r#"{{ "id": "h{e}", "on": {{ "base": {e} }}, "kind": {{ "single": {{ "gap": {}, "length": {} }} }}, "side": "{}" }}"#,
                r.range(0.5, 3.0),
                r.range(6.0, 15.0),
                r.pick(&["up", "down"]),
            ));
        } else if r.chance(0.2) {
            contours.push(format!(
                r#"{{ "id": "c{e}", "on": {{ "base": {e} }}, "profile": [[0,{t}],[0,{}],[{},{}]], "face": "bottom" }}"#,
                r.range(10.0, 40.0),
                if safe() { r.range(5.0, 20.0) } else { r.range(-20.0, 20.0) },
                r.range(10.0, 40.0),
            ));
        }
    }
    let mut holes = Vec::new();
    for _ in 0..(r.next() * 4.0) as usize {
        holes.push(format!(
            r#"{{ "center": [{}, {}], "diameter": {} }}"#,
            r.range(-10.0, w + 10.0),
            r.range(-10.0, h + 30.0),
            r.range(2.0, 12.0)
        ));
    }
    let relief = match (r.next() * 6.0) as usize {
        0 => r#""relief": { "width": 2, "depth": 1, "minimum_remnant": 4, "corner": { "shape": "round", "size": 12 } },"#.to_string(),
        // Small ones too: flat only, capped or refused depending on the bends.
        5 => format!(
            r#""relief": {{ "width": 2, "depth": 1, "minimum_remnant": 4, "corner": {{ "shape": "{}", "size": {} }} }},"#,
            r.pick(&["round", "square"]),
            r.range(1.0, 9.0)
        ),
        1 => r#""relief": { "width": 2, "depth": 1, "minimum_remnant": 4, "corner": { "shape": "linear_weld" } },"#.to_string(),
        2 => r#""relief": { "width": 2, "depth": 1, "minimum_remnant": 4, "seam": "overlap" },"#.to_string(),
        _ => String::new(),
    };
    format!(
        r#"{{ "id": "fuzz", "thickness": {t}, "bend_radius": {rad}, "k_factor": {}, {relief}
             "base": {{ "outline": {outline}, "holes": [{}] }},
             "flanges": [{}], "hems": [{}], "contours": [{}] }}"#,
        r.range(0.25, 0.5),
        holes.join(","),
        flanges.join(","),
        hems.join(","),
        contours.join(",")
    )
}

/// A profile part: a channel, hat or Z with random sizes, sometimes holes,
/// end flanges, hems on its ends or a fold across its base.
fn profile_part(r: &mut Rng) -> String {
    let t = [1.0, 1.5, 2.0, 3.0][(r.next() * 4.0) as usize % 4];
    let (w, h, b) = (
        r.range(40.0, 150.0),
        r.range(15.0, 50.0),
        r.range(30.0, 90.0),
    );
    let (points, material) = match (r.next() * 3.0) as usize {
        0 => (format!("[[0,{h}],[0,0],[{b},0],[{b},{h}]]"), "left"),
        1 => (
            format!(
                "[[0,0],[20,0],[20,{h}],[{},{h}],[{},0],[{},0]]",
                20.0 + b,
                20.0 + b,
                40.0 + b
            ),
            "left",
        ),
        _ => (
            format!("[[0,{h}],[0,{t}],[{b},{t}],[{b},{}]]", t - h),
            "right",
        ),
    };
    let mut features = Vec::new();
    if r.chance(0.3) {
        features.push(format!(
            r#""hems": [{{ "id": "h", "on": {{ "tip": "end" }}, "kind": {{ "single": {{ "gap": {}, "length": {} }} }}, "side": "{}" }}]"#,
            r.range(0.5, 3.0),
            r.range(6.0, 12.0),
            r.pick(&["up", "down"])
        ));
    }
    if r.chance(0.3) {
        features.push(format!(
            r#""flanges": [{{ "id": "e", "on": {{ "base": {} }}, "angle": 90, "length": {},
                             "datum": "outer", "position": "outside", "side": "up" }}]"#,
            r.pick(&["1", "3"]),
            r.range(8.0, h.min(30.0))
        ));
    }
    if r.chance(0.3) {
        let y = r.range(5.0, b - 5.0);
        features.push(format!(
            r#""base": {{ "profile": {{ "points": {points}, "material": "{material}", "width": {w} }},
                         "holes": [{{ "center": [{}, {y}], "diameter": {} }}] }}"#,
            r.range(5.0, w - 5.0),
            r.range(2.0, 8.0)
        ));
    } else {
        features.push(format!(
            r#""base": {{ "profile": {{ "points": {points}, "material": "{material}", "width": {w} }} }}"#
        ));
    }
    format!(
        r#"{{ "id": "fuzz", "thickness": {t}, "bend_radius": {}, "k_factor": {}, {} }}"#,
        r.range(0.5, 5.0),
        r.range(0.25, 0.5),
        features.join(",")
    )
}

/// A plate with one or two fold lines, maybe flanges and holes.
fn folded_part(r: &mut Rng) -> String {
    let t = [1.0, 1.5, 2.0, 3.0][(r.next() * 4.0) as usize % 4];
    let (w, h) = (r.range(60.0, 150.0), r.range(60.0, 120.0));
    let y1 = r.range(15.0, h / 2.0 - 5.0);
    let mut folds = vec![format!(
        r#"{{ "id": "a", "line": [[{w},{y1}],[0,{y1}]], "angle": {}, "location": "{}", "side": "{}" }}"#,
        r.range(30.0, 120.0),
        r.pick(&["centerline", "start", "end"]),
        r.pick(&["up", "down"])
    )];
    if r.chance(0.5) {
        let y2 = r.range(h / 2.0 + 5.0, h - 15.0);
        folds.push(format!(
            r#"{{ "id": "b", "line": [[0,{y2}],[{w},{y2}]], "angle": {}, "location": "centerline", "side": "{}" }}"#,
            r.range(30.0, 120.0),
            r.pick(&["up", "down"])
        ));
    }
    let flange = if r.chance(0.4) {
        format!(
            r#""flanges": [{{ "id": "f", "on": {{ "base": {} }}, "angle": 90, "length": {}, "datum": "outer", "position": "outside", "side": "{}" }}],"#,
            r.pick(&["0", "2"]),
            r.range(8.0, 25.0),
            r.pick(&["up", "down"])
        )
    } else {
        String::new()
    };
    let holes: Vec<String> = (0..(r.next() * 3.0) as usize)
        .map(|_| {
            format!(
                r#"{{ "center": [{}, {}], "diameter": {} }}"#,
                r.range(5.0, w - 5.0),
                r.range(5.0, h - 5.0),
                r.range(2.0, 8.0)
            )
        })
        .collect();
    format!(
        r#"{{ "id": "fuzz", "thickness": {t}, "bend_radius": {}, "k_factor": {}, {flange}
             "base": {{ "outline": [[0,0],[{w},0],[{w},{h}],[0,{h}]], "holes": [{}] }},
             "folds": [{}] }}"#,
        r.range(0.5, 5.0),
        r.range(0.25, 0.5),
        holes.join(","),
        folds.join(",")
    )
}

fn main() {
    let count: usize = std::env::args().nth(1).map_or(500, |s| s.parse().unwrap());
    let seed: u64 = std::env::args().nth(2).map_or(1, |s| s.parse().unwrap());
    panic::set_hook(Box::new(|_| {}));
    let mut rng = Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1);
    let (mut ok, mut refused, mut bad) = (0, 0, 0);
    for i in 0..count {
        let json = match std::env::var("KIND").as_deref() {
            Ok("profile") => profile_part(&mut rng),
            Ok("fold") => folded_part(&mut rng),
            _ => spec(&mut rng),
        };
        let part: gusset_sheet::spec::PartSpec = serde_json::from_str(&json).unwrap();
        let result = panic::catch_unwind(|| {
            let flat = gusset_sheet::flat::flat_pattern(&part)?;
            let brep = gusset_sheet::fold::fold(&part)?;
            let exact = gusset_sheet::fold::exact_mass(&part)?;
            let mesh = gusset_sheet::fold::mass(&brep, 0.02)?;
            let rel = (exact.volume - mesh.volume).abs() / exact.volume;
            if rel > 3e-3 {
                return Err(gusset_sheet::SheetError::Kernel(format!(
                    "volume: exact {} mesh {}",
                    exact.volume, mesh.volume
                )));
            }
            let _ = flat;
            Ok(())
        });
        match result {
            Ok(Ok(())) => ok += 1,
            Ok(Err(gusset_sheet::SheetError::Kernel(e))) => {
                bad += 1;
                println!("#{i} KERNEL {e}\n{json}\n");
            }
            Ok(Err(e)) => {
                refused += 1;
                if std::env::var("WHY").is_ok() {
                    println!(
                        "#{i} REFUSED {e}
{json}
"
                    );
                }
            }
            Err(p) => {
                bad += 1;
                let msg = p
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default();
                println!("#{i} PANIC {msg}\n{json}\n");
            }
        }
    }
    println!("{count} parts: {ok} built, {refused} refused, {bad} bad");
}
