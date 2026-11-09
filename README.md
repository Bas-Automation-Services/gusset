# Gusset

An open sheet metal kernel in Rust. Describe a part as JSON (a base plate and its flanges,
hems, contour flanges, folds, holes and cut-outs) and Gusset gives you the folded solid, the
flat pattern and its measurements, as STEP, DXF, SVG, a mesh or numbers.

Gusset is exact: the folded solid is an analytic B-rep, its mass properties are integrated, not
meshed. What it doesn't build yet it refuses with a clear message; it never approximates.

Made by [Bas Automation Services](https://www.basautomationservices.com).

## What it builds

- **Base plates** of lines and arcs (rounded corners), or a **profile** run along a width
  (a contour flange as the first feature: channels, hats, Z-sections).
- **Flanges** at any angle, every height datum, inside or outside, up or down, full or
  partial width with straight, round or tear reliefs; bend transitions (intersection, arc,
  trim to bend) where a flange ends next to a slanted edge; lips on their tips; rounded or
  chamfered tip corners.
- **Corners** between flanges at any angle: trim to bend with a miter gap, round, square,
  weld or tear reliefs, overlapping seams; flanges, hems and contours with different settings
  meeting by build order, or joined by a corner seam. Lips on both flanges at a corner are mitered there like one contour
  along both edges when flanges and lips are alike (one feature on two edges).
- **Hems:** single, double, teardrop, rolled, full or partial width.
- **Contour flanges** of lines and tangent arcs on a base edge or a flange's tip, full or
  partial width, or along several edges in one go, mitered at each corner (every leg and
  bend cut exactly, with corner reliefs).
- **Rolled cones and cylinders:** the sheet between two coaxial arcs (a lofted flange
  between them, or a straight line rolled around an axis), material on either side, unrolled
  on its neutral surface.
- **Folds** along a line across the plate or a flange, nested, at any angle to the outline; **jogs**
  (Z-steps by an offset).
- **Holes and cut-outs** anywhere in the flat, across bends, and across the outline as
  notches, also through a mitered bend's curved end; or placed on a flange's face, measured
  along its edge and up from the mold line.
- **Punches:** a shape of holes and cut-outs stamped at several points, turned by an angle;
  the flat pattern marks each stamp's centre and punch ID.
- **Bend order:** every bend numbered in the order it's made (outermost first, or as you
  list them), with a note at each bend line in the flat DXF and SVG.
- **Unfolding** by K-factor (per part or per feature) or by a bend table (DIN 6935 style: bend deduction by radius and angle).
- **Shop checks:** bend radius, leg length, hole distance to bends and edges, hole size.
- **Parts that fold into themselves** are refused, naming the features and where they meet;
  `clash::clashes` gives the overlapping areas on the flat pattern, `clash::preview` the
  folded sheet anyway with the clash marked.

## A part

```json
{
  "id": "bracket",
  "thickness": 2, "bend_radius": 2, "k_factor": 0.4,
  "base": {
    "outline": [[0, 0], [100, 0], [100, 60], [0, 60]],
    "holes": [{ "center": [50, 30], "diameter": 8 }]
  },
  "flanges": [{
    "id": "f", "on": { "base": 2 }, "angle": 90, "length": 30,
    "datum": "outer", "position": "outside", "side": "up"
  }]
}
```

Millimetres and degrees. The full description is `schema/part.schema.json`.

## Crates

| Crate | What |
|---|---|
| `gusset-sheet` | The kernel: `flat::flat_pattern`, `fold::fold` (an exact B-rep), `fold::exact_mass`, `check::check`, `clash::clashes` |
| `gusset-export` | Flat pattern as DXF and SVG (bend notes, punch marks); the folded part as STEP AP214 |
| `gusset-wasm` | All of it for the browser (wasm-bindgen): flat DXF/SVG, STEP, mesh, mass, checks, clashes |
| `gusset-capi` | A C library (`gusset.dll`): `include/gusset.h`, and `dotnet/Gusset.cs` for .NET 8 |

```rust
let part: PartSpec = serde_json::from_str(json)?;
let flat = gusset_sheet::flat::flat_pattern(&part)?;   // outline, holes, bend lines
let solid = gusset_sheet::fold::fold(&part)?;           // exact B-rep
let dxf = gusset_export::dxf::flat_dxf(&flat);
let step = gusset_export::step::part_step(&part)?;
```

```csharp
var (volume, area, centre) = Gusset.Mass(json);
File.WriteAllText("flat.dxf", Gusset.FlatDxf(json));
```

## Demo

**Try it: https://bas-automation-services.github.io/gusset/**

`demo/` is a static page that runs the kernel as WebAssembly: edit a part, see it folded
and flat, its mass, bend table and shop checks. A part that folds into itself shows the
clash in red, folded and flat.

    wasm-pack build crates/gusset-wasm --target web --release --out-dir ../../demo/pkg
    python -m http.server -d demo 8000

## Tests

`cargo test` checks flat patterns, folded solids, mass properties and exports of the main
cases against known values. `examples/measure.rs` prints flat size, flat area and folded
volume of any parts, for comparing with another CAD system; `examples/fuzz.rs` throws random
parts at the kernel.

The geometry kernel is a fork of [OpenGeometry](https://github.com/OpenGeometry-io/OpenGeometry)
(KERNEL.md: why, and what the fork changes).

## License

MPL-2.0.
