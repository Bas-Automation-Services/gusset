# Gusset

An open sheet metal kernel in Rust. Describe a part as JSON (a base plate and its flanges,
hems, contour flanges, folds, faces joined by bends, holes and cut-outs) and Gusset gives you the folded solid, the
flat pattern and its measurements, as STEP, DXF, SVG, a mesh or numbers.

Gusset is exact: the folded solid is an analytic B-rep, its mass properties are integrated, not
meshed. What it doesn't build yet it refuses with a clear message; it never approximates.

Made by [Bas Automation Services](https://www.basautomationservices.com).

## What it builds

- **Base plates** of lines and arcs (rounded corners), or a **profile** run along a width
  (a contour flange as the first feature: channels, hats, Z-sections).
- **Flanges** at any angle, every height datum, up or down, positioned inside or outside, the
  bend outside the edge, the outer virtual sharp on it or the bend tangent to its plane, full or
  partial width with straight, round or tear reliefs; bend transitions (intersection, arc,
  trim to bend) where a flange ends next to a slanted edge; lips on their tips; rounded or
  chamfered tip corners.
- **Corners** between flanges at any angle: trim to bend with a miter gap, round, square,
  weld or tear reliefs, obround reliefs (a slot on the bisector; the bends run on to their
  flanges' ends), constant width (the bends end on the miter, the gap apart, from their
  starts on), overlapping or butt seams; flanges, hems and contours with different settings
  meeting by build order, or joined by a corner seam. Lips on both flanges at a corner are mitered there like one contour
  along both edges when flanges and lips are alike (one feature on two edges).
- **Hems:** single, double, teardrop, rolled, full or partial width.
- **Contour flanges** of lines and tangent arcs on a base edge or a flange's tip, full or
  partial width, or along several edges in one go, mitered at each corner (every leg and
  bend cut exactly, with corner reliefs).
- **Rolled cones and cylinders:** the sheet between two coaxial arcs (a lofted flange
  between them, or a straight line rolled around an axis), material on either side, unrolled
  on its neutral surface; holes, cut-outs and punches cut in the flat, also across its edges
  as notches and however far round the axis, their walls square to the rolled sheet.
  Or made on a press brake: planar facets joined by bends, as many as a chord tolerance,
  facet angle or facet distance asks for, each bend's ends cut exactly where the arcs lie;
  holes, cut-outs and punches cut in the flat across facets and bends alike.
- **Lofts** between two open polylines (a rectangle to a smaller or shifted one, an open U, any
  convex polygon open at a gap), made on a press brake: a planar facet between each pair of
  matching segments, joined by bends along the lines between matching corners, each bend's
  ends cut exactly where the polylines lie; holes, cut-outs and punches in the flat.
- **Folds** along a line across the plate or a flange, nested, at any angle to the outline; **jogs**
  (Z-steps by an offset).
- **Faces joined by bends:** plates on planes of their own, anywhere in 3D, each joined to the
  base or an earlier face by a bend between two parallel edges. The faces are extended or
  trimmed to meet the bend, with bend reliefs where a face carries on beside it. Faces at an
  angle get one bend; parallel faces two (fixed edges or 45 degrees side by side, 90 degrees or
  a full radius facing each other). Holes and cut-outs on a face, flanges and hems on its edges.
- **Holes and cut-outs** anywhere in the flat, across bends, and across the outline as
  notches, also through a mitered bend's curved end; or placed on a flange's face, measured
  along its edge and up from the mold line.
- **Punches:** a shape of holes and cut-outs stamped at several points, turned by an angle,
  in the flat or on a flange's or face's own face; the flat pattern marks each stamp's centre
  and punch ID.
- **Embosses:** round or rectangular embosses (pads) pressed out of any flat plate (the base,
  a flange, a face), either way, at the sheet's thickness: a drafted or straight wall between
  two rounded corners, or none where the corners meet; domes (a spherical cap) and beads (a
  cylindrical crown with spherical ends) met straight from the corner. Exact (tori, cones,
  cylinders, spheres); the flat pattern marks the footprint.
- **Flanged holes:** a hole's edge pressed out of any flat plate into a collar, either way:
  a corner (sharp or rounded) and a straight or drafted wall, its end cut square to it;
  exact. The flat pattern holds the pre-hole, the collar developed along its neutral fibre.
- **Louvers:** a hood pressed out of any flat plate beside a lance, either way, at the sheet's
  thickness: a quarter round hood (cylinders) with spherical ends, leaving the plate over a
  rounded or sharp corner (tori), closed at the lance by two plane faces; exact. The flat
  pattern marks the footprint, the lance one of its sides.
- **Gussets** across a flange's bend: a rib pressed into the bend's inside, its inner face a
  plane across the corner (by depth from the bend or by its legs), a groove behind it, side
  walls square to the bend; any bend angle below 180, up or down, a wall of its own thickness;
  the flat pattern marks the footprint.
- **Tabs and slots** within one part: tabs on a flange's tip through slots where it lands on
  another flange's face (a self-locating box), the tip standing a gap off the face; count,
  length, height (through by default), clearances, offsets.
- **Cross breaks** on the base or a flange: marked in the flat (DXF layer CrossBreak, SVG),
  corner to corner of the face; no geometry.
- **Bend order:** every bend numbered in the order it's made (outermost first, or as you
  list them), with a note at each bend line in the flat DXF and SVG.
- **Unfolding** by K-factor (per part or per feature) or by a bend table (DIN 6935 style: bend deduction by radius and angle).
- **Shop checks:** bend radius, leg length, hole distance to bends and edges, hole size.
- **Quoting numbers:** the flat's smallest box at any rotation (arcs exact), its area with
  the holes out, cut length (outline and inner cuts apart), pierces and bend count (a fold
  broken by the outline counts once).
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
