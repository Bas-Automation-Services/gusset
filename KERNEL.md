# Gusset - OpenGeometry fit

Read on 2026-10-29 at upstream `46de326` (2026-10-24). That's about 67k lines of Rust, with
`analytic/booleans.rs` alone at 13.6k. The repo is active.

## Two B-rep systems: we use `analytic`

| Module | What | For Gusset |
|---|---|---|
| `brep::Brep` | Faceted: points + face normals, no surfaces. Used by the editor, freeform and the old exporter. | No |
| `analytic::BrepEnvelope` | Exact: surfaces, curves, pcurves, UV trims, per-entity tolerance, `Accuracy` budgets. Serialisable (JSON schema v2), `validate()`. | **Yes, only this** |

## What fits

- **Surfaces:** Plane, Cylinder, Cone, Sphere, Torus, and (fork) Wrapped. Plane and cylinder
  cover flats, bends, side faces and holes in flat regions; Wrapped covers cuts through bends.
- **Curves:** Line, Circle, Ellipse, Intersection (surface-surface), and (fork) Wrapped.
- **Pcurves:** Line2, Conic2, ProjectedCurve, IntersectionSide.
- **Face key + `FaceProvenance`:**
  - Every face has a `key: String`. Gusset's semantic reference goes there, for example
    `flange#3.bend.inner`.
  - Booleans keep provenance (Authored/Preserved/Split/Cut) and report a face mapping.
    References survive the kernel.
- **Builder (internal):** `analytic/primitives.rs` `pub(super) struct Builder` has
  `vertex`, `edge(curve, range, seam)`, `face`, `face_with_holes(key, surface, uv_bounds,
  outer_uses, hole_uses)` and `finish()`, which closes the shell and validates.
  - This is exactly how Gusset wants to deliver its faces: one closed shell built directly,
    with edges shared at the bend tangent lines, and no booleans.
- **Primitives matching sheet metal pieces** (useful as reference and for tests):
  - `annular_sector_extrusion(_with_openings)` is a bend region.
  - `linear_extrusion` and `arc_edged_extrusion_with_holes` are flat regions with holes.
- **STEP:** `analytic::exchange::export_step(&BrepEnvelope, unit)` writes `ADVANCED_FACE` on
  `PLANE` / `CYLINDRICAL_SURFACE` and a `MANIFOLD_SOLID_BREP`. Intersection curves are fitted
  as B-spline curves. It returns a `StepReport` with an error bound.
- **Mesh:** `analytic::tessellation::tessellate` + `TessellationCache` (viewer).
  `query::classify_point` (tests).
- **Errors:** a structured `GeometryError` enum, which is machine-readable.

## Gaps

**Rule:** whatever Gusset needs that the kernel lacks, we build ourselves in the fork, and it
goes upstream as a PR. Gaps never block us.

### No ruled / B-spline surfaces: solved with wrapped geometry (2026-10-30)

The kernel said so itself: "requires a ruled surface outside the current analytic surface
set" (`modeling.rs:231`). A cut through a bend needs:
- an edge on the inner and outer bend faces that is no line, circle or ellipse (its pcurve is
  an ellipse in the cylinder's UV, since a cylinder is developable);
- a **hole wall** that is ruled along radial lines between the two cylinders.

**Fork change 5 (below):** `Wrap` maps flat space onto a cylinder the way a bend does, and
`CurveGeometry::Wrapped` / `SurfaceGeometry::Wrapped` wrap any elementary curve or surface
exactly. A hole wall in a bend is a wrapped cylinder, its edges wrapped circles; a slanted
cut is a wrapped plane bounded by helices (wrapped lines). Exact, not approximated; only
STEP export fits B-splines (within the exchange budget).

**Lofted flange:** still "later"; it is not a wrap.

## Fork changes (each one also goes upstream as a PR)

**Upstream PRs (opened 2026-11-01).** Upstream released 3.0.0 on 2026-10-30, which moved
every file below; each change was ported to 3.0 and opened from a branch off upstream `main`.

**Gusset on 3.0 (2026-11-01).** Gusset uses fork branch `integration/3.0`: upstream 3.0 with
the seven PR branches merged, so the ports are tested by Gusset's tests before review. Fork `main`
stays the pre-3.0 history. Once a PR lands upstream, merge upstream `main` into the branch;
when all have landed, point Gusset back at upstream.
- **Paths:** `analytic::{BrepEnvelope, Frame3, CurveGeometry, ...}`, `analytic::topology::*`,
  `analytic::geometry::Wrap` and `analytic::primitives::{Builder, ...}` are all in `brep`;
  `analytic::tessellation` is `tessellation`, `analytic::exchange` is `exchange`, `Interval`
  is `math::Interval`. Units are strings (`export_step_in(&brep, "millimetre", "millimetre")`).
  `Builder::finish()` is `finish_solid()`.
- **Kept to the kernel in 3.0:** edge geometry (`Edge::geometry`, `GeometryStore::curves`) and
  `Frame3::local` / `vector`. Gusset keeps each line's direction on its own `Edge` when it
  creates it, and has `to_local` / `to_world` itself (`fold.rs`).
- **Fixes made only on the integration branch** (they only matter once PRs are combined):
  a wrapped curve is raised in a short non-planar loop like an intersection curve (#37 with
  #38); #39's orientation test turns whole bodies inside out, since #40's `validate()` already
  rejects a single flipped face. #38 itself got one fix (`8a22f57`): `shell.rs` stays at 600
  lines with #40.

| Change | Upstream PR | Fork branch |
|---|---|---|
| 1. Builder public | #36 | `feat/public-brep-builder` |
| 2 + 9. Orientation check (in `query::validity`, not in the builder: layer rules) | #39 | `fix/validate-face-orientation` |
| 3. Standard loop winding (not breaking: `from_json` turns old loops) | #40 | `fix/standard-loop-winding` |
| Boolean fixes (see Known kernel issues) | #41 | `fix/closed-loop-hole-booleans` |
| 4. STEP in the model's unit | #35 | `feat/step-model-unit` |
| 6 + 8. Curved trim tessellation | #37 | `fix/curved-trim-tessellation` |
| 5 + 7 (+ wrapped arcs of 8). Wrapped geometry | #38 | `feat/wrapped-geometry` |

License: **MPL-2.0** (confirmed by the maintainer). The `Cargo.toml` field still says
MIT/Apache; we can fix that along the way in a PR.

1. **Builder public: done** (`20b6c09` on fork `main`).
   - Public: `analytic::primitives::{Builder, Use, boundary, plane_boundary, uv_line}`.
   - `face()` returns the face id, plus `set_face_sense()` and `brep()`.
   - Test: `tests/public_brep_builder.rs`.
2. **Geometric orientation probe: done** (`f5a01bf`, branch
   `fix/validate-face-orientation`).
   - `query::validate_orientation`: per face, a point just in front of the face along its
     material normal must not be inside the solid while the point just behind it is
     outside. `Builder::finish` runs it.
   - It catches a whole solid that is inside-out, which the winding rule (3) can't see.
   - **Skipped** for bodies the tessellator refuses, and when the tolerance is below ~1e-12 of
     the model size.
   - **Finding:** at that scale `classify_point` is unreliable. On a 900 km box with
     tolerance 1e-9 it gives wrong results at the face centre.
3. **Standard loop winding: done** (`d5aa686`, branch `fix/standard-loop-winding`, breaking).
   - **Old convention:** loop winding did **not** say which side is material. Holes wound the
     same way as outer loops, and a `Reverse` face kept its loop. So a flipped face went
     unnoticed, and STEP export wrote hole bounds the wrong way per ISO 10303-42.
   - **Now (standard):** seen against the material normal, outer loops run counter-clockwise
     and holes clockwise. `validate()` enforces this.
   - **`BrepEnvelope::orient_loops()`** measures the winding per loop in UV and turns it if
     needed. `Builder::finish` and every Boolean result call it, so producers may build in
     either direction.
   - **Finding (upstream bug):** the cone ∩ sphere intersection
     (`generic_closed_loop_handles_periodic_winding`) already gave an inside-out cone patch
     on upstream `main`: a point inside the result classifies as Outside. The test is
     `#[ignore]` with that reason. Candidate for a separate fix plus upstream issue.
   - **Breaking:** old v2 JSON with the old winding is rejected by `from_json`.
   - **Gusset:** builds loops standard and checks them; STEP holes are now correct per spec.
4. **STEP in the model's unit: done** (`ea556c1`, branch `feat/step-model-unit`).
   - Analytic bodies are in metres by convention, so `export_step(.., "millimetre")` multiplied
     Gusset's millimetre coordinates by 1000.
   - `export_step_in(brep, model, file)` now takes the model's `LengthUnit`. Gusset exports
     millimetres unscaled.
   - STEP of all 32 test parts imports as one healthy solid with exactly its own volume.
5. **Wrapped geometry: done** (`f0927b7`, branch `feat/wrapped-geometry`).
   - `geometry::Wrap { frame, radius }`: flat `(x, y, z)` goes to
     `frame.point([z cos(x/r), z sin(x/r), y])`; x is arc length on radius r (Gusset uses the
     neutral layer, so flat lengths are kept there). Orientation-preserving.
   - `CurveGeometry::Wrapped` / `SurfaceGeometry::Wrapped`: point, derivatives (chain rule,
     second order), normal, projection (exact on the surface), interval bounds.
   - Tessellation: curvature from the exact second derivatives, sampled with a 1.5x margin
     (a strict bound is far too loose along the rulings).
   - Point classification: UV subdivision plus Newton (no implicit form), so the
     orientation probe works on wrapped faces.
   - STEP: bicubic B-spline surfaces over the face's own UV (pcurves stay valid), cubic
     B-spline edges; Hermite fits checked against the exact geometry.
   - Booleans, SSI and IFC refuse wrapped geometry (`UnsupportedGeometry`); placement only
     moves it rigidly.
   - Tests: wrapped families in the per-family geometry tests; `tests/wrapped_solid.rs`.
6. **Curved trim rings snapped: done** (`8c58c8d`, branch `feat/wrapped-geometry`).
   - Samples along an axis-parallel pcurve could differ from the neighbouring corner vertex
     by rounding; earcut then cut zero-area ears along the run and the mesh failed with
     "triangle orientation could not be resolved". Seen on a bend trimmed by a hole next to a
     corner. Coordinates within 64 eps of the previous sample are snapped.
   - The refinement-limit and orientation errors now name the face.
7. **Fitted wrapped surfaces reach past their faces: done** (`ce23aff`).
   - Some importers rebuild trims by intersecting neighbouring surfaces (ignoring pcurves). A
     B-spline fitted exactly over the face ended on the face's edges, so they trimmed short:
     up to 4e-5 off in volume. Fitted 10% past the faces (not around a full turn), STEP
     reads back within ~3e-7.
8. **Arcs get a segment per quarter turn: done** (`512cf48`).
   - Found by fuzzing: a small hole wholly inside a bend, under the orientation probe's
     coarse deflection, got one segment per half circle, a two-sample trim loop ("curved
     trim loop requires at least four samples"). Circles, ellipses and wrapped circles now
     get at least one segment per quarter turn.
9. **Orientation probe confirms with a second point: done** (`4e220eb`).
   - Found by fuzzing: the probe judged a face from one point (its largest triangle's
     centroid) and one ray; a ray grazing a hole's edge on a folded sheet miscounted and
     rejected sound parts (exact and mesh volumes agreed to 6e-6). A face is now reported
     only when its second largest triangle's centroid agrees.
10. **Conoid surfaces: done** (branch `feat/conoid`, merged into `integration/3.0` as
    `208a416`; not opened upstream yet).
   - `SurfaceGeometry::Conoid { frame, c }`: S(u, v) = O + v (cos u x + sin u y) +
     (c0 + c1 cos u + c2 sin u) z, UV = (angle, radius): radial rulings over a cylinder axis.
     A bend mitered at a corner ends on one, cut square to the sheet.
   - `conoid_ellipse(frame, c, radius)`: its edge at a radius as an exact ellipse (parameter
     = angle - phase). `projected_boundary(...)`: a projected pcurve for any edge on a face
     (here the ellipse on the bend's cylinder).
   - Projection, tessellation (grid sized for the twist), classification, placement;
     STEP as a bicubic B-spline fit like Wrapped; booleans, SSI and IFC refuse it.
11. Anything else we run into while building.

## Known kernel issues (2026-10-29)

**Where they are.** Almost all of them are in the **curved booleans**: sphere, cone and
torus with each other or with a cylinder. Upstream tests there mostly check `validate()`
and tessellation, not geometry (volume, inside/outside), so wrong solids passed.

**What Gusset uses** is mostly the builder, STEP export and tessellation; booleans hardly
at all. **Rule:** never trust kernel output without an independent check (exact volume,
probe).

**Fixed** (`04f16d6` on fork `main`, branch `fix/cone-sphere-intersection`):
- cone ∩ sphere kept the wrong sphere cap.
- cylinder ∩ sphere could come out inside-out.
- A hole ended up a period outside the face, so cone ∪ sphere was wrong by a whole cap.

**Refused, not computed:** a hole loop that crosses the seam of a kept face (u = 0/2π).
The kernel now returns `CoverageGap` instead of a wrong solid. Before, this gave wrong
solids that passed validation.

**Tessellation limit (2026-10-30).** A curved face that isn't a UV rectangle (a bend
trimmed by a hole) is meshed by earcut plus isotropic refinement. At 1e-5 mm deflection
that needs millions of triangles and hits the limit; 1e-4 mm takes about a second. Fine for
display and the orientation probe. For exact numbers Gusset integrates mass properties over
the layout instead (`fold::exact_mass`). Bisecting only each triangle's worst edge was tried
and was slower.

**Still open** (all pre-existing, not on Gusset's path):
- Duplicate face keys when both inputs have faces with the same key (e.g. `skin` on
  sphere and torus): "face IDs or keys are invalid".
- A sphere with its centre just inside a cylinder wall: some points are misclassified and
  the orientation probe fails.
- `classify_point` / `validate_orientation` sometimes throw "UV outside surface domain".
- The intersection solver fails on tilted or skewed sphere frames ("correction leaves
  branch tube").
- `classify_point` is unreliable when the tolerance is below about 1e-12 of the model size.

**Alternatives** (checked 2026-10-29):
- No mature pure-Rust kernel exists; truck and Fornjot are no further along.
- OpenCASCADE is mature, but as WASM it is about 30 MB.

**Decision: keep OpenGeometry.** Reassess once its STEP output has been tested on import.

## Practical

- Gusset uses a git dependency on the fork's `integration/3.0` branch (crate 3.0.0).
- **WASM:** OpenGeometry is a `cdylib` + `rlib` with `wasm-bindgen`. `gusset-wasm` links it as
  an rlib, so pin the same `wasm-bindgen` version.
- **Fork:** `github.com/Bas-Automation-Services/OpenGeometry`, locally at
  `D:\Repositories\BasAutomationServices\OpenGeometry` (remote `upstream` = OpenGeometry-io).
- **Built on this machine (2026-10-29, `46de326`):**
  - `cargo test --release --no-fail-fast`: 465 passed. Without `--no-fail-fast`, cargo stops
    at `occ_oracle` and skips the test binaries after it.
  - `occ_oracle` fails on Windows because `python3` resolves to the Microsoft Store stub,
    which errors instead of letting the test skip. Not a kernel problem, and a small
    upstream fix: treat "no stdout + failure" as a skip.
  - `wasm-pack build --target web --release` OK; `opengeometry_bg.wasm` is 3.0 MB.
