# Gusset Workbench

A static page that runs the kernel as WebAssembly: edit a part as JSON and see it folded,
its flat pattern, mass and bend table; copy the flat DXF or the STEP.

Build the package next to the page, then serve the folder:

    wasm-pack build crates/gusset-wasm --target web --release --out-dir ../../demo/pkg
    python -m http.server -d demo 8000

The examples in `index.html` are parts from Gusset's tests.
