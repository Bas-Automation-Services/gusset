# Releasing gusset-wasm

`crates/gusset-wasm` is published to npm as [`gusset-wasm`](https://www.npmjs.com/package/gusset-wasm)
by `.github/workflows/publish-wasm.yml` whenever a `v*` tag is pushed.

## Once

Create an npm access token that may publish `gusset-wasm` (a granular token, or an automation
token) and add it to the GitHub repository as the Actions secret `NPM_TOKEN`.

## Each release

1. Set `version` in `crates/gusset-wasm/Cargo.toml` to the new version and merge to `main`.
2. Tag that commit `v` + the same version and push the tag:

   ```sh
   git tag v0.1.0
   git push origin v0.1.0
   ```

## What the workflow does

1. Fails unless the tag is exactly `v` + the `gusset-wasm` crate version.
2. Runs `cargo test --release --workspace`.
3. Builds the package with `wasm-pack build crates/gusset-wasm --target web --release`.
   wasm-pack writes `package.json` from the crate's name, version, licence, description and
   repository, and copies `crates/gusset-wasm/LICENSE` (a copy of the root `LICENSE`; wasm-pack
   only looks in the crate folder) into the package.
4. Publishes `crates/gusset-wasm/pkg` with `npm publish --access public --provenance`.

`publish = false` in the crate only keeps it off crates.io; wasm-pack and npm do not read it.

## Using it from Vite

The package is an ES module that fetches `gusset_wasm_bg.wasm` next to itself. `vite build` handles
that, but the dev server pre-bundles dependencies and the default URL then misses. Either exclude
the package from pre-bundling:

```js
// vite.config.js
export default { optimizeDeps: { exclude: ["gusset-wasm"] } };
```

or pass the wasm URL yourself, which works in dev and build alike:

```js
import init, { step } from "gusset-wasm";
import wasmUrl from "gusset-wasm/gusset_wasm_bg.wasm?url";

await init({ module_or_path: wasmUrl });
```

Under Node, load the bytes and use `initSync`:

```js
import fs from "node:fs";
import { createRequire } from "node:module";
import { initSync } from "gusset-wasm";

const path = createRequire(import.meta.url).resolve("gusset-wasm/gusset_wasm_bg.wasm");
initSync({ module: fs.readFileSync(path) });
```
