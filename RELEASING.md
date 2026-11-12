# Releasing gusset-wasm

`crates/gusset-wasm` is published to npm as [`gusset-wasm`](https://www.npmjs.com/package/gusset-wasm)
by `.github/workflows/publish-wasm.yml` whenever a `v*` tag is pushed.

Publishing uses npm trusted publishing (OIDC): no npm token is stored anywhere. npm cannot set up a
trusted publisher for a package that does not exist yet, so the very first version goes out with a
short-lived token.

## Each release

1. Set `version` in `crates/gusset-wasm/Cargo.toml` to the new version and merge to `main`.
2. Tag that commit `v` + the same version and push the tag:

   ```sh
   git tag v0.2.0
   git push origin v0.2.0
   ```

That is all: no secrets, no login.

## First release (once)

1. On npmjs.com create a granular access token with read and write access to packages, expiring
   in 7 days. Add it to the GitHub repository as the Actions secret `NPM_TOKEN`.
2. Tag and push `v0.1.0` as above. The workflow publishes `gusset-wasm@0.1.0` with the token.
3. On npmjs.com open `gusset-wasm` → Settings → Trusted publishing, choose GitHub Actions and
   enter organisation `Bas-Automation-Services`, repository `gusset`, workflow filename
   `publish-wasm.yml` (no environment).
4. In the same settings, under Publishing access, choose "Require two-factor authentication and
   disallow tokens". Trusted publishing keeps working; tokens no longer can.
5. Delete the `NPM_TOKEN` secret from the repository and revoke the token on npmjs.com.

A new trusted publisher must complete its first publish within 2 days, or it expires and has to
be deleted and created again. If no release follows within 2 days of step 3, recreate it just
before the next release.

## What the workflow does

1. Fails unless the tag is exactly `v` + the `gusset-wasm` crate version.
2. Upgrades npm to 11.21.1 (trusted publishing needs 11.5.1 or later; Node 22 ships npm 10).
3. Runs `cargo test --release --workspace`.
4. Builds the package with `wasm-pack build crates/gusset-wasm --target web --release`.
   wasm-pack writes `package.json` from the crate's name, version, licence, description and
   repository, and copies `crates/gusset-wasm/LICENSE` (a copy of the root `LICENSE`; wasm-pack
   only looks in the crate folder) into the package.
5. Publishes `crates/gusset-wasm/pkg` with `npm publish --access public --provenance`. npm first
   exchanges the workflow's OIDC token for a publish token; that wins over `NPM_TOKEN`. Without a
   trusted publisher the exchange fails quietly and npm uses `NPM_TOKEN`. With neither, it stops
   with "This command requires you to be logged in" (add `--loglevel verbose` to see the `oidc`
   lines explaining why the exchange failed). `--provenance` is set explicitly so the token-published
   first version gets a provenance statement too; under trusted publishing npm would add one anyway.

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
