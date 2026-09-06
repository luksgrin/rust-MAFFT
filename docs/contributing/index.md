# Contributing

Issues, pull requests, and questions all land at
[github.com/luksgrin/rust-MAFFT](https://github.com/luksgrin/rust-MAFFT).
This page covers what you need to know to get a working development
environment and ship a change.

## Local setup

```sh
git clone --recurse-submodules https://github.com/luksgrin/rust-MAFFT
cd rust-MAFFT

# Sanity-build the Rust workspace
cargo build --workspace --release

# Full test suite (~430 tests; ~30s on a modern machine)
cargo test --workspace --release --tests
```

The `mafft-upstream/` submodule pins the reference C MAFFT version that
the [cross-validation harness](../architecture/cross-validation.md)
builds against. You need a C compiler on `$PATH` (`cc` is auto-detected
by the `cc` crate) — `apt install build-essential`, `xcode-select
--install`, or equivalent.

## Repository layout

```
rust-MAFFT/
├── crates/                  ← the 12-crate workspace (see Architecture)
├── docs/                    ← this site's source
├── mkdocs.yml
├── mafft-upstream/          ← submodule: reference C MAFFT
├── .github/workflows/
│   ├── ci.yml               ← push/PR (non-docs changes): cargo build + test
│   ├── docs.yml             ← push/PR touching docs inputs: build; main: deploy
│   ├── release.yml          ← release: 5 pre-built binaries to GH releases
│   ├── python.yml           ← push/PR: 5 wheels (py3.12); release: 25 wheels + sdist to PyPI
│   └── cargo-publish.yml    ← release: 10 crates to crates.io
├── README.md
└── TODO.md                  ← parity roadmap, working scratchpad
```

### CI on pull requests

`ci.yml` and `python.yml` run on pushes to `main`/`dev` and on pull
requests against *any* base branch (so stacked PRs get checks), but skip
changes that touch only `docs/**`, `**/*.md`, `mkdocs.yml`,
`CITATION.cff` or `LICENSE*`. `docs.yml` is the mirror image: it builds
only when something the site is generated from changes (`docs/**`,
`mkdocs.yml`, `crates/mafft-bin/src/**` for the CLI reference,
`crates/pymafft/**` for mkdocstrings, `scripts/gen-cli-reference.sh`, the
workflow itself) and always on release or manual dispatch.

On push/PR `python.yml` builds one wheel per target for Python 3.12 and
tests that wheel; the full 5 targets x Python 3.9-3.13 matrix runs on
`release: published` and `workflow_dispatch`. The regime is chosen by
the small `matrix` job at the top of the workflow, which emits the
Python-version lists that `build-wheels` / `test-wheels` read through
`fromJSON`.

## Building the docs locally

```sh
pip install mkdocs-material 'mkdocstrings[python]' biopython

# Build the pymafft extension so mkdocstrings can introspect it
cd crates/pymafft && pip install maturin && maturin develop --release
cd ../..

# Regenerate the CLI reference from the current binary
cargo build --release -p mafft-rs
./scripts/gen-cli-reference.sh > docs/cli-reference.md

# Serve at http://localhost:8000
mkdocs serve
```

## Coding guidelines

- **Don't break byte-identity.** Every change to the engine must keep
  the full cross-validation suite green. See
  [Byte-identity](../architecture/byte-identity.md) for the rationale.
- **No `unsafe` outside FFI.** The engine is `#![forbid(unsafe_code)]`
  at the module level. The only `unsafe` is in `mafft-c-bindings`
  (which is dev-only).
- **No `unwrap`/`expect` in library code.** Use proper error returns;
  the engine should never panic on real input.
- **Match C MAFFT structure where useful, not blindly.** Where the C
  code is genuinely confusing or relies on accidental globals,
  refactor to idiomatic Rust as long as the test suite passes.

## Adding a new flag

1. Add the field to `Args` in `crates/mafft-bin/src/lib.rs` (clap derive).
2. Wire it through to the engine via the appropriate `AlignmentMode` /
   knob.
3. Add an end-to-end test in `crates/mafft-core/tests/end_to_end.rs`.
4. If the flag changes alignment output, add a fixture-based
   cross-validation test that compares to the canonical C output
   (run upstream MAFFT with the same flag to generate the fixture).
5. Run `./scripts/gen-cli-reference.sh > docs/cli-reference.md`
   so the docs reflect the new flag.

## Release process

See [`docs/architecture/workspace.md`](../architecture/workspace.md#release-flow)
for the fan-out diagram. The button-press procedure:

1. **Bump version** in the root `Cargo.toml` `[workspace.package].version`
   field. (The `mafft-sys` stub is pinned at 0.0.1 and not driven by
   the workspace version.)
2. **Bump pymafft/pyproject.toml** `[project].version` to match.
3. **Commit**, tag, push:
   ```sh
   git commit -am "chore: bump version to vX.Y.Z"
   git tag vX.Y.Z
   git push --tags
   ```
4. **Create the GitHub release** from the tag. Publishing triggers:
    - `release.yml` → 5 pre-built binaries attached to the release page
    - `python.yml` → 25 wheels + sdist published to PyPI (OIDC)
    - `cargo-publish.yml` → 10 crates published to crates.io
    - `docs.yml` → docs site rebuilt on `main`
5. **Watch the workflows.** If a publish fails partway through, fix it
   and re-run the workflow (cargo and PyPI both tolerate re-pushes of
   the same version only if NONE of it landed; otherwise bump the
   patch).

## Required secrets

In the GitHub repo settings:

- **`CARGO_REGISTRY_TOKEN`** — for `cargo-publish.yml`. Generate at
  [crates.io → Account Settings → API Tokens](https://crates.io/me).
- **PyPI trusted publishing** — for `python.yml`. Register the repo +
  `python.yml` + environment `pypi` at
  [pypi.org → Manage publishing](https://pypi.org/manage/account/publishing/).
  No token needed.

Pages-related permissions are configured per-repo under Settings → Pages
(source: GitHub Actions).
