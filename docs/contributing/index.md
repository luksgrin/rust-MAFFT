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

# Quick check: the CLI crate's unit tests (no C compiler needed)
cargo test -p mafft-rs --lib
```

The `mafft-upstream/` submodule pins the reference C MAFFT version that
the [cross-validation harness](../architecture/cross-validation.md)
builds against. The full test suite needs it, plus a C compiler; see
[Reproduce the cross-validation suite](#reproduce-the-cross-validation-suite).

The pin only moves through the procedure in
[Upstream policy](../architecture/upstream-policy.md#upgrade-procedure):
`check-mafft-upstream.yml` opens an `upstream: …` issue when the pin
falls behind upstream, and the bump is a normal PR that re-derives the
C-canonical fixtures and reruns the BAliBASE sweep. Don't move the
submodule as a side effect of another change.

### Python bindings

`pymafft` is built with [maturin](https://www.maturin.rs/) into a virtual
environment. Only the standard library's `venv` and `pip` are needed:

```sh
cd crates/pymafft
python -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate
pip install maturin pytest biopython
cargo build --release -p mafft-rs   # the CLI the option-parity tests compare against
maturin develop --release
pytest tests/ -v
```

`maturin develop` installs the extension module but not the bundled CLI, so
the console-script tests skip. To include it, copy the binary in before
`maturin develop`:
`mkdir -p python/pymafft/_bin && cp ../../target/release/mafft-rs python/pymafft/_bin/`.
If you use [uv](https://docs.astral.sh/uv/), `uv venv` / `uv pip install`
work in place of the `venv` / `pip` steps.

## Reproduce the cross-validation suite

This re-checks rust-MAFFT against C MAFFT 7.526. It is not needed to use
the package; to check an install, see
[Verify your installation](../getting-started/install.md#verify-your-installation).

Prerequisites:

- the `mafft-upstream/` submodule (`git clone --recurse-submodules`, or
  `git submodule update --init` in an existing clone);
- a C compiler on `$PATH` (`cc` is auto-detected by the `cc` crate) —
  `apt install build-essential`, `xcode-select --install`, or
  equivalent. CI builds the C reference with GCC 13, 14 and 15 on Linux
  and Apple clang on macOS; any `CC=gcc-NN` works the same way. Do not
  put a `-std=` flag in `CFLAGS`: the `cc` crate appends `CFLAGS` after
  the shim's own `-std=c99`, and the upstream sources do not compile as
  C23 (GCC 15's default);
- time: the workspace suite takes a few minutes (about 1.5 minutes on an
  Apple M4 with the Rust dependencies already built). The BAliBASE sweep
  runs both implementations on every test set in every mode and takes
  much longer.

Results are relative to the C build on your platform; see
[Byte-identity is per platform build](../architecture/byte-identity.md#byte-identity-is-per-platform-build).

### Workspace suite

```sh
# 537 tests + 7 ignored across 48 suites; builds the C reference in-tree
# (mafft-c-bindings/build.rs) on the first run
cargo test --workspace --exclude pymafft --release
```

[Cross-validation](../architecture/cross-validation.md#reproducing-the-suite)
lists how to run a single layer, file or test.

### FP-policy variants

```sh
# Fixture tests under the other C compiler's contraction policy.
# On Apple Silicon the default fuses, so fp-contract-none matches gcc/MSVC.
# On Intel macOS, Linux (including aarch64) and Windows the default does not fuse,
# so fp-contract-fma matches Apple clang. FFI tests that compare
# against the in-tree C build skip when the policies do not match.
cargo test --workspace --exclude pymafft --release --features fp-contract-none

# Re-derive the policy-sensitive fixtures from the in-tree C build and diff
# them against the checked-in variant for this host (`write` overwrites)
make -C mafft-upstream/core
scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft check
```

### BAliBASE sweep

BAliBASE itself is not in the repository; the docstring of
`scripts/balibase_parity.py` says where to get it and the layout it
expects. The sweep compares `target/release/mafft-rs` against the in-tree
C build, mode by mode:

```sh
make -C mafft-upstream/core
export MAFFT_BINARIES="$PWD/mafft-upstream/binaries"
cargo build --release -p mafft-rs
scripts/balibase_parity.py /path/to/BAliBASE/ --mafft mafft-upstream/scripts/mafft \
    --modes "--maxiterate 1000" "--localpair --maxiterate 1000" \
            "--globalpair --maxiterate 1000" "--genafpair --maxiterate 1000" \
            "--globalpair --maxiterate 1000 --allowshift"
```

Add `--filter RV11 --limit 10` for a smoke run. Per-test results go to
`balibase_parity_results.tsv`; the exit code is non-zero if any test
diverges.

## Repository layout

```
rust-MAFFT/
├── crates/                  ← the 12-crate workspace (see Architecture)
├── docs/                    ← this site's source
├── mkdocs.yml
├── mafft-upstream/          ← submodule: reference C MAFFT
├── .github/workflows/
│   ├── ci.yml               ← push/PR (non-docs changes): build + test on x86-64 and arm64 (via build-test.yml)
│   ├── build-test.yml       ← reusable job: C in-tree, cargo test, FMA census, sentinel, fixture check
│   ├── check-mafft-upstream.yml ← weekly: open/update an issue when the submodule pin lags upstream
│   ├── docs.yml             ← push/PR touching docs inputs: build; main: deploy
│   ├── release.yml          ← release: 5 pre-built binaries to GH releases
│   ├── python.yml           ← push/PR: 5 abi3 wheels, tests on py3.12; release: same + sdist to PyPI
│   └── cargo-publish.yml    ← release: 10 crates to crates.io
├── README.md
└── TODO.md                  ← parity roadmap, working scratchpad
```

### CI on pull requests

`ci.yml` and `python.yml` run on pushes to `main` and on pull
requests against *any* base branch (so stacked PRs get checks), but skip
changes that touch only `docs/**`, `**/*.md`, `mkdocs.yml`,
`CITATION.cff` or `LICENSE*`. `docs.yml` is the mirror image: it builds
only when something the site is generated from changes (`docs/**`,
`mkdocs.yml`, `crates/mafft-bin/src/**` for the CLI reference,
`crates/pymafft/**` for mkdocstrings, `scripts/gen-cli-reference.sh`, the
workflow itself) and always on release or manual dispatch.

`python.yml` builds one stable-ABI (`abi3`, tagged `cp39-abi3`) wheel
per target on every run, so there is no Python axis on the build. The
tests install that same wheel on Python 3.12 on push/PR, and on 3.9, 3.12
and 3.14 on `release: published` and `workflow_dispatch`. The regime is
chosen by the small `matrix` job at the top of the workflow, which emits
the Python-version list that `test-wheels` reads through `fromJSON`. A
`lint` job runs `cargo check -p pymafft` with `RUSTFLAGS="-D warnings"`,
so PyO3 deprecations fail CI as soon as they appear.

## Building the docs locally

```sh
# `maturin develop` needs an active virtual environment
python -m venv .venv
source .venv/bin/activate           # Windows: .venv\Scripts\activate
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
    - `python.yml` → 5 abi3 wheels + sdist published to PyPI (OIDC)
    - `cargo-publish.yml` → 10 crates published to crates.io
    - `docs.yml` → docs site rebuilt on `main`
    - Zenodo (GitHub integration) → the release is archived and gets a
      version DOI under concept DOI 10.5281/zenodo.20623683; add the new
      version DOI to `CITATION.cff` (`doi:` and `identifiers:`) in a
      follow-up commit, together with the `version` / `date-released` bump
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
