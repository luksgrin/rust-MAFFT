# Cross-validation

The byte-identity claim is held in place by a test suite that imports
C MAFFT itself and runs both implementations on the same inputs,
comparing outputs at multiple granularities.

## The plumbing

The internal `mafft-c-bindings` crate (formerly `mafft-sys`) has a
`build.rs` that compiles MAFFT's C source from the `mafft-upstream/`
submodule into a static library:

```
mafft-upstream/core/*.c  ─►  mafft-c-bindings/build.rs  ─►  libmafft_c.a
                                                                 │
                                                                 ▼
                                                          unsafe extern "C"
                                                          bindings in Rust
```

The bindings are dev-deps of `mafft-scoring`, `mafft-tree`, `mafft-core`
(via a cargo-aliased name `mafft-sys` so test code can still write
`use mafft_sys::*`).

## What gets compared

Roughly ~140 tests across `crates/*/tests/cross_validate_*.rs`:

| Layer | Files | What's compared |
|--|--|--|
| Scoring | `cross_validate_weights.rs`, `cross_validate_calcw.rs` | Per-cell weight vectors, branch lengths, `calcW` SP-component |
| Distance | `cross_validate_parttree.rs`, `cross_validate_dist*.rs` | k-mer distances, PartTree pivots, UPGMA pair-merge order |
| Pairwise DP | `cross_validate_galign11.rs`, `cross_validate_lalign11.rs` | NW / SW pairwise alignments — final score AND traceback |
| Profile DP | `cross_validate_profile_align.rs`, `cross_validate_msalign.rs` | Profile DP outputs row-by-row |
| FFT | `cross_validate_falign.rs`, `cross_validate_localhom.rs` | FFT-anchored DP including localhom constraint integration |
| `--add` machinery | `cross_validate_insertnewgaps.rs` | New-merge-gap insertion, profilealignment |
| End-to-end | `crates/mafft-core/tests/end_to_end.rs` | Full alignment outputs vs C-canonical fixtures |

## Why per-function FFI tests, not just end-to-end

End-to-end tests catch *that* something differs. Per-function FFI tests
catch *which step* differs and *by how much* — which is everything when
debugging a 1-byte divergence on `BB30013`.

The pattern across cross-validate files:

1. Set up inputs identically on both sides (Rust types + C globals).
2. Call the Rust function.
3. Call the C function via FFI.
4. Compare outputs at the highest fidelity the function supports
   (often `assert_eq!` on byte slices; sometimes `f64::to_bits()`
   equality on per-cell DP scores).

When step 4 fails, you have a localized bug. Most of the residual
parity bugs that took byte-identity from 95% to 100% were found this
way.

## Running the suite

```sh
# All tests (workspace, release mode — the FFI shim is opt-level=3)
cargo test --workspace --release --tests

# Just one layer
cargo test -p mafft-fft  --release --tests
cargo test -p mafft-core --release --tests

# Just one file
cargo test -p mafft-core --release --test cross_validate_falign

# Just one test
cargo test -p mafft-core --release --test cross_validate_falign -- falign_matches_c_basic
```

## CI integration

The `ci.yml` workflow runs the whole suite on every push (`ubuntu-latest`),
alongside a parallel set of jobs that exercise the C MAFFT binaries on
the same `mafft-upstream/test/sample` input and compare to canonical
output files. Both must pass.

## Dump-mode debugging

For the trickier divergences (e.g., FFT cut selection differing on
flat-landscape inputs), several test files support env-var-gated
forensic dumps:

```sh
RS_DP_DUMP=/tmp/rs_dp.txt cargo run --release -p mafft-rs -- input.fa
RS_SUPPORT_DUMP=/tmp/supports.txt ...
```

The C source has matching `getenv()` blocks (locally patched in
`mafft-upstream/core/`) so the same dump can be produced on the C side.
Diffing the two `/tmp/rs_dp.txt` files pinpoints which merge step
diverges.

This is mostly archaeological at this point — every BAliBASE fixture
passes — but the infrastructure stays in place for future regressions.
