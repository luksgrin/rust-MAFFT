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

103 FFI test functions across 15 files (counted from `#[test]`
attributes in `crates/*/tests/cross_validate*.rs`):

| Layer | File (tests) | What's compared |
|--|--|--|
| Scoring | `mafft-scoring/tests/cross_validate.rs` (21) | BLOSUM45/50/62/80, JTT, TM and DNA `n_dis` matrices cell-by-cell (plain and FFT variants), penalty values, amino-acid index mapping |
| Weights | `mafft-tree/tests/cross_validate_weights.rs` (3) | Branch-length weights, per-group normalisation |
| Distance / trees | `mafft-tree/tests/cross_validate_parttree.rs` (9), `cross_validate_memsavetree.rs` (5) | k-mer (sextet) distances, PartTree pivots and UPGMA order, memsavetree / youngestlinkage topology step-by-step |
| Column profiles | `mafft-core/tests/cross_validate_cpmx.rs` (3), `cross_validate_profile.rs` (4) | `cpmx` column profile matrices, inter-group scores, unnormalised / unequal weights |
| Pairwise & constrained DP | `mafft-core/tests/cross_validate_constrained_align.rs` (12) | `G__align11` (incl. warp variants), `genL__align11` (E-INS-i params), constraint-mode `A__align` with gaps and multi-member groups |
| Profile DP | `mafft-core/tests/cross_validate_profile_align.rs` (5), `cross_validate_msalign.rs` (10) | `A__align` with `imp` matrix / non-zero `penalty_ex` / shifted matrices, `MSalignmm` recursion, base case, free tail, asymmetric lengths |
| FFT | `mafft-core/tests/cross_validate_fft.rs` (2), `cross_validate_bl50_fft.rs` (2) | Hand-ported Cooley-Tukey vs C `fft.c`; alignable-region detection and step-24 profile DP on the `--bl 50` sentinel |
| Refinement | `mafft-core/tests/cross_validate_counteff.rs` (2), `cross_validate_bb20027_dp.rs` (2), `cross_validate_refine_branch_dna.rs` (1) | Pass-1 widths and `counteff` weights on BB20027, segment-10 DNA refinement profile DP |
| `--add` machinery | `mafft-core/tests/cross_validate_insertnewgaps.rs` (1) | New-merge-gap insertion, `profilealignment` compression cell-by-cell |

Three further `mafft-core` test binaries also link the C shim but are
forensic rather than parity tests: `trace_refinement.rs` (15),
`exp_residual_dump_compare.rs` (1) and `exp_residual_step10.rs` (1).

The fixture-based end-to-end tests do **not** use the FFI; they compare
against committed C-canonical outputs: `mafft-core/tests/end_to_end.rs`
(120), `mafft-bin/tests/c_parity_dna.rs` (5), and the `imp` matrix
tests in `mafft-align/tests/` (3).

## Platform-relative comparison

C MAFFT 7.526 itself differs between clang/arm64 (which contracts
`a*b + c` into fused multiply-adds) and gcc/x86-64 (which does not) —
see [Byte-identity is per platform
build](byte-identity.md#byte-identity-is-per-platform-build). The FFI
tests compare Rust against the C that `mafft-c-bindings/build.rs`
compiled on the *same* machine, so they are meaningful only when Rust's
policy matches that C build. The default `mafft_types::fp::CONTRACTS_FMA`
(true on `aarch64`, false elsewhere) does exactly that on both CI
runners. Do not run the `cross_validate_*` binaries with
`--features fp-contract-none` on arm64 or `--features fp-contract-fma`
on x86-64: the C side would still be using the host compiler's policy
and the comparison would be between two different policies.

The fixture-based tests, by contrast, *can* be run under the foreign
policy: `fixture_path_fp` selects `<name>.fma` or `<name>.nofma` by
`CONTRACTS_FMA`, so `--features fp-contract-none` on an arm64 host
checks the x86-64 references. That is what the CI `policy-cross-check`
job does.

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
cargo test -p mafft-core --release --test cross_validate_fft

# Just one test
cargo test -p mafft-core --release --test cross_validate_fft -- fft_port_matches_c_forward_dc

# Fixture tests under the x86-64 policy on an arm64 host (non-FFI only;
# see "Platform-relative comparison" above)
cargo test -p mafft-core --release --features fp-contract-none --test end_to_end
```

## CI integration

The `ci.yml` workflow runs the whole suite on every push and pull
request (documentation-only changes excepted) on two
reference platforms via the reusable `build-test.yml` job: `Build
(x86-64)` on `ubuntu-latest` (gcc, no contraction) and `Build (arm64)`
on `macos-latest` (Apple clang, FMA contraction). Each leg builds the C
submodule in-tree, prints an FMA census of the binaries it just built,
runs `cargo test --workspace --lib` and `cargo test --workspace --release
--tests`, and finishes with an end-to-end sentinel that diffs
`mafft-rs --bl 50 --retree 2 --maxiterate 0` against the same-platform C
`mafft` wrapper (712 columns on arm64, 738 on x86-64).

A parallel set of x86-64 jobs exercises the C MAFFT binaries on the
same `mafft-upstream/test/sample` input and compares to canonical output
files, and a `policy-cross-check` job runs the non-FFI fixture tests on
arm64 under `--features fp-contract-none`. All must pass.

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
