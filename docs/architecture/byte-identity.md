# Byte-identity

rust-MAFFT's output matches C MAFFT 7.526 byte-for-byte across BAliBASE
3 (1930/1930 fixtures) for every supported mode. This page explains why
that's the bar and what it cost to hit it.

## Why hold to byte-identity

MAFFT is the *de facto* MSA reference. Every downstream tool and
pipeline — phylogenetics, structure prediction, comparative genomics —
treats its output as ground truth. A Rust port that's "close" forces
every downstream pipeline to choose:

1. **Re-validate the rewrite** on their own data, accepting that
   alignment differences may cascade through downstream analyses.
2. **Maintain a parallel C dependency** so they can verify rust-MAFFT
   isn't introducing drift.

Both are non-starters for adoption. Byte-identity removes the choice:
the Rust binary is a drop-in replacement, no warm-up trust period
required.

## Design decisions that fall out

### No global state

The C code uses ~400 `extern` globals. The Rust modules use owned
`ScoringContext`, `Topology`, `Profile` structs passed explicitly. The
engine is reentrant, embeddable, and trivially thread-safe at the
call-site boundary.

This wasn't done for ergonomic reasons — it was forced by byte-identity.
C MAFFT's globals get assigned in subtle orderings (e.g., `gapfaclocal`
depends on whether `--allowshift` is parsed before or after
`unalignlevel`). Replicating that in Rust required making the
"set-it-now" semantics explicit.

### Hand-ported Cooley-Tukey FFT

`mafft-fft/src/fft_c_compat.rs` is a hand port of MAFFT's `core/fft.c`,
matching the C rounding exactly. Off-the-shelf FFT libraries (rustfft,
realfft, FFTW) vary in butterfly grouping order, producing 1-ULP
correlation differences. On flat-landscape similarity matrices those
ULPs flip FFT-anchor selection — which propagates into different
alignment segmentation and ultimately different output bytes.

### FP order matters

In several hot loops we use `f64::mul_add` (FMA) explicitly to match
clang's `FP_CONTRACT`-on behavior. One concrete case
(`calcW` in `mafft-scoring/src/weighting.rs`):

```rust
// matches clang FP_CONTRACT=ON asm exactly; without this we drift
// 5 of 6 BAliBASE residuals
let s = a.mul_add(b, b.mul_add(c, a * c));
```

Doing this as a naive `a*b + b*c + a*c` produces visibly different
weights on `BB40043`, `BB30018`, `BB40010`, `BB30010`, `BB40004`.

### Tie-break ordering

Tied scores in Needleman-Wunsch / Smith-Waterman traceback are resolved
in C via `if (new >= best)` (note the `>=` — equality favors the new
candidate). Several constraint-mode scanners use a different inequality
(`if (new > best)`). Mixing them up flips one trace direction on one
column, which propagates downstream. We mirror each C usage site
exactly.

### Per-step weight normalization

C MAFFT's `Falign` (the FFT-anchored DP entry) expects the per-group
weight vectors to sum to 1.0. The Rust progressive-alignment caller
normalizes before invocation. Skipping this on a single merge step
produces an off-by-one column on `BB30013`.

## Cross-validation harness

The byte-identity bar is held in place by ~140 FFI tests
(`crates/*/tests/cross_validate_*.rs`) that compile MAFFT's C source
in-tree via `mafft-c-bindings` and call both implementations on the
same inputs, comparing outputs byte-for-byte.

This is documented in detail on the
[Cross-validation page](cross-validation.md).

## What this means for changes

Any change to `mafft-core` / `mafft-align` / `mafft-fft` / `mafft-tree`
/ `mafft-scoring` MUST keep the test suite green. New behaviour goes
behind a feature flag or a new function — not as a tweak to an existing
hot path. If you're tempted to "clean up" a slightly weird-looking
arithmetic in the engine, check the cross-validate test for that
function first; it's probably weird because the C output requires it.

## What this DOESN'T mean

- Algorithmic improvements are still possible — they just need their
  own modes / flags, and the byte-identical existing modes must keep
  working.
- Performance optimisation is fine wherever it doesn't change the
  observable output. We've shaved cycles via `mul_add` collapsing,
  pre-computed boundary tables, and LTO+fat codegen — none of which
  alter the alignment.
- The Python and CLI surfaces can evolve independently of the engine.
  Pretty-printing, progress callbacks, alternate output formats — all
  on the table.
