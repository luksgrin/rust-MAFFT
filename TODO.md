# Pending Work

## FFT-accelerated alignment in iterative refinement

**Status**: Not implemented (using plain `profile_align` instead of `Falign`)
**Priority**: High (last source of divergence in the refinement loop)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()`

C's `TreeDependentIteration()` calls `Falign()` (FFT-accelerated profile alignment with anchor detection) for re-alignment during refinement when FFT is enabled (`tditeration.c` line ~2153). Our refinement loop always uses plain `profile_align()`. This means FFT-NS-i diverges from C; NW-NS-i (which doesn't use FFT) is unaffected.

### How to fix

Wire `fft_profile_align()` (from `crates/mafft-align/src/fft_align.rs`) into `realign_all()` when `RefinementParams::use_fft` is true. The function already exists and is used in progressive alignment — it just needs to be called from the refinement path as well:

1. Thread the `use_fft` flag and scoring context's FFT-related fields into `realign_all()`.
2. When `use_fft` is true, call `fft_profile_align()` instead of `profile_align()`.
3. C also calls `commongappick_record()` (which records the gap map) before FFT alignment but skips it for non-FFT — match this behavior.

Note: `RefinementParams` already has a `use_fft` field, and the engine already passes the correct value. Only the `realign_all()` function needs updating.

## Per-group gap stripping in progressive alignment

**Status**: Not implemented (using global-all-gap stripping instead)
**Priority**: Medium (performance optimization — does not affect correctness)
**Location**: `crates/mafft-core/src/progressive.rs`, `merge_step()`

C's `commongappick()` strips columns that are all-gap within each group **independently** before profile alignment. Our code strips only columns that are all-gap across **all sequences globally**. The alignment result is the same, but the DP matrix is larger than necessary.

For a 500-column MSA where group1 has 100 residue-containing columns and group2 has 120:
- C's per-group: DP aligns 100 × 120 = 12,000 cells
- Our global: DP aligns up to 500 × 500 = 250,000 cells

### Why it's hard

Per-group stripping was attempted but fails on real data. The bug is in re-inserting stripped columns after alignment when the two groups have different kept-column sets. Six approaches were tried (two-pass interleaving, column classification, cursor models, template approaches) — all break when `kept1 != kept2` because "other" sequences can't follow both cursors simultaneously.

### How to fix

Port C's `insertnewgaps()` from `addfunctions.c` (lines 445-673). C's approach bypasses the cursor problem entirely: it operates on result character arrays with `=` markers for new gaps, using `gaplen[]`/`gapmap[]` arrays and independent per-group position counters (`posin12` for group1/group2, `j` for others). Steps:

1. After profile DP, reconstruct group1/group2 aligned sequences with `=` markers at new-gap positions.
2. Port `findnewgaps()` → `gaplen[]`.
3. Port `adjustgapmap()` → `gapmap[]`.
4. Port the main loop of `insertnewgaps()` using `j` and `posin12` counters.

The same optimization could also be applied to the refinement loop (`crates/mafft-core/src/refinement.rs`), though refinement currently avoids the problem by splitting ALL sequences into two groups (no "other" sequences exist).
