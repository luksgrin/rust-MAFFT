# Pending Work

## Non-stripped FFT alignment in refinement

**Status**: Partially addressed (FFT on gap-stripped profiles; C uses FFT on non-stripped profiles)
**Priority**: Low (FFT-NS-i width is 731 vs C's 721 on the sample dataset — close but not identical)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()`

C's Falign in refinement operates on non-stripped character arrays in a fixed-size buffer (`alloclen`), mutating sequences in-place. Our `fft_profile_align` builds new `Vec<u8>` from alignment operations (functional style). Without gap stripping, this causes unbounded width explosion: on non-stripped 717-column profiles, width grew to 26,454+ in a single iteration because the FFT's `profile_align` fallback (when anchors are poor) can produce width up to `prof1.len + prof2.len`.

Current state: `fft_profile_align` is called on gap-stripped profiles. This avoids the explosion and uses the correct algorithm (FFT correlation → anchors → segmented DP), but the different input (stripped vs non-stripped columns) produces different anchor positions.

### How to close the remaining gap

Port C's in-place mutation model: operate on `&mut [u8]` slices within a capacity-bounded buffer, overwriting sequences in-place. This requires restructuring `realign_all` to work with mutable character arrays instead of building new Vecs from operations. The `alloclen` bound prevents width explosion naturally.

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
