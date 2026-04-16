# Pending Work

## Refinement profile_align width growth (731 vs C's 721)

**Status**: Diagnosed — profile_align DP numerical divergence on non-stripped profiles
**Priority**: Low (FFT-NS-i width 731 vs C's 721 on the sample dataset — 1.4% divergence)
**Location**: `crates/mafft-align/src/profile.rs`, `profile_align()`

C's `MSalignmm` on two same-width non-stripped profiles (from the same alignment) produces near-zero width growth per call. Our `profile_align` produces ~7 columns growth per call, compounding to +484 across 69 branches in one iteration (717 → 1201). The profiles, gap penalties (ogcp/fgcp), scoring matrix, and DP algorithm are functionally identical. The divergence is a subtle numerical difference (likely in floating-point accumulation order or tie-breaking in the traceback) that only manifests on non-stripped profiles with many gap columns. Progressive alignment (which uses stripped profiles) is byte-identical to C.

### Current workaround

The no-anchor fallback (all refinement branches, since FFT segment detection fails on gap-diluted non-stripped profiles) uses `profile_align` on gap-stripped profiles, keeping width bounded by residue count.

### How to close the remaining gap

Investigate the per-position DP scores and traceback decisions for a single refinement branch, comparing C's MSalignmm output against ours. The divergence likely originates at gap-rich positions where the tie-breaking between Match and Insert/Delete differs. Adding per-cell RDBG output to both C and Rust would isolate the exact position and state.

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
