# Pending Work

## Refinement no-anchor fallback (width 731 vs C's 721)

**Status**: Partially addressed — anchored path matches C; no-anchor fallback diverges
**Priority**: Low (FFT-NS-i width 731 vs C's 721 on the sample dataset)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()`

When FFT finds anchors, both C and Rust run anchored DP on non-stripped profiles — this path matches. The divergence comes from branches where FFT finds **no anchors** (asymmetric splits like 1-vs-35):
- C creates a single synthetic segment spanning the full sequences and aligns it (bounded by `alloclen = nlenmax * 9`).
- Rust falls back to `profile_align` on gap-stripped profiles (bounded by residue count).

### How to close the remaining gap

Port C's single-segment fallback from Falign.c (lines 1377-1421): when count=0, set `cut1[0]=0, cut2[0]=0, cut1[1]=len1, cut2[1]=len2` and align the full sequences as one segment. This requires an `alloclen`-style width cap to prevent explosion. Could be implemented as `align_with_anchors` on non-stripped profiles with no anchors (equivalent to one big segment), plus a max-width check that rejects the result if `ops.len() > alloclen` and falls back to stripped alignment.

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
