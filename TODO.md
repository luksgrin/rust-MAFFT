# Pending Work

## Refinement alignment width (729 vs C's 721)

**Status**: 8-column gap remaining. Per-branch weights and all profile/DP components are FFI-validated byte-identical to C. Divergence comes from C's `Falign`-based segmented DP with FFT anchors, which our code does not fully replicate.
**Priority**: Low (1.1% width divergence, both produce valid alignments)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()`

### What was confirmed matching C (FFI cross-validation)

All individually-testable components FFI-validated byte-identical to C:
1. **Per-branch weights** (`weightFromABranch`): all 70 branches of 36-seq topology match to 2.22e-16.
2. **Per-group normalization** (`fastconjuction_noname`): exact match.
3. **Profile construction** (`cpmx_calc_new` + `st_OpeningGapCount` + `st_FinalGapCount` + `gapcountf`): exact match for any weights.
4. **Intergroup scoring** (`intergroup_score`): exact match.
5. **Profile alignment** (`MSalignmm` on stripped profiles): byte-identical at every branch of iter 0 on both the clean FFT-NS-2 input and the evolving state as refinement proceeds.

### Bugs fixed this session

- **MINLEN clamp** in `BranchWeights::new` (`crates/mafft-tree/src/weighting.rs`): C's `checkMinusLength` clamps branch lengths < 0.001. Without this, identical sequences produced 50x-smaller weights. After fix, weights match C to 3.47e-17.
- **Nongap-freq boundary** in `profile_align` (`crates/mafft-align/src/profile.rs`): C pads `gapfreq1pt[lgth1] = 1.0`; we were returning 0 for out-of-bounds. Fixed by explicit `if i < n { prof1.nongap_freq[i] } else { 1.0 }` at the DP boundary. Resolved 4 late-branch divergences that were visible in clean-state `profile_align` vs `MSalignmm` comparison.
- **`compute_split_score`** made sequential (no `par_iter`): deterministic accumulation order.

### Root cause of remaining 8-column gap

C's refinement calls `Falign` which:
1. Runs FFT anchor detection on **full** (unstripped) sequences.
2. Splits the sequences into segments at anchors.
3. Calls `MSalignmm` on each segment of the **full** sequences.

Our code:
1. Strips common-gap columns per group first.
2. Runs FFT anchor detection on full profiles.
3. If anchors found, runs anchored DP on full profiles.
4. If no anchors, runs `profile_align` on **stripped** profiles.

When our `profile_align(stripped)` and C's `MSalignmm(stripped)` are compared directly, they are byte-identical (0/69 divergences). But C's refinement actually calls `MSalignmm(full)` (inside `Falign`), which produces a different alignment width (746 vs our 717 on clean state). The full-vs-stripped alignment difference at specific branches cascades through the 16 iterations × 69 branches.

### How to fix

Port C's `Falign` fully — including its segment generation for the 0-anchor case — so that per-branch DP runs on the full sequences (not stripped). Naive attempts (just switching fallback to `profile_align` on full profiles) produce width explosion (829+) because our DP freely places gap columns when match scores are 0 at all-gap columns. Full Falign-parity requires replicating the segment-boundary logic and any internal constraints that bound gap placement.

## Other diverging modes

- **--bl 80** (BLOSUM80): diverges at retree 1 step 19.
- **--parttree/--dpparttree**: ~940-line diff.
- **--add/--addfragments**: ~900-line diff.

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
