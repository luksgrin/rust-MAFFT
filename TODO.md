# Pending Work

## Refinement alignment width (731 vs C's 721)

**Status**: Both C and Rust strip gaps and use profile_align/MSalignmm identically. Divergence is compounded epsilon-level float differences across 16 iterations × 69 branches.
**Priority**: Low (1.4% width divergence, both produce valid alignments)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()`

### What was confirmed matching C

Investigation of C's Falign.c revealed that **C DOES strip gap columns before calling MSalignmm** during refinement:
- `kobetsubunkatsu = 1` and `fftkeika = 1` in dvtditr.c (lines 51, 54)
- Falign.c line 1610: `if( kobetsubunkatsu && fftkeika ) commongappick( clus1, tmpres1 );`
- This strips per-group all-gap columns BEFORE MSalignmm, matching our `group_all_gap_columns` + stripped `profile_align`

Additionally confirmed matching:
- DP comparison operators (`>` for Insert/Delete vs wm, `>=` for running-max updates) — MSalignmm.c lines 841/847/860/866
- Gap cost formulas (`ogcp`/`fgcp`): `0.5 * (1.0 - count) * penalty * nongap_freq`
- Profile construction (`cpmx_calc_new` / `Profile::from_aligned`): same weighted accumulation
- Gap counting (`gapcountf` / `st_OpeningGapCount` / `st_FinalGapCount`): same logic
- Boundary gap handling (`sgap`/`egap` all `'o'` for single-segment case) equivalent to `st_*` functions
- `legacygapcost = 0` path with `headgapfreq = 1.0` matches our `hgf = 1.0`
- FFT segment detection fails on non-stripped profiles in both C and Rust (gap dilution)
- `currentw[m]` stale-value difference confirmed — does NOT affect the stripped-profile path

### Root cause: per-branch sequence weighting

C uses **per-branch weights** (`weightFromABranch` in treeOperation.c, weight=4 mode) that compute a different weight vector for each of the 69 branches per iteration. These weights are derived from the tree structure around each specific branch point, using synthetic branch lengths (harmonic mean recursion) and a 3-way branch weight formula. Our code uses **global weights** (`sequence_weights`, porting C's `counteff_simple_double`) that are the same for all branches.

The per-branch weighting changes which alignment decision is optimal at each branch, because profiles built with different weights produce different match scores and gap penalties. This produces the 14-column divergence visible from the very first refinement iteration (735 vs C's 721 after 1 iteration).

### Implementation status

`BranchWeights` in `crates/mafft-tree/src/weighting.rs` implements the framework (unrooted tree construction, synthetic length computation, recursive weight propagation), but the tree conversion from `Topology` to the unrooted Node structure has bugs — produced width 785 (worse than 731) when tested. The infrastructure is in place but disabled (refinement uses `_branch_weights` unused, falls back to `global_weights`).

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
