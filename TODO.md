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

### C's actual flow (discovered after inspecting dvtditr.c/Falign.c)

dvtditr.c sets `kobetsubunkatsu = 1` (line 54) and `use_fft = 1` (via `-F`). In Falign.c with kobetsubunkatsu=1:
1. The FFT correlation block (lines 1110–1282, guarded by `!kobetsubunkatsu`) is **SKIPPED**. C does **not** compute FFT for refinement.
2. `alignableReagion` runs ONCE at **lag=0** (maxk=1, kouho[0]=0): a deterministic sliding-window scan over the position-wise weighted match score `stra[i] = Σ n_disFFT[k][j]·prf1[k]·prf2[j]` / totaleff.
3. Collected segment centers become cut points: `[0, center_0, center_1, …, len]`.
4. For each segment `[cut1[i], cut1[i+1]]`: extract slice → `commongappick` strips per-group all-gap columns in the slice → `MSalignmm` aligns the stripped slice with boundary gap states `sgap/egap`.
5. Segment outputs are concatenated.

### What was ported (this session)

`realign_all()` in `crates/mafft-core/src/refinement.rs` was rewritten to match C's kobetsubunkatsu=1 flow:
- Compute per-position site scores at lag=0 via `full_prof1.match_score(i, full_prof2, i, substitution_matrix)`, divided by totaleff.
- `alignable_segments` (already ports C's `alignableReagion` sliding window) returns segment centers.
- Cut points assembled from `[0, centers…, width]` (deduped, sorted).
- For each segment slice, per-group `commongappick`-equivalent stripping (`seg_gap1`, `seg_gap2`, `seg_kept1`, `seg_kept2`).
- `profile_align` on the stripped segment with `head_gap=true, tail_gap=true` (matches outgap=1 in dvtditr).
- Per-segment outputs concatenated into `new_sequences`.

The previous FFT-anchor path (`find_fft_anchors` + `align_with_anchors`) is **no longer used** for refinement. C never runs FFT for refinement.

### What remains to match C exactly

- **Boundary gap state (`sgap/egap`) in segment DPs**. C's `MSalignmm_variousdist` switches from `st_OpeningGapCount` to `new_OpeningGapCount(sgappat)` when boundary gap info is supplied. This affects the opening/closing gap costs at segment boundaries. Our `Profile::from_aligned` always acts as `st_OpeningGapCount`. Impact: multi-segment branches only (15/207 ≈ 7% of branches in the 36-seq sample run 3 iterations). For the 0-segment single-segment case, `sgap1` is all `'o'` in C which is equivalent to `st_OpeningGapCount`, so our behavior matches.
- **Verification gap**: the 192/207 branches with 0 segments use effectively the same path as C (single segment, per-group strip, MSalignmm-equivalent DP), and `profile_align(stripped)` is FFI-validated byte-identical to `MSalignmm(stripped)`. Yet the full pipeline still diverges by 8 columns, suggesting the boundary-gap handling on the 15 multi-segment branches is the remaining source — or that `alignableReagion`'s segment output differs from ours (same algorithm but potentially different segment centers due to FP accumulation).

### Next investigation step

FFI-call C's `Falign` (or instrument C's MSalignmm_variousdist branch) for one of the specific divergent branches (e.g., iter 0 step 26 side 1 in our Rust trace, which C marks "identical") and compare segment counts, cut points, and per-segment alignments. The answer will either identify a boundary-gap bug or show that alignableReagion segment centers differ.

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
