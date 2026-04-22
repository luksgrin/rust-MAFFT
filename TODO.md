# Pending Work

## Refinement alignment width (720 vs C's 721)

**Status**: 1-column gap remaining (0.1%). Full Falign pipeline (alignableReagion + per-segment commongappick + profile_align) is FFI-validated byte-identical to C at every branch of iter 0, both on clean input and as state evolves.
**Priority**: Very low (0.1% width divergence, likely FP drift in accept decisions)
**Location**: `crates/mafft-core/src/refinement.rs`, `realign_all()` and `iterative_refine()`

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
- **`realign_all` rewritten to match C's `Falign` with `kobetsubunkatsu=1`**: dvtditr skips the FFT block and runs `alignableReagion` once at lag=0 for segment detection. Ported: site scores via `match_score(i,i)/totaleff`, `alignable_segments` for cut points, per-segment `commongappick`-equivalent stripping, per-segment `profile_align`.
- **Boundary gap state (`sgap`/`egap`)**: For multi-segment branches, adjust each segment's `ogcp[0]` and `fgcp[length-1]` to match C's `new_OpeningGapCount`/`new_FinalGapCount(sgappat)` semantics (a sequence already in a gap just before the segment does NOT count as "opening").
- **Identity check uses only 2 representative sequences** (tditeration.c:2184-2185): C picks `s1 = memlist1[0], s2 = memlist2[0]` and does `strcmp` on those two aligned strings only. Our code was checking **all** sequences, which flagged column rearrangements preserving the representatives as "changed" and triggered spurious accepts. **This change alone closed 9 of the 10 columns of divergence (729 → 720).**

### FFI validation tests added

- `compare_alignable_regions_with_c_iter0`: 0/69 divergences between our `alignable_segments` and C's `alignableReagion` for all iter-0 branches.
- `compare_segmented_falign_with_c_iter0`: 0/69 divergences between our segmented pipeline and C's `Falign` on clean state.
- `evolving_state_rust_vs_c_falign_iter0`: 0/69 divergences as state evolves through iter 0 accepts.

### C's actual flow (discovered after inspecting dvtditr.c/Falign.c)

dvtditr.c sets `kobetsubunkatsu = 1` (line 54) and `use_fft = 1` (via `-F`). In Falign.c with kobetsubunkatsu=1:
1. The FFT correlation block (lines 1110–1282, guarded by `!kobetsubunkatsu`) is **SKIPPED**. C does **not** compute FFT for refinement.
2. `alignableReagion` runs ONCE at **lag=0** (maxk=1, kouho[0]=0): a deterministic sliding-window scan over the position-wise weighted match score `stra[i] = Σ n_disFFT[k][j]·prf1[k]·prf2[j]` / totaleff.
3. Collected segment centers become cut points: `[0, center_0, center_1, …, len]`.
4. For each segment `[cut1[i], cut1[i+1]]`: extract slice → `commongappick` strips per-group all-gap columns in the slice → `MSalignmm` aligns the stripped slice with boundary gap states `sgap/egap` (enables `new_OpeningGapCount(sgappat)` / `new_FinalGapCount(egappat)`).
5. Segment outputs are concatenated, producing `aseq`.
6. Identity check: `strcmp(aseq[s1], bseq[s1]) == 0 && strcmp(aseq[s2], bseq[s2]) == 0` where s1 = memlist1[0], s2 = memlist2[0]. If identical, skip tscore and don't accept.
7. Otherwise compute `tscore = intergroup_score(aseq, ...)`. Accept if `tscore > mscore - cut/100*mscore`.

### What was ported (this session)

All steps 1–7 above are ported in `realign_all()` (steps 1–5) and `iterative_refine()` (steps 6–7) in `crates/mafft-core/src/refinement.rs`. FFI cross-validation tests confirm:

- `alignable_segments` = C's `alignableReagion` at lag=0, byte-identical segment centers across all iter-0 branches.
- Per-segment stripping + `profile_align` = C's `Falign(MSalignmm per segment)`, byte-identical outputs across all iter-0 branches (both clean and evolving state).
- Per-branch weights = C's `weightFromABranch` (exact to 2.22e-16).
- Per-group weight normalization = C's `fastconjuction_noname`.
- Intergroup scoring = C's `intergroup_score`.

The previous FFT-anchor path (`find_fft_anchors` + `align_with_anchors`) is **no longer used** for refinement. C never runs FFT for refinement.

### Remaining 1-column gap

Rust iter=4+ converges at width 720; C iter=4+ converges at 721. Every per-branch component matches C byte-identically via FFI — the 1-column gap comes from cumulative differences in accept decisions over 16 iterations × 69 branches. Likely sources:

- **Floating-point drift in `intergroup_score`**: our `compute_split_score` is FFI-validated against C's `intergroup_score` (exact match), but the accept boundary `tscore > mscore` is sensitive to 1-ULP differences that can propagate through 1104 branches.
- **Subtle non-determinism**: any remaining `par_iter` or HashMap-ordered accumulation — unlikely since `compute_split_score` was made sequential in this session.
- **Boundary gap correction completeness**: the `sgap`/`egap` corrections are applied but may not fully replicate C's `new_OpeningGapCount` for all edge cases (e.g., when stripping removes columns from the segment start).

### Next investigation step (if pursuing)

Add branch-by-branch diagnostic comparison between Rust's real refinement and C's MSalignmm-variousdist outputs across ALL iterations. The first branch where accept decisions diverge identifies whether the residual is FP drift in `tscore` or an algorithmic subtlety. Since the gap is only 0.1% and both outputs are valid alignments, this may not warrant further investment.

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
