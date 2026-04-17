# Pending Work

## Refinement profile_align width growth (731 vs C's 721)

**Status**: Requires cross-validation harness comparing C's MSalignmm vs Rust cell-by-cell
**Priority**: Medium (FFT-NS-i width 731 vs C's 721 — 1.4% divergence)
**Location**: `crates/mafft-align/src/profile.rs`, `profile_align()` line 465

C's MSalignmm on two same-width non-stripped profiles produces near-zero width growth per call. Our `profile_align` produces ~7 columns growth per call, compounding to +484 across 69 branches (717→1201) when used on non-stripped profiles. Progressive alignment (stripped profiles) is byte-identical to C.

### What was investigated

1. **Comparison operators**: All `>` vs `>=` in the DP inner loop match C exactly (MSalignmm.c lines 841/847/860/866).
2. **Gap cost arrays**: ogcp/fgcp computation matches C's `0.5 * (1.0 - count) * penalty * nongap_freq`.
3. **Scoring matrix**: Same matrix used (substitution_matrix = n_dis + offset).
4. **FFT anchor detection**: Fails on non-stripped profiles in both C and Rust due to gap dilution below segment threshold. C always falls through to single-segment alignment.
5. **`currentw[m]` stale value**: C's `match_calc` fills `[0..lgth2-1]`, leaving `currentw[lgth2]` as the stale value from the buffer swap. Our code zeroes it with `currentw[m] = 0.0`. Removing the zeroing matches C's buffer behavior for alternating rows BUT breaks scores for `tail_gap=true` (400→-1236 on a 4-position test) because the stale value from C's calloc'd init buffer differs from our vec swap pattern. The zeroing does NOT affect progressive alignment (which uses `tail_gap=false`).

### Current workaround

Refinement uses `profile_align` on gap-stripped profiles, keeping width bounded.

### Detailed trace data (from instrumented runs)

On the 36-sequence sample, first refinement iteration (69 branches, 717×717 non-stripped profiles):
- 7 branches produce 0 growth (M=717, I=0, D=0 — perfect diagonal match)
- Remaining branches produce +1 to +94 growth with **symmetric I=D** (paired Insert+Delete)
- The Insert+Delete pairs occur at gap-ambiguous positions where both alternatives score equally
- Growth compounds across branches: 717 → ~1201 after 69 branches

C's MSalignmm on the same 717×717 profiles produces 717 ops (0 growth) for all branches. This means C always chooses Match at positions where our DP chooses Insert+Delete. The tie-breaking divergence is in the DP fill, not the traceback — at gap-rich positions where `match_score ≈ 0` and `gap_cost ≈ 0`, C's accumulation order or rounding produces a value that makes Match win by epsilon, while ours produces a value that makes Insert+Delete win by epsilon.

### Next step

Add MSalignmm to `mafft-sys` FFI bindings and write a cross-validation test that calls both C and Rust on a single 717×717 profile pair from the first refinement branch. Dump `h[i][j]` and `ijp[i][j]` for both and find the first cell where they diverge. The divergence likely occurs at a gap-rich position near the N/C terminus where `nongap_freq ≈ 0`.

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
