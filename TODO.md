# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION,
> not just at convergence. A run that matches C's final width by coincidence
> while producing different intermediate alignments is NOT parity — it means
> our accept trajectory differs from C's, and the match is fragile
> (input-dependent).

## Current parity matrix (36-seq protein sample, mafft-upstream/test/sample)

Verified 2026-05-05 by running `target/release/mafft-rs <args> sample` against
`mafft-upstream/scripts/mafft <args> sample` on the upstream submodule
(MAFFT 7.526) and diffing the FASTA outputs.

| Mode / flags                               | C width | Rust width | diff lines | Status |
|--------------------------------------------|---------|------------|------------|--------|
| FFT-NS-2 (default)                         | 717     | 717        | 0          | byte-exact ✓ |
| NW-NS-2 (`--nofft`)                        | 717     | 717        | 0          | byte-exact ✓ |
| FFT-NS-i (`--maxiterate 100`)              | 721     | 721        | 0          | byte-exact ✓ |
| L-INS-1 (`--localpair --maxiterate 0`)     | 719     | 719        | 0          | byte-exact ✓ (closed 2026-05-05) |
| G-INS-1 (`--globalpair --maxiterate 0`)    | 746     | 746        | 0          | byte-exact ✓ (closed 2026-05-06) |
| E-INS-1 (`--genafpair --maxiterate 0`)     | 729     | 729        | 0          | byte-exact ✓ (closed 2026-05-06) |
| L-INS-i (`--localpair`)                    | 719     | 725        | 934        | divergent — see §3 |
| G-INS-i (`--globalpair`)                   | 746     | 714        | 939        | divergent — see §3 |
| E-INS-i (`--genafpair`)                    | 729     | 725        | 984        | divergent — see §3 |
| BL30 FFT (`--bl 30`)                       | 773     | 773        | 0          | byte-exact ✓ |
| BL45 FFT (`--bl 45`)                       | 729     | 729        | 0          | byte-exact ✓ |
| BL50 FFT (`--bl 50`)                       | 712     | 738        | 961        | divergent — see §4 |
| BL62 FFT (`--bl 62`, default)              | 717     | 717        | 0          | byte-exact ✓ |
| BL80 FFT (`--bl 80`)                       | 712     | 712        | 0          | byte-exact ✓ |
| BL80 NW  (`--bl 80 --nofft`)               | 712     | 712        | 0          | byte-exact ✓ |
| JTT 200 FFT (`--jtt 200`)                  | 729     | 729        | 0          | byte-exact ✓ |
| JTT 100 FFT (`--jtt 100`)                  | 732     | 732        | 4          | content tie-break — see §5 |
| TM 200 NW  (`--tm 200 --nofft`)            | 765     | 765        | 0          | byte-exact ✓ |
| TM 100 NW  (`--tm 100 --nofft`)            | 767     | 767        | 0          | byte-exact ✓ |
| TM 200 FFT (`--tm 200`)                    | 767     | 765        | 148        | divergent — see §5 |
| RNA NW  (`--nofft samplerna`)              | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |

Test suite: 223 Rust tests pass (`cargo test --workspace --exclude pymafft
--release`), 0 failed, 0 ignored. Plus 32 Python tests pass.

---

## §1. G-INS-1 — RESOLVED 2026-05-06

**Mode**: `--globalpair --maxiterate 0` now byte-identical to C.

**Two fixes** combined:

1. **`global_align` was a textbook 3-matrix Needleman-Wunsch**, while
   C's `G__align11` is a max-so-far DP with `>=` tie-break (mirroring
   `L__align11`'s structure but without the local-stop reset).
   Re-ported `global_align` to follow `Galign11.c:913-1474` exactly,
   including the running `mi` / `m_arr` trackers, `>=` for the
   prept-vs-mi/mjpt update, `fpenalty_ex_i` boundary handling at
   `i == lgth1`, and the `Atracking`-style traceback that emits the
   diagonal cell at the SOURCE `(ifi, jfi)` rather than the current
   cell.

2. **`outgap` was hardcoded to 0** (head/tail gap free) for all
   constraint-aware progressive merges. C's `scripts/mafft:2584` does
   NOT pass `$termgapopt = -O` for `--globalpair` (whereas L-INS-i and
   E-INS-i pass it via lines 2593/2601), so G-INS-i's `tbfast` runs
   with `outgap = 1` (head/tail gap penalised). Routed a
   `penalize_term_gaps: bool` flag through
   `progressive_align_with_constraints` → `merge_step_cached` →
   `profile_align_imp(head_gap, tail_gap, …)`. `engine.rs` sets it
   `true` for `GInsi` and `false` for L-INS-i / E-INS-i.

**Regression guard**: `ginsi_maxit0_byte_identical_to_c` in
`crates/mafft-core/tests/end_to_end.rs` (fixture
`tests/fixtures/sample.ginsi.maxit0`).

**FFI guard**: `rust_global_align_matches_c_g__align11` in
`crates/mafft-core/tests/cross_validate_constrained_align.rs` runs
`global_align` and `mafft_sys::G__align11` on the (M63632, K03494)
pair (which has a different N-terminal that triggers the head-gap
tie-break) and asserts byte-exact equality.

---

## §2. E-INS-1 — RESOLVED 2026-05-06

**Mode**: `--genafpair --maxiterate 0` now byte-identical to C.

**Three combined fixes**:

1. **Re-ported `genaffine_local_align`** (`crates/mafft-align/src/genaffine.rs`)
   to mirror C's `genL__align11` exactly (`genalign11.c:113-660`). The
   previous Rust implementation was a textbook 3-matrix DP; C uses the
   same max-so-far DP scheme as `L__align11` plus a separate "skip" gap
   state (running max `Mi` / per-column `largeM[j]`, fed to row-local
   `tbk`) with `penalty_OP` open and zero extension. Traceback uses two
   absolute-coordinate arrays `ijpi[i][j]` / `ijpj[i][j]` so a single
   step can change both i and j arbitrarily.

2. **Added `PairAligner::GeneralizedAffine`** in
   `crates/mafft-align/src/constraints.rs`. `build_homology_table` now
   takes an `op_penalty: f64` arg used only for this variant.
   `engine.rs` routes `AlignmentMode::EInsi` through it.

3. **Mirrored C's E-INS-i parameter overrides** (`scripts/mafft:1940-1948`):
   when `distance == "localgenaf"` the script resets `lexp = "0.0"` and
   `laof = "0.0"` so the regular gap-extension and matrix-offset are
   zeroed out, leaving only the skip-gap (`LGOP = -6.00`) as the
   long-range penalty. Without this our `genL__align11` was being driven
   with non-zero `penalty_ex` / `offset` and produced a different
   (wider) alignment than C.

**Regression guard**: `einsi_maxit0_byte_identical_to_c` in
`crates/mafft-core/tests/end_to_end.rs` (fixture
`tests/fixtures/sample.einsi.maxit0`).

**FFI guard**: `rust_genaffine_align_matches_c_gen_l__align11` in
`crates/mafft-core/tests/cross_validate_constrained_align.rs` runs
`genaffine_local_align` and `mafft_sys::genL__align11` on the
(M63632, U22180) pair and asserts byte-exact equality of score, width,
offsets, and aligned strings.

---

## §3. INS-i modes with iterative refinement — PARTIALLY RESOLVED

**Status as of 2026-05-06**:

| Mode    | C width | Rust width | Diff lines |
|---------|---------|------------|------------|
| L-INS-i | 719     | 725        | 934        |
| G-INS-i | 746     | 727        | 1000       |
| E-INS-i | 729     | 730        | 984        |

Small inputs match byte-exact:
- L-INS-i n=9 iter=2 ✓
- L-INS-i n=12 iter=5 ✓
- L-INS-i n=30 iter=1 ✓

Cascade divergence at n ∈ {13, 14, 15, 20, 36}.

**Three fixes landed 2026-05-06**:

1. **Accept/reject score now mirrors C's `mscore = oimpmatchdouble +
   tmpdouble`** (`tditeration.c:953,1094`). Added
   `compute_impmatch_diagonal` in `crates/mafft-core/src/refinement.rs`.
   For an existing alignment of width W, it builds the impmtx via
   `build_imp_matrix(..., width, width, ...)` (matching C's
   `part_imp_match_init_strict(..., length, length, ...)`) and sums
   `impmtx[i][i]` across the diagonal. `iterative_refine` combines this
   with the intergroup substitution score when computing both
   `old_score` and `tscore` for the accept threshold check.

2. **`Falign_localhom` port** (`Falign_localhom.c:163`,
   kobetsubunkatsu=1 path) in `realign_all_constrained_fft`. Mirrors
   the FFT-segmented loop of the unconstrained refinement but per-
   segment calls `profile_align_imp` with a local impmtx slice
   `local_imp[i][j] = global_imp[a + gapmap1[i]][a + gapmap2[j]]`
   where `a` = segment's parent-column start and `gapmap1`/`gapmap2`
   are the per-group `commongappick_record` mappings. This mirrors C's
   `part_imp_match_out_vead_gapmap`
   (`partSalignmm.c:71-83`).

3. **`partA__align` strict-`>` tie-break** for prept-vs-mi/mjpt update
   (partSalignmm.c:1218,1235; commented "// 2018/Apr"). C's progressive
   `A__align` uses `>=` (Salignmm.c:1926,1946). Added
   `profile_align_imp_with_tiebreak(..., strict_part_tiebreak: bool)`.
   The Falign_localhom path passes `true` (matches partA__align); the
   progressive constraint path keeps `false` (matches A__align).

**Regression guards**:
- `linsi_first9_iter2_byte_identical_to_c` — n=9, --maxiterate 2.
- `linsi_first12_iter5_byte_identical_to_c` — n=12, --maxiterate 5
  (exercises Falign_localhom's segment loop).

**Remaining gap (n ≥ 13 cascading divergence)**:

The Falign_localhom port and tie-break fix handle the structural
differences between progressive and refinement DP, but the 36-seq
output still differs from C by ~6-19 columns. Profiling shows:

- The first divergence at n=13 is a single-residue gap-shift in the
  tail (`EVS-T-S` vs `EVS--TS`) — a tie-break in either the DP or
  the score computation.
- For n=14+ the divergence cascades because each branch's accept/reject
  decision depends on the previous branch's output.

**Concrete next task** (effort: ~1 day): instrument C's
`partA__align` to dump per-cell DP state for the n=13 divergent
branch, compare against our `profile_align_imp_with_tiebreak` cell-by-
cell. The first cell where our `wm` / `mi` / `m[j]` / `ijp[i][j]`
differs from C's is the bug. Likely candidates:
- The gap-frequency multiplier `gf1va` / `gf1vapre` derivation for
  segment-stripped profiles (we may compute it from stripped seqs;
  C may compute from full + adjust via headgapfreq).
- The boundary `sgap[k]` / `egap[k]` per-sequence-status correction
  (we apply ogcp/fgcp tweaks only when neighbour columns are gap;
  C's `getkyokaigap` may behave subtly differently).

---

## §4. BLOSUM50 (`--bl 50`) — FFT path divergent despite cell-equal matrix

**Observed**: Rust width 738 vs C 712, diff 961 lines.

**Verified equal**:
- `cross_validate_blosum50_n_dis_cell_by_cell` (passes 0 mismatches)
- `cross_validate_blosum50_n_dis_fft_cell_by_cell` (passes 0 mismatches)
- BL50 raw 210-cell lower-triangle table matches `tmpmtx50` exactly
- Both Rust and C use the same gap penalty (`-1.53`), same `-h 0`, same
  `disttbfast` invocation flags (only `-b 50` vs `-b 62` differs)

**Pinpointed divergence (2026-04-28 instrumentation)**:

| Step | Clusters | C width / score | Rust width / score |
|------|----------|------------------|---------------------|
| 0–32 | various  | identical        | identical           |
| 33   | 13 × 8   | 624 / 47967.6    | 641 / 47954.6       |
| 34   | 21 × 15  | 721 / 56627.1    | 738 / 56582.8       |

The first 33 merges are identical in width AND score. Step 33 then
diverges by 17 columns and ~13 score points (0.027% — strongly
indicative of an FFT-anchor tie-break, not a scoring bug). Same input
profiles, same matrix, same gap penalty: our `fft_profile_align` picks
a slightly different anchor placement and the DP collapses into a wider
alignment.

**Diagnosis** (from prior investigation): C's `Falign`
(`Falign.c:1307-1372`) accumulates segments from all `NKOUHO=20`
candidate lags into a single flat list, sorts by start position
independently in seq1 and seq2, builds a sparse `crossscore[i][j]`
matrix where each segment writes its score at its `(rank-in-sort1+1,
rank-in-sort2+1)` cell, then runs `blockAlign2` over that combined
matrix. Our `find_fft_anchors` instead selects the single best lag and
runs `block_align` on its segments only.

For BL62 and 7 other modes the best lag dominates so the two
behaviours coincide; for BL50, JTT 100 (FFT), and TM 200 (FFT) multiple
lags carry comparable per-segment scores and C picks a different
optimal subset.

**Two prior structural rewrite attempts** (flat anchor list across all
candidates, fed to combined block_align) reached width 711 (vs 712) for
BL50 but broke BL62 (144-line diff). Reverted; partial improvements
landed in `block_align.rs` (permit-zero gating, traceback dedupe).

**Concrete next task** (effort: 3-4 h): instrument C's `Falign` to dump
the flat `(c1, c2, score, lag)` anchor-pair list pre-`blockAlign2` and
the `crossscore[][]` matrix pre-DP for the BL62 step-0 and BL50 step-33
merges. Diff against our equivalent dumps. The first divergence
identifies whether the issue is (a) sort tie-break, (b) corner
sentinel value, (c) reverse-mapping ambiguity, or (d) lag exclusion
bounds.

---

## §5. FFT tie-break residuals on `--jtt 100` and `--tm * (FFT)`

**Observed**:
- `--jtt 100`: width matches (732 / 732), 4 lines diff (one residue
  shifted by one column in seq 12: `MAA-W` vs `MA-AW`)
- `--tm 200 --nofft`: byte-exact ✓
- `--tm 100 --nofft`: byte-exact ✓
- `--tm 200` (FFT): width 765 vs 767, 148-line diff

Same class as §4 (BL50 FFT). Cell-by-cell matrix matches via
`cross_validate_jtt100_n_dis`, `cross_validate_tm_n_dis`,
`cross_validate_tm_n_dis_fft`. Float-precision tie-break in the FFT
correlation peak or anchor-selection sort.

**Concrete next task**: solve §4 first. The fix that closes BL50 will
likely close JTT 100 and TM 200 (FFT) too — same root cause.

---

## §6. PartTree (`--parttree`, `--dpparttree`)

**Observed**: `--parttree --nofft` has ~940-line diff from C.

**Priority**: Medium (separate algorithm, independent of refinement).

**Location**: `crates/mafft-tree/src/parttree.rs`.

The current port is functional but has a subtle bug in the
subtree-grouping / seed-selection logic. Needs a full pass against
C's `partSeedSelect`, `partTreeGrow`, and `getsuboptimal`-style
pruning paths in `mltaln9.c`.

**Concrete next task**: port C's parttree pipeline exhaustively, with
FFI-based cross-validation at each stage (seed set, subtree assignment,
distance matrix). Effort: 1-2 days.

---

## §7. Add / AddFragments (`--add`, `--addfragments`, `--keeplength`)

**Observed**: `--add --nofft` has ~900-line diff from C.

**Priority**: Medium-high (this feature also blocks the per-group-strip
perf fix in progressive alignment — see §8).

**Location**: `crates/mafft-core/src/add.rs`.

C's `addsingle.c` + `insertnewgaps()` from `addfunctions.c` (lines
445-673) is only partially ported. The `gaplen[]` / `gapmap[]` /
`posin12`-counter machinery is missing.

**Concrete next task**: port `insertnewgaps()` fully:
1. After profile DP, reconstruct group1/group2 aligned sequences with
   `=` markers at new-gap positions.
2. Port `findnewgaps()` → `gaplen[]`.
3. Port `adjustgapmap()` → `gapmap[]`.
4. Port the main loop of `insertnewgaps()` using `j` and `posin12`
   counters.

Effort: 2-3 days.

---

## §8. Per-group gap stripping in progressive alignment (perf, not correctness)

**Status**: Not implemented — we strip columns all-gap across all
sequences globally, not per-group like C's `commongappick()`.

**Priority**: Low for correctness (output matches C); medium for
performance.

**Location**: `crates/mafft-core/src/progressive.rs`, `merge_step()`.

For a 500-column MSA where group1 has 100 residue-containing columns
and group2 has 120:
- C's per-group: DP aligns 100 × 120 = 12,000 cells.
- Our global: DP aligns up to 500 × 500 = 250,000 cells.

**Why it's hard**: per-group stripping breaks when `kept1 != kept2`
because "other" sequences (not in either group) can't follow both
cursors simultaneously. Six prior approaches all failed.

**Concrete next task**: solve §7 first (port `insertnewgaps()`). Once
that exists, the progressive alignment can use it to reinsert stripped
columns correctly. SAME port unlocks both.

---

## §9. Other unimplemented / unvalidated modes

- RNA-aware alignment (`--xinsi`, `--qinsi`): wired to external tools
  (`mxscarnamod`, `contrafold`) but not validated end-to-end against C.
- `--allowshift`: gap-shift/warp DP path, partial implementation in
  `crates/mafft-align/src/shift.rs`.

**Concrete next task**: for each, run `mafft --xinsi sample` /
`mafft --qinsi sample` / `mafft --allowshift --globalpair sample`
against our binary and capture the first divergence. Effort scales
with the number of paths each mode activates.

---

## Recommended order of attack

1. **§3 (INS-i refinement loop)** — partial fix landed 2026-05-06 (added
   impmatch to accept/reject score; small inputs n≤12 byte-exact for
   L-INS-i). Closing fully needs a port of C's `Falign_localhom` (FFT-
   segmented DP with constraints) which is its own substantial chunk.
2. **§4 (BL50 FFT)** — structural FFT change; once this lands, §5
   likely closes for free.
3. **§7 (--add)** — unblocks §8 (perf) and is a real feature gap.
4. **§6 (parttree)** — algorithmic port, substantial work.
5. **§9 (RNA / allowshift)** — as needed.
