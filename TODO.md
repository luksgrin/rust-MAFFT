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
| L-INS-i (`--localpair`)                    | 731     | 731        | 0          | byte-exact ✓ (closed 2026-05-07) |
| G-INS-i (`--globalpair`)                   | 746     | 746        | 0          | byte-exact ✓ (closed 2026-05-07) |
| E-INS-i (`--genafpair`)                    | 729     | 729        | 0          | byte-exact ✓ (closed 2026-05-07) |
| BL30 FFT (`--bl 30`)                       | 773     | 773        | 0          | byte-exact ✓ |
| BL45 FFT (`--bl 45`)                       | 729     | 729        | 0          | byte-exact ✓ |
| BL50 FFT (`--bl 50`)                       | 712     | 712        | 0          | byte-exact ✓ (closed 2026-05-08) |
| BL62 FFT (`--bl 62`, default)              | 717     | 717        | 0          | byte-exact ✓ |
| BL80 FFT (`--bl 80`)                       | 712     | 712        | 0          | byte-exact ✓ |
| BL80 NW  (`--bl 80 --nofft`)               | 712     | 712        | 0          | byte-exact ✓ |
| JTT 200 FFT (`--jtt 200`)                  | 729     | 729        | 0          | byte-exact ✓ |
| JTT 100 FFT (`--jtt 100`)                  | 732     | 732        | 0          | byte-exact ✓ (closed 2026-05-09 by §4) |
| TM 200 NW  (`--tm 200 --nofft`)            | 765     | 765        | 0          | byte-exact ✓ |
| TM 100 NW  (`--tm 100 --nofft`)            | 767     | 767        | 0          | byte-exact ✓ |
| TM 100 FFT (`--tm 100`)                    | 767     | 767        | 0          | byte-exact ✓ (closed 2026-05-09 by §4) |
| TM 200 FFT (`--tm 200`)                    | 767     | 767        | 0          | byte-exact ✓ (closed 2026-05-10) |
| RNA NW  (`--nofft samplerna`)              | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |

Test suite: 250 Rust tests pass (`cargo test --workspace --exclude pymafft
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

## §3. L/G/E-INS-i with iterative refinement — RESOLVED 2026-05-07

**Status as of 2026-05-07** (hat2-distance fix + naivepairscore11 fix
fully close the L/G/E-INS-i cascades):

| Mode    | C width | Rust width | Diff lines |
|---------|---------|------------|------------|
| L-INS-i | 719     | 725        | 934        |
| G-INS-i | 746     | 727        | 1000       |
| E-INS-i | 729     | 730        | 984        |

L-INS-i is fully byte-identical to C across all tested input sizes:
- L-INS-i n=9 iter=2 ✓
- L-INS-i n=12 iter=5 ✓
- L-INS-i n=13 iter=2 ✓
- L-INS-i n=14 iter=1, iter=2 ✓ *(closed by 2026-05-07 hat2-distance fix)*
- L-INS-i n=15 iter=2 ✓ *(closed by hat2-distance fix)*
- L-INS-i n=30 iter=1 ✓
- L-INS-i n=36 iter=2 ✓ *(closed by hat2-distance fix)*
- maxit0 (progressive build) byte-identical for all sizes

**Six fixes landed 2026-05-06 / 2026-05-07**:

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

4. **Boundary nongap-freq threading** (partSalignmm.c:1058-1075). C's
   `partA__align` multiplies the row-0 / column-0 boundary corrections
   by `headgapfreq{1,2}` (sgap-derived nongap fractions of the column
   *just before* the segment in the full alignment) and uses
   `gapfreq{1,2}[lgth]` (egap-derived nongap fractions of the column
   *just after* the segment) at the inner DP's `i==lgth1` / `j==lgth2`
   boundary cells. Rust previously hardcoded 1.0 for these. Added
   `BoundaryFreqs { head1, head2, tail1, tail2 }` and
   `profile_align_imp_with_boundary`; `realign_all_constrained_fft`
   computes the four values from `sequences[idx][a-1]` / `[b]` weighted
   by `w1n` / `w2n`. Closes the n=13 divergence.

5. **`hat2`-rounded initial pairwise distances for refinement tree**
   (`scripts/mafft:2693`, `dvtditr.c:751`). The script invokes
   `dvtditr` *without* a preceding `dndpre` call for `--localpair`
   / `--globalpair` / `--genafpair`, so dvtditr reads the `hat2`
   file written by the initial pairlocalalign step — initial
   pairwise distances at 3-decimal precision (`%.3f`). Rust was
   recomputing distances from the progressive alignment via
   `compute_distance_matrix_scoring`, producing different edge
   lengths and hence different `weightFromABranch` weights, which
   tipped the accept/reject scoring at a few branches and
   cascade-diverged the iter≥1 alignment.

   `crates/mafft-core/src/engine.rs` now stashes the initial
   pairwise distance matrix from `build_homology_table` and, for
   modes that ran pairlocalalign (constrained), feeds it through
   `musclesupg` for the refinement tree after rounding each cell to
   `(d * 1000).round() / 1000`. For FFT-NS-i (no pairlocalalign;
   distance="ktuples"), the original `dndpre`-style recomputation
   path is preserved.

**Regression guards**:
- `linsi_first9_iter2_byte_identical_to_c` — n=9, --maxiterate 2.
- `linsi_first12_iter5_byte_identical_to_c` — n=12, --maxiterate 5
  (exercises Falign_localhom's segment loop).
- `linsi_first13_iter2_byte_identical_to_c` — n=13, --maxiterate 2
  (guards the boundary-freq fix).
- `linsi_first14_iter1_byte_identical_to_c` — pins the very first
  refinement iteration (the cleanest signal during the
  hat2-distance investigation).
- `linsi_first14_iter2_byte_identical_to_c`,
  `linsi_first15_iter2_byte_identical_to_c`,
  `linsi_first36_iter2_byte_identical_to_c` — guard the
  hat2-distance fix at sizes the prior cascade affected.
- `linsi_first14_maxit0_byte_identical_to_c`,
  `linsi_first15_maxit0_byte_identical_to_c`,
  `linsi_first36_maxit0_byte_identical_to_c` — pin progressive
  build byte-identity at the same sizes.

6. **`naivepairscore11` for E-INS-i pairwise distance** (closed 2026-05-07):
   For `--genafpair`, C's script passes `-N -Z` (per `scripts/mafft:1946`,
   line 2601) which sets `usenaivescoreinsteadofalignmentscore=1` in
   pairlocalalign. The relevant branch
   (`pairlocalalign.c:2225-2229`) runs `genL__align11` only to derive
   the alignment, then OVERRIDES `pscore` with `naivepairscore11(seq1,
   seq2, 0.0)` — a naive sum of substitution scores at aligned non-gap
   columns, with gap penalty zero — and feeds that into
   `score2dist`. Rust was using the genaffine score itself, which
   gave systematically larger distances (e.g. d[0][1] = 0.280 vs C's
   0.269 for the n=15 sample) and hence a different refinement tree.
   `crates/mafft-align/src/constraints.rs::build_homology_table` now
   computes `score_for_dist` separately for the
   `PairAligner::GeneralizedAffine` path: walk the alignment, add
   `matrix[i][j]` for each non-gap column, and use that for the
   distance conversion. The original `alignment.score` is still used
   for the homology table's region `opt`, mirroring C's split.

   Closes the n=15..n=36 E-INS-i divergences. Guarded by
   `einsi_first14_iter2_byte_identical_to_c`,
   `einsi_first36_iter2_byte_identical_to_c`.

**Investigation history (preserved for context)**:

The cascade was tracked down via per-branch dumps in both Rust
and C (`RUST_DUMP_BRANCHES` / `C_DUMP_BRANCHES`). Group splits
matched perfectly between Rust and C, but the per-branch
`weightFromABranch` weights diverged. Tracing through showed
the topology shape was correct but edge lengths differed —
specifically `step 0 ll=0.009` (C) vs `0.022` (Rust) for the
n=14 sample. That pointed to the distance matrix feeding
`musclesupg`. Dumping `eff` after `readhat2_pointer` in
`dvtditr.c:751-753` confirmed C uses a 3-decimal-rounded
distance matrix from the `hat2` file, while Rust was
recomputing distances from the progressive alignment.

---

## §4. BLOSUM50 (`--bl 50`) — RESOLVED 2026-05-08

**Status as of 2026-05-08**: byte-identical to C MAFFT 7.526. Width
712 = 712, 0-line diff. Regression test
`fftns2_bl50_byte_identical_to_c` added.

### Original symptom (pre-fix)

**Status as of 2026-05-08**: Rust width 738 vs C 712, diff ~74 lines.
Pre-rewrite (single-best-lag FFT): 961 lines diff.
Post-rewrite (full Falign-port FFT): 74 lines diff. **The structural
rewrite did NOT change BL50 width** — it merely brought our
implementation flow into structural alignment with C.

**Verified equal**:
- `cross_validate_blosum50_n_dis_cell_by_cell` (0 mismatches)
- `cross_validate_blosum50_n_dis_fft_cell_by_cell` (0 mismatches)
- BL50 raw 210-cell lower-triangle table matches `tmpmtx50` exactly
- **NEW** `bl50_alignable_reagion_matches_c` (FFI, 2026-05-08) — calls
  C's `alignableReagion` and Rust's `match_score` + `alignable_segments`
  on the same input. Segments produced are byte-identical, eliminating
  FFT segment scoring as the divergence source.

**Structural rewrite landed 2026-05-07**: `find_fft_anchors`
(`crates/mafft-align/src/fft_align.rs`) now mirrors C's `Falign` flow —
accumulates segments from ALL `NKOUHO=20` candidate lags into a flat
list, sorts independently in seq1 and seq2 dimensions, builds a sparse
`crossscore[rank1+1][rank2+1] = score` matrix with `1e7` corner
sentinels, runs `block_align`. Mirrors `Falign.c:1307-1465` including
the `if (tmpint == 0) break;` early termination on first-empty
candidate.

The rewrite preserves all FFT-NS-2 byte tests (BL62 / BL30 / BL45 /
BL80 / JTT200 still byte-identical) but does NOT change BL50 width
because at step 33 the first candidate (lag=39) returns 0 segments,
the `break` fires, and the function falls back to direct DP — same as
the pre-rewrite single-best-lag approach.

**Pinpointed divergence (2026-05-08 diagnostic dumps)**:

At BL50 step 33 (clusters 13×8):
- Per-step debug output (`MAFFT_DEBUG_STEPS=1`) shows steps 0-32 widths
  and scores match C's reported values; step 33 diverges:
  | Step | Clusters | C width / score | Rust width / score |
  |------|----------|------------------|---------------------|
  | 32   | 13 × 2   | identical        | 508 / 103322.3      |
  | 33   | 13 × 8   | 624 / 47967.6    | 641 / 47954.6       |
  | 34   | 21 × 15  | 721 / 56627.1    | 738 / 56582.8       |
- At step 33, the post-step-32 profiles for clusters [13seq] and [8seq]
  are very gappy (`prof1[0] nongap=0.228`, `prof2[0] nongap=0.220`).
- Top FFT candidate: lag=39, correlation score=29.14. Our
  `shift_and_score` at lag=39: max-window-20-sum=2204.6 <
  threshold=9600. **C's `alignableReagion` must produce the same
  result** since the segment-detection scoring is byte-identical
  (proved by `bl50_alignable_reagion_matches_c`). So C's
  `if (tmpint == 0) break;` fires too — both Rust and C fall back to
  direct full-profile DP.

**The bug is a tie-break divergence in profile DP at step 24.**
Established 2026-05-08 via end-to-end comparison:

The first divergent step is **step 24** (clus1=6, clus2=1, width=375
on both sides). Both Rust and C produce **identical scores
(119988.9)** but **different gap placements**:
- Rust `seq[0]`: `MNGTEG--DNFYVPF----SNK` (F early, 4 gaps after)
- C    `aseq[0]`: `MNGTEG--DNFYVP----FSNK` (4 gaps before F)

This 4-column gap-shift propagates: step 27 widens to 385 in both,
but the score diverges (Rust=111382.1 vs C=111516.6 — first divergent
score). By step 33 the cumulative effect is 17-column width gap.

**FFI tests confirm the DP/segment-detection layers match C
byte-for-byte**:
- `bl50_alignable_reagion_matches_c` — `match_score` +
  `alignable_segments` produce identical segments to C's
  `alignableReagion` (proves FFT scoring is correct).
- `bl50_step33_profile_dp_matches_c` (#[ignore]) — when called
  directly via FFI with the same prof1/prof2 input,
  `profile_align(prof1, prof2, ...)` and `A__align(...)` produce
  byte-identical output (width 638 step 33, width 385 step 26 with
  same scores). Whether sgap1/egap1 are NULL or all-'o', both
  produce the same alignment.

**Conclusion**: The divergence is NOT in score computation, NOT in
FFT, and NOT in the direct DP-with-clean-state. It's a tie-break in
the DP that surfaces only when running the full progressive
flow — likely due to **stale static-buffer state** in C's TLS
buffers carrying over from previous calls. C's static `ogcp1[lgth1]`
or `gapfreq1pt[lgth1]` slots may retain values from the previous
call, since `new_OpeningGapCount` only zeros indices `0..lgth-1`
while `st_OpeningGapCount` explicitly zeros `ogcp[len] = 0`. With BL62
the score landscape is wide enough that this stale data doesn't tip
DP cells; with BL50 (flatter scores) it does.

### Root cause and fix (closed 2026-05-08)

**Bug**: GCC `-O3` (and `clang`/`-O3` on aarch64) fuse `a + b * c` into
a fused multiply-add (FMA) instruction with a single rounding step.
Rust's `+=` followed by `*` produces two separate roundings. For most
of our DP cells the 1-ULP difference is irrelevant, but for matrices
with a flat score landscape (BL50 vs BL62/30/45/80), it surfaces as
tie-break divergences when two equal-score paths exist.

**Trace**: At BL50 step 24 (clus1=6, clus2=1), `match_calc_row`'s
`scarr[3] = sum_j matrix[j][3] * freq[j]` over j=10,12,13 with
matrix entries -226, -226, -303 and freqs 0.1196, 0.3237, 0.3096
produced:
- C (with FMA via `gcc -O3`): -1.93985622366210719747e+02
- Rust (without FMA): -1.93985622366210748169e+02

A 1-ULP difference (~2.84e-13). Through DP accumulation this
propagated to cell (22, 19) where `wm = previousw[j-1]` and
`g_jskip = mi + fgcp2[j-1] * gf1_i` were *equal* in C (both
1.41977729539934080094e+03 — bit-equal) but differed by 4.5e-13 in
Rust. With the strict `>` comparison, C kept the diagonal path
(ijp=0), Rust took the j-skip path (ijp=-5). The traceback then
diverged, placing F at column 28 vs column 32.

**Fix**: `crates/mafft-align/src/profile.rs::profile_align_imp_with_boundary`
now uses `f64::mul_add` for every accumulation involving a `+= a*b`
pattern in `match_calc_row` (substitution score) and the inner DP
loop (gap-frequency-modulated penalties). This produces bit-identical
results to C's gcc-O3-FMA output.

**Verified equal**:
- `bl50_alignable_reagion_matches_c` (FFT scoring, was already passing)
- `bl50_step33_profile_dp_matches_c` (#[ignore], diagnostic; now also
  produces byte-identical output for the §4 reproducer)
- `fftns2_bl50_byte_identical_to_c` (NEW regression test, BL50 width 712)
- All other FFT-NS-2 byte tests (BL30/45/62/80, JTT200) still pass.

### Original symptom (pre-fix, preserved for context)

---

## §5. FFT anchor / segment tie-break on `--tm * (FFT)` — RESOLVED 2026-05-10

**Mode**: `--tm 200` FFT now byte-identical to C. JTT 100 / TM 100 FFT
were already closed by §4's FMA fix.

**Four combined fixes** to the FFT pipeline (`fft_align.rs` + `profile.rs`):

1. **Bit-for-bit FFT C-port** (`crates/mafft-fft/src/fft_c_compat.rs`).
   Replaced `rustfft` with a hand port of MAFFT's Cooley-Tukey radix-2
   `fft()` (`fft.c:7-126`), including the `make_sintbl` /
   `make_bitrev` precomputed tables, butterfly multiply-add order, and
   forward-only `1/n` post-scale. `rustfft`'s scheduling produces
   1-ULP-different correlation values which were enough to flip
   anchor selection on flat-landscape matrices like `--tm 200`.
   Verified bit-identical via `cross_validate_fft.rs`.

2. **`shift_and_score` sign + length convention** (Falign.c:540-548 /
   alignableReagion.c:226-227). At lag `L`, the per-position score
   array is in the LATER-starting sequence's frame, with length
   `MIN(strlen aseq1, strlen aseq2)` after `zurasu2`:
   - L ≥ 0: `scores[i] = match(prof1[i], prof2[i+L])`, `len = min(n, m-L)`
   - L < 0: `scores[i] = match(prof1[i-L], prof2[i])`, `len = min(n+L, m)`
   Pre-fix Rust used `j = i - lag` with `len = max(n, m)` and
   zero-padding, which gave the OPPOSITE physical pairing — invisible
   at lag 0 but wrong everywhere else.

3. **Anchor coordinate mapping** in `find_fft_anchors`. With the new
   sign convention, segment center `c` (in the shifted frame) maps to
   absolute coords `(c1, c2)` where `c2 - c1 = lag`:
   - L ≥ 0: `c1 = center`, `c2 = center + lag`
   - L < 0: `c1 = center - lag`, `c2 = center`

4. **`align_with_anchors` as segment-boundary DP** (Falign.c:684-705).
   Pre-fix Rust treated anchors as forced match cells, emitting an
   `AlignOp::Match` at each anchor and DP'ing only the gaps between.
   C runs `count - 1` independent profile-DPs over segments
   `[cut1[i], cut1[i+1]) × [cut2[i], cut2[i+1])` with cuts bracketing
   the anchors — anchor positions are the START of the next segment,
   the DP decides match/gap freely. Now mirrored.

5. **2-channel polarity+volume FFT for protein**. Rust was always
   building 20-channel indicator vectors (one per residue). C's
   `Falign.c:342-348` takes the `seq_vec_2` path for protein
   (`fftscore && scoremtx != -1`), using just 2 channels filled with
   per-residue polarity and volume. Added
   `profile_to_property_channels` + a `property_channels` field on
   `FftAlignParams`; engine populates it from `scoring.polarity` /
   `scoring.volume` for non-nucleotide modes. This was the load-bearing
   fix for TM 200 — the 20-channel indicator correlation has different
   peak structure than C's 2-channel polarity+volume on flat matrices.

**Regression guards**: existing `cross_validate_fft.rs` (FFT FFI),
`cross_validate_bl50_fft.rs` (segment scoring + step-24 DP), plus the
TM 200 FFT byte-identity row in the parity matrix above.

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
