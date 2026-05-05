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
| G-INS-1 (`--globalpair --maxiterate 0`)    | 746     | 747        | 974        | divergent — see §1 |
| E-INS-1 (`--genafpair --maxiterate 0`)     | 729     | 719        | 948        | divergent — see §2 |
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

## §1. G-INS-1 — divergent without refinement

**Mode**: `--globalpair --maxiterate 0` (single progressive pass over
pairwise global-align-derived tree).

**Observed**: Rust width 747 vs C 746, diff 974 lines (massive content
divergence despite near-equal width). The first sequence's gap pattern
diverges at column 1 — this is not a refinement issue; it's a
progressive-pass issue.

**Suspected root**: our `PairAligner::Global` path through
`build_homology_table` (engine.rs:232) uses `global_align`. The
opt/importance computation may not match the L-INS-i analog after the
2026-05-05 fixes (full-precision distance, `iscore / sumoverlap` opt
scale, `recompute_importance` symmetrization). Verify each step by
reusing the FFI cell-by-cell test pattern from
`crates/mafft-core/tests/cross_validate_constrained_align.rs`.

**Concrete next task**:
1. Add an FFI test mirroring `constrained_align_real_2seq_matches_c` but
   driving `G__align11` instead of `L__align11` for the pairwise build.
2. Compare our `region.opt`, `region.importance`, and final impmtx
   to C's by running `mafft --globalpair --maxiterate 0 --debug` on a
   3-seq subset and extracting `hat3` + dumping intermediate state via
   the same `IMP_DUMP` instrumentation pattern used for L-INS-i.
3. The fix is likely small (importance/opt scaling, or a missing
   `outgap` setting on the global path). Effort: 0.5-1 day.

---

## §2. E-INS-1 — uses local pairwise instead of generalized-affine

**Mode**: `--genafpair --maxiterate 0` (E-INS-i without refinement).

**Observed**: Rust width 719 (10 narrower than C's 729), diff 948 lines.

**Root cause**: `engine.rs:229-230` selects `PairAligner::Local` for both
`LInsi` and `EInsi`. C's E-INS-i uses **generalized-affine** pairwise
(`genL__align11` / `pairlocalalign -E` / `pairlocalalign -K`), which has
distinct gap-extension and skip costs that produce a wider local window
than plain Smith-Waterman.

**Concrete next task**:
1. Add `PairAligner::GeneralizedAffine` variant in
   `crates/mafft-align/src/constraints.rs`.
2. Implement `genaff_local_align` (port of C's `genL__align11` from
   `genalign11.c`) and route it through `build_homology_table`.
3. Wire `EInsi { .. }` to `PairAligner::GeneralizedAffine` in
   `engine.rs::align`.
4. Mirror the same opt/importance pipeline (the 2026-05-05 fix should
   apply unchanged: opt = `iscore / sumoverlap`, full-precision
   distance, calcimportance_half symmetrization).

Effort: 1-2 days (the DP itself is ~150 lines; the integration mirrors
L-INS-i exactly).

---

## §3. INS-i modes with iterative refinement (`--localpair`,
   `--globalpair`, `--genafpair` defaults run 1000 iterations)

**Status as of 2026-05-05**:

| Mode    | C width | Rust width | Notes |
|---------|---------|------------|-------|
| L-INS-i | 719     | 725        | Was caused by `--maxiterate 0` parsing bug AND distance rounding; both fixed for `--maxiterate 0` (now byte-exact). With refinement enabled, refinement loop drifts. |
| G-INS-i | 746     | 714        | §1's progressive pass also diverges; refinement on top compounds. |
| E-INS-i | 729     | 725        | §2's pairwise mismatch + refinement drift. |

**What was fixed during 2026-05-05 session** (all in
`crates/mafft-core/src/engine.rs` and one in `crates/mafft-bin/src/main.rs`):

1. **`opt` scaling bug** (closed ~10 columns of L-INS-i gap):
   `crates/mafft-align/src/constraints.rs` was computing
   `opt = isumscore * 5.8 / (600 * sumoverlap)` (the value C's
   `putlocalhom2` writes to hat3), but C's `tbfast.c:2202` rescales by
   `* 600 / 5.8` immediately before calcimportance_half so the value
   actually consumed downstream is `isumscore / sumoverlap`. Our
   `recompute_importance` was producing impmtx values 100× smaller
   than C's. Now stores the post-rescale value directly.
2. **Distance rounding on the L-INS-i path** (closed ~6 columns):
   `engine.rs` was rounding distances to 3 decimals (mimicking the
   hat2 file round-trip via `%#6.3f`). That's correct for FFT-NS-i +
   dndpre but wrong for L-INS-i: `tbfast` with `callpairlocalalign=1`
   computes pairwise alignments in-memory and feeds full-precision
   `iscore[]` straight to `fixed_musclesupg_*`. Removed the rounding
   for the constraint-aware path.
3. **`--maxiterate 0` was treated as "use mode default 1000"**:
   `args.maxiterate` was a `usize` defaulting to 0, and
   `iters_for(1000)` returned 1000 whenever `args.maxiterate <= 0`.
   Changed to `Option<usize>` so explicit `Some(0)` disables refinement.
   This was the entire reason `--maxiterate 0` previously appeared to
   diverge by 6 columns — refinement was secretly running.

**Remaining gap (refinement loop)**: full default `--localpair` (1000
refinement iters) still produces width 725 vs C's 719. Without
refinement (`--maxiterate 0`) we're byte-exact, so the divergence is
purely in `iterative_refine` (`crates/mafft-core/src/refinement.rs`)
when constraints are present.

**Concrete next task** (effort: 1-2 days):
1. Capture C's per-iteration accept trajectory via
   `mafft --localpair --maxiterate 100 --debug` and extract the
   per-step DP score / width before/after each branch refinement
   from the trace files left in `--debug`'s tmpdir.
2. Compare against our `RUST_MAFFT_TRACE=1` output for the same input.
   The first iteration where Rust accepts a different rearrangement
   than C is the divergence point.
3. Likely candidates (in order of likelihood):
   - Constraint-aware refinement DP doesn't add the impmtx bonus per
     cell the same way the constrained progressive merge does. Compare
     `iterative_refine`'s DP call pattern to `progressive_align_with_
     constraints::merge_step_cached`.
   - Per-iteration tree weight recomputation: C may rebuild branch
     weights between iterations using the current alignment, we may
     not.
   - Random-tie-break order: C iterates branches in a specific order
     across iterations (alternating directions); verify ours matches.

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

1. **§3 (INS-i refinement loop)** — closing this completes the original
   user goal of "L/G/E-INS-i with refinement byte-exact". Most of the
   pieces are in place after the 2026-05-05 fixes; the work is now
   isolated to `iterative_refine` with constraints. Highest value per
   hour.
2. **§1 (G-INS-1)** — small contained bug in the pairwise/progressive
   pass. The L-INS-i fix from this session should largely apply; verify
   and extend.
3. **§2 (E-INS-1)** — adds a new pairwise variant. Independent of §1
   and §3 but exercises the same constraint pipeline.
4. **§4 (BL50 FFT)** — structural FFT change; once this lands, §5
   likely closes for free.
5. **§7 (--add)** — unblocks §8 (perf) and is a real feature gap.
6. **§6 (parttree)** — algorithmic port, substantial work.
7. **§9 (RNA / allowshift)** — as needed.
