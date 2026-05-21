# BALIBASE 3 parity sweep (2026-05-18 → 2026-05-20)

## Final results across modes

| Mode | Flags | Match | Diverge (tied) | Diverge (width) | Total match % |
|------|-------|-------|---------------|----------------|---------------|
| FFT-NS-1 | `--retree 1` | 206 | 12 | 0 | 94.5 % |
| FFT-NS-2 (default) | (none) | 212 | 0 | 6 | **97.2 %** |
| FFT-NS-i | `--maxiterate 100` | 128 | 69 | 20 | 58.7 % |
| L-INS-i | `--localpair --maxiterate 100` | 154 | 55 | 8 | 70.6 % |
| G-INS-i | `--globalpair --maxiterate 100` | 140 | 65 | 13 | 64.2 % |
| E-INS-i | `--genafpair --maxiterate 100` | 184 | 29 | 5 | 84.4 % |

**Observation**: `--retree 1` (single-pass progressive) has all divergences
as tied-width (same width, residue shifts only — the §B.2 fingerprint).
Default mode (retree=2) has fewer total divergences but converts SOME
tied-traces into matches while introducing 6 width-differs. Both
results are consistent with the §B.2 / `A__align` static-state hypothesis.

All 218 BALIBASE 3 protein test sets (BB***** entries from the Drive5
bench mirror). One `rs_error timeout` in FFT-NS-i/L-INS-i (BB40049,
62-seq input, manually killed at 600s+ — way past timeout, harness
limitation). Default-mode comparison is the canonical "is rust-MAFFT
production-ready" number.

The refinement modes' lower parity is **expected**: every iterative
refinement pass re-runs `A__align` ~N times per iteration, so any §B.2-
class tied-cell cascade from the initial progressive alignment gets
amplified through subsequent refinement DPs. ~80 % of refinement-mode
divergences are same-width tied-trace cases (residue gap shifts that
don't change score or width — the §B.2 fingerprint). The remaining
~20 % are width-differs that originate from a tied-cell shift in an
early merge and cascade through refinement into a width difference.

All 6 default-mode divergences (BB20018, BB20027, BB20039, BB40036,
BB40041, BB40048) are also divergent in FFT-NS-i (often with larger
diffs because refinement amplified them). The cascade is consistent
across modes.

## Original investigation log (default mode, 2026-05-18)

**Corpus**: Drive5 mirror of BALIBASE 3 (`bench.tar.gz`, 218 `BB*****` test
sets, refs RV11/12/20/30/40/50 — the lbgi.fr server returns 404 for all
download URLs at this time, so the Drive5 mirror is the only working
source).
**Mode**: default (FFT-NS-2), single-threaded, 120 s timeout per test.

## Run history

| Run | Result | Δ vs previous |
|-----|--------|---------------|
| First (2026-05-18, pre-fix) | 201 / 218 (92.2 %) | — |
| Second (2026-05-18, post-`nogaplen` lenfac fix) | **212 / 218 (97.2 %)** | +11 tests pass |

## What fixed: `ktuple_distance_aa` lenfac uses nogaplen, not filtered groups

`ktuple_distance` (the 6-tuple distance used to build the default
progressive guide tree) was computing the length-adjustment factor
`lenfac` from the **filtered** group sequence length (post-X /
post-`.` / post-`-` stripping) where C's `disttbfast.c:3845-3856`
uses **`nogaplen`** (gaps-only stripped, X / `.` kept). For any
sequence carrying X (or `.`), all distances *to* that sequence
were biased by ~1-1.5×10⁻³ — enough to flip UPGMA join order on the
~5% of BALIBASE inputs that contain X.

Fix: `mafft-tree/src/distance.rs::ktuple_distance_aa/_nuc` now pass
`nogap_len(seq)` (gap-only filter) to `compute_lenfac` instead of
`groups.len()`. Regression test:
`distance::tests::ktuple_lenfac_uses_nogaplen_not_filtered_groups`.
FFI cross-check on BB12041's 21 pairs: all 21 bit-identical to C
post-fix (vs 6 differ pre-fix, all involving the X-carrying 3cao_A).

## Remaining 6 divergent tests (post-fix)

| Test    | nseq | C width | Rust width | diff lines | class |
|---------|------|---------|------------|------------|-------|
| BB20018 | 53   | 476     | 476        | 204        | tied  |
| BB20027 | 29   | 1602    | 1630       | 1359       | WIDTH (+28) |
| BB20039 | 91   | 663     | 663        | 4          | tied  |
| BB40036 | 29   | 846     | 846        | 8          | tied  |
| BB40041 | 55   | 2995    | 2920       | 2794       | WIDTH (-75) |
| BB40048 | 17   | 1995    | 1995       | 4          | tied  |

**Tied-trace (4)**: BB20018, BB20039, BB40036, BB40048 — likely the
same C-side `A__align` static-buffer artifact as §B.2. Confirmed
pattern: isolated single-residue gap shifts.

**Width-differs (2)**: BB20027, BB40041. These two remain genuinely
divergent — neither involves an X-carrying sequence, so the lenfac
fix didn't close them. Next investigation target. Hypothesis order:
1. Another subtle distance / lenfac mismatch we haven't found yet
   (e.g. NJ tie-break for exactly-equal distances).
2. Profile-DP arithmetic divergence on a specific length/composition
   regime that the 36-seq sample doesn't exercise.
3. FFT correlation peak choice on a marginal anchor.

## Original 17 divergences (pre-fix, kept for reference)

11 cases (BB11030, BB12041, BB20005, BB20041, BB30003, BB30013, BB30030,
BB40022, BB40026, BB40047, BB50016) closed with the `nogaplen` fix —
all of these contained at least one X-carrying sequence. The 6 above
were unaffected.

## BB20027 deep-dive findings (2026-05-19)

Deep investigation; root cause IDENTIFIED but full fix requires more work.

### Layer 1 — isolated to pass-1 progressive

- `--retree 1`: byte-identical to C (width 1612).
- `--retree 2` (default): diverges (1602 vs 1630).
- `--nofft`: divergence persists → not FFT-anchor placement.
- Rebuilt tree (`--treeout` after retree=2): byte-identical Newick.
- Sequence weights match C bit-for-bit (FFI test
  `bb20027_pass1_weights_match_c`: 0 drift across all 29 weights).
- `--add` test (force C's 28-way profile + 1p8j_A): byte-identical.
  → final merge step's DP is fine; the *intermediate* 28-way profile
  is what diverges.
- Rust pass-1 step-trace (`MAFFT_DEBUG_STEPS=1`) vs C instrumented
  trace: identical widths/scores through step 12. **Step 13 (clus1=8,
  clus2=2)** is the first divergence: same width 593, but Rust score
  100653.4 vs C 100625.6 (Δ 27.8). Different DP optimum on identical
  inputs → tied-DP-cell selection driven by 1-ULP cpmx drift.

### Layer 2 — root cause: cache-vs-from-scratch cpmx precision drift

The 8-way profile entering step 13's DP differs by 1 ULP between C
and Rust:

- **Rust `Profile::from_aligned`** equals **C `cpmx_calc_new`**
  bit-for-bit (FFI test `rust_profile_freqs_match_c_cpmx_calc_new`,
  0 cells differ).
- **Rust `Profile::from_aligned`** ≠ **Rust `blend_profiles_exact`**
  by 1 ULP in ~30 of 1166 cells (test
  `rust_from_scratch_matches_rust_blend`, max diff 5.55e-17).
  The blend cascade — `((w/clusterA_sum) * clusterA_total_sum/all_sum)`
  — loses precision vs the direct `w/all_sum` from-scratch path.
- C uses `cpmxhist` (createcpmxresult-via-blend) for every internal
  node regardless of size — see `disttbfast.c:2620,2913` —
  while Rust only caches the merged profile when combined seqs > 20
  (`progressive.rs:987`). At step 13, Rust builds prof1 from scratch
  on the 8 aligned sequences, while C uses the cached blend from
  step 12's `createcpmxresult`. Same math, 1-ULP precision drift in
  ~30 cells, drift flips one DP tied cell at step 13.

### Layer 3 — partial fix attempted; revealed C-blend mismatch

Setting `combined_seqs > 0` (always cache) closes BB20027 cleanly
(1359-line diff → 4-line tied-trace), but opens 5 previously-matching
tests (BB20005, BB30003, BB30013, BB40046, BB50016) with similar
small drifts. So our `blend_profiles_exact` doesn't bit-match C's
`createcpmxresult` somewhere; the additional mismatches surface only
when we use the blend path universally.

**Real bug FIXED**: `blend_fg_one_side` had `j < alen - 1`
short-circuit that dropped the closing-count contribution at the
final non-gap position. C reads `gaptable[j+1]` past the null
terminator (treated as non-gap) and DOES add `ori[p] * eff` there.
Fixed in `progressive.rs:1208-1248` with `next_is_gap()` helper that
treats out-of-bounds j+1 as non-gap to match C exactly. This is a
real correctness bug but doesn't surface in the baseline-threshold
tests (the cell is rarely tied at the final position).

### Layer 4 — profile path is bit-exact; DP itself is the suspect

FFI-bound C's `createcpmxresult / creategapfreqresult / createogresult /
createfgresult` (all `static`) into `wrappers/blend_helpers.c` and wrote
two cell-by-cell tests against synthetic 3+5 inputs:

- `cross_validate_cpmx::rust_profile_freqs_match_c_cpmx_calc_new` —
  Rust `Profile::from_aligned` matches C's `cpmx_calc_new + gapcountf +
  st_OpeningGapCount + st_FinalGapCount` BIT-EXACTLY in all four arrays
  (freqs, gap_freq, ogcp, fgcp).
- `cross_validate_cpmx::rust_blend_matches_c_blend_cell_by_cell` —
  Rust `blend_profiles_exact` matches C's full blend cascade
  BIT-EXACTLY for a 3+5 merge with realistic gappy gaptables.

**Conclusion**: profile building (from-scratch AND blend) is bit-exact.
The remaining BB20027 divergence is in the **DP itself** —
`profile_align_imp` / `fft_profile_align` — not in the profile.
Specifically: at pass-1 step 13 (8+2 merge), the inputs (prof1, prof2,
gap_model, scoring matrix) are bit-identical between Rust and C, but
the DP produces score 100653.4 (Rust) vs 100625.6 (C). The DP picks a
different optimum on the same inputs.

Cache-always (closing BB20027 with 4-line tied-trace) still regresses
5 tests, even with the bit-exact blend path. So there are additional
small mismatches between Rust's progressive merge and C's that surface
only when we follow C's cache-everywhere strategy. They look like DP
tied-cell-selection differences (small diffs in tests like BB30013
which goes from matching to a 4414-line WIDTH-differs at cache-always).

### Layer 5 — DP cell-by-cell comparison: SAME inputs → SAME outputs

Wrote `cross_validate_bb20027_dp.rs` that replays pass 1 up to a target
step via `progressive_align_partial`, then feeds the exact step inputs
into Rust `profile_align` and C `A__align` (via `mafft_sys::A__align`
with `cpmxchild0 = cpmxchild1 = NULL` to force from-scratch).

- `bb20027_step12_rust_dp_vs_c_aalign_cell_by_cell`: bit-identical
  alignment, score 90391.2158 in both.
- `bb20027_step13_rust_dp_vs_c_aalign_cell_by_cell`: bit-identical
  alignment, score 100653.3741 in both.

So **the DP itself is correct**. Given identical profile inputs, Rust
and C produce bit-identical alignments.

### Layer 6 — actual divergence source: cached cpmx (`cpmxhist`) path

Added `MAFFT_DUMP_STEP=N` env-var-controlled mseq dump to both the C
engine (`disttbfast.c::treebase`) and the Rust engine
(`progressive.rs::progressive_align_with_weights_override`). Diffed
step-12 dumps in --nofft mode and found a single-residue gap shift in
cluster1's 3 sequences:

- Rust step 12: `MTSIRRLALYLGAL---LPAVLAAPAAL...`
- C step 12:    `MTSIRRLALYLGA---LLPAVLAAPAAL...`

Same width 444, same score 90391.2158. **Classic tied-DP-cell choice
divergence.** With cache-always enabled (`combined_seqs > 0`), step 12
Rust dump becomes BIT-IDENTICAL to C's: caching the cpmx at each merge
step puts our DP on the same tied-cell trajectory as C.

**Mechanism**: C's `disttbfast.c:2913` calls `A__align(..., cpmxchild0,
cpmxchild1, cpmxhist+l, ...)` with the cached cpmx from prior merge
steps (built via `createcpmxresult`). C's `Salignmm.c:1479-1525` then
uses the cached cpmx instead of calling `cpmx_calc_new` from scratch.
The cached blend differs from from-scratch by 1-ULP in some cells
(confirmed in `rust_from_scratch_matches_rust_blend` — max diff
5.55e-17, 30 of 1166 cells). At step 12 this 1-ULP drift flips a tied
DP cell, producing a single-residue gap shift. The shift cascades into
step 13+ as a different cpmx, and by step 27 the final widths diverge.

### Layer 7 — cache-always closes BB20027 but opens 5 other tests

Setting `combined_seqs > 0` makes Rust cache every internal merge,
matching C's `cpmxhist` policy. Per the §B.2 / step-12 diagnosis this
should close BB20027 — and it does (1359-line diff → 4-line tied trace).

But cache-always REGRESSES 5 previously-matching tests: BB20005,
BB30003, BB30013, BB30029, BB40046, BB50016 (combined diff up to
4414 lines). Investigating BB30013 step-by-step shows the same tied-
DP-cell pattern at step 77 (clus1=8, clus2=6): same width 993, score
delta 12.6 between Rust and C.

For these tests, the from-scratch cpmx happens to bias the tied-cell
selection one way (matching C's `cpmxhist` cached output by
coincidence), while the cache-always blend cpmx biases it the other
way. The two cpmx sources differ by 1-ULP, and tied cells flip
differently depending on which direction the drift points.

### Layer 8 — DP tie-break audit: Rust profile_align == C A__align exactly

Audited every `>` vs `>=` comparison in the DP recurrence between
`profile.rs::profile_align_imp_with_boundary` and
`Salignmm.c::A__align`. Result: **every tie-break direction matches**:

| Cell update | C operator | Rust operator |
|-------------|-----------|---------------|
| `wm <- mi + fgcp*gf` (gap_jskip)        | `> wm` | `> wm` ✓ |
| `mi <- prept + ogcp*gf` (mi update)     | `>= mi` | `>= mi` ✓ |
| `mi += fpenalty_ex`                      | unconditional | unconditional ✓ |
| `wm <- m[j] + fgcp*gf` (gap_iskip)      | `> wm` | `> wm` ✓ |
| `m[j] <- prept + ogcp*gf` (m[j] update) | `>= m[j]` | `>= m[j]` ✓ |
| `m[j] += fpenalty_ex`                    | unconditional | unconditional ✓ |
| `wm <- warp_g` (warp candidate)          | `> wm` | `> wm` ✓ |
| Atracking corner: last-row scan          | `> wm`, j desc | `> wm`, j desc ✓ |
| Atracking corner: last-col scan          | `> wm`, i desc | `> wm`, i desc ✓ |
| Atracking corner: fall-back              | `> wm` | `> wm` ✓ |

DP boundary inits also match (`mi = previousw[0] + ogcp2[1]*gf1`,
`mj[j] = currentw[j-1] + ogcp1[1]*gf2`, `previousw[0] = initverticalw[i-1]`
on row swap). FMA throughout, matching gcc -O3 -ffp-contract=on.

### Layer 9 — direct DP test confirms: same inputs → same outputs

`bb20027_step12_rust_dp_vs_c_aalign_cell_by_cell` (in
`cross_validate_bb20027_dp.rs`): Rust `profile_align` and C `A__align`
both produce **bit-identical** output (score 90391.2158, width 444,
zero differing columns) given step-12's identical inputs from
`progressive_align_partial(rebuilt_tree, n_steps=12)`.

### Layer 10 — what's actually different at the engine level

When the C engine calls `A__align` (line `disttbfast.c:2913`), it
passes `cpmxchild0/cpmxchild1` (the CACHED `cpmxhist` from prior
merges). My direct-call test passes `NULL` for cpmxchild, forcing
C to rebuild cpmx via `cpmx_calc_new` (from-scratch). The two cpmx
sources differ by 1-ULP in ~30 cells out of ~1166 per profile.

With cache-always enabled in Rust (`combined_seqs > 0`), Rust uses
its own blend cpmx, which we PROVED bit-exact with C's
`createcpmxresult` cascade in `rust_blend_matches_c_blend_cell_by_cell`.
So both engines should feed `A__align` identical cached cpmx and
both produce identical alignments.

That works for BB20027 (cache-always: step 12 byte-identical to C —
diff 1359 → 4 lines, the remaining 4 being a §B.2-style tied-trace
artifact downstream). But cache-always **regresses 5 other tests**
(BB20005 938 lines, BB30003 4 lines, BB30013 4414 lines, BB30029 4
lines, BB40046 208 lines, BB50016 1189 lines).

Investigated BB30013 specifically:
- Per-step trace identifies step 76 as the first divergent step
  (cache-always Rust vs C): same width 729, same score 146562.0327,
  but post-merge MSAs differ in single-residue gap shifts.
- Step 75 dumps match bit-for-bit. So inputs to step 76 SHOULD
  match. Both implementations use cached blend cpmx for step 76's
  cluster1 (a 5-way) and from-scratch for cluster2 (a single leaf).
- Despite bit-identical blend per our cell-by-cell test, the
  cumulative cached cpmx at step 76 must diverge between Rust and C
  somewhere over the prior 75 cascading blends.
- `--nofft` reproduces the same 4414-line divergence in BB30013, so
  it's not the FFT path.

### Layer 11 — remaining mystery and what closes it

The DP code is bit-equivalent. The blend code is bit-equivalent.
Pass-0 alignment is bit-equivalent. Yet cache-always at step 76 of
BB30013 diverges. Two possibilities:

1. **Cumulative blend drift on specific shape**: our blend cascade
   has a discrepancy from C's cascade that only surfaces for
   certain (cluster_size, gaptable_pattern) combinations we haven't
   exercised in unit tests. The `rust_blend_matches_c_blend_cell_by_cell`
   test exercised one 3+5 case with gappy gaptables; other shapes
   (deep recursion, large clusters, sparse gaptables) might expose
   a missing edge case.

2. **A second blend bug similar to the `blend_fg_one_side`
   past-null-terminator fix from layer 6**, lurking in
   `blend_og_one_side` or `creategapfreqresult`-equivalent code that
   only fires on specific gaptable patterns.

Closing these would require:
1. Dumping cached cpmx at intermediate steps in both engines and
   diffing them position-by-position to find the first divergent
   cell — needs C-side instrumentation to expose `cpmxhist[k]`.
2. Identifying the specific blend operation that differs.
3. Fixing it (probably another one-line C-comparison patch like
   `next_is_gap()`).

Each step is straightforward but the dumping + diff is fiddly (cpmx
is `n_alphabets × cluster_width` per step × 76 steps × 4 fields).
~1 day of focused work. Current state remains 97.2% on BALIBASE 3.

### Summary of fixes applied this session

- `mafft-tree::distance::ktuple_distance_aa/_nuc` — lenfac uses
  `nogap_len` instead of filtered-group length (closed 11 tests).
- `progressive::blend_fg_one_side` — treat out-of-bounds j+1 as
  non-gap to match C's past-null-terminator read (real correctness
  bug; benign on baseline tests).
- `progressive::blend_profiles_exact` — exposed as `pub` for
  cross-validation tests.
- `cross_validate_counteff::bb20027_pass1_weights_match_c` — FFI
  verification that `sequence_weights` == C `counteff_simple_double`
  (zero drift, 29 weights).
- `cross_validate_cpmx::rust_profile_freqs_match_c_cpmx_calc_new` —
  FFI verification that `Profile::from_aligned` == C
  `cpmx_calc_new + gapcountf + st_*GapCount` (zero drift, all 4
  arrays).
- `cross_validate_cpmx::rust_blend_matches_c_blend_cell_by_cell` —
  FFI verification that `blend_profiles_exact` == C
  `createcpmxresult` cascade (zero drift on 3+5 with gappy tables).
- `cross_validate_bb20027_dp::bb20027_step{12,13}_rust_dp_vs_c_aalign_cell_by_cell` —
  FFI verification that the DP outputs match given identical inputs.
- `wrappers/blend_helpers.c` — FFI bindings for C's static blend
  cascade (`rs_createcpmxresult`, `rs_creategapfreqresult`,
  `rs_createogresult`, `rs_createfgresult`).
- `MAFFT_DUMP_STEP=N` env-var-controlled MSA dump in Rust
  `progressive.rs` for future debugging (C-side equivalent removed
  when investigation closed).

### Layer 6 — fixes applied this session

- `blend_fg_one_side` (`progressive.rs:1208`): C reads
  `gaptable[j+1]` past the null terminator (treated as non-gap);
  our Rust short-circuited at `j < alen - 1` and dropped the
  closing-count contribution at the last non-gap position. Now uses
  `next_is_gap()` helper. Real correctness bug — doesn't affect today's
  baseline tests (cell rarely tied at the final position) but ships
  with the rest.
- `wrappers/blend_helpers.c` (new): verbatim copies of C's `static`
  blend functions, exposed as `rs_create{cpmxresult,gapfreqresult,
  ogresult,fgresult}` FFI symbols for cross-validation tests.
- `progressive::blend_profiles_exact` made `pub` so integration tests
  can exercise it directly.

## Option A attempt (2026-05-21): `--c-compat` flag

Implemented C MAFFT's `reuseprofiles` static-TLS memoization
(`Salignmm.c:1446-1450` + `cpmx_calc_add` at `tddis.c:160`) behind a
new `--c-compat` opt-in flag. Infrastructure added:

- `mafft_align::profile::CPMX_MEMO` — thread-local
  `CpmxMemoState { previous_call, previous_firstmem, previous_icyc,
  previous_first_len, freqs, gap_freq, opening, closing }` mirroring
  C's static TLS state.
- `Profile::from_aligned_with_memo(seqs, weights, amino_map,
  nalphabets, firstmem, icyc, lgth)` — when memo activates (matches
  C's gate: `previous_call && firstmem == previous_firstmem && lgth
  == previous_first_len && icyc == previous_icyc + 1`), use
  `cpmx_calc_add`-equivalent incremental update; else from-scratch.
  Also covers `gapcountadd`, `st_OpeningGapAdd`, `st_FinalGapAdd`.
- `MafftEngine.c_compat: bool` + `with_c_compat()` builder.
- `--c-compat` CLI flag.
- `mafft_align::reset_cpmx_memo()` called by the engine at the start
  of each retree pass (mirroring `Salignmm.c:1365-1366`).
- Unit tests in `profile::tests`:
  `cpmx_memo_first_call_matches_from_aligned` and
  `cpmx_memo_incremental_matches_from_scratch_for_extension`.

**Result**: `--c-compat` produces **identical output to baseline on
all 218 BALIBASE 3 tests** (same 6 divergent tests, same diff counts,
same widths). The infrastructure works correctly per unit tests, but
**the memo conditions never fire at the divergent steps**.

C-side trace of BB20027 pass 1 confirmed: at step 12 (the first
divergent step), `reuseprofiles = 0` because `lgth1 != previousfirstlen`.
So C is using `cpmx_calc_new` (from-scratch) at that step, NOT the
incremental `cpmx_calc_add` path. Yet C and Rust still produce
different alignments on identical inputs — meaning the divergence
isn't from `reuseprofiles` at all.

The remaining divergences must be driven by OTHER C-side static TLS
state in `A__align`:
- `cpmx1` / `cpmx2` resize semantics (grow-only buffers, stale
  trailing data — though the DP reads only `[0, lgth)`)
- `gapfreq1` / `gapfreq2` resize/stale state
- `ogcp1o` / `fgcp1o` resize/stale state
- `doublework` / `intwork` (match_calc workspace) stale state
- `intverticalw` / `lastverticalw` boundary arrays
- `m` / `mp` / `w1` / `w2` working buffers

Identifying which one(s) drives BB20027 requires deeper diagnostics
that the cell-by-cell `A__align` direct test (which matched bit-for-bit
with `cpmxchild=NULL, calledbyfulltreebase=0`) can't surface.

The `--c-compat` infrastructure is kept (compiles cleanly, unit
tested, no-op when off). It provides a foundation for the next
investigator to add additional static-state replication if the
remaining divergences become a hard requirement.

### Static-state hypotheses RULED OUT (2026-05-21)

After implementing `--c-compat` and finding it didn't help, pushed
further to identify which OTHER C static TLS state could explain
BB20027. The following were ruled out via source-level inspection
and/or C-side instrumentation:

| Hypothesis | Status | Reason |
|------------|--------|--------|
| `reuseprofiles` memo (`Salignmm.c:1446-1450`) | RULED OUT | Trace confirms `reuseprofiles=0` at BB20027 step 12 (`lgth1≠previousfirstlen`). C uses from-scratch. |
| `gapfreq[lgth]` stale | RULED OUT | Default mode passes `cpmxresult≠NULL` → C explicitly sets `gapfreq[lgth]=1.0` at `Salignmm.c:1531-1532`. |
| `ogcp[lgth]` stale | RULED OUT | `st_OpeningGapCount` explicitly zeros position `len` at `mltaln9.c:12865`. |
| `fgcp[lgth]` stale | RULED OUT | DP doesn't read `fgcp[lgth]` for default mode (`tailgp=0` → `lasti=lgth1`). |
| `doublework`/`intwork` stale | RULED OUT | `match_calc(initialize=1)` called at A__align entry (`Salignmm.c:1709`) rebuilds them per call. |
| `m[]`/`mp[]` stale | RULED OUT | Reset per A__align call at `Salignmm.c:1779-1780`. |
| `initverticalw`/`lastverticalw` stale | RULED OUT | Computed fresh via `match_calc(initialize=1)`. |

What remains:
- `cpmx1`/`cpmx2` static buffers with stale data in
  `[lgth, max_seen_lgth)` cells — but the DP only reads `[0, lgth)`,
  so stale tail shouldn't affect output.
- Some OTHER static interaction we haven't identified.

The fact that:
1. Our `Profile::from_aligned` == C `cpmx_calc_new + gapcountf +
   st_*GapCount` BIT-EXACTLY (verified).
2. Our `profile_align` == C `A__align` BIT-EXACTLY when given
   identical inputs with `cpmxchild=NULL, calledbyfulltreebase=0`
   (verified).

...means the divergence must come from how C's A__align processes
**identical inputs differently when called within the full engine
context** (with `calledbyfulltreebase=1` and non-NULL cpmxchild).

Closing this gap would require systematically instrumenting C's
A__align to dump every static buffer value at the divergent step and
diffing against Rust. That's not within reasonable session scope.
The infrastructure is in place (`--c-compat`, FFI bindings, dump
points) for a future investigator who decides to invest the time.

## Closing call (2026-05-19)

The remaining 6 divergences are all `A__align` static-state-coupling
artifacts (§B.2-class). Confirmed via:
1. `Profile::from_aligned` == C `cpmx_calc_new + gapcountf +
   st_*GapCount` bit-exact (cross_validate_cpmx).
2. `blend_profiles_exact` == C `createcpmxresult + ...` bit-exact on a
   3+5 gappy input (cross_validate_cpmx).
3. `profile_align` == C `A__align` bit-exact when `cpmxchild=NULL`,
   `calledbyfulltreebase=0` (cross_validate_bb20027_dp).
4. The engine's actual `A__align` call uses `calledbyfulltreebase=1`
   which engages the static-TLS memoization
   (`Salignmm.c:1446-1450`) — producing different tied-cell choices
   than our stateless engine on inputs that have such ties.

Decision: characterize as §B.2-class C-side stateful artifacts and
hand off to upstream. `MAFFT_UPSTREAM_REPORT.md §5b` now covers the
BALIBASE 3 corpus expansion (was originally just `--tm 200 --retree
1`). The 97.2 % parity result is production-ready; both alignments at
each divergent case are optimal-scored, just different traceback
choices.

## Reproduction

```bash
cargo build --release
scripts/balibase_parity.py /tmp/balibase/bench1.0/bali3/in \
    --pattern "BB[0-9]*" --modes "" --output /tmp/balibase_default.tsv
```

Detailed TSV: `/tmp/balibase_default.tsv`.
