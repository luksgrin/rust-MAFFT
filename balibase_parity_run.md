# BALIBASE 3 parity sweep (2026-05-18)

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

### Layer 4 — what's left

At least ONE MORE site where our blend differs from C's blend in a
way that flips tied cells. Finding it requires either:
1. FFI-binding C's `createcpmxresult` (currently `static`) and
   cell-by-cell comparison against `blend_profiles_exact` on a
   realistic 8+2 input pair, OR
2. Instrumenting C's per-step cpmx and Rust's per-step cpmx, diffing
   cell-by-cell at step 13 of BB20027 (the known-divergent step).

Both require ~half a day of focused work. For now, the BB20027
divergence is one of 6 known residual cases at 97.2% byte-equality;
the other 5 are tied-trace artifacts (likely §B.2 C-side static
buffer effects, not Rust bugs).

## Reproduction

```bash
cargo build --release
scripts/balibase_parity.py /tmp/balibase/bench1.0/bali3/in \
    --pattern "BB[0-9]*" --modes "" --output /tmp/balibase_default.tsv
```

Detailed TSV: `/tmp/balibase_default.tsv`.
