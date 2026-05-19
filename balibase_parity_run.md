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

Spent a long session narrowing BB20027 (29 seqs, +28 column delta).
Did NOT close it but isolated the divergence to a specific code path:

- `--retree 1` (single-pass progressive): **byte-identical** to C
  (both width 1612).
- `--retree 2` (default, second progressive pass on rebuilt tree):
  diverges (C width 1602 vs Rust 1630).
- `--nofft` does not affect the divergence (FFT anchor placement is
  not the cause).
- C and Rust's rebuilt tree (`--treeout` after `--retree 2`) is
  **byte-identical** Newick — same topology, same branch lengths.
- Sequence weights derived from the rebuilt tree therefore match
  (since `sequence_weights` is deterministic given the tree).
- Inputs to pass 1: raw input sequences, byte-identical in both
  implementations (we use `sequences.clone()` each pass, matching
  C's `gappick0(bseq, seq)`).
- `--add` test (give both implementations C's 28-way intermediate
  alignment + the 1p8j_A sequence, ask them to add it): byte-identical
  output. So the *final* merge step's profile-vs-sequence DP is fine.
- Therefore the divergence is in some earlier merge step of pass 1.
  The 28-way intermediate profile that Rust pass-1 builds is subtly
  different from C's, despite (raw input, tree, weights) being
  byte-identical.
- BB20027 has NO non-standard residues (no X / `.` / `J`), so it's
  not another `nogaplen` lenfac case.
- DP TLS pool reset between passes (mafft-align/profile.rs's
  `reset_dp_pools`) does not change the output. Cross-pass buffer
  pollution is not the cause.
- The result is deterministic across 3 runs (md5sum identical).

**What's left to check (would require multi-hour C instrumentation)**:
1. Instrument `disttbfast.c::treebase` to log per-merge-step input
   group widths, output width, and DP score. Diff Rust pass-1's
   step trace against C's step-by-step.
2. Suspect: the per-step `eff1`/`eff2` weight blending in
   `progressive.rs::merge_step_cached` may diverge from C's
   `fastconjuction_noname` accumulation for non-uniform group
   sizes. Pass 0 happens to not trigger it; pass 1's tree does.
3. Suspect: the gap-stripping inside `Profile::from_aligned`
   when rebuilding a profile from aligned-with-gaps sequences
   (pass 0 starts from gap-free input; pass 1 may re-encounter
   pre-built profiles with internal gaps from earlier merges).

## Next investigation step

Pick BB40041 next — the other remaining WIDTH-differs case. Same
diagnostic flow: `--retree 1` to see if it's also a pass-1-only
divergence, `--nofft` to rule out FFT. If it shares the BB20027
profile, both are likely a single root cause in pass-1 progressive
merge accumulation.

## Reproduction

```bash
cargo build --release
scripts/balibase_parity.py /tmp/balibase/bench1.0/bali3/in \
    --pattern "BB[0-9]*" --modes "" --output /tmp/balibase_default.tsv
```

Detailed TSV: `/tmp/balibase_default.tsv`.
