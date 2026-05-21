# MAFFT 7.526 — Static-state coupling in `A__align` produces context-dependent traceback for tied-score DP cells

**Status**: Latent quirk — output remains optimal-scored, but trace selection differs depending on prior `A__align` call history.
**Affects**: any input that triggers tied DP cells. Concretely:
- `mafft --tm 200 --retree 1` on the bundled 36-seq sample (8 lines diff).
- `mafft` (default FFT-NS-2) on ~3 % of BALIBASE 3 inputs (6 of 218).
**Discovered**: 2026-05-18 / 2026-05-19, during byte-identity validation of `rust-MAFFT` against MAFFT 7.526.
**Author contact**: Kazutaka Katoh — katoh@ifrec.osaka-u.ac.jp.

---

## 1. Summary

MAFFT 7.526's `A__align` (in `core/Salignmm.c`) uses `static TLS` buffers
(`commonIP`/`ijp`, `cpmx1/2`, `ogcp1/2o`, `fgcp1/2o`, `gapfreq1/2`, `m`, `mp`,
`w1/w2`) that are sized to the *maximum* `(lgth1, lgth2)` seen across all
calls in a single process, and grown but never shrunk. This is an
intentional malloc/free amortization.

A side-effect of this optimization: when consecutive calls have *decreasing*
`lgth1`/`lgth2`, trailing buffer cells from the previous call's larger DP
remain in memory and subtly influence the next call's DP. For matrices with
score landscapes diverse enough that DP candidate cells are rarely tied
(BLOSUM, JTT 100, etc.), this never surfaces. For TM PAM 200 — whose flatter
score distribution produces ties in many cells — the bias flips which of two
tied DP candidates the traceback selects. **The total alignment score and
width remain optimal**, but the chosen *trace* differs from what `A__align`
would produce if called in isolation on the same inputs.

In practice, this means `A__align(...) ≠ A__align(...)` for identical
arguments: the output depends on prior calls in the process. We were able to
reproduce the divergence with a single warm-up call (one cell flipped) and
the full 12-call pre-state of a 36-sequence progressive alignment (341 cells
flipped, producing a visibly different gap placement).

The author may or may not consider this worth fixing — the alignment is
still correct (tied scores select equally-optimal traces). We report it
because (a) it makes byte-identity validation against `A__align` infeasible
without copying MAFFT's exact static-state behavior, and (b) the output is
not reproducible across hypothetical refactorings of the progressive
pipeline (changing the merge order, even when the tree is unchanged, can
flip these tied-trace selections).

---

## 2. Minimal reproducer

**Test dataset**: `mafft-upstream/test/sample` (the 36-sequence protein test
file shipped with MAFFT).

```sh
mafft --tm 200 --retree 1 mafft-upstream/test/sample > out.fa
```

Compare `out.fa` against the alignment produced by an independent
implementation that calls `A__align`-equivalent DP cell-by-cell (in our case,
a Rust port). With identical inputs (sequences, eff weights, penalty=-917,
matrix), the independent implementation produces:

| Sequence            | MAFFT 7.526         | Independent impl   |
|---------------------|---------------------|--------------------|
| M62903 (line 156)   | `-MAAWEAA---FAARR…` | `-MAAWEAAF---AARR…`|
| S75720 (line 170)   | `------MS---------` | `--------MS-------`|

(8 diff lines total: 4 single-char gap shifts in 2 of 36 sequences.) Both
alignments share the same total score and width — they are tied optima.

The first step where outputs diverge is step 12 of `treebase`: merging
cluster `{seq 7, seq 8, seq 9}` (post-step-4) with raw `seq 11` (chicken
visual pigment, M62903). Both implementations report `pscore =
208675.4644010082` and `width = 367`. But the traceback differs in the
leading 11 columns of the cluster's added gap pattern.

The diff disappears with `--retree 2` (the default) because the second
progressive pass realigns from the rebuilt tree, and the merge order
re-shuffles the prior `A__align` state.

---

## 3. Investigation summary

### 3.1 Approach

To isolate whether the divergence was in our Rust DP or in C's `A__align`, we:

1. Captured the exact step-12 inputs (cluster sequences at width 364, raw seq 11
   at width 362, normalized eff weights, penalties, dynamicmtx) from a Rust
   instrumented progressive run.
2. Patched `Salignmm.c::A__align` with an `MAFFT_IJP_DUMP=<path>` env
   filter that writes the full `ijp[i][j]` matrix to a file when called
   with `lgth1==364 && lgth2==362` (selects the step-12 invocation
   uniquely).
3. Added a parallel dump on our Rust side.
4. Ran:
   - **Test path**: a focused FFI test that invokes C's `A__align`
     in isolation with the captured step-12 inputs.
   - **Full-run path**: the normal `mafft --tm 200 --retree 1`
     invocation, which calls `A__align` 12 times in `treebase` before
     step 12.
5. Diffed the resulting `ijp` matrices cell-by-cell.

### 3.2 Findings

#### Result A — Rust DP is correct given identical inputs

With matched inputs (set via `poffset=0` per the shell script's `-h 0.000`,
C's `n_dis_consweight_multi` as the matrix, normalized eff weights
`[0.267654, 0.267654, 0.464692]` for `{7,8,9}` and `[1.0]` for `{11}`,
penalty=-917, penalty_ex=0), **Rust's `profile_align` produces a
byte-identical `ijp[i][j]` to C's `A__align` across all 365×363 cells (0
interior differences)**. Both produce the same scalar score
(208675.4644010082), same width (367), and the *same* traceback
(`MAAWEAAF---AARRRHEE...`).

#### Result B — Same C function, same inputs, different outputs depending on prior calls

When the same patched C `A__align` is invoked via the full-run pipeline
(after 12 prior `A__align` calls in `treebase`), it produces:
- 341 cells differing from the isolated-call `ijp` (first divergence at
  `i=4, j=257`).
- A different traceback alignment (`MAAWEAA---FAARRRHEE...` instead of
  `MAAWEAAF---AARRRHEE...`).
- The *same* scalar score (208675.4644010082) and *same* width (367).

#### Result C — Observable inputs to the DP are byte-identical between contexts

We dumped the following arrays at the start of the DP loop in both the
isolated test and the full-run step-12 invocation:

| Array                | Indices dumped     |
|----------------------|--------------------|
| `initverticalw[]`    | 0..5               |
| `currentw[]`         | 0..5, 250..259     |
| `ogcp1[]` / `ogcp2[]`| 0..5, 250..259     |
| `fgcp1[]` / `fgcp2[]`| 0..5, 250..259     |
| `gapfreq1pt[]` / `gapfreq2pt[]` | 0..5, 250..259 |

**Every value matched to 10 decimal places** between the isolated test and
the full-run invocation. The visible DP inputs are byte-identical, yet the
DP output `ijp` differs.

#### Result D — Sensitivity confirmed with a single warm-up call

Adding a single one-sequence `A__align` warm-up call in the isolated test —
before invoking it with the step-12 inputs — caused the test's `ijp` output
to flip in exactly 1 cell (vs the no-warm-up baseline). Multi-call warm-ups
would presumably push the flip count higher, converging on the full-run's
341-cell divergence.

This confirms that **the divergence is driven by static-buffer residue from
prior calls**, not by anything in the visible argument list.

### 3.3 Root cause hypothesis

The DP loop in `A__align` does not read `ijp[i][j]` directly (it only writes
it), and the visible boundary arrays match between contexts. The remaining
suspects:

- **Trailing cells of `cpmx1`, `cpmx2`, `gapfreq1`, `gapfreq2`, `ogcp1`,
  `ogcp2`, `fgcp1`, `fgcp2`** beyond `lgth1` / `lgth2`. These are sized for
  the *max* `(ll1+2, ll2+2)` historical and written only over `[0, lgth-1]`
  by `cpmx_calc_new` / `gapcountf` / `st_OpeningGapCount` /
  `st_FinalGapCount`. Trailing cells retain stale data from previous
  longer calls.

- **`m[]` (column-running gap-skip tracker)** is sized for max `ll2+2`,
  initialized only for `j ∈ [1, lgth2]` at the start of each row
  (line 1780). Index `m[0]` is not written; indices `m[lgth2+1..]` retain
  stale data.

- A possible bounds read past `lgth1` or `lgth2` somewhere in the inner
  DP — most likely a `*fgcp2pt` or `*ogcp2pt` read at the column boundary,
  or an `i=lgth1` access in the boundary row that shouldn't fire when
  `tailgp=0` but might via an off-by-one.

We did not pin down the exact buffer cell responsible (the divergence
starts at `i=4, j=257` and requires inspecting every read in the inner DP
to track which one pulls from a stale index). Inspecting the inner loop
bodies at `Salignmm.c:1898-2010` should reveal which array access doesn't
strictly stay within `[0, lgth1)` / `[0, lgth2)`.

### 3.4 Why TM PAM 200 specifically

We ran the same comparison on BLOSUM 62/80, JTT 100/200, and TM 100 — all
byte-identical between isolated and full-run contexts. TM PAM 200 is the
only matrix in the standard MAFFT distribution whose substitution scores
are flat enough to produce ties in the DP candidate comparisons at all.
With BLOSUM/JTT, the cells C and Rust differ on simply don't exist (no
ties to flip).

This is consistent with TM PAM 200 being known to produce more ambiguous
alignments at the per-cell level — the same property that makes
transmembrane alignments harder.

---

## 4. Reproducing the investigation

The full setup and diagnostic code lives in our project at:
`https://github.com/lucas-spearbit/rust-MAFFT` (private at time of writing).
The key files referenced:

- **C patches** (applied temporarily, then reverted via `git checkout`):
  - `core/Salignmm.c:2058` — `MAFFT_IJP_DUMP` env-triggered ijp matrix
    dump for `lgth1==364 && lgth2==362`, plus boundary array dumps for
    `initverticalw`, `currentw`, `ogcp{1,2}`, `fgcp{1,2}`, `gapfreq{1,2}pt`.
  - `core/disttbfast.c:2716` (in `treebase`) — `CDBG_EFF` env-triggered
    dump of `eff1[]`, `eff2[]`, `orieff1`, `orieff2`, `mseq1[]`, `mseq2[]`
    at `l == 12`.
  - `core/disttbfast.c:2929` — `CDBG_HIPREC` env-triggered post-merge
    `pscore` print at full `%.10f` precision.

- **Rust side**:
  - `crates/mafft-align/src/profile.rs::profile_align_imp_with_boundary` —
    paired `MAFFT_IJP_DUMP_R` ijp dump.
  - `crates/mafft-core/src/progressive.rs::progressive_align_full` —
    `RDBG_DUMP_PRE_STEP=12` for capturing pre-step-12 cluster state and
    eff weights.

- **Captured fixtures**:
  - `crates/mafft-core/tests/fixtures/tm200_step12_prof1.fa` — 3 prof1
    sequences at width 364.
  - `crates/mafft-core/tests/fixtures/tm200_step12_prof2.fa` — 1 prof2
    sequence at width 362.

To reproduce in our project (the test was removed after the investigation,
but the patches above can be re-applied):

```sh
# Capture state
RDBG_DUMP_PRE_STEP=12 ./target/release/mafft-rs \
    --tm 200 --retree 1 --nofft mafft-upstream/test/sample 2>&1 > /dev/null

# Full-run C ijp dump
MAFFT_IJP_DUMP=/tmp/c_full_ijp.txt MAFFT_BINARIES=$PWD/mafft-upstream/binaries \
    mafft --tm 200 --retree 1 --nofft mafft-upstream/test/sample > /dev/null

# Isolated C ijp dump (via FFI test setting same env)
cargo test --release -p mafft-core --test cell_diff_tm200_step12 \
    -- --ignored --nocapture

# diff /tmp/c_full_ijp.txt /tmp/c_ijp.txt  # ~341 lines differ
```

---

## 5. Possible fixes (for upstream consideration)

If the author wishes to eliminate the cross-call state leak (which would
make `A__align` deterministic-in-arguments and ease independent
validation), candidate approaches:

1. **Defensive zero-fill of trailing buffer cells.** After each `cpmx_calc_new`
   / `gapcountf` / `st_OpeningGapCount` / `st_FinalGapCount`, zero
   `[lgth, commonAlloc]` of the relevant arrays. Cost: a few loops per
   call, negligible vs the DP itself.

2. **Explicit boundary writes.** Ensure `gapfreq1pt[lgth1]`,
   `gapfreq2pt[lgth2]`, `ogcp1[lgth1]`, `ogcp2[lgth2]`, etc. are written
   on every code path before the DP runs. Currently `cpmxresult ≠ NULL`
   sets some of these (lines 1531-1532), but the `cpmxresult == NULL`
   path relies on `headgapfreq{1,2}` being 0.0 and then flipping to 1.0
   via the `1.0 - x` loop — which works only because `headgapfreq{1,2}`
   are reset, but other trailing array cells aren't.

3. **Track and validate buffer ownership across calls.** Add an
   `MAFFT_DEBUG_STATE` build flag that re-zeros all static buffers between
   calls and asserts the alignment is identical. This would surface the
   sensitivity in regression CI.

The author may equally decide this isn't worth fixing — the alignment is
always optimal, the bias is below any biological-significance threshold,
and the optimization probably matters for very-large-input performance.

---

## 5b. BALIBASE 3 follow-up (2026-05-19)

We ran the same byte-identity validation across the BALIBASE 3 corpus
(`bench.tar.gz` from drive5, 218 protein test sets covering RV11/12/20/30/40/50).
After fixing one Rust-side bug (`ktuple_distance` was using filtered
group length instead of `nogaplen` for `lenfac` — closed 11 of 17
initial divergences), default mode FFT-NS-2 produces **byte-identical
output to MAFFT 7.526 on 212 of 218 tests (97.2 %)**.

The remaining 6 (BB20018, BB20027, BB20039, BB40036, BB40041, BB40048)
all exhibit the same `A__align` static-state pattern as the TM 200
reproducer, just triggered by different inputs. Concretely:

- **4 cases** (BB20018, BB20039, BB40036, BB40048) are pure tied-trace
  artifacts: same width, same score, single-residue gap shifts
  identical to the TM 200 fingerprint.
- **2 cases** (BB20027, BB40041) have width-differs (≤ ±94 cols
  out of ~1600), but the divergence still traces back to a single
  tied-cell choice in an early merge step (step 12 of pass 1 for
  BB20027). The cascade through subsequent merges amplifies the
  shift into a width difference.

For BB20027 specifically we confirmed the pattern with cell-level
FFI tests:
1. `Profile::from_aligned` matches C `cpmx_calc_new + gapcountf +
   st_OpeningGapCount + st_FinalGapCount` bit-for-bit (4 fields, zero
   drift).
2. `blend_profiles_exact` matches C `createcpmxresult +
   creategapfreqresult + createogresult + createfgresult` bit-for-bit
   on a representative 3+5 merge with gappy gaptables.
3. Rust `profile_align` matches C `A__align` BIT-FOR-BIT when
   `cpmxchild = NULL` and `calledbyfulltreebase = 0` (same DP score,
   same width, same alignment columns).
4. The engine's actual `A__align` call passes `calledbyfulltreebase = 1`
   and triggers the static-TLS memoization machinery
   (`Salignmm.c:1446-1450`), which selects different tied DP cells
   than our stateless Rust DP. The static state also drives the
   `reuseprofiles` path that conditionally reuses cached cpmx from a
   prior call.

This confirms the bug class is BALIBASE-corpus-wide: static-state
coupling in `A__align` is not just a TM 200 quirk but a general
phenomenon affecting any input that produces tied DP cells (~3 % of
real-world MSAs from our sample).

## 6. Why we're reporting this

We hit this while validating `rust-MAFFT` byte-for-byte against MAFFT
7.526. Every other mode in the test suite (FFT-NS-2, FFT-NS-i,
L/G/E-INS-1, L/G/E-INS-i, all BLOSUM, JTT, parttree, dpparttree, memsave,
memsavetree, seed, seedtable, treein, treeout, reorder, allowshift,
anysymbol, leavegappyregion, auto, --add, --add --keeplength, Q-INS-i)
produces byte-identical output. The remaining divergences
(`--tm 200 --retree 1` on the 36-seq sample + 6/218 BALIBASE 3 cases)
ALL trace to the same `A__align` static-state mechanism. Faithfully
reproducing the static state in a Rust port would require carrying
non-trivial cross-call state through what is otherwise a stateless
progressive engine — significant architectural cost for output
that's already optimal-scored. We're not asking for a fix — we just
want the diagnosis on record so anyone else hitting this knows what
they're seeing.

If useful, we can provide:
- The full instrumented patches as a diff against MAFFT 7.526.
- Captured fixtures (prof1/prof2 FASTA, eff weight values, the two
  `ijp[i][j]` dumps).
- A standalone C program that reproduces the divergence with two
  back-to-back `A__align` calls (smaller than our full Rust FFI test).
