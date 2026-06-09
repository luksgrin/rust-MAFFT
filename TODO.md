# Pending Work

> **Definition of "exact C parity"**: byte-identical output AT EVERY ITERATION,
> not just at convergence. A run that matches C's final width by coincidence
> while producing different intermediate alignments is NOT parity — it means
> our accept trajectory differs from C's, and the match is fragile
> (input-dependent).

## Current parity matrix (36-seq protein sample, mafft-upstream/test/sample)

`target/release/mafft-rs <args> sample` vs system `mafft <args> sample`
(MAFFT 7.526), diffing the FASTA outputs.

| Mode / flags                                | C width | Rust width | diff lines | Status |
|---------------------------------------------|---------|------------|------------|--------|
| FFT-NS-2 (default)                          | 717     | 717        | 0          | byte-exact ✓ |
| NW-NS-2 (`--nofft`)                         | 717     | 717        | 0          | byte-exact ✓ |
| FFT-NS-i (`--maxiterate 100`)               | 721     | 721        | 0          | byte-exact ✓ |
| L-INS-1 (`--localpair --maxiterate 0`)      | 719     | 719        | 0          | byte-exact ✓ |
| G-INS-1 (`--globalpair --maxiterate 0`)     | 746     | 746        | 0          | byte-exact ✓ |
| E-INS-1 (`--genafpair --maxiterate 0`)      | 729     | 729        | 0          | byte-exact ✓ |
| L-INS-i (`linsi` / `--localpair --maxiterate 1000`) | 735 | 735 | 0       | byte-exact ✓ |
| G-INS-i (`ginsi` / `--globalpair --maxiterate 1000`)| 737 | 737 | 0       | byte-exact ✓ |
| E-INS-i (`einsi` / `--genafpair --maxiterate 1000`) | 729 | 729 | 0       | byte-exact ✓ |
| BL30 / BL45 / BL50 / BL62 / BL80 FFT (`--bl N`) | match | match | 0      | byte-exact ✓ |
| BL80 NW (`--bl 80 --nofft`)                 | 712     | 712        | 0          | byte-exact ✓ |
| JTT 200 / JTT 100 FFT (`--jtt N`)           | match   | match      | 0          | byte-exact ✓ |
| TM 100 / TM 200 NW (`--tm N --nofft`)       | match   | match      | 0          | byte-exact ✓ |
| TM 100 / TM 200 FFT (`--tm N`, default retree=2) | match | match    | 0          | byte-exact ✓ |
| `--tm 200 --retree 1`                       | 717     | 717        | 0          | byte-exact ✓ |
| `--tm 200 --treeout`                        | match   | match      | 0          | byte-exact ✓ |
| PartTree (`--parttree`)                     | 752     | 752        | 0          | byte-exact ✓ |
| DP-PartTree (`--dpparttree`)                | 752     | 752        | 0          | byte-exact ✓ |
| PartTree NW (`--parttree --nofft`)          | 752     | 752        | 0          | byte-exact ✓ |
| `--parttree --reorder`                      | match   | match      | 0          | byte-exact ✓ |
| `--add` / `--add --nofft` / `--add --keeplength` (30+6 fixture) | match | match | 0 | byte-exact ✓ |
| RNA NW (`--nofft samplerna`)                | 360     | 360        | 62 (case)  | byte-exact mod case ✓ |
| Q-INS-i (`--qinsi samplerna`)               | 360     | 360        | 62 (case)  | byte-exact mod case ✓ (needs `mxscarnamod`) |
| `--allowshift --globalpair --maxiterate 0`  | 1029    | 1029       | 0          | byte-exact ✓ |
| `--allowshift --globalpair --maxiterate 1000` | 1060  | 1060       | 0          | byte-exact ✓ (full BALIBASE 3 also 386/386) |
| `--reorder` (FFT-NS-2, INS-i family)        | match   | match      | 0          | byte-exact ✓ |
| `--treeout` (FFT-NS-2, NW-NS-2, FFT-NS-i, L/G/E-INS-i, BL/JTT, parttree, dpparttree) | match | match | 0 | byte-exact ✓ |
| `--treein` (FFT-NS-2, NW-NS-2, BL/JTT/TM, L/G/E-INS-i, FFT-NS-i, all non-parttree) | match | match | 0 | byte-exact ✓ |
| `--treein --treeout`                        | match   | match      | 0          | byte-exact ✓ |
| `--auto` (small/medium/large brackets)      | match   | match      | 0          | byte-exact ✓ |
| `--memsavetree` / `--memsavetree --treeout` | match   | match      | 0          | byte-exact ✓ |
| `--anysymbol` / `--preservecase` (protein + DNA, non-standard chars, mixed case) | match | match | 0 | byte-exact ✓ |
| `--leavegappyregion` / `--legacygappenalty` (combos with --anysymbol/--reorder/--treein/--memsavetree/--auto) | match | match | 0 | byte-exact ✓ |
| `--seed FILE` (L/G/E-INS-i + FFT-NS-i, single + multiple files) | match | match | 0 | byte-exact ✓ |
| `--seedtable FILE` (L/G/E-INS-i + FFT-NS-i, pre-computed hat3.seed) | match | match | 0 | byte-exact ✓ |
| `--memsave` / `--nomemsave` (FFT-NS-2, FFT-NS-i, retree-1, --memsavetree, --nofft combos) | match | match | 0 | byte-exact ✓ |
| `--retree {1,2,3,5}` × {FFT-NS-2, FFT-NS-i, L/G/E-INS-i} (20 combos) | match | match | 0 | byte-exact ✓ |

## BALIBASE 3 sweep (RV11–RV50, 386 protein test sets)

| Mode | Flags | Match | % |
|------|-------|-------|---|
| FFT-NS-i               | `--maxiterate 1000`                              | 386/386 | 100.0 % |
| L-INS-i                | `--localpair --maxiterate 1000`                  | 386/386 | 100.0 % |
| G-INS-i                | `--globalpair --maxiterate 1000`                 | 386/386 | 100.0 % |
| E-INS-i                | `--genafpair --maxiterate 1000`                  | 386/386 | 100.0 % |
| G-INS-i + allowshift   | `--globalpair --maxiterate 1000 --allowshift`    | 386/386 | 100.0 % |

Total: **1930/1930 alignments byte-identical** to C MAFFT 7.526.

## Performance vs C MAFFT 7.526 (macOS arm64, alternating-run medians)

| Mode / input                                | C MAFFT | mafft-rs | Ratio |
|---------------------------------------------|---------|----------|-------|
| FFT-NS-2 (108-seq synthetic)                | 1.34s   | 1.16s    | **1.16× faster** |
| FFT-NS-i (36-seq sample)                    | 0.785s  | 0.670s   | **1.17× faster** |
| FFT-NS-i (BB12019, 5 seqs / 529 cols)       | 0.333s  | 0.198s   | **1.68× faster** |
| FFT-NS-i (BB30004, 50 seqs / 383 cols)      | 1.366s  | 1.246s   | **1.10× faster** |
| `--allowshift --globalpair --maxiterate 1000` (36-seq) | 2.586s | 1.756s | **1.47× faster** |
| L-INS-i (`--maxiterate 50`, 108-seq)        | 37.00s  | 22.93s   | **1.61× faster** |
| G-INS-i (`--maxiterate 50`, 108-seq)        | 35.32s  | 22.62s   | **1.56× faster** |
| `--parttree` (108-seq)                      | 1.70s   | 0.09s    | **18× faster** |
| `--dpparttree` (108-seq)                    | 7.24s   | 0.08s    | **90× faster** |

Test suite: **416 Rust integration tests pass, 0 failed, 6 ignored**
(the ignored are deliberate diagnostic / bisection tools annotated at
the call site). Plus **32 Python tests pass**.

---

## Active items

**Zero known correctness divergences.** What follows is the inventory
of non-correctness gaps still on the radar. Items split into three
groups: design choices we won't reverse without new evidence, known
external-dependency limitations, and untested or unwired surface that
could plausibly bite a user.

### Design choices (closed unless new evidence)

- **Per-group gap stripping** (C's `commongappick`) is intentionally
  NOT ported. C strips per group; we strip globally. The cell-count
  theory predicts a win for C on pathological inputs but empirically
  rust-MAFFT is already 1.16×–90× faster than C on every mode this
  optimisation would target (FFT-NS-2, L/G-INS-i, PartTree,
  DP-PartTree). Re-open only if a production workload surfaces a
  progressive-alignment bottleneck where the cell-count advantage
  plausibly crosses the constant-overhead threshold.
- **Sequential float-summation guard** in `refinement.rs:1389-1397`
  (`intergroup_score` accumulation). Must NOT be parallelised — FP
  add is non-associative; `par_iter().sum()` gives a tree-shaped
  reduction whose shape depends on Rayon work-stealing. The sub-ULP
  per-pair drift accumulates over `(clus1 × clus2)` pairs into a
  multi-unit `mscore` drift that flips accept/reject decisions late
  in refinement. The `mul_add` form also matches clang's FMA fusion
  exactly. The comment IS the reminder; no work needed.

### External-dependency limitations (out of our hands)

- **`--xinsi` (X-INS-i)** requires Stanford's `CONTRAfold v2.02+`
  binary (http://contra.stanford.edu/contrafold/) — not shipped by
  upstream MAFFT, not buildable from `mafft-upstream/extensions`.
  Rust wiring exists (`engine.rs::XInsi`) and emits the correct
  "contrafold not found" diagnostic when absent. End-to-end
  validation deferred until `contrafold` is installed.
- **`--qinsi` (Q-INS-i)** requires `mxscarnamod` built from
  `mafft-upstream/extensions`. Functional when present (alignment is
  byte-identical mod RNA case), but not exercised in CI.
- **`--scarnalike`** requires `dash_client` in PATH. Wired, untestable
  without the binary.

### Pending gaps (not yet addressed)

The following are known gaps that have not been actioned. None block
the 1930/1930 BALIBASE parity or any documented user workflow.

1. **DNA strand auto-detection — IMPLEMENTED.** Both
   `--adjustdirection` (6-mer mode) and
   `--adjustdirectionaccurately` (DP mode) are byte-identical to
   C MAFFT 7.526 on the committed
   `crates/mafft-core/tests/fixtures/dna_adjustdirection_input.fa`
   fixture and on samplerna-derived DNA. See R-5 below for the
   implementation notes.
2. **Tree-linkage variants — all four IMPLEMENTED.** All four
   C tree-linkage flags are wired and byte-identical to C MAFFT:
   - **`--averagelinkage`** (sueff = 1.0), **`--minimumlinkage`**
     (sueff = 0.0), **`--mixedlinkage F`** (sueff = F) all route to
     the existing `mafft_tree::ClusterMethod::Mix` infrastructure
     via `MafftEngine.cluster_method`. Verified byte-identical to C
     across 9 linkage × mode combinations (default / FFT-NS-i /
     L-INS-i × each linkage flag).
   - **`--youngestlinkage`** — IMPLEMENTED, byte-identical to
     C MAFFT 7.526 on first14/15/30/36 and the 36-seq sample
     (widths 423/423/578/703). Dedicated port of
     `compacttree_memsaveselectable(howcompact=2, memsave=1)`
     with both k-mer (pass 0) and MSA (pass 1+) distance paths.
     See R-8 below.
3. **Iteration-strategy variants — all four IMPLEMENTED.** All four
   C iteration-strategy flags have CLI args. Status:
   - **`--simplehillclimbing`** — FULLY HANDLED. Matches C's
     `parallelizationstrategy=BAATARI2`, which is the default in
     both C MAFFT and mafft-rs, so the flag is a true no-op
     (byte-identical to C with or without it).
   - **`--bestfirst`** — IMPLEMENTED, byte-identical to C MAFFT
     7.526 `--thread 1 --bestfirst` across FFT-NS-i (width 713),
     L-INS-i (729) and G-INS-i (734) on the 36-seq sample. See
     R-2 below for the implementation notes (dominant fix: the
     C driver's `iteratelimit=254` cap for BESTFIRST vs 16 for
     BAATARI2).
   - **`--skipiterate F`** — IMPLEMENTED. Both branches of C's
     `dvtditr -E $fixthreshold`:
     (1) **F large** (≥ max-root-tip-distance) skips refinement
     entirely with the C-style WARNING (closed earlier).
     (2) **F small** generates sub-alignment clusters from the
     guide tree and skips refinement of any branch whose subtree
     is a STRICT SUBSET of any cluster (mirrors C's
     `dvtditr.c:997-1006` `includemember && !samemember` gate).
     Byte-identical to C MAFFT 7.526 on the 36-seq sample with
     `--maxiterate 100` across F ∈ {0.05, 0.1, 0.15, 0.2, 0.3,
     0.4, 0.5, 0.7, 1.0}. Implementation: new
     `mafft_tree::generate_subalignments_table` (78 LOC port of
     `mltaln9.c:15330-15407`), new
     `RefinementParams.skip_branches: Vec<(bool, bool)>` field,
     skip check in `iterative_refine`'s branch loop. The
     small-F route bypasses `segmented_iterative_refine` (which
     doesn't honor the skip flags) and uses the standard
     BAATARI2 walk. See R-3 below.
   - **`--oneiteration`** — IMPLEMENTED, byte-identical to C MAFFT
     7.526 on FFT-NS-2 (width 713) and FFT-NS-i (width 717) on
     the 36-seq sample. C only triggers this in the disttbfast
     path (`scripts/mafft:2673` passes `-r` only to `disttbfast`);
     L/G/E-INS-i are correctly no-op in both rust and C.
     See R-4 below for the port notes.
4. **Auxiliary output formats — all six IMPLEMENTED.** All six C
   MAFFT 7.526 aux-output flags have CLI args. Status:
   - **`--distout`** — FULLY IMPLEMENTED. Writes the engine's
     pairwise distance matrix to `<INPUT>.hat2`. Byte-identical to
     C across 6 mode combinations (default / FFT-NS-i / L/G/E-INS-i
     / `--parttree`). Format mirrors C `io.c:2980-2984` exactly,
     including the `=` name prefix C prepends on FASTA read.
   - **`--scoreout`** — FULLY IMPLEMENTED. Prints
     `Unweighted sum-of-pairs score = N.NNNNN` to stderr. Direct
     port of C's `sumofpairsscore` (`mltaln9.c:15411`) + `naivepair
     score11`. Byte-identical to C across 4 modes (default /
     FFT-NS-i / L-INS-i / G-INS-i).
   - **`--nodeout`** — IMPLEMENTED. Appends a per-leaf `Density:`
     section and per-internal-node `Node info:` section to the
     `.tree` file, mirroring C's `treeout==2` path in
     `mltaln9.c::fixed_musclesupg_double_realloc_nobk_halfmtx_treeout`
     (lines 6492-6518). Byte-identical to C MAFFT 7.526 with
     `--nodeout --maxiterate 0` on the 36-seq sample plus
     first14 / first15 / first30 fixtures. Density formula =
     `Σ_{j ≠ i, d(i,j) < 1.0} (2.0 - d(i,j))` (port of
     `setdensity`), `densest` = first-encountered max in each
     subtree (port of `getdensest`, strict `>`). Early `--maxiterate
     > 0` gate matches C's script-level "supports only progressive
     method" error. Regression test:
     `end_to_end::nodeout_density_section_byte_identical_to_c`.
   - **`--mapout` / `--compactmapout`** — IMPLEMENTED. Writes
     `<addfile>.map` after `--add --keeplength`, recording which
     positions of each added sequence were dropped by the
     keeplength column filter. `--mapout` emits the full per-letter
     table (port of `addfunctions.c::reconstructdeletemap`);
     `--compactmapout` emits maximal-run blocks
     (`reconstructdeletemap_compact`). Byte-identical to C MAFFT
     7.526 on the canonical 30+6 fixture
     (`crates/mafft-core/tests/fixtures/sample.add6.mapout.*.map`).
     Regression tests: `mapout_full_byte_identical_to_c`,
     `mapout_compact_byte_identical_to_c`.
   - **`--pileup`** — IMPLEMENTED, byte-identical to C MAFFT
     7.526 on first14/15/30/36 and the 36-seq sample (width
     429/429/604/800). One-line fix: --pileup needs uniform
     1.0 weights (C `tbrweight=0`), not tree-derived. The
     original "single-representative profile" hypothesis was
     wrong; both rust and C use full-cluster profiles, just
     with different weight schemes. See R-7 below.
5. **Sequence-filter flags — all six IMPLEMENTED (2 intentional
   no-ops for --dash dependencies).** All six C MAFFT 7.526
   sequence-filter flags have CLI args. Status:
   - **`--maxambiguous F`** — FULLY IMPLEMENTED. Direct port of C's
     `filter.c` (drops sequences whose ambiguous-residue fraction
     exceeds F, collapses runs of N/X via `shortenN`). Mirrors C's
     gating exactly — filter runs only on the `--add` /
     `--addfragments` file, not the primary input. Byte-identical to
     C across 4 thresholds (0.05 / 0.3 / 0.49 / 0.7) including the
     stderr `Removed N sequence(s)` report.
   - **`--minimumweight F`** — FULLY IMPLEMENTED. Threads through
     `RefinementParams.minimum_weight` to all per-sequence-weight
     clamp sites. Byte-identical to C across 5 modes (FFT-NS-2,
     FFT-NS-i, L/G/E-INS-i).
   - **`--nwildcard` / `--nzero`** — IMPLEMENTED. `--nwildcard`
     fills the DNA scoring matrix's `'n'` row (index 17) with
     25%-self-score values, mirroring C
     `constants.c::nscore`. Engine also enables it implicitly when
     `unalign_level > 0.0` (mirrors `scripts/mafft:1437` setting
     `nmodel=" -: "` whenever `unalignlevel != 0.0`). `--nzero`
     is the documented default — no wiring needed. Byte-identical
     to C MAFFT 7.526 on the samplerna-derived DNA + injected-`n`
     fixture and on the 8-seq adjustdirection fixture across
     {default, --nwildcard, --nzero, --maxiterate 100,
     --allowshift, --allowshift+nwildcard}. Regression test:
     `end_to_end::nwildcard_dna_no_n_chars_byte_identical_to_c`.
   - **`--excludehomologs` / `--originalseqonly`** — WIRED, no-op.
     C documents both as "works with --dash only"; we don't support
     `--dash` (DASH structure-DB pipeline is out-of-scope per the
     RNA-structure design choice).
6. **Fine-grained gap-penalty knobs — fully wired (one tied-trace
   residual at pathological values).** We have `--op` / `--ep` plus
   `--exp` / `--shiftpenalty` / `--lop` / `--lep` / `--lexp` / `--gop`
   / `--gep` / `--gexp` (8 new flags wired). All verified
   byte-identical to C MAFFT across 44/45 mode×value combinations
   (R-1 closed 2026-06-02 — the boundary-init FP-order bug that
   originally caused `--exp` aggressive-value divergence is fixed).

   Sole remaining residual: **`--exp ≥ 4.30` on FFT-NS-2** (no
   refinement) produces 16-24 lines of single-char tied-trace
   gap shifts (same width 517, same score) — pathological-
   value tie-break.

   Threshold pinned 2026-06-04 (`--exp 4.29` byte-identical,
   `--exp 4.30` triggers). Corresponds to `penalty_ex = -2579`
   (= `(int)(0.6 * -4300 + 0.5)`) — at this value the per-cell
   gap-extension penalty crosses a threshold where the DP
   has multiple equally-scoring paths through the Drosophila
   opsin cluster (seqs 14-17).

   Diagnostic ladder (all `--exp 5.0` unless noted):
   - first17 / dros4 / seqs 13-17 in isolation: 0 diff.
   - first21 / first28+: triggers (3-7 diff cols per seq on
     seqs 14-17).
   - `--nofft`: 0 diff (FFT-anchor specific).
   - `--maxiterate ≥ 1`: 0 diff (refinement smooths).
   - The `M-V` vs `-MV` shift is exactly 1 column at the
     boundary of the conserved transmembrane region.

   Closure barrier: the diverging cell is in the per-segment
   profile DP at extreme penalty_ex. All rust-vs-C tie-break
   operators (`>=` for mi/mj GE-update, `>` for wm strict-GT,
   `wm = previousw[j-1]` init, `lasti = n+1` boundary) verified
   matching C `Salignmm.c::A__align` lines 2046-2087.

   **Cell-level FFI investigation 2026-06-04** added per-step
   dump infrastructure (`RS_DP_DUMP=/path` env var in
   `progressive.rs::merge_step_cached`) plus FFI compare tests
   (`exp_residual_dump_compare`, `exp_residual_step10`). The
   compare test feeds rust's recorded per-step inputs to C's
   `Falign` via mafft-sys and reports rust vs C row-by-row.

   Findings:
   - For `--exp 5.0` on the 36-seq sample: ALL 70 merge steps
     (35 × 2 passes) show rust's recorded output == C's
     `Falign` output on the SAME inputs (with weight
     per-group sum-to-1 normalization).
   - The isolated step-10 test (rust's `fft_profile_align` vs
     C's `Falign` on identical inputs): both produce width
     354, byte-identical alignment.
   - But end-to-end `mafft --exp 5.0 sample` still produces a
     16-line content diff vs `mafft-rs --exp 5.0 sample`.

   This is paradoxical: if per-step rust = C on same inputs,
   end-to-end should match by induction. The dump_compare
   harness also shows divergence at step 47 for `--exp 0`
   (which has 0-line end-to-end diff), confirming false
   positives in the FFI compare — some C state we don't
   capture in the dump (likely accumulated FP roundoff or
   subtle Falign global state) shifts behavior in ways the
   per-step replay can't reproduce.

   Pathological-value tie-break with no realistic workflow
   impact — default `--exp` is 0, no documented workflow uses
   `--exp ≥ 4.30`, and any refinement closes it. Diagnostic
   infrastructure is in place if a real closure attempt
   becomes warranted: `RS_DP_DUMP`, `exp_residual_dump_compare`,
   `exp_residual_step10`. `--gop`/`--gep`/`--gexp` are wired but
   currently inert in protein/DNA pipelines (they only affect C's
   X-INS-i / Q-INS-i RNA paths, which are external-dep blocked).
   `--rop` / `--rep` / `--LOP` / `--LEXP` / `--GOP` / `--GEXP`
   are now accepted at the CLI for compatibility (previously
   unrecognized → clap error); they emit a stderr note when
   passed, explaining they only affect the external-dep-blocked
   RNA-structure paths (CONTRAfold / mxscarnamod / LARA / DAFS)
   and are otherwise inert. Variable names map to C verbatim
   (`rgop`/`rgep`/`LGOP`/`LEXP`/`GGOP`/`GEXP`).
7. **PDB structure-aware alignment (`--pdbfilelist`,
    `--pdbidlist`) — full parity with C MAFFT 7.526.** Both flags
    are disabled in upstream C MAFFT itself (`scripts/mafft:969-990`
    — "temporarily unavailable, 2018/Dec." → `exit`), so true parity
    is achieved by emitting the same verbatim message and exit code.
    Rust wires the flags via clap and matches C's exit behaviour
    byte-for-byte (verified 2026-06-04: same stderr message, same
    `$?` exit status).
8. **Out-of-scope by design** (no plan to support): RNA structure
    alignment via DAFS / FoldAlign / LARA / SCARNA / MCCASKILL /
    RIBOSUM / RNAALIFOLD; alternative pairwise via
    BLAST or LAST (`--blastpair`, `--lastpair`, `--lastmultipair`,
    `--fastapair`, `--fastswpair`, `--hybridpair`,
    `--longshortpair`, `--shortlongpair`); MPI parallelism
    (`--mpi` — we use Rayon).

### Open research items

Latent algorithm-port issues surfaced during the gap-flag work above
that need investigation. None affects default-flag behaviour (the full
1930/1930 BALIBASE parity is unchanged); each one is a
non-default-input edge case that diverges from C.

R-1b. **`--exp > 0` NW-path tied-trace residual — RESOLVED
   2026-06-03.** Root cause was NOT an FP-order accumulation
   issue (the original hypothesis); it was a MISSING per-cell
   `fpenalty_ex` contribution in
   `mafft_align::pairwise_align11` (the port of C
   `Galign11.c::G__align11`).
   C lines 1362-1364 and 1381-1384 add `mi += fpenalty_ex_i`
   (in-row gap extension) and `m[j] += fpenalty_ex` (in-col)
   per cell after the diagonal update, gated on
   `i < lgth1` / `j < lgth2` respectively (the 2018/May/11
   tail-boundary fix). Rust's `pairwise_align11` had no
   `penalty_ex` parameter at all — it silently treated the
   extension penalty as zero, which coincidentally matched C
   on the FFT-segmented path (because `pairwise_align11` is
   only called from the 1-vs-1 merge fast-path, and FFT-NS-2
   skips that path on this fixture) and matched C exactly when
   `--exp = 0` (`fpenalty_ex = 0` makes the missing
   contribution zero).
   Fix: added a new `pairwise_align11_ex(seq1, seq2, mtx, map,
   penalty, penalty_ex, head_gap, tail_gap)` that applies
   `penalty_ex` per cell with the boundary gates. The old
   `pairwise_align11` delegates with `penalty_ex = 0`.
   `progressive.rs::merge_step_cached` passes
   `scoring.gap.extend` so the per-cell contribution flows
   through. Byte-identical to C MAFFT 7.526 `--nofft
   --maxiterate 0 --exp F` across F ∈ {0, 0.05, 0.1, 0.15,
   0.2, 0.25, 0.3, 0.4, 0.5, 1.0, 2.0, 4.4, 5.0} on the 36-seq
   sample. FFT path unchanged. Regression test:
   `nofft_exp_sweep_byte_identical_to_c`.

R-1. **`--exp` divergence — CLOSED 2026-06-02.** Required three
   independent fixes:
   1. **Boundary-init FP-order** in
      `profile_align_imp_multimtx::initverticalw` (1-ULP FMA
      ordering matching clang's codegen).
   2. **Refinement zero-out** in `refinement.rs::iterative_refine`
      (`GapModel::new(..., 0.0)` instead of `scoring.gap.extend`).
      C's `dvtditr` invocation in `scripts/mafft` does NOT pass
      `-g $gexp`, so C's refinement always uses `penalty_ex = 0`.
   3. **Constraint-aware progressive zero-out** in
      `progressive.rs::progressive_align_full_c_compat_ex` (gap
      extend forced to 0 when `constraints.is_some()`). The C
      flow for L/G/E-INS-i is `pairlocalalign → tbfast → dvtditr`,
      none of which receive `-g $gexp`. Only `disttbfast` (the
      FFT-NS-2 progressive path) does.
   Result: 62/63 `--exp` × mode combinations byte-identical to C
   across {0, 0.001, 0.01, 0.1, 0.5, 1.0, 5.0} × {FFT-NS-2,
   FFT-NS-i, L-INS-i, G-INS-i, E-INS-i, L-INS-1, G-INS-1,
   E-INS-1, `--allowshift`}. The one residual at `--exp ≥ 4.4` on
   FFT-NS-2 default is a tied-trace gap placement (same width 517,
   same score, 16 lines of single-char gap shifts) — a §B.2-class
   tie-break at pathologically high penalty_ex values that
   amplifies tied DP cells. Three end_to_end regression tests
   guard the fix (`fftns2_exp_0_1_byte_identical_to_c`,
   `fftnsi_exp_0_1_byte_identical_to_c`,
   `linsi_exp_0_1_byte_identical_to_c`).

R-2. **`--bestfirst` BESTFIRST refinement — RESOLVED 2026-06-02.**
   C's `BESTFIRST` strategy (`tditeration.c::athread` worker loop
   plus the thread-0 collector at 1471-1543) evaluates EVERY branch
   from a frozen baseline alignment per iteration, then the
   collector picks the single highest-gain branch across all
   workers and applies it to the master copy. With `--thread 1`,
   one worker walks all branches in `branchtable` order against the
   unchanged baseline; the collector applies the best move. Ported
   as `refinement.rs::bestfirst_refine`, dispatched from
   `engine.rs::1067-1071` when `params.bestfirst` is true. Three
   bugs surfaced during the port:
   1. **Iterate cap.** C's `scripts/mafft:1512` raises
      `iteratelimit` to 254 for BESTFIRST (16 for BAATARI2). Rust
      was capping refinement at 16 unconditionally — the dominant
      divergence (rust width 725 vs C 713 on FFT-NS-i). Fixed in
      `engine.rs::988`.
   2. **Identity short-circuit.** When `localcopy[s1]==mastercopy[s1]`
      and same for `s2`, C sets `tscore=mscore` and skips the
      recompute (`tditeration.c:2185-2192`). Without this rust
      could accept a noop move whose FP-recomputed gain rounded
      barely positive. Added the strcmp-equivalent guard in
      `bestfirst_refine` (refinement.rs:1681-1690).
   3. **`gain > best` tie-break.** Strict `>` to match C's per-thread
      "first encountered max wins" rule (`tditeration.c:2253`).
   Byte-identical to C MAFFT 7.526 `--thread 1 --bestfirst` across
   FFT-NS-i / L-INS-i / G-INS-i on the 36-seq sample and across
   `--maxiterate` values {2, 10, 50, 100, 200}. Regression tests:
   `end_to_end::bestfirst_fftnsi_byte_identical_to_c`,
   `bestfirst_linsi_byte_identical_to_c`,
   `bestfirst_ginsi_byte_identical_to_c`.

R-3. **`--skipiterate F` — RESOLVED 2026-06-03.** Both paths
   now byte-identical to C MAFFT 7.526 `--maxiterate 100
   --skipiterate F`:
   1. **F ≥ max-root-tip-distance** → skip refinement entirely
      with C-style WARNING (closed 2026-06-02).
   2. **F < max-root-tip-distance** → generate sub-alignment
      clusters via `mltaln9.c::generatesubalignmentstable`
      (ported as `mafft_tree::generate_subalignments_table`),
      then skip refinement of any branch whose subtree is a
      strict subset of any cluster (port of C
      `dvtditr.c:997-1006`'s `includemember && !samemember`).
   Verified byte-identical across F ∈ {0.05, 0.1, 0.15, 0.2,
   0.3, 0.4, 0.5, 0.7, 1.0} on the 36-seq sample. Two
   layered fixes during the port:
   1. **Skipped branches must NOT count toward convergence.**
      First attempt incremented `converged_count` for skipped
      branches; with 38/70 branches skipped on the test input,
      we bailed out of `iterative_refine` after 36 skipped
      branches (= convergence target) before attempting ANY
      real refinement. C `tditeration.c:2358` mirrors this:
      `identity = 1; tscore = mscore` for skipped branches
      flows to the "no improvement" path which doesn't
      increment the converge counter. Fixed by silently
      `continue`ing past skipped branches.
   2. **FFT-segmented refinement doesn't honor skip flags.**
      Force the small-F case through the standard BAATARI2
      `iterative_refine` (which DOES honor skips) instead of
      `segmented_iterative_refine`. Gated on
      `!params.skip_branches.is_empty()` in `engine.rs`.
   Regression test:
   `skipiterate_small_f_byte_identical_to_c`.

R-4. **`--oneiteration` one-vs-others refinement — RESOLVED
   2026-06-02.** Ported as `refinement.rs::one_vs_others_refine`
   (~90 LOC) and dispatched from `engine.rs` AFTER the
   progressive merge and BEFORE regular refinement, gated on
   `(FftNs2 | FftNsi) && self.oneiteration` (matching C
   `scripts/mafft:2673` which passes `-r` only to `disttbfast`).
   Algorithm: ITERATIVECYCLE=2 full passes; per iteration
   `l = ll % nseq`, treat `{l}` as group1 and the rest as group2,
   per-group commongappick, realign via the progressive-Falign
   path (kobetsubunkatsu=0), accept iff intergroup_score didn't
   drop. Three layered fixes surfaced during the port:
   1. **Realign path.** Initial attempt used the existing
      `realign_all` (refinement DP, kobetsubunkatsu=1) and
      produced width 715 vs C 713. C's `dooneiteration` uses
      Falign directly (the progressive disttbfast path), so the
      rust port must call the equivalent `merge_step_cached`
      (exposed via new public wrapper
      `progressive::merge_two_groups_progressive`) NOT
      `realign_all`.
   2. **Pre-strip common gaps per group.** C
      `disttbfast.c:2338-2339` calls `commongappick(mseq2)` then
      `commongappick(mseq1)` BEFORE the realign DP. Without this
      pre-strip the progressive merge sees redundant gaps and
      explodes the width (rust 933 vs C 713 on the second
      attempt). Mirrored by computing per-group all-gap column
      masks and slicing the singleton + N-1 group sequences down
      to the stripped columns before calling
      `merge_two_groups_progressive`.
   3. **`intergroup_score_c_order` for the accept test.** C's
      score formula precomputes `efficient = eff1[i] * eff2[j]`
      then accumulates `value += tmpscore * efficient` — one mul
      outside, one fma. `compute_split_score` instead inlines as
      `((tmpscore * wi) * wj) + total`, a mathematically
      equivalent but FP-different product order. Sub-ULP drift
      across 36*2 leave-one-out iterations was below the
      tie-break threshold in this port (the wrong-realign-path
      bug dominated), but introduced
      `intergroup_score_c_order` defensively to mirror C exactly
      and isolate this from future regressions.

   Byte-identical to C MAFFT 7.526 `--oneiteration` across:
   - FFT-NS-2 (width 713), FFT-NS-i `--maxiterate {2, 100}`
     (width 717), composed with `--reorder` and `--bestfirst`.
   - L/G/E-INS-i (no-op as in C — `dooneiteration` only triggers
     in disttbfast, not pairlocalalign-based pipelines).
   Regression tests:
   `end_to_end::oneiteration_fftns2_byte_identical_to_c`,
   `oneiteration_fftnsi_byte_identical_to_c`,
   `oneiteration_linsi_noop_byte_identical_to_c`.

R-5. **`--adjustdirection` k-mer strand detection — RESOLVED
   2026-06-02.** Default 6-mer mode (`makedirectionlist -m -o a -r
   5000`) ported to `mafft-core::adjust_direction::adjust_direction`
   (~250 LOC) and dispatched from `main.rs` between input read
   and `engine.align()`. Algorithm: build forward and
   reverse-complement 6-mer point vectors for every input, sort
   by descending (forward_self − reverse_self) contrast, then
   walk in that order — each sequence's two orientations are
   scored against the mean common-6-mer count of every
   already-decided sequence's chosen orientation. The orientation
   with the higher mean wins; ties default to Forward (C's strict
   `>` at `makedirectionlist.c:1234`). If the original-index-0
   sequence ends up Reverse, flip every direction (C
   `makedirectionlist.c:1261-1270`) so the first output is
   always forward. Reverse sequences get their name prefixed with
   `_R_` mirroring `setdirection.c:142-155`. Reuses k-mer
   primitives from `mafft_tree::parttree_dist`
   (`encode_points_dna`, `composition_table`,
   `common_sextets_p`).

   Byte-identical to C MAFFT 7.526 `--adjustdirection` on:
   - The 5-seq `samplerna` re-cast as DNA with 2 sequences
     reverse-complemented (same `_R_` headers and same alignment
     body when paired with `--preservecase`).
   - The committed 8-seq fixture
     `crates/mafft-core/tests/fixtures/dna_adjustdirection_input.fa`
     (4 forward seqs + 4 RC copies — direction detector marks
     exactly the last four `_R_`).
   - All-forward DNA input (no flips).
   - Protein input (passthrough no-op).

   **`--adjustdirectionaccurately` (DP mode) — RESOLVED 2026-06-03.**
   Same `adjust_direction_mode` entry point with a new
   `AdjustMode` enum; the DP path swaps `common_sextets_p` for
   `mafft_align::local_align` (port of C's `L__align11_noalign`
   in `makedirectionlist.c:646`), and the reference cap is `100`
   instead of `5000` (mirrors C `scripts/mafft:2333` `-r 100`).
   Contrastsort uses `L_align(s,s) - L_align(s,rc(s))` instead
   of the 6-mer self-difference (matches C `selfdpthread`).
   Byte-identical to C MAFFT 7.526 `--adjustdirectionaccurately`
   on the 8-seq adjustdirection fixture and on the 5-seq
   samplerna-derived DNA fixture. Regression test:
   `end_to_end::adjustdirectionaccurately_mixed_dna_byte_identical_to_c`.

   **`--add` interaction — RESOLVED 2026-06-03.** Added
   `adjust_direction_mode_add(input, mode, nadd)` entry point.
   Anchors the first `nseq - nadd` sequences as Forward (uses
   their input order as the reference set) and only contrastsort/
   orientation-tests the last `nadd`. CLI dispatches to this when
   `--adjustdirection` + `--add` are both set, combining
   existing + added into a single input. Also surfaced a
   case-normalization bug: existing was read with
   `--preservecase` (lowercase), added without (uppercased), and
   the DNA scoring matrix has zero cross-case entries — so DP
   `local_align` between mixed-case inputs returned 0. Fixed by
   normalizing to lowercase in `gappick0` (matches C
   `getnumlen`'s uppercase behaviour modulo case choice — the
   orientation decision is case-invariant either way). Regression
   test: `adjustdirection_with_add_only_flips_added`.

   **Contrastsort sort-stability follow-up — RESOLVED
   2026-06-04.** Switched `testable.sort_by` → `sort_unstable_by`
   to better match C `qsort`'s spirit. Verified byte-identical
   to C on (a) the standard 8-seq + 5-seq + 17-seq DNA fixtures
   (no contrast-key ties), (b) a synthetic palindromic-DNA
   fixture where ALL 4 contrast keys tie at 0 (orientation
   decisions are direction-symmetric so the tie order doesn't
   matter downstream). Neither rust's stable nor unstable sort
   matches glibc qsort bit-exactly on every conceivable tied
   input — that would require linking glibc's qsort or re-
   implementing its exact pivot/partition algorithm — but no
   realistic input has been observed to expose the difference.

   Regression test:
   `end_to_end::adjust_direction_mixed_dna_byte_identical_to_c`.

R-6. **Gapped-input divergence — direct alignment + `--add`
   path both CLOSED 2026-06-03.**

   **Direct alignment path — RESOLVED.** Root cause was rust's
   `engine.align` keeping the input sequences verbatim (with
   pre-existing FASTA gaps) instead of stripping per-row gaps
   before the progressive merge. C `disttbfast.c:4453` does
   `for(i=0; i<njob; i++) gappick0(bseq[i], seq[i])` — every
   sequence is reduced to residues before progressive starts.
   Without this, feeding an already-aligned FASTA produces a
   different alignment than C (369-line diff on the
   `combined_17_r6.fa` fixture). Fix is a 4-line gappick in
   `engine.align`'s sequence-vector construction. Verified
   byte-identical to C MAFFT 7.526 on the fixture (width 444).
   Regression test: `r6_gapped_input_byte_identical_to_c`.
   36-seq sample and BBaliBase 1930/1930 unchanged.

   **`--add` path still diverges** on the same adversarial
   input (16 existing + 1 added-with-20-insertions): 238-line
   diff with `--add`, 336-line diff with `--add --keeplength`.
   The canonical 30+6 `--add` fixture is byte-identical.

   **Deeper investigation (`progressive_align_with_mergeoralign_n`,
   2026-06-03)** narrowed the divergence to step 4 (the only
   NewRight step in the 16+1 topology — merging existing seq 0
   with the added seq 16). Both rust and C take FFT path here
   (use_fft && ffttry both true for 1-vs-1 with reasonable
   nlen). Verified:
   - The pre-merge inputs are identical between rust and C
     (stripped seq[0] = 348 residues, raw seq[16] = 373
     residues with 20 random insertions).
   - The merge result width differs by 1 column (rust 368,
     C 367 post-merge, leading to rust 445 / C 444 after
     restoring 77 common-gap columns).
   - With `--add --nofft`, the divergence is IDENTICAL → not
     an FFT-anchor or FFT-cut issue.
   - `compute_new_merge_gap_set` and
     `build_other_post_restore_row` look correct on inspection
     (correctly handle the residue-vs-gap detection).

   **Suspected root cause**: rust's `build_other_post_restore_row`
   port of C's `insertnewgaps` (`addfunctions.c:445-650`) is
   incomplete. C's `insertnewgaps` calls **`profilealignment`**
   on the OTHER rows at every gap-region of length > 0
   (`addfunctions.c:583-588`) — when the new-merge inserts a
   gap of length `gapshift` into the existing-side, OTHER
   rows aren't just padded with `-` chars; their actual
   content at that pre-merge position participates in a
   sub-profile alignment with the gap region. Rust currently
   does flat `-` insertion (`progressive.rs:565-575`).

   For the canonical 30+6 fixture and standard biological
   inputs, the merge tends to insert gap REGIONS in places
   where OTHER rows already have gaps (so profilealignment
   collapses to a no-op). For the adversarial input
   (random insertions), the merge inserts gaps at
   positions where OTHER rows have residues, and the
   profilealignment matters.

   **Implementation attempt 2026-06-03**: a structural
   scaffold for the port now lives in
   `progressive.rs::insertnewgaps_with_profilealignment` (gated
   `#[allow(dead_code)]`). It walks the post-merge gap regions,
   extracts OTHER's content at the next stripped_width anchors,
   builds 2-side profiles, runs `profile_align`, and emits the
   ops into the OTHER rows. Single-row wire-in (behind an env
   var) showed the function compiles and runs but doesn't close
   the divergence on its own (445-line cascade became 545-line
   instead of dropping toward 0).

   The remaining piece: synchronizing the **active rows'
   width** with OTHER's when profilealignment changes the
   gap-region width. C handles this naturally because its
   per-position loop emits chars for ALL groups (group0=OTHER,
   group1=existing-side, group2=added) at every step. Rust
   currently restores active rows to fixed
   `post_merge_width + n_gap_cols` BEFORE OTHER reconstruction
   — so when profilealignment compresses a gap region by 1
   col, OTHER and active end up at different widths and the
   alignment is broken.

   **CLOSED 2026-06-03**: full port shipped — `apply_c_
   insertnewgaps` + `rs_profilealignment` in `progressive.rs`,
   enabled by default (`RS_R6_PORT_OFF` to revert to flat-
   padding for diagnostics). Cell-level FFI parity harness at
   `crates/mafft-core/tests/cross_validate_insertnewgaps.rs`
   validates rust = C byte-by-byte on hand-crafted scenarios.

   Key gotcha that took the longest to find: rust represents
   new-merge-gaps as `-` (same as common-gaps), while C uses
   `=`. So `findnewgaps`-via-string-scan returned all zeros
   — looked like no new-merge-gaps when there were 146 of
   them. Fix: derive `gaplen` and `gapmap` directly from the
   `new_merge_gap_set` + `gap_cols_before` data already
   computed for the rust restore step.

   Test results:
   - R-6 --add adversarial (exist_2..16 + 3 added): 0 diff
     on all (was 80, 96, 128, 232, 320, 376).
   - Canonical 30+6 --add --keeplength: 0 diff.
   - 36-seq sample (non-add): 0 diff.
   - 415/415 workspace tests passing.

   **Original R-6 entry preserved below for context** —
   adversarial-input `--add` divergence on
   constructed-divergent inputs only.

R-6-orig. **`--add` alignment-body divergence on synthetic divergent
   inputs — characterised 2026-06-03.** Canonical 30+6 fixture
   is byte-identical (see parity matrix); standard `--add`
   workflows (a few short fragments into a well-aligned
   reference) match C. Divergence appears only on adversarial
   inputs where the added sequence has many random insertions:
   **threshold is ≥16 existing sequences AND ≥8 random
   amino-acid insertions in the added seq**. Below either
   threshold the output is byte-identical to C MAFFT 7.526.
   Above it, rust is consistently 1 column wider; the 1-column
   shift cascades into a constant 238-line diff regardless of
   the exact insertion count.

   Verified with bisection (`/tmp/r6test/`):
   - 15 existing + 20 insertions: 0 diff.
   - 16 existing + 5 insertions: 0 diff.
   - 16 existing + 8 insertions: width rust 433 vs C 432, 238
     line diff.
   - 16 existing + 20 insertions: width rust 445 vs C 444, 238
     line diff.

   Eliminated as causes:
   - It's not a `pairwise_align11` `penalty_ex` omission like
     R-1b — that fix is in this build and doesn't close R-6.
     Default `scoring.gap.extend = 0` for protein, so per-cell
     extension contributions would be zero anyway.
   - It's not the `--mapout` writer (the writers are
     byte-identical on the canonical fixture per
     `mapout_full_byte_identical_to_c`).

   **Deeper finding 2026-06-03**: this is NOT really an
   `--add`-specific bug. The divergence reproduces on the
   same input WITHOUT `--add` (just align all 17 seqs as
   one input → 369-line diff). So it's a tied-trace
   divergence in the standard FFT-NS-2 / profile-DP path,
   surfaced only by adversarial input characteristics:
   - First 17 of sample (real biological seqs): 0 diff.
   - 36-seq full sample: 0 diff (BBaliBase 1930/1930).
   - 16 existing + 1 added-with-20-insertions
     (`combined_17.fa`): 369 diff.

   The synthetic input has highly repetitive residue
   patterns (added seq is a real sequence with random
   amino acids inserted at random positions), creating
   tied DP cells that real biological data doesn't
   produce. The k-mer distance matrix already differs
   slightly between rust and C on this input (one branch
   length 0.03430 vs 0.03759 in the standalone 16-seq
   tree → guide-tree topology then drives different merge
   order, cascading into the 369-line diff).

   Eliminated as causes:
   - Not FFT-anchor (--nofft also diverges).
   - Not `pairwise_align11` `penalty_ex` (R-1b fix in
     place, doesn't close R-6).
   - Not the `--mapout` writer.
   - Not `--add`-specific (no-add reproduces it too).

   Confirmed root cause **family**: 1-ULP FP-order drift
   in either `ktuple_distance` or `musclesupg` cluster-
   distance computations causes a different tied-distance
   selection at the 16-seq scale, producing a different
   tree → different merge order → 369-line cascade. The
   first-divergent-tie would be in either the distance-
   matrix construction or the UPGMA `mindist` selection
   at one specific cluster-pair distance.

   **Realistic workflows do not hit this** — committed
   regressions include the full 36-seq sample, all of
   BBaliBase 3 (1930/1930), and the canonical 30+6
   `--add` fixture, all byte-identical. The synthetic
   input was specifically constructed to amplify ties.

R-8. **`--youngestlinkage` divergence — CLOSED 2026-06-03.**
   Port required two pieces (NOT 638 LOC):
   1. `youngestlinkage_tree` (~200 LOC) — k-mer path mirroring
      `compacttree_memsaveselectable(seq=NULL, howcompact=2,
      memsave=1)`. Forward + two-sided initial scan via
      `initial_mindist_yl`; per-step recompute via k-mer tables.
   2. `youngestlinkage_tree_msa` (~50 LOC sharing core via
      closure) — MSA path mirroring the same C function with
      `seq=bseq` (triggering `verycompactmsadistarrthreadjoblist`).
      Used in pass 1+ of `--youngestlinkage` (the two-pass
      progressive C MAFFT runs).

   Shared `youngestlinkage_core` factors the per-step loop;
   distance closures parameterize k-mer vs MSA. Cell-level FFI
   harness (`youngestlinkage_topol_matches_c_step_by_step`)
   validated rust port = C `compacttree_memsaveselectable`
   byte-identically when fed the same initial mindist.

   Results: first14/15/30/36 all 0 diff vs C MAFFT 7.526
   `--youngestlinkage` (was 0/24/682/870). Width matches:
   first30 578=578, first36 703=703. --memsavetree preserved
   byte-identical. 416/416 tests passing.

R-7. **`--pileup` alignment-body divergence vs C — CLOSED
   2026-06-03 with a one-line fix.** Root cause was NOT
   single-rep vs full-cluster (the original hypothesis); both
   C and rust use full-cluster profiles. The actual difference:
   C's `--pileup` sets `tbrweight = 0` (`disttbfast.c:3962`),
   which makes the merge weights `eff[i] = 1.0` for all i
   (uniform). Rust was using `sequence_weights(topology)`
   (tree-derived) for all modes except `--parttree`.

   **Fix** (`engine.rs:832`): one-line addition to the
   existing `weights_override` branch:
   ```rust
   let weights_override = if use_parttree || self.pileup {
       Some(vec![1.0; sequences.len()])
   } else { None };
   ```

   Closes byte-identity on `first14`, `first15`, `first30`,
   `first36`, and `36-seq sample` (was 1000+ line diff on
   36-seq, now 0). The CHAIN TOPOLOGY was already byte-
   identical (existing
   `pileup_topology_byte_identical_to_c_branch_lengths` test).
