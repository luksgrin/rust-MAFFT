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

Test suite: **364 Rust integration tests pass, 0 failed, 4 ignored** (the
4 ignored are deliberate diagnostic / bisection tools annotated at the
call site). Plus **32 Python tests pass**.

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

1. **DNA strand auto-detection — wired as no-op stub.** Both
   `--adjustdirection` and `--adjustdirectionaccurately` accept at
   the CLI with an honest "not yet implemented" stderr note. The
   k-mer-based detection algorithm itself (port of C's
   `makedirectionlist.c`, ~1300 LOC) is tracked as research item
   R-5 below. Affects DNA users only; protein workflows
   unaffected.
2. **Tree-linkage variants — partially exposed.** Three of the four
   C tree-linkage flags are now wired and byte-identical to C MAFFT:
   - **`--averagelinkage`** (sueff = 1.0), **`--minimumlinkage`**
     (sueff = 0.0), **`--mixedlinkage F`** (sueff = F) all route to
     the existing `mafft_tree::ClusterMethod::Mix` infrastructure
     via `MafftEngine.cluster_method`. Verified byte-identical to C
     across 9 linkage × mode combinations (default / FFT-NS-i /
     L-INS-i × each linkage flag).
   - **`--youngestlinkage`** is a separate algorithm in C
     (`treeext="youngestlinkage"`, not a sueff value); accepted at
     the CLI as a no-op with a warning. Porting the youngest-
     linkage algorithm is an open item.
3. **Iteration-strategy variants — partially exposed.** All four
   C iteration-strategy flags now have CLI args. Status:
   - **`--simplehillclimbing`** — FULLY HANDLED. Matches C's
     `parallelizationstrategy=BAATARI2`, which is the default in
     both C MAFFT and mafft-rs, so the flag is a true no-op
     (byte-identical to C with or without it).
   - **`--bestfirst`** — WIRED, no-op. C's `BESTFIRST` strategy
     refines branches in score-order (best-first) rather than the
     BAATARI2 hill-climbing default. Porting requires restructuring
     the refinement loop to pre-score all branches per iteration.
     See R-2 below.
   - **`--skipiterate F`** — WIRED, no-op. C's `dvtditr -E
     $fixthreshold` skips branches whose distance-from-tip exceeds
     F. Porting requires adding the branch-skip gate to the
     refinement loop. See R-3 below.
   - **`--oneiteration`** — WIRED, no-op. C's `disttbfast -r`
     forces a single iteration of distance recomputation. Our
     distance recomputation is already gated by `--retree`; the
     interaction needs investigation before mapping. See R-4 below.
4. **Auxiliary output formats — partially exposed.** All six C
   MAFFT 7.526 aux-output flags now have CLI args. Status:
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
   - **`--nodeout`** — WIRED, partial. Currently behaves as
     `--treeout`. The C `--nodeout` additionally appends a per-leaf
     `Density:` section to the tree file; that section is not yet
     emitted. Open item.
   - **`--mapout` / `--compactmapout`** — WIRED, no-op. The flag is
     accepted but no `.map` file is written. C tracks
     input-to-output column mapping during the `--add` /
     `--addfragments` pipeline; that plumbing is not yet exposed.
     Open item.
   - **`--pileup`** — WIRED, routes through standard FASTA. The C
     `--pileup` invokes a DIFFERENT alignment strategy
     (`treeext="pileup"` in `scripts/mafft`), not just a different
     format. Porting the pileup strategy is a separate item; out of
     scope here.
5. **Sequence-filter flags — partially exposed.** All six C MAFFT
   7.526 sequence-filter flags now have CLI args. Status:
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
   - **`--nwildcard` / `--nzero`** — WIRED, no-op. Accepted at the
     CLI but the N-row substitution-matrix tweak (C's `-:` / blank
     `nmodel` in `constants()`) is not yet applied. Affects DNA
     workflows only. Open item.
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
   Sole remaining residual: `--exp ≥ 4.5` on FFT-NS-2 (no refinement)
   produces 16 lines of single-char tied-trace gap shifts (same
   width, same score) — pathological-value tie-break, out-of-scope
   for realistic workflows. `--gop`/`--gep`/`--gexp` are wired but
   currently inert in protein/DNA pipelines (they only affect C's
   X-INS-i / Q-INS-i RNA paths, which are external-dep blocked).
   Not implemented: `--rop` / `--rep` (RNA-only); `--LOP` /
   `--LEXP` / `--GOP` / `--GEXP` (LARA RNA only).
7. **Out-of-scope by design** (no plan to support): RNA structure
    alignment via DAFS / FoldAlign / LARA / SCARNA / MCCASKILL /
    RIBOSUM / RNAALIFOLD; PDB structure-aware alignment
    (`--pdbfilelist`, `--pdbidlist`); alternative pairwise via
    BLAST or LAST (`--blastpair`, `--lastpair`, `--lastmultipair`,
    `--fastapair`, `--fastswpair`, `--hybridpair`,
    `--longshortpair`, `--shortlongpair`); MPI parallelism
    (`--mpi` — we use Rayon).

### Open research items

Latent algorithm-port issues surfaced during the gap-flag work above
that need investigation. None affects default-flag behaviour (the full
1930/1930 BALIBASE parity is unchanged); each one is a
non-default-input edge case that diverges from C.

R-1. **`--exp` boundary-init FP-order bug — CLOSED 2026-06-02.**
   Surfaced while wiring `--exp` (gap #6). The `initverticalw`
   boundary init in `profile_align_imp_multimtx` used a nested
   `mul_add` (`fma(C,D, fma(A,B, init))`) instead of clang's
   `init + (A*B + C*D)` order — a 1-ULP FP non-associativity
   mismatch invisible at `--exp 0` (the `+= 0 * i` second term
   absorbs the drift) but amplified at any non-zero `--exp` by the
   per-row `+= fpenalty_ex * i`. Fix: rewrote to match clang's
   FMA codegen (same pattern `currentw` init was already using —
   a 3-line change to `profile.rs:1218-1230`). Result: 44/45
   `--exp` × mode combinations now byte-identical to C
   (`--exp ∈ {0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0}`
   × {FFT-NS-2, FFT-NS-i, L-INS-i, G-INS-i}). The one residual at
   `--exp ≥ 4.5` on FFT-NS-2 (no refinement) is a tied-trace gap
   placement (same width, same score, 16 lines of single-char gap
   shifts) — a §B.2-class tie-break at pathologically high
   penalty_ex values. Out-of-scope for realistic workflows.

R-2. **`--bestfirst` refinement strategy not implemented.** Surfaced
   while wiring iteration-strategy flags. C's `BESTFIRST`
   (`parallelizationstrategy=BESTFIRST` →
   `dvtditr -p BESTFIRST`) restructures the iterative-refinement
   loop: instead of BAATARI2's hill-climbing (refine each branch
   sequentially, accept improvements as found), BESTFIRST scores
   every branch first, sorts by score, then refines them in
   score-order until no improvement. Porting requires (a)
   surfacing per-branch baseline scores BEFORE attempting any
   merge, (b) sorting branches by descending score per iteration,
   (c) a second pass to actually refine in that order. Existing
   refinement loop in `refinement.rs::iterative_refine` walks
   branches in topology order (matching BAATARI2). **To
   investigate:** add a `RefinementStrategy` enum to
   `RefinementParams`; implement the BESTFIRST path as a separate
   loop that pre-scores then re-orders; verify byte-identity
   against C with `--bestfirst` on the 36-seq sample.

R-3. **`--skipiterate F` branch-skip gate not implemented.**
   Surfaced while wiring iteration-strategy flags. C's `dvtditr
   -E $fixthreshold` skips refining branches where
   `distFromABranch > F`. Conceptually a single conditional
   inside the per-branch refinement loop. **To investigate:** add
   `fix_threshold: Option<f64>` to `RefinementParams`; in
   `iterative_refine`, compute `dist_from_a_branch` (already
   exists for `--allowshift`) and skip if `> fix_threshold`;
   verify byte-identity against C with `--skipiterate 0.5`.

R-4. **`--oneiteration` interaction with `--retree` not yet
   mapped.** Surfaced while wiring iteration-strategy flags. C's
   `disttbfast -r` forces a single distance-recomputation
   iteration. We already gate distance-recompute on `--retree`,
   so the relationship between `--oneiteration` and `--retree` in
   our engine needs investigation before wiring. **To investigate:**
   build a small fixture where C `--oneiteration` and C
   `--retree 1` produce different outputs; if they're identical
   for the inputs we care about, document `--oneiteration` as a
   `--retree 1` alias; otherwise port the `-r` semantics
   specifically.

R-5. **DNA strand auto-detection (`--adjustdirection` /
   `--adjustdirectionaccurately`) not yet implemented.** Surfaced
   while wiring gap #1. C MAFFT's strand-detection lives in
   `mafft-upstream/core/makedirectionlist.c` (~1300 LOC) +
   `setdirection.c` (~213 LOC). Algorithm sketch:
   - For each input DNA sequence, build 6-mer count tables for
     both forward orientation and reverse-complement.
   - Compare each sequence against a reference (the longest or
     most-conserved among the first `reflim` sequences) using
     k-mer-derived similarity in both orientations.
   - If the reverse-complement scores higher, mark the sequence
     to be flipped and prepend `_R_` to its name in the output.
   - `--adjustdirectionaccurately` (mode=2) uses `-r 100` (only
     100 reference sequences) and `-d` (more accurate DP-based
     scoring); `--adjustdirection` (mode=1) uses `-r 5000` and
     the faster k-mer-only path.
   - `setdirection` then walks the `_direction` file and applies
     the orientation flip in place (reverse + complement the
     residues).
   **To investigate / port:**
   1. Implement `reverse_complement` and `kmer_count` helpers in
      a new `mafft-tree/src/direction.rs` (k-mer infrastructure
      already exists for the existing ktuple distance code).
   2. Pick a reference: longest sequence is a fine first cut;
      C's `contrastsort` reference selection is a refinement
      we can defer.
   3. For each non-reference sequence, compute
      `cosine_similarity(kmer_fwd, kmer_ref)` and
      `cosine_similarity(kmer_rc, kmer_ref)`. If RC > FWD by a
      threshold (`-t 0.00` in C = always flip when RC better),
      mark for reversal.
   4. In `main.rs` (after read, before align), if any sequence
      is flagged, replace it with its reverse-complement and
      prepend `_R_` to the name.
   5. Validate byte-identity against C on a synthetic input
      with mixed orientations. Note that C's reference-selection
      and per-pair scoring is intricate enough that exact byte
      identity will likely require iterative refinement of the
      port; a "produces same final orientation choice" guarantee
      is the realistic first milestone.
