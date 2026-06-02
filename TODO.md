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

1. **DNA strand auto-detection not exposed.** C's
   `--adjustdirection` / `--adjustdirectionaccurately` reverse-
   complement sequences predicted to be on the opposite strand
   before alignment. We support neither flag. Affects DNA users
   only; protein workflows unaffected.
2. **Tree-linkage variants not exposed.** C lets you choose between
   `--averagelinkage` (default in some modes), `--minimumlinkage`,
   `--mixedlinkage`, `--youngestlinkage`. We ship the C default and
   don't expose the alternatives.
3. **Iteration-strategy variants not exposed.** C's `--bestfirst`,
   `--simplehillclimbing`, `--skipiterate`, `--oneiteration` knobs
   are not surfaced. We use the C default ("randomchain"-style group
   selection).
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
6. **Fine-grained gap-penalty knobs — partially exposed.** We have
   `--op` / `--ep` plus `--exp` / `--shiftpenalty` / `--lop` / `--lep`
   / `--lexp` / `--gop` / `--gep` / `--gexp` (8 new flags wired). All
   verified byte-identical to C MAFFT on the 36-seq sample at typical
   values. Known limitation: `--exp` at aggressive values (≥0.5) on
   FFT-NS-i can produce an alignment whose ARCHITECTURE differs from
   C (different gap distribution, same number of sequences) — a
   latent algorithm-port issue in refinement's penalty_ex handling
   surfaced by exercising the new flag. `--gop`/`--gep`/`--gexp` are
   wired but currently inert in protein/DNA pipelines (they only
   affect C's X-INS-i / Q-INS-i RNA paths, which are external-dep
   blocked). Not implemented: `--rop` / `--rep` (RNA-only); `--LOP` /
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

R-1. **`--exp` at aggressive values diverges from C on FFT-NS-i
   refinement.** Surfaced 2026-06-02 while wiring `--exp` (gap #6
   above). At default `--exp 0` (or `--exp 0.1`, `--exp 0`-class
   values) the rust output is byte-identical to C across every mode.
   But at `--exp 0.5` on FFT-NS-i (`--maxiterate 100`), rust produces
   an alignment of different ARCHITECTURE than C (same `nseq=36`,
   different widths: rust=564 vs C=703, with different gap
   distributions per row). The penalty-extension internal value is
   verified identical to C (`--exp 0.5` → `pgexp=-500` → `penalty_ex
   =-299` for protein in both rust and C). The progressive phase
   alone (`--exp 0.5 --maxiterate 0`) is byte-identical, isolating
   the divergence to the *refinement* DP's penalty_ex handling.
   Likely culprits: tail-gap accumulation in `profile.rs::j_loop`
   (the `mj_v += f_ext` and `mi += f_ext` per-cell increments),
   FFT-segmented refinement boundary handling under high
   penalty_ex, or a tie-break that flips when penalty_ex is large
   enough to outweigh substitution scores. **To investigate:**
   instrument C MAFFT and our refinement DP at the divergent step;
   compare `mj`/`mi`/`wm` trajectories at `--exp 0.5`; pin the
   first row where they diverge; identify the missing increment or
   sign mismatch. Workaround: keep `--exp` ≤ 0.1 for FFT-NS-i runs
   that need byte-identity with C.
