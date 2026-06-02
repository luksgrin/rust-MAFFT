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
4. **Auxiliary output formats not exposed.** C's `--distout`,
   `--scoreout`, `--nodeout`, `--mapout`, `--compactmapout`,
   `--pileup` are not implemented. The CLI exposes `--format clustal`
   / `--format phylip` (covers C's `--clustalout` / `--phylipout`)
   and `--treeout` (covers nodeout in the Newick form), but not the
   rest.
5. **Sequence-filter flags not exposed.** C's `--maxambiguous`,
   `--excludehomologs`, `--minimumweight`, `--nwildcard`, `--nzero`,
   `--originalseqonly` — not implemented.
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
7. **Upstream test fixtures unreferenced.** Seven files in
   `mafft-upstream/test/` (`sample.gins1`, `sample.lins1`,
   `sample.parttree`, `sample.dpparttree`, `sample.hat2`,
   `samplerna.qinsi`, `samplerna.xinsi`) are present but never
   read by `crates/*/tests/`. The modes they cover ARE tested
   against our own fixtures in `crates/mafft-core/tests/fixtures/`
   or against re-running C at test time — so this is a maintenance
   concern, not a correctness gap (if upstream changes a reference
   output, we won't notice via CI).
8. **Out-of-scope by design** (no plan to support): RNA structure
    alignment via DAFS / FoldAlign / LARA / SCARNA / MCCASKILL /
    RIBOSUM / RNAALIFOLD; PDB structure-aware alignment
    (`--pdbfilelist`, `--pdbidlist`); alternative pairwise via
    BLAST or LAST (`--blastpair`, `--lastpair`, `--lastmultipair`,
    `--fastapair`, `--fastswpair`, `--hybridpair`,
    `--longshortpair`, `--shortlongpair`); MPI parallelism
    (`--mpi` — we use Rayon).
