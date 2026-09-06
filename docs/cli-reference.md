# CLI Reference

Auto-generated from `mafft-rs --help` on every docs build. The
authoritative source for every flag is the CLI itself — run
`mafft-rs --help` locally for the version you have installed.

**Generated for**: `mafft-rs 0.2.0`

## Usage

```text
MAFFT multiple sequence alignment tool (Rust implementation)

Usage: mafft-rs [OPTIONS] [INPUT]

Arguments:
  [INPUT]  Input FASTA file (reads from stdin if omitted)

Options:
  -o, --output <FILE>                Output file (writes to stdout if omitted)
      --localpair                    Use L-INS-i (local, iterative; most accurate for <200 seqs)
      --globalpair                   Use G-INS-i (global, iterative; for globally alignable seqs)
      --genafpair                    Use E-INS-i (generalized affine, iterative; for seqs with large gaps)
      --add <FILE>                   Add new sequences to an existing alignment (provide aligned FASTA as INPUT, new sequences as --add FILE)
      --addfragments <FILE>          Add fragment sequences to an existing alignment (same as --add but for short fragments)
      --keeplength                   Preserve existing alignment column structure when adding sequences
      --allowshift                   Enable per-step dynamic matrix scaling (sets unalignlevel=0.8). Allows divergent regions to stay unaligned at shallow merges. Requires `--globalpair`
      --unalignlevel <F>             Per-step substitution-score offset = (distfromtip - unalignlevel) * 600 when distfromtip < unalignlevel, else 0. Default 0 = no scaling. `--allowshift` sets this to 0.8 if not explicitly given
      --maxiterate <MAXITERATE>      Maximum number of iterative refinement cycles. Unset → mode default (1000 for INS-i modes, 0 for FFT-NS-2). Explicit 0 disables refinement
      --retree <RETREE>              Number of guide tree rebuilds [default: 2] [default: 2]
      --nofft                        Disable FFT: force pure DP for all alignment steps (NW-NS-2 mode)
      --parttree                     Use PartTree guide tree for large datasets (10K+ sequences)
      --dpparttree                   Use DP-based PartTree (more accurate than --parttree, slower)
      --groupsize <GROUPSIZE>        Group size for PartTree partitioning [default: 150]
      --qinsi                        Use Q-INS-i: RNA secondary structure from McCaskill base-pair probabilities
      --xinsi                        Use X-INS-i: RNA secondary structure from CONTRAfold predictions
      --scarnalike                   Use SCARNA-like structural alignment via DASH
      --pdbidlist <FILE>             PDB ID list file for structure-aware alignment. **Non-functional**: matches C MAFFT 7.526's behaviour — the upstream `scripts/mafft:969-983` disables this flag with "temporarily unavailable, 2018/Dec." and `exit`s before any structural alignment runs. Rust accepts the flag and exits with the same message rather than reporting an unknown arg
      --pdbfilelist <FILE>           PDB file list for structure-aware alignment. **Non-functional**: same status as `--pdbidlist` — C MAFFT 7.526 (`scripts/mafft:980-990`) disables this with "temporarily unavailable, 2018/Dec." and `exit`s. Rust matches that behaviour
      --format <FORMAT>              Output format: fasta (default), clustal, phylip [default: fasta]
      --linewidth <LINEWIDTH>        FASTA line width (0 for unlimited) [default: 60] [default: 60]
      --namelength <N>               Name field width in CLUSTAL/PHYLIP output. Default: 15 for CLUSTAL, 10 for PHYLIP (matches C MAFFT's `clustalout_pointer` / `phylipout_pointer`). Names longer than the field are truncated
      --op <OP>                      Gap opening penalty (positive float, e.g. 1.53) [default: 1.53]
      --ep <EP>                      Offset (gap extension-like penalty, positive float, e.g. 0.123) [default: 0.123]
      --exp <EXP>                    Gap extension penalty (`--exp`). Positive float (e.g. 0.1); negated internally to match C's `gexp = -1.0 * arg` convention. Default 0 (no per-residue extension cost)
      --lop <LOP>                    L-INS-i pairwise gap-open (`--lop`). Signed float, no negation (matches C: `lgop=-2.00` default). `allow_hyphen_values` so negative numbers like `-3.0` are parsed as the value, not a flag
      --lep <LEP>                    L-INS-i pairwise offset (`--lep`). Signed float, no negation (matches C: `laof=0.100` default)
      --lexp <LEXP>                  L-INS-i pairwise gap-extend (`--lexp`). Signed float, no negation (matches C: `lexp=-0.100` default)
      --gop <N>                      X-INS-i / Q-INS-i generalized-affine pair gap-open (`--gop`). Inert in protein/DNA pipelines (only used by RNA structure modes — see `TODO.md` external-dep limitations). Default `pggop=-1.53` in C
      --gep <N>                      X-INS-i / Q-INS-i generalized-affine pair offset (`--gep`). Inert in protein/DNA pipelines. Default `pgaof=0.10`
      --gexp <N>                     X-INS-i / Q-INS-i generalized-affine pair gap-extend (`--gexp`). Inert in protein/DNA pipelines. Default `pgexp=-0.10`
      --rop <N>                      RNA-only: ribosum gap-open penalty (`--rop`, C variable `rgop`, default `-1.530`). Forwarded to `rnaopt` in C MAFFT only for the RNA-structure paths (mccaskill / contrafold / dafs / rnaalifold), all of which depend on external binaries we don't ship. Accepted at the CLI for compatibility but emits a "no-op without --xinsi/--qinsi" warning if used
      --rep <N>                      RNA-only: ribosum gap-extend penalty (`--rep`, C variable `rgep`, default `-0.000`). Same RNA-structure-only gating as `--rop`
      --LOP <N>                      LARA-only: gap-open penalty for the LARA RNA path (`--LOP`, C variable `LGOP`, default `-6.00`). LARA mode requires the `lara` binary which is not shipped here. Accepted at the CLI with a warning if used outside a LARA-mode flag
      --LEXP <N>                     LARA-only: gap-extend penalty for the LARA RNA path (`--LEXP`, C variable `LEXP`). Same gating as `--LOP`
      --GOP <N>                      LARA-only: alternative-format gap-open penalty (`--GOP`, C variable `GGOP`, default `-6.00`). Same gating as `--LOP`
      --GEXP <N>                     LARA-only: alternative-format gap-extend penalty (`--GEXP`, C variable `GEXP`). Same gating as `--LOP`
      --shiftpenalty <SHIFTPENALTY>  Shift penalty factor for `--allowshift` (`--shiftpenalty`). Multiplied by the gap-open penalty to get the per-cell shift cost. Default 2.0 (matches C `spfactor=2.0` when `--allowshift` is on)
      --bl <BL>                      BLOSUM matrix number (30, 45, 50, 62, 80). Only used with --localpair/--globalpair when scoring with BLOSUM
      --jtt <JTT>                    JTT PAM number for substitution scoring (e.g. 100, 200). Mirrors `mafft --jtt N` (default PAM 200)
      --tm <TM>                      Transmembrane (TM) PAM number for substitution scoring (e.g. 100, 200). Mirrors `mafft --tm N` (default PAM 200)
      --kimura <KIMURA>              Kimura R parameter for DNA distance model [default: 2]
      --thread <THREAD>              Number of threads (0 = use all available cores) [default: 0] [default: 0]
  -q, --quiet                        Quiet mode: suppress progress messages
      --cite                         Print the citation for MAFFT (and rust-MAFFT, once published) and exit. rust-MAFFT is a port of MAFFT by Kazutaka Katoh et al.; the scientific contribution is theirs and must be cited in any published work
      --reorder                      Output sequences in guide-tree DFS order (matching C MAFFT `--reorder`)
      --inputorder                   Output sequences in input order (default; matches C MAFFT `--inputorder`)
      --treeout                      Write the guide tree to `<INPUT>.tree` in Newick format (matches C MAFFT `--treeout`). Ignored when input is read from stdin
      --distout                      Write the pairwise distance matrix used by the guide-tree construction to `<INPUT>.hat2` (matches C MAFFT `--distout`). Ignored when reading from stdin (no path to derive the output name from). The matrix is whatever the engine actually used for tree construction — k-mer for FFT-NS-2/FFT-NS-i, pairwise alignment-score-derived for L-INS-i / G-INS-i / E-INS-i
      --scoreout                     Print the unweighted sum-of-pairs score of the final alignment to stderr (matches C MAFFT `--scoreout`'s `Unweighted sum-of-pairs score = N.NNNNN` line)
      --nodeout                      Write the guide tree to `<INPUT>.tree` followed by a per-leaf `Density:` section and a `Node info:` section (matches C MAFFT `--nodeout` → `treeout==2` path in `mltaln9.c:6492-6518`). Implies `--treeout`
      --pileup                       Pileup alignment strategy: build a comb-tree guide (sequence 0 joins 1, then that pair joins 2, etc.) instead of UPGMA, and run a single progressive pass with no refinement. Mirrors C MAFFT's `--pileup` ("Pileup-NS-1" strategy, `scripts/mafft:2169` + `mltaln9.c::createchain`)
      --mapout                       Write the per-input-column mapping (gap-insertion positions) to `<ADDFILE>.map`. Requires `--add` / `--addfragments` plus `--keeplength` (matches C MAFFT `--mapout` → internal `-Z -Y`). Mirrors `reconstructdeletemap` (addfunctions.c:1985)
      --compactmapout                Compact form of `--mapout`. Same requirements. Mirrors C's `reconstructdeletemap_compact` (`addfunctions.c:2047`)
      --maxambiguous <F>             Filter input sequences whose ambiguous-residue fraction exceeds N (range 0.0–1.0). Mirrors C MAFFT's `--maxambiguous` (`filter.c`): for protein, ambiguous = anything outside `ARNDCQEGHILKMFPSTWYV` (case-insensitive). For DNA/RNA, ambiguous = anything outside `ATGCU`. Sequences exceeding the threshold are removed before alignment; runs of N/X are collapsed to a single character (`shortenN`). Default 1.0 (no filtering)
      --minimumweight <F>            Floor for per-sequence weights used in refinement and progressive merge. Mirrors C's `tbfast -W $minimumweight` (`scripts/mafft:1029`). Sequences with weight below this floor are clamped up to it. Default 0.00001
      --nwildcard                    Treat N (DNA/RNA ambiguous) as a wildcard that matches anything positively (mirrors C MAFFT's `--nwildcard`, internal flag `-:`): fills the DNA matrix's `n` row with 25%-self-score values (C `constants.c::nscore`). Affects only DNA workflows; also enabled implicitly when `--unalignlevel` > 0, as in C
      --nzero                        Treat N (DNA/RNA ambiguous) as scoring 0 against everything (mirrors C MAFFT's `--nzero`). This is the default in both C MAFFT and mafft-rs, so the flag is a no-op unless it overrides a prior `--nwildcard`
      --excludehomologs              Exclude near-identical sequences during alignment (mirrors `--excludehomologs`). Documented in C as "works with --dash only"; we don't support `--dash`, so this flag is a no-op for now
      --originalseqonly              Output only the original input sequences (mirrors C's `--originalseqonly`). Documented as "works with --dash only"; no-op without `--dash`
      --averagelinkage               Use pure average linkage for UPGMA cluster joining (sueff = 1.0, matches C `--averagelinkage` / `tbfast -X 1.0`). Mutually exclusive with `--minimumlinkage` and `--mixedlinkage`
      --minimumlinkage               Use pure single-linkage (minimum) for UPGMA cluster joining (sueff = 0.0, matches C `--minimumlinkage` / `tbfast -X 0.0`)
      --mixedlinkage <F>             Use a weighted mix of minimum and average linkage for UPGMA cluster joining (`--mixedlinkage F` → sueff = F, matches C `tbfast -X F`). Range 0.0–1.0. The C default is 0.1; this flag is for explicit overrides
      --youngestlinkage              Use the "youngest" linkage scheme (matches C `--youngestlinkage`). C's `youngestlinkage` is a separate algorithm, not just a `sueff` value; ported as `mafft_tree::youngestlinkage_tree` (k-mer and MSA passes of `compacttree_memsaveselectable`), byte-identical to C
      --bestfirst                    Use C MAFFT's `BESTFIRST` parallelisation strategy for iterative refinement — evaluate every branch against a frozen baseline each cycle, then apply the single best move. Default is `BAATARI2` (simple hill climbing). Byte-identical to C `--thread 1 --bestfirst` (iterate cap 254, as C's script sets)
      --simplehillclimbing           Use the simple hill-climbing refinement strategy (matches C `--simplehillclimbing` → `parallelizationstrategy=BAATARI2`). This IS the default in both C MAFFT and mafft-rs, so the flag is a true no-op except when overriding a prior `--bestfirst`
      --skipiterate <F>              Skip refinement of branches whose distance-from-tip exceeds F (matches C `--skipiterate F` → `dvtditr -E $fixthreshold`). Large F skips refinement entirely with C's WARNING; small F skips every branch whose subtree is a strict subset of a sub-alignment cluster (`dvtditr.c:997-1006`). Byte-identical to C
      --oneiteration                 Run C's one-vs-others refinement after the progressive merge (matches C `--oneiteration` → `disttbfast -r`, `disttbfast.c::dooneiteration`). Applies to FFT-NS-2 / FFT-NS-i only; a no-op for L/G/E-INS-i, as in C. Byte-identical to C
      --adjustdirection              Auto-detect input DNA strand orientation and reverse-complement sequences on the wrong strand before alignment (matches C `--adjustdirection`). Implemented: the k-mer-based detection algorithm is a port of `mafft-upstream/core/makedirectionlist.c` (see `mafft_core::adjust_direction`, `TODO.md` R-5, resolved 2026-06-03). Affects DNA workflows only; protein inputs silently bypass it
      --adjustdirectionaccurately    Slower, more accurate variant of `--adjustdirection` (matches C `--adjustdirectionaccurately`, internally `adjustdirection=2`). Also implemented — it selects the DP scorer (`AdjustMode::Dp`) instead of the k-mer one
      --treein <FILE>                Use a user-supplied guide tree (matches C MAFFT `--treein FILE`). FILE must be in MAFFT's internal tree format: nseq-1 lines of `im jm len0 len1` (1-indexed sequence numbers, im < jm). Convert a standard Newick file with `mafft-upstream/core/newick2mafft.rb`
      --auto                         Automatically select alignment strategy based on input size, matching C MAFFT `--auto` (`scripts/mafft:1290-1343`). Picks L-INS-i, FFT-NS-i, FFT-NS-2, FFT-NS-1, --dpparttree, or --parttree depending on the number of sequences and the longest sequence length. Overrides any other algorithm-selection flag
      --memsavetree                  Use the memory-saving guide-tree algorithm (matches C MAFFT `--memsavetree`). Builds a UPGMA-like tree using k-mer distances computed on the fly, avoiding the O(N²) memory cost of a full distance matrix. Recommended for very large inputs (100k+ seqs); also enabled automatically by `--auto` in that bracket
      --anysymbol                    Allow any non-standard characters in input (matches C MAFFT `--anysymbol`). Before alignment, non-standard residues are replaced with `X` (protein) or `n` (DNA) so the alignment DP can score them; the alignment is then post-processed to restore the original characters (and case). Mirrors C's `replaceu` + `restoreu` external steps
      --preservecase                 Alias for `--anysymbol` (matches C MAFFT `--preservecase`). C maps both flags to the same internal `anysymbol=1` variable
      --leavegappyregion             Restore the pre-7.110 gap-cost behaviour (matches C MAFFT `--leavegappyregion` / `--legacygappenalty`). The profile DP stops down-weighting columns by their gap fraction (`legacygapcost = 1` — `Salignmm.c:1604-1610`), which lets alignments leave heavily-gapped regions untouched instead of inserting more gaps to align around them
      --seed <FILE>                  Use a pre-aligned seed alignment as a strong-importance constraint (matches C MAFFT `--seed FILE`). The flag is repeatable — every seed file's sequences are prepended to the user input with a `_seed_` name prefix, and all in-group pairs are written to a `hat3.seed`-style local-homology table with `opt` multiplied by `tsuyosa = user_nseq² * 100` so the refinement DP follows them tightly. Forces `maxiterate ≥ 2` (`scripts/mafft:1911-1923`)
      --seedtable <FILE>             Pre-computed seed local-homology table (matches C MAFFT `--seedtable FILE`, `scripts/mafft:1021-1024`, `scripts/mafft:2437-2438`). The file follows the same 9-field text format `multi2hat3s.c:214` writes for `--seed`: `i j overlapaa opt start1 end1 start2 end2 k` (one record per line, 0-based sequence indices and inclusive residue positions; `opt` is the already-`tsuyosa`-boosted score). Unlike `--seed`, no sequences are prepended — the file's `i`/`j` reference whatever indices the user's input FASTA contains. Mutually exclusive with `--seed`, `--add`/`--addfragments`, `--parttree`/`--dpparttree`, and `--memsave`. Forces `maxiterate ≥ 2` (`scripts/mafft:1911-1923`)
      --memsave                      Memory-saving mode (matches C MAFFT `--memsave` → `tbfast -M -B`, `scripts/mafft:543-544`). In C, this routes the profile DP through `MSalignmm` (Hirschberg-style linear-space divide-and-conquer) for the group merge, avoiding the O(N×M) allocation. For sequences that fit in memory (≤ 30000 residues) the alignment is byte-identical to default-mode output — C's auto-switch at `len > 30000` (`tbfast.c:1096`) makes the two paths converge for typical inputs. Our engine uses full-memory DP regardless; the flag is accepted (for CLI parity) and validated against the same script-level gating as C. NOTE: the actual Hirschberg DP is not yet ported — long sequences that would auto-trigger C's memsave path may OOM here. Tracked in TODO §B.3
      --nomemsave                    Disable the auto-switch to memory-saving DP for long sequences (matches C MAFFT `--nomemsave` → `tbfast -N`, `scripts/mafft:545-546`). In C this sets `nevermemsave = 1` so long-sequence inputs use the full-memory DP. We always use the full DP, so the flag is accepted for CLI parity and has no runtime effect
      --nuc                          Force the input to be treated as nucleotide, overriding the ATGC-frequency auto-detection (matches C MAFFT `--nuc` → `seqtype="-D"`, `scripts/mafft:547-548`). Mutually exclusive with `--amino`
      --amino                        Force the input to be treated as protein, overriding the ATGC-frequency auto-detection (matches C MAFFT `--amino` → `seqtype="-P"`, `scripts/mafft:549-550`). Mutually exclusive with `--nuc`
      --c-compat                     Replicate C MAFFT's static-TLS `reuseprofiles` memoization (`Salignmm.c:1446-1450`) so tied-DP-cell choices match C byte-for-byte. Off by default — the stateless progressive engine is the design goal. Enable when downstream byte-equality with C MAFFT 7.526 is hard-required on inputs that surface the `A__align` static-state artifact (BB20027-class cases, see the deep-dive in `balibase_parity_run.md`; every BAliBASE fixture is byte-identical without this flag today)
  -h, --help                         Print help
  -V, --version                      Print version
```

## Notes on flag semantics

- `--maxiterate 0` disables iterative refinement (FFT-NS-2 default).
- `--maxiterate N` (N > 0) enables refinement up to N cycles per
  strategy. INS-i strategies default to 1000 if `--maxiterate` is
  unset.
- `--quiet` suppresses the per-step progress noise to stderr; the
  alignment itself still writes to stdout.
- `--nofft` forces NW/SW pairwise alignment for every merge step —
  exact behavior, but slower than FFT-anchored DP.
- `--pdbidlist` and `--pdbfilelist` print a "temporarily unavailable,
  2018/Dec." message and exit 0, matching the C MAFFT 7.526 behavior
  verbatim (these flags were disabled upstream).

## Adding new flags

CLI flags live in `crates/mafft-bin/src/lib.rs` as fields on the
`Args` struct (parsed via `clap` derive macros). After adding a flag,
re-run this script to regenerate the reference; the next docs build
will pick it up automatically.
