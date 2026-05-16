use std::io::{self, BufReader, Write};
use std::path::PathBuf;

use clap::Parser;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::{read_fasta, read_fasta_from_reader, read_fasta_casepreserve, read_fasta_from_reader_casepreserve};
use mafft_types::{Sequence, SequenceSet, ScoringModel};

/// MAFFT-rs: Multiple sequence alignment (Rust implementation)
#[derive(Parser, Debug)]
#[command(name = "mafft-rs", version, about)]
struct Args {
    /// Input FASTA file (reads from stdin if omitted)
    #[arg(value_name = "INPUT")]
    input: Option<PathBuf>,

    /// Output file (writes to stdout if omitted)
    #[arg(short, long, value_name = "FILE")]
    output: Option<PathBuf>,

    // --- Algorithm selection ---
    /// Use L-INS-i (local, iterative; most accurate for <200 seqs)
    #[arg(long)]
    localpair: bool,

    /// Use G-INS-i (global, iterative; for globally alignable seqs)
    #[arg(long)]
    globalpair: bool,

    /// Use E-INS-i (generalized affine, iterative; for seqs with large gaps)
    #[arg(long)]
    genafpair: bool,

    /// Add new sequences to an existing alignment (provide aligned FASTA as INPUT,
    /// new sequences as --add FILE)
    #[arg(long, value_name = "FILE")]
    add: Option<PathBuf>,

    /// Add fragment sequences to an existing alignment (same as --add but for short fragments)
    #[arg(long, value_name = "FILE")]
    addfragments: Option<PathBuf>,

    /// Preserve existing alignment column structure when adding sequences
    #[arg(long)]
    keeplength: bool,

    /// Enable per-step dynamic matrix scaling (sets unalignlevel=0.8). Allows
    /// divergent regions to stay unaligned at shallow merges. Requires
    /// `--globalpair`.
    #[arg(long)]
    allowshift: bool,

    /// Per-step substitution-score offset = (distfromtip - unalignlevel) * 600
    /// when distfromtip < unalignlevel, else 0. Default 0 = no scaling.
    /// `--allowshift` sets this to 0.8 if not explicitly given.
    #[arg(long, value_name = "F")]
    unalignlevel: Option<f64>,

    /// Maximum number of iterative refinement cycles. Unset → mode default
    /// (1000 for INS-i modes, 0 for FFT-NS-2). Explicit 0 disables refinement.
    #[arg(long)]
    maxiterate: Option<usize>,

    /// Number of guide tree rebuilds [default: 2]
    #[arg(long, default_value_t = 2)]
    retree: usize,

    /// Disable FFT: force pure DP for all alignment steps (NW-NS-2 mode)
    #[arg(long)]
    nofft: bool,

    /// Use PartTree guide tree for large datasets (10K+ sequences)
    #[arg(long)]
    parttree: bool,

    /// Use DP-based PartTree (more accurate than --parttree, slower)
    #[arg(long)]
    dpparttree: bool,

    /// Group size for PartTree partitioning [default: 150]
    #[arg(long)]
    groupsize: Option<usize>,

    /// Use Q-INS-i: RNA secondary structure from McCaskill base-pair probabilities
    #[arg(long)]
    qinsi: bool,

    /// Use X-INS-i: RNA secondary structure from CONTRAfold predictions
    #[arg(long)]
    xinsi: bool,

    /// Use SCARNA-like structural alignment via DASH
    #[arg(long)]
    scarnalike: bool,

    /// Output format: fasta (default), clustal, phylip
    #[arg(long, default_value = "fasta")]
    format: String,

    /// FASTA line width (0 for unlimited) [default: 60]
    #[arg(long, default_value_t = 60)]
    linewidth: usize,

    // --- Scoring parameters ---
    /// Gap opening penalty (positive float, e.g. 1.53) [default: 1.53]
    #[arg(long)]
    op: Option<f64>,

    /// Offset (gap extension-like penalty, positive float, e.g. 0.123) [default: 0.123]
    #[arg(long)]
    ep: Option<f64>,

    /// BLOSUM matrix number (30, 45, 50, 62, 80). Only used with --localpair/--globalpair when scoring with BLOSUM.
    #[arg(long)]
    bl: Option<i32>,

    /// JTT PAM number for substitution scoring (e.g. 100, 200). Mirrors `mafft --jtt N` (default PAM 200).
    #[arg(long, conflicts_with_all = ["bl", "tm"])]
    jtt: Option<i32>,

    /// Transmembrane (TM) PAM number for substitution scoring (e.g. 100, 200). Mirrors `mafft --tm N` (default PAM 200).
    #[arg(long, conflicts_with_all = ["bl", "jtt"])]
    tm: Option<i32>,

    /// Kimura R parameter for DNA distance model [default: 2]
    #[arg(long)]
    kimura: Option<i32>,

    /// Number of threads (0 = use all available cores) [default: 0]
    #[arg(long, default_value_t = 0)]
    thread: usize,

    /// Quiet mode: suppress progress messages
    #[arg(long, short)]
    quiet: bool,

    /// Output sequences in guide-tree DFS order (matching C MAFFT `--reorder`).
    #[arg(long, conflicts_with = "inputorder")]
    reorder: bool,

    /// Output sequences in input order (default; matches C MAFFT `--inputorder`).
    #[arg(long)]
    inputorder: bool,

    /// Write the guide tree to `<INPUT>.tree` in Newick format (matches
    /// C MAFFT `--treeout`). Ignored when input is read from stdin.
    #[arg(long)]
    treeout: bool,

    /// Use a user-supplied guide tree (matches C MAFFT `--treein FILE`).
    /// FILE must be in MAFFT's internal tree format: nseq-1 lines of
    /// `im jm len0 len1` (1-indexed sequence numbers, im < jm). Convert
    /// a standard Newick file with `mafft-upstream/core/newick2mafft.rb`.
    #[arg(long, value_name = "FILE")]
    treein: Option<std::path::PathBuf>,

    /// Automatically select alignment strategy based on input size, matching
    /// C MAFFT `--auto` (`scripts/mafft:1290-1343`). Picks L-INS-i, FFT-NS-i,
    /// FFT-NS-2, FFT-NS-1, --dpparttree, or --parttree depending on the
    /// number of sequences and the longest sequence length. Overrides any
    /// other algorithm-selection flag.
    #[arg(long)]
    auto: bool,

    /// Use the memory-saving guide-tree algorithm (matches C MAFFT
    /// `--memsavetree`). Builds a UPGMA-like tree using k-mer distances
    /// computed on the fly, avoiding the O(N²) memory cost of a full
    /// distance matrix. Recommended for very large inputs (100k+ seqs);
    /// also enabled automatically by `--auto` in that bracket.
    #[arg(long)]
    memsavetree: bool,

    /// Allow any non-standard characters in input (matches C MAFFT
    /// `--anysymbol`). Before alignment, non-standard residues are
    /// replaced with `X` (protein) or `n` (DNA) so the alignment DP
    /// can score them; the alignment is then post-processed to
    /// restore the original characters (and case). Mirrors C's
    /// `replaceu` + `restoreu` external steps.
    #[arg(long)]
    anysymbol: bool,

    /// Alias for `--anysymbol` (matches C MAFFT `--preservecase`).
    /// C maps both flags to the same internal `anysymbol=1` variable.
    #[arg(long)]
    preservecase: bool,

    /// Restore the pre-7.110 gap-cost behaviour (matches C MAFFT
    /// `--leavegappyregion` / `--legacygappenalty`). The profile DP
    /// stops down-weighting columns by their gap fraction
    /// (`legacygapcost = 1` — `Salignmm.c:1604-1610`), which lets
    /// alignments leave heavily-gapped regions untouched instead of
    /// inserting more gaps to align around them.
    #[arg(long, alias = "legacygappenalty")]
    leavegappyregion: bool,

    /// Use a pre-aligned seed alignment as a strong-importance
    /// constraint (matches C MAFFT `--seed FILE`). The flag is
    /// repeatable — every seed file's sequences are prepended to the
    /// user input with a `_seed_` name prefix, and all in-group
    /// pairs are written to a `hat3.seed`-style local-homology table
    /// with `opt` multiplied by `tsuyosa = user_nseq² * 100` so the
    /// refinement DP follows them tightly. Forces
    /// `maxiterate ≥ 2` (`scripts/mafft:1911-1923`).
    #[arg(long = "seed", value_name = "FILE")]
    seed_files: Vec<PathBuf>,

    /// Memory-saving mode (matches C MAFFT `--memsave` →
    /// `tbfast -M -B`, `scripts/mafft:543-544`). In C, this routes the
    /// profile DP through `MSalignmm` (Hirschberg-style linear-space
    /// divide-and-conquer) for the group merge, avoiding the O(N×M)
    /// allocation. For sequences that fit in memory (≤ 30000 residues)
    /// the alignment is byte-identical to default-mode output — C's
    /// auto-switch at `len > 30000` (`tbfast.c:1096`) makes the two
    /// paths converge for typical inputs. Our engine uses full-memory
    /// DP regardless; the flag is accepted (for CLI parity) and
    /// validated against the same script-level gating as C. NOTE: the
    /// actual Hirschberg DP is not yet ported — long sequences that
    /// would auto-trigger C's memsave path may OOM here. Tracked in
    /// TODO §B.3.
    #[arg(long)]
    memsave: bool,

    /// Disable the auto-switch to memory-saving DP for long sequences
    /// (matches C MAFFT `--nomemsave` → `tbfast -N`,
    /// `scripts/mafft:545-546`). In C this sets `nevermemsave = 1` so
    /// long-sequence inputs use the full-memory DP. We always use the
    /// full DP, so the flag is accepted for CLI parity and has no
    /// runtime effect.
    #[arg(long)]
    nomemsave: bool,
}

fn main() {
    let args = Args::parse();

    // Configure thread pool
    if args.thread > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.thread)
            .build_global()
            .ok(); // ignore error if pool already initialized
    }

    // Read input. `--anysymbol`/`--preservecase` need every original
    // character preserved (case + non-standard residues) so the
    // post-alignment restore pass can put them back; the default
    // reader normalizes (`* → -`, drops non-alpha) which would lose
    // exactly the chars we need.
    let anysymbol_read = args.anysymbol || args.preservecase;
    let input = match &args.input {
        Some(path) => {
            let result = if anysymbol_read {
                read_fasta_casepreserve(path)
            } else {
                read_fasta(path)
            };
            result.unwrap_or_else(|e| {
                eprintln!("Error reading {}: {e}", path.display());
                std::process::exit(1);
            })
        }
        None => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            let result = if anysymbol_read {
                read_fasta_from_reader_casepreserve(reader)
            } else {
                read_fasta_from_reader(reader)
            };
            result.unwrap_or_else(|e| {
                eprintln!("Error reading stdin: {e}");
                std::process::exit(1);
            })
        }
    };

    let user_nseq = input.nseq();
    if user_nseq == 0 {
        eprintln!("Error: no sequences found in input");
        std::process::exit(1);
    }

    // `--memsave` gating: C MAFFT rejects `--memsave` for every
    // non-ktuples distance mode (`scripts/mafft:1866-1869`). That
    // covers `--localpair`, `--globalpair`, `--genafpair`,
    // `--lastpair`, `--multipair`, etc. — the seed `Impossible`
    // diagnostic. Mirror that here so users get the same error
    // surface.
    if args.memsave
        && (args.localpair || args.globalpair || args.genafpair
            || args.qinsi || args.xinsi || args.scarnalike)
    {
        eprintln!("Impossible");
        std::process::exit(1);
    }
    if args.memsave && !args.seed_files.is_empty() {
        // C rejects --seed + --memsave: MSalignmm doesn't accept
        // local-homology constraints (`tbfast.c:1117-1118`).
        eprintln!("Impossible");
        std::process::exit(1);
    }

    // `--seed FILE` (repeatable): read each pre-aligned seed file with
    // gaps preserved (the seed-pair LH extraction needs the gap
    // pattern), prepend the gap-stripped seed sequences to the user
    // input with `_seed_` name prefixes, and build the seed local-
    // homology table. Mirrors C MAFFT `scripts/mafft:2400-2436` +
    // `multi2hat3s.c`.
    //
    // The combined sequence list (seeds then user input) is what the
    // engine sees; the original user_nseq is preserved here so we can
    // restore the input subset and report mode names accurately.
    let mut input = input;
    let seed_groups_aligned: Vec<Vec<Vec<u8>>>;
    let mut seed_seq_count: usize = 0;
    if !args.seed_files.is_empty() {
        let mut groups: Vec<Vec<Vec<u8>>> = Vec::with_capacity(args.seed_files.len());
        for path in &args.seed_files {
            let seed_set = read_fasta_casepreserve(path).unwrap_or_else(|e| {
                eprintln!("Error reading {}: {e}", path.display());
                std::process::exit(1);
            });
            groups.push(seed_set.sequences.iter().map(|s| s.data.clone()).collect());
            // Prepend renamed (gap-stripped) seed sequences to the input
            // ahead of the user data — matching C's `multi2hat3s` output
            // followed by `cat infile2 >> infile` (`scripts/mafft:2435`).
            for s in &seed_set.sequences {
                let ungapped: Vec<u8> = s.data.iter().copied()
                    .filter(|&c| c != b'-' && c != b'.').collect();
                let renamed = Sequence {
                    name: format!("_seed_{}", s.name),
                    data: ungapped,
                };
                input.sequences.insert(seed_seq_count, renamed);
                seed_seq_count += 1;
            }
        }
        seed_groups_aligned = groups;
    } else {
        seed_groups_aligned = Vec::new();
    }
    let total_nseq = input.nseq();
    if !args.quiet && seed_seq_count > 0 {
        eprintln!("--seed: {} seed sequences across {} file(s)",
                  seed_seq_count, args.seed_files.len());
    }

    // `--anysymbol` / `--preservecase`: snapshot the originals (case and
    // non-standard chars intact) and substitute X (protein) / n (DNA)
    // for any character outside the alignment alphabet before passing
    // the sequences to the DP. After alignment we restore the original
    // characters via name-keyed lookup. Mirrors C `replaceu` +
    // `restoreu` (`mafft-upstream/core/replaceu.c`, `restoreu.c`).
    let anysymbol = args.anysymbol || args.preservecase;
    let originals: Option<std::collections::HashMap<String, Vec<u8>>> = if anysymbol {
        let is_dna = input.seq_type.is_nucleotide();
        let map: std::collections::HashMap<String, Vec<u8>> = input.sequences.iter()
            .map(|s| (s.name.clone(), s.data.clone())).collect();
        for s in input.sequences.iter_mut() {
            replace_unusual(&mut s.data, is_dna);
        }
        Some(map)
    } else {
        None
    };

    // Check SCARNA-like mode (requires DASH client — network service)
    if args.scarnalike {
        eprintln!("SCARNA-like mode requires the DASH structural alignment client.");
        eprintln!("Install dash_client and ensure it is in your PATH or set MAFFT_BINARIES.");
        eprintln!("See: https://mafft.cbrc.jp/alignment/software/source.html");
        std::process::exit(1);
    }

    // `--auto`: pick mode + retree based on input size, mirroring
    // `scripts/mafft:1290-1343`. Overrides --localpair / --globalpair /
    // --genafpair / --parttree / --dpparttree / --maxiterate / --retree.
    let auto_choice = if args.auto {
        let nlen = input.sequences.iter().map(|s| s.data.len()).max().unwrap_or(0);
        Some(decide_auto(total_nseq, nlen))
    } else {
        None
    };

    // Determine alignment mode
    let mut mode = if let Some(ref a) = auto_choice {
        a.mode.clone()
    } else {
        determine_mode(&args)
    };

    // `--seed`: C MAFFT forces `iterate ≥ 2` when seed alignments are
    // present (`scripts/mafft:1911-1923`) — the seed-pair `hat3.seed`
    // constraints only fire during refinement. Lift `0`/`1` iteration
    // counts to 2, and promote progressive-only FFT-NS-2 to FFT-NS-i
    // with 2 iterations.
    if !args.seed_files.is_empty() {
        mode = match mode {
            AlignmentMode::FftNs2 => AlignmentMode::FftNsi { iterations: 2 },
            AlignmentMode::FftNsi { iterations } => {
                AlignmentMode::FftNsi { iterations: iterations.max(2) }
            }
            AlignmentMode::LInsi { iterations } => {
                AlignmentMode::LInsi { iterations: iterations.max(2) }
            }
            AlignmentMode::GInsi { iterations } => {
                AlignmentMode::GInsi { iterations: iterations.max(2) }
            }
            AlignmentMode::EInsi { iterations } => {
                AlignmentMode::EInsi { iterations: iterations.max(2) }
            }
            AlignmentMode::QInsi { iterations } => {
                AlignmentMode::QInsi { iterations: iterations.max(2) }
            }
            AlignmentMode::XInsi { iterations } => {
                AlignmentMode::XInsi { iterations: iterations.max(2) }
            }
        };
    }

    if !args.quiet {
        let mode_name = match &mode {
            AlignmentMode::FftNs2 => "FFT-NS-2",
            AlignmentMode::FftNsi { .. } => "FFT-NS-i",
            AlignmentMode::GInsi { .. } => "G-INS-i",
            AlignmentMode::LInsi { .. } => "L-INS-i",
            AlignmentMode::EInsi { .. } => "E-INS-i",
            AlignmentMode::QInsi { .. } => "Q-INS-i",
            AlignmentMode::XInsi { .. } => "X-INS-i",
        };
        let seq_type = if input.seq_type.is_nucleotide() { "nuc" } else { "aa" };
        eprintln!("mafft-rs v{}", env!("CARGO_PKG_VERSION"));
        eprintln!("{total_nseq} sequences ({seq_type}), strategy: {mode_name}");
    }

    // Build engine. With `--auto`, the retree count comes from the size
    // heuristic; otherwise the CLI `--retree` value (default 2) wins.
    let retree = auto_choice.as_ref().map(|a| a.retree).unwrap_or(args.retree);
    let mut engine = MafftEngine::new(mode).with_retree(retree);
    if let Some(op) = args.op {
        engine = engine.with_gap_open(op);
    }
    if let Some(ep) = args.ep {
        engine = engine.with_gap_offset(ep);
    }
    if let Some(bl) = args.bl {
        engine = engine.with_scoring_model(ScoringModel::Blosum(bl));
    }
    if let Some(pam) = args.jtt {
        engine = engine.with_scoring_model(ScoringModel::Jtt(pam));
    }
    if let Some(pam) = args.tm {
        engine = engine.with_scoring_model(ScoringModel::Tm(pam));
    }
    if let Some(kr) = args.kimura {
        engine = engine.with_kimura(kr);
    }
    if args.nofft {
        engine = engine.with_nofft(true);
    }
    // `--allowshift` sets unalignlevel=0.8 unless `--unalignlevel` is also
    // given (mirrors `scripts/mafft:1424-1428`).
    let unalign_level = match (args.unalignlevel, args.allowshift) {
        (Some(v), _) => v,
        (None, true) => 0.8,
        (None, false) => 0.0,
    };
    if args.allowshift {
        engine = engine.with_allowshift(true);
    }
    if unalign_level > 0.0 {
        engine = engine.with_unalign_level(unalign_level);
    }
    // `--auto` may override parttree/dpparttree based on the size heuristic.
    let parttree = auto_choice.as_ref().map(|a| a.parttree).unwrap_or(args.parttree);
    let dpparttree = auto_choice.as_ref().map(|a| a.dpparttree).unwrap_or(args.dpparttree);
    if parttree {
        engine = engine.with_parttree(true);
    }
    if dpparttree {
        engine = engine.with_dpparttree(true);
    }
    if let Some(gs) = args.groupsize {
        engine = engine.with_groupsize(gs);
    }
    if args.reorder {
        engine = engine.with_reorder(true);
    }
    if let Some(ref tree_path) = args.treein {
        if !tree_path.exists() {
            eprintln!("Cannot open {}", tree_path.display());
            std::process::exit(1);
        }
        engine.treein_path = Some(tree_path.clone());
    }
    if args.leavegappyregion {
        engine.legacy_gap_cost = true;
    }
    // `--memsave`: route the non-FFT progressive merge through
    // `mafft_align::msalignmm` (Hirschberg DP, linear-space). Verified
    // byte-identical to C MSalignmm via FFI cross-validation including
    // the asymmetric-length case closed 2026-05-16 (midw indexing fix
    // — `midw[j] += wm`, not `j+1`).
    if args.memsave {
        engine.memsave_dp = true;
    }
    // `--memsavetree` overrides distance-based UPGMA tree construction with
    // C MAFFT's compacttree_memsaveselectable algorithm. `--auto` may also
    // request memsavetree in the 100k+ bracket — pass that through too.
    let memsavetree_active = args.memsavetree
        || auto_choice.as_ref().map(|a| a.memsavetree).unwrap_or(false);
    if memsavetree_active {
        engine.memsavetree = true;
    }

    // `--seed`: build the seed local-homology table. The engine will
    // (a) merge it into the L-INS-i/G-INS-i/E-INS-i pairwise table
    //     before `recompute_importance`, or
    // (b) use it directly when the chosen mode has no pairwise step
    //     (FFT-NS-i with `--seed`).
    if seed_seq_count > 0 {
        let scoring_model = if input.seq_type.is_nucleotide() {
            ScoringModel::Dna
        } else {
            engine.scoring_model
        };
        let scoring = mafft_scoring::build_context(scoring_model, input.seq_type);
        let mut seed_groups: Vec<mafft_align::SeedGroup> = Vec::new();
        let mut next_idx = 0usize;
        for group in &seed_groups_aligned {
            let n = group.len();
            seed_groups.push(mafft_align::SeedGroup {
                aligned: group.iter().map(|s| s.as_slice()).collect(),
                global_indices: (next_idx..next_idx + n).collect(),
            });
            next_idx += n;
        }
        let seed_table = mafft_align::build_seed_homology_table(
            &seed_groups,
            total_nseq,
            user_nseq,
            &scoring.consweight_matrix,
            &scoring.amino_map,
        );
        engine.seed_homology = Some(seed_table);
    }

    // Handle --add / --addfragments
    let add_file = args.add.as_ref().or(args.addfragments.as_ref());
    let mut msa = if let Some(add_path) = add_file {
        let new_input = read_fasta(add_path).unwrap_or_else(|e| {
            eprintln!("Error reading {}: {e}", add_path.display());
            std::process::exit(1);
        });
        if !args.quiet {
            eprintln!("Adding {} sequences to existing alignment", new_input.nseq());
        }
        engine.add_to_alignment(&input, &new_input, args.keeplength)
    } else {
        engine.align(&input)
    };

    // `--anysymbol`: restore each aligned row to its original characters
    // (case and non-standard residues intact). Mirrors C `restoreu`
    // (`mafft-upstream/core/restoreu.c::fillorichar`): for each aligned
    // sequence, walk every non-gap position and copy the next character
    // from the gap-stripped original.
    if let Some(orig_map) = originals {
        for i in 0..msa.sequences.len() {
            let Some(orig) = orig_map.get(&msa.names[i]) else { continue };
            let orig_no_gaps: Vec<u8> = orig.iter().copied()
                .filter(|c| *c != b'-' && *c != b'.').collect();
            let mut k = 0;
            for c in msa.sequences[i].iter_mut() {
                if *c != b'-' && *c != b'.' && k < orig_no_gaps.len() {
                    *c = orig_no_gaps[k];
                    k += 1;
                }
            }
        }
    }

    if !args.quiet {
        eprintln!("Alignment: {} columns", msa.width());
    }

    // --treeout: write the guide tree to `<INPUT>.tree` in Newick format,
    // mirroring C MAFFT's `cp $TMPFILE/infile.tree $infilename.tree`
    // (`scripts/mafft:2817-2819`). Requires a file-backed input — when
    // reading from stdin we have no path to derive the output name.
    //
    // PartTree uses a distinct tree format: numeric leaves only, no
    // branch lengths (`splittbfast.c:1275-1301,2532-2553`).
    if args.treeout {
        if let Some(input_path) = &args.input {
            let tree_path = {
                let mut p = input_path.clone();
                p.as_mut_os_string().push(".tree");
                p
            };
            let seq_type = input.seq_type;
            let scoring_model = if seq_type.is_nucleotide() {
                mafft_types::ScoringModel::Dna
            } else {
                match args.bl {
                    Some(n) => mafft_types::ScoringModel::Blosum(n),
                    None => match args.jtt {
                        Some(p) => mafft_types::ScoringModel::Jtt(p),
                        None => match args.tm {
                            Some(p) => mafft_types::ScoringModel::Tm(p),
                            None => mafft_types::ScoringModel::Blosum(62),
                        },
                    },
                }
            };
            let scoring = mafft_scoring::build_context(scoring_model, seq_type);

            let newick_opt: Option<String> = if args.dpparttree {
                // `--dpparttree` uses cycle=1 (one `splittbfast` call) with
                // `-U` (`doalign=1`) so distances are computed via
                // `G__align11_noalign( n_disLN, -1200, -60, ... )` on the
                // RAW sequences (`splittbfast.c:1700`). `n_disLN` is the
                // base substitution matrix shifted by `offset - offsetLN`
                // (`constants.c:1431-1437`): for protein default with
                // `poffset = 0` this is `n_dis - 60` (offsetLN = 60).
                //
                // Selfscore uses the BASE matrix diagonal (no offsetLN shift)
                // per `splittbfast.c:3011-3017`.
                let raw_seqs: Vec<Vec<u8>> = input.sequences.iter()
                    .map(|s| s.data.clone()).collect();
                let base_matrix = &scoring.consweight_matrix;
                let amino_map = &scoring.amino_map;
                let nalpha = base_matrix.len();
                let offset_ln = 60.0f64;
                // n_disLN-equivalent: shift residue×residue cells by
                // `-offsetLN`. Cells involving non-residue indices stay 0.
                let nscored = scoring.nscoredalphabets;
                let mut dist_matrix: Vec<Vec<f64>> = (0..nalpha).map(|i| (0..nalpha).map(|j| {
                    if i < nscored && j < nscored {
                        base_matrix[i][j] - offset_ln
                    } else {
                        0.0
                    }
                }).collect()).collect();
                // Mirror C's `makedynamicmtx` '−' row/col skip
                // (`mltaln9.c:15197-15203`): the gap-index row/col is NOT
                // shifted so it stays zero, matching C's amino_dynamicmtx.
                let gap_idx = amino_map[b'-' as usize] as usize;
                if gap_idx < nalpha {
                    for j in 0..nalpha { dist_matrix[gap_idx][j] = 0.0; }
                    for i in 0..nalpha { dist_matrix[i][gap_idx] = 0.0; }
                }
                let selfscore_diag = |i: usize| -> i64 {
                    let mut s = 0.0f64;
                    for &c in &raw_seqs[i] {
                        let idx = amino_map[c as usize] as usize;
                        if idx < nalpha { s += base_matrix[idx][idx]; }
                    }
                    s as i64
                };
                let gap = mafft_align::GapModel::new(-1200.0, -60.0);
                let dist_matrix_ref = &dist_matrix;
                // `splittbfast.c:560` sets `outgap = 1` by default, so
                // `G__align11_noalign` penalizes BOTH terminal gaps —
                // equivalent to head_gap=true, tail_gap=true.
                let pair_dp = |i: usize, j: usize| -> f64 {
                    if i == j { return selfscore_diag(i) as f64; }
                    let aln = mafft_align::global_align(
                        &raw_seqs[i], &raw_seqs[j], dist_matrix_ref, amino_map, &gap, true, true,
                    );
                    aln.score
                };
                let seqs_equal = |i: usize, j: usize| -> bool {
                    raw_seqs[i] == raw_seqs[j]
                };
                let orilen = |i: usize| -> usize { raw_seqs[i].len() };
                mafft_tree::parttree_split::run_parttree_pipeline_with_scorer(
                    raw_seqs.len(), selfscore_diag, orilen, pair_dp, seqs_equal, 50,
                ).map(|r| mafft_tree::parttree_split::parttree_result_to_newick(&r))
            } else if args.parttree {
                // PartTree (cycle=2): C overwrites `infile.tree` with CALL 2's
                // (`fromaln=1`) tree, so we use the same fromaln scoring
                // on the FIRST-pass aligned MSA.
                let source_msa: &Vec<Vec<u8>> = msa.first_pass_sequences
                    .as_ref().unwrap_or(&msa.sequences);
                Some(mafft_tree::parttree_split::compute_parttree_newick_fromaln(
                    source_msa,
                    &scoring.consweight_matrix,
                    &scoring.amino_map,
                    scoring.gap.open as f64,
                ))
            } else {
                msa.guide_tree.as_ref().map(|t|
                    mafft_tree::topology_to_newick(t, &msa.names))
            };
            if let Some(mut newick) = newick_opt {
                // C `mltaln9.c:2818` appends `#by loadtree\n` to the
                // tree file when `--treein` was used. This is a comment
                // line (not part of the Newick string itself) but we
                // mirror it for byte-identical `--treeout` parity.
                if args.treein.is_some() {
                    newick.push_str("#by loadtree\n");
                }
                match std::fs::write(&tree_path, newick) {
                    Ok(_) if !args.quiet => eprintln!("Wrote guide tree to {}", tree_path.display()),
                    Ok(_) => {}
                    Err(e) => eprintln!("Warning: could not write {}: {e}", tree_path.display()),
                }
            }
        } else {
            eprintln!("Warning: --treeout requires a file input (stdin not supported)");
        }
    }

    // Build output SequenceSet (with gaps)
    let output_seqs = SequenceSet {
        sequences: msa.sequences.iter().zip(msa.names.iter()).map(|(seq, name)| {
            Sequence {
                name: name.clone(),
                data: seq.clone(),
            }
        }).collect(),
        seq_type: input.seq_type,
    };

    // Write output
    let write_result = match &args.output {
        Some(path) => {
            let file = std::fs::File::create(path).unwrap_or_else(|e| {
                eprintln!("Error creating {}: {e}", path.display());
                std::process::exit(1);
            });
            let mut writer = io::BufWriter::new(file);
            write_output(&output_seqs, &mut writer, &args)
        }
        None => {
            let stdout = io::stdout();
            let mut writer = io::BufWriter::new(stdout.lock());
            write_output(&output_seqs, &mut writer, &args)
        }
    };

    if let Err(e) = write_result {
        eprintln!("Error writing output: {e}");
        std::process::exit(1);
    }
}

/// `--anysymbol` preprocessor — substitute every character outside
/// the alignment alphabet with the appropriate "unknown" symbol, then
/// canonicalize case to match `replaceu.c::replace_unusual`:
/// - Protein: usual = "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";
///   unknown = 'X', case = `toupper`.
/// - DNA:     usual = "ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";
///   unknown = 'n', case = `tolower`.
fn replace_unusual(seq: &mut [u8], is_dna: bool) {
    let usual_protein: &[u8] = b"ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";
    let usual_dna: &[u8] = b"ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";
    let (usual, unknown) = if is_dna {
        (usual_dna, b'n')
    } else {
        (usual_protein, b'X')
    };
    for c in seq.iter_mut() {
        if !usual.contains(c) {
            *c = unknown;
        } else if is_dna {
            *c = c.to_ascii_lowercase();
        } else {
            *c = c.to_ascii_uppercase();
        }
    }
}

/// Resolved alignment strategy for `--auto`.
#[derive(Debug, Clone)]
struct AutoChoice {
    mode: AlignmentMode,
    retree: usize,
    parttree: bool,
    dpparttree: bool,
    /// C MAFFT sets `treeext="memsavetree"` (and thus `compacttree=2`) for
    /// the 100k-200k brackets — see `scripts/mafft:1319-1328`. Surface it
    /// so the CLI can flip on `engine.memsavetree`.
    memsavetree: bool,
}

/// Mirror C `scripts/mafft:1290-1343` `--auto` heuristic. Picks mode and
/// retree count from `nseq` (sequence count) and `nlen` (longest input
/// sequence length, ungapped).
///
/// C uses `memsavetree` (large-N tree algorithm) at nseq ≥ 100k; we
/// don't have that yet, so for the 100k-200k bracket we fall through to
/// FFT-NS-2 / FFT-NS-1 (the alignment phase is the same; only the tree
/// construction differs). Output for those sizes may diverge from C.
fn decide_auto(nseq: usize, nlen: usize) -> AutoChoice {
    if nlen < 3000 && nseq < 100 {
        AutoChoice { mode: AlignmentMode::LInsi { iterations: 1000 }, retree: 1, parttree: false, dpparttree: false, memsavetree: false }
    } else if nlen < 1000 && nseq < 200 {
        AutoChoice { mode: AlignmentMode::LInsi { iterations: 2 }, retree: 1, parttree: false, dpparttree: false, memsavetree: false }
    } else if nlen < 10000 && nseq < 500 {
        AutoChoice { mode: AlignmentMode::FftNsi { iterations: 2 }, retree: 2, parttree: false, dpparttree: false, memsavetree: false }
    } else if nseq < 20000 {
        AutoChoice { mode: AlignmentMode::FftNs2, retree: 2, parttree: false, dpparttree: false, memsavetree: false }
    } else if nseq < 100000 {
        // C: cycle=2, memsavetree on. See `scripts/mafft:1315-1321`.
        AutoChoice { mode: AlignmentMode::FftNs2, retree: 2, parttree: false, dpparttree: false, memsavetree: true }
    } else if nseq < 200000 {
        // C: cycle=1, memsavetree on. See `scripts/mafft:1322-1328`.
        AutoChoice { mode: AlignmentMode::FftNs2, retree: 1, parttree: false, dpparttree: false, memsavetree: true }
    } else if nlen < 3000 {
        // PartTree + localalign distance (= --dpparttree).
        AutoChoice { mode: AlignmentMode::FftNs2, retree: 1, parttree: true, dpparttree: true, memsavetree: false }
    } else {
        // PartTree + ktuple distance.
        AutoChoice { mode: AlignmentMode::FftNs2, retree: 1, parttree: true, dpparttree: false, memsavetree: false }
    }
}

fn determine_mode(args: &Args) -> AlignmentMode {
    // C's `defaultiterate` (`scripts/mafft:86`) is 0 when invoked as
    // `mafft` — the bare flags `--localpair`/`--globalpair`/`--genafpair`
    // do NOT change it. Only the script-name aliases `linsi` / `ginsi` /
    // `einsi` (lines 142-156) set `defaultiterate=1000`. So `mafft
    // --localpair` runs progressive-only; `linsi` runs 1000-iter refinement.
    //
    // Explicit `--maxiterate N` always wins. `args.maxiterate` is
    // `Option<usize>`, so `None` means unset and `Some(n)` is the user's
    // choice (including `Some(0)`).
    let iters_for = |default: usize| args.maxiterate.unwrap_or(default);
    if args.qinsi {
        AlignmentMode::QInsi { iterations: iters_for(1000) }
    } else if args.xinsi {
        AlignmentMode::XInsi { iterations: iters_for(1000) }
    } else if args.localpair {
        AlignmentMode::LInsi { iterations: iters_for(0) }
    } else if args.globalpair {
        AlignmentMode::GInsi { iterations: iters_for(0) }
    } else if args.genafpair {
        AlignmentMode::EInsi { iterations: iters_for(0) }
    } else if let Some(n) = args.maxiterate.filter(|&n| n > 0) {
        AlignmentMode::FftNsi { iterations: n }
    } else {
        AlignmentMode::FftNs2
    }
}

fn write_output<W: Write>(
    seqs: &SequenceSet,
    writer: &mut W,
    args: &Args,
) -> Result<(), mafft_io::IoError> {
    match args.format.as_str() {
        "clustal" | "clw" => {
            mafft_io::write_clustal(seqs, writer, None, None, None)
        }
        "phylip" | "phy" => {
            mafft_io::write_phylip(seqs, writer, None, None)
        }
        _ => {
            mafft_io::write_fasta_to_writer_with_width(seqs, writer, args.linewidth)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mode_name(m: &AlignmentMode) -> &'static str {
        match m {
            AlignmentMode::FftNs2 => "FFT-NS-2",
            AlignmentMode::FftNsi { .. } => "FFT-NS-i",
            AlignmentMode::LInsi { .. } => "L-INS-i",
            AlignmentMode::GInsi { .. } => "G-INS-i",
            AlignmentMode::EInsi { .. } => "E-INS-i",
            AlignmentMode::QInsi { .. } => "Q-INS-i",
            AlignmentMode::XInsi { .. } => "X-INS-i",
        }
    }

    #[test]
    fn auto_small_picks_linsi_1000() {
        let a = decide_auto(50, 500);
        assert_eq!(mode_name(&a.mode), "L-INS-i");
        if let AlignmentMode::LInsi { iterations } = a.mode {
            assert_eq!(iterations, 1000);
        }
        assert_eq!(a.retree, 1);
        assert!(!a.parttree && !a.dpparttree);
    }

    #[test]
    fn auto_medium_picks_linsi_2() {
        // nlen<1000, nseq<200 but not <100 → L-INS-i iterate=2.
        let a = decide_auto(150, 800);
        assert_eq!(mode_name(&a.mode), "L-INS-i");
        if let AlignmentMode::LInsi { iterations } = a.mode {
            assert_eq!(iterations, 2);
        }
    }

    #[test]
    fn auto_large_picks_fft_nsi() {
        // nlen<10000, nseq<500 but not the LInsi brackets → FFT-NS-i iter=2.
        let a = decide_auto(300, 5000);
        assert_eq!(mode_name(&a.mode), "FFT-NS-i");
        if let AlignmentMode::FftNsi { iterations } = a.mode {
            assert_eq!(iterations, 2);
        }
        assert_eq!(a.retree, 2);
    }

    #[test]
    fn auto_xlarge_picks_fft_ns2() {
        // nseq<20000 → FFT-NS-2 retree=2.
        let a = decide_auto(5000, 5000);
        assert_eq!(mode_name(&a.mode), "FFT-NS-2");
        assert_eq!(a.retree, 2);
        assert!(!a.parttree);
    }

    #[test]
    fn auto_huge_picks_fft_ns1() {
        // nseq>=100000 (but <200000) → FFT-NS-2 retree=1.
        let a = decide_auto(150_000, 500);
        assert_eq!(mode_name(&a.mode), "FFT-NS-2");
        assert_eq!(a.retree, 1);
    }

    #[test]
    fn auto_giant_short_picks_dpparttree() {
        // nseq>=200000, nlen<3000 → --parttree --dpparttree retree=1.
        let a = decide_auto(250_000, 500);
        assert!(a.parttree);
        assert!(a.dpparttree);
        assert_eq!(a.retree, 1);
    }

    #[test]
    fn auto_giant_long_picks_parttree() {
        // nseq>=200000, nlen>=3000 → --parttree (ktuple) retree=1.
        let a = decide_auto(250_000, 5000);
        assert!(a.parttree);
        assert!(!a.dpparttree);
        assert_eq!(a.retree, 1);
    }

    #[test]
    fn replace_unusual_protein_canonicalizes() {
        // Protein usual set: 20 AA + lowercase + '-' '.'.
        // Lowercase gets uppercased; '*', '@', 'U' (selenocys), 'x',
        // 'B', 'J', 'Z' (extended set absent from `usual_protein`) →
        // 'X'.
        let mut seq = b"MKAUlsgVPxxBJL@&*fkdgna-.".to_vec();
        replace_unusual(&mut seq, false);
        assert_eq!(&seq, b"MKAXLSGVPXXXXLXXXFKDGNA-.");
    }

    #[test]
    fn replace_unusual_dna_canonicalizes() {
        // DNA usual set: ATGCU + lowercase + IUPAC ambig + 'X' + '-'.
        // Uppercase ATGC gets lowercased; '@' / '*' / digits-equivalents
        // → 'n'. IUPAC `R`, `Y`, `N` stay (lowercased).
        let mut seq = b"ATGCUNnxx@*RYBDHKMSVW-".to_vec();
        replace_unusual(&mut seq, true);
        // Every alphabetic char in `usual_dna` lowercased; unknowns → 'n'.
        assert_eq!(&seq, b"atgcunnxxnnrybdhkmsvw-");
    }

    #[test]
    fn replace_unusual_preserves_length() {
        // The function operates in place and must not change length.
        let original = b"AaBbCc@*xx.-".to_vec();
        let mut seq = original.clone();
        replace_unusual(&mut seq, false);
        assert_eq!(seq.len(), original.len());
        let mut seq2 = original.clone();
        replace_unusual(&mut seq2, true);
        assert_eq!(seq2.len(), original.len());
    }
}
