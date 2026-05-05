use std::io::{self, BufReader, Write};
use std::path::PathBuf;

use clap::Parser;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::{read_fasta, read_fasta_from_reader};
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

    /// Enable long-range gap shift penalty (warp)
    #[arg(long)]
    allowshift: bool,

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

    // Read input
    let input = match &args.input {
        Some(path) => {
            read_fasta(path).unwrap_or_else(|e| {
                eprintln!("Error reading {}: {e}", path.display());
                std::process::exit(1);
            })
        }
        None => {
            let stdin = io::stdin();
            let reader = BufReader::new(stdin.lock());
            read_fasta_from_reader(reader).unwrap_or_else(|e| {
                eprintln!("Error reading stdin: {e}");
                std::process::exit(1);
            })
        }
    };

    let nseq = input.nseq();
    if nseq == 0 {
        eprintln!("Error: no sequences found in input");
        std::process::exit(1);
    }

    // Check SCARNA-like mode (requires DASH client — network service)
    if args.scarnalike {
        eprintln!("SCARNA-like mode requires the DASH structural alignment client.");
        eprintln!("Install dash_client and ensure it is in your PATH or set MAFFT_BINARIES.");
        eprintln!("See: https://mafft.cbrc.jp/alignment/software/source.html");
        std::process::exit(1);
    }

    // Determine alignment mode
    let mode = determine_mode(&args);

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
        eprintln!("{nseq} sequences ({seq_type}), strategy: {mode_name}");
    }

    // Build engine
    let mut engine = MafftEngine::new(mode).with_retree(args.retree);
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
    if args.allowshift {
        engine = engine.with_allowshift(true);
    }
    if args.parttree {
        engine = engine.with_parttree(true);
    }
    if args.dpparttree {
        engine = engine.with_dpparttree(true);
    }
    if let Some(gs) = args.groupsize {
        engine = engine.with_groupsize(gs);
    }

    // Handle --add / --addfragments
    let add_file = args.add.as_ref().or(args.addfragments.as_ref());
    let msa = if let Some(add_path) = add_file {
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

    if !args.quiet {
        eprintln!("Alignment: {} columns", msa.width());
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

fn determine_mode(args: &Args) -> AlignmentMode {
    // C's `defaultiterate=1000` for INS-i modes (`scripts/mafft:145,150,155`)
    // applies only when `--maxiterate` is not given on the command line.
    // Explicit `--maxiterate 0` disables refinement; `args.maxiterate` is
    // `Option<usize>`, so `None` means unset and `Some(n)` is the user's
    // choice (including `Some(0)`).
    let iters_for = |default: usize| args.maxiterate.unwrap_or(default);
    if args.qinsi {
        AlignmentMode::QInsi { iterations: iters_for(1000) }
    } else if args.xinsi {
        AlignmentMode::XInsi { iterations: iters_for(1000) }
    } else if args.localpair {
        AlignmentMode::LInsi { iterations: iters_for(1000) }
    } else if args.globalpair {
        AlignmentMode::GInsi { iterations: iters_for(1000) }
    } else if args.genafpair {
        AlignmentMode::EInsi { iterations: iters_for(1000) }
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
