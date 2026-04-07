use std::io::{self, BufReader, Write};
use std::path::PathBuf;

use clap::Parser;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::{read_fasta, read_fasta_from_reader, write_fasta_to_writer_with_width};
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

    /// Maximum number of iterative refinement cycles [default: 0]
    #[arg(long, default_value_t = 0)]
    maxiterate: usize,

    /// Number of guide tree rebuilds [default: 2]
    #[arg(long, default_value_t = 2)]
    retree: usize,

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

    // Determine alignment mode
    let mode = determine_mode(&args);

    if !args.quiet {
        let mode_name = match &mode {
            AlignmentMode::FftNs2 => "FFT-NS-2",
            AlignmentMode::FftNsi { .. } => "FFT-NS-i",
            AlignmentMode::GInsi { .. } => "G-INS-i",
            AlignmentMode::LInsi { .. } => "L-INS-i",
            AlignmentMode::EInsi { .. } => "E-INS-i",
        };
        let seq_type = if input.seq_type.is_nucleotide() { "nuc" } else { "aa" };
        eprintln!("mafft-rs v{}", env!("CARGO_PKG_VERSION"));
        eprintln!("{nseq} sequences ({seq_type}), strategy: {mode_name}");
    }

    // Align
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
    // --kimura is stored for future use when DNA PAM generation is implemented
    let msa = engine.align(&input);

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
    if args.localpair {
        let iters = if args.maxiterate > 0 { args.maxiterate } else { 1000 };
        AlignmentMode::LInsi { iterations: iters }
    } else if args.globalpair {
        let iters = if args.maxiterate > 0 { args.maxiterate } else { 1000 };
        AlignmentMode::GInsi { iterations: iters }
    } else if args.genafpair {
        let iters = if args.maxiterate > 0 { args.maxiterate } else { 1000 };
        AlignmentMode::EInsi { iterations: iters }
    } else if args.maxiterate > 0 {
        AlignmentMode::FftNsi { iterations: args.maxiterate }
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
