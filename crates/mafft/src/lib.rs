//! MAFFT multiple sequence alignment — high-level Rust API.
//!
//! This is the ergonomic entry point. Add it to your project with:
//!
//! ```text
//! cargo add mafft
//! ```
//!
//! and write:
//!
//! ```no_run
//! use mafft::{MafftEngine, AlignmentMode, SequenceSet, read_fasta};
//!
//! let input: SequenceSet = read_fasta("input.fasta").unwrap();
//! let engine = MafftEngine::new(AlignmentMode::FftNs2);
//! let msa = engine.align(&input);
//! for (name, seq) in msa.names.iter().zip(msa.sequences.iter()) {
//!     println!(">{name}");
//!     println!("{}", std::str::from_utf8(seq).unwrap());
//! }
//! ```
//!
//! # What this crate re-exports
//!
//! * [`mafft_core`] — alignment engine, modes, MSA result types
//! * [`mafft_types`] — `Sequence`, `SequenceSet`, scoring models, segment types
//! * [`mafft_io`] — FASTA / hat2 / localhom readers and writers
//!
//! Items from these crates are flattened into the root namespace below,
//! so most callers never need to write `mafft::core::` / `mafft::io::`
//! paths — `use mafft::*` is enough.
//!
//! # When to depend on the sub-crates directly
//!
//! Reach for `mafft-core` / `mafft-align` / `mafft-tree` / `mafft-scoring`
//! / `mafft-fft` directly only if you need to:
//!
//! * cut compile time by avoiding the I/O layer,
//! * pin a sub-crate to a specific version independently of the rest, or
//! * extend internals (e.g. custom guide trees, custom scoring matrices).
//!
//! For everything else, depend on `mafft`.
//!
//! # The command line's flag layer, in-process (`cli` feature)
//!
//! `MafftEngine` takes an [`AlignmentMode`]; the command line *chooses* one
//! (`--auto` from sequence count and length) and does more before the engine
//! runs (`--adjustdirection` strand detection, `--nuc` / `--amino` type
//! forcing, the residue case fold). Re-deriving any of that in a caller is
//! how it silently diverges from C MAFFT. With `features = ["cli"]` the
//! [`cli`] module re-exports the `mafft-rs` crate, whose entry points take
//! the same argv the shell would and go through the same code as the
//! binary:
//!
//! ```no_run
//! # #[cfg(feature = "cli")] {
//! use mafft::cli::{run_from_seqs, SilentProgress};
//! use mafft::{Sequence, SequenceSet, SeqType};
//!
//! let input = SequenceSet {
//!     sequences: vec![
//!         Sequence { name: "a".into(), data: b"atggctagcttggacc".to_vec() },
//!         Sequence { name: "b".into(), data: b"atggctagcttgcacc".to_vec() },
//!     ],
//!     seq_type: SeqType::Dna,
//! };
//! // Same flags as `mafft --auto --adjustdirection --thread 1 --nuc FILE`,
//! // but the sequences stay in memory and the rows come back as a value.
//! let msa = run_from_seqs(
//!     ["mafft", "--auto", "--adjustdirection", "--thread", "1", "--nuc"],
//!     &input,
//!     &SilentProgress,
//! )?;
//! assert_eq!(msa.names.len(), 2);
//! # }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! `run_from(argv, &mut out)` is the same thing with a FASTA file in and
//! formatted text out; `Mafft::new().auto().nuc().run_seqs(&input)` is the
//! typed builder over either.
//!
//! # Related crates
//!
//! * [`mafft-rs`](https://crates.io/crates/mafft-rs) — standalone CLI:
//!   `cargo install mafft-rs`; also the crate behind the `cli` feature
//! * [`pymafft`](https://pypi.org/project/pymafft/) — Python bindings
//!   (`pip install pymafft`)

pub use mafft_core::*;
pub use mafft_types::*;
pub use mafft_io::*;

/// The `mafft-rs` crate's argv-driven entry points (`run_from`,
/// `run_from_seqs`, `Mafft`, `MafftError`, `Progress`). Needs the `cli`
/// feature.
#[cfg(feature = "cli")]
pub mod cli {
    pub use mafft_rs::*;
}

/// Sub-crate re-exports under explicit names, for callers who prefer
/// disambiguation over the flattened root namespace.
pub mod core {
    pub use mafft_core::*;
}

/// Sequence / scoring / segment types.
pub mod types {
    pub use mafft_types::*;
}

/// FASTA / hat2 / localhom I/O.
pub mod io {
    pub use mafft_io::*;
}
