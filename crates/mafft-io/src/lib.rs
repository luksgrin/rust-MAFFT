//! I/O routines for MAFFT: FASTA, Clustal, PHYLIP, hat2, and localhom formats.
//!
//! FASTA reading/writing delegates to `noodles-fasta`. The other formats are
//! MAFFT-specific or simple enough that we implement them directly.

mod error;
mod fasta;
mod clustal;
mod phylip;
mod hat2;
mod localhom;
mod detect;

pub use error::IoError;
pub use fasta::{read_fasta, read_fasta_from_reader, write_fasta, write_fasta_to_writer};
pub use clustal::write_clustal;
pub use phylip::write_phylip;
pub use hat2::{read_hat2, write_hat2, Hat2Matrix};
pub use localhom::{read_localhom_table, write_localhom_table};
pub use detect::detect_seq_type;
