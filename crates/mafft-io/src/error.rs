use std::io;

#[derive(Debug, thiserror::Error)]
pub enum IoError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    #[error("FASTA parse error: {0}")]
    FastaParse(String),

    #[error("hat2 format error: {0}")]
    Hat2Format(String),

    #[error("localhom format error: {0}")]
    LocalHomFormat(String),

    #[error("no sequences found in input")]
    EmptyInput,

    #[error("sequence count mismatch: expected {expected}, got {got}")]
    SeqCountMismatch { expected: usize, got: usize },
}
