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

    /// A description line preceded by blanks (`^[[:blank:]]+>`). C MAFFT
    /// rejects the whole input (`scripts/mafft:1827-1834`, exit 1); a
    /// lenient parser would silently fold the line into the previous
    /// sequence. `line` is 1-based, `text` the offending line verbatim.
    #[error("The first character of a description line must be \nthe greater-than (>) symbol, not a blank.\nPlease check the format around the following line(s):\n{line}:{text}")]
    BlankBeforeHeader { line: usize, text: String },

    /// `=`, `<` or `>` inside a sequence line on the case-preserving
    /// (`--anysymbol` / `--preservecase`) path. C's `charfilter`
    /// (`io.c:1329-1352`) exits 1 on these because `readData_pointer`
    /// uses them as record markers.
    #[error("Characters '= < >' can be used only in the title lines in the --anysymbol or --text mode.")]
    IllegalTitleCharInSequence,

    #[error("sequence count mismatch: expected {expected}, got {got}")]
    SeqCountMismatch { expected: usize, got: usize },
}
