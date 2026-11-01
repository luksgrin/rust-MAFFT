/// The type of biological sequences being aligned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqType {
    Protein,
    Dna,
    Rna,
    Text,
    Unknown,
}

impl SeqType {
    /// Convert from the C `dorp` global variable convention.
    /// 'd' = DNA/RNA, 'p' = Protein, NOTSPECIFIED = Unknown.
    pub fn from_dorp(dorp: i32) -> Self {
        match dorp as u8 {
            b'd' => Self::Dna,
            b'p' => Self::Protein,
            _ => Self::Unknown,
        }
    }

    pub fn is_nucleotide(self) -> bool {
        matches!(self, Self::Dna | Self::Rna)
    }
}

/// A FASTA header after the leading `>`, kept as the bytes on the line.
///
/// C MAFFT copies those bytes through. They are not required to be
/// UTF-8, and a NUL does not end the name (C's line reader would stop
/// there; this type keeps the bytes after it). [`Display`] writes valid
/// UTF-8 unchanged and uses the replacement character only when the
/// bytes are not UTF-8, which is for diagnostics. Writers that must
/// match C emit [`SeqName::as_bytes`].
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SeqName(Vec<u8>);

impl SeqName {
    pub fn as_bytes(&self) -> &[u8] {
        &self.0
    }

    pub fn into_bytes(self) -> Vec<u8> {
        self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `prefix` plus these bytes. Same result as formatting the prefix
    /// in front of a UTF-8 name, and it does not recode a name that is
    /// not UTF-8.
    pub fn prefixed(&self, prefix: &[u8]) -> Self {
        let mut v = Vec::with_capacity(prefix.len() + self.0.len());
        v.extend_from_slice(prefix);
        v.extend_from_slice(&self.0);
        Self(v)
    }

    pub fn starts_with(&self, prefix: &str) -> bool {
        self.0.starts_with(prefix.as_bytes())
    }

    pub fn contains(&self, needle: &str) -> bool {
        if needle.is_empty() {
            return true;
        }
        let n = needle.as_bytes();
        self.0.windows(n.len()).any(|w| w == n)
    }

    /// Split on the first `delim` byte. The `|` in a compact-map name.
    pub fn split_once_byte(&self, delim: u8) -> Option<(&[u8], &[u8])> {
        let i = self.0.iter().position(|&c| c == delim)?;
        Some((&self.0[..i], &self.0[i + 1..]))
    }

    /// First whitespace-delimited word, matching `str::split_whitespace`
    /// when the name is UTF-8 (so a non-breaking space still splits) and
    /// ASCII whitespace otherwise.
    pub fn first_word(&self) -> &[u8] {
        match std::str::from_utf8(&self.0) {
            Ok(s) => s.split_whitespace().next().unwrap_or(s).as_bytes(),
            Err(_) => {
                let start = self
                    .0
                    .iter()
                    .position(|c| !c.is_ascii_whitespace())
                    .unwrap_or(self.0.len());
                let rest = &self.0[start..];
                let end = rest
                    .iter()
                    .position(|c| c.is_ascii_whitespace())
                    .unwrap_or(rest.len());
                &rest[..end]
            }
        }
    }

    /// UTF-8 view. `None` when the header is not UTF-8; callers that
    /// must keep the bytes use [`as_bytes`](Self::as_bytes).
    pub fn as_str(&self) -> Option<&str> {
        std::str::from_utf8(&self.0).ok()
    }
}

impl std::fmt::Debug for SeqName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match std::str::from_utf8(&self.0) {
            Ok(s) => std::fmt::Debug::fmt(s, f),
            Err(_) => std::fmt::Debug::fmt(&self.0, f),
        }
    }
}

impl std::fmt::Display for SeqName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match std::str::from_utf8(&self.0) {
            Ok(s) => f.write_str(s),
            Err(_) => f.write_str(&String::from_utf8_lossy(&self.0)),
        }
    }
}

impl From<&str> for SeqName {
    fn from(s: &str) -> Self {
        Self(s.as_bytes().to_vec())
    }
}

impl From<String> for SeqName {
    fn from(s: String) -> Self {
        Self(s.into_bytes())
    }
}

impl From<Vec<u8>> for SeqName {
    fn from(v: Vec<u8>) -> Self {
        Self(v)
    }
}

impl From<&[u8]> for SeqName {
    fn from(v: &[u8]) -> Self {
        Self(v.to_vec())
    }
}

impl AsRef<[u8]> for SeqName {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl std::borrow::Borrow<[u8]> for SeqName {
    fn borrow(&self) -> &[u8] {
        &self.0
    }
}

impl PartialEq<&str> for SeqName {
    fn eq(&self, other: &&str) -> bool {
        self.0 == other.as_bytes()
    }
}

impl PartialEq<str> for SeqName {
    fn eq(&self, other: &str) -> bool {
        self.0 == other.as_bytes()
    }
}

impl PartialEq<String> for SeqName {
    fn eq(&self, other: &String) -> bool {
        self.0 == other.as_bytes()
    }
}

impl PartialEq<SeqName> for str {
    fn eq(&self, other: &SeqName) -> bool {
        other.0 == self.as_bytes()
    }
}

impl PartialEq<SeqName> for &str {
    fn eq(&self, other: &SeqName) -> bool {
        other.0 == self.as_bytes()
    }
}

impl PartialEq<SeqName> for String {
    fn eq(&self, other: &SeqName) -> bool {
        other.0 == self.as_bytes()
    }
}

/// A named biological sequence.
///
/// `name` is the FASTA header bytes ([`SeqName`]), not a Unicode string.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: SeqName,
    pub data: Vec<u8>,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

/// A collection of sequences to be aligned.
#[derive(Debug, Clone)]
pub struct SequenceSet {
    pub sequences: Vec<Sequence>,
    pub seq_type: SeqType,
}

impl SequenceSet {
    pub fn new(seq_type: SeqType) -> Self {
        Self {
            sequences: Vec::new(),
            seq_type,
        }
    }

    pub fn nseq(&self) -> usize {
        self.sequences.len()
    }

    pub fn max_len(&self) -> usize {
        self.sequences.iter().map(|s| s.len()).max().unwrap_or(0)
    }
}

/// An RNA base pair probability.
///
/// Replaces the C `RNApair` struct.
#[derive(Debug, Clone, Copy)]
pub struct RnaBasePair {
    pub up_pos: i32,
    pub up_score: f64,
    pub down_pos: i32,
    pub down_score: f64,
    pub best_pos: i32,
    pub best_score: f64,
}

impl Default for RnaBasePair {
    fn default() -> Self {
        Self {
            up_pos: -1,
            up_score: 0.0,
            down_pos: -1,
            down_score: 0.0,
            best_pos: -1,
            best_score: 0.0,
        }
    }
}
