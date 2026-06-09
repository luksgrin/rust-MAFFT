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

/// A named biological sequence.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: String,
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
