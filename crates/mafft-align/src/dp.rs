/// Shared types for dynamic programming alignment.

/// An individual alignment operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignOp {
    /// Match or mismatch (consume one residue from each sequence).
    Match,
    /// Gap in sequence 1 (insertion in seq2).
    Insert,
    /// Gap in sequence 2 (deletion from seq1).
    Delete,
}

/// A pairwise alignment result.
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Aligned sequence 1 (with gap characters '-' inserted).
    pub seq1: Vec<u8>,
    /// Aligned sequence 2 (with gap characters '-' inserted).
    pub seq2: Vec<u8>,
    /// The alignment score.
    pub score: f64,
    /// The sequence of alignment operations.
    pub operations: Vec<AlignOp>,
}

impl Alignment {
    /// Length of the alignment (including gaps).
    pub fn len(&self) -> usize {
        self.seq1.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq1.is_empty()
    }

    /// Count the number of identical positions.
    pub fn identity_count(&self) -> usize {
        self.seq1
            .iter()
            .zip(self.seq2.iter())
            .filter(|(a, b)| a == b && **a != b'-')
            .count()
    }

    /// Fractional identity (identical positions / aligned length excluding gaps).
    pub fn identity(&self) -> f64 {
        let aligned = self
            .seq1
            .iter()
            .zip(self.seq2.iter())
            .filter(|(a, b)| **a != b'-' && **b != b'-')
            .count();
        if aligned == 0 {
            0.0
        } else {
            self.identity_count() as f64 / aligned as f64
        }
    }
}

/// Gap penalty model.
#[derive(Debug, Clone)]
pub struct GapModel {
    /// Gap opening penalty (negative value).
    pub open: f64,
    /// Gap extension penalty (negative value).
    pub extend: f64,
}

impl GapModel {
    pub fn new(open: f64, extend: f64) -> Self {
        Self { open, extend }
    }
}

impl Default for GapModel {
    fn default() -> Self {
        Self {
            open: -918.0,  // MAFFT default: (int)(600/1000 * -1530 + 0.5)
            extend: 0.0,   // MAFFT default
        }
    }
}

/// Score a substitution using a matrix indexed by internal residue codes.
///
/// `amino_map` converts an ASCII character to an internal index.
/// Returns the score from `matrix[idx1][idx2]`.
pub fn score_pair(
    a: u8,
    b: u8,
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
) -> f64 {
    let i = amino_map[a as usize] as usize;
    let j = amino_map[b as usize] as usize;
    if i < matrix.len() && j < matrix[0].len() {
        matrix[i][j] as f64
    } else {
        0.0
    }
}
