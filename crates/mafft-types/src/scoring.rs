use crate::SeqType;

/// Gap penalty parameters.
#[derive(Debug, Clone)]
pub struct GapPenalties {
    /// Gap opening penalty.
    pub open: i32,
    /// Gap extension penalty.
    pub extend: i32,
    /// Offset added to all pairs.
    pub offset: i32,
}

impl Default for GapPenalties {
    fn default() -> Self {
        Self {
            open: -1530,
            extend: -100,
            offset: 0,
        }
    }
}

/// Which substitution matrix model to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringModel {
    Blosum(i32), // 30, 45, 50, 62, 80
    Jtt,
    Tm,
    Dna,
    UserDefined,
}

/// The complete scoring context needed by alignment algorithms.
///
/// Owns the substitution matrix and all derived scoring data.
/// This replaces the C globals `amino_dis`, `n_dis`, `polarity[]`, etc.
#[derive(Debug, Clone)]
pub struct ScoringContext {
    /// Substitution score matrix (indexed by internal residue codes).
    /// For proteins: 20x20 (or 26x26 with ambiguity codes).
    /// For DNA: 4x4 (or larger with ambiguity).
    pub substitution_matrix: Vec<Vec<i32>>,
    /// Consistency-weighted version of the substitution matrix.
    pub consweight_matrix: Vec<Vec<f64>>,
    /// Polarity values for each residue (amino acid property).
    pub polarity: [f64; 256],
    /// Volume values for each residue (amino acid property).
    pub volume: [f64; 256],
    /// Character-to-internal-index mapping.
    pub amino_map: [u8; 256],
    /// Scoring model in use.
    pub model: ScoringModel,
    /// Sequence type.
    pub seq_type: SeqType,
    /// Gap penalties.
    pub gap: GapPenalties,
    /// Number of alphabets (26 for protein, 4-16 for DNA).
    pub nalphabets: usize,
    /// Number of scored alphabets (20 for protein standard AAs).
    pub nscoredalphabets: usize,
    /// FFT-specific scoring matrix (n_dis + offset adjustment).
    /// Used by `alignableReagion()` for segment detection scoring.
    pub fft_matrix: Vec<Vec<i32>>,
}

impl ScoringContext {
    pub fn matrix_size(&self) -> usize {
        self.substitution_matrix.len()
    }
}
