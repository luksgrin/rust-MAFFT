//! Scoring matrices and initialization for MAFFT alignment.
//!
//! Ports the C `constants.c`, `blosum.c`, `JTT.c`, `DNA.h`, and `miyata.h`
//! into pure Rust. Provides owned `ScoringContext` structs that replace the
//! C global variables.

mod alphabet;
mod blosum;
mod jtt;
mod dna;
mod properties;
mod penalties;
mod normalize;

pub use alphabet::{ProteinAlphabet, DnaAlphabet, Alphabet, PROTEIN_ALPHABET, DNA_ALPHABET};
pub use blosum::blosum_matrix;
pub use jtt::{jtt_matrix, jtt_frequencies, tm_frequencies};
pub use dna::{default_dna_matrix, DNA_RIBOSUM4, DNA_RIBOSUM16};
pub use properties::{POLARITY, VOLUME, normalized_polarity, normalized_volume};
pub use penalties::{GapParams, default_dna_gap_params, default_protein_gap_params};
pub use normalize::{normalize_matrix, build_scoring_matrix};

use mafft_types::{SeqType, ScoringModel, ScoringContext, GapPenalties};

/// Build a complete `ScoringContext` for the given model and sequence type.
///
/// This is the main entry point, replacing the C `constants()` function.
/// Unlike the C version, this returns an owned struct instead of setting
/// globals.
pub fn build_context(model: ScoringModel, seq_type: SeqType) -> ScoringContext {
    match seq_type {
        SeqType::Dna | SeqType::Rna => build_dna_context(),
        SeqType::Protein => build_protein_context(model),
        _ => build_protein_context(model), // default to protein
    }
}

fn build_dna_context() -> ScoringContext {
    let raw = default_dna_matrix();
    let gap = default_dna_gap_params();
    let n = DNA_ALPHABET.len();

    // Build 26x26 integer matrix (the DNA alphabet has 26 entries)
    let mut matrix = vec![vec![0i32; n]; n];
    for i in 0..10 {
        for j in 0..10 {
            matrix[i][j] = raw[i][j];
        }
    }

    // Copy u -> t mappings (index 4 copies from index 3)
    for i in 0..5 {
        matrix[4][i] = matrix[3][i];
        matrix[i][4] = matrix[i][3];
    }
    // Upper-case copies (indices 5-9 mirror 0-4)
    for i in 5..10 {
        for j in 5..10 {
            matrix[i][j] = matrix[i - 5][j - 5];
        }
    }

    let consweight = matrix
        .iter()
        .map(|row| row.iter().map(|&v| v as f64).collect())
        .collect();

    ScoringContext {
        substitution_matrix: matrix,
        consweight_matrix: consweight,
        polarity: [0.0; 256],
        volume: [0.0; 256],
        amino_map: DNA_ALPHABET.build_amino_map(),
        model: ScoringModel::Dna,
        seq_type: SeqType::Dna,
        gap: GapPenalties {
            open: gap.penalty,
            extend: gap.penalty_ex,
            offset: gap.offset,
        },
        nalphabets: 26,
        nscoredalphabets: 10,
    }
}

fn build_protein_context(model: ScoringModel) -> ScoringContext {
    let gap_params = default_protein_gap_params();

    let (raw_20x20, freq) = match model {
        ScoringModel::Blosum(n) => {
            let mat = blosum_matrix(n);
            let freq = blosum::blosum_frequencies();
            (mat, freq)
        }
        ScoringModel::Jtt => {
            let mat = jtt::build_jtt_pam_matrix(false, 200);
            let freq = jtt_frequencies();
            (mat, freq)
        }
        ScoringModel::Tm => {
            let mat = jtt::build_jtt_pam_matrix(true, 200);
            let freq = tm_frequencies();
            (mat, freq)
        }
        _ => {
            // Default: JTT with PAM 200
            let mat = jtt::build_jtt_pam_matrix(false, 200);
            let freq = jtt_frequencies();
            (mat, freq)
        }
    };

    let needs_rescale = matches!(model, ScoringModel::Blosum(_));
    let normalized = build_scoring_matrix(&raw_20x20, &freq, gap_params.offset, needs_rescale);

    // Expand 20x20 -> 26x26 (adding B, Z, X, '.', '-', J)
    let n = PROTEIN_ALPHABET.len();
    let mut matrix = vec![vec![0i32; n]; n];
    for i in 0..20 {
        for j in 0..20 {
            matrix[i][j] = normalized[i][j];
        }
    }
    // B = average of N, D (indices 2, 3)
    for i in 0..20 {
        matrix[20][i] = round_half_away((normalized[2][i] + normalized[3][i]) as f64 / 2.0);
        matrix[i][20] = matrix[20][i];
    }
    // Z = average of Q, E (indices 5, 6)
    for i in 0..20 {
        matrix[21][i] = round_half_away((normalized[5][i] + normalized[6][i]) as f64 / 2.0);
        matrix[i][21] = matrix[21][i];
    }
    // X = average of all 20
    for i in 0..20 {
        let avg: f64 = (0..20).map(|j| normalized[i][j] as f64).sum::<f64>() / 20.0;
        matrix[22][i] = round_half_away(avg);
        matrix[i][22] = matrix[22][i];
    }

    let consweight = matrix
        .iter()
        .map(|row| row.iter().map(|&v| v as f64).collect())
        .collect();

    let pol = normalized_polarity();
    let vol = normalized_volume();
    let mut polarity_arr = [0.0f64; 256];
    let mut volume_arr = [0.0f64; 256];
    let alpha = &PROTEIN_ALPHABET;
    for (i, &ch) in alpha.chars.iter().enumerate() {
        if i < 20 {
            polarity_arr[ch as usize] = pol[i];
            volume_arr[ch as usize] = vol[i];
        }
    }

    ScoringContext {
        substitution_matrix: matrix,
        consweight_matrix: consweight,
        polarity: polarity_arr,
        volume: volume_arr,
        amino_map: PROTEIN_ALPHABET.build_amino_map(),
        model,
        seq_type: SeqType::Protein,
        gap: GapPenalties {
            open: gap_params.penalty,
            extend: gap_params.penalty_ex,
            offset: gap_params.offset,
        },
        nalphabets: 26,
        nscoredalphabets: 20,
    }
}

/// Round to nearest integer, rounding 0.5 away from zero (C's shishagonyuu).
fn round_half_away(x: f64) -> i32 {
    if x > 0.0 {
        (x + 0.5) as i32
    } else if x < 0.0 {
        (x - 0.5) as i32
    } else {
        0
    }
}
