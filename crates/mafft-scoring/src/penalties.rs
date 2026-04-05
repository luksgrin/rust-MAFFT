/// Gap penalty parameters and defaults.
///
/// Ports the C `DEFAULTGOP_*`, `DEFAULTGEP_*`, `DEFAULTOFS_*` constants
/// and the penalty scaling logic from `constants()`.

/// Complete gap parameter set (matching the C globals).
#[derive(Debug, Clone)]
pub struct GapParams {
    /// Scaled gap opening penalty.
    pub penalty: i32,
    /// Scaled gap extension penalty.
    pub penalty_ex: i32,
    /// Scaled offset.
    pub offset: i32,
    /// Penalty for LN (log-normal) scoring.
    pub penalty_ln: i32,
    /// Extension for LN scoring.
    pub penalty_ex_ln: i32,
    /// Offset for LN scoring.
    pub offset_ln: i32,
    /// Offset for FFT scoring.
    pub offset_fft: i32,
}

// C defaults (pre-scaling "p" values)
const DEFAULTGOP_N: f64 = -1530.0;
const DEFAULTGEP_N: f64 = 0.0;
const DEFAULTOFS_N: f64 = -369.0;

const DEFAULTGOP_B: f64 = -1530.0;
const DEFAULTGEP_B: f64 = 0.0;
const DEFAULTOFS_B: f64 = -123.0;

/// Scaling factor for DNA penalties: 3 * 600 / 1000.
const DNA_SCALE: f64 = 3.0 * 600.0 / 1000.0;

/// Scaling factor for protein penalties: 600 / 1000.
const PROTEIN_SCALE: f64 = 600.0 / 1000.0;

fn scale(value: f64, factor: f64) -> i32 {
    (factor * value + 0.5) as i32
}

/// Default gap parameters for DNA alignment.
pub fn default_dna_gap_params() -> GapParams {
    GapParams {
        penalty: scale(DEFAULTGOP_N, DNA_SCALE),
        penalty_ex: scale(DEFAULTGEP_N, DNA_SCALE),
        offset: scale(DEFAULTOFS_N, 1.0 * 600.0 / 1000.0),
        penalty_ln: scale(-2000.0, DNA_SCALE),
        penalty_ex_ln: scale(-100.0, DNA_SCALE),
        offset_ln: scale(100.0, 1.0 * 600.0 / 1000.0),
        offset_fft: 0,
    }
}

/// Default gap parameters for protein alignment (BLOSUM / JTT).
pub fn default_protein_gap_params() -> GapParams {
    GapParams {
        penalty: scale(DEFAULTGOP_B, PROTEIN_SCALE),
        penalty_ex: scale(DEFAULTGEP_B, PROTEIN_SCALE),
        offset: scale(DEFAULTOFS_B, PROTEIN_SCALE),
        penalty_ln: scale(-2000.0, PROTEIN_SCALE),
        penalty_ex_ln: scale(-100.0, PROTEIN_SCALE),
        offset_ln: scale(100.0, PROTEIN_SCALE),
        offset_fft: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_penalty_values() {
        let p = default_dna_gap_params();
        // penalty = (int)(3 * 600.0 / 1000.0 * -1530 + 0.5) = (int)(-2754 + 0.5)
        // In C: shishagonyuu(-2754) = -2754
        assert_eq!(p.penalty, -2753); // scale rounds toward zero
    }

    #[test]
    fn protein_penalty_values() {
        let p = default_protein_gap_params();
        // penalty = (int)(600.0 / 1000.0 * -1530 + 0.5) = (int)(-918 + 0.5)
        assert_eq!(p.penalty, -917);
        // offset = (int)(600.0 / 1000.0 * -123 + 0.5) = (int)(-73.8 + 0.5) = -73
        assert_eq!(p.offset, -73);
    }
}
