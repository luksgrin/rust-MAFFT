/// DNA scoring matrices from DNA.h.

use crate::round_half_away;

/// Generate a DNA scoring matrix via Kimura PAM model.
///
/// Ports C's `generatenuc1pam()` + PAM exponentiation + normalization.
/// This is what C uses by default (kimuraR=2, pamN=200).
///
/// Returns a 4x4 integer scoring matrix for a,g,c,t.
pub fn generate_dna_pam(kimura_r: i32, pam_n: usize, offset: i32) -> [[i32; 4]; 4] {
    let freq = [0.25f64; 4]; // uniform frequencies (C default for fmodel=0)
    let kr = kimura_r as f64;

    // Kimura rate matrix
    let rate = [
        [0.0, kr,  1.0, 1.0],
        [kr,  0.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, kr ],
        [1.0, 1.0, kr,  0.0],
    ];

    // Compute mutability and delta
    let mut mutability = [0.0f64; 4];
    let mut total = 0.0;
    for i in 0..4 {
        let mut m = 0.0;
        for j in 0..4 { m += rate[i][j] * freq[j]; }
        mutability[i] = m;
        total += m * freq[i];
    }
    let delta = 0.01 / total;

    // Build PAM1
    let mut pam1 = [[0.0f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                pam1[i][j] = delta * rate[i][j] * freq[j];
            } else {
                pam1[i][j] = 1.0 - delta * mutability[i];
            }
        }
    }

    // Exponentiate: pamx = pam1^pam_n
    let mut pamx = [[0.0f64; 4]; 4];
    for i in 0..4 { pamx[i][i] = 1.0; } // identity
    for _ in 0..pam_n {
        let prev = pamx;
        for i in 0..4 {
            for j in 0..4 {
                let mut s = 0.0;
                for k in 0..4 { s += prev[i][k] * pam1[k][j]; }
                pamx[i][j] = s;
            }
        }
    }

    // Divide by background frequency
    for i in 0..4 {
        for j in 0..4 {
            pamx[i][j] /= freq[j];
        }
    }

    // Log transform
    for i in 0..4 {
        for j in 0..4 {
            if pamx[i][j] <= 0.0 { pamx[i][j] = 0.00001; }
            pamx[i][j] = pamx[i][j].log10() * 1000.0;
        }
    }

    // Normalize: subtract weighted average
    let mut average = 0.0;
    for i in 0..4 { for j in 0..4 { average += pamx[i][j] * freq[i] * freq[j]; } }
    for i in 0..4 { for j in 0..4 { pamx[i][j] -= average; } }

    // Scale by 600 / diagonal average (C uses uniform 1/4 weighting)
    average = 0.0;
    for i in 0..4 { average += pamx[i][i] * 0.25; }
    for i in 0..4 { for j in 0..4 { pamx[i][j] *= 600.0 / average; } }

    // Subtract offset
    let ofs = offset as f64;
    for i in 0..4 { for j in 0..4 { pamx[i][j] -= ofs; } }

    // Round
    let mut result = [[0i32; 4]; 4];
    for i in 0..4 { for j in 0..4 { result[i][j] = round_half_away(pamx[i][j]); } }
    result
}

/// Default DNA scoring matrix (26x26, mostly zeros).
/// Uses the static table (only used when kimuraR == 9999 in C).
pub fn default_dna_matrix() -> [[i32; 26]; 26] {
    let mut m = [[0i32; 26]; 26];

    // Core 4x4: a(0), g(1), c(2), t(3)
    // Matches score 1000, transitions 600, transversions 0
    m[0][0] = 1000; m[0][1] = 600;
    m[1][0] = 600;  m[1][1] = 1000;
    m[2][2] = 1000; m[2][3] = 600;
    m[3][2] = 600;  m[3][3] = 1000;

    // Gap penalties in the matrix (index 25 = 'O')
    for i in 0..25 {
        m[i][25] = -500;
        m[25][i] = -500;
    }
    m[25][25] = 500;
    m[24][25] = 0; // gap-O
    m[25][24] = 0;

    m
}

/// RIBOSUM 4x4 matrix for RNA stem scoring.
#[rustfmt::skip]
pub static DNA_RIBOSUM4: [[f64; 4]; 4] = [
    //   a       g       c       t
    [  2.22,  -1.46,  -1.86,  -1.39], // a
    [ -1.46,   1.03,  -2.48,  -1.74], // g
    [ -1.86,  -2.48,   1.16,  -1.05], // c
    [ -1.39,  -1.74,  -1.05,   1.65], // t
];

/// RIBOSUM 16x16 matrix for RNA base-pair interactions.
#[rustfmt::skip]
pub static DNA_RIBOSUM16: [[f64; 16]; 16] = [
    [ -2.49, -8.24, -7.04, -4.32, -6.86, -8.39, -5.03, -5.84, -8.84, -4.68,-14.37,-12.64, -4.01, -6.16,-11.32, -9.05],
    [ -8.24, -0.80, -8.89, -5.13, -8.61, -5.38, -5.77, -6.60,-10.41, -4.57,-14.53,-10.14, -5.43, -5.94, -8.87,-11.07],
    [ -7.04, -8.89, -2.11, -2.04, -9.73,-11.05, -3.81, -4.72, -9.37, -5.86, -9.08,-10.45, -5.33, -6.93, -8.67, -7.83],
    [ -4.32, -5.13, -2.04,  4.49, -5.33, -5.61,  2.70,  0.59, -5.56,  1.67, -6.71, -5.17,  1.61, -0.51, -4.81, -2.98],
    [ -6.86, -8.61, -9.73, -5.33, -1.05, -8.67, -4.88, -6.10, -7.98, -6.00,-12.43, -7.71, -5.85, -7.55, -6.63,-11.54],
    [ -8.39, -5.38,-11.05, -5.61, -8.67, -1.98, -4.13, -5.77,-11.36, -4.66,-12.58,-13.69, -5.75, -4.27,-12.01,-10.79],
    [ -5.03, -5.77, -3.81,  2.70, -4.88, -4.13,  5.62,  1.21, -5.95,  2.11, -3.70, -5.84,  1.60, -0.08, -4.49, -3.90],
    [ -5.84, -6.60, -4.72,  0.59, -6.10, -5.77,  1.21,  3.47, -7.93, -0.27, -7.88, -5.61, -0.57, -2.09, -5.30, -4.45],
    [ -8.84,-10.41, -9.37, -5.56, -7.98,-11.36, -5.95, -7.93, -5.13, -3.57,-10.45, -8.49, -2.42, -5.63, -7.08, -8.39],
    [ -4.68, -4.57, -5.86,  1.67, -6.00, -4.66,  2.11, -0.27, -3.57,  5.36, -5.71, -4.96,  2.75,  1.32, -4.91, -3.67],
    [-14.37,-14.53, -9.08, -6.71,-12.43,-12.58, -3.70, -7.88,-10.45, -5.71, -3.59, -5.77, -6.88, -8.41, -7.40, -5.41],
    [-12.64,-10.14,-10.45, -5.17, -7.71,-13.69, -5.84, -5.61, -8.49, -4.96, -5.77, -2.28, -4.72, -7.36, -3.83, -5.21],
    [ -4.01, -5.43, -5.33,  1.61, -5.85, -5.75,  1.60, -0.57, -2.42,  2.75, -6.88, -4.72,  4.97,  1.14, -2.98, -3.39],
    [ -6.16, -5.94, -6.93, -0.51, -7.55, -4.27, -0.08, -2.09, -5.63,  1.32, -8.41, -7.36,  1.14,  3.36, -4.76, -4.28],
    [-11.32, -8.87, -8.67, -4.81, -6.63,-12.01, -4.49, -5.30, -7.08, -4.91, -7.40, -3.83, -2.98, -4.76, -3.21, -5.97],
    [ -9.05,-11.07, -7.83, -2.98,-11.54,-10.79, -3.90, -4.45, -8.39, -3.67, -5.41, -5.21, -3.39, -4.28, -5.97, -0.02],
];
