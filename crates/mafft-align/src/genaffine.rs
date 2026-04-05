/// Generalized affine gap local alignment (for E-INS-i).
///
/// Ports the C `genL__align11()` from genalign11.c.
///
/// The generalized affine model adds a second gap state that allows
/// "jumping" across long unaligned regions with a separate penalty
/// (`penalty_op`), distinct from the standard affine open penalty.
/// This is essential for sequences with large internal insertions.

use crate::dp::{score_pair, AlignOp, Alignment, GapModel};
use crate::local::LocalAlignment;

/// Extended gap model with a separate "generalized" opening penalty.
#[derive(Debug, Clone)]
pub struct GenAffineGapModel {
    /// Standard affine gap model.
    pub affine: GapModel,
    /// Generalized gap opening penalty (for long-range jumps).
    pub open_generalized: f64,
}

impl Default for GenAffineGapModel {
    fn default() -> Self {
        Self {
            affine: GapModel::default(),
            open_generalized: -918.0,
        }
    }
}

/// Perform generalized affine gap local alignment.
///
/// Adds a second pair of DP states beyond standard affine, allowing
/// long-range jumps with a separate penalty. This makes E-INS-i
/// effective for multi-domain proteins with large internal gaps.
pub fn genaffine_local_align(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap_model: &GenAffineGapModel,
) -> LocalAlignment {
    let n = seq1.len();
    let m = seq2.len();

    if n == 0 || m == 0 {
        return LocalAlignment {
            alignment: Alignment {
                seq1: Vec::new(),
                seq2: Vec::new(),
                score: 0.0,
                operations: Vec::new(),
            },
            offset1: 0,
            offset2: 0,
        };
    }

    let gap_open = gap_model.affine.open;
    let gap_ext = gap_model.affine.extend;
    let gap_op = gap_model.open_generalized;

    // Standard DP states
    let mut h = vec![vec![0.0f64; m + 1]; n + 1];
    let mut d = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut ins = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];

    // Traceback: 0=diag, 1=del, 2=ins, 3=gen_del, 4=gen_ins, 5=stop
    let mut traceback = vec![vec![5u8; m + 1]; n + 1];
    // For generalized gap traceback: store jump source (i, j)
    let mut gen_src = vec![vec![(0usize, 0usize); m + 1]; n + 1];

    let mut best_score = 0.0f64;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    // Running max for generalized gaps
    let mut tbk_col: Vec<(f64, usize)> = vec![(f64::NEG_INFINITY, 0); m + 1]; // (score, src_i)

    for i in 1..=n {
        let mut tbk_row: (f64, usize) = (f64::NEG_INFINITY, 0); // (score, src_j)

        for j in 1..=m {
            let sub = score_pair(seq1[i - 1], seq2[j - 1], matrix, amino_map);

            let diag = h[i - 1][j - 1] + sub;

            d[i][j] = (d[i - 1][j] + gap_ext).max(h[i - 1][j] + gap_open);
            ins[i][j] = (ins[i][j - 1] + gap_ext).max(h[i][j - 1] + gap_open);

            // Generalized: jump from best previous in same column (del) or row (ins)
            let gd_score = tbk_col[j].0 + gap_op;
            let gi_score = tbk_row.0 + gap_op;

            let mut wm = diag;
            let mut tb = 0u8;

            if d[i][j] > wm { wm = d[i][j]; tb = 1; }
            if ins[i][j] > wm { wm = ins[i][j]; tb = 2; }
            if gd_score > wm {
                wm = gd_score;
                tb = 3;
                gen_src[i][j] = (tbk_col[j].1, j); // jump from (src_i, j)
            }
            if gi_score > wm {
                wm = gi_score;
                tb = 4;
                gen_src[i][j] = (i, tbk_row.1); // jump from (i, src_j)
            }

            if wm < 0.0 {
                wm = 0.0;
                tb = 5;
                d[i][j] = f64::NEG_INFINITY;
                ins[i][j] = f64::NEG_INFINITY;
            }

            h[i][j] = wm;
            traceback[i][j] = tb;

            if wm > best_score {
                best_score = wm;
                best_i = i;
                best_j = j;
            }

            // Update running max (store the diagonal predecessor's position)
            if h[i - 1][j - 1] > tbk_col[j].0 {
                tbk_col[j] = (h[i - 1][j - 1], i - 1);
            }
            if h[i - 1][j - 1] > tbk_row.0 {
                tbk_row = (h[i - 1][j - 1], j - 1);
            }
        }
    }

    // Traceback
    let mut ops = Vec::new();
    let mut i = best_i;
    let mut j = best_j;

    while i > 0 && j > 0 {
        match traceback[i][j] {
            0 => {
                ops.push(AlignOp::Match);
                i -= 1;
                j -= 1;
            }
            1 => {
                ops.push(AlignOp::Delete);
                i -= 1;
            }
            2 => {
                ops.push(AlignOp::Insert);
                j -= 1;
            }
            3 => {
                // Generalized deletion: emit gaps in seq2 from i down to src_i
                let (src_i, _) = gen_src[i][j];
                while i > src_i + 1 {
                    ops.push(AlignOp::Delete);
                    i -= 1;
                }
                // Land on diagonal of source
                ops.push(AlignOp::Match);
                i -= 1;
                j -= 1;
            }
            4 => {
                // Generalized insertion: emit gaps in seq1 from j down to src_j
                let (_, src_j) = gen_src[i][j];
                while j > src_j + 1 {
                    ops.push(AlignOp::Insert);
                    j -= 1;
                }
                ops.push(AlignOp::Match);
                i -= 1;
                j -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();

    let offset1 = i;
    let offset2 = j;

    let mut a1 = Vec::with_capacity(ops.len());
    let mut a2 = Vec::with_capacity(ops.len());
    let mut si = offset1;
    let mut sj = offset2;

    for &op in &ops {
        match op {
            AlignOp::Match => {
                a1.push(seq1[si]);
                a2.push(seq2[sj]);
                si += 1;
                sj += 1;
            }
            AlignOp::Delete => {
                a1.push(seq1[si]);
                a2.push(b'-');
                si += 1;
            }
            AlignOp::Insert => {
                a1.push(b'-');
                a2.push(seq2[sj]);
                sj += 1;
            }
        }
    }

    LocalAlignment {
        alignment: Alignment {
            seq1: a1,
            seq2: a2,
            score: best_score,
            operations: ops,
        },
        offset1,
        offset2,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_matrix() -> (Vec<Vec<i32>>, [u8; 256]) {
        let mut mtx = vec![vec![-100i32; 5]; 5];
        for i in 0..4 {
            mtx[i][i] = 100;
        }
        let mut map = [0xFFu8; 256];
        map[b'A' as usize] = 0;
        map[b'C' as usize] = 1;
        map[b'G' as usize] = 2;
        map[b'T' as usize] = 3;
        map[b'-' as usize] = 4;
        (mtx, map)
    }

    #[test]
    fn identical_sequences() {
        let (mtx, map) = simple_matrix();
        let gap = GenAffineGapModel {
            affine: GapModel::new(-200.0, -10.0),
            open_generalized: -300.0,
        };
        let result = genaffine_local_align(b"ACGT", b"ACGT", &mtx, &map, &gap);
        assert!((result.alignment.score - 400.0).abs() < 1e-6);
    }

    #[test]
    fn handles_long_insertion() {
        let (mtx, map) = simple_matrix();
        let gap = GenAffineGapModel {
            affine: GapModel::new(-200.0, -50.0),
            open_generalized: -150.0,
        };
        let s1 = b"ACGTACGT";
        let s2 = b"ACGTTTTTTTTTTTTTACGT";
        let result = genaffine_local_align(s1, s2, &mtx, &map, &gap);
        assert!(result.alignment.score > 200.0);
    }

    #[test]
    fn score_non_negative() {
        let (mtx, map) = simple_matrix();
        let gap = GenAffineGapModel::default();
        let result = genaffine_local_align(b"AAAA", b"CCCC", &mtx, &map, &gap);
        assert!(result.alignment.score >= 0.0);
    }
}
