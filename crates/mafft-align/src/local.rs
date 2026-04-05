/// Local alignment (Smith-Waterman) with affine gap penalties.
///
/// Ports the C `L__align11()` from Lalign11.c.

use crate::dp::{score_pair, AlignOp, Alignment, GapModel};

/// Result of a local alignment, including the offsets into original sequences.
#[derive(Debug, Clone)]
pub struct LocalAlignment {
    /// The alignment itself.
    pub alignment: Alignment,
    /// Start position in sequence 1 (0-based).
    pub offset1: usize,
    /// Start position in sequence 2 (0-based).
    pub offset2: usize,
}

/// Perform local alignment of two sequences.
///
/// Returns the highest-scoring local alignment. Unlike global alignment,
/// the alignment can start and end at any position in either sequence.
///
/// - `score_offset`: additional offset subtracted from all scores (shifts
///   the effective zero threshold for local alignment termination).
///   Corresponds to the C `scoreoffset` parameter.
pub fn local_align(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    score_offset: f64,
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

    // The local threshold: scores below this reset to 0
    let local_thr = score_offset * 600.0;

    // DP matrices
    let mut h = vec![vec![0.0f64; m + 1]; n + 1];
    let mut d = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut ins = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];

    // 0=diag, 1=del, 2=ins, 3=local_stop (reset to 0)
    let mut traceback = vec![vec![3u8; m + 1]; n + 1];

    let mut best_score = 0.0f64;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=n {
        for j in 1..=m {
            let sub = score_pair(seq1[i - 1], seq2[j - 1], matrix, amino_map);

            // Diagonal
            let diag = h[i - 1][j - 1] + sub;

            // Deletion
            d[i][j] = (d[i - 1][j] + gap.extend).max(h[i - 1][j] + gap.open);

            // Insertion
            ins[i][j] = (ins[i][j - 1] + gap.extend).max(h[i][j - 1] + gap.open);

            // Best path (local: floor at local_thr)
            let mut wm = diag;
            let mut tb = 0u8;

            if d[i][j] > wm {
                wm = d[i][j];
                tb = 1;
            }
            if ins[i][j] > wm {
                wm = ins[i][j];
                tb = 2;
            }

            if wm < local_thr {
                // Local stop: reset
                wm = 0.0;
                tb = 3;
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
        }
    }

    // Traceback from best position
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
            _ => break, // local stop
        }
    }

    ops.reverse();

    let offset1 = i;
    let offset2 = j;

    // Build aligned sequences
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
    fn finds_local_match() {
        let (mtx, map) = simple_matrix();
        let gap = GapModel::new(-200.0, -10.0);
        // Embed a matching region in mismatching flanks
        let s1 = b"TTTTACGTTTTT";
        let s2 = b"GGGGACGTGGGG";
        let result = local_align(s1, s2, &mtx, &map, &gap, 0.0);

        // Should find the ACGT match
        let ungapped1: Vec<u8> = result.alignment.seq1.iter().filter(|&&c| c != b'-').cloned().collect();
        let ungapped2: Vec<u8> = result.alignment.seq2.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(ungapped1, b"ACGT");
        assert_eq!(ungapped2, b"ACGT");
        assert_eq!(result.offset1, 4);
        assert_eq!(result.offset2, 4);
    }

    #[test]
    fn local_score_non_negative() {
        let (mtx, map) = simple_matrix();
        let gap = GapModel::new(-200.0, -10.0);
        let result = local_align(b"AAAA", b"CCCC", &mtx, &map, &gap, 0.0);
        assert!(result.alignment.score >= 0.0);
    }

    #[test]
    fn identical_sequences_full_match() {
        let (mtx, map) = simple_matrix();
        let gap = GapModel::new(-200.0, -10.0);
        let result = local_align(b"ACGT", b"ACGT", &mtx, &map, &gap, 0.0);
        assert!((result.alignment.score - 400.0).abs() < 1e-6);
        assert_eq!(result.offset1, 0);
        assert_eq!(result.offset2, 0);
    }
}
