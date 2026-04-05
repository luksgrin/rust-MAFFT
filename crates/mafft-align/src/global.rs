/// Global alignment (Needleman-Wunsch) with affine gap penalties.
///
/// Ports the C `G__align11()` from Galign11.c.

use crate::dp::{score_pair, AlignOp, Alignment, GapModel};

/// Perform global alignment of two sequences.
///
/// - `seq1`, `seq2`: raw residue sequences (uppercase ASCII).
/// - `matrix`: substitution score matrix (alphabet × alphabet).
/// - `amino_map`: ASCII char → internal index (256-element lookup).
/// - `gap`: affine gap model (open + extend).
/// - `head_gap`: if true, penalize terminal gaps at the start.
/// - `tail_gap`: if true, penalize terminal gaps at the end.
pub fn global_align(
    seq1: &[u8],
    seq2: &[u8],
    matrix: &[Vec<i32>],
    amino_map: &[u8; 256],
    gap: &GapModel,
    head_gap: bool,
    tail_gap: bool,
) -> Alignment {
    let n = seq1.len();
    let m = seq2.len();

    if n == 0 || m == 0 {
        return empty_alignment(seq1, seq2);
    }

    let head_open = if head_gap { gap.open } else { 0.0 };
    let head_ext = if head_gap { gap.extend } else { 0.0 };

    // DP matrices (n+1) x (m+1)
    // H[i][j] = best score aligning seq1[0..i] and seq2[0..j]
    // D[i][j] = best score ending with gap in seq2 (deletion from seq1)
    // I[i][j] = best score ending with gap in seq1 (insertion)
    let mut h = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut d = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];
    let mut ins = vec![vec![f64::NEG_INFINITY; m + 1]; n + 1];

    // Traceback: 0=diag, 1=from D (del), 2=from I (ins)
    let mut traceback = vec![vec![0u8; m + 1]; n + 1];

    // Initialization
    h[0][0] = 0.0;
    for i in 1..=n {
        d[i][0] = head_open + head_ext * i as f64;
        h[i][0] = d[i][0];
        traceback[i][0] = 1;
    }
    for j in 1..=m {
        ins[0][j] = head_open + head_ext * j as f64;
        h[0][j] = ins[0][j];
        traceback[0][j] = 2;
    }

    // Fill
    for i in 1..=n {
        for j in 1..=m {
            let sub = score_pair(seq1[i - 1], seq2[j - 1], matrix, amino_map);

            // Diagonal (match/mismatch)
            let diag = h[i - 1][j - 1] + sub;

            // Deletion: gap in seq2
            let d_extend = d[i - 1][j] + gap.extend;
            let d_open = h[i - 1][j] + gap.open;
            d[i][j] = d_extend.max(d_open);

            // Insertion: gap in seq1
            let i_extend = ins[i][j - 1] + gap.extend;
            let i_open = h[i][j - 1] + gap.open;
            ins[i][j] = i_extend.max(i_open);

            // Best path
            h[i][j] = diag;
            traceback[i][j] = 0;

            if d[i][j] > h[i][j] {
                h[i][j] = d[i][j];
                traceback[i][j] = 1;
            }
            if ins[i][j] > h[i][j] {
                h[i][j] = ins[i][j];
                traceback[i][j] = 2;
            }
        }
    }

    // Find best endpoint (considering tail gap penalties)
    let (mut ei, mut ej) = (n, m);
    let mut best_score = h[n][m];

    if !tail_gap {
        // Check last row (free end gaps in seq2)
        for j in 0..m {
            if h[n][j] > best_score {
                best_score = h[n][j];
                ei = n;
                ej = j;
            }
        }
        // Check last column (free end gaps in seq1)
        for i in 0..n {
            if h[i][m] > best_score {
                best_score = h[i][m];
                ei = i;
                ej = m;
            }
        }
    }

    // Traceback
    let mut ops = Vec::new();
    let mut i = ei;
    let mut j = ej;

    // Add trailing gaps if we didn't end at (n, m)
    while i < n {
        ops.push(AlignOp::Delete);
        i += 1;
    }
    while j < m {
        ops.push(AlignOp::Insert);
        j += 1;
    }

    i = ei;
    j = ej;
    while i > 0 || j > 0 {
        match traceback[i][j] {
            0 => {
                // Diagonal
                if i == 0 || j == 0 {
                    break;
                }
                ops.push(AlignOp::Match);
                i -= 1;
                j -= 1;
            }
            1 => {
                // Deletion (gap in seq2)
                ops.push(AlignOp::Delete);
                i -= 1;
            }
            2 => {
                // Insertion (gap in seq1)
                ops.push(AlignOp::Insert);
                j -= 1;
            }
            _ => break,
        }
    }

    ops.reverse();

    // Build aligned sequences
    let (aligned1, aligned2) = build_aligned_seqs(seq1, seq2, &ops);

    Alignment {
        seq1: aligned1,
        seq2: aligned2,
        score: best_score,
        operations: ops,
    }
}

fn build_aligned_seqs(seq1: &[u8], seq2: &[u8], ops: &[AlignOp]) -> (Vec<u8>, Vec<u8>) {
    let mut a1 = Vec::with_capacity(ops.len());
    let mut a2 = Vec::with_capacity(ops.len());
    let mut i = 0;
    let mut j = 0;

    for &op in ops {
        match op {
            AlignOp::Match => {
                a1.push(seq1[i]);
                a2.push(seq2[j]);
                i += 1;
                j += 1;
            }
            AlignOp::Delete => {
                a1.push(seq1[i]);
                a2.push(b'-');
                i += 1;
            }
            AlignOp::Insert => {
                a1.push(b'-');
                a2.push(seq2[j]);
                j += 1;
            }
        }
    }
    (a1, a2)
}

fn empty_alignment(seq1: &[u8], seq2: &[u8]) -> Alignment {
    let mut a1 = Vec::new();
    let mut a2 = Vec::new();
    let mut ops = Vec::new();

    for &c in seq1 {
        a1.push(c);
        a2.push(b'-');
        ops.push(AlignOp::Delete);
    }
    for &c in seq2 {
        a1.push(b'-');
        a2.push(c);
        ops.push(AlignOp::Insert);
    }
    Alignment {
        seq1: a1,
        seq2: a2,
        score: 0.0,
        operations: ops,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn simple_matrix() -> (Vec<Vec<i32>>, [u8; 256]) {
        // Simple 5x5 matrix: A=0, C=1, G=2, T=3, -=4
        // Match=100, mismatch=-100
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
        let gap = GapModel::new(-200.0, -10.0);
        let aln = global_align(b"ACGT", b"ACGT", &mtx, &map, &gap, true, true);

        assert_eq!(aln.seq1, b"ACGT");
        assert_eq!(aln.seq2, b"ACGT");
        assert!((aln.score - 400.0).abs() < 1e-6); // 4 * 100
        assert!((aln.identity() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn simple_gap() {
        let (mtx, map) = simple_matrix();
        let gap = GapModel::new(-150.0, -10.0);
        // ACGT vs AGT: best alignment inserts gap in seq2
        let aln = global_align(b"ACGT", b"AGT", &mtx, &map, &gap, true, true);

        assert_eq!(aln.seq1.len(), aln.seq2.len());
        // Should have 3 matches and 1 gap
        let gaps: usize = aln.seq2.iter().filter(|&&c| c == b'-').count()
            + aln.seq1.iter().filter(|&&c| c == b'-').count();
        assert!(gaps >= 1);
    }

    #[test]
    fn alignment_is_valid() {
        let (mtx, map) = simple_matrix();
        let gap = GapModel::new(-200.0, -10.0);
        let aln = global_align(b"ACGTACGT", b"ACGACG", &mtx, &map, &gap, true, true);

        // Aligned sequences must have equal length
        assert_eq!(aln.seq1.len(), aln.seq2.len());

        // Removing gaps from aligned seq1 should give original
        let ungapped1: Vec<u8> = aln.seq1.iter().filter(|&&c| c != b'-').cloned().collect();
        let ungapped2: Vec<u8> = aln.seq2.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(ungapped1, b"ACGTACGT");
        assert_eq!(ungapped2, b"ACGACG");
    }
}
