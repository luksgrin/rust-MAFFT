/// Optimal anchor pair selection via DP.
///
/// Ports the C `blockAlign2()` from fftFunctions.c.
///
/// Given cross-scores between segment pairs from two sequence groups,
/// selects the optimal non-overlapping subset of anchor pairs that
/// maximizes the total alignment score.

/// Select optimal anchor pairs from a cross-score matrix.
///
/// - `cross_scores[i][j]`: score of pairing segment `i` from group1 with
///   segment `j` from group2.
/// - `gap_penalty`: cost for skipping segments (non-diagonal moves).
///
/// Returns the selected pairs as `(cut1_indices, cut2_indices)`.
pub fn block_align(
    cross_scores: &[Vec<f64>],
    gap_penalty: f64,
) -> (Vec<usize>, Vec<usize>) {
    let ncut = cross_scores.len();
    if ncut == 0 {
        return (Vec::new(), Vec::new());
    }
    if ncut == 1 {
        return (vec![0], vec![0]);
    }

    // DP: score[i][j] = best cumulative score up to pairing segment i with j
    let mut dp = vec![vec![0.0f64; ncut]; ncut];
    // Track: 0 = diagonal, positive = skip in j, negative = skip in i
    let mut track = vec![vec![0i32; ncut]; ncut];

    // Initialize first row and column
    for i in 0..ncut {
        for j in 0..ncut {
            dp[i][j] = cross_scores[i][j];
        }
    }

    // Fill DP
    for i in 1..ncut {
        for j in 1..ncut {
            // Diagonal: continue from (i-1, j-1)
            let mut best = dp[i - 1][j - 1];
            track[i][j] = 0;

            // Skip segments in j (gap in group2's segments)
            for k in 0..j.saturating_sub(1) {
                let score = dp[i - 1][k] + gap_penalty;
                if score > best {
                    best = score;
                    track[i][j] = (j - k) as i32;
                }
            }

            // Skip segments in i (gap in group1's segments)
            for k in 0..i.saturating_sub(1) {
                let score = dp[k][j - 1] + gap_penalty;
                if score > best {
                    best = score;
                    track[i][j] = -((i - k) as i32);
                }
            }

            dp[i][j] = cross_scores[i][j] + best;
        }
    }

    // Traceback from (ncut-1, ncut-1)
    let mut result_i = Vec::new();
    let mut result_j = Vec::new();

    let mut i = ncut - 1;
    let mut j = ncut - 1;

    loop {
        if cross_scores[i][j] > 0.0 {
            result_i.push(i);
            result_j.push(j);
        }

        let shift = track[i][j];
        if shift == 0 {
            if i == 0 || j == 0 { break; }
            i -= 1;
            j -= 1;
        } else if shift > 0 {
            if i == 0 { break; }
            j -= shift as usize;
            i -= 1;
            if j == 0 { break; }
        } else {
            if j == 0 { break; }
            i -= (-shift) as usize;
            j -= 1;
            if i == 0 { break; }
        }
    }

    result_i.reverse();
    result_j.reverse();
    (result_i, result_j)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diagonal_scores_selected() {
        // Strong diagonal, weak off-diagonal
        let scores = vec![
            vec![10.0, 0.0, 0.0],
            vec![0.0, 10.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let (ci, cj) = block_align(&scores, -5.0);
        assert_eq!(ci, vec![0, 1, 2]);
        assert_eq!(cj, vec![0, 1, 2]);
    }

    #[test]
    fn skips_zero_score_pairs() {
        let scores = vec![
            vec![10.0, 0.0],
            vec![0.0, 0.0],
        ];
        let (ci, cj) = block_align(&scores, -5.0);
        // Should only include the non-zero pair
        assert!(ci.contains(&0));
        assert!(!ci.iter().zip(cj.iter()).any(|(&i, &j)| scores[i][j] == 0.0 && i > 0));
    }

    #[test]
    fn empty_input() {
        let (ci, cj) = block_align(&[], -5.0);
        assert!(ci.is_empty());
        assert!(cj.is_empty());
    }
}
