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

    // Fill DP. Mirrors C's `blockAlign2` (`fftFunctions.c:455-509`).
    //
    // C calls `permit(seg_a, seg_b)` to gate non-corner skip transitions —
    // but the production C code defines `permit` as `return( 0 );` followed
    // by dead code (`fftFunctions.c:378-384`). With `permit == 0` the
    // condition `if( k && k<ncut-1 && j<ncut-1 && !permit(...) ) continue;`
    // simplifies to: skip iff `k != 0 AND k < ncut-1 AND j < ncut-1`. So
    // for interior cells only k=0 (corner-jump) is allowed; only boundary
    // cells (last row or column) consider full skips.
    for i in 1..ncut {
        for j in 1..ncut {
            // Diagonal: continue from (i-1, j-1)
            let mut best = dp[i - 1][j - 1];
            track[i][j] = 0;

            // Skip segments in j (gap in group2's segments).
            // C: klim = j-2; for(k=0; k<klim; k++).
            // Gating: skip iff `k != 0 && k < ncut-1 && j < ncut-1`.
            // i.e. allowed when k=0, or k=ncut-1, or j=ncut-1.
            for k in 0..j.saturating_sub(2) {
                if k != 0 && k < ncut - 1 && j < ncut - 1 { continue; }
                let score = dp[i - 1][k] + gap_penalty;
                if score > best {
                    best = score;
                    track[i][j] = (j - k) as i32;
                }
            }

            // Skip segments in i (gap in group1's segments).
            // C: klim = i-2; same gating with i in place of j.
            for k in 0..i.saturating_sub(2) {
                if k != 0 && k < ncut - 1 && i < ncut - 1 { continue; }
                let score = dp[k][j - 1] + gap_penalty;
                if score > best {
                    best = score;
                    track[i][j] = -((i - k) as i32);
                }
            }

            dp[i][j] = cross_scores[i][j] + best;
        }
    }

    // Traceback from (ncut-1, ncut-1). Mirrors `fftFunctions.c:520-545`:
    // walk back via `track` until either coordinate hits 0.
    let mut path_i: Vec<usize> = Vec::new();
    let mut path_j: Vec<usize> = Vec::new();

    let mut i = ncut - 1;
    let mut j = ncut - 1;
    path_i.push(i);
    path_j.push(j);

    loop {
        if i == 0 || j == 0 { break; }
        let shift = track[i][j];
        if shift == 0 {
            i -= 1;
            j -= 1;
        } else if shift > 0 {
            // gap in group2: came from (i-1, j-shift)
            let new_j = j.checked_sub(shift as usize);
            i -= 1;
            j = match new_j { Some(v) => v, None => break };
        } else {
            // gap in group1: came from (i+shift, j-1)
            let new_i = i.checked_sub((-shift) as usize);
            j -= 1;
            i = match new_i { Some(v) => v, None => break };
        }
        path_i.push(i);
        path_j.push(j);
    }

    path_i.reverse();
    path_j.reverse();

    // Filter + dedupe: mirrors `fftFunctions.c:547-560`. Drop cells whose
    // `cross_scores` is zero (corner sentinels stripped, unfilled cells
    // skipped), and when a kept cell shares a row OR column with the
    // previous kept cell, keep only the higher-scoring one.
    let mut result_i: Vec<usize> = Vec::new();
    let mut result_j: Vec<usize> = Vec::new();
    for (&ci, &cj) in path_i.iter().zip(path_j.iter()) {
        if cross_scores[ci][cj] == 0.0 { continue; }
        if let (Some(&prev_i), Some(&prev_j)) = (result_i.last(), result_j.last()) {
            if (ci == prev_i || cj == prev_j)
                && cross_scores[ci][cj] > cross_scores[prev_i][prev_j]
            {
                result_i.pop();
                result_j.pop();
            } else if ci == prev_i || cj == prev_j {
                // tie or worse — keep the previous, drop current
                continue;
            }
        }
        result_i.push(ci);
        result_j.push(cj);
    }

    (result_i, result_j)
}

/// Variant of block_align using running-max tracking (O(n²) time).
///
/// Ports C's `blockAlign3()`. Uses `jumpscore`/`jumppos` arrays to track
/// the best previous score in each row/column, avoiding the inner loop
/// scan of `blockAlign2`.
pub fn block_align3(
    cross_scores: &[Vec<f64>],
    gap_penalty: f64,
) -> (Vec<usize>, Vec<usize>) {
    let ncut = cross_scores.len();
    if ncut == 0 { return (Vec::new(), Vec::new()); }
    if ncut == 1 { return (vec![0], vec![0]); }

    let mut dp = vec![vec![0.0f64; ncut]; ncut];
    let mut track = vec![vec![0i32; ncut]; ncut];

    for i in 0..ncut {
        for j in 0..ncut {
            dp[i][j] = cross_scores[i][j];
        }
    }

    // Running max trackers
    let mut jumpscore_col = vec![f64::NEG_INFINITY; ncut]; // best score in column j
    let mut jumppos_col = vec![0usize; ncut];

    for i in 1..ncut {
        let mut jumpscore_row = f64::NEG_INFINITY;
        let mut jumppos_row = 0usize;

        for j in 1..ncut {
            let mut best = dp[i - 1][j - 1];
            track[i][j] = 0;

            // Jump from best previous in row (skip columns)
            let row_score = jumpscore_row + gap_penalty;
            if row_score > best {
                best = row_score;
                track[i][j] = (j - jumppos_row) as i32;
            }

            // Jump from best previous in column (skip rows)
            let col_score = jumpscore_col[j] + gap_penalty;
            if col_score > best {
                best = col_score;
                track[i][j] = -((i - jumppos_col[j]) as i32);
            }

            dp[i][j] = cross_scores[i][j] + best;

            // Update running max for row
            if dp[i - 1][j] > jumpscore_row {
                jumpscore_row = dp[i - 1][j];
                jumppos_row = j;
            }

            // Update running max for column
            if dp[i][j - 1] > jumpscore_col[j] {
                jumpscore_col[j] = dp[i][j - 1];
                jumppos_col[j] = i;
            }
        }
    }

    // Traceback (same as blockAlign2)
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
    fn block_align3_diagonal() {
        let scores = vec![
            vec![10.0, 0.0, 0.0],
            vec![0.0, 10.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let (ci, cj) = block_align3(&scores, -5.0);
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
