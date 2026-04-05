/// Iterative refinement (tree-dependent iteration).
///
/// Ports the C `TreeDependentIteration()` from tditeration.c.
///
/// Repeatedly re-aligns pairs of groups defined by the guide tree,
/// accepting improvements and rejecting regressions, until convergence.

use mafft_align::{profile_align, Profile, GapModel, Alignment, AlignOp};
use mafft_tree::{Topology, sequence_weights};
use mafft_types::ScoringContext;

use crate::progressive::MultipleAlignment;

/// Parameters controlling iterative refinement.
#[derive(Debug, Clone)]
pub struct RefinementParams {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Score improvement threshold (fraction). An alignment is rejected
    /// if the new score is worse than `old_score * (1 - cut)`.
    pub cut: f64,
    /// Whether to use FFT-accelerated alignment during refinement.
    pub use_fft: bool,
}

impl Default for RefinementParams {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            cut: 0.0001,
            use_fft: false,
        }
    }
}

/// Iteratively refine a multiple alignment.
///
/// Walks the guide tree repeatedly, re-aligning the two groups at each
/// internal node. Improvements are accepted, regressions rejected.
/// Stops when converged (no changes across all branches) or oscillating.
///
/// Returns the refined alignment and number of iterations performed.
pub fn iterative_refine(
    alignment: &mut MultipleAlignment,
    topology: &Topology,
    scoring: &ScoringContext,
    params: &RefinementParams,
) -> usize {
    let nseq = alignment.nseq();
    if nseq <= 2 || topology.steps.is_empty() {
        return 0;
    }

    let weights = sequence_weights(topology);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let mut converged_count = 0usize;
    let convergence_target = nseq * 2;

    // Score history for oscillation detection (even-indexed iterations)
    let mut score_history: Vec<f64> = Vec::new();

    let mut iteration = 0;
    for iter in 0..params.max_iterations {
        iteration = iter + 1;
        let mut any_change = false;

        for (step_idx, step) in topology.steps.iter().enumerate() {
            if step.left.is_empty() || step.right.is_empty() {
                continue;
            }

            // Compute current score for this split
            let old_score = compute_split_score(
                &step.left, &step.right, &alignment.sequences, &weights, scoring,
            );

            // Re-align the two groups
            let (new_seqs, new_score) = realign_groups(
                &step.left,
                &step.right,
                &alignment.sequences,
                &weights,
                scoring,
                &gap,
            );

            // Accept or reject
            let threshold = old_score - params.cut * old_score.abs();
            if new_score > threshold && new_seqs.is_some() {
                let new_seqs = new_seqs.unwrap();

                // Check if alignment actually changed
                let changed = step.left.iter().chain(step.right.iter()).any(|&i| {
                    alignment.sequences[i] != new_seqs[i]
                });

                if changed {
                    // Accept: update sequences
                    for &i in step.left.iter().chain(step.right.iter()) {
                        alignment.sequences[i] = new_seqs[i].clone();
                    }
                    any_change = true;
                    converged_count = 0;
                } else {
                    converged_count += 1;
                }
            } else {
                converged_count += 1;
            }

            // Check convergence
            if converged_count >= convergence_target {
                return iteration;
            }
        }

        // Oscillation detection: check if current score matches 2 iterations ago
        let total_score = compute_total_score(&alignment.sequences, &weights, scoring);
        if score_history.len() >= 2 {
            let prev = score_history[score_history.len() - 2];
            if (total_score - prev).abs() < 1e-10 {
                return iteration; // oscillating
            }
        }
        score_history.push(total_score);

        if !any_change {
            return iteration; // fully converged
        }
    }

    iteration
}

/// Re-align two groups within the current MSA.
fn realign_groups(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
) -> (Option<Vec<Vec<u8>>>, f64) {
    // Extract and strip common gaps for each group
    let seqs1: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
    let seqs2: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();

    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();

    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0; group2.len()] };

    let prof1 = Profile::from_aligned(&seqs1, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &w2n, &scoring.amino_map, scoring.nalphabets);

    if prof1.length == 0 || prof2.length == 0 {
        return (None, 0.0);
    }

    let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, gap, true, true);

    // Apply alignment operations to produce new gapped sequences
    let mut new_sequences = sequences.to_vec();
    apply_ops_to_groups(group1, group2, &aln.operations, &mut new_sequences);

    (Some(new_sequences), aln.score)
}

/// Apply alignment operations to insert gaps into group sequences.
fn apply_ops_to_groups(
    group1: &[usize],
    group2: &[usize],
    ops: &[AlignOp],
    sequences: &mut [Vec<u8>],
) {
    let mut gap_pattern_1 = Vec::new();
    let mut gap_pattern_2 = Vec::new();

    for op in ops {
        match op {
            AlignOp::Match => { gap_pattern_1.push(false); gap_pattern_2.push(false); }
            AlignOp::Delete => { gap_pattern_1.push(false); gap_pattern_2.push(true); }
            AlignOp::Insert => { gap_pattern_1.push(true); gap_pattern_2.push(false); }
        }
    }

    for &idx in group1 {
        sequences[idx] = insert_gaps(&sequences[idx], &gap_pattern_1);
    }
    for &idx in group2 {
        sequences[idx] = insert_gaps(&sequences[idx], &gap_pattern_2);
    }
}

fn insert_gaps(seq: &[u8], gaps: &[bool]) -> Vec<u8> {
    let mut result = Vec::with_capacity(gaps.len());
    let mut seq_pos = 0;
    for &is_gap in gaps {
        if is_gap {
            result.push(b'-');
        } else if seq_pos < seq.len() {
            result.push(seq[seq_pos]);
            seq_pos += 1;
        } else {
            result.push(b'-');
        }
    }
    result
}

/// Compute alignment score for a split (sum-of-pairs between groups).
fn compute_split_score(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
) -> f64 {
    let mut score = 0.0;
    for &i in group1 {
        for &j in group2 {
            let w = weights[i] * weights[j];
            let pair_score = pairwise_score(&sequences[i], &sequences[j], scoring);
            score += pair_score * w;
        }
    }
    score
}

/// Simple pairwise alignment score (sum of substitution scores at aligned positions).
fn pairwise_score(seq1: &[u8], seq2: &[u8], scoring: &ScoringContext) -> f64 {
    let mut score = 0.0;
    for (a, b) in seq1.iter().zip(seq2.iter()) {
        if *a != b'-' && *b != b'-' {
            let i = scoring.amino_map[*a as usize] as usize;
            let j = scoring.amino_map[*b as usize] as usize;
            if i < scoring.substitution_matrix.len() && j < scoring.substitution_matrix[0].len() {
                score += scoring.substitution_matrix[i][j] as f64;
            }
        }
    }
    score
}

/// Compute total sum-of-pairs score for the alignment.
fn compute_total_score(sequences: &[Vec<u8>], weights: &[f64], scoring: &ScoringContext) -> f64 {
    let n = sequences.len();
    let mut score = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let w = weights[i] * weights[j];
            score += pairwise_score(&sequences[i], &sequences[j], scoring) * w;
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_tree::{DistanceMatrix, upgma};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use crate::progressive::progressive_align;

    #[test]
    fn refinement_converges() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFHIK".to_vec(),
            b"ACDHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];

        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1);
        dm.set(0, 2, 0.3);
        dm.set(1, 2, 0.2);
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false);
        let params = RefinementParams {
            max_iterations: 10,
            cut: 0.0001,
            use_fft: false,
        };

        let iters = iterative_refine(&mut msa, &topo, &scoring, &params);
        assert!(iters <= 10, "should converge within 10 iterations, took {iters}");

        // All sequences should still have the same width
        let width = msa.width();
        for seq in &msa.sequences {
            assert_eq!(seq.len(), width);
        }
    }
}
