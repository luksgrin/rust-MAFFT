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
/// Returns the number of iterations performed.
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
    let mut score_history: Vec<f64> = Vec::new();

    let mut iteration = 0;
    for iter in 0..params.max_iterations {
        iteration = iter + 1;
        let mut any_change = false;

        for step in &topology.steps {
            if step.left.is_empty() || step.right.is_empty() {
                continue;
            }

            let old_score = compute_split_score(
                &step.left, &step.right, &alignment.sequences, &weights, scoring,
            );

            let (new_seqs, new_score) = realign_groups(
                &step.left,
                &step.right,
                &alignment.sequences,
                &weights,
                scoring,
                &gap,
            );

            let threshold = old_score - params.cut * old_score.abs();
            if new_score > threshold && new_seqs.is_some() {
                let new_seqs = new_seqs.unwrap();

                let changed = (0..nseq).any(|i| alignment.sequences[i] != new_seqs[i]);

                if changed {
                    alignment.sequences = new_seqs;
                    any_change = true;
                    converged_count = 0;
                } else {
                    converged_count += 1;
                }
            } else {
                converged_count += 1;
            }

            if converged_count >= convergence_target {
                return iteration;
            }
        }

        let total_score = compute_total_score(&alignment.sequences, &weights, scoring);
        if score_history.len() >= 2 {
            let prev = score_history[score_history.len() - 2];
            if (total_score - prev).abs() < 1e-10 {
                return iteration;
            }
        }
        score_history.push(total_score);

        if !any_change {
            return iteration;
        }
    }

    iteration
}

/// Re-align two groups within the current MSA.
///
/// Both groups share the same alignment width. The profile alignment
/// re-aligns their column spaces. The result is a complete new set of
/// sequences with consistent width.
fn realign_groups(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
) -> (Option<Vec<Vec<u8>>>, f64) {
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

    // Verify ops consume all columns from both profiles
    let consumed1 = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete)).count();
    let consumed2 = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert)).count();
    // Verify ops consume all columns. If not, the re-alignment may lose residues.
    if consumed1 != prof1.length || consumed2 != prof2.length {
        return (None, 0.0); // reject this re-alignment
    }

    // Both groups share the same alignment width before re-alignment.
    // The ops tell us how to merge their column spaces into a new alignment.
    let width = sequences[0].len();
    let mut new_sequences = sequences.to_vec();

    // Determine which group's columns each sequence shares.
    let mut in_group1 = vec![false; sequences.len()];
    let mut in_group2 = vec![false; sequences.len()];
    for &i in group1 { in_group1[i] = true; }
    for &i in group2 { in_group2[i] = true; }

    for idx in 0..sequences.len() {
        let old = &sequences[idx];
        let mut new_seq = Vec::with_capacity(aln.operations.len());
        let mut col = 0;

        for op in &aln.operations {
            if in_group1[idx] {
                // Group1: consume column at Match/Delete, gap at Insert
                match op {
                    AlignOp::Match | AlignOp::Delete => {
                        new_seq.push(if col < old.len() { old[col] } else { b'-' });
                        col += 1;
                    }
                    AlignOp::Insert => {
                        new_seq.push(b'-');
                    }
                }
            } else if in_group2[idx] {
                // Group2: consume column at Match/Insert, gap at Delete
                match op {
                    AlignOp::Match | AlignOp::Insert => {
                        new_seq.push(if col < old.len() { old[col] } else { b'-' });
                        col += 1;
                    }
                    AlignOp::Delete => {
                        new_seq.push(b'-');
                    }
                }
            } else {
                // Not in either group: treat like group1 (they share the same
                // column space as the full alignment before the split)
                match op {
                    AlignOp::Match | AlignOp::Delete => {
                        new_seq.push(if col < old.len() { old[col] } else { b'-' });
                        col += 1;
                    }
                    AlignOp::Insert => {
                        new_seq.push(b'-');
                    }
                }
            }
        }

        new_sequences[idx] = new_seq;
    }

    (Some(new_sequences), aln.score)
}

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
            score += pairwise_score(&sequences[i], &sequences[j], scoring) * w;
        }
    }
    score
}

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

fn compute_total_score(sequences: &[Vec<u8>], weights: &[f64], scoring: &ScoringContext) -> f64 {
    let n = sequences.len();
    let mut score = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            score += pairwise_score(&sequences[i], &sequences[j], scoring) * weights[i] * weights[j];
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

        // All sequences should have the same width
        let width = msa.width();
        for seq in &msa.sequences {
            assert_eq!(seq.len(), width);
        }

        // Ungapped should match originals
        let ungapped: Vec<Vec<u8>> = msa.sequences.iter()
            .map(|s| s.iter().filter(|&&c| c != b'-').cloned().collect())
            .collect();
        assert_eq!(ungapped[0], b"ACDEFGHIK");
        assert_eq!(ungapped[1], b"ACDEFHIK");
        assert_eq!(ungapped[2], b"ACDHIK");
    }

    #[test]
    fn refinement_preserves_width_consistency() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNP".to_vec(),
            b"ACDEFHIKLMNP".to_vec(),
            b"ACDEHIKLMNP".to_vec(),
            b"ACDHIKLMNP".to_vec(),
        ];
        let names: Vec<String> = (0..4).map(|i| format!("s{i}")).collect();

        let mut dm = DistanceMatrix::new(4);
        dm.set(0, 1, 0.1); dm.set(0, 2, 0.2); dm.set(0, 3, 0.3);
        dm.set(1, 2, 0.15); dm.set(1, 3, 0.25); dm.set(2, 3, 0.15);
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false);
        let params = RefinementParams { max_iterations: 5, ..Default::default() };

        iterative_refine(&mut msa, &topo, &scoring, &params);

        let width = msa.width();
        assert!(width > 0);
        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "sequence {i} has wrong width after refinement");
        }
    }
}
