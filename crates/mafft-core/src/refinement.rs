/// Iterative refinement (tree-dependent iteration).
///
/// Ports the C `TreeDependentIteration()` from tditeration.c.
///
/// Repeatedly re-aligns pairs of groups defined by the guide tree,
/// accepting improvements and rejecting regressions, until convergence.
///
/// Key insight from the C code: at each tree branch, ALL sequences are
/// split into two groups (subtree vs everything else). There are never
/// "uninvolved" sequences — every sequence is in one group or the other.
///
/// Branch enumeration matches C exactly:
/// - For each topology step, both sides (k=0: left, k=1: right) are
///   processed, EXCEPT the root step (last step) where only k=1 (right)
///   is used (since at the root, left-vs-complement and right-vs-complement
///   produce the same split, just flipped).
/// - Even iterations traverse steps forward (0 → N-1), odd iterations
///   traverse backward (N-1 → 0). Within each step, k always goes 0→1.
/// - Total branches per iteration: (nseq-1)*2 - 1.

use rayon::prelude::*;

use mafft_align::{
    profile_align, constrained_profile_align, ConstrainedAlignParams,
    Profile, GapModel, AlignOp,
};
use mafft_fft::SegmentParams;
use mafft_tree::{Topology, sequence_weights};
use mafft_types::{ScoringContext, LocalHomologyTable};

use crate::progressive::MultipleAlignment;

/// Parameters controlling iterative refinement.
#[derive(Debug, Clone)]
pub struct RefinementParams {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Score improvement threshold (fraction of old score).
    /// C default is 0.0 (accept only strict improvements).
    pub cut: f64,
    /// Whether to use FFT-accelerated alignment during refinement.
    pub use_fft: bool,
}

impl Default for RefinementParams {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            cut: 0.0,
            use_fft: false,
        }
    }
}

/// A branch identifier for oscillation tracking: (step_index, side).
/// side 0 = left, side 1 = right.
type BranchId = (usize, usize);

/// Iteratively refine a multiple alignment.
///
/// At each tree branch, splits ALL sequences into two groups (subtree vs
/// rest), re-aligns the two groups, and accepts improvements.
pub fn iterative_refine(
    alignment: &mut MultipleAlignment,
    topology: &Topology,
    scoring: &ScoringContext,
    params: &RefinementParams,
    constraints: Option<&LocalHomologyTable>,
) -> usize {
    let nseq = alignment.nseq();
    if nseq <= 2 || topology.steps.is_empty() {
        return 0;
    }

    let weights = sequence_weights(topology);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let mut converged_count = 0usize;
    let convergence_target = nseq * 2;

    let nsteps = topology.steps.len();
    let root_idx = nsteps - 1;
    let all_indices: Vec<usize> = (0..nseq).collect();

    // Pre-compute branch splits for all (step, side) pairs.
    // For each step, side 0 = left vs complement, side 1 = right vs complement.
    // Root step only has side 1 (side 0 is redundant at the root).
    let mut branch_map: Vec<Vec<(usize, Vec<usize>, Vec<usize>)>> = Vec::with_capacity(nsteps);
    for (step_idx, step) in topology.steps.iter().enumerate() {
        let is_root = step_idx == root_idx;
        let mut sides = Vec::new();

        if !is_root {
            // Side 0: step.left vs complement
            let complement: Vec<usize> = all_indices
                .iter()
                .filter(|i| !step.left.contains(i))
                .copied()
                .collect();
            if !step.left.is_empty() && !complement.is_empty() {
                sides.push((0, step.left.clone(), complement));
            }
        }

        // Side 1: step.right vs complement
        let complement: Vec<usize> = all_indices
            .iter()
            .filter(|i| !step.right.contains(i))
            .copied()
            .collect();
        if !step.right.is_empty() && !complement.is_empty() {
            sides.push((1, step.right.clone(), complement));
        }

        branch_map.push(sides);
    }

    // Per-branch score history for oscillation detection.
    // history[iteration][(step_idx, side)] = score after processing that branch.
    let mut history: Vec<std::collections::HashMap<BranchId, f64>> = Vec::new();

    let mut iteration = 0;
    for iter in 0..params.max_iterations {
        iteration = iter + 1;
        let mut any_change = false;
        let mut iter_scores: std::collections::HashMap<BranchId, f64> = std::collections::HashMap::new();

        // C alternates step traversal direction: even → forward, odd → reverse.
        let step_order: Vec<usize> = if iter % 2 == 0 {
            (0..nsteps).collect()
        } else {
            (0..nsteps).rev().collect()
        };

        for &step_idx in &step_order {
            for (side, group1, group2) in &branch_map[step_idx] {
                let branch_id: BranchId = (step_idx, *side);

                let old_score = compute_split_score(
                    group1, group2, &alignment.sequences, &weights, scoring,
                );

                let new_seqs = realign_all(
                    group1, group2, &alignment.sequences, &weights, scoring, &gap,
                    constraints,
                );

                if let Some((new_seqs, _new_score)) = new_seqs {
                    // C's identity check: compare representative sequences
                    let changed = (0..nseq).any(|i| alignment.sequences[i] != new_seqs[i]);

                    if !changed {
                        // Identical — no change, count toward convergence
                        let tscore = old_score;
                        iter_scores.insert(branch_id, tscore);
                        converged_count += 1;
                    } else {
                        // Compute score of the new alignment for this split
                        let tscore = compute_split_score(
                            group1, group2, &new_seqs, &weights, scoring,
                        );

                        let threshold = old_score - params.cut / 100.0 * old_score;
                        if tscore > threshold {
                            // Accept
                            alignment.sequences = new_seqs;
                            any_change = true;
                            converged_count = 0;
                        } else {
                            // Reject
                            converged_count += 1;
                        }
                        iter_scores.insert(branch_id, tscore);
                    }
                } else {
                    iter_scores.insert(branch_id, old_score);
                    converged_count += 1;
                }

                if converged_count >= convergence_target {
                    return iteration;
                }

                // Oscillation detection: check if this branch's score matches
                // the score from 2, 4, 6... iterations ago (same branch).
                if iter >= 2 {
                    let tscore = iter_scores[&branch_id];
                    let mut oscillating = false;
                    let mut ii = history.len() as isize - 2; // iterate-2
                    while ii >= 0 {
                        if let Some(&prev_score) = history[ii as usize].get(&branch_id) {
                            if tscore == prev_score {
                                oscillating = true;
                                break;
                            }
                        }
                        ii -= 2;
                    }
                    if oscillating {
                        return iteration;
                    }
                }
            }
        }

        history.push(iter_scores);

        if !any_change {
            return iteration;
        }
    }

    iteration
}

/// Re-align all sequences split into two groups.
///
/// Since group1 + group2 = ALL sequences, there are no "other" sequences
/// to worry about. Every sequence is in exactly one group.
fn realign_all(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    constraints: Option<&LocalHomologyTable>,
) -> Option<(Vec<Vec<u8>>, f64)> {
    let width = sequences[0].len();

    // Per-group gap stripping (matching C's `commongappick` during refinement).
    // Since group1 + group2 = ALL sequences, stripping both is safe:
    // no "other" sequences need column re-insertion.
    let gap1 = group_all_gap_columns(group1, sequences, width);
    let gap2 = group_all_gap_columns(group2, sequences, width);
    let kept1: Vec<usize> = (0..width).filter(|&c| !gap1[c]).collect();
    let kept2: Vec<usize> = (0..width).filter(|&c| !gap2[c]).collect();

    // Build stripped sequences
    let stripped1: Vec<Vec<u8>> = group1.iter()
        .map(|&i| kept1.iter().map(|&c| sequences[i][c]).collect())
        .collect();
    let stripped2: Vec<Vec<u8>> = group2.iter()
        .map(|&i| kept2.iter().map(|&c| sequences[i][c]).collect())
        .collect();

    let s1_refs: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
    let s2_refs: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();

    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();
    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0; group2.len()] };

    let prof1 = Profile::from_aligned(&s1_refs, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&s2_refs, &w2n, &scoring.amino_map, scoring.nalphabets);

    if prof1.length == 0 || prof2.length == 0 {
        return None;
    }

    // Use constrained alignment if local homology table is available
    let aln = if let Some(lh_table) = constraints {
        let params = ConstrainedAlignParams {
            gap: gap.clone(),
            segment_params: SegmentParams::protein(),
            constraint_weight: 1.0,
        };
        constrained_profile_align(
            &prof1, &prof2,
            &scoring.substitution_matrix,
            lh_table,
            group1, group2,
            &params,
        )
    } else {
        profile_align(&prof1, &prof2, &scoring.substitution_matrix, gap, true, true)
    };

    // Verify ops consume all columns from both profiles.
    let consumed1 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Delete))
        .count();
    let consumed2 = aln.operations.iter()
        .filter(|op| matches!(op, AlignOp::Match | AlignOp::Insert))
        .count();

    if consumed1 != prof1.length || consumed2 != prof2.length {
        eprintln!(
            "WARN: ops mismatch: consumed1={} prof1.len={} consumed2={} prof2.len={}",
            consumed1, prof1.length, consumed2, prof2.length
        );
        return None;
    }

    // Build new sequences using column mappings from stripped profiles.
    // Since ALL sequences are in one of the two groups, no re-insertion needed.
    let mut new_sequences = vec![Vec::with_capacity(aln.operations.len()); sequences.len()];
    let mut cursor1 = 0usize;
    let mut cursor2 = 0usize;

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                let oc1 = kept1[cursor1];
                let oc2 = kept2[cursor2];
                for &i in group1 { new_sequences[i].push(sequences[i][oc1]); }
                for &i in group2 { new_sequences[i].push(sequences[i][oc2]); }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                let oc1 = kept1[cursor1];
                for &i in group1 { new_sequences[i].push(sequences[i][oc1]); }
                for &i in group2 { new_sequences[i].push(b'-'); }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                let oc2 = kept2[cursor2];
                for &i in group1 { new_sequences[i].push(b'-'); }
                for &i in group2 { new_sequences[i].push(sequences[i][oc2]); }
                cursor2 += 1;
            }
        }
    }

    Some((new_sequences, aln.score))
}

fn group_all_gap_columns(group: &[usize], sequences: &[Vec<u8>], width: usize) -> Vec<bool> {
    let mut all_gap = vec![true; width];
    for &idx in group {
        for (col, &ch) in sequences[idx].iter().enumerate() {
            if ch != b'-' && ch != b'.' {
                all_gap[col] = false;
            }
        }
    }
    all_gap
}

fn compute_split_score(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
) -> f64 {
    group1
        .par_iter()
        .map(|&i| {
            group2
                .iter()
                .map(|&j| pairwise_score(&sequences[i], &sequences[j], scoring) * weights[i] * weights[j])
                .sum::<f64>()
        })
        .sum()
}

/// Branchless pairwise scoring for auto-vectorization.
///
/// The gap check is converted to a mask multiply: if either residue is '-',
/// the score contribution is 0. This eliminates branches that prevent SIMD.
#[inline]
fn pairwise_score(seq1: &[u8], seq2: &[u8], scoring: &ScoringContext) -> f64 {
    let map = &scoring.amino_map;
    let mtx = &scoring.substitution_matrix;
    let mtx_size = mtx.len();

    // Accumulate in i64 to avoid f64 conversion per position.
    // The substitution matrix contains i32 values; summing as i64 is exact.
    let mut acc = 0i64;

    for k in 0..seq1.len().min(seq2.len()) {
        let a = seq1[k];
        let b = seq2[k];
        // Branchless: non_gap is 1 if both are not '-', else 0
        let non_gap = ((a != b'-') & (b != b'-')) as i64;
        let i = map[a as usize] as usize;
        let j = map[b as usize] as usize;
        // Bounds check is predictable (almost always true for valid sequences)
        let s = if i < mtx_size && j < mtx_size { mtx[i][j] as i64 } else { 0 };
        acc += s * non_gap;
    }

    acc as f64
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

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams {
            max_iterations: 10,
            ..Default::default()
        };

        let iters = iterative_refine(&mut msa, &topo, &scoring, &params, None);
        assert!(iters <= 10);

        let width = msa.width();
        for seq in &msa.sequences {
            assert_eq!(seq.len(), width);
        }

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

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams { max_iterations: 5, ..Default::default() };

        iterative_refine(&mut msa, &topo, &scoring, &params, None);

        let width = msa.width();
        assert!(width > 0);
        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "sequence {i} has wrong width");
        }

        // Verify residue preservation
        for (i, seq) in msa.sequences.iter().enumerate() {
            let residue_count = seq.iter().filter(|&&c| c != b'-').count();
            assert_eq!(residue_count, seqs[i].len(),
                "sequence {i} lost residues: {} vs {}", residue_count, seqs[i].len());
        }
    }

    #[test]
    fn refinement_six_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIKLMNPQR".to_vec(),
            b"ACDEFHIKLMNPQR".to_vec(),
            b"ACDEHIKLMNPQR".to_vec(),
            b"ACDHIKLMNPQR".to_vec(),
            b"ACDHIKLMNP".to_vec(),
            b"ACDHIKLM".to_vec(),
        ];
        let names: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();

        let mut dm = DistanceMatrix::new(6);
        for i in 0..6 {
            for j in (i + 1)..6 {
                dm.set(i, j, (j - i) as f64 * 0.1);
            }
        }
        let topo = upgma(&dm);

        let mut msa = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        let params = RefinementParams { max_iterations: 3, ..Default::default() };

        iterative_refine(&mut msa, &topo, &scoring, &params, None);

        let width = msa.width();
        for (i, seq) in msa.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "sequence {i} has wrong width after refinement");
            let residues = seq.iter().filter(|&&c| c != b'-').count();
            assert_eq!(residues, seqs[i].len(),
                "sequence {i} lost residues during refinement");
        }
    }
}
