/// Progressive alignment following a guide tree.
///
/// Walks the tree bottom-up. At each internal node, the two child clusters
/// are aligned using profile DP (or FFT-accelerated alignment), producing
/// a merged alignment for the parent node.

use mafft_align::{profile_align, fft_profile_align, Profile, GapModel, Alignment, AlignOp, FftAlignParams};
use mafft_tree::{Topology, JoinStep, sequence_weights};
use mafft_types::ScoringContext;

/// Result of progressive alignment: aligned sequences with gap characters.
#[derive(Debug, Clone)]
pub struct MultipleAlignment {
    /// Aligned sequences (each has the same length, with '-' for gaps).
    pub sequences: Vec<Vec<u8>>,
    /// Sequence names (parallel to `sequences`).
    pub names: Vec<String>,
    /// Alignment score from the final merge step.
    pub score: f64,
}

impl MultipleAlignment {
    /// Alignment width (number of columns).
    pub fn width(&self) -> usize {
        self.sequences.first().map_or(0, |s| s.len())
    }

    /// Number of sequences.
    pub fn nseq(&self) -> usize {
        self.sequences.len()
    }
}

/// Perform progressive alignment of sequences following a guide tree.
///
/// - `sequences`: unaligned input sequences (raw residues, no gaps).
/// - `names`: sequence names.
/// - `topology`: guide tree (sequence of join steps).
/// - `scoring`: scoring context (matrix, gap penalties, etc.).
/// - `use_fft`: if true, use FFT-accelerated alignment for large groups.
pub fn progressive_align(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
) -> MultipleAlignment {
    let nseq = sequences.len();
    if nseq == 0 {
        return MultipleAlignment {
            sequences: Vec::new(),
            names: Vec::new(),
            score: 0.0,
        };
    }
    if nseq == 1 {
        return MultipleAlignment {
            sequences: sequences.to_vec(),
            names: names.to_vec(),
            score: 0.0,
        };
    }

    // Compute sequence weights from tree
    let weights = sequence_weights(topology);

    // Current aligned sequences (start as unaligned, gaps inserted progressively)
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();
    let mut last_score = 0.0;

    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Process each join step bottom-up
    for step in &topology.steps {
        let (group1_ops, group2_ops, score) = align_two_groups(
            &step.left,
            &step.right,
            &aligned,
            &weights,
            scoring,
            &gap,
            use_fft,
        );

        // Apply alignment operations to insert gaps
        apply_alignment_to_groups(
            &step.left,
            &step.right,
            &group1_ops,
            &group2_ops,
            &mut aligned,
        );

        last_score = score;
    }

    MultipleAlignment {
        sequences: aligned,
        names: names.to_vec(),
        score: last_score,
    }
}

/// Align two groups of sequences and return the alignment operations for each.
fn align_two_groups(
    group1: &[usize],
    group2: &[usize],
    aligned: &[Vec<u8>],
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    use_fft: bool,
) -> (Vec<AlignOp>, Vec<AlignOp>, f64) {
    // Build profiles from current aligned sequences
    let seqs1: Vec<&[u8]> = group1.iter().map(|&i| aligned[i].as_slice()).collect();
    let seqs2: Vec<&[u8]> = group2.iter().map(|&i| aligned[i].as_slice()).collect();

    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();

    // Normalize weights within each group
    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1_norm: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0 / group1.len() as f64; group1.len()] };
    let w2_norm: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0 / group2.len() as f64; group2.len()] };

    let prof1 = Profile::from_aligned(&seqs1, &w1_norm, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &w2_norm, &scoring.amino_map, scoring.nalphabets);

    let aln = if use_fft && prof1.length > 80 && prof2.length > 80 {
        let fft_params = FftAlignParams {
            num_candidates: 20,
            segment_params: if scoring.seq_type.is_nucleotide() {
                mafft_fft::SegmentParams::dna()
            } else {
                mafft_fft::SegmentParams::protein()
            },
            gap: gap.clone(),
            head_gap: true,
            tail_gap: true,
            num_channels: scoring.nscoredalphabets,
        };
        fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &fft_params)
    } else {
        profile_align(&prof1, &prof2, &scoring.substitution_matrix, gap, true, true)
    };

    // Split operations into group1 and group2 gap patterns
    let mut g1_ops = Vec::new();
    let mut g2_ops = Vec::new();

    for &op in &aln.operations {
        match op {
            AlignOp::Match => {
                g1_ops.push(false); // no gap in group1
                g2_ops.push(false); // no gap in group2
            }
            AlignOp::Delete => {
                g1_ops.push(false); // group1 has a residue
                g2_ops.push(true);  // group2 gets a gap
            }
            AlignOp::Insert => {
                g1_ops.push(true);  // group1 gets a gap
                g2_ops.push(false); // group2 has a residue
            }
        }
    }

    (aln.operations.clone(), aln.operations, aln.score)
}

/// Apply alignment operations to merge two groups into a consistent MSA.
///
/// The key challenge: group1 and group2 may already have gaps from previous
/// merges, and the profile alignment operates on column-level positions
/// (including gap columns). The ops describe how to merge the two groups'
/// column spaces into one.
///
/// - Match: both groups contribute a column (1:1)
/// - Delete: group1 contributes a column, group2 gets a gap column
/// - Insert: group2 contributes a column, group1 gets a gap column
///
/// Sequences NOT in either group also need gap columns inserted at the
/// appropriate positions.
fn apply_alignment_to_groups(
    group1: &[usize],
    group2: &[usize],
    _group1_ops: &[AlignOp],
    group2_ops: &[AlignOp],
    aligned: &mut [Vec<u8>],
) {
    let ops = group2_ops;

    // Get current widths of each group (should be consistent within group)
    let width1 = group1.first().map(|&i| aligned[i].len()).unwrap_or(0);
    let width2 = group2.first().map(|&i| aligned[i].len()).unwrap_or(0);

    // Build merged columns for each group
    // For group1: walk through its columns, inserting gap columns at Insert ops
    // For group2: walk through its columns, inserting gap columns at Delete ops
    let new_width = ops.len();

    let mut involved = vec![false; aligned.len()];
    for &idx in group1 { involved[idx] = true; }
    for &idx in group2 { involved[idx] = true; }

    // Rebuild group1 sequences
    for &idx in group1 {
        let old = &aligned[idx];
        let mut new_seq = Vec::with_capacity(new_width);
        let mut col = 0;
        for op in ops {
            match op {
                AlignOp::Match | AlignOp::Delete => {
                    if col < old.len() {
                        new_seq.push(old[col]);
                    } else {
                        new_seq.push(b'-');
                    }
                    col += 1;
                }
                AlignOp::Insert => {
                    new_seq.push(b'-'); // gap column for group1
                }
            }
        }
        aligned[idx] = new_seq;
    }

    // Rebuild group2 sequences
    for &idx in group2 {
        let old = &aligned[idx];
        let mut new_seq = Vec::with_capacity(new_width);
        let mut col = 0;
        for op in ops {
            match op {
                AlignOp::Match | AlignOp::Insert => {
                    if col < old.len() {
                        new_seq.push(old[col]);
                    } else {
                        new_seq.push(b'-');
                    }
                    col += 1;
                }
                AlignOp::Delete => {
                    new_seq.push(b'-'); // gap column for group2
                }
            }
        }
        aligned[idx] = new_seq;
    }

    // Sequences not in either group: they share the column space of whichever
    // group they were previously merged with. We need to insert gap columns
    // at positions where new columns were introduced.
    // Since both groups had the same alignment width before this step
    // (they share the MSA), we insert gaps at Delete (new from group1's
    // perspective beyond group2) and Insert positions.
    //
    // Actually, before this merge, group1 and group2 may have DIFFERENT
    // widths if they haven't been merged before. Other sequences share
    // width with whichever group they were last merged with.
    //
    // The simplest correct approach: other sequences have the same width
    // as group1 (or group2). We insert gap columns where that group got gaps.
    for idx in 0..aligned.len() {
        if involved[idx] { continue; }
        let old = &aligned[idx];
        let old_len = old.len();

        // Determine which group's column space this sequence shares.
        // If it has the same width as group1, use group1's pattern.
        // Otherwise use group2's pattern.
        let mut new_seq = Vec::with_capacity(new_width);
        if old_len == width1 {
            // Same as group1: consume at Match/Delete, gap at Insert
            let mut col = 0;
            for op in ops {
                match op {
                    AlignOp::Match | AlignOp::Delete => {
                        if col < old.len() { new_seq.push(old[col]); }
                        else { new_seq.push(b'-'); }
                        col += 1;
                    }
                    AlignOp::Insert => { new_seq.push(b'-'); }
                }
            }
        } else if old_len == width2 {
            // Same as group2: consume at Match/Insert, gap at Delete
            let mut col = 0;
            for op in ops {
                match op {
                    AlignOp::Match | AlignOp::Insert => {
                        if col < old.len() { new_seq.push(old[col]); }
                        else { new_seq.push(b'-'); }
                        col += 1;
                    }
                    AlignOp::Delete => { new_seq.push(b'-'); }
                }
            }
        } else {
            // Width doesn't match either group — this can happen if the
            // sequence hasn't been involved in any merge yet (raw sequence).
            // Treat it like group1 pattern and pad.
            let mut col = 0;
            for op in ops {
                match op {
                    AlignOp::Match | AlignOp::Delete => {
                        if col < old.len() { new_seq.push(old[col]); }
                        else { new_seq.push(b'-'); }
                        col += 1;
                    }
                    AlignOp::Insert => { new_seq.push(b'-'); }
                }
            }
        }
        aligned[idx] = new_seq;
    }
}

/// Insert gaps into a sequence according to a boolean pattern.
/// `gaps[i] == true` means insert a gap at position i in the output.
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

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_tree::{DistanceMatrix, upgma};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};

    #[test]
    fn progressive_two_identical() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![b"ACDEFGHIK".to_vec(), b"ACDEFGHIK".to_vec()];
        let names = vec!["s1".into(), "s2".into()];

        let mut dm = DistanceMatrix::new(2);
        dm.set(0, 1, 0.0);
        let topo = upgma(&dm);

        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
        assert_eq!(result.nseq(), 2);
        assert_eq!(result.sequences[0], result.sequences[1]);
    }

    #[test]
    fn progressive_three_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFHIK".to_vec(),  // missing G
            b"ACDHIK".to_vec(),    // missing EFG
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];

        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1);
        dm.set(0, 2, 0.3);
        dm.set(1, 2, 0.2);
        let topo = upgma(&dm);

        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
        assert_eq!(result.nseq(), 3);

        // All aligned sequences should have the same length
        let width = result.width();
        for seq in &result.sequences {
            assert_eq!(seq.len(), width, "all sequences should have same length");
        }
    }
}
