/// Progressive alignment following a guide tree.
///
/// Matches C's `treebase()` from `disttbfast.c`: at each merge step,
/// only group1 and group2 sequences are modified. "Other" sequences
/// (not in either group) are left untouched, matching C's behavior
/// where `mergeoralign = 'a'` means no special gap handling for others.

use mafft_align::{profile_align, fft_profile_align, Profile, GapModel, Alignment, AlignOp, FftAlignParams};
use mafft_tree::{Topology, sequence_weights};
use mafft_types::ScoringContext;

#[derive(Debug, Clone)]
pub struct MultipleAlignment {
    pub sequences: Vec<Vec<u8>>,
    pub names: Vec<String>,
    pub score: f64,
}

impl MultipleAlignment {
    pub fn width(&self) -> usize {
        self.sequences.first().map_or(0, |s| s.len())
    }
    pub fn nseq(&self) -> usize {
        self.sequences.len()
    }
}

pub fn progressive_align(
    sequences: &[Vec<u8>],
    names: &[String],
    topology: &Topology,
    scoring: &ScoringContext,
    use_fft: bool,
    shift_penalty: Option<f64>,
) -> MultipleAlignment {
    let nseq = sequences.len();
    if nseq == 0 {
        return MultipleAlignment { sequences: Vec::new(), names: Vec::new(), score: 0.0 };
    }
    if nseq == 1 {
        return MultipleAlignment { sequences: sequences.to_vec(), names: names.to_vec(), score: 0.0 };
    }

    let weights = sequence_weights(topology);
    // Start with original sequences (varying lengths, like C).
    // No padding — sequences within a group will have matching widths
    // because they were merged together at a previous step.
    let mut aligned: Vec<Vec<u8>> = sequences.to_vec();

    let mut last_score = 0.0;
    let mut gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    if let Some(shift) = shift_penalty {
        gap = gap.with_shift(shift);
    }

    for step in &topology.steps {
        last_score = merge_step(
            &step.left, &step.right, &mut aligned, &weights, scoring, &gap, use_fft,
        );
    }

    // After all merges, pad all sequences to the same width for output.
    // At this point, all sequences should have been involved in at least
    // the final merge, so they should all have the same width. But pad
    // just in case of any edge cases.
    let max_width = aligned.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut aligned {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment { sequences: aligned, names: names.to_vec(), score: last_score }
}

/// Merge two groups by profile alignment.
///
/// Matches C's behavior: only group1 and group2 sequences are modified.
/// "Other" sequences are NOT touched — they keep their current content
/// and width. This prevents spurious gap inflation.
fn merge_step(
    group1: &[usize],
    group2: &[usize],
    aligned: &mut Vec<Vec<u8>>,
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    use_fft: bool,
) -> f64 {
    // Group1 sequences should all have the same width (from their last merge).
    // Group2 sequences should all have the same width (possibly different from group1).
    let width1 = aligned[group1[0]].len();
    let width2 = aligned[group2[0]].len();

    // Build profiles directly from the group sequences.
    // No global gap stripping needed — each group's sequences are already
    // internally consistent from their previous merge.
    let seqs1: Vec<&[u8]> = group1.iter().map(|&i| aligned[i].as_slice()).collect();
    let seqs2: Vec<&[u8]> = group2.iter().map(|&i| aligned[i].as_slice()).collect();

    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();
    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0; group2.len()] };

    let prof1 = Profile::from_aligned(&seqs1, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &w2n, &scoring.amino_map, scoring.nalphabets);

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

    // Build new sequences for group1 and group2 ONLY.
    // "Other" sequences are not touched (matching C's behavior).
    let new_width = aln.operations.len();
    let mut cursor1 = 0usize;
    let mut cursor2 = 0usize;

    // Pre-build new sequences for both groups
    let mut new_seqs_g1: Vec<Vec<u8>> = vec![Vec::with_capacity(new_width); group1.len()];
    let mut new_seqs_g2: Vec<Vec<u8>> = vec![Vec::with_capacity(new_width); group2.len()];

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                for (gi, &idx) in group1.iter().enumerate() {
                    new_seqs_g1[gi].push(
                        if cursor1 < width1 { aligned[idx][cursor1] } else { b'-' }
                    );
                }
                for (gi, &idx) in group2.iter().enumerate() {
                    new_seqs_g2[gi].push(
                        if cursor2 < width2 { aligned[idx][cursor2] } else { b'-' }
                    );
                }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                for (gi, &idx) in group1.iter().enumerate() {
                    new_seqs_g1[gi].push(
                        if cursor1 < width1 { aligned[idx][cursor1] } else { b'-' }
                    );
                }
                for gi in 0..group2.len() {
                    new_seqs_g2[gi].push(b'-');
                }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                for gi in 0..group1.len() {
                    new_seqs_g1[gi].push(b'-');
                }
                for (gi, &idx) in group2.iter().enumerate() {
                    new_seqs_g2[gi].push(
                        if cursor2 < width2 { aligned[idx][cursor2] } else { b'-' }
                    );
                }
                cursor2 += 1;
            }
        }
    }

    // Write back to aligned — only group1 and group2 are modified
    for (gi, &idx) in group1.iter().enumerate() {
        aligned[idx] = new_seqs_g1[gi].clone();
    }
    for (gi, &idx) in group2.iter().enumerate() {
        aligned[idx] = new_seqs_g2[gi].clone();
    }

    aln.score
}

fn group_all_gap_columns(group: &[usize], aligned: &[Vec<u8>], width: usize) -> Vec<bool> {
    let mut all_gap = vec![true; width];
    for &idx in group {
        for (col, &ch) in aligned[idx].iter().enumerate() {
            if ch != b'-' && ch != b'.' {
                all_gap[col] = false;
            }
        }
    }
    all_gap
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_tree::{DistanceMatrix, upgma};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};

    fn check_alignment(result: &MultipleAlignment, original: &[Vec<u8>]) {
        let width = result.width();
        assert!(width > 0);
        for (i, seq) in result.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "seq {i} wrong width: {} vs {width}", seq.len());
            let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
            assert_eq!(ungapped, original[i], "seq {i} residues not preserved");
        }
    }

    #[test]
    fn progressive_two_identical() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![b"ACDEFGHIK".to_vec(), b"ACDEFGHIK".to_vec()];
        let names = vec!["s1".into(), "s2".into()];
        let mut dm = DistanceMatrix::new(2);
        dm.set(0, 1, 0.0);
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        assert_eq!(result.sequences[0], result.sequences[1]);
        check_alignment(&result, &seqs);
    }

    #[test]
    fn progressive_three_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let seqs = vec![
            b"ACDEFGHIK".to_vec(),
            b"ACDEFHIK".to_vec(),
            b"ACDHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];
        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1); dm.set(0, 2, 0.3); dm.set(1, 2, 0.2);
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        check_alignment(&result, &seqs);
    }

    #[test]
    fn progressive_six_sequences_preserves_residues() {
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
        for i in 0..6 { for j in (i+1)..6 { dm.set(i, j, (j-i) as f64 * 0.1); } }
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false, None);
        check_alignment(&result, &seqs);
    }
}
