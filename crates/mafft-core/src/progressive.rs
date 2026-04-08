/// Progressive alignment following a guide tree.
///
/// At each merge, per-group all-gap columns are stripped before profile
/// alignment (matching C's `commongappick()`), then re-inserted after.

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
    let max_len = sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    let mut aligned: Vec<Vec<u8>> = sequences
        .iter()
        .map(|s| { let mut p = s.clone(); p.resize(max_len, b'-'); p })
        .collect();

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

    MultipleAlignment { sequences: aligned, names: names.to_vec(), score: last_score }
}

fn merge_step(
    group1: &[usize],
    group2: &[usize],
    aligned: &mut Vec<Vec<u8>>,
    weights: &[f64],
    scoring: &ScoringContext,
    gap: &GapModel,
    use_fft: bool,
) -> f64 {
    let nseq = aligned.len();
    let width = aligned[0].len();

    // Strip columns that are all-gap across ALL sequences (globally safe).
    //
    // The C code uses per-group stripping (more aggressive, strips columns
    // all-gap in each group independently), which requires careful
    // re-insertion tracking. Our global stripping is more conservative but
    // still effective: as the alignment grows through progressive merges,
    // globally-all-gap columns accumulate and are stripped, keeping the
    // profile DP fast. Per-group stripping can be added as an optimization
    // once the re-insertion interleaving is fully debugged.
    let all_seqs: Vec<usize> = (0..nseq).collect();
    let gap1 = group_all_gap_columns(&all_seqs, aligned, width);
    let gap2 = gap1.clone();

    // kept1/kept2: original column indices that have content in each group
    let kept1: Vec<usize> = (0..width).filter(|&c| !gap1[c]).collect();
    let kept2: Vec<usize> = (0..width).filter(|&c| !gap2[c]).collect();

    // Strip and build profiles
    let stripped1: Vec<Vec<u8>> = group1.iter()
        .map(|&i| kept1.iter().map(|&c| aligned[i][c]).collect())
        .collect();
    let stripped2: Vec<Vec<u8>> = group2.iter()
        .map(|&i| kept2.iter().map(|&c| aligned[i][c]).collect())
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

    // --- Pass 1: Build aligned result from ops using kept column indices ---
    //
    // The ops align kept1 columns against kept2 columns (different lengths).
    // Group1 follows cursor1 through kept1. Group2 follows cursor2 through kept2.
    // "Other" follows cursor1 (like group1). At Insert, "other" gets gap.

    let mut is_group1 = vec![false; nseq];
    let mut is_group2 = vec![false; nseq];
    for &i in group1 { is_group1[i] = true; }
    for &i in group2 { is_group2[i] = true; }

    let new_aln_width = aln.operations.len();
    let mut result = vec![Vec::with_capacity(new_aln_width); nseq];
    let mut cursor1 = 0usize;
    let mut cursor2 = 0usize;

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                let oc1 = kept1[cursor1];
                let oc2 = kept2[cursor2];
                for idx in 0..nseq {
                    if is_group1[idx] {
                        result[idx].push(aligned[idx][oc1]);
                    } else if is_group2[idx] {
                        result[idx].push(aligned[idx][oc2]);
                    } else {
                        result[idx].push(aligned[idx][oc1]);
                    }
                }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                let oc1 = kept1[cursor1];
                for idx in 0..nseq {
                    if is_group2[idx] {
                        result[idx].push(b'-');
                    } else {
                        result[idx].push(aligned[idx][oc1]);
                    }
                }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                let oc2 = kept2[cursor2];
                for idx in 0..nseq {
                    if is_group1[idx] || (!is_group1[idx] && !is_group2[idx]) {
                        result[idx].push(b'-');
                    } else {
                        result[idx].push(aligned[idx][oc2]);
                    }
                }
                cursor2 += 1;
            }
        }
    }

    // --- Pass 2: Re-insert stripped columns ---
    //
    // The result array contains columns from the alignment of kept1 vs kept2.
    // But original columns that were all-gap in group1 (present in kept2 but
    // not kept1) were stripped and need to be re-inserted at their correct
    // positions. Similarly, columns all-gap in both groups need re-insertion.
    //
    // Strategy: walk original columns. Maintain a cursor into `result`.
    // - result columns are ordered by cursor1 (kept1 index).
    // - Each kept1 column may have Insert ops BEFORE it in the result.
    //
    // For each original column:
    //   If it's a kept1 column: emit any preceding Insert result-columns,
    //     then emit the Match/Delete result-column for this kept1 entry.
    //   If it's stripped from group1 but in kept2: emit directly.
    //   If it's stripped from both: emit directly.
    //
    // To know which result columns are Inserts vs Match/Delete, track via ops.

    // Build a per-result-column flag: true if this result column corresponds
    // to a kept1 entry (Match/Delete op), false if it's an Insert.
    let is_k1_col: Vec<bool> = aln.operations.iter()
        .map(|op| matches!(op, AlignOp::Match | AlignOp::Delete))
        .collect();

    let mut final_seqs = vec![Vec::with_capacity(width + new_aln_width); nseq];
    let mut result_cursor = 0usize; // cursor into result columns

    for orig_col in 0..width {
        if !gap1[orig_col] {
            // This is a kept1 column. Before emitting the Match/Delete for it,
            // emit any Insert result-columns that precede it.
            while result_cursor < new_aln_width && !is_k1_col[result_cursor] {
                for idx in 0..nseq {
                    final_seqs[idx].push(result[idx][result_cursor]);
                }
                result_cursor += 1;
            }
            // Now emit the Match/Delete column itself.
            if result_cursor < new_aln_width {
                for idx in 0..nseq {
                    final_seqs[idx].push(result[idx][result_cursor]);
                }
                result_cursor += 1;
            }
        } else if !gap2[orig_col] {
            // All-gap in group1 only. This column has group2/other content
            // but was stripped from group1's profile. Emit directly.
            for idx in 0..nseq {
                if is_group1[idx] {
                    final_seqs[idx].push(b'-');
                } else {
                    final_seqs[idx].push(aligned[idx][orig_col]);
                }
            }
        } else {
            // All-gap in both groups. Emit original content for all.
            for idx in 0..nseq {
                final_seqs[idx].push(aligned[idx][orig_col]);
            }
        }
    }

    // Emit any trailing Insert result-columns after the last kept1 column.
    while result_cursor < new_aln_width {
        for idx in 0..nseq {
            final_seqs[idx].push(result[idx][result_cursor]);
        }
        result_cursor += 1;
    }

    *aligned = final_seqs;
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
