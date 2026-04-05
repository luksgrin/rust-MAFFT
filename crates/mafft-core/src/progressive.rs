/// Progressive alignment following a guide tree.
///
/// Walks the tree bottom-up. At each internal node, the two child clusters
/// are aligned, and the result is applied to ALL sequences.
///
/// Key design (matching C code):
/// - All sequences share the same alignment width at all times.
/// - At each merge, common gap columns are stripped before profile alignment.
/// - After alignment, gap columns are re-inserted and all sequences expanded.

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
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    for step in &topology.steps {
        last_score = merge_step(
            &step.left, &step.right, &mut aligned, &weights, scoring, &gap, use_fft,
        );
    }

    MultipleAlignment { sequences: aligned, names: names.to_vec(), score: last_score }
}

/// Perform one merge step: strip common gaps, align profiles, expand all.
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

    // Step 1: Find common-gap columns (all-gap in BOTH groups).
    // These carry no information and waste DP time.
    let is_common_gap = find_common_gap_columns(group1, group2, aligned, width);
    let kept_cols: Vec<usize> = (0..width).filter(|&c| !is_common_gap[c]).collect();
    let num_common_gaps = width - kept_cols.len();

    // Step 2: Strip common-gap columns for profile construction.
    let stripped1: Vec<Vec<u8>> = group1.iter()
        .map(|&i| kept_cols.iter().map(|&c| aligned[i][c]).collect())
        .collect();
    let stripped2: Vec<Vec<u8>> = group2.iter()
        .map(|&i| kept_cols.iter().map(|&c| aligned[i][c]).collect())
        .collect();

    let s1_refs: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
    let s2_refs: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();

    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();
    let sum1: f64 = w1.iter().sum();
    let sum2: f64 = w2.iter().sum();
    let w1n: Vec<f64> = if sum1 > 0.0 { w1.iter().map(|w| w / sum1).collect() } else { vec![1.0; group1.len()] };
    let w2n: Vec<f64> = if sum2 > 0.0 { w2.iter().map(|w| w / sum2).collect() } else { vec![1.0; group2.len()] };

    // Step 3: Build profiles from stripped sequences and align.
    let prof1 = Profile::from_aligned(&s1_refs, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&s2_refs, &w2n, &scoring.amino_map, scoring.nalphabets);

    let stripped_width = kept_cols.len();

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

    // Step 4: Build new sequences from alignment ops.
    //
    // The ops align the stripped columns. Both stripped profiles have the
    // same length (stripped_width). cursor1 and cursor2 each walk through
    // kept_cols indices 0..stripped_width-1.
    //
    // For each op:
    // - Match: cursor1 advances (maps to kept_cols[cursor1] in original),
    //          cursor2 advances (maps to kept_cols[cursor2] in original).
    //   Group1 takes from original[kept_cols[cursor1]].
    //   Group2 takes from original[kept_cols[cursor2]].
    //   Others follow cursor1.
    //
    // - Delete: cursor1 advances. Group1 takes content, group2 gets gap.
    //   Others follow cursor1.
    //
    // - Insert: cursor2 advances. Group2 takes content, group1 gets gap.
    //   Others get gap (cursor2's column will be covered by cursor1 later).
    //
    // After processing all ops, re-insert common-gap columns at their
    // original relative positions.

    let mut is_group1 = vec![false; nseq];
    let mut is_group2 = vec![false; nseq];
    for &i in group1 { is_group1[i] = true; }
    for &i in group2 { is_group2[i] = true; }

    // Build the aligned-but-no-common-gaps result
    let new_stripped_width = aln.operations.len();
    let mut new_stripped = vec![Vec::with_capacity(new_stripped_width); nseq];

    let mut cursor1 = 0usize;
    let mut cursor2 = 0usize;

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                let old_col1 = kept_cols[cursor1];
                let old_col2 = kept_cols[cursor2];
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_stripped[idx].push(aligned[idx][old_col1]);
                    } else if is_group2[idx] {
                        new_stripped[idx].push(aligned[idx][old_col2]);
                    } else {
                        new_stripped[idx].push(aligned[idx][old_col1]);
                    }
                }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                let old_col1 = kept_cols[cursor1];
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_stripped[idx].push(aligned[idx][old_col1]);
                    } else if is_group2[idx] {
                        new_stripped[idx].push(b'-');
                    } else {
                        new_stripped[idx].push(aligned[idx][old_col1]);
                    }
                }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                let old_col2 = kept_cols[cursor2];
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_stripped[idx].push(b'-');
                    } else if is_group2[idx] {
                        new_stripped[idx].push(aligned[idx][old_col2]);
                    } else {
                        new_stripped[idx].push(b'-');
                    }
                }
                cursor2 += 1;
            }
        }
    }

    // Step 5: Re-insert common-gap columns at their original positions.
    //
    // We need to interleave the new_stripped result with the common-gap
    // columns. The common-gap columns were at specific positions in the
    // original alignment. We place them back relative to the kept columns.
    //
    // Build the final sequences by walking through original column indices.
    // For each original column index c:
    // - If c was a common gap: insert a gap column for all sequences
    // - If c was a kept column: take the next column from new_stripped
    //
    // But the new_stripped may be WIDER than the number of kept columns
    // (because the alignment introduced new columns via Delete/Insert ops).
    // So we can't do a 1:1 mapping with original positions.
    //
    // Instead, we interleave by tracking the original-column order:
    // - Walk through original columns 0..width-1
    // - For each kept column, emit the corresponding new_stripped columns
    //   (there may be multiple if the alignment expanded at that point)
    // - For each common-gap column, emit a gap column
    //
    // To do this properly, we need to know which new_stripped columns
    // correspond to which kept_cols entries. The mapping is:
    // - new_stripped columns corresponding to cursor1 positions map to
    //   kept_cols[0], kept_cols[1], ... (Match/Delete ops advance cursor1)
    // - Insert ops don't correspond to any kept_cols entry (they're new)
    //
    // So we build a list of (kept_col_index_or_none, new_stripped_col):
    let mut col_origins: Vec<Option<usize>> = Vec::with_capacity(new_stripped_width);
    {
        let mut c1 = 0usize;
        let mut c2 = 0usize;
        for op in &aln.operations {
            match op {
                AlignOp::Match => {
                    col_origins.push(Some(c1)); // maps to kept_cols[c1]
                    c1 += 1;
                    c2 += 1;
                }
                AlignOp::Delete => {
                    col_origins.push(Some(c1)); // maps to kept_cols[c1]
                    c1 += 1;
                }
                AlignOp::Insert => {
                    col_origins.push(None); // new column, no original position
                    c2 += 1;
                }
            }
        }
    }

    // Now build final sequences by interleaving:
    let final_width = new_stripped_width + num_common_gaps;
    let mut final_aligned = vec![Vec::with_capacity(final_width); nseq];

    let mut stripped_cursor = 0usize; // cursor into new_stripped columns
    let mut kept_cursor = 0usize; // which kept_col we expect next

    for orig_col in 0..width {
        if is_common_gap[orig_col] {
            // Re-insert common-gap column: all sequences get gap
            for idx in 0..nseq {
                final_aligned[idx].push(b'-');
            }
        } else {
            // This original column corresponds to kept_cols[kept_cursor].
            // Emit all new_stripped columns that map to this kept position,
            // plus any Insert columns (None) that precede the next kept position.
            while stripped_cursor < new_stripped_width {
                match col_origins[stripped_cursor] {
                    Some(k) if k == kept_cursor => {
                        // This new column maps to the current kept position
                        for idx in 0..nseq {
                            final_aligned[idx].push(new_stripped[idx][stripped_cursor]);
                        }
                        stripped_cursor += 1;
                        break; // move to next original column
                    }
                    None => {
                        // Insert column (new, no original position) — emit it
                        for idx in 0..nseq {
                            final_aligned[idx].push(new_stripped[idx][stripped_cursor]);
                        }
                        stripped_cursor += 1;
                    }
                    Some(k) => {
                        // Maps to a different kept position — we've moved past
                        // the current one. This shouldn't happen if cursor1
                        // advances monotonically, but handle gracefully.
                        break;
                    }
                }
            }
            kept_cursor += 1;
        }
    }

    // Emit any remaining new_stripped columns (Insert ops after the last kept column)
    while stripped_cursor < new_stripped_width {
        for idx in 0..nseq {
            final_aligned[idx].push(new_stripped[idx][stripped_cursor]);
        }
        stripped_cursor += 1;
    }

    *aligned = final_aligned;
    aln.score
}

/// Find columns that are all-gap in ALL sequences (true common gap columns).
///
/// Only columns where every sequence has a gap can be safely stripped,
/// since "other" sequences (not in either group) may have real content
/// at columns that are all-gap within the two merge groups.
fn find_common_gap_columns(
    _group1: &[usize],
    _group2: &[usize],
    aligned: &[Vec<u8>],
    width: usize,
) -> Vec<bool> {
    let mut all_gap = vec![true; width];

    for seq in aligned.iter() {
        for (col, &ch) in seq.iter().enumerate() {
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
        assert!(width > 0, "alignment width should be > 0");
        for (i, seq) in result.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "seq {i} has wrong width: {} vs {width}", seq.len());
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
        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
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
        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
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
        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
        check_alignment(&result, &seqs);
    }
}
