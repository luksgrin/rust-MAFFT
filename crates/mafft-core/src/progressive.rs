/// Progressive alignment following a guide tree.
///
/// Walks the tree bottom-up. At each internal node, the two child clusters
/// are aligned, and the result is applied to ALL sequences.
///
/// Key design (matching C code):
/// - All sequences share the same alignment width at all times.
/// - At each merge, profiles are built from the full-width sequences.
/// - The profile DP aligns column-by-column (both profiles have width W).
/// - After alignment, ALL sequences are expanded to the new width by
///   inserting gap columns at the appropriate positions.

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

/// Perform one merge step: align two groups and expand all sequences.
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

    // Build profiles from full-width aligned sequences.
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

    // Both profiles have the same length (= width, the current alignment width).
    // The ops align these W columns against each other.
    //
    // Now we need to expand ALL sequences to the new width.
    //
    // The ops tell us:
    // - Match at (i,j): old column i (from group1 cursor) is aligned with
    //   old column j (from group2 cursor). New alignment has ONE column for this.
    // - Delete: old column i (group1) appears alone. Group2 gets gap.
    // - Insert: old column j (group2) appears alone. Group1 gets gap.
    //
    // For "other" sequences: they have content at ALL old columns. At a Match
    // position, old columns i and j are being merged into one. If i == j,
    // "other" sequence contributes its column i content. If i != j, we must
    // choose — we use column i (from cursor1) since that's the "main" timeline.
    // Column j's content from "other" will be lost... UNLESS it appears
    // separately in a Delete/Insert op.
    //
    // This is correct because: the DP aligns ALL W columns from group1 with
    // ALL W columns from group2. Every old column appears exactly once:
    // either consumed by cursor1 (at Match/Delete) or by cursor2 (at
    // Match/Insert). For cursor1, columns appear in order 0,1,2,...,W-1.
    // Same for cursor2. So every old column index is visited exactly once
    // by each cursor. At a Match op, cursor1's column and cursor2's column
    // are merged. For "other" sequences, we include cursor1's column at
    // Match/Delete, and cursor2's column at Insert. This covers every old
    // column exactly once (cursor1 consumes W columns across Match+Delete,
    // cursor2 consumes W columns across Match+Insert, but Match overlaps
    // both, so total distinct old columns = W).
    //
    // Wait, that means "other" sequences contribute W columns total:
    // (Match+Delete columns from cursor1) + (Insert columns from cursor2)
    // = W + (W - Match_count) = 2W - Match_count = new width.
    // But "other" sequences only HAVE W columns of content. We'd be pulling
    // W columns from cursor1 and (W-Match) from cursor2 = 2W-Match total.
    // That's more than W columns pulled from a sequence of length W!
    //
    // The issue: cursor1 and cursor2 both index the SAME W columns (0..W-1)
    // in the same order. At a Match op with cursor1=3 and cursor2=5,
    // "other" sequence would use column 3 here. At a later Insert with
    // cursor2=6 (having skipped 4,5 because they were used in Delete ops),
    // wait no — cursor2 advances monotonically through 0..W-1 just like
    // cursor1. So at Insert, "other" uses cursor2's column, which is a
    // different column from what cursor1 consumed.
    //
    // So total columns consumed from "other":
    // - At Match+Delete: cursor1 values (= indices 0..W-1, all W columns)
    // - At Insert: cursor2 values that aren't already at Match positions
    //
    // But cursor2 also goes 0..W-1. At Match ops, cursor2's value is used
    // for group2 but "other" uses cursor1's value instead. So "other"
    // doesn't use cursor2 at Match ops. At Insert ops, "other" uses cursor2.
    // Insert count = W - Match count.
    // So "other" pulls W columns from cursor1 + (W - Match) from cursor2.
    // That's 2W - Match total. But "other" only has W columns!
    //
    // This means we'd be reading some columns TWICE (once from cursor1 at
    // Match/Delete, and again from cursor2 at Insert). That causes
    // duplication of residues for "other" sequences.
    //
    // THE FIX: "other" sequences should use a SINGLE cursor that advances
    // through their own columns 0..W-1 in order. At each op (Match, Delete,
    // or Insert), advance this cursor. The "other" sequence contributes its
    // current column at EVERY op position. This means: at new positions
    // where group1 has a gap (Insert), "other" still contributes content.
    // At new positions where group2 has a gap (Delete), "other" still
    // contributes content. The other sequence's column order is preserved,
    // and exactly W columns are consumed (one per op that advances the
    // "other" cursor)... but wait, there are MORE than W ops (new width > W).
    //
    // Ugh. The new width is W + Delete_count + Insert_count - Match_count?
    // No: consumed_from_1 = Match + Delete = W. consumed_from_2 = Match + Insert = W.
    // New width = Match + Delete + Insert = W - Match + W = 2W - Match.
    //
    // So new width > W (unless Match = W, meaning identity alignment).
    // "Other" sequences have only W columns. They need 2W-Match new columns.
    // So (W-Match) new columns must be gaps for "other" sequences.
    //
    // Which ops should be gaps for "other"? The ops that bring in "new"
    // columns that aren't part of the shared column space. Since both
    // cursor1 and cursor2 go through columns 0..W-1, the "new" columns
    // are when cursor1 and cursor2 are at different positions and both
    // consuming simultaneously. But that only happens at Match ops.
    //
    // OK, I think the correct approach is: "other" sequences use cursor1
    // at Match/Delete, and get gaps at Insert. This consumes exactly W
    // columns (= Match + Delete). At Insert positions (new width - W
    // positions), "other" gets gaps. Total = W real + (new_width - W) gaps
    // = new_width. ✓
    //
    // But then "other" sequences' columns that correspond to cursor2
    // positions (which may have real content) are IGNORED at Insert positions.
    // Those columns get "used up" by cursor2 at Match and Insert ops, but
    // "other" doesn't see them. This means "other" sequences lose content
    // at columns that are in cursor2 but not cursor1.
    //
    // EXCEPT: cursor1 ALSO goes through 0..W-1 in order. So cursor1
    // visits every column 0..W-1. Therefore "other" sequences (following
    // cursor1) also visit every column 0..W-1. No content is lost.
    //
    // At Insert positions, cursor2 visits a column that cursor1 ALSO
    // visits (later, at a Delete or Match op). So "other"'s content from
    // that column will appear at the cursor1 position, not the Insert
    // position. The Insert position gets a gap for "other", and the actual
    // content appears elsewhere. This is correct!

    let mut new_aligned = vec![Vec::with_capacity(aln.operations.len()); nseq];

    let mut cursor1 = 0usize; // group1's column cursor through 0..W-1
    let mut cursor2 = 0usize; // group2's column cursor through 0..W-1

    // Pre-compute membership for speed
    let mut is_group1 = vec![false; nseq];
    let mut is_group2 = vec![false; nseq];
    for &i in group1 { is_group1[i] = true; }
    for &i in group2 { is_group2[i] = true; }

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_aligned[idx].push(aligned[idx][cursor1]);
                    } else if is_group2[idx] {
                        new_aligned[idx].push(aligned[idx][cursor2]);
                    } else {
                        // "Other": follow cursor1 (covers all W old columns)
                        new_aligned[idx].push(aligned[idx][cursor1]);
                    }
                }
                cursor1 += 1;
                cursor2 += 1;
            }
            AlignOp::Delete => {
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_aligned[idx].push(aligned[idx][cursor1]);
                    } else if is_group2[idx] {
                        new_aligned[idx].push(b'-');
                    } else {
                        new_aligned[idx].push(aligned[idx][cursor1]);
                    }
                }
                cursor1 += 1;
            }
            AlignOp::Insert => {
                for idx in 0..nseq {
                    if is_group1[idx] {
                        new_aligned[idx].push(b'-');
                    } else if is_group2[idx] {
                        new_aligned[idx].push(aligned[idx][cursor2]);
                    } else {
                        // Gap for "other" — cursor2's column will be
                        // covered when cursor1 reaches it later.
                        new_aligned[idx].push(b'-');
                    }
                }
                cursor2 += 1;
            }
        }
    }

    *aligned = new_aligned;
    aln.score
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
            b"ACDEFHIK".to_vec(),
            b"ACDHIK".to_vec(),
        ];
        let names = vec!["s1".into(), "s2".into(), "s3".into()];
        let mut dm = DistanceMatrix::new(3);
        dm.set(0, 1, 0.1); dm.set(0, 2, 0.3); dm.set(1, 2, 0.2);
        let topo = upgma(&dm);
        let result = progressive_align(&seqs, &names, &topo, &scoring, false);
        assert_eq!(result.nseq(), 3);
        let width = result.width();
        for seq in &result.sequences { assert_eq!(seq.len(), width); }
        for (i, seq) in result.sequences.iter().enumerate() {
            let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
            assert_eq!(ungapped, seqs[i], "seq {i} residues not preserved");
        }
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
        let width = result.width();
        for (i, seq) in result.sequences.iter().enumerate() {
            assert_eq!(seq.len(), width, "seq {i} wrong width");
            let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
            assert_eq!(ungapped, seqs[i], "seq {i} residues not preserved");
        }
    }
}
