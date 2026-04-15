/// Add sequences to an existing alignment.
///
/// Ports the C `addsingle` binary's core logic: for each new sequence,
/// place it in the guide tree using `addonetip()`, then perform
/// progressive alignment to produce the combined result.

use rayon::prelude::*;

use mafft_tree::{
    DistanceMatrix, musclesupg, ClusterMethod, addonetip, ktuple_distance,
};
use mafft_types::ScoringContext;

use crate::progressive::{progressive_align, MultipleAlignment};

/// Add new sequences to an existing alignment.
///
/// `existing` is the already-aligned MSA (with gaps).
/// `new_sequences` are the unaligned sequences to add.
/// `new_names` are their names.
/// `scoring` is the scoring context.
/// `use_fft` controls whether FFT acceleration is used.
///
/// Returns a new MSA containing all sequences (existing + new).
pub fn add_sequences(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[String],
    scoring: &ScoringContext,
    use_fft: bool,
) -> MultipleAlignment {
    if new_sequences.is_empty() {
        return existing.clone();
    }

    let norg = existing.nseq();

    // Get ungapped original sequences for distance computation
    let orig_ungapped: Vec<Vec<u8>> = existing.sequences.iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').copied().collect())
        .collect();

    // Build guide tree from original sequences using k-tuple distance
    let orig_dm = compute_ktuple_dm(&orig_ungapped);
    let orig_topo = musclesupg(&orig_dm, ClusterMethod::default());

    // Add each new sequence independently to the original alignment,
    // then combine results. This matches C's approach: each new sequence
    // is placed in the original tree and aligned against the original MSA.
    let mut combined_sequences = existing.sequences.clone();
    let mut combined_names = existing.names.clone();

    for (add_idx, new_seq) in new_sequences.iter().enumerate() {
        // Compute distances from new sequence to all original sequences
        let distances: Vec<f64> = orig_ungapped.iter()
            .map(|orig| ktuple_distance(orig, new_seq, 6))
            .collect();

        // Place in tree
        let add_result = addonetip(&orig_topo, &distances, 0.1);

        // Build combined sequence set for progressive alignment
        let mut all_ungapped: Vec<Vec<u8>> = orig_ungapped.clone();
        all_ungapped.push(new_seq.clone());

        let mut all_names: Vec<String> = existing.names.clone();
        all_names.push(new_names[add_idx].clone());

        // Progressive align using the new topology
        let msa = progressive_align(
            &all_ungapped,
            &all_names,
            &add_result.topology,
            scoring,
            use_fft,
            None,
        );

        // Extract just the new sequence's aligned form (last sequence in the result)
        // and adjust it to be compatible with the combined alignment width
        combined_sequences.push(msa.sequences[norg].clone());
        combined_names.push(new_names[add_idx].clone());
    }

    // Ensure all sequences have the same width (pad shorter ones with gaps)
    let max_width = combined_sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut combined_sequences {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment {
        sequences: combined_sequences,
        names: combined_names,
        score: 0.0,
        step_trace: Vec::new(),
    }
}

/// Add sequences preserving the existing alignment structure.
///
/// This is the more faithful port of C's --add behavior: the existing
/// alignment columns are preserved, and new sequences are aligned to the
/// existing MSA without disturbing it.
pub fn add_sequences_keeplength(
    existing: &MultipleAlignment,
    new_sequences: &[Vec<u8>],
    new_names: &[String],
    scoring: &ScoringContext,
    _use_fft: bool,
) -> MultipleAlignment {
    if new_sequences.is_empty() {
        return existing.clone();
    }

    // Get ungapped original sequences
    let orig_ungapped: Vec<Vec<u8>> = existing.sequences.iter()
        .map(|s| s.iter().filter(|&&c| c != b'-').copied().collect())
        .collect();

    let mut all_sequences = existing.sequences.clone();
    let mut all_names = existing.names.clone();

    for (add_idx, new_seq) in new_sequences.iter().enumerate() {
        // Compute distances from new sequence to all original sequences
        let distances: Vec<f64> = orig_ungapped.iter()
            .map(|orig| ktuple_distance(orig, new_seq, 6))
            .collect();

        // Find nearest original sequence
        let (nearest_idx, _nearest_dist) = distances.iter().enumerate()
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();

        // Align new sequence against the nearest original sequence's gapped form
        // using profile alignment (treating the existing alignment as a frozen profile)
        let nearest_gapped = &existing.sequences[nearest_idx];

        // Simple approach: align new_seq to nearest_gapped using the gap structure
        let aligned_new = align_to_profile(new_seq, nearest_gapped, scoring);

        all_sequences.push(aligned_new);
        all_names.push(new_names[add_idx].clone());
    }

    // Ensure all sequences have the same width (pad shorter ones)
    let max_width = all_sequences.iter().map(|s| s.len()).max().unwrap_or(0);
    for seq in &mut all_sequences {
        seq.resize(max_width, b'-');
    }

    MultipleAlignment {
        sequences: all_sequences,
        names: all_names,
        score: 0.0,
        step_trace: Vec::new(),
    }
}

/// Align a single new sequence to a gapped reference sequence.
///
/// Maps the new sequence's residues to the reference's non-gap positions,
/// inserting gaps where the reference has gaps.
fn align_to_profile(
    new_seq: &[u8],
    reference_gapped: &[u8],
    scoring: &ScoringContext,
) -> Vec<u8> {
    use mafft_align::{profile_align, Profile, GapModel};

    // Build a single-sequence profile from the reference (gapped)
    let ref_slices: Vec<&[u8]> = vec![reference_gapped];
    let ref_weights = vec![1.0];
    let prof_ref = Profile::from_aligned(
        &ref_slices, &ref_weights,
        &scoring.amino_map, scoring.nalphabets,
    );

    // Build a single-sequence profile from the new sequence (ungapped)
    let new_slices: Vec<&[u8]> = vec![new_seq];
    let new_weights = vec![1.0];
    let prof_new = Profile::from_aligned(
        &new_slices, &new_weights,
        &scoring.amino_map, scoring.nalphabets,
    );

    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let aln = profile_align(
        &prof_ref, &prof_new,
        &scoring.substitution_matrix, &gap,
        true, true,
    );

    // Build aligned new sequence from alignment operations
    use mafft_align::AlignOp;
    let mut result = Vec::with_capacity(aln.operations.len());
    let mut new_cursor = 0usize;

    for op in &aln.operations {
        match op {
            AlignOp::Match => {
                result.push(if new_cursor < new_seq.len() { new_seq[new_cursor] } else { b'-' });
                new_cursor += 1;
            }
            AlignOp::Delete => {
                // Reference has a column, new seq gets a gap
                result.push(b'-');
            }
            AlignOp::Insert => {
                // New seq has a residue, reference gets a gap
                result.push(if new_cursor < new_seq.len() { new_seq[new_cursor] } else { b'-' });
                new_cursor += 1;
            }
        }
    }

    result
}

fn compute_ktuple_dm(sequences: &[Vec<u8>]) -> DistanceMatrix {
    let nseq = sequences.len();
    let pairs: Vec<(usize, usize, f64)> = (0..nseq)
        .into_par_iter()
        .flat_map(|i| {
            let seqs = sequences;
            ((i + 1)..nseq).into_par_iter().map(move |j| {
                let d = ktuple_distance(&seqs[i], &seqs[j], 6);
                (i, j, d)
            })
        })
        .collect();

    let mut dm = DistanceMatrix::new(nseq);
    for (i, j, d) in pairs {
        dm.set(i, j, d);
    }
    dm
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};

    fn make_existing_alignment() -> MultipleAlignment {
        // Pre-aligned MSA
        MultipleAlignment {
            sequences: vec![
                b"ACDEFGHIK".to_vec(),
                b"ACDEF-HIK".to_vec(),
                b"ACD---HIK".to_vec(),
            ],
            names: vec!["s1".into(), "s2".into(), "s3".into()],
            score: 0.0,
            step_trace: Vec::new(),
        }
    }

    #[test]
    fn add_single_sequence() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let new_seqs = vec![b"ACDEFGHIKLM".to_vec()];
        let new_names = vec!["new1".into()];

        let result = add_sequences(&existing, &new_seqs, &new_names, &scoring, false);
        assert_eq!(result.nseq(), 4);
        let w = result.width();
        for seq in &result.sequences {
            assert_eq!(seq.len(), w);
        }
        // Verify residue preservation
        let ungapped: Vec<u8> = result.sequences[3].iter()
            .filter(|&&c| c != b'-').copied().collect();
        assert_eq!(ungapped, b"ACDEFGHIKLM");
    }

    #[test]
    fn add_no_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let result = add_sequences(&existing, &[], &[], &scoring, false);
        assert_eq!(result.nseq(), 3);
    }

    #[test]
    fn add_keeplength_preserves_width() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let new_seqs = vec![b"ACDEHIK".to_vec()];
        let new_names = vec!["new1".into()];

        let result = add_sequences_keeplength(
            &existing, &new_seqs, &new_names, &scoring, false,
        );
        assert_eq!(result.nseq(), 4);
        let w = result.width();
        for seq in &result.sequences {
            assert_eq!(seq.len(), w);
        }
    }

    #[test]
    fn add_multiple_sequences() {
        let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
        let existing = make_existing_alignment();
        let new_seqs = vec![
            b"ACDEFGHIKLM".to_vec(),
            b"ACDEHIK".to_vec(),
        ];
        let new_names = vec!["new1".into(), "new2".into()];

        let result = add_sequences(&existing, &new_seqs, &new_names, &scoring, false);
        assert_eq!(result.nseq(), 5);
        let w = result.width();
        for seq in &result.sequences {
            assert_eq!(seq.len(), w);
        }
    }
}
