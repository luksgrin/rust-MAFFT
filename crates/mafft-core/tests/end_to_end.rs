/// End-to-end integration tests: read real test data, align, verify output.

use std::path::PathBuf;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::read_fasta;

fn test_data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../mafft-upstream/test")
        .join(name)
}

#[test]
fn align_sample_fasta() {
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let input = read_fasta(test_data_path("sample")).unwrap();
    assert_eq!(input.nseq(), 36);

    let msa = engine.align(&input);

    // All sequences should have the same width
    let width = msa.width();
    assert!(width > 0, "alignment should have non-zero width");
    for (i, seq) in msa.sequences.iter().enumerate() {
        assert_eq!(
            seq.len(),
            width,
            "sequence {} has length {} but expected {}",
            msa.names[i],
            seq.len(),
            width
        );
    }

    // Removing gaps should recover original sequences
    for (i, seq) in msa.sequences.iter().enumerate() {
        let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(
            ungapped, input.sequences[i].data,
            "ungapped sequence {} doesn't match original",
            msa.names[i]
        );
    }
}

#[test]
fn align_sample_with_refinement() {
    // Use only the first 6 sequences to keep refinement fast
    let full_input = read_fasta(test_data_path("sample")).unwrap();
    let input = mafft_types::SequenceSet {
        sequences: full_input.sequences[..6].to_vec(),
        seq_type: full_input.seq_type,
    };

    // First verify progressive-only works
    let engine_prog = MafftEngine::new(AlignmentMode::FftNs2);
    let msa_prog = engine_prog.align(&input);
    for (i, seq) in msa_prog.sequences.iter().enumerate() {
        let residues = seq.iter().filter(|&&c| c != b'-').count();
        assert_eq!(residues, input.sequences[i].data.len(),
            "progressive lost residues for seq {i}: {} vs {}", residues, input.sequences[i].data.len());
    }

    // Now test with refinement
    let engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 });
    let msa = engine.align(&input);

    let width = msa.width();
    assert!(width > 0);
    for (i, seq) in msa.sequences.iter().enumerate() {
        assert_eq!(seq.len(), width, "sequence {i} width mismatch after refinement");
    }

    // Ungapped residue counts should be preserved
    for (i, seq) in msa.sequences.iter().enumerate() {
        let residue_count = seq.iter().filter(|&&c| c != b'-').count();
        assert_eq!(residue_count, input.sequences[i].data.len(),
            "sequence {i} lost residues during refinement: {} vs {}",
            residue_count, input.sequences[i].data.len());
    }
}

#[test]
fn compare_against_c_reference() {
    // Read the C reference alignment (test/sample.fftns2)
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();

    // Run Rust engine on same input
    let input = read_fasta(test_data_path("sample")).unwrap();
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let msa = engine.align(&input);

    // 1. Same number of sequences
    assert_eq!(msa.nseq(), c_ref.nseq(), "different number of sequences");

    // 2. Same residue content per sequence (ungapped)
    for i in 0..msa.nseq() {
        let rust_ungapped: Vec<u8> = msa.sequences[i].iter().filter(|&&c| c != b'-').cloned().collect();
        let c_ungapped: Vec<u8> = c_ref.sequences[i].data.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(
            rust_ungapped, c_ungapped,
            "sequence {i} has different residue content between Rust and C"
        );
    }

    // 3. Compare alignment quality via sum-of-pairs identity
    let rust_sp = sum_of_pairs_identity(&msa.sequences);
    let c_seqs: Vec<Vec<u8>> = c_ref.sequences.iter().map(|s| s.data.clone()).collect();
    let c_sp = sum_of_pairs_identity(&c_seqs);

    // Report quality comparison (not a hard failure — different algorithms
    // produce different alignments, but quality should be in the same ballpark)
    let ratio = if c_sp > 0.0 { rust_sp / c_sp } else { 1.0 };
    eprintln!(
        "Alignment quality: Rust SP={rust_sp:.4}, C SP={c_sp:.4}, ratio={ratio:.4}"
    );
    // Rust alignment should be at least 50% as good as C's
    // (a loose bound — we're not matching C's exact algorithm)
    assert!(
        ratio > 0.5,
        "Rust alignment quality too low: {rust_sp:.4} vs C's {c_sp:.4} (ratio {ratio:.4})"
    );
}

/// Compute sum-of-pairs identity score for an alignment.
fn sum_of_pairs_identity(sequences: &[Vec<u8>]) -> f64 {
    let n = sequences.len();
    if n < 2 { return 0.0; }
    let mut total_match = 0u64;
    let mut total_aligned = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let len = sequences[i].len().min(sequences[j].len());
            for k in 0..len {
                let a = sequences[i][k];
                let b = sequences[j][k];
                if a != b'-' && b != b'-' {
                    total_aligned += 1;
                    if a == b {
                        total_match += 1;
                    }
                }
            }
        }
    }
    if total_aligned == 0 { 0.0 } else { total_match as f64 / total_aligned as f64 }
}

#[test]
fn align_rna_sample() {
    let engine = MafftEngine::default();
    let input = read_fasta(test_data_path("samplerna")).unwrap();
    assert!(input.seq_type.is_nucleotide());

    let msa = engine.align(&input);
    let width = msa.width();
    assert!(width > 0);
    for seq in &msa.sequences {
        assert_eq!(seq.len(), width);
    }
}
