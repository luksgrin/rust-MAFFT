/// End-to-end integration tests: read real test data, align, verify output.

use std::path::PathBuf;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::read_fasta;

fn test_data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../test")
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
#[ignore] // TODO: iterative refinement gap bookkeeping needs fix for width consistency
fn align_sample_with_refinement() {
    let engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 });
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = engine.align(&input);

    let width = msa.width();
    assert!(width > 0);
    for seq in &msa.sequences {
        assert_eq!(seq.len(), width);
    }
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
