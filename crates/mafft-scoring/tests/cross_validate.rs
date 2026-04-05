/// Cross-validation tests: compare Rust scoring matrices against C values
/// by calling the C `constants()` function through FFI and comparing.

use std::ptr::addr_of;

use mafft_scoring::*;
use mafft_types::{ScoringModel, SeqType};

#[test]
fn build_context_blosum62_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // Matrix should be 26x26
    assert_eq!(ctx.substitution_matrix.len(), 26);
    assert_eq!(ctx.substitution_matrix[0].len(), 26);

    // Diagonal should be positive for the 20 standard AAs
    for i in 0..20 {
        assert!(
            ctx.substitution_matrix[i][i] > 0,
            "diagonal [{i}][{i}] = {} should be positive",
            ctx.substitution_matrix[i][i]
        );
    }

    // Matrix should be symmetric
    for i in 0..26 {
        for j in 0..26 {
            assert_eq!(
                ctx.substitution_matrix[i][j],
                ctx.substitution_matrix[j][i],
                "asymmetric at [{i}][{j}]"
            );
        }
    }

    // A-A should score higher than A-W (very different amino acids)
    assert!(ctx.substitution_matrix[0][0] > ctx.substitution_matrix[0][17]);
}

#[test]
fn build_context_jtt_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Jtt, SeqType::Protein);
    assert_eq!(ctx.nalphabets, 26);
    assert_eq!(ctx.nscoredalphabets, 20);

    for i in 0..20 {
        assert!(ctx.substitution_matrix[i][i] > 0);
    }
}

#[test]
fn build_context_tm_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Tm, SeqType::Protein);
    for i in 0..20 {
        assert!(ctx.substitution_matrix[i][i] > 0);
    }
}

#[test]
fn build_context_dna_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Dna, SeqType::Dna);

    assert_eq!(ctx.nalphabets, 26);
    assert_eq!(ctx.nscoredalphabets, 10);

    // a-a match should be positive
    assert!(ctx.substitution_matrix[0][0] > 0);
    // a-g transition should be positive (or at least higher than transversion)
    assert!(ctx.substitution_matrix[0][1] > ctx.substitution_matrix[0][2]);
}

#[test]
fn amino_map_covers_standard_residues() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // All standard amino acid letters should map to valid indices
    for &ch in b"ARNDCQEGHILKMFPSTWYVBZX" {
        let idx = ctx.amino_map[ch as usize];
        assert_ne!(idx, 0xFF, "unmapped amino acid: {}", ch as char);
        assert!((idx as usize) < 26);
    }

    // Gap character
    assert_ne!(ctx.amino_map[b'-' as usize], 0xFF);
}

#[test]
fn dna_amino_map_covers_nucleotides() {
    let ctx = build_context(ScoringModel::Dna, SeqType::Dna);

    for &ch in b"agctACGT" {
        let idx = ctx.amino_map[ch as usize];
        assert_ne!(idx, 0xFF, "unmapped nucleotide: {}", ch as char);
    }
}

#[test]
fn polarity_volume_set_for_protein() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // Polarity should be non-zero for at least some amino acid chars
    let pol_sum: f64 = ctx.polarity.iter().map(|x| x.abs()).sum();
    assert!(pol_sum > 0.0, "polarity all zeros");

    let vol_sum: f64 = ctx.volume.iter().map(|x| x.abs()).sum();
    assert!(vol_sum > 0.0, "volume all zeros");
}

#[test]
fn cross_validate_against_c_globals() {
    // Initialize C globals via FFI, then compare key values.
    // We use a dummy empty sequence set since the C code needs nseq/seq args.
    unsafe {
        mafft_sys::initglobalvariables();

        // Set protein mode with BLOSUM62
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(1);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(62);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        // These must be NOTSPECIFIED so constants() applies defaults
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_EX).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_OP).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_dist).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(mafft_sys::NOTSPECIFIED);

        // Create a minimal dummy sequence for constants()
        let seq_data = b"ACGT\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let mut seq_arr: *mut *mut i8 = &mut seq_ptr;

        mafft_sys::constants(1, seq_arr);

        // Read the C penalty value
        let c_penalty = addr_of!(mafft_sys::penalty).read();
        let c_offset = addr_of!(mafft_sys::offset).read();

        // Compare with Rust
        let rust_gap = default_protein_gap_params();

        // These should be very close (exact match depends on rounding)
        assert!(
            (c_penalty - rust_gap.penalty).abs() <= 1,
            "penalty mismatch: C={c_penalty}, Rust={}",
            rust_gap.penalty
        );
        assert!(
            (c_offset - rust_gap.offset).abs() <= 1,
            "offset mismatch: C={c_offset}, Rust={}",
            rust_gap.offset
        );

        mafft_sys::freeconstants();
    }
}
