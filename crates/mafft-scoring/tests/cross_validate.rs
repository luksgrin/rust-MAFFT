/// Cross-validation tests: compare Rust scoring matrices against C values
/// cell-by-cell by calling the C `constants()` function through FFI.
///
/// C globals are shared mutable state, so all FFI tests must be serialized.

use std::ptr::addr_of;
use std::sync::Mutex;

use mafft_scoring::*;
use mafft_types::{ScoringModel, SeqType};

/// Mutex to serialize all tests that touch C globals.
static C_MUTEX: Mutex<()> = Mutex::new(());

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Initialize C globals to NOTSPECIFIED so constants() applies defaults.
unsafe fn init_c_globals() {
    mafft_sys::initglobalvariables();
    std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::ppenalty_EX).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::ppenalty_OP).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::ppenalty_dist).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::poffset).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
    std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
}

/// Call C constants() with a dummy sequence.
unsafe fn call_c_constants(dorp: u8, scoremtx: i32, nblosum: i32) {
    std::ptr::addr_of_mut!(mafft_sys::dorp).write(dorp as i32);
    std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(scoremtx);
    std::ptr::addr_of_mut!(mafft_sys::nblosum).write(nblosum);
    std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);

    let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
    let mut seq_ptr = seq_data.as_ptr() as *mut i8;
    let seq_arr: *mut *mut i8 = &mut seq_ptr;
    mafft_sys::constants(1, seq_arr);
}

/// Read the C n_dis matrix into a Vec<Vec<i32>>.
unsafe fn read_c_n_dis() -> Vec<Vec<i32>> {
    let nalpha = addr_of!(mafft_sys::nalphabets).read() as usize;
    let n_dis_ptr = addr_of!(mafft_sys::n_dis).read();

    let mut matrix = vec![vec![0i32; nalpha]; nalpha];
    for i in 0..nalpha {
        let row_ptr = *n_dis_ptr.add(i);
        for j in 0..nalpha {
            matrix[i][j] = *row_ptr.add(j);
        }
    }
    matrix
}

/// Read the C n_disFFT matrix into a Vec<Vec<i32>>.
unsafe fn read_c_n_dis_fft() -> Vec<Vec<i32>> {
    let nalpha = addr_of!(mafft_sys::nalphabets).read() as usize;
    let ptr = addr_of!(mafft_sys::n_disFFT).read();
    if ptr.is_null() {
        return vec![vec![0; nalpha]; nalpha];
    }

    let mut matrix = vec![vec![0i32; nalpha]; nalpha];
    for i in 0..nalpha {
        let row_ptr = *ptr.add(i);
        for j in 0..nalpha {
            matrix[i][j] = *row_ptr.add(j);
        }
    }
    matrix
}

// ---------------------------------------------------------------------------
// Basic validation tests (no FFI)
// ---------------------------------------------------------------------------

#[test]
fn build_context_blosum62_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    assert_eq!(ctx.substitution_matrix.len(), 26);
    assert_eq!(ctx.substitution_matrix[0].len(), 26);
    for i in 0..20 {
        assert!(ctx.substitution_matrix[i][i] > 0, "diagonal [{i}][{i}] = {}", ctx.substitution_matrix[i][i]);
    }
    for i in 0..26 {
        for j in 0..26 {
            assert_eq!(ctx.substitution_matrix[i][j], ctx.substitution_matrix[j][i], "asymmetric at [{i}][{j}]");
        }
    }
}

#[test]
fn build_context_jtt_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Jtt, SeqType::Protein);
    assert_eq!(ctx.nalphabets, 26);
    assert_eq!(ctx.nscoredalphabets, 20);
    for i in 0..20 { assert!(ctx.substitution_matrix[i][i] > 0); }
}

#[test]
fn build_context_tm_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Tm, SeqType::Protein);
    for i in 0..20 { assert!(ctx.substitution_matrix[i][i] > 0); }
}

#[test]
fn build_context_dna_produces_valid_matrix() {
    let ctx = build_context(ScoringModel::Dna, SeqType::Dna);
    assert_eq!(ctx.nalphabets, 26);
    assert_eq!(ctx.nscoredalphabets, 10);
    assert!(ctx.substitution_matrix[0][0] > 0);
    assert!(ctx.substitution_matrix[0][1] > ctx.substitution_matrix[0][2]);
}

#[test]
fn amino_map_covers_standard_residues() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    for &ch in b"ARNDCQEGHILKMFPSTWYVBZX" {
        assert_ne!(ctx.amino_map[ch as usize], 0xFF, "unmapped: {}", ch as char);
    }
    assert_ne!(ctx.amino_map[b'-' as usize], 0xFF);
}

#[test]
fn dna_amino_map_covers_nucleotides() {
    let ctx = build_context(ScoringModel::Dna, SeqType::Dna);
    for &ch in b"agctACGT" {
        assert_ne!(ctx.amino_map[ch as usize], 0xFF, "unmapped: {}", ch as char);
    }
}

#[test]
fn polarity_volume_set_for_protein() {
    let ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    assert!(ctx.polarity.iter().map(|x| x.abs()).sum::<f64>() > 0.0);
    assert!(ctx.volume.iter().map(|x| x.abs()).sum::<f64>() > 0.0);
}

// ---------------------------------------------------------------------------
// Cell-by-cell FFI cross-validation
// ---------------------------------------------------------------------------

#[test]
fn cross_validate_blosum62_n_dis_cell_by_cell() {
    let _lock = C_MUTEX.lock().unwrap();
    let rust_ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    let c_matrix = unsafe {
        init_c_globals();
        call_c_constants(b'p', 1, 62);
        let m = read_c_n_dis();
        mafft_sys::freeconstants();
        m
    };

    let nalpha = 26;
    assert_eq!(c_matrix.len(), nalpha);
    assert_eq!(rust_ctx.substitution_matrix.len(), nalpha);

    let mut mismatches = Vec::new();
    for i in 0..nalpha {
        for j in 0..nalpha {
            let c_val = c_matrix[i][j];
            let r_val = rust_ctx.substitution_matrix[i][j];
            if c_val != r_val {
                mismatches.push((i, j, c_val, r_val));
            }
        }
    }

    if !mismatches.is_empty() {
        let total_cells = nalpha * nalpha;
        let n_mismatch = mismatches.len();
        let max_diff = mismatches.iter().map(|(_, _, c, r)| (c - r).abs()).max().unwrap_or(0);

        eprintln!(
            "BLOSUM62 n_dis: {n_mismatch}/{total_cells} cells differ (max diff = {max_diff})"
        );
        // Print first 10 mismatches
        for (i, j, c, r) in mismatches.iter().take(10) {
            eprintln!("  n_dis[{i}][{j}]: C={c}, Rust={r} (diff={})", c - r);
        }

        // Allow small differences (rounding) but flag large ones
        assert!(
            max_diff <= 2,
            "BLOSUM62 n_dis has cells differing by more than 2: max_diff={max_diff}, {n_mismatch} mismatches"
        );
    } else {
        eprintln!("BLOSUM62 n_dis: all {}/{} cells match exactly", nalpha * nalpha, nalpha * nalpha);
    }
}

#[test]
fn cross_validate_blosum62_n_dis_fft_cell_by_cell() {
    let _lock = C_MUTEX.lock().unwrap();
    let rust_ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    let c_matrix = unsafe {
        init_c_globals();
        call_c_constants(b'p', 1, 62);
        let m = read_c_n_dis_fft();
        mafft_sys::freeconstants();
        m
    };

    let nalpha = 26;
    let mut mismatches = Vec::new();
    for i in 0..nalpha {
        for j in 0..nalpha {
            let c_val = c_matrix[i][j];
            let r_val = rust_ctx.fft_matrix[i][j];
            if c_val != r_val {
                mismatches.push((i, j, c_val, r_val));
            }
        }
    }

    if !mismatches.is_empty() {
        let n_mismatch = mismatches.len();
        let max_diff = mismatches.iter().map(|(_, _, c, r)| (c - r).abs()).max().unwrap_or(0);
        eprintln!("BLOSUM62 n_disFFT: {n_mismatch}/{} cells differ (max diff = {max_diff})", nalpha * nalpha);
        for (i, j, c, r) in mismatches.iter().take(10) {
            eprintln!("  n_disFFT[{i}][{j}]: C={c}, Rust={r} (diff={})", c - r);
        }
        assert!(max_diff <= 2, "BLOSUM62 n_disFFT max_diff={max_diff}");
    } else {
        eprintln!("BLOSUM62 n_disFFT: all cells match exactly");
    }
}

#[test]
fn cross_validate_jtt_n_dis_cell_by_cell() {
    let _lock = C_MUTEX.lock().unwrap();
    let rust_ctx = build_context(ScoringModel::Jtt, SeqType::Protein);

    let c_matrix = unsafe {
        init_c_globals();
        call_c_constants(b'p', 0, 0);
        let m = read_c_n_dis();
        mafft_sys::freeconstants();
        m
    };

    let nalpha = 26;
    let mut mismatches = Vec::new();
    for i in 0..nalpha {
        for j in 0..nalpha {
            let c_val = c_matrix[i][j];
            let r_val = rust_ctx.substitution_matrix[i][j];
            if c_val != r_val {
                mismatches.push((i, j, c_val, r_val));
            }
        }
    }

    if !mismatches.is_empty() {
        let n_mismatch = mismatches.len();
        let max_diff = mismatches.iter().map(|(_, _, c, r)| (c - r).abs()).max().unwrap_or(0);
        eprintln!("JTT n_dis: {n_mismatch}/{} cells differ (max diff = {max_diff})", nalpha * nalpha);
        for (i, j, c, r) in mismatches.iter().take(10) {
            eprintln!("  n_dis[{i}][{j}]: C={c}, Rust={r} (diff={})", c - r);
        }
        assert!(max_diff <= 2, "JTT n_dis max_diff={max_diff}");
    } else {
        eprintln!("JTT n_dis: all cells match exactly");
    }
}

#[test]
fn cross_validate_dna_n_dis_cell_by_cell() {
    let _lock = C_MUTEX.lock().unwrap();
    let rust_ctx = build_context(ScoringModel::Dna, SeqType::Dna);

    let c_matrix = unsafe {
        init_c_globals();
        call_c_constants(b'd', -1, 0);
        let m = read_c_n_dis();
        mafft_sys::freeconstants();
        m
    };

    let nalpha = 26;
    let mut mismatches = Vec::new();
    for i in 0..nalpha {
        for j in 0..nalpha {
            let c_val = c_matrix[i][j];
            let r_val = rust_ctx.substitution_matrix[i][j];
            if c_val != r_val {
                mismatches.push((i, j, c_val, r_val));
            }
        }
    }

    if !mismatches.is_empty() {
        let n_mismatch = mismatches.len();
        let max_diff = mismatches.iter().map(|(_, _, c, r)| (c - r).abs()).max().unwrap_or(0);
        eprintln!("DNA n_dis: {n_mismatch}/{} cells differ (max diff = {max_diff})", nalpha * nalpha);
        for (i, j, c, r) in mismatches.iter().take(10) {
            eprintln!("  n_dis[{i}][{j}]: C={c}, Rust={r} (diff={})", c - r);
        }
        assert!(max_diff <= 2, "DNA n_dis max_diff={max_diff}");
    } else {
        eprintln!("DNA n_dis: all cells match exactly");
    }
}

#[test]
fn cross_validate_penalty_values() {
    let _lock = C_MUTEX.lock().unwrap();
    unsafe {
        init_c_globals();
        call_c_constants(b'p', 1, 62);

        let c_penalty = addr_of!(mafft_sys::penalty).read();
        let c_offset = addr_of!(mafft_sys::offset).read();

        let rust_gap = default_protein_gap_params();

        assert!(
            (c_penalty - rust_gap.penalty).abs() <= 1,
            "penalty mismatch: C={c_penalty}, Rust={}", rust_gap.penalty
        );
        assert!(
            (c_offset - rust_gap.offset).abs() <= 1,
            "offset mismatch: C={c_offset}, Rust={}", rust_gap.offset
        );

        mafft_sys::freeconstants();
    }
}

#[test]
fn cross_validate_amino_mapping() {
    let _lock = C_MUTEX.lock().unwrap();
    let rust_ctx = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    unsafe {
        init_c_globals();
        call_c_constants(b'p', 1, 62);

        // Compare amino[] (index → char) for the 26 alphabet entries
        let c_nalpha = addr_of!(mafft_sys::nalphabets).read() as usize;
        assert_eq!(c_nalpha, rust_ctx.nalphabets, "nalphabets mismatch");

        for i in 0..c_nalpha {
            let c_char = addr_of!(mafft_sys::amino).read()[i];
            // Find which Rust index maps to this char
            let r_idx = rust_ctx.amino_map[c_char as usize];
            assert_eq!(
                r_idx as usize, i,
                "amino mapping mismatch: C amino[{i}]='{}' but Rust maps '{}' to index {}",
                c_char as char, c_char as char, r_idx
            );
        }

        mafft_sys::freeconstants();
    }
}

#[test]
fn debug_jtt_pam1_vs_c() {
    let _lock = C_MUTEX.lock().unwrap();
    
    // Build Rust PAM1 (just 1 iteration to isolate the issue)
    let rust_pam1 = mafft_scoring::jtt::build_jtt_pam_matrix(false, 1);
    
    // Build C's with pamN=1
    unsafe {
        init_c_globals();
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(0);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(0);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(1);
        std::ptr::addr_of_mut!(mafft_sys::TMorJTT).write(201); // JTT

        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);

        let c_matrix = read_c_n_dis();
        mafft_sys::freeconstants();

        // Apply same normalization to Rust PAM1
        let freq = mafft_scoring::jtt::jtt_frequencies();
        let gap = mafft_scoring::default_protein_gap_params();
        let rust_norm = mafft_scoring::build_scoring_matrix(&rust_pam1, &freq, gap.offset, true);
        
        let mut mismatches = 0;
        let mut max_diff = 0i32;
        for i in 0..20 {
            for j in 0..20 {
                let diff = (c_matrix[i][j] - rust_norm[i][j]).abs();
                if diff > 0 { mismatches += 1; }
                if diff > max_diff { max_diff = diff; }
            }
        }
        eprintln!("JTT PAM1: {mismatches}/400 cells differ, max_diff={max_diff}");
        if mismatches > 0 {
            for i in 0..3 {
                for j in 0..3 {
                    eprintln!("  [{i}][{j}]: C={}, Rust={}", c_matrix[i][j], rust_norm[i][j]);
                }
            }
        }
    }
}
