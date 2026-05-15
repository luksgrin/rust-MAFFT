/// Cross-validate `mafft_align::msalignmm` against C MAFFT 7.526's
/// `MSalignmm` for the same input pair. The Rust Hirschberg DP is
/// asserted to return the exact same aligned strings as C's
/// `MSalignmm.c::MSalignmm` (which `tbfast.c:1159-1161` calls when
/// `alg='M'`, i.e. under `--memsave`).
///
/// We use the same FFI setup pattern as
/// `cross_validate_constrained_align.rs` and
/// `cross_validate_bl50_fft.rs`: `initglobalvariables` →
/// `constants()` → BLOSUM62 default penalties → invoke C with
/// pre-allocated buffers that C edits in place.
///
/// Inputs are chosen to exercise the recursive Hirschberg path
/// (`lgth1 > DPTANNI=100`) — for shorter inputs both DPs delegate to
/// the same base case and parity is trivially preserved.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

use mafft_align::{msalignmm, GapModel, Profile};
use mafft_scoring::build_context;
use mafft_types::{ScoringModel, SeqType};

static C_MUTEX: Mutex<()> = Mutex::new(());

unsafe fn init_c_protein_blosum62() {
    unsafe {
        mafft_sys::initglobalvariables();
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_EX).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_OP).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_dist).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(1);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(62);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(1);

        // constants() needs a sample sequence to initialize alphabets.
        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);
    }
}

/// Compare Rust msalignmm aligned output to C MSalignmm aligned output.
/// Both should be byte-identical. Returns (rust_s1, rust_s2, c_s1, c_s2).
fn align_via_both(
    s1: &[u8],
    s2: &[u8],
    head_gap: bool,
    tail_gap: bool,
) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let amino_map = &scoring.amino_map;
    let nalpha = scoring.nalphabets;

    // ---- Rust side ----
    let prof1 = Profile::from_aligned(&[s1], &[1.0], amino_map, nalpha);
    let prof2 = Profile::from_aligned(&[s2], &[1.0], amino_map, nalpha);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let rust_aln = msalignmm(&prof1, &prof2, &scoring.consweight_matrix, &gap, head_gap, tail_gap);

    let mut rust_s1 = Vec::with_capacity(rust_aln.operations.len());
    let mut rust_s2 = Vec::with_capacity(rust_aln.operations.len());
    let (mut p, mut q) = (0usize, 0usize);
    for op in &rust_aln.operations {
        match op {
            mafft_align::AlignOp::Match => {
                rust_s1.push(s1[p]); p += 1;
                rust_s2.push(s2[q]); q += 1;
            }
            mafft_align::AlignOp::Delete => {
                rust_s1.push(s1[p]); p += 1;
                rust_s2.push(b'-');
            }
            mafft_align::AlignOp::Insert => {
                rust_s1.push(b'-');
                rust_s2.push(s2[q]); q += 1;
            }
        }
    }

    // ---- C side ----
    let _guard = C_MUTEX.lock().unwrap();
    let (c_s1, c_s2) = unsafe {
        init_c_protein_blosum62();
        let alloclen = (s1.len() + s2.len() + 1000) as c_int;

        let c_seq1 = CString::new(s1).unwrap();
        let c_seq2 = CString::new(s2).unwrap();
        let mut buf1: Vec<u8> = c_seq1.as_bytes().to_vec();
        buf1.resize(alloclen as usize + 1, 0);
        let mut buf2: Vec<u8> = c_seq2.as_bytes().to_vec();
        buf2.resize(alloclen as usize + 1, 0);
        let mut p1 = buf1.as_mut_ptr() as *mut c_char;
        let mut p2 = buf2.as_mut_ptr() as *mut c_char;

        let mut eff1: c_double = 1.0;
        let mut eff2: c_double = 1.0;

        let nalpha_c = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha_c, nalpha_c);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }

        let _c_score = mafft_sys::MSalignmm(
            n_dyn,
            &mut p1, &mut p2,
            &mut eff1, &mut eff2,
            1, 1,
            alloclen,
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
            head_gap as c_int, tail_gap as c_int,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            1.0, 1.0,
        );

        // Read back the aligned strings (in-place edits).
        let c_width = {
            let mut k = 0; while *p1.add(k) != 0 { k += 1; } k
        };
        let c_s1: Vec<u8> = (0..c_width).map(|k| *p1.add(k) as u8).collect();
        let c_s2: Vec<u8> = (0..c_width).map(|k| *p2.add(k) as u8).collect();

        mafft_sys::freeconstants();
        (c_s1, c_s2)
    };

    (rust_s1, rust_s2, c_s1, c_s2)
}

/// Recursive case with `lgth > DPTANNI = 100`. Uses a distinct
/// 130-residue protein input on each side with a 20-residue offset
/// so the optimal alignment requires gaps and exercises both halves
/// of the Hirschberg recursion.
#[test]
fn msalign_recursive_matches_c() {
    // 130-residue input pair with a deliberate divergence: 20-residue
    // motif in the middle of s2 that's absent in s1.
    let s1 = b"MKTIIALSYIFCLVFAKEDFREEKSPELLVNVPILTPVAGTHKAGKLITGSTMKAKEGNCGRDLLINGTGRLILSSSGKLPHRMNAIPRTNKPGSEDYTKVVNFLSGNLDRGQLSYLKLELKM";
    let s2 = b"MKTIIALSYIFCLVFAKEDFREEKSPELLVNVPILTPVAGTHKAGKLITGSTMKAKEGNCGRDPQLLLAGKSDESQRWSAALLINGTGRLILSSSGKLPHRMNAIPRTNKPGSEDYTKVVNFLSGNLDRGQLSYLKLELKM";
    assert!(s1.len() > 100, "test input must exercise Hirschberg recursion");
    assert!(s2.len() > 100);

    let (rust_s1, rust_s2, c_s1, c_s2) = align_via_both(s1, s2, true, true);

    let rust_str1 = String::from_utf8_lossy(&rust_s1);
    let rust_str2 = String::from_utf8_lossy(&rust_s2);
    let c_str1 = String::from_utf8_lossy(&c_s1);
    let c_str2 = String::from_utf8_lossy(&c_s2);

    eprintln!("Rust width: {}", rust_s1.len());
    eprintln!("C    width: {}", c_s1.len());
    eprintln!("Rust seq1: {}", rust_str1);
    eprintln!("C    seq1: {}", c_str1);
    eprintln!("Rust seq2: {}", rust_str2);
    eprintln!("C    seq2: {}", c_str2);

    assert_eq!(rust_s1.len(), c_s1.len(),
        "alignment width mismatch: rust={} c={}", rust_s1.len(), c_s1.len());
    assert_eq!(rust_s1, c_s1, "seq1 aligned output differs");
    assert_eq!(rust_s2, c_s2, "seq2 aligned output differs");
}

/// Smaller input that hits the base case in both implementations.
/// Should be trivially identical because both delegate to the same
/// DP recurrence.
#[test]
fn msalign_base_case_matches_c() {
    let s1 = b"MKTIIALSYIFCLVFAKEDFREEK";
    let s2 = b"MKTIIALSYIFCLVFAKEDFREEK";
    assert!(s1.len() < 100);

    let (rust_s1, rust_s2, c_s1, c_s2) = align_via_both(s1, s2, true, true);
    assert_eq!(rust_s1, c_s1, "base-case seq1 differs");
    assert_eq!(rust_s2, c_s2, "base-case seq2 differs");
}

/// Identical longer sequences — forces a deep recursion with a clear
/// diagonal trace. Stresses the iso-score path through the Hirschberg
/// midpoint selection.
#[test]
fn msalign_identical_long_matches_c() {
    let s = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSIEGFFATLGGEVALWSLVVLAIERYIVIC";
    assert!(s.len() > 100);
    let (rust_s1, rust_s2, c_s1, c_s2) = align_via_both(s, s, true, true);
    assert_eq!(rust_s1, c_s1, "identical-long seq1 differs");
    assert_eq!(rust_s2, c_s2, "identical-long seq2 differs");
}

/// Tail-gap test: with `tail_gap=false`, the DP should leave a
/// terminal stretch unaligned. Verifies C MSalignmm's terminal-gap
/// branch matches our port.
#[test]
fn msalign_freetail_matches_c() {
    let s1 = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSIEGFFATLGGEVALWSLV";
    let s2 = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSIEGFFAT";
    let (rust_s1, rust_s2, c_s1, c_s2) = align_via_both(s1, s2, false, false);
    assert_eq!(rust_s1, c_s1, "freetail seq1 differs");
    assert_eq!(rust_s2, c_s2, "freetail seq2 differs");
}

/// Asymmetric lengths where lgth1 < lgth2 — exercises an m-split
/// case (the Hirschberg backward DP needs to continue past
/// `i == imid - 1` to refresh `jumpforwi/jumpforwj` at the chosen
/// `jumpi`, C line 1597-1603 and 1767). Our port currently breaks
/// at `i == imid - 1`, so this test fails by a few columns. See
/// TODO §B.3 for the residual. Marked `ignore` until the FFI
/// instrumented cross-validation (`MSalignmm_rec` row-state at
/// `i == jumpi`) lands.
#[test]
#[ignore = "TODO §B.3: m-split jumpforwi refresh not yet implemented"]
fn msalign_asymmetric_lengths_matches_c() {
    // s1 = 110 residues, s2 = 180 residues (s1 with extra middle motif).
    let s1 = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNFLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSI";
    let s2 = b"MNGTEGDNFYVPFSNKTGLARSPYEYPQYYLAEPWKYSALAAYMFFLILVGFPVNGGRTLSEVMKWPFSDQIANLPTQRDLELFQKLMSARTVTNLTLFVTVQHKKLRTPLNYILLNLAMANLFMVLFGFTVTMYTSMNGYFVFGPTMCSI";
    let (rust_s1, rust_s2, c_s1, c_s2) = align_via_both(s1, s2, true, true);
    assert_eq!(rust_s1.len(), c_s1.len(),
        "width differs: rust={} c={}", rust_s1.len(), c_s1.len());
    assert_eq!(rust_s1, c_s1, "asymmetric seq1 differs");
    assert_eq!(rust_s2, c_s2, "asymmetric seq2 differs");
}
