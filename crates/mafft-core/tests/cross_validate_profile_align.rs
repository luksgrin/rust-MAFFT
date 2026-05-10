/// Cross-validate profile_align against C's MSalignmm.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

use mafft_align::{Profile, profile_align, GapModel};
use mafft_scoring::build_context;
use mafft_types::{ScoringModel, SeqType};

static C_MUTEX: Mutex<()> = Mutex::new(());

unsafe fn alloc_zeroed(size: usize) -> *mut u8 {
    let layout = std::alloc::Layout::from_size_align(size.max(8), 8).unwrap();
    unsafe { std::alloc::alloc_zeroed(layout) }
}

unsafe fn init_c_protein() {
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

        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);
    }
}

/// Build a C-style n_dynamicmtx (double** indexed by char codes).
unsafe fn build_c_dynamicmtx(scoring_matrix: &[Vec<i32>]) -> *mut *mut c_double {
    // C's n_dynamicmtx is indexed by [0..nalphabets-1][0..nalphabets-1] like n_dis.
    let nalpha = scoring_matrix.len() as c_int;
    let mtx = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
    for i in 0..scoring_matrix.len() {
        for j in 0..scoring_matrix[i].len() {
            *(*mtx.add(i)).add(j) = scoring_matrix[i][j] as f64;
        }
    }
    mtx
}

#[test]
fn profile_align_matches_c_msalignmm() {
    let _guard = C_MUTEX.lock().unwrap();

    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // Two small groups with gaps (refinement-style input).
    let group1: Vec<&[u8]> = vec![
        b"ACDEFGHIKLM",
        b"ACDE-GHIKLM",
    ];
    let group2: Vec<&[u8]> = vec![
        b"ACDE-GHIKL-",
        b"ACDEFGHIKLM",
        b"A-DEFGHIK-M",
    ];
    let w1 = vec![0.5, 0.5];
    let w2 = vec![0.33333, 0.33334, 0.33333];

    let len1 = group1[0].len();
    let len2 = group2[0].len();

    // Rust side
    let prof1 = Profile::from_aligned(&group1, &w1, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&group2, &w2, &scoring.amino_map, scoring.nalphabets);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let rust_aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    eprintln!("Rust: ops.len()={}, score={:.2}", rust_aln.operations.len(), rust_aln.score);

    // C side
    unsafe {
        init_c_protein();

        // Build sequences in writable buffers with enough alloclen.
        let alloclen = (len1 + len2) * 10;
        let c_seqs1_buf: Vec<Vec<u8>> = group1.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen, 0);
            v
        }).collect();
        let c_seqs2_buf: Vec<Vec<u8>> = group2.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen, 0);
            v
        }).collect();
        let mut c_seq1_ptrs: Vec<*mut c_char> = c_seqs1_buf.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let mut c_seq2_ptrs: Vec<*mut c_char> = c_seqs2_buf.iter().map(|v| v.as_ptr() as *mut c_char).collect();

        // Actually we need writable buffers since MSalignmm modifies them.
        let c_seq1_boxed: Vec<Box<[u8]>> = group1.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();
        let c_seq2_boxed: Vec<Box<[u8]>> = group2.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();

        // Build mutable pointer arrays for the boxed slices.
        let mut c_seq1_ptrs: Vec<*mut c_char> = c_seq1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let mut c_seq2_ptrs: Vec<*mut c_char> = c_seq2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();

        // eff arrays
        let eff1: *mut c_double = alloc_zeroed(w1.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in w1.iter().enumerate() { *eff1.add(i) = w; }
        let eff2: *mut c_double = alloc_zeroed(w2.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in w2.iter().enumerate() { *eff2.add(i) = w; }

        // n_dynamicmtx
        let n_dyn = build_c_dynamicmtx(&scoring.substitution_matrix);

        // sgap/egap = null (triggers st_* gap count path which matches our Profile::from_aligned)
        let score = mafft_sys::MSalignmm(
            n_dyn,
            c_seq1_ptrs.as_mut_ptr(),
            c_seq2_ptrs.as_mut_ptr(),
            eff1, eff2,
            w1.len() as c_int,
            w2.len() as c_int,
            alloclen as c_int,
            std::ptr::null_mut(),  // sgap1
            std::ptr::null_mut(),  // sgap2
            std::ptr::null_mut(),  // egap1
            std::ptr::null_mut(),  // egap2
            std::ptr::null_mut(),  // chudanpt
            0,
            std::ptr::null_mut(),
            1,  // headgp
            1,  // tailgp
            std::ptr::null_mut(),  // cpmxchild0
            std::ptr::null_mut(),  // cpmxchild1
            std::ptr::null_mut(),  // cpmxresult
            1.0,  // orieff1
            1.0,  // orieff2
        );

        // Read the aligned sequences back
        let c_len = {
            let s = c_seq1_ptrs[0];
            let mut n = 0;
            while *s.add(n) != 0 { n += 1; }
            n
        };
        eprintln!("C:    aligned_len={}, score={:.2}", c_len, score);

        // Compare widths
        let rust_width = rust_aln.operations.len();
        eprintln!("  Rust width: {}, C width: {}", rust_width, c_len);

        // Print both aligned group1 sequences for comparison
        for (i, s) in c_seq1_boxed.iter().enumerate() {
            let c_str = std::str::from_utf8(&s[..c_len]).unwrap_or("???");
            eprintln!("  C  g1[{i}]: {c_str}");
        }
        for (i, s) in c_seq2_boxed.iter().enumerate() {
            let c_str = std::str::from_utf8(&s[..c_len]).unwrap_or("???");
            eprintln!("  C  g2[{i}]: {c_str}");
        }

        // Reconstruct Rust aligned sequences from ops
        let mut rust_g1 = vec![Vec::<u8>::new(); group1.len()];
        let mut rust_g2 = vec![Vec::<u8>::new(); group2.len()];
        let mut c1 = 0;
        let mut c2 = 0;
        use mafft_align::AlignOp;
        for op in &rust_aln.operations {
            match op {
                AlignOp::Match => {
                    for (si, s) in group1.iter().enumerate() { rust_g1[si].push(s[c1]); }
                    for (si, s) in group2.iter().enumerate() { rust_g2[si].push(s[c2]); }
                    c1 += 1; c2 += 1;
                }
                AlignOp::Delete => {
                    for (si, s) in group1.iter().enumerate() { rust_g1[si].push(s[c1]); }
                    for si in 0..group2.len() { rust_g2[si].push(b'-'); }
                    c1 += 1;
                }
                AlignOp::Insert => {
                    for si in 0..group1.len() { rust_g1[si].push(b'-'); }
                    for (si, s) in group2.iter().enumerate() { rust_g2[si].push(s[c2]); }
                    c2 += 1;
                }
            }
        }
        for (i, s) in rust_g1.iter().enumerate() {
            eprintln!("  R  g1[{i}]: {}", String::from_utf8_lossy(s));
        }
        for (i, s) in rust_g2.iter().enumerate() {
            eprintln!("  R  g2[{i}]: {}", String::from_utf8_lossy(s));
        }

        mafft_sys::freeconstants();

        // Require exact width match
        assert_eq!(rust_width, c_len, "width mismatch: rust={rust_width} c={c_len}");
    }
}

/// Test profile_align with asymmetric groups (1 vs many) — common refinement case.
#[test]
fn profile_align_1_vs_many_matches_c() {
    let _guard = C_MUTEX.lock().unwrap();

    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // 1 vs 5 split — typical refinement 1-vs-N case.
    let group1: Vec<&[u8]> = vec![
        b"----MNGTEGDNF-YVPFSNK-TGL-ARSPYEYPQY----",
    ];
    let group2: Vec<&[u8]> = vec![
        b"MN--GTEGDNFYVPFSNKTGLARSPYE-------YPQYAE",
        b"MNGTEGDNFYVPFS----NKTGLARSPYEYPQ---Y--AE",
        b"-MN-GTEGDNFYVPFSNKTGLARSPYEYPQY---------",
        b"--MNGTEGDNFYVPFSNKTGL--ARSPYEYPQY-------",
        b"-M-NGTEGDNFYVPFSNKTGLARSPYEYPQY---------",
    ];
    let w1 = vec![1.0];
    let w2 = vec![0.2; 5];

    let len1 = group1[0].len();
    let len2 = group2[0].len();
    assert_eq!(len1, len2);

    // Rust side
    let prof1 = Profile::from_aligned(&group1, &w1, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&group2, &w2, &scoring.amino_map, scoring.nalphabets);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let rust_aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);

    // C side
    unsafe {
        init_c_protein();

        let alloclen = (len1 + len2) * 10;
        let c_seq1_boxed: Vec<Box<[u8]>> = group1.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();
        let c_seq2_boxed: Vec<Box<[u8]>> = group2.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();
        let mut c_seq1_ptrs: Vec<*mut c_char> = c_seq1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let mut c_seq2_ptrs: Vec<*mut c_char> = c_seq2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();

        let eff1: *mut c_double = alloc_zeroed(w1.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in w1.iter().enumerate() { *eff1.add(i) = w; }
        let eff2: *mut c_double = alloc_zeroed(w2.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in w2.iter().enumerate() { *eff2.add(i) = w; }

        let n_dyn = build_c_dynamicmtx(&scoring.substitution_matrix);

        let c_score = mafft_sys::MSalignmm(
            n_dyn,
            c_seq1_ptrs.as_mut_ptr(),
            c_seq2_ptrs.as_mut_ptr(),
            eff1, eff2,
            w1.len() as c_int,
            w2.len() as c_int,
            alloclen as c_int,
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
            1, 1,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            1.0, 1.0,
        );

        let c_len = {
            let s = c_seq1_ptrs[0];
            let mut n = 0;
            while *s.add(n) != 0 { n += 1; }
            n
        };

        eprintln!("Rust: ops={} score={:.2}", rust_aln.operations.len(), rust_aln.score);
        eprintln!("C:    len={} score={:.2}", c_len, c_score);

        for (i, s) in c_seq1_boxed.iter().enumerate() {
            eprintln!("  C  g1[{i}]: {}", String::from_utf8_lossy(&s[..c_len]));
        }
        for (i, s) in c_seq2_boxed.iter().enumerate() {
            eprintln!("  C  g2[{i}]: {}", String::from_utf8_lossy(&s[..c_len]));
        }

        mafft_sys::freeconstants();

        assert_eq!(rust_aln.operations.len(), c_len, "width mismatch");
        assert!((rust_aln.score - c_score).abs() < 1e-3, "score mismatch: rust={} c={}", rust_aln.score, c_score);
    }
}
