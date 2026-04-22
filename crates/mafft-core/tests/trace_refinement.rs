/// Detailed trace of a single refinement branch: Rust vs C, every step.
///
/// Starts from the byte-identical FFT-NS-2 progressive output and traces
/// the FIRST refinement branch through both implementations.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

use mafft_align::{Profile, profile_align, GapModel};
use mafft_io::read_fasta;
use mafft_scoring::build_context;
use mafft_tree::{BranchWeights, DistanceMatrix, ClusterMethod, musclesupg, scoring_matrix_distance};
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
        // MAFFT script: -h 0.000 → poffset=0
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        // Ensure FFT window/threshold are filled in by constants() defaults.
        std::ptr::addr_of_mut!(mafft_sys::fftWinSize).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftThreshold).write(mafft_sys::NOTSPECIFIED);
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

/// Build a C topology from a Rust Topology.
/// Returns a 3D int array suitable for C's weightFromABranch.
unsafe fn build_c_topol(topo: &mafft_tree::Topology) -> *mut *mut *mut c_int {
    let nsteps = topo.steps.len();
    let topol: *mut *mut *mut c_int = alloc_zeroed(nsteps * std::mem::size_of::<*mut *mut c_int>()) as _;
    for k in 0..nsteps {
        let row: *mut *mut c_int = alloc_zeroed(2 * std::mem::size_of::<*mut c_int>()) as _;
        let left = &topo.steps[k].left;
        let left_arr: *mut c_int = alloc_zeroed((left.len() + 1) * std::mem::size_of::<c_int>()) as _;
        for (i, &s) in left.iter().enumerate() {
            *left_arr.add(i) = s as c_int;
        }
        *left_arr.add(left.len()) = -1;
        *row.add(0) = left_arr;

        let right = &topo.steps[k].right;
        let right_arr: *mut c_int = alloc_zeroed((right.len() + 1) * std::mem::size_of::<c_int>()) as _;
        for (i, &s) in right.iter().enumerate() {
            *right_arr.add(i) = s as c_int;
        }
        *right_arr.add(right.len()) = -1;
        *row.add(1) = right_arr;
        *topol.add(k) = row;
    }
    topol
}

unsafe fn build_c_len(topo: &mafft_tree::Topology) -> *mut *mut c_double {
    let nsteps = topo.steps.len();
    let len: *mut *mut c_double = alloc_zeroed(nsteps * std::mem::size_of::<*mut c_double>()) as _;
    for k in 0..nsteps {
        let row: *mut c_double = alloc_zeroed(2 * std::mem::size_of::<c_double>()) as _;
        *row.add(0) = topo.steps[k].left_length;
        *row.add(1) = topo.steps[k].right_length;
        *len.add(k) = row;
    }
    len
}

#[test]
fn trace_first_branch_vs_c() {
    let _guard = C_MUTEX.lock().unwrap();

    // Load byte-identical progressive output.
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();
    let width = sequences[0].len();
    eprintln!("=== TRACE: {nseq} seqs, width {width} ===");

    // Build distance matrix via our scoring_matrix_distance
    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(
                &sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist,
            );
            dm.set(i, j, d);
        }
    }

    // Build refinement tree (musclesupg)
    let topo = musclesupg(&dm, ClusterMethod::default());
    eprintln!("Topology: {} steps", topo.steps.len());
    for k in 0..5.min(topo.steps.len()) {
        eprintln!("  step {k}: left={:?} right={:?} ll={} rl={}",
            topo.steps[k].left, topo.steps[k].right,
            topo.steps[k].left_length, topo.steps[k].right_length);
    }

    // Rust branch weights
    let bw = BranchWeights::new(&topo);
    let rust_weights = bw.weights_for_branch(&topo, 0, 0);
    eprintln!("Rust weights for step=0 side=0: first 6 = {:?}", &rust_weights[..6.min(rust_weights.len())]);

    // C branch weights via FFI
    unsafe {
        init_c_protein();

        let topol = build_c_topol(&topo);
        let len = build_c_len(&topo);
        let bw_mtx: *mut *mut c_double = alloc_zeroed(topo.steps.len() * std::mem::size_of::<*mut c_double>()) as _;
        for k in 0..topo.steps.len() {
            let row: *mut c_double = alloc_zeroed(2 * std::mem::size_of::<c_double>()) as _;
            *row.add(0) = 1.0;
            *row.add(1) = 1.0;
            *bw_mtx.add(k) = row;
        }
        std::ptr::addr_of_mut!(mafft_sys::sueff_global).write(0.1);
        std::ptr::addr_of_mut!(mafft_sys::treemethod).write(b'X' as c_int);
        let total_nodes = 2 * nseq as usize;
        let stopol: *mut mafft_sys::Node = alloc_zeroed(total_nodes * std::mem::size_of::<mafft_sys::Node>()) as _;
        mafft_sys::treeCnv(stopol, nseq as c_int, topol, len, bw_mtx);
        mafft_sys::calcBranchWeight(bw_mtx, nseq as c_int, stopol, topol, len);

        let mut c_weights = vec![0.0f64; nseq];
        mafft_sys::weightFromABranch(nseq as c_int, c_weights.as_mut_ptr(), stopol, topol, 0, 0);
        eprintln!("C    weights for step=0 side=0: first 6 = {:?}", &c_weights[..6.min(c_weights.len())]);

        // Dump C's bw matrix (per-edge branch weights)
        eprintln!("C bw matrix (step, lor -> weight):");
        for k in 0..topo.steps.len().min(5) {
            for lor in 0..2 {
                let w = *(*bw_mtx.add(k)).add(lor);
                eprintln!("  bw[{k}][{lor}] = {}", w);
            }
        }

        // Dump Rust's branch_weight per-node (first 5 internal nodes)
        eprintln!("Rust branch_weight per internal node:");
        for k in 0..5.min(topo.steps.len()) {
            let (children, weights) = bw.debug_node(k);
            eprintln!("  node[{k}]: children={:?} bw={:?}", children, weights);
        }

        // Compare C's stopol[0] with Rust nodes[0]
        eprintln!("C stopol[0]:");
        for d in 0..3 {
            let ch = (*stopol.add(0)).children[d];
            let len = (*stopol.add(0)).length[d];
            let wp = (*stopol.add(0)).weightptr[d];
            let w = if wp.is_null() { 0.0 } else { *wp };
            let node_idx = if ch.is_null() { -1i64 } else { (ch as i64 - stopol as i64) / std::mem::size_of::<mafft_sys::Node>() as i64 };
            eprintln!("  children[{d}]: node={}, len={:.4}, bw={:.4}", node_idx, len, w);
        }

        // Compare weight vectors
        let max_w_diff = rust_weights.iter().zip(c_weights.iter())
            .map(|(r, c)| (r - c).abs())
            .fold(0.0f64, f64::max);
        eprintln!("Weight max diff: {:.2e}", max_w_diff);
        assert!(max_w_diff < 1e-6, "weights diverge");

        // Determine group1 and group2 for this branch.
        let group1: Vec<usize> = topo.steps[0].left.clone();
        let mut group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
        eprintln!("group1 ({}): {:?}", group1.len(), group1);
        eprintln!("group2 ({}): first 5 = {:?}", group2.len(), &group2[..5.min(group2.len())]);

        // ==============================================
        // STEP: Per-group normalized weights (Rust side)
        // ==============================================
        const MIN_W: f64 = 0.00001;
        let rw1: Vec<f64> = group1.iter().map(|&i| rust_weights[i].max(MIN_W)).collect();
        let rw2: Vec<f64> = group2.iter().map(|&i| rust_weights[i].max(MIN_W)).collect();
        let rs1: f64 = rw1.iter().sum();
        let rs2: f64 = rw2.iter().sum();
        let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
        let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

        // Per-group normalized weights (C side via fastconjuction_noname)
        let eff: *mut c_double = alloc_zeroed(nseq * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in c_weights.iter().enumerate() { *eff.add(i) = w; }

        let mut memlist1: Vec<c_int> = group1.iter().map(|&i| i as c_int).chain(std::iter::once(-1)).collect();
        let mut memlist2: Vec<c_int> = group2.iter().map(|&i| i as c_int).chain(std::iter::once(-1)).collect();

        // Need dummy seq arrays
        let cseqs_buf: Vec<CString> = sequences.iter().map(|s| CString::new(s.clone()).unwrap()).collect();
        let mut cseqs_ptrs: Vec<*mut c_char> = cseqs_buf.iter().map(|c| c.as_ptr() as *mut c_char).collect();
        let mut aseq1: Vec<*mut c_char> = vec![std::ptr::null_mut(); nseq];
        let mut aseq2: Vec<*mut c_char> = vec![std::ptr::null_mut(); nseq];
        let mut peff1 = vec![0.0f64; nseq];
        let mut peff2 = vec![0.0f64; nseq];
        let d1: *mut c_char = alloc_zeroed(1000) as _;
        let d2: *mut c_char = alloc_zeroed(1000) as _;

        let cc1 = mafft_sys::fastconjuction_noname(
            memlist1.as_mut_ptr(), cseqs_ptrs.as_mut_ptr(), aseq1.as_mut_ptr(),
            peff1.as_mut_ptr(), eff, d1, 0.00001, std::ptr::null_mut(),
        );
        let cc2 = mafft_sys::fastconjuction_noname(
            memlist2.as_mut_ptr(), cseqs_ptrs.as_mut_ptr(), aseq2.as_mut_ptr(),
            peff2.as_mut_ptr(), eff, d2, 0.00001, std::ptr::null_mut(),
        );

        assert_eq!(cc1 as usize, group1.len());
        assert_eq!(cc2 as usize, group2.len());

        eprintln!("Rust w1n: {:?}", &rw1n[..3.min(rw1n.len())]);
        eprintln!("C    w1n: {:?}", &peff1[..3.min(cc1 as usize)]);
        eprintln!("Rust w2n: first 3 = {:?}", &rw2n[..3.min(rw2n.len())]);
        eprintln!("C    w2n: first 3 = {:?}", &peff2[..3.min(cc2 as usize)]);

        let max_w1n_diff = rw1n.iter().zip(peff1.iter())
            .map(|(r, c)| (r - c).abs()).fold(0.0f64, f64::max);
        let max_w2n_diff = rw2n.iter().zip(peff2.iter())
            .map(|(r, c)| (r - c).abs()).fold(0.0f64, f64::max);
        eprintln!("w1n max diff: {max_w1n_diff:.2e}");
        eprintln!("w2n max diff: {max_w2n_diff:.2e}");
        assert!(max_w1n_diff < 1e-6);
        assert!(max_w2n_diff < 1e-6);

        // ==============================================
        // STEP: old_score via intergroup_score
        // ==============================================
        // Rust: compute via pairwise_score (port of intergroup_score)
        let rust_old_score = compute_intergroup_score_rust(
            &group1, &group2, &sequences, &rw1n, &rw2n, &scoring,
        );

        // C: call C's intergroup_score directly
        // Need char** arrays for group1 and group2 sequences
        let g1_seqs: Vec<CString> = group1.iter().map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();
        let g2_seqs: Vec<CString> = group2.iter().map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();
        let mut g1_ptrs: Vec<*mut c_char> = g1_seqs.iter().map(|c| c.as_ptr() as *mut c_char).collect();
        let mut g2_ptrs: Vec<*mut c_char> = g2_seqs.iter().map(|c| c.as_ptr() as *mut c_char).collect();

        let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in peff1[..cc1 as usize].iter().enumerate() { *e1.add(i) = w; }
        let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in peff2[..cc2 as usize].iter().enumerate() { *e2.add(i) = w; }

        let mut c_old_score = 0.0f64;
        mafft_sys::intergroup_score(
            g1_ptrs.as_mut_ptr(), g2_ptrs.as_mut_ptr(),
            e1, e2,
            group1.len() as c_int, group2.len() as c_int,
            width as c_int,
            &mut c_old_score as *mut c_double,
        );

        eprintln!("Rust old_score: {:.4}", rust_old_score);
        eprintln!("C    old_score: {:.4}", c_old_score);
        let old_score_diff = (rust_old_score - c_old_score).abs();
        eprintln!("old_score diff: {:.2e}", old_score_diff);
        assert!(old_score_diff < 1e-2, "old_score diverges");

        // ==============================================
        // STEP: Realignment (Falign in C / realign_all in Rust)
        // We'll compare just the profile_align part on stripped profiles
        // ==============================================

        // Gap strip both groups per-group
        let gap1_cols: Vec<bool> = (0..width).map(|c| group1.iter().all(|&i| sequences[i][c] == b'-')).collect();
        let gap2_cols: Vec<bool> = (0..width).map(|c| group2.iter().all(|&i| sequences[i][c] == b'-')).collect();
        let kept1: Vec<usize> = (0..width).filter(|&c| !gap1_cols[c]).collect();
        let kept2: Vec<usize> = (0..width).filter(|&c| !gap2_cols[c]).collect();
        eprintln!("kept1: {} cols, kept2: {} cols", kept1.len(), kept2.len());

        let stripped1: Vec<Vec<u8>> = group1.iter().map(|&i| kept1.iter().map(|&c| sequences[i][c]).collect()).collect();
        let stripped2: Vec<Vec<u8>> = group2.iter().map(|&i| kept2.iter().map(|&c| sequences[i][c]).collect()).collect();
        let s1r: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
        let s2r: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();

        // Rust profile_align
        let prof1 = Profile::from_aligned(&s1r, &rw1n, &scoring.amino_map, scoring.nalphabets);
        let prof2 = Profile::from_aligned(&s2r, &rw2n, &scoring.amino_map, scoring.nalphabets);
        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let rust_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
        eprintln!("Rust profile_align: ops={} score={:.2}", rust_aln.operations.len(), rust_aln.score);

        // C MSalignmm
        let alloclen = width * 10;
        let c_s1_boxed: Vec<Box<[u8]>> = stripped1.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();
        let c_s2_boxed: Vec<Box<[u8]>> = stripped2.iter().map(|s| {
            let mut v = s.to_vec();
            v.resize(alloclen + 1, 0);
            v.into_boxed_slice()
        }).collect();
        let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();

        let n_dyn = {
            let nalpha = scoring.substitution_matrix.len() as c_int;
            let m = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
            for i in 0..scoring.substitution_matrix.len() {
                for j in 0..scoring.substitution_matrix[i].len() {
                    *(*m.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
                }
            }
            m
        };

        let c_msa_score = mafft_sys::MSalignmm(
            n_dyn,
            c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
            e1, e2,
            group1.len() as c_int, group2.len() as c_int,
            alloclen as c_int,
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
            1, 1,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            1.0, 1.0,
        );

        let c_aln_len = {
            let s = c_s1_ptrs[0];
            let mut n = 0;
            while *s.add(n) != 0 { n += 1; }
            n
        };

        eprintln!("C    MSalignmm: ops={} score={:.2}", c_aln_len, c_msa_score);
        eprintln!("aln ops diff: rust={}, c={}, diff={}",
            rust_aln.operations.len(), c_aln_len,
            rust_aln.operations.len() as i32 - c_aln_len as i32);
        eprintln!("aln score diff: {:.2e}", (rust_aln.score - c_msa_score).abs());

        // If widths differ, compare sequences
        if rust_aln.operations.len() != c_aln_len {
            eprintln!("*** WIDTH DIVERGENCE in profile_align ***");
            // Print first aligned seq of each group from both
            for (i, s) in c_s1_boxed.iter().enumerate().take(2) {
                eprintln!("  C  g1[{i}]: {}", String::from_utf8_lossy(&s[..c_aln_len]));
            }
        }

        // Reconstruct aligned sequences from rust_aln.operations
        let rust_aln_len = rust_aln.operations.len();
        let mut rust_aligned_g1: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln_len); group1.len()];
        let mut rust_aligned_g2: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln_len); group2.len()];
        let mut i1 = 0usize;
        let mut i2 = 0usize;
        for op in &rust_aln.operations {
            match op {
                mafft_align::AlignOp::Match => {
                    for (k, s) in stripped1.iter().enumerate() { rust_aligned_g1[k].push(s[i1]); }
                    for (k, s) in stripped2.iter().enumerate() { rust_aligned_g2[k].push(s[i2]); }
                    i1 += 1; i2 += 1;
                }
                mafft_align::AlignOp::Insert => {
                    for (k, _) in stripped1.iter().enumerate() { rust_aligned_g1[k].push(b'-'); }
                    for (k, s) in stripped2.iter().enumerate() { rust_aligned_g2[k].push(s[i2]); }
                    i2 += 1;
                }
                mafft_align::AlignOp::Delete => {
                    for (k, s) in stripped1.iter().enumerate() { rust_aligned_g1[k].push(s[i1]); }
                    for (k, _) in stripped2.iter().enumerate() { rust_aligned_g2[k].push(b'-'); }
                    i1 += 1;
                }
            }
        }

        let mut seq_mismatches = 0usize;
        let mut first_col_diff: Option<usize> = None;
        if rust_aln_len == c_aln_len {
            for i in 0..group1.len() {
                let rust_s = &rust_aligned_g1[i];
                let c_s = &c_s1_boxed[i][..c_aln_len];
                for c in 0..rust_s.len() {
                    if rust_s[c] != c_s[c] {
                        if first_col_diff.is_none() { first_col_diff = Some(c); }
                        seq_mismatches += 1;
                    }
                }
            }
            for i in 0..group2.len() {
                let rust_s = &rust_aligned_g2[i];
                let c_s = &c_s2_boxed[i][..c_aln_len];
                for c in 0..rust_s.len() {
                    if rust_s[c] != c_s[c] {
                        if first_col_diff.is_none() { first_col_diff = Some(c); }
                        seq_mismatches += 1;
                    }
                }
            }
        }
        eprintln!("Sequence byte mismatches: {} (first col diff: {:?})", seq_mismatches, first_col_diff);
        if let Some(c) = first_col_diff {
            let show_start = c.saturating_sub(10);
            let show_end = (c + 10).min(rust_aln_len);
            eprintln!("  rust g1[0][{}..{}]: {}", show_start, show_end,
                String::from_utf8_lossy(&rust_aligned_g1[0][show_start..show_end]));
            eprintln!("  c    g1[0][{}..{}]: {}", show_start, show_end,
                String::from_utf8_lossy(&c_s1_boxed[0][show_start..show_end]));
        }

        mafft_sys::freeconstants();
    }
}

#[test]
fn all_branches_weights_match_c_36seq() {
    let _guard = C_MUTEX.lock().unwrap();

    // Load byte-identical progressive output.
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(
                &sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist,
            );
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let topol = build_c_topol(&topo);
        let len = build_c_len(&topo);
        let bw_mtx: *mut *mut c_double = alloc_zeroed(topo.steps.len() * std::mem::size_of::<*mut c_double>()) as _;
        for k in 0..topo.steps.len() {
            let row: *mut c_double = alloc_zeroed(2 * std::mem::size_of::<c_double>()) as _;
            *row.add(0) = 1.0;
            *row.add(1) = 1.0;
            *bw_mtx.add(k) = row;
        }
        std::ptr::addr_of_mut!(mafft_sys::sueff_global).write(0.1);
        std::ptr::addr_of_mut!(mafft_sys::treemethod).write(b'X' as c_int);
        let total_nodes = 2 * nseq as usize;
        let stopol: *mut mafft_sys::Node = alloc_zeroed(total_nodes * std::mem::size_of::<mafft_sys::Node>()) as _;
        mafft_sys::treeCnv(stopol, nseq as c_int, topol, len, bw_mtx);
        mafft_sys::calcBranchWeight(bw_mtx, nseq as c_int, stopol, topol, len);

        let mut max_diff = 0.0f64;
        let mut first_mismatch: Option<(usize, usize, usize, f64, f64)> = None;
        let mut mismatches = 0usize;
        for k in 0..topo.steps.len() {
            for side in 0..2u32 {
                let rust_w = bw.weights_for_branch(&topo, k, side as usize);
                let mut c_w = vec![0.0f64; nseq];
                mafft_sys::weightFromABranch(
                    nseq as c_int, c_w.as_mut_ptr(), stopol, topol,
                    k as c_int, side as c_int,
                );
                for seq_i in 0..nseq {
                    let d = (rust_w[seq_i] - c_w[seq_i]).abs();
                    if d > max_diff { max_diff = d; }
                    if d > 1e-6 {
                        mismatches += 1;
                        if first_mismatch.is_none() {
                            first_mismatch = Some((k, side as usize, seq_i, rust_w[seq_i], c_w[seq_i]));
                        }
                    }
                }
            }
        }
        eprintln!("Total 36-seq branch weight max diff: {:.2e}", max_diff);
        eprintln!("Mismatches > 1e-6: {}", mismatches);
        if let Some((k, s, si, r, c)) = first_mismatch {
            eprintln!("First mismatch: step={k} side={s} seq={si} rust={r} c={c}");
        }
        assert!(max_diff < 1e-6, "weights diverge for some branch");
        mafft_sys::freeconstants();
    }
}

/// Run the Rust refinement loop for iter 0 branch-by-branch, and at each
/// branch, compare our profile_align output against C's MSalignmm on the
/// current (evolving) state. This catches state-dependent divergences that
/// a clean-state sweep misses.
#[test]
fn evolving_state_profile_align_vs_msalignmm_iter0() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let mut sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(
                &sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist,
            );
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let topol = build_c_topol(&topo);
        let len = build_c_len(&topo);
        let bw_mtx: *mut *mut c_double = alloc_zeroed(topo.steps.len() * std::mem::size_of::<*mut c_double>()) as _;
        for k in 0..topo.steps.len() {
            let row: *mut c_double = alloc_zeroed(2 * std::mem::size_of::<c_double>()) as _;
            *row.add(0) = 1.0;
            *row.add(1) = 1.0;
            *bw_mtx.add(k) = row;
        }
        std::ptr::addr_of_mut!(mafft_sys::sueff_global).write(0.1);
        std::ptr::addr_of_mut!(mafft_sys::treemethod).write(b'X' as c_int);
        let stopol: *mut mafft_sys::Node = alloc_zeroed(2 * nseq * std::mem::size_of::<mafft_sys::Node>()) as _;
        mafft_sys::treeCnv(stopol, nseq as c_int, topol, len, bw_mtx);
        mafft_sys::calcBranchWeight(bw_mtx, nseq as c_int, stopol, topol, len);

        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }
        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let nsteps = topo.steps.len();

        let mut branches: Vec<(usize, usize)> = Vec::new();
        for step in 0..nsteps {
            if step == nsteps - 1 { branches.push((step, 1)); }
            else { branches.push((step, 0)); branches.push((step, 1)); }
        }

        let mut divergences = 0usize;
        let mut first_diverge: Option<(usize, usize)> = None;
        for (step, side) in &branches {
            let group1: Vec<usize> = if *side == 0 { topo.steps[*step].left.clone() } else { topo.steps[*step].right.clone() };
            let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
            if group1.is_empty() || group2.is_empty() { continue; }

            let weights = bw.weights_for_branch(&topo, *step, *side);
            const MIN_W: f64 = 0.00001;
            let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rs1: f64 = rw1.iter().sum();
            let rs2: f64 = rw2.iter().sum();
            let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
            let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

            let width = sequences[0].len();
            let gap1: Vec<bool> = (0..width).map(|c| group1.iter().all(|&i| sequences[i][c] == b'-')).collect();
            let gap2: Vec<bool> = (0..width).map(|c| group2.iter().all(|&i| sequences[i][c] == b'-')).collect();
            let kept1: Vec<usize> = (0..width).filter(|&c| !gap1[c]).collect();
            let kept2: Vec<usize> = (0..width).filter(|&c| !gap2[c]).collect();
            let stripped1: Vec<Vec<u8>> = group1.iter().map(|&i| kept1.iter().map(|&c| sequences[i][c]).collect()).collect();
            let stripped2: Vec<Vec<u8>> = group2.iter().map(|&i| kept2.iter().map(|&c| sequences[i][c]).collect()).collect();
            if stripped1.is_empty() || stripped1[0].is_empty() || stripped2.is_empty() || stripped2[0].is_empty() { continue; }
            let s1r: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
            let s2r: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();

            let prof1 = Profile::from_aligned(&s1r, &rw1n, &scoring.amino_map, scoring.nalphabets);
            let prof2 = Profile::from_aligned(&s2r, &rw2n, &scoring.amino_map, scoring.nalphabets);
            let rust_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

            let mut r_g1: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln.operations.len()); group1.len()];
            let mut r_g2: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln.operations.len()); group2.len()];
            let mut i1 = 0usize; let mut i2 = 0usize;
            for op in &rust_aln.operations {
                match op {
                    mafft_align::AlignOp::Match => {
                        for (k, s) in stripped1.iter().enumerate() { r_g1[k].push(s[i1]); }
                        for (k, s) in stripped2.iter().enumerate() { r_g2[k].push(s[i2]); }
                        i1 += 1; i2 += 1;
                    }
                    mafft_align::AlignOp::Insert => {
                        for r in r_g1.iter_mut() { r.push(b'-'); }
                        for (k, s) in stripped2.iter().enumerate() { r_g2[k].push(s[i2]); }
                        i2 += 1;
                    }
                    mafft_align::AlignOp::Delete => {
                        for (k, s) in stripped1.iter().enumerate() { r_g1[k].push(s[i1]); }
                        for r in r_g2.iter_mut() { r.push(b'-'); }
                        i1 += 1;
                    }
                }
            }

            let alloclen = width * 10;
            let c_s1_boxed: Vec<Box<[u8]>> = stripped1.iter().map(|s| {
                let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let c_s2_boxed: Vec<Box<[u8]>> = stripped2.iter().map(|s| {
                let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
            let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }
            let _ = mafft_sys::MSalignmm(
                n_dyn, c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
                e1, e2, group1.len() as c_int, group2.len() as c_int,
                alloclen as c_int, std::ptr::null_mut(), std::ptr::null_mut(),
                std::ptr::null_mut(), std::ptr::null_mut(),
                std::ptr::null_mut(), 0, std::ptr::null_mut(), 1, 1,
                std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
                1.0, 1.0,
            );

            let c_aln_len = {
                let s = c_s1_ptrs[0];
                let mut n = 0; while *s.add(n) != 0 { n += 1; } n
            };
            let rust_len = rust_aln.operations.len();

            let mut divergent = false;
            if rust_len != c_aln_len { divergent = true; }
            else {
                'outer: for i in 0..group1.len() {
                    for c in 0..rust_len {
                        if r_g1[i][c] != c_s1_boxed[i][c] { divergent = true; break 'outer; }
                    }
                }
                if !divergent {
                    'outer2: for i in 0..group2.len() {
                        for c in 0..rust_len {
                            if r_g2[i][c] != c_s2_boxed[i][c] { divergent = true; break 'outer2; }
                        }
                    }
                }
            }
            if divergent {
                divergences += 1;
                if first_diverge.is_none() { first_diverge = Some((*step, *side)); }
            }

            // Apply Rust's accept/reject decision to evolve state.
            // We simulate Rust's behavior: compute tscore and old_score, accept if tscore > old_score.
            // Reconstruct the new full-width sequences from the alignment.
            let new_width = rust_aln.operations.len();
            let mut new_seqs = vec![vec![b'-'; new_width]; nseq];
            let mut i1 = 0usize; let mut i2 = 0usize;
            for (col, op) in rust_aln.operations.iter().enumerate() {
                match op {
                    mafft_align::AlignOp::Match => {
                        for (k, &idx) in group1.iter().enumerate() { new_seqs[idx][col] = stripped1[k][i1]; }
                        for (k, &idx) in group2.iter().enumerate() { new_seqs[idx][col] = stripped2[k][i2]; }
                        i1 += 1; i2 += 1;
                    }
                    mafft_align::AlignOp::Insert => {
                        for &idx in group1.iter() { new_seqs[idx][col] = b'-'; }
                        for (k, &idx) in group2.iter().enumerate() { new_seqs[idx][col] = stripped2[k][i2]; }
                        i2 += 1;
                    }
                    mafft_align::AlignOp::Delete => {
                        for (k, &idx) in group1.iter().enumerate() { new_seqs[idx][col] = stripped1[k][i1]; }
                        for &idx in group2.iter() { new_seqs[idx][col] = b'-'; }
                        i1 += 1;
                    }
                }
            }

            // Compute simple intergroup scores for accept decision
            let old_score = compute_intergroup_score_rust(&group1, &group2, &sequences, &rw1n, &rw2n, &scoring);
            let tscore = compute_intergroup_score_rust(&group1, &group2, &new_seqs, &rw1n, &rw2n, &scoring);
            let changed = (0..nseq).any(|i| sequences[i] != new_seqs[i]);
            let decision = if !changed { "identical" } else if tscore > old_score { "accepted" } else { "rejected" };
            let common_gap_cols = (0..width).filter(|&c| gap1[c] && gap2[c]).count();
            eprintln!("step={} side={} {} old={:.3} tscore={:.3} diff={:.3} div={} cgc={} w={} k1={} k2={} newW={}",
                step, side, decision, old_score, tscore, tscore - old_score, divergent,
                common_gap_cols, width, kept1.len(), kept2.len(), rust_aln.operations.len());
            if changed && tscore > old_score {
                sequences = new_seqs;
            }
        }
        eprintln!("Evolving-state divergent branches (iter 0): {}/{}", divergences, branches.len());
        if let Some((s, sd)) = first_diverge {
            eprintln!("First divergence: step={} side={}", s, sd);
        }
        mafft_sys::freeconstants();
    }
}

/// Test profile_align on FULL (non-stripped) sequences vs C MSalignmm on FULL.
#[test]
fn full_sequence_profile_align_vs_msalignmm() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(&sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist);
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }
        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let nsteps = topo.steps.len();

        let mut divergences = 0;
        let mut total = 0;
        for step in 0..nsteps {
            let sides: Vec<usize> = if step == nsteps - 1 { vec![1] } else { vec![0, 1] };
            for side in sides {
                let group1: Vec<usize> = if side == 0 { topo.steps[step].left.clone() } else { topo.steps[step].right.clone() };
                let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
                if group1.is_empty() || group2.is_empty() { continue; }

                let weights = bw.weights_for_branch(&topo, step, side);
                const MIN_W: f64 = 0.00001;
                let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
                let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
                let rs1: f64 = rw1.iter().sum();
                let rs2: f64 = rw2.iter().sum();
                let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
                let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

                // Use FULL (non-stripped) sequences.
                let full1: Vec<Vec<u8>> = group1.iter().map(|&i| sequences[i].clone()).collect();
                let full2: Vec<Vec<u8>> = group2.iter().map(|&i| sequences[i].clone()).collect();
                let f1r: Vec<&[u8]> = full1.iter().map(|s| s.as_slice()).collect();
                let f2r: Vec<&[u8]> = full2.iter().map(|s| s.as_slice()).collect();

                let prof1 = Profile::from_aligned(&f1r, &rw1n, &scoring.amino_map, scoring.nalphabets);
                let prof2 = Profile::from_aligned(&f2r, &rw2n, &scoring.amino_map, scoring.nalphabets);
                let rust_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

                let width = sequences[0].len();
                let alloclen = width * 10;
                let c_s1_boxed: Vec<Box<[u8]>> = full1.iter().map(|s| {
                    let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
                }).collect();
                let c_s2_boxed: Vec<Box<[u8]>> = full2.iter().map(|s| {
                    let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
                }).collect();
                let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
                let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
                let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
                for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
                let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
                for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }
                let _ = mafft_sys::MSalignmm(
                    n_dyn, c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
                    e1, e2, group1.len() as c_int, group2.len() as c_int,
                    alloclen as c_int, std::ptr::null_mut(), std::ptr::null_mut(),
                    std::ptr::null_mut(), std::ptr::null_mut(),
                    std::ptr::null_mut(), 0, std::ptr::null_mut(), 1, 1,
                    std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
                    1.0, 1.0,
                );
                let c_aln_len = {
                    let s = c_s1_ptrs[0];
                    let mut n = 0; while *s.add(n) != 0 { n += 1; } n
                };
                let rust_len = rust_aln.operations.len();
                total += 1;
                if rust_len != c_aln_len {
                    divergences += 1;
                    if divergences <= 3 {
                        eprintln!("step={} side={} width-diff rust={} c={}", step, side, rust_len, c_aln_len);
                    }
                }
            }
        }
        eprintln!("Full-sequence divergences: {}/{}", divergences, total);
        mafft_sys::freeconstants();
    }
}

/// Run C's real refinement pipeline for 1 iteration via FFI — compare against
/// our Rust refinement after 1 iteration to see if states diverge.
#[test]
fn rust_refine_vs_c_direct_iter1() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let initial: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    // Run our real refinement for iter=1
    let mut msa = mafft_core::MultipleAlignment {
        sequences: initial.clone(),
        names: vec![String::new(); nseq],
        score: 0.0,
        step_trace: Vec::new(),
    };
    // Build topology
    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(&initial[i], &initial[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist);
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());

    let params = mafft_core::RefinementParams {
        max_iterations: 1, use_fft: true, ..Default::default()
    };
    mafft_core::iterative_refine(&mut msa, &topo, &scoring, &params, None);

    let rust_after: Vec<Vec<u8>> = msa.sequences.clone();
    eprintln!("Rust iter 1 width: {}", rust_after[0].len());

    unsafe {
        init_c_protein();
        // Just run Falign on branch 6 side 0 with the CLEAN state to compare
        // its output to what real Rust refinement produced.
        let bw = BranchWeights::new(&topo);
        let (step, side) = (6usize, 0usize);
        let group1: Vec<usize> = topo.steps[step].left.clone();
        let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
        let weights = bw.weights_for_branch(&topo, step, side);
        const MIN_W: f64 = 0.00001;
        let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
        let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
        let rs1: f64 = rw1.iter().sum(); let rs2: f64 = rw2.iter().sum();
        let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
        let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

        let width = initial[0].len();
        let alloclen = width * 3;
        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }
        std::ptr::addr_of_mut!(mafft_sys::alg).write(b'M' as i8);
        std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(1);
        std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(1);
        std::ptr::addr_of_mut!(mafft_sys::use_fft).write(1);
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(1);

        let mut fftlog: c_int = 0;
        let c_s1_boxed: Vec<Box<[u8]>> = group1.iter().map(|&i| {
            let mut v = initial[i].clone(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
        }).collect();
        let c_s2_boxed: Vec<Box<[u8]>> = group2.iter().map(|&i| {
            let mut v = initial[i].clone(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
        }).collect();
        let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
        let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
        let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
        for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }
        let _ = mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
            c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
            e1, e2, std::ptr::null_mut(), std::ptr::null_mut(),
            group1.len() as c_int, group2.len() as c_int,
            alloclen as c_int, &mut fftlog as *mut c_int,
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
        let c_aln_len = {
            let s = c_s1_ptrs[0];
            let mut n = 0; while *s.add(n) != 0 { n += 1; } n
        };
        eprintln!("C Falign(step=6 side=0, clean): output width = {}", c_aln_len);

        // Compare C's output for group1 idx=0 to Rust's real refinement state at same idx
        let c_first = &c_s1_boxed[0][..c_aln_len];
        let rust_first = &rust_after[group1[0]];
        if c_first != rust_first.as_slice() {
            // Find first difference
            let min = c_first.len().min(rust_first.len());
            let mut first_diff = None;
            for i in 0..min {
                if c_first[i] != rust_first[i] { first_diff = Some(i); break; }
            }
            eprintln!("C vs Rust real refinement seq[group1[0]]: c_len={}, r_len={}, first_diff={:?}",
                c_first.len(), rust_first.len(), first_diff);
            if let Some(d) = first_diff {
                let s = d.saturating_sub(5);
                let e = (d + 10).min(min);
                eprintln!("  c[{s}..{e}] = {}", String::from_utf8_lossy(&c_first[s..e]));
                eprintln!("  r[{s}..{e}] = {}", String::from_utf8_lossy(&rust_first[s..e]));
            }
        } else {
            eprintln!("C Falign output matches Rust real refinement at seq[group1[0]]");
        }

        mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            0, 0, 0, std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
        mafft_sys::alignableReagion(0, 0, std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }
}

/// Simulate a full Rust iter 0 refinement while also running C's Falign at each
/// branch with the SAME evolving state, and compare sequences / accept decisions.
#[test]
fn evolving_state_rust_vs_c_falign_iter0() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let mut sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(&sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist);
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }
        std::ptr::addr_of_mut!(mafft_sys::alg).write(b'M' as i8);
        std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(1);
        std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(1);
        std::ptr::addr_of_mut!(mafft_sys::use_fft).write(1);
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(1);

        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let nsteps = topo.steps.len();
        let mut branches: Vec<(usize, usize)> = Vec::new();
        for step in 0..nsteps {
            if step == nsteps - 1 { branches.push((step, 1)); }
            else { branches.push((step, 0)); branches.push((step, 1)); }
        }

        let mut divergences = 0usize;
        let mut first_diverge: Option<(usize, usize)> = None;
        for (step, side) in &branches {
            let group1: Vec<usize> = if *side == 0 { topo.steps[*step].left.clone() } else { topo.steps[*step].right.clone() };
            let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
            if group1.is_empty() || group2.is_empty() { continue; }

            let weights = bw.weights_for_branch(&topo, *step, *side);
            const MIN_W: f64 = 0.00001;
            let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rs1: f64 = rw1.iter().sum(); let rs2: f64 = rw2.iter().sum();
            let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
            let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

            let width = sequences[0].len();

            // ---- Rust segmented pipeline ----
            let full1: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
            let full2: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();
            let full_prof1 = Profile::from_aligned(&full1, &rw1n, &scoring.amino_map, scoring.nalphabets);
            let full_prof2 = Profile::from_aligned(&full2, &rw2n, &scoring.amino_map, scoring.nalphabets);
            let totaleff = rw1n.iter().sum::<f64>() * rw2n.iter().sum::<f64>();
            let len = full_prof1.length.min(full_prof2.length);
            let mut scores = vec![0.0f64; len];
            for i in 0..len {
                scores[i] = full_prof1.match_score(i, &full_prof2, i, &scoring.substitution_matrix) / totaleff;
            }
            let segs = mafft_fft::alignable_segments(&scores, &mafft_fft::SegmentParams::protein());
            let mut cuts: Vec<usize> = Vec::with_capacity(segs.len() + 2);
            cuts.push(0);
            for s in &segs { cuts.push(s.center.min(width)); }
            cuts.push(width);
            cuts.sort(); cuts.dedup();

            let mut rust_new: Vec<Vec<u8>> = vec![Vec::new(); nseq];
            for w in cuts.windows(2) {
                let (a, b) = (w[0], w[1]);
                if a >= b { continue; }
                let seg1: Vec<Vec<u8>> = group1.iter().map(|&i| sequences[i][a..b].to_vec()).collect();
                let seg2: Vec<Vec<u8>> = group2.iter().map(|&i| sequences[i][a..b].to_vec()).collect();
                let sw = b - a;
                let gap1: Vec<bool> = (0..sw).map(|c| seg1.iter().all(|s| s[c] == b'-')).collect();
                let gap2: Vec<bool> = (0..sw).map(|c| seg2.iter().all(|s| s[c] == b'-')).collect();
                let kept1: Vec<usize> = (0..sw).filter(|&c| !gap1[c]).collect();
                let kept2: Vec<usize> = (0..sw).filter(|&c| !gap2[c]).collect();
                let s1: Vec<Vec<u8>> = seg1.iter().map(|s| kept1.iter().map(|&c| s[c]).collect()).collect();
                let s2: Vec<Vec<u8>> = seg2.iter().map(|s| kept2.iter().map(|&c| s[c]).collect()).collect();
                if s1.is_empty() || s1[0].is_empty() || s2.is_empty() || s2[0].is_empty() { continue; }
                let r1: Vec<&[u8]> = s1.iter().map(|s| s.as_slice()).collect();
                let r2: Vec<&[u8]> = s2.iter().map(|s| s.as_slice()).collect();
                let p1 = Profile::from_aligned(&r1, &rw1n, &scoring.amino_map, scoring.nalphabets);
                let p2 = Profile::from_aligned(&r2, &rw2n, &scoring.amino_map, scoring.nalphabets);
                let aln = profile_align(&p1, &p2, &scoring.substitution_matrix, &gap, true, true);
                let mut i1 = 0usize; let mut i2 = 0usize;
                for op in &aln.operations {
                    match op {
                        mafft_align::AlignOp::Match => {
                            for (k, &idx) in group1.iter().enumerate() { rust_new[idx].push(s1[k][i1]); }
                            for (k, &idx) in group2.iter().enumerate() { rust_new[idx].push(s2[k][i2]); }
                            i1 += 1; i2 += 1;
                        }
                        mafft_align::AlignOp::Delete => {
                            for (k, &idx) in group1.iter().enumerate() { rust_new[idx].push(s1[k][i1]); }
                            for &idx in &group2 { rust_new[idx].push(b'-'); }
                            i1 += 1;
                        }
                        mafft_align::AlignOp::Insert => {
                            for &idx in &group1 { rust_new[idx].push(b'-'); }
                            for (k, &idx) in group2.iter().enumerate() { rust_new[idx].push(s2[k][i2]); }
                            i2 += 1;
                        }
                    }
                }
            }

            // ---- C Falign on same inputs ----
            let mut fftlog: c_int = 0;
            let alloclen = width * 3;
            let c_s1_boxed: Vec<Box<[u8]>> = group1.iter().map(|&i| {
                let mut v = sequences[i].clone(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let c_s2_boxed: Vec<Box<[u8]>> = group2.iter().map(|&i| {
                let mut v = sequences[i].clone(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
            let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }
            let _ = mafft_sys::Falign(
                std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
                c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
                e1, e2, std::ptr::null_mut(), std::ptr::null_mut(),
                group1.len() as c_int, group2.len() as c_int,
                alloclen as c_int, &mut fftlog as *mut c_int,
                std::ptr::null_mut(), 0, std::ptr::null_mut(),
            );
            let c_aln_len = {
                let s = c_s1_ptrs[0];
                let mut n = 0; while *s.add(n) != 0 { n += 1; } n
            };

            let rust_len = rust_new[group1[0]].len();
            let mut byte_mismatches = 0usize;
            if rust_len == c_aln_len {
                for (k, &idx) in group1.iter().enumerate() {
                    if rust_new[idx].as_slice() != &c_s1_boxed[k][..c_aln_len] { byte_mismatches += 1; }
                }
                for (k, &idx) in group2.iter().enumerate() {
                    if rust_new[idx].as_slice() != &c_s2_boxed[k][..c_aln_len] { byte_mismatches += 1; }
                }
            }
            if rust_len != c_aln_len || byte_mismatches > 0 {
                divergences += 1;
                if first_diverge.is_none() { first_diverge = Some((*step, *side)); }
                if divergences <= 3 {
                    eprintln!("step={} side={} DIFF: rust_len={} c_len={} byte_mm={}",
                        step, side, rust_len, c_aln_len, byte_mismatches);
                }
            }

            // Check if C's output equals input (C's identity check on representatives).
            // `OneClusterAndTheOther_fast` picks s1, s2 as the first members of each group.
            let s1_idx = group1[0];
            let s2_idx = group2[0];
            let c_s1_out = &c_s1_boxed[0][..c_aln_len];
            let c_s2_out = &c_s2_boxed[0][..c_aln_len];
            let c_identity = c_aln_len == width
                && c_s1_out == sequences[s1_idx].as_slice()
                && c_s2_out == sequences[s2_idx].as_slice();

            // Apply C's accept decision to evolve state.
            let old_score = compute_intergroup_score_rust(&group1, &group2, &sequences, &rw1n, &rw2n, &scoring);
            // Reconstruct C's output into a full new_sequences to compute tscore
            let mut c_new: Vec<Vec<u8>> = vec![Vec::new(); nseq];
            for (k, &idx) in group1.iter().enumerate() { c_new[idx] = c_s1_boxed[k][..c_aln_len].to_vec(); }
            for (k, &idx) in group2.iter().enumerate() { c_new[idx] = c_s2_boxed[k][..c_aln_len].to_vec(); }
            let tscore = compute_intergroup_score_rust(&group1, &group2, &c_new, &rw1n, &rw2n, &scoring);

            let rust_changed_all = (0..nseq).any(|i| sequences[i] != rust_new[i]);
            let rust_changed_reps = sequences[s1_idx] != rust_new[s1_idx] || sequences[s2_idx] != rust_new[s2_idx];
            if *step == 26 && *side == 1 {
                eprintln!("step=26 side=1: width={} c_len={} c_identity={} rust_changed_all={} rust_changed_reps={} old={:.3} tscore={:.3} delta={:.3}",
                    width, c_aln_len, c_identity, rust_changed_all, rust_changed_reps, old_score, tscore, tscore - old_score);
            }
            if !c_identity && tscore > old_score {
                sequences = c_new;
            }
        }
        eprintln!("Evolving-state Rust vs C Falign divergences (iter 0): {}/{}", divergences, branches.len());
        if let Some((s, sd)) = first_diverge {
            eprintln!("First divergence at step={} side={}", s, sd);
        }

        mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            0, 0, 0, std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
        mafft_sys::alignableReagion(0, 0, std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }
}

/// Compare our per-segment Falign-equivalent pipeline output against C's Falign
/// output at every branch of iter 0. Uses MSalignmm per segment (no sgap/egap
/// — same as our Rust code) and concatenates.
#[test]
fn compare_segmented_falign_with_c_iter0() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(&sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist);
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }
        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let nsteps = topo.steps.len();
        let mut branches: Vec<(usize, usize)> = Vec::new();
        for step in 0..nsteps {
            if step == nsteps - 1 { branches.push((step, 1)); }
            else { branches.push((step, 0)); branches.push((step, 1)); }
        }

        let width = sequences[0].len();
        let mut divergences = 0usize;
        for (step, side) in &branches {
            let group1: Vec<usize> = if *side == 0 { topo.steps[*step].left.clone() } else { topo.steps[*step].right.clone() };
            let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
            if group1.is_empty() || group2.is_empty() { continue; }

            let weights = bw.weights_for_branch(&topo, *step, *side);
            const MIN_W: f64 = 0.00001;
            let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rs1: f64 = rw1.iter().sum(); let rs2: f64 = rw2.iter().sum();
            let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
            let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

            // --- Run our segmented Rust pipeline (mirroring realign_all use_fft path) ---
            let full1: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
            let full2: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();
            let full_prof1 = Profile::from_aligned(&full1, &rw1n, &scoring.amino_map, scoring.nalphabets);
            let full_prof2 = Profile::from_aligned(&full2, &rw2n, &scoring.amino_map, scoring.nalphabets);
            let totaleff = rw1n.iter().sum::<f64>() * rw2n.iter().sum::<f64>();
            let len = full_prof1.length.min(full_prof2.length);
            let mut scores = vec![0.0f64; len];
            for i in 0..len {
                scores[i] = full_prof1.match_score(i, &full_prof2, i, &scoring.substitution_matrix) / totaleff;
            }
            let segs = mafft_fft::alignable_segments(&scores, &mafft_fft::SegmentParams::protein());
            let mut cuts: Vec<usize> = Vec::with_capacity(segs.len() + 2);
            cuts.push(0);
            for s in &segs { cuts.push(s.center.min(width)); }
            cuts.push(width);
            cuts.sort(); cuts.dedup();

            let mut rust_new: Vec<Vec<u8>> = vec![Vec::new(); nseq];
            for w in cuts.windows(2) {
                let (a, b) = (w[0], w[1]);
                if a >= b { continue; }
                let seg1: Vec<Vec<u8>> = group1.iter().map(|&i| sequences[i][a..b].to_vec()).collect();
                let seg2: Vec<Vec<u8>> = group2.iter().map(|&i| sequences[i][a..b].to_vec()).collect();
                let sw = b - a;
                let gap1: Vec<bool> = (0..sw).map(|c| seg1.iter().all(|s| s[c] == b'-')).collect();
                let gap2: Vec<bool> = (0..sw).map(|c| seg2.iter().all(|s| s[c] == b'-')).collect();
                let kept1: Vec<usize> = (0..sw).filter(|&c| !gap1[c]).collect();
                let kept2: Vec<usize> = (0..sw).filter(|&c| !gap2[c]).collect();
                let s1: Vec<Vec<u8>> = seg1.iter().map(|s| kept1.iter().map(|&c| s[c]).collect()).collect();
                let s2: Vec<Vec<u8>> = seg2.iter().map(|s| kept2.iter().map(|&c| s[c]).collect()).collect();
                if s1.is_empty() || s1[0].is_empty() || s2.is_empty() || s2[0].is_empty() { continue; }
                let r1: Vec<&[u8]> = s1.iter().map(|s| s.as_slice()).collect();
                let r2: Vec<&[u8]> = s2.iter().map(|s| s.as_slice()).collect();
                let p1 = Profile::from_aligned(&r1, &rw1n, &scoring.amino_map, scoring.nalphabets);
                let p2 = Profile::from_aligned(&r2, &rw2n, &scoring.amino_map, scoring.nalphabets);
                let aln = profile_align(&p1, &p2, &scoring.substitution_matrix, &gap, true, true);
                let mut i1 = 0usize; let mut i2 = 0usize;
                for op in &aln.operations {
                    match op {
                        mafft_align::AlignOp::Match => {
                            for (k, &idx) in group1.iter().enumerate() { rust_new[idx].push(s1[k][i1]); }
                            for (k, &idx) in group2.iter().enumerate() { rust_new[idx].push(s2[k][i2]); }
                            i1 += 1; i2 += 1;
                        }
                        mafft_align::AlignOp::Delete => {
                            for (k, &idx) in group1.iter().enumerate() { rust_new[idx].push(s1[k][i1]); }
                            for &idx in &group2 { rust_new[idx].push(b'-'); }
                            i1 += 1;
                        }
                        mafft_align::AlignOp::Insert => {
                            for &idx in &group1 { rust_new[idx].push(b'-'); }
                            for (k, &idx) in group2.iter().enumerate() { rust_new[idx].push(s2[k][i2]); }
                            i2 += 1;
                        }
                    }
                }
            }

            // --- Run C's Falign on the same inputs ---
            let mut fftlog: c_int = 0;
            let alloclen = width * 3;
            let c_s1_boxed: Vec<Box<[u8]>> = group1.iter().map(|&i| {
                let mut v = sequences[i].clone();
                v.resize(alloclen + 1, 0);
                v.into_boxed_slice()
            }).collect();
            let c_s2_boxed: Vec<Box<[u8]>> = group2.iter().map(|&i| {
                let mut v = sequences[i].clone();
                v.resize(alloclen + 1, 0);
                v.into_boxed_slice()
            }).collect();
            let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();

            let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
            let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }

            // Ensure C's alg is set to 'M' (MSalignmm) and fftkeika=1 so Falign
            // uses the commongappick-per-segment path we want.
            std::ptr::addr_of_mut!(mafft_sys::alg).write(b'M' as i8);
            std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(1);
            std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(1);
            std::ptr::addr_of_mut!(mafft_sys::use_fft).write(1);
            std::ptr::addr_of_mut!(mafft_sys::outgap).write(1);

            let _ = mafft_sys::Falign(
                std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
                c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
                e1, e2, std::ptr::null_mut(), std::ptr::null_mut(),
                group1.len() as c_int, group2.len() as c_int,
                alloclen as c_int, &mut fftlog as *mut c_int,
                std::ptr::null_mut(), 0, std::ptr::null_mut(),
            );

            // Read back C's output sequences.
            let c_aln_len = {
                let s = c_s1_ptrs[0];
                let mut n = 0; while *s.add(n) != 0 { n += 1; } n
            };

            let rust_len = rust_new[group1[0]].len();
            if rust_len != c_aln_len {
                divergences += 1;
                if divergences <= 3 {
                    eprintln!("step={} side={} DIFF: rust_len={}, c_len={}", step, side, rust_len, c_aln_len);
                }
                continue;
            }
            let mut byte_mismatches = 0usize;
            for (k, &idx) in group1.iter().enumerate() {
                let r = &rust_new[idx];
                let c = &c_s1_boxed[k][..c_aln_len];
                if r.as_slice() != c { byte_mismatches += 1; }
            }
            for (k, &idx) in group2.iter().enumerate() {
                let r = &rust_new[idx];
                let c = &c_s2_boxed[k][..c_aln_len];
                if r.as_slice() != c { byte_mismatches += 1; }
            }
            if byte_mismatches > 0 {
                divergences += 1;
                if divergences <= 3 {
                    eprintln!("step={} side={} SAME LEN {} but {} seq rows differ",
                        step, side, rust_len, byte_mismatches);
                }
            }
        }
        eprintln!("Segmented-Falign divergences (iter 0, clean state): {}/{}", divergences, branches.len());

        mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            0, 0, 0, std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
        mafft_sys::alignableReagion(0, 0, std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }
}

/// Compare Rust `alignable_segments` output to C's `alignableReagion` at lag=0
/// for every branch of iter 0 on the clean FFT-NS-2 input.
#[test]
fn compare_alignable_regions_with_c_iter0() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(&sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist);
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();

        // fftWinSize/fftThreshold are set by constants() → defaults 20 / 80.
        // Ensure they match our Rust defaults.
        let c_fftwinsize = std::ptr::addr_of!(mafft_sys::fftWinSize).read();
        let c_fftthreshold = std::ptr::addr_of!(mafft_sys::fftThreshold).read();
        eprintln!("C fftWinSize={c_fftwinsize}, fftThreshold={c_fftthreshold}");

        let nsteps = topo.steps.len();
        let mut branches: Vec<(usize, usize)> = Vec::new();
        for step in 0..nsteps {
            if step == nsteps - 1 { branches.push((step, 1)); }
            else { branches.push((step, 0)); branches.push((step, 1)); }
        }

        let mut mismatches: Vec<(usize, usize, usize, usize)> = Vec::new();
        for (step, side) in &branches {
            let group1: Vec<usize> = if *side == 0 { topo.steps[*step].left.clone() } else { topo.steps[*step].right.clone() };
            let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
            if group1.is_empty() || group2.is_empty() { continue; }

            let weights = bw.weights_for_branch(&topo, *step, *side);
            const MIN_W: f64 = 0.00001;
            let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rs1: f64 = rw1.iter().sum();
            let rs2: f64 = rw2.iter().sum();
            let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
            let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

            // Rust side: compute site scores at lag=0 and run alignable_segments.
            let full1: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
            let full2: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();
            let prof1 = Profile::from_aligned(&full1, &rw1n, &scoring.amino_map, scoring.nalphabets);
            let prof2 = Profile::from_aligned(&full2, &rw2n, &scoring.amino_map, scoring.nalphabets);
            let totaleff = rw1n.iter().sum::<f64>() * rw2n.iter().sum::<f64>();
            let len = prof1.length.min(prof2.length);
            let mut scores = vec![0.0f64; len];
            for i in 0..len {
                scores[i] = prof1.match_score(i, &prof2, i, &scoring.substitution_matrix) / totaleff;
            }
            let rust_segs = mafft_fft::alignable_segments(&scores, &mafft_fft::SegmentParams::protein());
            let rust_centers: Vec<usize> = rust_segs.iter().map(|s| s.center).collect();

            // C side: call alignableReagion with the same inputs at lag=0.
            let g1_seqs: Vec<CString> = group1.iter().map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();
            let g2_seqs: Vec<CString> = group2.iter().map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();
            let mut g1_ptrs: Vec<*mut c_char> = g1_seqs.iter().map(|c| c.as_ptr() as *mut c_char).collect();
            let mut g2_ptrs: Vec<*mut c_char> = g2_seqs.iter().map(|c| c.as_ptr() as *mut c_char).collect();
            let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
            let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }

            // Allocate MAXSEG Segment slots (ballpark — 1000 is plenty)
            const MAX_SEG: usize = 1000;
            let seg_buf: *mut mafft_sys::Segment = alloc_zeroed(MAX_SEG * std::mem::size_of::<mafft_sys::Segment>()) as _;
            let c_count = mafft_sys::alignableReagion(
                group1.len() as c_int,
                group2.len() as c_int,
                g1_ptrs.as_mut_ptr(),
                g2_ptrs.as_mut_ptr(),
                e1,
                e2,
                seg_buf,
            ) as usize;
            let c_centers: Vec<usize> = (0..c_count).map(|i| (*seg_buf.add(i)).center as usize).collect();

            if rust_centers != c_centers {
                mismatches.push((*step, *side, rust_centers.len(), c_count));
                if mismatches.len() <= 5 {
                    eprintln!("step={} side={}: rust centers={:?}, c centers={:?}", step, side, rust_centers, c_centers);
                }
            }
        }
        eprintln!("Segment mismatch branches (iter 0, clean state): {}/{}", mismatches.len(), branches.len());

        // Deallocate C state
        mafft_sys::alignableReagion(0, 0, std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }
}

/// Compare Rust profile_align vs C MSalignmm at EVERY branch of iter 0.
/// Reports the first branch where the produced alignments differ.
#[test]
fn full_iter0_profile_align_vs_msalignmm() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample.fftns2"))
        .expect("load sample.fftns2");
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let nseq = input.sequences.len();
    let mut sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let penalty_dist = scoring.gap.open;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = scoring_matrix_distance(
                &sequences[i], &sequences[j],
                &scoring.substitution_matrix, &scoring.amino_map, penalty_dist,
            );
            dm.set(i, j, d);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let bw = BranchWeights::new(&topo);

    unsafe {
        init_c_protein();
        let topol = build_c_topol(&topo);
        let len = build_c_len(&topo);
        let bw_mtx: *mut *mut c_double = alloc_zeroed(topo.steps.len() * std::mem::size_of::<*mut c_double>()) as _;
        for k in 0..topo.steps.len() {
            let row: *mut c_double = alloc_zeroed(2 * std::mem::size_of::<c_double>()) as _;
            *row.add(0) = 1.0;
            *row.add(1) = 1.0;
            *bw_mtx.add(k) = row;
        }
        std::ptr::addr_of_mut!(mafft_sys::sueff_global).write(0.1);
        std::ptr::addr_of_mut!(mafft_sys::treemethod).write(b'X' as c_int);
        let stopol: *mut mafft_sys::Node = alloc_zeroed(2 * nseq * std::mem::size_of::<mafft_sys::Node>()) as _;
        mafft_sys::treeCnv(stopol, nseq as c_int, topol, len, bw_mtx);
        mafft_sys::calcBranchWeight(bw_mtx, nseq as c_int, stopol, topol, len);

        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }

        let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
        let nsteps = topo.steps.len();

        // Enumerate branches in iter 0 order (forward, with root skip)
        let mut branches: Vec<(usize, usize)> = Vec::new();
        for step in 0..nsteps {
            if step == nsteps - 1 {
                branches.push((step, 1));
            } else {
                branches.push((step, 0));
                branches.push((step, 1));
            }
        }

        let mut divergences: Vec<(usize, usize, String)> = Vec::new();
        for (step, side) in &branches {
            let group1: Vec<usize> = if *side == 0 {
                topo.steps[*step].left.clone()
            } else {
                topo.steps[*step].right.clone()
            };
            let group2: Vec<usize> = (0..nseq).filter(|i| !group1.contains(i)).collect();
            if group1.is_empty() || group2.is_empty() { continue; }

            // Rust weights
            let weights = bw.weights_for_branch(&topo, *step, *side);
            const MIN_W: f64 = 0.00001;
            let rw1: Vec<f64> = group1.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rw2: Vec<f64> = group2.iter().map(|&i| weights[i].max(MIN_W)).collect();
            let rs1: f64 = rw1.iter().sum();
            let rs2: f64 = rw2.iter().sum();
            let rw1n: Vec<f64> = rw1.iter().map(|w| w / rs1).collect();
            let rw2n: Vec<f64> = rw2.iter().map(|w| w / rs2).collect();

            // Strip per-group gap columns
            let width = sequences[0].len();
            let gap1: Vec<bool> = (0..width).map(|c| group1.iter().all(|&i| sequences[i][c] == b'-')).collect();
            let gap2: Vec<bool> = (0..width).map(|c| group2.iter().all(|&i| sequences[i][c] == b'-')).collect();
            let kept1: Vec<usize> = (0..width).filter(|&c| !gap1[c]).collect();
            let kept2: Vec<usize> = (0..width).filter(|&c| !gap2[c]).collect();

            let stripped1: Vec<Vec<u8>> = group1.iter().map(|&i| kept1.iter().map(|&c| sequences[i][c]).collect()).collect();
            let stripped2: Vec<Vec<u8>> = group2.iter().map(|&i| kept2.iter().map(|&c| sequences[i][c]).collect()).collect();
            if stripped1.is_empty() || stripped1[0].is_empty() || stripped2.is_empty() || stripped2[0].is_empty() {
                continue;
            }
            let s1r: Vec<&[u8]> = stripped1.iter().map(|s| s.as_slice()).collect();
            let s2r: Vec<&[u8]> = stripped2.iter().map(|s| s.as_slice()).collect();

            // Rust profile_align
            let prof1 = Profile::from_aligned(&s1r, &rw1n, &scoring.amino_map, scoring.nalphabets);
            let prof2 = Profile::from_aligned(&s2r, &rw2n, &scoring.amino_map, scoring.nalphabets);
            let rust_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

            // Reconstruct rust aligned
            let mut r_g1: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln.operations.len()); group1.len()];
            let mut r_g2: Vec<Vec<u8>> = vec![Vec::with_capacity(rust_aln.operations.len()); group2.len()];
            let mut i1 = 0usize; let mut i2 = 0usize;
            for op in &rust_aln.operations {
                match op {
                    mafft_align::AlignOp::Match => {
                        for (k, s) in stripped1.iter().enumerate() { r_g1[k].push(s[i1]); }
                        for (k, s) in stripped2.iter().enumerate() { r_g2[k].push(s[i2]); }
                        i1 += 1; i2 += 1;
                    }
                    mafft_align::AlignOp::Insert => {
                        for r in r_g1.iter_mut() { r.push(b'-'); }
                        for (k, s) in stripped2.iter().enumerate() { r_g2[k].push(s[i2]); }
                        i2 += 1;
                    }
                    mafft_align::AlignOp::Delete => {
                        for (k, s) in stripped1.iter().enumerate() { r_g1[k].push(s[i1]); }
                        for r in r_g2.iter_mut() { r.push(b'-'); }
                        i1 += 1;
                    }
                }
            }

            // C MSalignmm
            let alloclen = width * 10;
            let c_s1_boxed: Vec<Box<[u8]>> = stripped1.iter().map(|s| {
                let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let c_s2_boxed: Vec<Box<[u8]>> = stripped2.iter().map(|s| {
                let mut v = s.to_vec(); v.resize(alloclen + 1, 0); v.into_boxed_slice()
            }).collect();
            let mut c_s1_ptrs: Vec<*mut c_char> = c_s1_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();
            let mut c_s2_ptrs: Vec<*mut c_char> = c_s2_boxed.iter().map(|v| v.as_ptr() as *mut c_char).collect();

            let e1: *mut c_double = alloc_zeroed(group1.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw1n.iter().enumerate() { *e1.add(i) = w; }
            let e2: *mut c_double = alloc_zeroed(group2.len() * std::mem::size_of::<c_double>()) as _;
            for (i, &w) in rw2n.iter().enumerate() { *e2.add(i) = w; }

            let _c_score = mafft_sys::MSalignmm(
                n_dyn,
                c_s1_ptrs.as_mut_ptr(), c_s2_ptrs.as_mut_ptr(),
                e1, e2,
                group1.len() as c_int, group2.len() as c_int,
                alloclen as c_int,
                std::ptr::null_mut(), std::ptr::null_mut(),
                std::ptr::null_mut(), std::ptr::null_mut(),
                std::ptr::null_mut(), 0, std::ptr::null_mut(),
                1, 1,
                std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
                1.0, 1.0,
            );

            let c_aln_len = {
                let s = c_s1_ptrs[0];
                let mut n = 0;
                while *s.add(n) != 0 { n += 1; }
                n
            };

            // Compare
            let rust_len = rust_aln.operations.len();
            if rust_len != c_aln_len {
                divergences.push((*step, *side, format!("len: rust={} c={}", rust_len, c_aln_len)));
                if divergences.len() <= 5 {
                    eprintln!("DIFF step={} side={} rust_len={} c_len={}", step, side, rust_len, c_aln_len);
                }
                continue;
            }
            let mut byte_diff = 0;
            for i in 0..group1.len() {
                for c in 0..rust_len {
                    if r_g1[i][c] != c_s1_boxed[i][c] { byte_diff += 1; break; }
                }
            }
            for i in 0..group2.len() {
                for c in 0..rust_len {
                    if r_g2[i][c] != c_s2_boxed[i][c] { byte_diff += 1; break; }
                }
            }
            if byte_diff > 0 {
                divergences.push((*step, *side, format!("bytes: {} seq rows differ", byte_diff)));
                if divergences.len() <= 5 {
                    eprintln!("DIFF step={} side={} (same len {}, {} seq rows differ) g1={} g2={}",
                        step, side, rust_len, byte_diff, group1.len(), group2.len());
                    eprintln!("  group1[{}]={:?}", group1.len(), &group1);
                    eprintln!("  kept1.len={} kept2.len={}", kept1.len(), kept2.len());
                    eprintln!("  rust score={:.3} c score={:.3}", rust_aln.score, _c_score);
                    // Print first row that differs
                    for i in 0..group1.len() {
                        for c in 0..rust_len {
                            if r_g1[i][c] != c_s1_boxed[i][c] {
                                let start = c.saturating_sub(8);
                                let end = (c + 8).min(rust_len);
                                eprintln!("  g1[{}] col {}: rust={} c={}",
                                    i, c,
                                    String::from_utf8_lossy(&r_g1[i][start..end]),
                                    String::from_utf8_lossy(&c_s1_boxed[i][start..end]));
                                break;
                            }
                        }
                    }
                }
            }
        }
        eprintln!("Total divergent branches: {}/{}", divergences.len(), branches.len());
        for (step, side, msg) in divergences.iter().take(10) {
            eprintln!("  step={} side={} {}", step, side, msg);
        }
        mafft_sys::freeconstants();

        // Assert that at least no divergences on non-updated state (iter 0).
        // We haven't updated `sequences` — we're testing against the initial FFT-NS-2 state.
    }
}

fn compute_intergroup_score_rust(
    group1: &[usize],
    group2: &[usize],
    sequences: &[Vec<u8>],
    w1n: &[f64],
    w2n: &[f64],
    scoring: &mafft_types::ScoringContext,
) -> f64 {
    let penalty = scoring.gap.open as f64;
    let mtx = &scoring.substitution_matrix;
    let map = &scoring.amino_map;
    let mtx_size = mtx.len();

    let pairwise = |a: &[u8], b: &[u8]| -> f64 {
        let l = a.len().min(b.len());
        let mut score = 0.0f64;
        let mut k = 0;
        while k < l {
            let ca = a[k];
            let cb = b[k];
            if ca == b'-' && cb == b'-' { k += 1; continue; }
            if ca == b'-' {
                score += penalty;
                k += 1;
                while k < l && a[k] == b'-' { k += 1; }
                continue;
            }
            if cb == b'-' {
                score += penalty;
                k += 1;
                while k < l && b[k] == b'-' { k += 1; }
                continue;
            }
            let i = map[ca as usize] as usize;
            let j = map[cb as usize] as usize;
            if i < mtx_size && j < mtx_size { score += mtx[i][j] as f64; }
            k += 1;
        }
        score
    };

    let mut total = 0.0f64;
    for (i_local, &i) in group1.iter().enumerate() {
        for (j_local, &j) in group2.iter().enumerate() {
            total += pairwise(&sequences[i], &sequences[j]) * w1n[i_local] * w2n[j_local];
        }
    }
    total
}
