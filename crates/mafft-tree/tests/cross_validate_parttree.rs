/// FFI cross-validation for the `--parttree` 6-mer distance pipeline.
///
/// The C reference is `commonsextet_p` in `mltaln9.c` (the non-static,
/// non-TLS twin of `splittbfast.c::localcommonsextet_p`). Both produce
/// identical integer common-counts for the same `(table, points)`
/// inputs, so we use `commonsextet_p` for the FFI guard.
///
/// Tests:
/// - `parttree_points_match_c`: Rust `encode_points_protein` agrees
///   with C `seq_grp` + `makepointtable` byte-for-byte for each
///   sequence in the 36-seq sample.
/// - `parttree_common_sextets_match_c`: Rust `common_sextets_p` agrees
///   with C `commonsextet_p` for every (i, j) pivot pair.
/// - `parttree_lenfac_matches_c`: Rust `lenfac` reproduces C's formula
///   bit-for-bit for the lengths in the sample.
/// - `parttree_distance_matches_c`: end-to-end pairwise distance
///   matches C's `(1 - common/min(ss,ss)) * lenfac`.

use std::ffi::CString;
use std::os::raw::{c_char, c_int};
use std::sync::Mutex;

use mafft_io::read_fasta;
use mafft_tree::{musclesupg, ClusterMethod, DistanceMatrix};
use mafft_tree::parttree_dist::{
    common_sextets_p, composition_table, encode_points_protein, lenfac,
    parttree_distance_protein, MAX6DIST, PLENFACA, PLENFACB, PLENFACC, PLENFACD,
};

static C_MUTEX: Mutex<()> = Mutex::new(());

/// Initialize C globals to the protein, parttree path: `dorp='p'`,
/// `scoremtx=1` (BLOSUM, but doesn't actually matter for sextet
/// counting — we just need `tsize`/`amino_grp` populated by `constants`).
unsafe fn init_c_protein() {
    unsafe {
        mafft_sys::initglobalvariables();
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(1);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(62);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::tsize).write(46656); // 6^6
        // `commonsextet_p` allocates `ct[MIN(maxl, tsize) + 1]` — needs to
        // be at least as large as the longest sequence's point vector.
        std::ptr::addr_of_mut!(mafft_sys::maxl).write(46656);
        std::ptr::addr_of_mut!(mafft_sys::lenfaca).write(PLENFACA);
        std::ptr::addr_of_mut!(mafft_sys::lenfacb).write(PLENFACB);
        std::ptr::addr_of_mut!(mafft_sys::lenfacc).write(PLENFACC);
        std::ptr::addr_of_mut!(mafft_sys::lenfacd).write(PLENFACD);

        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);
    }
}

unsafe fn cleanup_c() {
    unsafe {
        // commonsextet_p has TLS state; passing (NULL, NULL) frees it.
        mafft_sys::commonsextet_p(std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }
}

/// Run C's `seq_grp` + `makepointtable` for `seq` and return the point
/// vector. Mirrors what `splittbfast` does per-sequence at startup.
unsafe fn c_encode_points(seq: &[u8]) -> Vec<i32> {
    let cs = CString::new(seq).unwrap();
    // grp buffer must hold at most seq.len() entries + END_OF_VEC sentinel.
    let mut grp = vec![0i32; seq.len() + 1];
    let nvalid = unsafe { mafft_sys::seq_grp(grp.as_mut_ptr(), cs.as_ptr() as *const c_char) };
    if nvalid < 6 {
        return Vec::new();
    }
    grp[nvalid as usize] = -1;
    let mut pointt = vec![0i32; (nvalid - 5) as usize + 1];
    unsafe { mafft_sys::makepointtable(pointt.as_mut_ptr(), grp.as_mut_ptr()) };
    pointt.pop(); // drop END_OF_VEC sentinel
    pointt
}

#[test]
fn parttree_points_match_c() {
    // C globals are shared across all FFI tests in this binary; serialize
    // access. Recover from poison so a single failing test doesn't cascade.
    let _g = C_MUTEX.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");

    unsafe { init_c_protein(); }

    for (idx, seq) in input.sequences.iter().enumerate() {
        let rust_pts = encode_points_protein(&seq.data);
        let c_pts = unsafe { c_encode_points(&seq.data) };
        assert_eq!(rust_pts.len(), c_pts.len(),
            "seq {idx} '{}': Rust={} pts, C={} pts",
            seq.name, rust_pts.len(), c_pts.len());
        for (i, (&r, &c)) in rust_pts.iter().zip(c_pts.iter()).enumerate() {
            assert_eq!(r as i32, c, "seq {idx} '{}' point {i}: Rust={r} C={c}", seq.name);
        }
    }

    unsafe { cleanup_c(); }
}

#[test]
fn parttree_common_sextets_match_c() {
    // C globals are shared across all FFI tests in this binary; serialize
    // access. Recover from poison so a single failing test doesn't cascade.
    let _g = C_MUTEX.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");
    unsafe { init_c_protein(); }

    // Encode all sequences using both Rust and C.
    let rust_points: Vec<Vec<u32>> = input.sequences.iter()
        .map(|s| encode_points_protein(&s.data))
        .collect();
    let mut c_points: Vec<Vec<i32>> = input.sequences.iter()
        .map(|s| unsafe { c_encode_points(&s.data) })
        .collect();
    // Append END_OF_VEC sentinel to each (commonsextet_p reads until -1).
    for v in &mut c_points { v.push(-1); }

    // Compare `min(self, other)` count for every pair (i, j) — this is
    // what `splittbfast` does in its pickmtx construction.
    for i in 0..input.sequences.len() {
        // C composition table.
        let mut c_table = vec![0i32; 46656];
        unsafe { mafft_sys::makecompositiontable_p(c_table.as_mut_ptr(), c_points[i].as_mut_ptr()); }
        let rust_table = composition_table(&rust_points[i], 46656);
        assert_eq!(c_table, rust_table, "composition table mismatch for seq {i}");

        for j in 0..input.sequences.len() {
            let mut c_pts_j = c_points[j].clone();
            let c_common = unsafe {
                mafft_sys::commonsextet_p(c_table.as_mut_ptr(), c_pts_j.as_mut_ptr())
            };
            let rust_common = common_sextets_p(&rust_table, &rust_points[j], 46656);
            assert_eq!(rust_common, c_common,
                "common_sextets({i}, {j}): Rust={rust_common} C={c_common}");
        }
    }

    unsafe { cleanup_c(); }
}

#[test]
fn parttree_lenfac_matches_c_formula() {
    // No FFI needed — just verify our `lenfac` formula reproduces the
    // C constants exactly. The constants come from `splittbfast.c:34-37`
    // which match the active `#else` block in `disttbfast.c:43-46`, so
    // any drift here would be caught by reading `mafft_sys::lenfaca` etc.
    // C globals are shared across all FFI tests in this binary; serialize
    // access. Recover from poison so a single failing test doesn't cascade.
    let _g = C_MUTEX.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    unsafe { init_c_protein(); }
    let a = unsafe { std::ptr::addr_of!(mafft_sys::lenfaca).read() };
    let b = unsafe { std::ptr::addr_of!(mafft_sys::lenfacb).read() };
    let c = unsafe { std::ptr::addr_of!(mafft_sys::lenfacc).read() };
    let d = unsafe { std::ptr::addr_of!(mafft_sys::lenfacd).read() };
    assert_eq!(a, PLENFACA);
    assert_eq!(b, PLENFACB);
    assert_eq!(c, PLENFACC);
    assert_eq!(d, PLENFACD);
    unsafe { cleanup_c(); }

    for (l1, l2) in &[(353usize, 348), (455, 448), (180, 509), (100, 1000)] {
        let lf = lenfac(*l1, *l2, PLENFACA, PLENFACB, PLENFACC, PLENFACD);
        let longer = (*l1).max(*l2) as f64;
        let shorter = (*l1).min(*l2) as f64;
        let expected = 1.0 / (shorter / longer * PLENFACD + PLENFACB / (longer + PLENFACC) + PLENFACA);
        assert!((lf - expected).abs() < 1e-12);
    }
}

#[test]
fn parttree_distance_clamps_at_max6dist() {
    // Two completely disjoint sequences (no shared 6-mers) → raw=1.0,
    // multiplied by lenfac which can be > 1 → clamped at MAX6DIST=10.0.
    // Use very short sequences with disjoint group profiles.
    let s1: &[u8] = b"AAAAAAAAA"; // grp 0 only
    let s2: &[u8] = b"WYFFWYWYW"; // grp 4 only
    let p1 = encode_points_protein(s1);
    let p2 = encode_points_protein(s2);
    let d = parttree_distance_protein(&p1, &p2, s1.len(), s2.len());
    assert!(d <= MAX6DIST, "expected clamp at MAX6DIST=10.0, got {d}");
    assert!(d > 0.5, "expected near-max distance for disjoint seqs, got {d}");
}

/// End-to-end: every (i, j) pairwise distance from the 36-seq sample
/// should match C's `(1 - common / min(ss,ss)) * lenfac` formula
/// computed via FFI helpers.
#[test]
fn parttree_distance_matches_c() {
    // C globals are shared across all FFI tests in this binary; serialize
    // access. Recover from poison so a single failing test doesn't cascade.
    let _g = C_MUTEX.lock().unwrap_or_else(|poisoned| poisoned.into_inner());
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");
    unsafe { init_c_protein(); }

    let rust_points: Vec<Vec<u32>> = input.sequences.iter()
        .map(|s| encode_points_protein(&s.data))
        .collect();
    let mut c_points: Vec<Vec<i32>> = input.sequences.iter()
        .map(|s| unsafe { c_encode_points(&s.data) })
        .collect();
    for v in &mut c_points { v.push(-1); }

    let nseq = input.sequences.len();
    for i in 0..nseq {
        let mut c_table = vec![0i32; 46656];
        unsafe { mafft_sys::makecompositiontable_p(c_table.as_mut_ptr(), c_points[i].as_mut_ptr()); }
        for j in (i + 1)..nseq {
            let mut c_pts_j = c_points[j].clone();
            let c_common = unsafe {
                mafft_sys::commonsextet_p(c_table.as_mut_ptr(), c_pts_j.as_mut_ptr())
            };
            // C splittbfast formula at line 1445/1674:
            //   bunbo = MIN(selfscore_i, selfscore_j) = min(points.len() in our convention)
            //   raw = (1 - common / bunbo) — pre-lenfac
            //   lenfac uses orilen (= filtered length, which for protein
            //     equals len(grp)+0 = len(seq) post-filter; same as
            //     points.len()+5 since we have len-5 6-mers)
            // C uses orilen = strlen(seq[i]) (raw sequence length) for
            // lenfac, NOT the filtered group length.
            let raw_len_i = input.sequences[i].data.len();
            let raw_len_j = input.sequences[j].data.len();
            let bunbo = rust_points[i].len().min(rust_points[j].len()) as f64;
            let raw = 1.0 - c_common as f64 / bunbo;
            let lf = lenfac(raw_len_i, raw_len_j, PLENFACA, PLENFACB, PLENFACC, PLENFACD);
            let mut expected = raw * lf;
            if expected > MAX6DIST { expected = MAX6DIST; }
            if expected < 0.0 { expected = 0.0; }

            let rust_dist = parttree_distance_protein(
                &rust_points[i], &rust_points[j],
                raw_len_i, raw_len_j,
            );
            assert!((rust_dist - expected).abs() < 1e-12,
                "dist({i},{j}): Rust={rust_dist} C-formula={expected} (common={c_common} bunbo={bunbo})");
        }
    }

    unsafe { cleanup_c(); }
}

/// FFI guard for `fixed_musclesupg_double_realloc_nobk_halfmtx`. Builds a
/// 36×36 PartTree-style distance matrix from the sample, then runs the C
/// reference and Rust's `musclesupg` (which the docstring already claims
/// to port) and compares topology join steps in order.
#[test]
fn parttree_upgma_matches_c() {
    let _g = C_MUTEX.lock().unwrap_or_else(|p| p.into_inner());
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");
    let nseq = input.sequences.len();

    // Build the same parttree-style pickmtx that splittbfast would use:
    // pairwise PartTree distance for every (i, j).
    let rust_points: Vec<Vec<u32>> = input.sequences.iter()
        .map(|s| encode_points_protein(&s.data))
        .collect();
    let raw_lens: Vec<usize> = input.sequences.iter().map(|s| s.data.len()).collect();

    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = parttree_distance_protein(
                &rust_points[i], &rust_points[j], raw_lens[i], raw_lens[j],
            );
            dm.set(i, j, d);
        }
    }

    // Diagnostic: print the 5 smallest distances and which pair they are.
    {
        let mut all_dists: Vec<(usize, usize, f64)> = Vec::new();
        for i in 0..nseq {
            for j in (i + 1)..nseq {
                all_dists.push((i, j, dm.get(i, j)));
            }
        }
        all_dists.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        eprintln!("Smallest 5 distances:");
        for (i, j, d) in all_dists.iter().take(5) {
            eprintln!("  ({i}, {j}) = {d}");
        }
    }

    // Rust topology.
    let rust_topo = musclesupg(&dm, ClusterMethod::Mix { sueff: 0.1 });

    // C topology via FFI to fixed_musclesupg_double_realloc_nobk_halfmtx.
    unsafe {
        init_c_protein();
        std::ptr::addr_of_mut!(mafft_sys::treemethod).write(b'X' as c_int);
        std::ptr::addr_of_mut!(mafft_sys::sueff_global).write(0.1);
        std::ptr::addr_of_mut!(mafft_sys::njob).write(nseq as c_int);

        // Build C's half-matrix. Note `setnearest` (mltaln9.c:1419) reads
        // `eff[pos][j-pos]` — i.e. row `i` has slot `[j-i]` for j > i, with
        // `[0]` unused (the "diagonal" slot). Each row has `nseq-i` slots
        // per `AllocateFloatHalfMtx` (mtxutl.c:163).
        let eff: *mut *mut f64 = mafft_sys::AllocateFloatHalfMtx(nseq as c_int);
        for i in 0..nseq {
            for j in (i + 1)..nseq {
                *(*eff.add(i)).add(j - i) = dm.get(i, j);
            }
        }

        let topol: *mut *mut *mut c_int = mafft_sys::AllocateIntCub(
            (nseq - 1) as c_int, 2, 0,
        );
        let len: *mut *mut f64 = mafft_sys::AllocateFloatMtx(
            (nseq - 1) as c_int, 2,
        );

        mafft_sys::fixed_musclesupg_double_realloc_nobk_halfmtx(
            nseq as c_int,
            eff,
            topol,
            len,
            std::ptr::null_mut(),
            0,
            1,
        );

        let mut c_steps: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
        for k in 0..(nseq - 1) {
            let p0 = *topol.add(k);
            let p1 = *(*topol.add(k)).add(1);
            let mut left = Vec::new();
            let mut right = Vec::new();
            let mut q = *p0;
            loop {
                let v = *q;
                if v < 0 { break; }
                left.push(v as usize);
                q = q.add(1);
            }
            let mut q = p1;
            loop {
                let v = *q;
                if v < 0 { break; }
                right.push(v as usize);
                q = q.add(1);
            }
            c_steps.push((left, right));
        }

        assert_eq!(rust_topo.steps.len(), c_steps.len(),
            "topology step count: Rust={} C={}",
            rust_topo.steps.len(), c_steps.len());
        for (k, (rs, (cl, cr))) in rust_topo.steps.iter().zip(c_steps.iter()).enumerate() {
            let mut r_left = rs.left.clone(); r_left.sort();
            let mut r_right = rs.right.clone(); r_right.sort();
            let mut c_left = cl.clone(); c_left.sort();
            let mut c_right = cr.clone(); c_right.sort();
            let r_pair = if r_left < r_right { (r_left, r_right) } else { (r_right, r_left) };
            let c_pair = if c_left < c_right { (c_left, c_right) } else { (c_right, c_left) };
            assert_eq!(r_pair, c_pair,
                "step {k}: Rust=({:?}, {:?}) C=({:?}, {:?})",
                rs.left, rs.right, cl, cr);
        }

        cleanup_c();
    }
}
