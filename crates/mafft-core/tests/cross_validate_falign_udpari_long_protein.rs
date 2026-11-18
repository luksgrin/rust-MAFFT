//! Cross-validate the protein merges of `disttbfast` after `alg='M'`
//! against C: `Falign_udpari_long` (`Falign.c:2583`) for the FFT branch,
//! and `MSalignmm` with `treebase`'s profile cache (`cpmxhist`,
//! `disttbfast.c:2889`) for the DP branch.
//!
//! The DP branch is where the 477-against-12 merge of the 500-record
//! PF00069 set diverged: `MSalignmm` stores the merged profile it blended
//! from the children (`createcpmxresult`, `MSalignmm.c:283`) for the parent
//! merge, and blends the cached child's opening/closing gap counts from its
//! zero-filled scratch buffers (`MSalignmm.c:2459-2460`). The synthetic
//! groups below reproduce that chain: a 21-row group with a subfamily-wide
//! deletion (so the counts are far from zero) merged with one sequence,
//! then that cached result merged with another.
//!
//! `UDPARI_CASE_DIR` points the last test at a directory with `g1.txt`,
//! `g2.txt` (aligned rows) and `w1.txt`, `w2.txt` (raw weights) from the
//! engine's merge dump, for investigating a real merge: `UDPARI_DP=1`
//! compares the DP branch, `UDPARI_CHAIN=<prefix,...>` the cache chain,
//! otherwise the FFT branch.

#![allow(non_snake_case)]

use std::os::raw::{c_char, c_double, c_int};
use std::path::PathBuf;
use std::sync::Mutex;

use mafft_align::{fft_profile_align, fft_total_eff, AlignOp, FftAlignParams, GapModel, Profile};
use mafft_core::progressive::{blend_profiles_exact, blend_profiles_with_counts};
use mafft_types::{ScoringContext, ScoringModel, SeqType};

static C_MUTEX: Mutex<()> = Mutex::new(());

unsafe extern "C" {
    fn Falign_udpari_long(
        whichmtx: *mut *mut c_int,
        scoringmatrices: *mut *mut *mut c_double,
        n_dynamicmtx: *mut *mut c_double,
        seq1: *mut *mut c_char,
        seq2: *mut *mut c_char,
        eff1: *mut c_double,
        eff2: *mut c_double,
        eff1s: *mut *mut c_double,
        eff2s: *mut *mut c_double,
        clus1: c_int,
        clus2: c_int,
        alloclen: c_int,
        fftlog: *mut c_int,
    ) -> c_double;
}

unsafe fn init_c_protein(njob: usize, nlenmax: usize) {
    unsafe {
        mafft_sys::initglobalvariables();
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_EX).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_OP).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_dist).write(mafft_sys::NOTSPECIFIED);
        // The script passes `-h 0.000` for protein (`defaultaof`), so
        // `poffset = 0` (`disttbfast.c:260`) and the DP matrix is unshifted;
        // `NOTSPECIFIED` would shift it by the built-in default.
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(1);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(62);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        // `Falign_udpari_long` only has a `case 'M'` (`Falign.c:3350`).
        std::ptr::addr_of_mut!(mafft_sys::alg).write(b'M' as c_char);
        std::ptr::addr_of_mut!(mafft_sys::fftscore).write(1);
        std::ptr::addr_of_mut!(mafft_sys::fftWinSize).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftThreshold).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(0);
        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut c_char;
        mafft_sys::constants(1, &mut seq_ptr);
        std::ptr::addr_of_mut!(mafft_sys::njob).write(njob as c_int);
        std::ptr::addr_of_mut!(mafft_sys::nlenmax).write(nlenmax as c_int);
        // FFT-NS-2 gets `-O` from the script: outgap = 0.
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(0);
    }
}

fn scoring() -> ScoringContext {
    mafft_scoring::build_context(ScoringModel::Blosum(62), SeqType::Protein)
}

fn apply_ops(ops: &[AlignOp], g1: &[Vec<u8>], g2: &[Vec<u8>]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut r1 = vec![Vec::new(); g1.len()];
    let mut r2 = vec![Vec::new(); g2.len()];
    let (mut c1, mut c2) = (0usize, 0usize);
    for op in ops {
        let (t1, t2) = match op {
            AlignOp::Match => (true, true),
            AlignOp::Delete => (true, false),
            AlignOp::Insert => (false, true),
        };
        for k in 0..g1.len() { r1[k].push(if t1 { g1[k][c1] } else { b'-' }); }
        for k in 0..g2.len() { r2[k].push(if t2 { g2[k][c2] } else { b'-' }); }
        if t1 { c1 += 1; }
        if t2 { c2 += 1; }
    }
    (r1, r2)
}

/// The merge's gap tables (`mgt1`/`mgt2` in `MSalignmm`): `'o'` where the
/// group has a column, `'-'` where the merge inserted a gap into it.
fn gap_tables(ops: &[AlignOp]) -> (Vec<u8>, Vec<u8>) {
    let mut gt1 = Vec::with_capacity(ops.len());
    let mut gt2 = Vec::with_capacity(ops.len());
    for op in ops {
        match op {
            AlignOp::Match => { gt1.push(b'o'); gt2.push(b'o'); }
            AlignOp::Delete => { gt1.push(b'o'); gt2.push(b'-'); }
            AlignOp::Insert => { gt1.push(b'-'); gt2.push(b'o'); }
        }
    }
    (gt1, gt2)
}

/// `fastconjuction_noname` (`tddis.c:543-556`): clamp, sum in order, divide.
fn normalize(w: &[f64]) -> Vec<f64> {
    let w: Vec<f64> = w.iter().map(|&x| x.max(0.00001)).collect();
    let total: f64 = w.iter().fold(0.0, |a, b| a + b);
    w.iter().map(|x| x / total).collect()
}

fn raw_sum(w: &[f64]) -> f64 { w.iter().map(|&x| x.max(0.00001)).fold(0.0, |a, b| a + b) }

fn profile(scoring: &ScoringContext, rows: &[Vec<u8>], eff: &[f64]) -> Profile {
    let r: Vec<&[u8]> = rows.iter().map(|s| s.as_slice()).collect();
    Profile::from_aligned(&r, eff, &scoring.amino_map, scoring.nalphabets)
}

fn rust_udpari(scoring: &ScoringContext, g1: &[Vec<u8>], e1: &[f64], g2: &[Vec<u8>], e2: &[f64]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let p1 = profile(scoring, g1, e1);
    let p2 = profile(scoring, g2, e2);
    // Same property channels as `progressive::merge_step_cached`.
    let nscored = scoring.nscoredalphabets;
    let mut polarity_by_idx = vec![0.0f64; nscored];
    let mut volume_by_idx = vec![0.0f64; nscored];
    for ch in 0u16..256 {
        let idx = scoring.amino_map[ch as usize] as usize;
        if idx < nscored {
            polarity_by_idx[idx] = scoring.polarity[ch as usize];
            volume_by_idx[idx] = scoring.volume[ch as usize];
        }
    }
    let mut params = FftAlignParams::protein();
    params.gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    params.head_gap = false;
    params.tail_gap = false;
    params.num_channels = nscored;
    params.property_channels = Some((polarity_by_idx, volume_by_idx));
    params.total_eff = fft_total_eff(e1, e2);
    params.udpari_long = true;
    params.site_matrix = Some(mafft_align::matrix_i32_to_f64(&scoring.fft_matrix));
    let aln = fft_profile_align(&p1, &p2, &scoring.consweight_matrix, &params);
    apply_ops(&aln.operations, g1, g2)
}

/// NUL-padded copies of the rows and their pointer table for C.
struct CRows { bufs: Vec<Vec<u8>>, ptrs: Vec<*mut c_char> }

fn c_rows(g: &[Vec<u8>], alloclen: usize) -> CRows {
    let mut bufs: Vec<Vec<u8>> = g.iter().map(|s| { let mut v = s.clone(); v.resize(alloclen + 1, 0); v }).collect();
    let ptrs = bufs.iter_mut().map(|v| v.as_mut_ptr() as *mut c_char).collect();
    CRows { bufs, ptrs }
}

fn rows_back(c: &CRows) -> Vec<Vec<u8>> {
    c.bufs.iter().map(|v| v[..v.iter().position(|&c| c == 0).unwrap()].to_vec()).collect()
}

fn c_udpari(g1: &[Vec<u8>], e1: &[f64], g2: &[Vec<u8>], e2: &[f64]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let alloclen = (g1[0].len() + g2[0].len()) * 2 + 100;
    let mut b1 = c_rows(g1, alloclen);
    let mut b2 = c_rows(g2, alloclen);
    let mut fe1 = e1.to_vec();
    let mut fe2 = e2.to_vec();
    let mut fftlog: c_int = 0;
    unsafe {
        let n_dyn = std::ptr::addr_of!(mafft_sys::n_dis_consweight_multi).read();
        Falign_udpari_long(
            std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
            b1.ptrs.as_mut_ptr(), b2.ptrs.as_mut_ptr(),
            fe1.as_mut_ptr(), fe2.as_mut_ptr(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            g1.len() as c_int, g2.len() as c_int, alloclen as c_int,
            &mut fftlog,
        );
    }
    (rows_back(&b1), rows_back(&b2))
}

fn release_c() {
    unsafe {
        Falign_udpari_long(
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), 0, 0, 0, std::ptr::null_mut(),
        );
    }
}

/// C `MSalignmm` as `treebase` calls it (`disttbfast.c:2889`): no
/// `sgap`/`egap`, `outgap=0`, the children's cached cpmx slots and the
/// slot the merged cpmx is stored into (`cpmxhist`).
unsafe fn c_msalignmm_cached(
    g1: &[Vec<u8>], e1: &[f64], ori1: f64,
    g2: &[Vec<u8>], e2: &[f64], ori2: f64,
    child0: *mut *mut *mut c_double,
    child1: *mut *mut *mut c_double,
    result: *mut *mut *mut c_double,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let alloclen = (g1[0].len() + g2[0].len()) * 2 + 100;
    let mut b1 = c_rows(g1, alloclen);
    let mut b2 = c_rows(g2, alloclen);
    let mut fe1 = e1.to_vec();
    let mut fe2 = e2.to_vec();
    unsafe {
        let n_dyn = std::ptr::addr_of!(mafft_sys::n_dis_consweight_multi).read();
        mafft_sys::MSalignmm(
            n_dyn, b1.ptrs.as_mut_ptr(), b2.ptrs.as_mut_ptr(), fe1.as_mut_ptr(), fe2.as_mut_ptr(),
            g1.len() as c_int, g2.len() as c_int, alloclen as c_int,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
            0, 0, child0, child1, result, ori1, ori2,
        );
    }
    (rows_back(&b1), rows_back(&b2))
}

/// Strip the columns that are gaps in every row of the group.
fn strip_common_gaps(rows: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let w = rows[0].len();
    let keep: Vec<bool> = (0..w).map(|j| rows.iter().any(|r| r[j] != b'-')).collect();
    rows.iter().map(|r| r.iter().zip(&keep).filter(|(_, k)| **k).map(|(c, _)| *c).collect()).collect()
}

/// Mismatch counts between C's cached `cpmxhist` entry (`nalphabets + 3`
/// rows: cpmx, gap frequency, opening and closing counts) and a Rust
/// profile, over `alen` columns.
unsafe fn cache_diffs(slot: *mut *mut c_double, prof: &Profile, nalphabets: usize, alen: usize) -> [usize; 4] {
    let mut d = [0usize; 4];
    unsafe {
        for k in 0..nalphabets {
            let row = *slot.add(k);
            for j in 0..alen {
                if *row.add(j) != prof.freqs[j][k] { d[0] += 1; }
            }
        }
        let gf = *slot.add(nalphabets);
        let og = *slot.add(nalphabets + 1);
        let fg = *slot.add(nalphabets + 2);
        for j in 0..alen {
            if *gf.add(j) != prof.nongap_freq[j] { d[1] += 1; }
            if *og.add(j) != prof.ogcp[j] { d[2] += 1; }
            if *fg.add(j) != prof.fgcp[j] { d[3] += 1; }
        }
    }
    d
}

/// One level of the DP-branch chain on the Rust side: profile (fresh or
/// cached), `MSalignmm`, and the blend `MSalignmm` would store. Also
/// returns what `A__align` would store (`Salignmm.c:2119-2122`, the
/// counts it aligned with), to check that the fixture tells them apart.
fn rust_chain_level(
    scoring: &ScoringContext,
    cached1: Option<Profile>,
    a: &[Vec<u8>], wa: &[f64],
    b: &[Vec<u8>], wb: &[f64],
) -> (Vec<Vec<u8>>, Profile, Profile) {
    let (na, nb) = (normalize(wa), normalize(wb));
    let from_cache = cached1.is_some();
    let pa = cached1.unwrap_or_else(|| profile(scoring, a, &na));
    let pb = profile(scoring, b, &nb);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let aln = mafft_align::msalignmm(&pa, &pb, &scoring.consweight_matrix, &gap, false, false);
    let (q1, q2) = apply_ops(&aln.operations, a, b);
    let mut rows = q1; rows.extend(q2);
    let (gt1, gt2) = gap_tables(&aln.operations);
    // `MSalignmm.c:2436-2448`.
    let (sa, sb) = (raw_sum(wa), raw_sum(wb));
    let ns = |w: &[f64]| w.iter().fold(0.0f64, |x, y| x + y);
    let e1 = ns(&na) * sa / (sa + sb);
    let e2 = ns(&nb) * sb / (sa + sb);
    // `MSalignmm.c:2459-2460`: a cached child's counts come from the
    // zero-filled scratch buffers.
    let zeros = vec![0.0f64; pa.ogcp.len().max(pa.fgcp.len())];
    let counts: [&[f64]; 4] = [
        if from_cache { &zeros } else { &pa.ogcp },
        if from_cache { &zeros } else { &pa.fgcp },
        &pb.ogcp, &pb.fgcp,
    ];
    let blend = blend_profiles_with_counts(&pa, &pb, e1, e2, &gt1, &gt2, scoring.nalphabets, counts);
    let naive = blend_profiles_exact(&pa, &pb, e1, e2, &gt1, &gt2, scoring.nalphabets);
    (rows, blend, naive)
}

/// A chain of DP-branch merges: `splits` are the prefix lengths of `g1`
/// at which the earlier merges joined, oldest first. Each level merges the
/// prefix (cached cpmx from the previous level, fresh for the first) with
/// the next slice (fresh) and stores its cpmx; the final step merges `g1`
/// (cached) with `g2`. Every level is checked against C's rows and C's
/// cached arrays. Returns the final rows from C and from Rust.
fn check_chain(scoring: &ScoringContext, g1: &[Vec<u8>], w1: &[f64], splits: &[usize], g2: &[Vec<u8>], w2: &[f64]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let maxlen = g1[0].len().max(g2[0].len());
    unsafe { init_c_protein(g1.len() + g2.len(), maxlen * 2 + 100); }
    let mut slot_prev: *mut *mut c_double = std::ptr::null_mut();
    let mut rust_prev: Option<Profile> = None;
    let mut bounds: Vec<usize> = splits.to_vec();
    bounds.push(g1.len());
    for k in 1..bounds.len() {
        let (lo, hi) = (bounds[k - 1], bounds[k]);
        let a = strip_common_gaps(&g1[..lo]);
        let b = strip_common_gaps(&g1[lo..hi]);
        let (wa, wb) = (&w1[..lo], &w1[lo..hi]);
        let mut slot_null: *mut *mut c_double = std::ptr::null_mut();
        let mut slot_this: *mut *mut c_double = std::ptr::null_mut();
        let (pa, pb) = unsafe {
            c_msalignmm_cached(&a, &normalize(wa), raw_sum(wa), &b, &normalize(wb), raw_sum(wb),
                &mut slot_prev, &mut slot_null, &mut slot_this)
        };
        let mut c_level: Vec<Vec<u8>> = pa; c_level.extend(pb);
        let known = strip_common_gaps(&g1[..hi]);
        assert!(c_level == known, "chain level {k}: C does not reproduce the known merge of {lo}+{}", hi - lo);
        assert!(!slot_this.is_null(), "chain level {k}: C stored no cpmx (icyc + jcyc <= 20?)");
        let was_cached = rust_prev.is_some();
        let (r_level, blend, naive) = rust_chain_level(scoring, rust_prev.take(), &a, wa, &b, wb);
        assert!(r_level == known, "chain level {k}: rust MSalignmm differs from the known merge");
        let alen = known[0].len();
        let d = unsafe { cache_diffs(slot_this, &blend, scoring.nalphabets, alen) };
        assert_eq!(d, [0, 0, 0, 0], "chain level {k}: blend differs from C's cpmxhist (cpmx, gapfreq, ogcp, fgcp mismatches)");
        if was_cached {
            let dn = unsafe { cache_diffs(slot_this, &naive, scoring.nalphabets, alen) };
            assert!(dn[2] > 0 || dn[3] > 0, "chain level {k}: the fixture does not exercise the scratch-buffer counts");
        }
        rust_prev = Some(blend);
        slot_prev = slot_this;
    }
    let mut slot_null2: *mut *mut c_double = std::ptr::null_mut();
    let mut slot_this: *mut *mut c_double = std::ptr::null_mut();
    let (c1, c2) = unsafe {
        c_msalignmm_cached(g1, &normalize(w1), raw_sum(w1), g2, &normalize(w2), raw_sum(w2),
            &mut slot_prev, &mut slot_null2, &mut slot_this)
    };
    let mut c_final = c1; c_final.extend(c2);
    let (r_final, _, _) = rust_chain_level(scoring, rust_prev, g1, w1, g2, w2);
    (c_final, r_final)
}

fn first_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    a.iter().zip(b).position(|(x, y)| x != y).or(if a.len() != b.len() { Some(a.len().min(b.len())) } else { None })
}

/// The FFT branch: `Falign_udpari_long` against `fft_profile_align`.
fn check_groups(g1: &[Vec<u8>], w1: &[f64], g2: &[Vec<u8>], w2: &[f64], label: &str) {
    let scoring = scoring();
    let maxlen = g1[0].len().max(g2[0].len());
    unsafe { init_c_protein(g1.len() + g2.len(), maxlen * 2 + 100); }
    let e1 = normalize(w1);
    let e2 = normalize(w2);
    let (c1, c2) = c_udpari(g1, &e1, g2, &e2);
    release_c();
    let (r1, r2) = rust_udpari(&scoring, g1, &e1, g2, &e2);
    let col = first_diff(&c1[0], &r1[0]).or_else(|| first_diff(&c2[0], &r2[0]));
    assert!(
        c1 == r1 && c2 == r2,
        "{label}: C width {} rust width {}, first difference at column {col:?}",
        c1[0].len(),
        r1[0].len(),
    );
}

/// The DP branch of `treebase` after the `'M'` switch (`disttbfast.c:2889`)
/// on fresh profiles: `MSalignmm` on the whole profiles, no `sgap`/`egap`.
fn check_groups_dp(g1: &[Vec<u8>], w1: &[f64], g2: &[Vec<u8>], w2: &[f64], label: &str) {
    let scoring = scoring();
    let maxlen = g1[0].len().max(g2[0].len());
    unsafe { init_c_protein(g1.len() + g2.len(), maxlen * 2 + 100); }
    let e1 = normalize(w1);
    let e2 = normalize(w2);
    let mut slot: *mut *mut c_double = std::ptr::null_mut();
    let (c1, c2) = unsafe {
        c_msalignmm_cached(g1, &e1, raw_sum(w1), g2, &e2, raw_sum(w2),
            std::ptr::null_mut(), std::ptr::null_mut(), &mut slot)
    };
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let aln = mafft_align::msalignmm(&profile(&scoring, g1, &e1), &profile(&scoring, g2, &e2), &scoring.consweight_matrix, &gap, false, false);
    let (q1, q2) = apply_ops(&aln.operations, g1, g2);
    let col = first_diff(&c1[0], &q1[0]).or_else(|| first_diff(&c2[0], &q2[0]));
    assert!(c1 == q1 && c2 == q2, "{label} DP branch: C width {} rust width {}, first difference {col:?}", c1[0].len(), q1[0].len());
}

/// Deterministic synthetic family: 21 aligned rows of 44 columns from one
/// ancestor, rows 11-20 with a 6-column deletion block, plus two
/// unaligned relatives with their own short deletions. The weights are
/// those of a 40-record tree with equal branches.
fn synthetic_family() -> (Vec<Vec<u8>>, Vec<f64>, Vec<Vec<u8>>, Vec<Vec<u8>>) {
    const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
    let mut state: u64 = 0x2545F4914F6CDD1D;
    let mut next = move || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
    let base: Vec<u8> = (0..44).map(|_| AA[(next() % 20) as usize]).collect();
    let mut rows: Vec<Vec<u8>> = Vec::new();
    for r in 0..21 {
        let mut s = base.clone();
        for j in 0..44 {
            if next() % 100 < 12 { s[j] = AA[(next() % 20) as usize]; }
        }
        if r >= 11 {
            for j in 18..24 { s[j] = b'-'; }
        }
        rows.push(s);
    }
    let mut single = |drop: std::ops::Range<usize>| -> Vec<u8> {
        let mut s = base.clone();
        for j in 0..44 {
            if next() % 100 < 20 { s[j] = AA[(next() % 20) as usize]; }
        }
        s.iter().enumerate().filter(|(j, _)| !drop.contains(j)).map(|(_, c)| *c).collect()
    };
    let b = single(30..34);
    let c = single(5..8);
    let w = vec![1.0 / 40.0; 21];
    (rows, w, vec![b], vec![c])
}

#[test]
fn udpari_long_protein_groups_match_c() {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let (rows, w, b, _c) = synthetic_family();
    check_groups(&rows, &w, &b, &[1.0 / 40.0], "synthetic 21 vs 1");
}

#[test]
fn msalignmm_dp_branch_protein_groups_match_c() {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let (rows, w, b, _c) = synthetic_family();
    check_groups_dp(&rows, &w, &b, &[1.0 / 40.0], "synthetic 21 vs 1");
}

/// Bug: the port rebuilt every profile; `treebase` feeds a DP merge the
/// profile the previous DP merge stored, and `MSalignmm` blends the cached
/// child's gap counts from zeroed scratch buffers. 21+1 stores a cpmx
/// (`icyc + jcyc > 20`), 22 (cached) + 1 consumes and stores it, and the
/// final merge consumes that.
#[test]
fn msalignmm_cpmx_chain_matches_c() {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let scoring = scoring();
    let (rows, w, b, c) = synthetic_family();
    // The 22-row group as C builds it: 21 + 1, fresh.
    unsafe { init_c_protein(23, 44 * 2 + 100); }
    let (p1, p2) = unsafe {
        let mut slot: *mut *mut c_double = std::ptr::null_mut();
        c_msalignmm_cached(&rows, &normalize(&w), raw_sum(&w), &b, &[1.0], 1.0 / 40.0,
            std::ptr::null_mut(), std::ptr::null_mut(), &mut slot)
    };
    let mut g1 = p1; g1.extend(p2);
    let mut w1 = w.clone(); w1.push(1.0 / 40.0);
    let (c_final, r_final) = check_chain(&scoring, &g1, &w1, &[21], &c, &[1.0 / 40.0]);
    assert!(c_final == r_final, "final merge: C width {} rust width {}", c_final[0].len(), r_final[0].len());
}

fn read_lines(p: &std::path::Path) -> Vec<String> {
    std::fs::read_to_string(p).expect("read").lines().map(|s| s.to_string()).collect()
}

/// Investigation entry point, see the module comment.
#[test]
fn udpari_long_protein_case_dir_matches_c() {
    let Ok(dir) = std::env::var("UDPARI_CASE_DIR") else { return };
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let dir = PathBuf::from(dir);
    let g1: Vec<Vec<u8>> = read_lines(&dir.join("g1.txt")).into_iter().map(|s| s.into_bytes()).collect();
    let g2: Vec<Vec<u8>> = read_lines(&dir.join("g2.txt")).into_iter().map(|s| s.into_bytes()).collect();
    let w1: Vec<f64> = read_lines(&dir.join("w1.txt")).iter().map(|s| s.parse().unwrap()).collect();
    let w2: Vec<f64> = read_lines(&dir.join("w2.txt")).iter().map(|s| s.parse().unwrap()).collect();
    if let Ok(spec) = std::env::var("UDPARI_CHAIN") {
        let splits: Vec<usize> = spec.split(',').map(|x| x.parse().unwrap()).collect();
        let (c_final, r_final) = check_chain(&scoring(), &g1, &w1, &splits, &g2, &w2);
        eprintln!("chain final: C width {} rust width {}", c_final[0].len(), r_final[0].len());
        assert!(c_final == r_final, "case dir chain: final merge differs");
    } else if std::env::var_os("UDPARI_DP").is_some() {
        check_groups_dp(&g1, &w1, &g2, &w2, "case dir");
    } else {
        check_groups(&g1, &w1, &g2, &w2, "case dir");
    }
}
