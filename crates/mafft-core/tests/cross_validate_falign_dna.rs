//! Cross-validate the nucleotide FFT anchor search and `fft_profile_align`
//! against C's `Falign` on long coding sequences (issue #43).
//!
//! Input: `dna_01715.fa`, three *S. pneumoniae* CDS of 3.6–5.0 kb on which
//! the FFT progressive modes diverged from C while `--nofft`, `--localpair`
//! and `--globalpair` matched.
//!
//! * `soukan_matches_c_bitwise` rebuilds C's correlation array from C's own
//!   `fft` and `calcNaiseki` with `Falign`'s nucleotide channel layout (one
//!   channel, `n20or4or2 = 1`, `Falign.c:923`) and requires
//!   `mafft_fft::soukan_c_compat` to equal it bit for bit, and the 20
//!   candidate lags of `getKouho` to agree.
//! * `falign_replay_matches_c` runs every pairwise merge (1 vs 1) and the
//!   2-vs-1 merge on C's own pair output through both C `Falign` and Rust
//!   `fft_profile_align`, and requires identical rows.

#![allow(non_snake_case)]

use std::os::raw::{c_char, c_double, c_int};
use std::path::PathBuf;
use std::sync::Mutex;

use mafft_align::{fft_profile_align, fft_total_eff, AlignOp, FftAlignParams, GapModel, Profile};
use mafft_fft::soukan_c_compat;
use num_complex::Complex64;

static C_MUTEX: Mutex<()> = Mutex::new(());

unsafe extern "C" {
    fn calcNaiseki(value: *mut mafft_sys::Fukusosuu, x: *mut mafft_sys::Fukusosuu, y: *mut mafft_sys::Fukusosuu);
    fn getKouho(kouho: *mut c_int, nkouho: c_int, soukan: *mut c_double, nlen2: c_int) -> c_int;
}

fn load_input() -> Vec<Vec<u8>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../mafft-bin/tests/fixtures/dna_01715.fa");
    let set = mafft_io::read_fasta(&path).expect("read dna_01715.fa");
    set.sequences.into_iter().map(|s| s.data.to_ascii_lowercase()).collect()
}

unsafe fn init_c_dna(njob: usize, nlenmax: usize) {
    unsafe {
        mafft_sys::initglobalvariables();
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'd' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(-1);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::alg).write(b'A' as c_char);
        std::ptr::addr_of_mut!(mafft_sys::fftWinSize).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftThreshold).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(0);
        let seq_data = b"acgt\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut c_char;
        mafft_sys::constants(1, &mut seq_ptr);
        std::ptr::addr_of_mut!(mafft_sys::njob).write(njob as c_int);
        std::ptr::addr_of_mut!(mafft_sys::nlenmax).write(nlenmax as c_int);
        // The mafft script passes `-O` to disttbfast for FFT-NS-2: outgap = 0.
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(0);
    }
}

/// C `Falign`'s FFT length: smallest power of two strictly above `len`.
fn c_nlen(len: usize) -> usize {
    let mut nlen = 1;
    while len >= nlen { nlen <<= 1; }
    nlen
}

/// `seq_vec_3` with `n20or4or2 = 1`: the `a` indicator, weighted.
fn a_channel(seqs: &[&[u8]], effs: &[f64], nlen: usize) -> Vec<Complex64> {
    let mut v = vec![Complex64::new(0.0, 0.0); nlen];
    for (s, &e) in seqs.iter().zip(effs) {
        for (i, &c) in s.iter().enumerate() {
            if c == b'a' { v[i].re += e; }
        }
    }
    v
}

fn c_soukan(a: &[Complex64], b: &[Complex64]) -> Vec<f64> {
    let nlen = a.len();
    let mut x: Vec<mafft_sys::Fukusosuu> = a.iter().map(|c| mafft_sys::Fukusosuu { R: c.re, I: c.im }).collect();
    let mut y: Vec<mafft_sys::Fukusosuu> = b.iter().map(|c| mafft_sys::Fukusosuu { R: c.re, I: c.im }).collect();
    let mut w = vec![mafft_sys::Fukusosuu { R: 0.0, I: 0.0 }; nlen];
    unsafe {
        mafft_sys::fft(nlen as c_int, y.as_mut_ptr(), 0);
        mafft_sys::fft(nlen as c_int, x.as_mut_ptr(), 0);
        for l in 0..nlen {
            let mut v = mafft_sys::Fukusosuu { R: 0.0, I: 0.0 };
            calcNaiseki(&mut v, x.as_mut_ptr().add(l), y.as_mut_ptr().add(l));
            w[l].R += v.R;
            w[l].I += v.I;
        }
        mafft_sys::fft(-(nlen as c_int), w.as_mut_ptr(), 0);
    }
    let nlen2 = nlen / 2;
    let mut s = vec![0.0; nlen];
    for m in 0..=nlen2 { s[m] = w[nlen2 - m].R; }
    for m in nlen2 + 1..nlen { s[m] = w[nlen + nlen2 - m].R; }
    s
}

fn c_kouho(soukan: &[f64]) -> Vec<i32> {
    let mut s = soukan.to_vec();
    let mut k = vec![0 as c_int; 20];
    unsafe { getKouho(k.as_mut_ptr(), 20, s.as_mut_ptr(), s.len() as c_int); }
    k
}

fn rust_kouho(soukan: &[f64]) -> Vec<i32> {
    mafft_fft::get_top_candidates(soukan, 20).into_iter().map(|c| c.lag).collect()
}

#[test]
fn soukan_matches_c_bitwise() {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let seqs = load_input();
    for i in 0..seqs.len() {
        for j in (i + 1)..seqs.len() {
            let nlen = c_nlen(seqs[i].len().max(seqs[j].len()));
            let a = a_channel(&[&seqs[i]], &[1.0], nlen);
            let b = a_channel(&[&seqs[j]], &[1.0], nlen);
            let c = c_soukan(&a, &b);
            let r = soukan_c_compat(&[a.clone()], &[b.clone()]);
            let ndiff = c.iter().zip(&r).filter(|(x, y)| x.to_bits() != y.to_bits()).count();
            assert_eq!(ndiff, 0, "pair ({i},{j}): {ndiff}/{nlen} soukan cells differ from C");
            assert_eq!(rust_kouho(&r), c_kouho(&c), "pair ({i},{j}): candidate lags differ");
        }
    }
    unsafe { mafft_sys::fft(0, std::ptr::null_mut(), 1); }
}

/// Rows of an alignment given as operations over `g1` / `g2`.
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

fn rust_falign(
    scoring: &mafft_types::ScoringContext,
    g1: &[Vec<u8>], e1: &[f64],
    g2: &[Vec<u8>], e2: &[f64],
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let r1: Vec<&[u8]> = g1.iter().map(|s| s.as_slice()).collect();
    let r2: Vec<&[u8]> = g2.iter().map(|s| s.as_slice()).collect();
    let p1 = Profile::from_aligned(&r1, e1, &scoring.amino_map, scoring.nalphabets);
    let p2 = Profile::from_aligned(&r2, e2, &scoring.amino_map, scoring.nalphabets);
    let mut params = FftAlignParams::dna();
    params.gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    params.head_gap = false;
    params.tail_gap = false;
    params.total_eff = fft_total_eff(e1, e2);
    let aln = fft_profile_align(&p1, &p2, &scoring.consweight_matrix, &params);
    apply_ops(&aln.operations, g1, g2)
}

fn c_falign(g1: &[Vec<u8>], e1: &[f64], g2: &[Vec<u8>], e2: &[f64]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let alloclen = (g1[0].len() + g2[0].len()) * 2 + 100;
    let mk = |g: &[Vec<u8>]| -> Vec<Vec<u8>> {
        g.iter().map(|s| { let mut v = s.clone(); v.resize(alloclen + 1, 0); v }).collect()
    };
    let mut b1 = mk(g1);
    let mut b2 = mk(g2);
    let mut p1: Vec<*mut c_char> = b1.iter_mut().map(|v| v.as_mut_ptr() as *mut c_char).collect();
    let mut p2: Vec<*mut c_char> = b2.iter_mut().map(|v| v.as_mut_ptr() as *mut c_char).collect();
    let mut fe1 = e1.to_vec();
    let mut fe2 = e2.to_vec();
    let mut fftlog: c_int = 0;
    unsafe {
        let n_dyn = std::ptr::addr_of!(mafft_sys::n_dis_consweight_multi).read();
        mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
            p1.as_mut_ptr(), p2.as_mut_ptr(),
            fe1.as_mut_ptr(), fe2.as_mut_ptr(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            g1.len() as c_int, g2.len() as c_int, alloclen as c_int,
            &mut fftlog, std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
    }
    let rows = |b: &[Vec<u8>]| -> Vec<Vec<u8>> {
        b.iter().map(|v| v[..v.iter().position(|&c| c == 0).unwrap()].to_vec()).collect()
    };
    (rows(&b1), rows(&b2))
}

fn first_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    a.iter().zip(b).position(|(x, y)| x != y).or(if a.len() != b.len() { Some(a.len().min(b.len())) } else { None })
}

#[test]
fn falign_replay_matches_c() {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let seqs = load_input();
    let scoring = mafft_scoring::build_context(mafft_types::ScoringModel::Dna, mafft_types::SeqType::Dna);
    let maxlen = seqs.iter().map(|s| s.len()).max().unwrap();
    unsafe { init_c_dna(seqs.len(), maxlen * 2 + 100); }
    assert_eq!(unsafe { std::ptr::addr_of!(mafft_sys::penalty).read() }, scoring.gap.open as c_int);
    assert_eq!(unsafe { std::ptr::addr_of!(mafft_sys::penalty_ex).read() }, scoring.gap.extend as c_int);

    let mut failures = Vec::new();
    let mut check = |label: String, g1: &[Vec<u8>], e1: &[f64], g2: &[Vec<u8>], e2: &[f64]| {
        let (c1, c2) = c_falign(g1, e1, g2, e2);
        let (r1, r2) = rust_falign(&scoring, g1, e1, g2, e2);
        if c1 != r1 || c2 != r2 {
            let col = first_diff(&c1[0], &r1[0]).or_else(|| first_diff(&c2[0], &r2[0]));
            eprintln!("{label}: C width {} rust width {}, first diff col {col:?}", c1[0].len(), r1[0].len());
            failures.push(label);
        }
        (c1, c2)
    };

    for i in 0..seqs.len() {
        for j in (i + 1)..seqs.len() {
            let (c1, c2) = check(format!("pair ({i},{j})"), &seqs[i..=i], &[1.0], &seqs[j..=j], &[1.0]);
            // Then the remaining sequence against C's pair alignment, with
            // unequal weights so `totaleff` and the 2-sequence profile are
            // exercised.
            let k = (0..seqs.len()).find(|&k| k != i && k != j).unwrap();
            let grp: Vec<Vec<u8>> = vec![c1[0].clone(), c2[0].clone()];
            let w = [0.3711, 0.6289];
            check(format!("({i},{j}) vs {k}"), &grp, &w, &seqs[k..=k], &[1.0]);
            check(format!("{k} vs ({i},{j})"), &seqs[k..=k], &[1.0], &grp, &w);
        }
    }
    unsafe {
        mafft_sys::Falign(
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), 0, 0, 0,
            std::ptr::null_mut(), std::ptr::null_mut(), 0, std::ptr::null_mut(),
        );
    }
    assert!(failures.is_empty(), "Rust fft_profile_align differs from C Falign on: {failures:?}");
}
