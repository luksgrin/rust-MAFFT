//! Cross-validate `fft_profile_align` with `udpari_long` against C's
//! `Falign_udpari_long` (`Falign.c:2583`), the FFT aligner `disttbfast`
//! uses once `alg = 'M'` (`--memsave`, or a profile past 30000 columns).
//!
//! Inputs are 1-vs-1 nucleotide pairs: the `homopolymer_dna` memory-safety
//! fixture (300 × `A` + 300 nt against 300 × `A` + 295 nt), whose A-run
//! yields segments at hundreds of lags; the `long_seq_memsave_2402` pair;
//! and two low-complexity sequences that take the many-anchor chain.
//! On nucleotides the correlation is `seq_vec_4` (`Falign.c:83-99`), one
//! complex channel, not `Falign`'s `a` indicator.

#![allow(non_snake_case)]

use std::os::raw::{c_char, c_double, c_int};
use std::path::PathBuf;
use std::sync::Mutex;

use mafft_align::{fft_profile_align, fft_total_eff, AlignOp, FftAlignParams, GapModel, NucRows, Profile};

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

fn load(rel: &str) -> Vec<Vec<u8>> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../mafft-bin/tests/fixtures").join(rel);
    let set = mafft_io::read_fasta(&path).expect("read fixture");
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
        // `Falign_udpari_long` only has a `case 'M'` (`Falign.c:3350`).
        std::ptr::addr_of_mut!(mafft_sys::alg).write(b'M' as c_char);
        std::ptr::addr_of_mut!(mafft_sys::fftWinSize).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftThreshold).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftkeika).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kobetsubunkatsu).write(0);
        let seq_data = b"acgt\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut c_char;
        mafft_sys::constants(1, &mut seq_ptr);
        std::ptr::addr_of_mut!(mafft_sys::njob).write(njob as c_int);
        std::ptr::addr_of_mut!(mafft_sys::nlenmax).write(nlenmax as c_int);
        // FFT-NS-2 gets `-O` from the script: outgap = 0.
        std::ptr::addr_of_mut!(mafft_sys::outgap).write(0);
    }
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

fn rust_udpari(
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
    params.udpari_long = true;
    params.nuc_rows = Some(NucRows { rows1: g1.to_vec(), eff1: e1.to_vec(), rows2: g2.to_vec(), eff2: e2.to_vec() });
    let aln = fft_profile_align(&p1, &p2, &scoring.consweight_matrix, &params);
    apply_ops(&aln.operations, g1, g2)
}

fn c_udpari(g1: &[Vec<u8>], e1: &[f64], g2: &[Vec<u8>], e2: &[f64]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
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
        Falign_udpari_long(
            std::ptr::null_mut(), std::ptr::null_mut(), n_dyn,
            p1.as_mut_ptr(), p2.as_mut_ptr(),
            fe1.as_mut_ptr(), fe2.as_mut_ptr(),
            std::ptr::null_mut(), std::ptr::null_mut(),
            g1.len() as c_int, g2.len() as c_int, alloclen as c_int,
            &mut fftlog,
        );
    }
    let rows = |b: &[Vec<u8>]| -> Vec<Vec<u8>> {
        b.iter().map(|v| v[..v.iter().position(|&c| c == 0).unwrap()].to_vec()).collect()
    };
    (rows(&b1), rows(&b2))
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

fn first_diff(a: &[u8], b: &[u8]) -> Option<usize> {
    a.iter().zip(b).position(|(x, y)| x != y).or(if a.len() != b.len() { Some(a.len().min(b.len())) } else { None })
}

fn check_pair(rel: &str, i: usize, j: usize) {
    let _g = C_MUTEX.lock().unwrap_or_else(std::sync::PoisonError::into_inner);
    let seqs = load(rel);
    let scoring = mafft_scoring::build_context(mafft_types::ScoringModel::Dna, mafft_types::SeqType::Dna);
    let maxlen = seqs.iter().map(|s| s.len()).max().unwrap();
    unsafe { init_c_dna(seqs.len(), maxlen * 2 + 100); }
    let (c1, c2) = c_udpari(&seqs[i..=i], &[1.0], &seqs[j..=j], &[1.0]);
    release_c();
    let (r1, r2) = rust_udpari(&scoring, &seqs[i..=i], &[1.0], &seqs[j..=j], &[1.0]);
    let col = first_diff(&c1[0], &r1[0]).or_else(|| first_diff(&c2[0], &r2[0]));
    assert!(
        c1 == r1 && c2 == r2,
        "{rel} ({i}, {j}): C width {} rust width {}, first difference at column {col:?}",
        c1[0].len(),
        r1[0].len(),
    );
}

#[test]
fn udpari_long_homopolymer_pair_matches_c() {
    check_pair("memsafety/homopolymer_dna.fa", 0, 1);
}

#[test]
fn udpari_long_long_seq_2402_pair_matches_c() {
    check_pair("long_seq_memsave_2402.fa", 0, 1);
}



/// Two low-complexity sequences of about 6.3 kb. They yield more than 5000
/// segments over the 500 lags, so `Falign_udpari_long` skips `blockAlign2`
/// and keeps a monotonic chain of centers in `mymergesort` order
/// (`Falign.c:3164-3202`).
#[test]
fn udpari_long_many_anchor_chain_matches_c() {
    check_pair("memsafety/adjustdirection_lowcomplexity.fa", 0, 2);
}
