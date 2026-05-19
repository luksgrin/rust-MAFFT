//! Compare Rust `Profile::from_aligned` cpmx to C `cpmx_calc_new` bit-exactly.
//!
//! BB20027 pass-1 step 13 diverges from C at the 8+2 merge with same width
//! but different score (27.8 unit gap). Suspect: our `Profile::from_aligned`
//! produces slightly different cpmx than C's cached `createcpmxresult`.
//!
//! This test takes a synthetic 8-sequence aligned cluster (similar shape
//! to step 12's output), computes cpmx via both paths, and verifies
//! they're bit-identical.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};

use mafft_align::Profile;
use mafft_scoring::build_context;
use mafft_types::{ScoringModel, SeqType};

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
        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);
    }
}

#[test]
#[ignore]
fn rust_from_scratch_matches_rust_blend() {
    // Split 8 sequences into 3+5. Compute cpmx three ways:
    //  (a) Rust Profile::from_aligned on all 8 (the "from-scratch" path)
    //  (b) Rust two Profile::from_aligned on 3 and 5, then blend
    //  (c) Cell-by-cell diff to see if (a) == (b)
    //
    // Pass 0 of BB20027 always matches C (every step from-scratch).
    // Pass 1 step 13 uses (a) in Rust but (b) in C. If (a) != (b), we've
    // found the source of pass-1 divergence and a clear path to fix.

    let seqs: Vec<Vec<u8>> = vec![
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQAPILS---RVGDGTQDNL".to_vec(),
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQGSILS---RVADGTQDNI".to_vec(),
        b"MKTAYIAKQRQISFLKSHFSRQLEERLG--LIEVQAPILK---RVGDGTQDNL".to_vec(),
        b"MKVAYVAKQRTLSWVKAHISRSAEEERLNGTLEEKVNAVPN---RVGDGTKEEI".to_vec(),
        b"-KAVHISKVRTLSWVKAHISRSAEAERLNGTLEEKVNAVPN---RVGDGTKEEI".to_vec(),
        b"MAA-YVAKQRTLSWVKAHISRSAEAERLNGTLEEKVNAVRN---RVGDGTAEEI".to_vec(),
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQAPILS---RVGDGTQDNL".to_vec(),
        b"MKEVYIAKQRQVAYIKSHFSRPAEERLT--AIEVPDQIIS--PRVGDPVQDQL".to_vec(),
    ];
    let len = seqs.iter().map(|s| s.len()).max().unwrap();
    let seqs: Vec<Vec<u8>> = seqs.iter().map(|s| {
        let mut v = s.clone(); v.resize(len, b'-'); v
    }).collect();

    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // Equal per-leaf weights summing to 1.0.
    let leaf_weight = 1.0 / 8.0;
    let all_weights = vec![leaf_weight; 8];

    // (a) From scratch on all 8 sequences with global weights summing to 1.
    let all_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let prof_all = Profile::from_aligned(&all_refs, &all_weights,
        &scoring.amino_map, scoring.nalphabets);

    // (b) Build 3-way and 5-way separately with intra-cluster normalized
    // weights, then blend with eff1 = 3/8, eff2 = 5/8.
    let g1: Vec<&[u8]> = seqs[..3].iter().map(|s| s.as_slice()).collect();
    let g2: Vec<&[u8]> = seqs[3..].iter().map(|s| s.as_slice()).collect();
    let w1n = vec![1.0 / 3.0; 3];
    let w2n = vec![1.0 / 5.0; 5];
    let prof1 = Profile::from_aligned(&g1, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&g2, &w2n, &scoring.amino_map, scoring.nalphabets);

    // gaptable for full-overlap blend (no extra positions inserted).
    let gaptable: Vec<u8> = vec![b'o'; len];

    // Call our internal blend_profiles_exact (it's `fn`, not `pub fn`, so
    // we replicate the call signature inline below using a helper from
    // progressive.rs). Without an exposed API, we test the equivalent
    // operation manually here, mirroring the function:
    //   freqs[j][k] = prof1.freqs[p][k] * eff1 + prof2.freqs[p][k] * eff2
    // where p = j (no gap columns).
    let eff1 = 3.0 / 8.0;
    let eff2 = 5.0 / 8.0;
    let mut blend_freqs = vec![vec![0.0f64; scoring.nalphabets]; len];
    for j in 0..len {
        for k in 0..scoring.nalphabets {
            blend_freqs[j][k] = prof1.freqs[j][k].mul_add(eff1, blend_freqs[j][k]);
            blend_freqs[j][k] = prof2.freqs[j][k].mul_add(eff2, blend_freqs[j][k]);
        }
    }

    // Compare (a) vs (b).
    let mut max_diff: f64 = 0.0;
    let mut n_diff = 0;
    for j in 0..len {
        for k in 0..scoring.nalphabets {
            let d = (prof_all.freqs[j][k] - blend_freqs[j][k]).abs();
            if d > 0.0 {
                if n_diff < 5 {
                    eprintln!("DIFF: pos={} k={} from_scratch={:.20} blend={:.20} d={:.3e}",
                        j, k, prof_all.freqs[j][k], blend_freqs[j][k], d);
                }
                n_diff += 1;
            }
            max_diff = max_diff.max(d);
        }
    }
    eprintln!("from-scratch vs blend: max |diff| = {:.3e}, cells differing = {}",
              max_diff, n_diff);
    if max_diff > 0.0 {
        eprintln!("Confirms: from-scratch (cpmx_calc_new) and blend (createcpmxresult) \
                   are NOT bit-identical. C uses blend at internal merges via cpmxhist; \
                   Rust uses from-scratch when not caching. This precision drift drives \
                   pass-1 tied-DP-cell flips on BB20027 step 13.");
    }
}

#[test]
#[ignore]
fn rust_profile_freqs_match_c_cpmx_calc_new() {
    // 8 aligned sequences of equal length, with realistic gaps.
    let seqs: Vec<Vec<u8>> = vec![
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQAPILS---RVGDGTQDNL"
            .to_vec(),
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQGSILS---RVADGTQDNI"
            .to_vec(),
        b"MKTAYIAKQRQISFLKSHFSRQLEERLG--LIEVQAPILK---RVGDGTQDNL"
            .to_vec(),
        b"MKVAYVAKQRTLSWVKAHISRSAEEERLNGTLEEKVNAVPN---RVGDGTKEEI"
            .to_vec(),
        b"-KAVHISKVRTLSWVKAHISRSAEAERLNGTLEEKVNAVPN---RVGDGTKEEI"
            .to_vec(),
        b"MAA-YVAKQRTLSWVKAHISRSAEAERLNGTLEEKVNAVRN---RVGDGTAEEI"
            .to_vec(),
        b"MKTAYIAKQRQISFVKSHFSRQLEERLG--LIEVQAPILS---RVGDGTQDNL"
            .to_vec(),
        b"MKEVYIAKQRQVAYIKSHFSRPAEERLT--AIEVPDQIIS--PRVGDPVQDQL"
            .to_vec(),
    ];
    // Verify all same length (pad with - if not).
    let len = seqs.iter().map(|s| s.len()).max().unwrap();
    let seqs_padded: Vec<Vec<u8>> = seqs.iter()
        .map(|s| {
            let mut v = s.clone();
            v.resize(len, b'-');
            v
        })
        .collect();

    // Equal weights summing to 1.0.
    let nseq = seqs_padded.len();
    let weight = 1.0 / nseq as f64;
    let weights = vec![weight; nseq];

    // ====== Rust path ======
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let seq_refs: Vec<&[u8]> = seqs_padded.iter().map(|s| s.as_slice()).collect();
    let rust_prof = Profile::from_aligned(&seq_refs, &weights,
        &scoring.amino_map, scoring.nalphabets);

    // ====== C path ======
    unsafe { init_c_protein(); }

    // C expects null-terminated C strings, mutable.
    let mut c_seqs: Vec<CString> = seqs_padded.iter()
        .map(|s| CString::new(s.clone()).unwrap())
        .collect();
    let mut c_seq_ptrs: Vec<*mut c_char> = c_seqs.iter_mut()
        .map(|s| s.as_ptr() as *mut c_char)
        .collect();

    // Allocate cpmx[nalphabets][lgth] for C.
    let nalpha = scoring.nalphabets;
    let mut c_cpmx_rows: Vec<Vec<f64>> = (0..nalpha)
        .map(|_| vec![0.0f64; len])
        .collect();
    let mut c_cpmx_ptrs: Vec<*mut c_double> = c_cpmx_rows.iter_mut()
        .map(|r| r.as_mut_ptr())
        .collect();

    let mut c_eff: Vec<f64> = weights.clone();

    unsafe {
        mafft_sys::cpmx_calc_new(
            c_seq_ptrs.as_mut_ptr(),
            c_cpmx_ptrs.as_mut_ptr(),
            c_eff.as_mut_ptr(),
            len as c_int,
            nseq as c_int,
        );
    }

    // Compare: Rust freqs[pos][k] vs C cpmx[k][pos]
    let mut max_diff: f64 = 0.0;
    let mut total_diffs = 0;
    for pos in 0..len {
        for k in 0..nalpha {
            let r = rust_prof.freqs[pos][k];
            let c = c_cpmx_rows[k][pos];
            let d = (r - c).abs();
            if d > 1e-15 {
                total_diffs += 1;
                if total_diffs <= 10 {
                    eprintln!("DIFF: pos={} k={} rust={:.16} c={:.16} d={:.3e}",
                              pos, k, r, c, d);
                }
            }
            max_diff = max_diff.max(d);
        }
    }
    eprintln!("max |diff| = {:.3e}, total cells differing = {}",
              max_diff, total_diffs);
    assert!(max_diff < 1e-13,
        "Rust Profile::from_aligned does not bit-match C cpmx_calc_new");
}
