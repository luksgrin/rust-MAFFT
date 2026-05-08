/// Cross-validate Rust's FFT segment-detection against C's
/// `alignableReagion` on a BL50 input. Confirms our `match_score` /
/// `alignable_segments` produce byte-identical scores to C — eliminates
/// FFT scoring as the source of the BL50 §4 divergence.
///
/// Also includes a step-33 profile-DP comparison test that drives the
/// progressive merge to step 32, then compares Rust's `profile_align`
/// fallback to C's `A__align` on the divergent input.

use std::ffi::CString;
use std::os::raw::{c_char, c_double, c_int};
use std::sync::Mutex;

use mafft_align::{Profile, profile_align, GapModel};
use mafft_io::read_fasta;
use mafft_scoring::build_context;
use mafft_tree::{DistanceMatrix, ClusterMethod, ktuple_distance, musclesupg};
use mafft_types::{ScoringModel, SeqType};

static C_MUTEX: Mutex<()> = Mutex::new(());

unsafe fn alloc_zeroed(size: usize) -> *mut u8 {
    let layout = std::alloc::Layout::from_size_align(size.max(8), 8).unwrap();
    unsafe { std::alloc::alloc_zeroed(layout) }
}

unsafe fn init_c_protein_bl50() {
    unsafe {
        mafft_sys::initglobalvariables();
        std::ptr::addr_of_mut!(mafft_sys::ppenalty).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_ex).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_EX).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_OP).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::ppenalty_dist).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::poffset).write(0);
        std::ptr::addr_of_mut!(mafft_sys::kimuraR).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::pamN).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::dorp).write(b'p' as i32);
        std::ptr::addr_of_mut!(mafft_sys::scoremtx).write(1);
        std::ptr::addr_of_mut!(mafft_sys::nblosum).write(50);
        std::ptr::addr_of_mut!(mafft_sys::fmodel).write(0);
        std::ptr::addr_of_mut!(mafft_sys::fftWinSize).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::fftThreshold).write(mafft_sys::NOTSPECIFIED);
        std::ptr::addr_of_mut!(mafft_sys::consweight_multi).write(1.0);

        let seq_data = b"ACDEFGHIKLMNPQRSTVWY\0";
        let mut seq_ptr = seq_data.as_ptr() as *mut i8;
        let seq_arr: *mut *mut i8 = &mut seq_ptr;
        mafft_sys::constants(1, seq_arr);
    }
}

#[test]
fn bl50_alignable_reagion_matches_c() {
    let _guard = C_MUTEX.lock().unwrap();

    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");
    let scoring = build_context(ScoringModel::Blosum(50), SeqType::Protein);
    let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    // Pick the first two sequences. These are 353 and 348 residues long,
    // exercising the same per-position scoring as the post-step-32 13×8
    // cluster pair (which contains them) without needing to drive the
    // full progressive merge.
    let s1 = &sequences[0];
    let s2 = &sequences[1];

    // Rust: build single-sequence profiles and compute per-position
    // scores at lag=0.
    let prof1 = Profile::from_aligned(
        &[s1.as_slice()],
        &[1.0],
        &scoring.amino_map,
        scoring.nalphabets,
    );
    let prof2 = Profile::from_aligned(
        &[s2.as_slice()],
        &[1.0],
        &scoring.amino_map,
        scoring.nalphabets,
    );

    let len = prof1.length.min(prof2.length);
    let mut rust_scores = vec![0.0f64; len];
    for i in 0..len {
        rust_scores[i] = prof1.match_score(i, &prof2, i, &scoring.substitution_matrix);
    }

    // Run Rust segment detection (matches the BL50 default protein params).
    let rust_segs = mafft_fft::alignable_segments(
        &rust_scores,
        &mafft_fft::SegmentParams::protein(),
    );

    // C: call alignableReagion with the same input.
    let mut c_seg_data: Vec<(i32, i32, i32, f64)> = Vec::new();
    unsafe {
        init_c_protein_bl50();

        let cs1 = CString::new(s1.clone()).unwrap();
        let cs2 = CString::new(s2.clone()).unwrap();
        let mut cs1_ptr: *mut c_char = cs1.as_ptr() as *mut c_char;
        let mut cs2_ptr: *mut c_char = cs2.as_ptr() as *mut c_char;

        let mut eff1 = vec![1.0f64];
        let mut eff2 = vec![1.0f64];

        const MAX_SEG: usize = 1000;
        let seg_buf: *mut mafft_sys::Segment = alloc_zeroed(
            MAX_SEG * std::mem::size_of::<mafft_sys::Segment>()
        ) as _;

        let s1_arr: *mut *mut c_char = &mut cs1_ptr;
        let s2_arr: *mut *mut c_char = &mut cs2_ptr;

        let count = mafft_sys::alignableReagion(
            1, 1,
            s1_arr, s2_arr,
            eff1.as_mut_ptr(), eff2.as_mut_ptr(),
            seg_buf,
        ) as usize;

        for i in 0..count {
            let s = &*seg_buf.add(i);
            c_seg_data.push((s.start, s.end, s.center, s.score));
        }

        // Cleanup
        mafft_sys::alignableReagion(0, 0, std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
        mafft_sys::freeconstants();
    }

    eprintln!("Rust: {} segments, C: {} segments", rust_segs.len(), c_seg_data.len());
    assert_eq!(rust_segs.len(), c_seg_data.len(),
        "BL50 segment count differs: Rust={} C={}",
        rust_segs.len(), c_seg_data.len());

    for (i, (rs, (cs, ce, cc, cscore))) in rust_segs.iter().zip(c_seg_data.iter()).enumerate() {
        assert_eq!(rs.start, *cs as usize, "seg[{i}] start: rust={} c={}", rs.start, cs);
        assert_eq!(rs.end, *ce as usize, "seg[{i}] end: rust={} c={}", rs.end, ce);
        assert_eq!(rs.center, *cc as usize, "seg[{i}] center: rust={} c={}", rs.center, cc);
        assert!((rs.score - cscore).abs() < 1e-6,
            "seg[{i}] score: rust={} c={}", rs.score, cscore);
    }
}

/// Drive Rust progressive merge through step 32 with `MAFFT_STOP_AT_STEP=33`,
/// then compare Rust's `profile_align` and C's `A__align` outputs on the
/// step-33 13×8 cluster pair. Identifies the first divergent residue
/// position to nail down the BL50 direct-DP bug.
#[test]
#[ignore = "TODO §4 — diagnostic, not a regression test (depends on dump file)"]
fn bl50_step33_profile_dp_matches_c() {
    let _guard = C_MUTEX.lock().unwrap();

    // Build the same topology the engine uses.
    let input = read_fasta(std::path::Path::new("../../mafft-upstream/test/sample"))
        .expect("load mafft-upstream/test/sample");
    let scoring = build_context(ScoringModel::Blosum(50), SeqType::Protein);
    let nseq = input.sequences.len();
    let raw_seqs: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();

    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            dm.set(i, j, ktuple_distance(&raw_seqs[i], &raw_seqs[j], 6));
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());

    // Read step-33 input from the dump file (run mafft-rs with
    // MAFFT_STOP_AT_STEP=33 first to populate it).
    let dump_path = "/tmp/bl50.step33.txt";
    let dump = match std::fs::read_to_string(dump_path) {
        Ok(s) => s,
        Err(_) => {
            eprintln!("WARN: {dump_path} not found; run mafft-rs with MAFFT_STOP_AT_STEP=33 first");
            return;
        }
    };
    let mut sequences: Vec<Vec<u8>> = vec![Vec::new(); nseq];
    for line in dump.lines() {
        if let Some(rest) = line.strip_prefix("seq") {
            if let Some((idx_s, seq)) = rest.split_once(' ') {
                if let Ok(idx) = idx_s.parse::<usize>() {
                    sequences[idx] = seq.as_bytes().to_vec();
                }
            }
        }
    }
    eprintln!("Loaded step-33 dump: nseq={}, width={}", sequences.len(), sequences[0].len());

    let target_step: usize = std::env::var("BL50_TEST_STEP")
        .ok().and_then(|s| s.parse().ok()).unwrap_or(33);
    let step33 = &topo.steps[target_step];
    let group1: Vec<usize> = step33.left.clone();
    let group2: Vec<usize> = step33.right.clone();
    eprintln!("Step 33: |g1|={} g1={:?}", group1.len(), group1);
    eprintln!("Step 33: |g2|={} g2={:?}", group2.len(), group2);
    eprintln!("Step 33: g1[0].len={} g2[0].len={}",
        sequences[group1[0]].len(), sequences[group2[0]].len());

    let weights = mafft_tree::sequence_weights(&topo);
    let w1: Vec<f64> = group1.iter().map(|&i| weights[i]).collect();
    let w2: Vec<f64> = group2.iter().map(|&i| weights[i]).collect();
    let s1w: f64 = w1.iter().sum();
    let s2w: f64 = w2.iter().sum();
    let w1n: Vec<f64> = w1.iter().map(|w| w / s1w).collect();
    let w2n: Vec<f64> = w2.iter().map(|w| w / s2w).collect();
    eprintln!("Rust w1n: {:?}", w1n.iter().map(|x| format!("{:.10}", x)).collect::<Vec<_>>());
    eprintln!("Rust w2n: {:?}", w2n.iter().map(|x| format!("{:.10}", x)).collect::<Vec<_>>());

    let s1_refs: Vec<&[u8]> = group1.iter().map(|&i| sequences[i].as_slice()).collect();
    let s2_refs: Vec<&[u8]> = group2.iter().map(|&i| sequences[i].as_slice()).collect();

    let prof1 = Profile::from_aligned(&s1_refs, &w1n, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&s2_refs, &w2n, &scoring.amino_map, scoring.nalphabets);
    eprintln!("Profiles: prof1.length={} prof2.length={}", prof1.length, prof2.length);

    // Dump test profiles for diff against engine dump.
    if let Ok(path) = std::env::var("MAFFT_TEST_PROFILE_PATH") {
        let mut out = String::new();
        out.push_str("# test step 33\n");
        out.push_str(&format!("# prof1.length={} prof2.length={}\n",
            prof1.length, prof2.length));
        for k in 0..prof1.length {
            out.push_str(&format!("p1[{k}] gap={:.6} nong={:.6}",
                prof1.gap_freq[k], prof1.nongap_freq[k]));
            for a in 0..prof1.nalphabets {
                if prof1.freqs[k][a] != 0.0 {
                    out.push_str(&format!(" f{a}={:.6}", prof1.freqs[k][a]));
                }
            }
            out.push_str(&format!(" o={:.6} f={:.6}\n",
                prof1.ogcp[k], prof1.fgcp[k]));
        }
        for k in 0..prof2.length {
            out.push_str(&format!("p2[{k}] gap={:.6} nong={:.6}",
                prof2.gap_freq[k], prof2.nongap_freq[k]));
            for a in 0..prof2.nalphabets {
                if prof2.freqs[k][a] != 0.0 {
                    out.push_str(&format!(" f{a}={:.6}", prof2.freqs[k][a]));
                }
            }
            out.push_str(&format!(" o={:.6} f={:.6}\n",
                prof2.ogcp[k], prof2.fgcp[k]));
        }
        std::fs::write(&path, out).ok();
        eprintln!("[TEST_PROFILE_DUMP] dumped to {path}");
    }

    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Rust DP: mirror the FFT fallback path (head_gap=false, tail_gap=false
    // because the engine's FftAlignParams sets both to false).
    let rust_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, false, false);
    let rust_width = rust_aln.operations.len();
    eprintln!("Rust profile_align: width={} score={:.2}", rust_width, rust_aln.score);

    // C DP: A__align on the same input.
    let mut c_seq1: Vec<CString> = group1.iter()
        .map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();
    let mut c_seq2: Vec<CString> = group2.iter()
        .map(|&i| CString::new(sequences[i].clone()).unwrap()).collect();

    unsafe {
        init_c_protein_bl50();
        std::ptr::addr_of_mut!(mafft_sys::njob).write(nseq as c_int);
        // Pull penalty values from C's globals after constants() built them.
        let c_penalty = std::ptr::addr_of!(mafft_sys::penalty).read();
        let c_penalty_ex = std::ptr::addr_of!(mafft_sys::penalty_ex).read();
        eprintln!("C penalty={} penalty_ex={}", c_penalty, c_penalty_ex);

        let alloclen = (sequences[0].len() + 1000) as c_int;
        let icyc = group1.len() as c_int;
        let jcyc = group2.len() as c_int;

        // Each C string needs a buffer of at least alloclen+1 bytes for in-place editing.
        let mut buf1: Vec<Vec<u8>> = c_seq1.iter().map(|s| {
            let mut v = s.as_bytes().to_vec();
            v.resize(alloclen as usize + 1, 0);
            v
        }).collect();
        let mut buf2: Vec<Vec<u8>> = c_seq2.iter().map(|s| {
            let mut v = s.as_bytes().to_vec();
            v.resize(alloclen as usize + 1, 0);
            v
        }).collect();
        let mut p1: Vec<*mut c_char> = buf1.iter_mut().map(|v| v.as_mut_ptr() as *mut c_char).collect();
        let mut p2: Vec<*mut c_char> = buf2.iter_mut().map(|v| v.as_mut_ptr() as *mut c_char).collect();
        let mut e1: Vec<c_double> = w1n.clone();
        let mut e2: Vec<c_double> = w2n.clone();

        // n_dynamicmtx as f64 matrix from scoring.substitution_matrix
        let nalpha = scoring.substitution_matrix.len() as c_int;
        let n_dyn = mafft_sys::AllocateDoubleMtx(nalpha, nalpha);
        for i in 0..scoring.substitution_matrix.len() {
            for j in 0..scoring.substitution_matrix[i].len() {
                *(*n_dyn.add(i)).add(j) = scoring.substitution_matrix[i][j] as f64;
            }
        }

        let mut impmatch = 0.0f64;

        // Pass all-'o' sgap/egap to mirror C's Falign single-segment path.
        let mut sgap1: Vec<u8> = vec![b'o'; group1.len() + 1];
        let mut sgap2: Vec<u8> = vec![b'o'; group2.len() + 1];
        let mut egap1: Vec<u8> = vec![b'o'; group1.len() + 1];
        let mut egap2: Vec<u8> = vec![b'o'; group2.len() + 1];
        let use_sgap = std::env::var("BL50_TEST_USE_SGAP").is_ok();
        let (sgap1_p, sgap2_p, egap1_p, egap2_p) = if use_sgap {
            (sgap1.as_mut_ptr() as *mut c_char,
             sgap2.as_mut_ptr() as *mut c_char,
             egap1.as_mut_ptr() as *mut c_char,
             egap2.as_mut_ptr() as *mut c_char)
        } else {
            (std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut())
        };

        let c_score = mafft_sys::A__align(
            n_dyn, c_penalty, c_penalty_ex,
            p1.as_mut_ptr(), p2.as_mut_ptr(),
            e1.as_mut_ptr(), e2.as_mut_ptr(),
            icyc, jcyc, alloclen,
            0, &mut impmatch,
            sgap1_p, sgap2_p, egap1_p, egap2_p,
            std::ptr::null_mut(), 0, std::ptr::null_mut(),
            0, 0,  // headgp, tailgp = 0 (matching engine.rs FftAlignParams head_gap/tail_gap = false)
            -1, -1,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            0.0, 0.0,
        );

        // Find result width by reading back p1[0]
        let c_width = {
            let s = p1[0];
            let mut k = 0; while *s.add(k) != 0 { k += 1; } k
        };
        eprintln!("C A__align: width={} score={:.2}", c_width, c_score);

        // Save first 80 chars of C's first sequence
        let c_s1: Vec<u8> = (0..c_width.min(80)).map(|k| *p1[0].add(k) as u8).collect();
        let c_s1_str = std::str::from_utf8(&c_s1).unwrap_or("");
        eprintln!("C    s1[..80] = {c_s1_str}");

        // Reconstruct Rust's first sequence's first 80 chars from operations
        // (Rust's profile_align doesn't materialize sequences; need to apply ops)
        let mut rust_s1 = Vec::new();
        let mut p = 0;
        for op in rust_aln.operations.iter().take(80) {
            match op {
                mafft_align::AlignOp::Match | mafft_align::AlignOp::Delete => {
                    if p < sequences[group1[0]].len() {
                        rust_s1.push(sequences[group1[0]][p]);
                        p += 1;
                    }
                }
                mafft_align::AlignOp::Insert => rust_s1.push(b'-'),
            }
        }
        let rust_s1_str = std::str::from_utf8(&rust_s1).unwrap_or("");
        eprintln!("Rust s1[..80] = {rust_s1_str}");

        mafft_sys::A__align(std::ptr::null_mut(), 0, 0, std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), 0, 0, 0, 0, std::ptr::null_mut(),
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(),
            std::ptr::null_mut(), 0, std::ptr::null_mut(), 0, 0, -1, -1,
            std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), 0.0, 0.0);
        mafft_sys::freeconstants();
    }
}
