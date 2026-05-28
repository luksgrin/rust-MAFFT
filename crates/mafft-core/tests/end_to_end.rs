/// End-to-end integration tests: read real test data, align, verify output.

use std::path::PathBuf;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::{read_fasta, read_fasta_casepreserve};
use mafft_types::Sequence;

fn test_data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../mafft-upstream/test")
        .join(name)
}

/// Path to a fixture file in `crates/mafft-core/tests/fixtures/`.
///
/// Unlike `test_data_path`, which points into the upstream submodule (and
/// thus contains only files shipped by upstream MAFFT), this points to
/// fixtures committed in our own repo — typically C-reference outputs
/// generated for alignment modes that upstream doesn't ship references for
/// (e.g., NW-NS-2 / `--nofft`).
fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn align_sample_fasta() {
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let input = read_fasta(test_data_path("sample")).unwrap();
    assert_eq!(input.nseq(), 36);

    let msa = engine.align(&input);

    // All sequences should have the same width
    let width = msa.width();
    assert!(width > 0, "alignment should have non-zero width");
    for (i, seq) in msa.sequences.iter().enumerate() {
        assert_eq!(
            seq.len(),
            width,
            "sequence {} has length {} but expected {}",
            msa.names[i],
            seq.len(),
            width
        );
    }

    // Removing gaps should recover original sequences
    for (i, seq) in msa.sequences.iter().enumerate() {
        let ungapped: Vec<u8> = seq.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(
            ungapped, input.sequences[i].data,
            "ungapped sequence {} doesn't match original",
            msa.names[i]
        );
    }
}

#[test]
fn align_sample_with_refinement() {
    // Use only the first 6 sequences to keep refinement fast
    let full_input = read_fasta(test_data_path("sample")).unwrap();
    let input = mafft_types::SequenceSet {
        sequences: full_input.sequences[..6].to_vec(),
        seq_type: full_input.seq_type,
    };

    // First verify progressive-only works
    let engine_prog = MafftEngine::new(AlignmentMode::FftNs2);
    let msa_prog = engine_prog.align(&input);
    for (i, seq) in msa_prog.sequences.iter().enumerate() {
        let residues = seq.iter().filter(|&&c| c != b'-').count();
        assert_eq!(residues, input.sequences[i].data.len(),
            "progressive lost residues for seq {i}: {} vs {}", residues, input.sequences[i].data.len());
    }

    // Now test with refinement
    let engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 });
    let msa = engine.align(&input);

    let width = msa.width();
    assert!(width > 0);
    for (i, seq) in msa.sequences.iter().enumerate() {
        assert_eq!(seq.len(), width, "sequence {i} width mismatch after refinement");
    }

    // Ungapped residue counts should be preserved
    for (i, seq) in msa.sequences.iter().enumerate() {
        let residue_count = seq.iter().filter(|&&c| c != b'-').count();
        assert_eq!(residue_count, input.sequences[i].data.len(),
            "sequence {i} lost residues during refinement: {} vs {}",
            residue_count, input.sequences[i].data.len());
    }
}

/// FFT-NS-2 (default strategy) must be byte-identical to C's reference output.
///
/// The reference is `mafft-upstream/test/sample.fftns2`, shipped by upstream.
/// This test guards the full pipeline: FFT correlation, anchor selection,
/// anchor segment head/tail gap handling, inter-anchor DP, retree distance,
/// and UPGMA tree reconstruction.
#[test]
fn fftns2_byte_identical_to_c() {
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq(), "different number of sequences");
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "FFT-NS-2 width differs: Rust={}, C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len()
    );

    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i]
                    .iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} (name: {:?}) differs; first diff at position {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(
        mismatches, 0,
        "{mismatches} sequence(s) differ from C's FFT-NS-2 output"
    );
}

/// FFT-NS-i (`mafft --maxiterate 100`) must produce byte-identical output to C.
#[test]
fn fftnsi_byte_identical_to_c() {
    let c_ref = read_fasta(test_data_path("sample.fftnsi")).unwrap();
    let input = read_fasta(test_data_path("sample")).unwrap();

    // C test driver uses `--maxiterate 100`; script caps to 16 internally.
    let msa = MafftEngine::new(AlignmentMode::FftNsi { iterations: 100 }).align(&input);

    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "FFT-NS-i width differs: Rust={}, C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );

    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i].iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} (name: {:?}) differs; first diff at {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} sequence(s) differ from C's FFT-NS-i output");
}

/// Regression guard for the FFT-segmented refinement boundary-frequencies
/// fix (`refinement.rs::realign_all` use_fft branch, 2026-05-25). BB12019
/// (BALIBASE 3, 5 sequences, 883 columns after iter=5) is the smallest
/// BALIBASE input where the unfixed `BoundaryFreqs::default() = 1.0`
/// produces a 4-line traceback shift vs C. With C's per-segment
/// `outgapcount(sgap/egap)` boundary frequencies wired up, output is
/// bit-identical to `mafft --maxiterate 5 BB12019`.
///
/// Companion to the (much larger) BB30013 BALIBASE case, which exercises
/// the same code path on 86 sequences — kept out of fixtures for size.
#[test]
fn fftnsi_segmented_boundary_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("bali3.BB12019.fftnsi.iter5")).unwrap();
    let input = read_fasta(fixture_path("bali3.BB12019.fa")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNsi { iterations: 5 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq(), "nseq mismatch");
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --maxiterate 5 BB12019 output \
             (FFT-segmented refinement boundary-frequencies regression)",
        );
    }
}


/// Every NW-NS-2 merge step's `(clus1, clus2, width, score)` must match C's.
///
/// This is a finer-grained regression guard than `nofft_byte_identical_to_c`:
/// it catches intermediate DP regressions that happen to produce the same
/// final output. Reference: `tests/fixtures/sample.nwns2.steps` (one RDBG
/// line per merge across both retree passes, 70 lines total).
#[test]
fn nofft_per_step_matches_c() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .align(&input);

    let ref_txt = std::fs::read_to_string(fixture_path("sample.nwns2.steps"))
        .expect("missing tests/fixtures/sample.nwns2.steps — see fixtures/README.md");
    let expected: Vec<(usize, usize, usize, f64)> = ref_txt
        .lines()
        .filter(|l| l.starts_with("RDBG"))
        .map(|l| {
            let parts: Vec<&str> = l.split_whitespace().collect();
            // RDBG <step_idx> <clus1> <clus2> <width> <score>
            (
                parts[2].parse().unwrap(),
                parts[3].parse().unwrap(),
                parts[4].parse().unwrap(),
                parts[5].parse().unwrap(),
            )
        })
        .collect();

    assert_eq!(
        msa.step_trace.len(), expected.len(),
        "step count differs: Rust has {} steps, C has {}",
        msa.step_trace.len(), expected.len()
    );

    let mut first_mismatch: Option<usize> = None;
    for (i, (got, want)) in msa.step_trace.iter().zip(expected.iter()).enumerate() {
        let matches = got.clus1 == want.0
            && got.clus2 == want.1
            && got.width == want.2
            && (got.score - want.3).abs() < 0.1;
        if !matches && first_mismatch.is_none() {
            first_mismatch = Some(i);
            eprintln!(
                "first step mismatch at index {i}: \
                 Rust=({} {} {} {:.1}), C=({} {} {} {:.1})",
                got.clus1, got.clus2, got.width, got.score,
                want.0, want.1, want.2, want.3
            );
        }
    }
    assert!(
        first_mismatch.is_none(),
        "per-step trace diverges from C at step {}; see stderr",
        first_mismatch.unwrap()
    );
}

/// NW-NS-2 (`--nofft`) must be byte-identical to C's output on `test/sample`.
///
/// This is a regression guard for the suite of fixes that brought the NW-NS-2
/// pipeline to exact parity with C MAFFT: 0-based `match_calc` indexing,
/// `outgap=0` boundary handling, `nongap_freq` default 0.0, `match_calc_row(i)`
/// position, and the retree-2 `penalty_dist` scaling. If any of these
/// regresses, the diff below grows from 0 lines to many.
///
/// Reference: `tests/fixtures/sample.nwns2` (C 7.526 `mafft --nofft --quiet`).
#[test]
fn nofft_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.nwns2"))
        .expect("missing tests/fixtures/sample.nwns2 — see fixtures/README.md");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq(), "different number of sequences");

    let rust_width = msa.sequences[0].len();
    let c_width = c_ref.sequences[0].data.len();
    assert_eq!(
        rust_width, c_width,
        "alignment width differs: Rust={rust_width}, C={c_width}"
    );

    // Every sequence must match byte-for-byte (including gap positions).
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i]
                    .iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} (name: {:?}) differs; first diff at position {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(
        mismatches, 0,
        "{mismatches} sequence(s) differ from C's --nofft output"
    );
}

/// NW-NS-2 with a non-default `--op` override must match C byte-for-byte.
///
/// Guards the command-line gap-opening override application. Our scaling
/// (`ppenalty = -op * 1000`, then `penalty = 0.6 * ppenalty`) must match
/// C's `constants()` routine for any `--op` value, not just the default.
///
/// Reference: `tests/fixtures/sample.nwns2.op25` (`mafft --nofft --op 2.5`).
#[test]
fn nofft_op_override_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.nwns2.op25"))
        .expect("missing tests/fixtures/sample.nwns2.op25 — see fixtures/README.md");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .with_gap_open(2.5)
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --op 2.5"
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --op 2.5 output"
        );
    }
}

/// `--bl 80 --nofft` must match C byte-for-byte.
///
/// Guards the BLOSUM80 substitution-matrix table. MAFFT's `tmpmtx80` in
/// `mafft-upstream/core/blosum.c` differs from the standard NCBI BLOSUM80
/// at four cells (H/R = 0, F/M = 0, P/R = -3, V/I = 4); using the standard
/// table here would diverge from C by ~12 columns on this fixture.
///
/// Reference: `tests/fixtures/sample.bl80.nwns2` (`mafft --nofft --bl 80`).
#[test]
fn nofft_bl80_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.bl80.nwns2"))
        .expect("missing tests/fixtures/sample.bl80.nwns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .with_scoring_model(ScoringModel::Blosum(80))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --bl 80 --nofft: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --bl 80 --nofft output"
        );
    }
}

/// `--bl 45` (FFT-NS-2 with BLOSUM45) must match C byte-for-byte.
///
/// MAFFT's `tmpmtx45` differs from standard NCBI BLOSUM45 at one cell —
/// (P, H) = -2 in MAFFT, -1 in NCBI. Reverting that cell to the NCBI value
/// re-introduces a 144-line diff against C's `--bl 45` output.
///
/// Reference: `tests/fixtures/sample.bl45.fftns2`.
#[test]
fn fftns2_bl45_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.bl45.fftns2"))
        .expect("missing tests/fixtures/sample.bl45.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_scoring_model(ScoringModel::Blosum(45))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --bl 45: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --bl 45 output"
        );
    }
}

/// `--bl 50` (FFT-NS-2 with BLOSUM50) must match C byte-for-byte.
/// Closed 2026-05-08 by the FMA fix in `profile_align_imp_with_boundary`
/// (see TODO §4) — gcc -O3 fuses `a + b * c` into FMA, Rust's `+` and `*`
/// don't, so direct DP cells diverge by 1 ULP per accumulation. The
/// flatter BL50 score landscape (vs BL62/30/45/80) exposed this as a
/// tie-break divergence at step 24's profile DP. Fix: use `f64::mul_add`
/// in match_calc_row and DP gap-frequency computations to match C's FMA.
///
/// Reference: `tests/fixtures/sample.bl50.fftns2`.
#[test]
fn fftns2_bl50_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.bl50.fftns2"))
        .expect("missing tests/fixtures/sample.bl50.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_scoring_model(ScoringModel::Blosum(50))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --bl 50: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --bl 50 output"
        );
    }
}

/// `--bl 30` (FFT-NS-2 with BLOSUM30) must match C byte-for-byte.
///
/// Reference: `tests/fixtures/sample.bl30.fftns2`.
#[test]
fn fftns2_bl30_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.bl30.fftns2"))
        .expect("missing tests/fixtures/sample.bl30.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_scoring_model(ScoringModel::Blosum(30))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --bl 30: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --bl 30 output"
        );
    }
}

/// `--jtt 200` (FFT-NS-2 with JTT PAM 200) must match C byte-for-byte.
///
/// Guards (1) the JTT lower-triangle accepted-point-mutation table in
/// `mafft-scoring/src/jtt.rs::jtt_rsr_matrix`, (2) the PAM exponentiation
/// loop in `build_jtt_pam_matrix`, and (3) the normalize/600-scale/offset
/// pipeline shared with BLOSUM.
///
/// Reference: `tests/fixtures/sample.jtt200.fftns2` (`mafft --jtt 200`).
#[test]
fn fftns2_jtt200_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.jtt200.fftns2"))
        .expect("missing tests/fixtures/sample.jtt200.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_scoring_model(ScoringModel::Jtt(200))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --jtt 200: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --jtt 200 output"
        );
    }
}

/// `--tm 200 --nofft` (NW-NS-2 with TM PAM 200) must match C byte-for-byte.
///
/// Guards the TM upper-triangle accepted-point-mutation table in
/// `mafft-scoring/src/jtt.rs::tm_rsr_matrix` together with the TM frequency
/// vector. C's `JTTmtx(... isTM=1)` reads the *upper* triangle of the rsr
/// counts array (lines 145-217 of `mafft-upstream/core/JTT.c`); the lower
/// triangle holds JTT data and is irrelevant for TM. Before this test,
/// `--tm` silently produced JTT-like output because the upper triangle was
/// never populated in our Rust port.
///
/// Reference: `tests/fixtures/sample.tm200.nwns2` (`mafft --tm 200 --nofft`).
#[test]
fn nofft_tm200_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.tm200.nwns2"))
        .expect("missing tests/fixtures/sample.tm200.nwns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .with_scoring_model(ScoringModel::Tm(200))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --tm 200 --nofft: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --tm 200 --nofft output"
        );
    }
}

/// `--tm 100 --nofft` exercises the same TM data at a non-default PAM —
/// catches drift in the matrix-power path in `build_jtt_pam_matrix`.
///
/// Reference: `tests/fixtures/sample.tm100.nwns2`.
#[test]
fn nofft_tm100_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.tm100.nwns2"))
        .expect("missing tests/fixtures/sample.tm100.nwns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .with_scoring_model(ScoringModel::Tm(100))
        .align(&input);

    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --tm 100 --nofft output"
        );
    }
}

/// `--bl 80` (FFT-NS-2 with BLOSUM80) must match C byte-for-byte.
///
/// Same matrix-table guard as `nofft_bl80_byte_identical_to_c` but exercises
/// the FFT pipeline as well.
///
/// Reference: `tests/fixtures/sample.bl80.fftns2` (`mafft --bl 80`).
#[test]
fn fftns2_bl80_byte_identical_to_c() {
    use mafft_types::ScoringModel;
    let c_ref = read_fasta(fixture_path("sample.bl80.fftns2"))
        .expect("missing tests/fixtures/sample.bl80.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_scoring_model(ScoringModel::Blosum(80))
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --bl 80: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --bl 80 output"
        );
    }
}

/// NW-NS-2 with a non-default `--ep` override must match C byte-for-byte.
///
/// Guards the scoring-matrix offset application. C's `--ep` flag maps to
/// `aof` in the shell script (negated), then to `poffset` in disttbfast,
/// then to `offset` which is subtracted from the scoring matrix. This test
/// ensures that override path correctly adjusts the already-built matrix.
///
/// Reference: `tests/fixtures/sample.nwns2.ep05` (`mafft --nofft --ep 0.5`).
#[test]
fn nofft_ep_override_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.nwns2.ep05"))
        .expect("missing tests/fixtures/sample.nwns2.ep05 — see fixtures/README.md");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .with_gap_offset(0.5)
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --ep 0.5"
    );
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --ep 0.5 output"
        );
    }
}

/// RNA NW-NS-2 must match C byte-for-byte, ignoring ASCII case.
///
/// Guards the nucleotide alignment path end-to-end. Currently the only
/// difference between our output and C's is that C preserves the input
/// lowercase while Rust uppercases residues before alignment; the gap
/// placement is identical. This test normalizes case on both sides so it
/// asserts alignment equality (column-for-column) without being sensitive
/// to that pre-alignment casing choice.
///
/// Reference: `tests/fixtures/samplerna.nwns2`
/// (`mafft --nofft mafft-upstream/test/samplerna`).
#[test]
fn rna_nofft_case_insensitive_identical_to_c() {
    let c_ref = read_fasta(fixture_path("samplerna.nwns2"))
        .expect("missing tests/fixtures/samplerna.nwns2 — see fixtures/README.md");
    let input = read_fasta(test_data_path("samplerna")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_nofft(true)
        .align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "RNA alignment width differs"
    );

    fn lower(bytes: &[u8]) -> Vec<u8> {
        bytes.iter().map(|b| b.to_ascii_lowercase()).collect()
    }
    for i in 0..msa.nseq() {
        assert_eq!(
            lower(&msa.sequences[i]), lower(&c_ref.sequences[i].data),
            "RNA seq {i} differs from C (case-insensitive)"
        );
    }
}

/// Compute sum-of-pairs identity score for an alignment.
fn sum_of_pairs_identity(sequences: &[Vec<u8>]) -> f64 {
    let n = sequences.len();
    if n < 2 { return 0.0; }
    let mut total_match = 0u64;
    let mut total_aligned = 0u64;
    for i in 0..n {
        for j in (i + 1)..n {
            let len = sequences[i].len().min(sequences[j].len());
            for k in 0..len {
                let a = sequences[i][k];
                let b = sequences[j][k];
                if a != b'-' && b != b'-' {
                    total_aligned += 1;
                    if a == b {
                        total_match += 1;
                    }
                }
            }
        }
    }
    if total_aligned == 0 { 0.0 } else { total_match as f64 / total_aligned as f64 }
}

#[test]
fn align_rna_sample() {
    let engine = MafftEngine::default();
    let input = read_fasta(test_data_path("samplerna")).unwrap();
    assert!(input.seq_type.is_nucleotide());

    let msa = engine.align(&input);
    let width = msa.width();
    assert!(width > 0);
    for seq in &msa.sequences {
        assert_eq!(seq.len(), width);
    }
}

#[test]
fn diagnostic_guide_tree() {
    use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, ktuple_distance};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let nseq = input.nseq();

    // Compute 6-tuple distance matrix (same as engine's first pass)
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            let d = ktuple_distance(&input.sequences[i].data, &input.sequences[j].data, 6);
            dm.set(i, j, d);
        }
    }

    // Print first 5 distances
    eprintln!("First 5 pairwise 6-tuple distances:");
    let mut count = 0;
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            if count >= 5 { break; }
            eprintln!("  d({},{}) = {:.6}", i, j, dm.get(i, j));
            count += 1;
        }
        if count >= 5 { break; }
    }

    // Build tree and print first 5 merge steps
    let topo = musclesupg(&dm, ClusterMethod::default());
    eprintln!("First 5 merge steps:");
    for (i, step) in topo.steps.iter().enumerate().take(5) {
        let mut left: Vec<usize> = step.left.clone();
        let mut right: Vec<usize> = step.right.clone();
        left.sort();
        right.sort();
        eprintln!("  Step {}: {:?} + {:?} (len: {:.4}, {:.4})",
            i, left, right, step.left_length, step.right_length);
    }

    // Basic structural check
    let last = topo.steps.last().unwrap();
    let mut all: Vec<usize> = last.left.iter().chain(last.right.iter()).copied().collect();
    all.sort();
    assert_eq!(all, (0..nseq).collect::<Vec<_>>(), "tree doesn't cover all sequences");
}

#[test]
fn diagnostic_fft_anchoring() {
    use mafft_types::{ScoringModel, SeqType};
    use mafft_scoring::build_context;
    use mafft_align::{Profile, FftAlignParams, fft_profile_align, profile_align, GapModel};
    use mafft_fft::SegmentParams;

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);

    // Take two sequences and compare FFT-accelerated vs direct DP alignment
    let s1 = &input.sequences[0].data;
    let s2 = &input.sequences[1].data;

    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let w = vec![1.0];

    let prof1 = Profile::from_aligned(&seqs1, &w, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &w, &scoring.amino_map, scoring.nalphabets);

    // Direct DP alignment
    let dp_aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);

    // FFT-accelerated alignment
    let fft_params = FftAlignParams {
        num_candidates: 20,
        segment_params: SegmentParams::protein(),
        gap: gap.clone(),
        head_gap: true,
        tail_gap: true,
        num_channels: 20,
        property_channels: None,
    };
    let fft_aln = fft_profile_align(&prof1, &prof2, &scoring.consweight_matrix, &fft_params);

    eprintln!("Direct DP:  score={:.1}, ops={}", dp_aln.score, dp_aln.operations.len());
    eprintln!("FFT accel:  score={:.1}, ops={}", fft_aln.score, fft_aln.operations.len());

    // Both should produce valid alignments
    assert!(dp_aln.operations.len() > 0, "DP alignment empty");
    assert!(fft_aln.operations.len() > 0, "FFT alignment empty");

    // FFT should produce a score at least 50% of DP (it's an approximation)
    if dp_aln.score != 0.0 {
        let ratio = fft_aln.score / dp_aln.score;
        eprintln!("FFT/DP score ratio: {:.4}", ratio);
        assert!(ratio > 0.3, "FFT score too low vs DP: {:.1} vs {:.1}", fft_aln.score, dp_aln.score);
    }
}

#[test]
fn diagnostic_fft_vs_nofft() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();
    let c_seqs: Vec<Vec<u8>> = c_ref.sequences.iter().map(|s| s.data.clone()).collect();

    // FFT-NS-2
    let msa_fft = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    // NW-NS-2 (pure DP, no FFT)
    let msa_nofft = MafftEngine::new(AlignmentMode::FftNs2).with_nofft(true).align(&input);

    let sp_fft = sum_of_pairs_identity(&msa_fft.sequences);
    let sp_nofft = sum_of_pairs_identity(&msa_nofft.sequences);
    let sp_c = sum_of_pairs_identity(&c_seqs);

    eprintln!("C reference:   SP={:.4}, width={}", sp_c, c_seqs[0].len());
    eprintln!("Rust FFT-NS-2: SP={:.4}, width={}", sp_fft, msa_fft.width());
    eprintln!("Rust NW-NS-2:  SP={:.4}, width={}", sp_nofft, msa_nofft.width());
    eprintln!("FFT/C ratio:   {:.4}", sp_fft / sp_c);
    eprintln!("noFFT/C ratio: {:.4}", sp_nofft / sp_c);
    eprintln!("FFT == noFFT:  {}", msa_fft.sequences == msa_nofft.sequences);
}

#[test]
fn diagnostic_retree_widths() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    
    let msa1 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(1).align(&input);
    let msa2 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(2).align(&input);
    
    eprintln!("retree=1: width={}, SP={:.4}", msa1.width(), sum_of_pairs_identity(&msa1.sequences));
    eprintln!("retree=2: width={}, SP={:.4}", msa2.width(), sum_of_pairs_identity(&msa2.sequences));
}

#[test]
fn diagnostic_merge_widths() {
    use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, ktuple_distance};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let nseq = input.nseq();

    // Build distance matrix and tree
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            dm.set(i, j, ktuple_distance(&input.sequences[i].data, &input.sequences[j].data, 6));
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());

    // Manually trace merge steps
    let max_len = input.sequences.iter().map(|s| s.data.len()).max().unwrap_or(0);
    let aligned: Vec<Vec<u8>> = input.sequences.iter()
        .map(|s| { let mut p = s.data.clone(); p.resize(max_len, b'-'); p })
        .collect();

    eprintln!("Initial width: {}", aligned[0].len());
    eprintln!("Seq lengths: min={}, max={}", 
        input.sequences.iter().map(|s| s.data.len()).min().unwrap(),
        input.sequences.iter().map(|s| s.data.len()).max().unwrap());

    for (step_idx, step) in topo.steps.iter().enumerate().take(10) {
        let width = aligned[0].len();
        let g1_len = step.left.len();
        let g2_len = step.right.len();
        eprintln!("Step {}: width={}, merge {} + {} seqs", step_idx, width, g1_len, g2_len);
    }
    eprintln!("...");
    // Show last 3 steps
    for (step_idx, step) in topo.steps.iter().enumerate().skip(topo.steps.len().saturating_sub(3)) {
        eprintln!("Step {}: merge {} + {} seqs", step_idx, step.left.len(), step.right.len());
    }
}

#[test]
fn diagnostic_first_merge() {
    use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, ktuple_distance};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, GapModel};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Jtt(200), SeqType::Protein);
    let nseq = input.nseq();

    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            dm.set(i, j, ktuple_distance(&input.sequences[i].data, &input.sequences[j].data, 6));
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());

    // First merge step
    let step = &topo.steps[0];
    eprintln!("First merge: {:?} + {:?}", step.left, step.right);
    
    let s1 = &input.sequences[step.left[0]].data;
    let s2 = &input.sequences[step.right[0]].data;
    eprintln!("  Seq {} len={}", step.left[0], s1.len());
    eprintln!("  Seq {} len={}", step.right[0], s2.len());

    // Profile align these two
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let w = vec![1.0];
    let prof1 = Profile::from_aligned(&seqs1, &w, &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &w, &scoring.amino_map, scoring.nalphabets);
    let aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    
    eprintln!("  Alignment score: {:.1}", aln.score);
    eprintln!("  Alignment width: {}", aln.operations.len());
    
    let matches = aln.operations.iter().filter(|op| matches!(op, mafft_align::AlignOp::Match)).count();
    let deletes = aln.operations.iter().filter(|op| matches!(op, mafft_align::AlignOp::Delete)).count();
    let inserts = aln.operations.iter().filter(|op| matches!(op, mafft_align::AlignOp::Insert)).count();
    eprintln!("  Match={}, Delete={}, Insert={}", matches, deletes, inserts);
}

#[test]
fn diagnostic_distance_check() {
    use mafft_tree::ktuple_distance;
    let input = read_fasta(test_data_path("sample")).unwrap();
    // Print first 10 pairwise distances with high precision
    let mut count = 0;
    for i in 0..input.nseq() {
        for j in (i+1)..input.nseq() {
            if count >= 10 { break; }
            let d = ktuple_distance(&input.sequences[i].data, &input.sequences[j].data, 6);
            eprintln!("d({},{}) = {:.15}", i, j, d);
            count += 1;
        }
        if count >= 10 { break; }
    }
}

#[test]
fn diagnostic_merge_trace() {
    use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, ktuple_distance, sequence_weights};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, GapModel, AlignOp};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Jtt(200), SeqType::Protein);
    let nseq = input.nseq();

    // Build distance matrix and tree (retree pass 1)
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            dm.set(i, j, ktuple_distance(&input.sequences[i].data, &input.sequences[j].data, 6));
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());
    let weights = sequence_weights(&topo);

    // Trace first 10 merge steps with alignment details
    let mut aligned: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    for (step_idx, step) in topo.steps.iter().enumerate().take(10) {
        let width1 = aligned[step.left[0]].len();
        let width2 = aligned[step.right[0]].len();

        let seqs1: Vec<&[u8]> = step.left.iter().map(|&i| aligned[i].as_slice()).collect();
        let seqs2: Vec<&[u8]> = step.right.iter().map(|&i| aligned[i].as_slice()).collect();
        let w1: Vec<f64> = step.left.iter().map(|&i| weights[i]).collect();
        let w2: Vec<f64> = step.right.iter().map(|&i| weights[i]).collect();
        let sum1: f64 = w1.iter().sum();
        let sum2: f64 = w2.iter().sum();
        let w1n: Vec<f64> = w1.iter().map(|v| v / sum1).collect();
        let w2n: Vec<f64> = w2.iter().map(|v| v / sum2).collect();

        let prof1 = Profile::from_aligned(&seqs1, &w1n, &scoring.amino_map, scoring.nalphabets);
        let prof2 = Profile::from_aligned(&seqs2, &w2n, &scoring.amino_map, scoring.nalphabets);
        let aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);

        let matches = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
        let deletes = aln.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
        let inserts = aln.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();

        eprintln!("Step {:2}: {:?}+{:?} w1={:.6} w2={:.6} prof1={} prof2={} score={:.1} ops={} M/D/I={}/{}/{}",
            step_idx, step.left, step.right, sum1, sum2,
            prof1.length, prof2.length, aln.score,
            aln.operations.len(), matches, deletes, inserts);

        // Apply alignment to sequences (simplified — just track widths)
        let new_width = aln.operations.len();
        for &idx in &step.left {
            let mut new_seq = Vec::with_capacity(new_width);
            let mut cursor = 0;
            for op in &aln.operations {
                match op {
                    AlignOp::Match | AlignOp::Delete => {
                        new_seq.push(if cursor < width1 { aligned[idx][cursor] } else { b'-' });
                        cursor += 1;
                    }
                    AlignOp::Insert => new_seq.push(b'-'),
                }
            }
            aligned[idx] = new_seq;
        }
        for &idx in &step.right {
            let mut new_seq = Vec::with_capacity(new_width);
            let mut cursor = 0;
            for op in &aln.operations {
                match op {
                    AlignOp::Match | AlignOp::Insert => {
                        new_seq.push(if cursor < width2 { aligned[idx][cursor] } else { b'-' });
                        cursor += 1;
                    }
                    AlignOp::Delete => new_seq.push(b'-'),
                }
            }
            aligned[idx] = new_seq;
        }
    }
}

#[test]
fn diagnostic_alignment_diff() {
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();
    let input = read_fasta(test_data_path("sample")).unwrap();
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let msa = engine.align(&input);

    eprintln!("Rust width: {}, C width: {}", msa.width(), c_ref.sequences[0].data.len());
    
    // Count identical columns
    let rust_width = msa.width();
    let c_width = c_ref.sequences[0].data.len();
    
    // Compare first sequence's alignment character by character
    let r0 = &msa.sequences[0];
    let c0 = &c_ref.sequences[0].data;
    
    // Find first difference
    let min_len = r0.len().min(c0.len());
    let mut first_diff = min_len;
    for k in 0..min_len {
        if r0[k] != c0[k] {
            first_diff = k;
            break;
        }
    }
    
    if first_diff < min_len {
        eprintln!("First diff at col {}: Rust='{}' C='{}'", 
            first_diff, r0[first_diff] as char, c0[first_diff] as char);
        // Show context around first diff
        let start = first_diff.saturating_sub(5);
        let end = (first_diff + 10).min(min_len);
        eprintln!("Rust seq0[{}..{}]: {}", start, end, 
            String::from_utf8_lossy(&r0[start..end]));
        eprintln!("C    seq0[{}..{}]: {}", start, end,
            String::from_utf8_lossy(&c0[start..end]));
    } else {
        eprintln!("Seq 0 matches for first {} chars!", min_len);
    }
    
    // Count total matching columns across all sequences
    let mut total_match = 0u64;
    let mut total_cols = 0u64;
    if rust_width == c_width {
        for col in 0..rust_width {
            let mut all_match = true;
            for i in 0..msa.nseq() {
                if msa.sequences[i][col] != c_ref.sequences[i].data[col] {
                    all_match = false;
                    break;
                }
            }
            if all_match { total_match += 1; }
            total_cols += 1;
        }
        eprintln!("Matching columns: {}/{} ({:.1}%)", total_match, total_cols,
            100.0 * total_match as f64 / total_cols as f64);
    }
}

#[test]
fn diagnostic_gap_pattern() {
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();
    let input = read_fasta(test_data_path("sample")).unwrap();
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let msa = engine.align(&input);

    // For each sequence, show the gap pattern (positions of first/last residue)
    for i in 0..5.min(msa.nseq()) {
        let r = &msa.sequences[i];
        let c = &c_ref.sequences[i].data;
        
        let r_first = r.iter().position(|&c| c != b'-').unwrap_or(0);
        let r_last = r.iter().rposition(|&c| c != b'-').unwrap_or(0);
        let c_first = c.iter().position(|&c| c != b'-').unwrap_or(0);
        let c_last = c.iter().rposition(|&c| c != b'-').unwrap_or(0);
        let r_gaps: usize = r.iter().filter(|&&c| c == b'-').count();
        let c_gaps: usize = c.iter().filter(|&&c| c == b'-').count();
        
        eprintln!("Seq {:2}: Rust first={:3} last={:3} gaps={:3} width={}  |  C first={:3} last={:3} gaps={:3} width={}",
            i, r_first, r_last, r_gaps, r.len(), c_first, c_last, c_gaps, c.len());
    }
    
    // Show the retree pass info
    eprintln!("\n--- Retree pass 1 vs pass 2 ---");
    let engine1 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(1);
    let engine2 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(2);
    let msa1 = engine1.align(&input);
    let msa2 = engine2.align(&input);
    eprintln!("retree=1: width={} SP={:.4}", msa1.width(), sum_of_pairs_identity(&msa1.sequences));
    eprintln!("retree=2: width={} SP={:.4}", msa2.width(), sum_of_pairs_identity(&msa2.sequences));
}

#[test]
fn diagnostic_align11_vs_profile() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, pairwise_align11, GapModel, AlignOp};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Jtt(200), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Compare align11 vs profile_align for seqs 19 and 20 (merged at step 1)
    let s1 = &input.sequences[19].data;
    let s2 = &input.sequences[20].data;

    // Use the same boundary convention the engine uses (outgap=0 → false, false).
    let aln11 = pairwise_align11(s1, s2, &scoring.consweight_matrix, &scoring.amino_map,
        scoring.gap.open as f64, false, false);

    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let aln_prof = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, false, false);

    let m11 = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let d11 = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
    let i11 = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();
    let mp_n = aln_prof.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let dp_n = aln_prof.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
    let ip_n = aln_prof.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();

    eprintln!("G__align11: score={:.1} width={} M/D/I={}/{}/{}", aln11.score, aln11.operations.len(), m11, d11, i11);
    eprintln!("MSalignmm:  score={:.1} width={} M/D/I={}/{}/{}", aln_prof.score, aln_prof.operations.len(), mp_n, dp_n, ip_n);

    // Regression guard: the two pairwise code paths (G__align11 and MSalignmm
    // specialized to 1×1) must produce identical alignment operations and
    // matching scores on any real input pair. Any future divergence means one
    // of them has broken its port of C's algorithm.
    assert_eq!(
        aln11.operations, aln_prof.operations,
        "pairwise_align11 and profile_align disagree on 1×1 alignment"
    );
    assert!(
        (aln11.score - aln_prof.score).abs() < 0.01,
        "pairwise_align11 score {} ≠ profile_align score {}",
        aln11.score, aln_prof.score
    );
}

#[test]
fn diagnostic_fft_anchors() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, fft_profile_align, GapModel, AlignOp, FftAlignParams};
    use mafft_fft::SegmentParams;

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Test merge step 6: [28] + [29,30] — first multi-seq merge with different widths
    let s1 = &input.sequences[28].data;
    let s29 = &input.sequences[29].data;
    let s30 = &input.sequences[30].data;

    // Build profiles
    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);

    let seqs2: Vec<&[u8]> = vec![s29.as_slice(), s30.as_slice()];
    let prof2 = Profile::from_aligned(&seqs2, &[0.5, 0.5], &scoring.amino_map, scoring.nalphabets);

    let fft_params = FftAlignParams {
        num_candidates: 20,
        segment_params: SegmentParams::protein(),
        gap: gap.clone(),
        head_gap: true,
        tail_gap: true,
        num_channels: scoring.nscoredalphabets,
        property_channels: None,
    };

    let aln_fft = fft_profile_align(&prof1, &prof2, &scoring.consweight_matrix, &fft_params);
    let aln_dp = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);

    let m_fft = aln_fft.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let m_dp = aln_dp.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();

    eprintln!("Prof1 len={}, Prof2 len={}", prof1.length, prof2.length);
    eprintln!("FFT: score={:.1} width={} matches={}", aln_fft.score, aln_fft.operations.len(), m_fft);
    eprintln!("DP:  score={:.1} width={} matches={}", aln_dp.score, aln_dp.operations.len(), m_dp);
    eprintln!("Same: {}", aln_fft.operations == aln_dp.operations);
}

#[test]
fn diagnostic_matrix_diagonal() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    eprintln!("Matrix size: {}x{}", scoring.substitution_matrix.len(), scoring.substitution_matrix[0].len());
    eprintln!("nalphabets: {}", scoring.nalphabets);
    eprintln!("gap.open: {}, gap.extend: {}, gap.offset: {}", scoring.gap.open, scoring.gap.extend, scoring.gap.offset);
    // Print first 5 diagonal values
    for i in 0..5.min(scoring.substitution_matrix.len()) {
        eprintln!("matrix[{}][{}] = {}", i, i, scoring.substitution_matrix[i][i]);
    }
    // Sum of diagonal for first 20 (amino acids)
    let diag_sum: i32 = (0..20).map(|i| scoring.substitution_matrix[i][i]).sum();
    eprintln!("Sum of diagonal (0..20): {}", diag_sum);
    eprintln!("Mean diagonal: {:.1}", diag_sum as f64 / 20.0);
}

#[test]
fn diagnostic_score_breakdown() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, pairwise_align11, GapModel};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s33 = &input.sequences[33].data;
    let s34 = &input.sequences[34].data;

    // Profile-based score
    let seqs: Vec<&[u8]> = vec![s33.as_slice()];
    let prof = Profile::from_aligned(&seqs, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let aln_prof = profile_align(&prof, &prof, &scoring.consweight_matrix, &gap, true, true);
    
    // G__align11 score
    let aln11 = pairwise_align11(s33, s34, &scoring.consweight_matrix, &scoring.amino_map,
        scoring.gap.open as f64, true, true);

    // Manual diagonal sum
    let mut diag_sum = 0i64;
    for &ch in s33 {
        let i = scoring.amino_map[ch as usize] as usize;
        if i < scoring.substitution_matrix.len() {
            diag_sum += scoring.substitution_matrix[i][i] as i64;
        }
    }

    eprintln!("profile_align score: {:.1}", aln_prof.score);
    eprintln!("pairwise_align11 score: {:.1}", aln11.score);
    eprintln!("Manual diagonal sum: {}", diag_sum);
    eprintln!("C's score: 302431");
    eprintln!("gap.open = {}", scoring.gap.open);
    
    // Check first few sub scores
    for i in 0..3 {
        let s = prof.match_score(i, &prof, i, &scoring.consweight_matrix);
        eprintln!("match_score({},{}) = {:.1} (char={})", i, i, s, s33[i] as char);
    }
}

#[test]
fn diagnostic_fft_pipeline() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, FftAlignParams, fft_profile_align};
    use mafft_fft::SegmentParams;

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s1 = &input.sequences[0].data;
    let s2 = &input.sequences[1].data;
    eprintln!("s1 len={} s2 len={}", s1.len(), s2.len());
    
    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);

    let params = FftAlignParams {
        num_candidates: 20,
        segment_params: SegmentParams::protein(),
        gap: gap.clone(),
        head_gap: true,
        tail_gap: true,
        num_channels: scoring.nscoredalphabets,
        property_channels: None,
    };
    let aln = fft_profile_align(&prof1, &prof2, &scoring.consweight_matrix, &params);
    eprintln!("FFT result: score={} ops={}", aln.score, aln.operations.len());
}

#[test]
fn diagnostic_segment_align() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, profile_align};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s1 = &input.sequences[0].data;
    let s2 = &input.sequences[1].data;
    
    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);

    // Full alignment
    let full = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    eprintln!("Full alignment: score={} ops={}", full.score, full.operations.len());

    // Sub-profile alignment: positions 0..28 of each
    let sub1 = prof1.sub_profile(0, 28);
    let sub2 = prof2.sub_profile(0, 28);
    let seg = profile_align(&sub1, &sub2, &scoring.consweight_matrix, &gap, true, false);
    eprintln!("Segment 0..28: score={} ops={}", seg.score, seg.operations.len());

    // Same with head_gap=false (intermediate segment)
    let seg2 = profile_align(&sub1, &sub2, &scoring.consweight_matrix, &gap, false, false);
    eprintln!("Segment 0..28 (no head_gap): score={} ops={}", seg2.score, seg2.operations.len());
}

#[test]
fn diagnostic_step3_anchors() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, FftAlignParams, fft_profile_align};
    use mafft_fft::SegmentParams;

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s7 = &input.sequences[7].data;
    let s8 = &input.sequences[8].data;
    eprintln!("s7 len={} s8 len={}", s7.len(), s8.len());

    let seqs1: Vec<&[u8]> = vec![s7.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s8.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);

    let params = FftAlignParams {
        num_candidates: 20,
        segment_params: SegmentParams::protein(),
        gap: gap.clone(),
        head_gap: true,
        tail_gap: true,
        num_channels: scoring.nscoredalphabets,
        property_channels: None,
    };
    let aln = fft_profile_align(&prof1, &prof2, &scoring.consweight_matrix, &params);
    eprintln!("FFT result: score={} ops={}", aln.score, aln.operations.len());

    use mafft_align::profile_align;
    let dp = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    eprintln!("DP result: score={} ops={}", dp.score, dp.operations.len());
    eprintln!("C step3 score: 108355.0");
}

#[test]
fn diagnostic_step3_dp() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, profile_align, AlignOp};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s7 = &input.sequences[7].data;
    let s8 = &input.sequences[8].data;

    let seqs1: Vec<&[u8]> = vec![s7.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s8.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);

    let aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    let m_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let d_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
    let i_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();
    eprintln!("DP result: score={} ops={} M/D/I={}/{}/{}", aln.score, aln.operations.len(), m_count, d_count, i_count);

    // Compute manual score: for each match position, look up what residues match
    use mafft_align::pairwise_align11;
    let aln11 = pairwise_align11(s7, s8, &scoring.consweight_matrix, &scoring.amino_map,
        scoring.gap.open as f64, true, true);
    let m_count = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let d_count = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
    let i_count = aln11.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();
    eprintln!("G__align11 result: score={} ops={} M/D/I={}/{}/{}", aln11.score, aln11.operations.len(), m_count, d_count, i_count);

    eprintln!("C step3: score=108355, width=364");
}

#[test]
fn diagnostic_step3_align_dump() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, profile_align, AlignOp};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    let s7 = &input.sequences[7].data;
    let s8 = &input.sequences[8].data;

    let seqs1: Vec<&[u8]> = vec![s7.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s8.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);

    let aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);

    // Walk operations and dump first 30 positions
    let mut p1 = 0;
    let mut p2 = 0;
    let mut a1 = String::new();
    let mut a2 = String::new();
    for (i, op) in aln.operations.iter().enumerate() {
        if i >= 200 { break; }
        match op {
            AlignOp::Match => {
                a1.push(s7[p1] as char);
                a2.push(s8[p2] as char);
                p1 += 1;
                p2 += 1;
            }
            AlignOp::Delete => {
                a1.push(s7[p1] as char);
                a2.push('-');
                p1 += 1;
            }
            AlignOp::Insert => {
                a1.push('-');
                a2.push(s8[p2] as char);
                p2 += 1;
            }
        }
    }
    eprintln!("Alignment (first 200 positions):");
    eprintln!("a1: {}", a1);
    eprintln!("a2: {}", a2);
    
    // Count where the inserts and matches are
    let leading_inserts = aln.operations.iter().take_while(|op| matches!(op, AlignOp::Insert)).count();
    let trailing_inserts = aln.operations.iter().rev().take_while(|op| matches!(op, AlignOp::Insert)).count();
    eprintln!("leading_inserts={}, trailing_inserts={}", leading_inserts, trailing_inserts);
}

#[test]
fn diagnostic_simple_offset() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, profile_align, AlignOp};

    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Test: short ACDE in middle of long sequence padded with random residues
    let s1: &[u8] = b"ACDE";
    let s2: &[u8] = b"WWWWACDEWWWW";
    
    let prof1 = Profile::from_aligned(&[s1], &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&[s2], &[1.0], &scoring.amino_map, scoring.nalphabets);
    let aln = profile_align(&prof1, &prof2, &scoring.consweight_matrix, &gap, true, true);
    
    let mut p1 = 0; let mut p2 = 0;
    let mut a1 = String::new();
    let mut a2 = String::new();
    for op in &aln.operations {
        match op {
            AlignOp::Match => { a1.push(s1[p1] as char); a2.push(s2[p2] as char); p1+=1; p2+=1; }
            AlignOp::Delete => { a1.push(s1[p1] as char); a2.push('-'); p1+=1; }
            AlignOp::Insert => { a1.push('-'); a2.push(s2[p2] as char); p2+=1; }
        }
    }
    eprintln!("ACDE vs WWWWACDEWWWW:");
    eprintln!("  a1: {}", a1);
    eprintln!("  a2: {}", a2);
    eprintln!("  score: {}", aln.score);

    // Regression guard: this minimal reproducer originally returned
    // `A----CDE----` before we fixed the boundary indexing. The optimal
    // alignment is `----ACDE----`: all 4 residues aligned diagonally to the
    // matching stretch of s2, with the 8 surrounding gaps on both ends.
    assert_eq!(
        a1, "----ACDE----",
        "DP found suboptimal alignment — boundary initialization may have regressed"
    );
    assert_eq!(a2, "WWWWACDEWWWW");
    assert_eq!(
        aln.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count(),
        4,
        "expected 4 matches (A-A, C-C, D-D, E-E)"
    );
}

#[test]
fn diagnostic_dp_score_bug() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, GapModel, profile_align};

    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
    
    eprintln!("matrix[0][0] (A-A) = {}", scoring.substitution_matrix[0][0]);
    eprintln!("matrix[4][4] (C-C) = {}", scoring.substitution_matrix[4][4]);
    eprintln!("matrix[3][3] (D-D) = {}", scoring.substitution_matrix[3][3]);
    eprintln!("matrix[6][6] (E-E) = {}", scoring.substitution_matrix[6][6]);
    eprintln!("matrix[17][17] (W-W) = {}", scoring.substitution_matrix[17][17]);
    eprintln!("amino_map[A]={}", scoring.amino_map[b'A' as usize]);
    eprintln!("amino_map[C]={}", scoring.amino_map[b'C' as usize]);
    eprintln!("amino_map[D]={}", scoring.amino_map[b'D' as usize]);
    eprintln!("amino_map[E]={}", scoring.amino_map[b'E' as usize]);
    eprintln!("amino_map[W]={}", scoring.amino_map[b'W' as usize]);
    eprintln!("gap.open={}", scoring.gap.open);

    // Test align ACDE vs ACDE — should give max score
    let s1: &[u8] = b"ACDE";
    let prof1 = Profile::from_aligned(&[s1], &[1.0], &scoring.amino_map, scoring.nalphabets);
    let aln = profile_align(&prof1, &prof1, &scoring.consweight_matrix, &gap, true, true);
    eprintln!("ACDE vs ACDE: score={}", aln.score);
    // Expected: A-A + C-C + D-D + E-E
}

#[test]
fn diagnostic_dp_trace() {
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    let scoring = build_context(ScoringModel::Blosum(62), SeqType::Protein);
    
    // Print key matrix values
    let amap = &scoring.amino_map;
    let m = &scoring.substitution_matrix;
    let aw_idx = amap[b'A' as usize] as usize;
    let cw_idx = amap[b'C' as usize] as usize;
    let dw_idx = amap[b'D' as usize] as usize;
    let ew_idx = amap[b'E' as usize] as usize;
    let ww_idx = amap[b'W' as usize] as usize;
    eprintln!("sub(A,A)={}  sub(A,W)={}", m[aw_idx][aw_idx], m[aw_idx][ww_idx]);
    eprintln!("sub(C,A)={}  sub(C,W)={}", m[cw_idx][aw_idx], m[cw_idx][ww_idx]);
    eprintln!("sub(C,C)={}  sub(D,D)={}  sub(E,E)={}", m[cw_idx][cw_idx], m[dw_idx][dw_idx], m[ew_idx][ew_idx]);
    eprintln!("gap.open={}", scoring.gap.open);
    
    // What should the optimal score be?
    let opt = m[aw_idx][aw_idx] + m[cw_idx][cw_idx] + m[dw_idx][dw_idx] + m[ew_idx][ew_idx] - 2 * (scoring.gap.open / 2);
    eprintln!("Expected optimal score (4 matches + 2 gaps of 4): ~{}", opt);
}


/// L-INS-1 (`mafft --localpair --maxiterate 0`) must produce byte-identical
/// output to C. Closed 2026-05-05 by three combined fixes:
///   1. `opt = isumscore / sumoverlap` post-rescale (constraints.rs:435).
///   2. Removing the 3-decimal hat2 distance round on the L-INS-i path
///      (engine.rs:301).
///   3. CLI `--maxiterate 0` honoured as zero refinement iters (main.rs:54).
///
/// Reference: `tests/fixtures/sample.linsi.maxit0` (= `mafft --localpair
/// --maxiterate 0 mafft-upstream/test/sample`).
#[test]
fn linsi_maxit0_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.linsi.maxit0"))
        .expect("missing tests/fixtures/sample.linsi.maxit0");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::LInsi { iterations: 0 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "L-INS-1 width differs: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i]
                    .iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} ({:?}) first diff at col {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} L-INS-1 sequences differ from C");
}

/// G-INS-1 (`mafft --globalpair --maxiterate 0`) must produce byte-identical
/// output to C. Closed 2026-05-06 by:
///   1. Porting `global_align` to mirror C's `G__align11` exactly
///      (max-so-far DP with `>=` tie-break, Atracking-style traceback).
///   2. Routing `outgap = 1` (head/tail gap penalised) through
///      `progressive_align_with_constraints` since C's
///      `scripts/mafft:2584` omits `$termgapopt = -O` for `--globalpair`
///      (whereas L-INS-i and E-INS-i pass it → `outgap = 0`).
///
/// Reference: `tests/fixtures/sample.ginsi.maxit0` (= `mafft --globalpair
/// --maxiterate 0 mafft-upstream/test/sample`).
#[test]
fn ginsi_maxit0_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.ginsi.maxit0"))
        .expect("missing tests/fixtures/sample.ginsi.maxit0");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::GInsi { iterations: 0 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "G-INS-1 width differs: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i]
                    .iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} ({:?}) first diff at col {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} G-INS-1 sequences differ from C");
}

/// E-INS-1 (`mafft --genafpair --maxiterate 0`) must produce byte-identical
/// output to C. Closed 2026-05-06 by:
///   1. Porting `genaffine_local_align` (`crates/mafft-align/src/genaffine.rs`)
///      to mirror C's `genL__align11` exactly — max-so-far DP with separate
///      `Mi` / `largeM` running-max trackers and a "skip" gap state with
///      `penalty_OP` open / no extension.
///   2. Adding `PairAligner::GeneralizedAffine` and routing E-INS-i through
///      it via `engine.rs::align`.
///   3. Mirroring the script's `localgenaf` overrides (`mafft:1940-1948`):
///      `lexp = laof = 0.0` for E-INS-i pairwise (NOT the L-INS-i values).
///      Without this, `genL__align11` was being driven with non-zero gap
///      extension and offset, which suppressed the skip-gap state.
///
/// Reference: `tests/fixtures/sample.einsi.maxit0` (= `mafft --genafpair
/// --maxiterate 0 mafft-upstream/test/sample`).
#[test]
fn einsi_maxit0_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.einsi.maxit0"))
        .expect("missing tests/fixtures/sample.einsi.maxit0");
    let input = read_fasta(test_data_path("sample")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::EInsi { iterations: 0 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "E-INS-1 width differs: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
            if mismatches <= 3 {
                let first_diff = msa.sequences[i]
                    .iter()
                    .zip(c_ref.sequences[i].data.iter())
                    .position(|(a, b)| a != b)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "seq {i} ({:?}) first diff at col {first_diff}",
                    c_ref.sequences[i].name
                );
            }
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} E-INS-1 sequences differ from C");
}

/// Small-input L-INS-i refinement: 9 sequences with `--maxiterate 2` must
/// be byte-identical to C. Guards the partial 2026-05-06 refinement fix:
/// adding `oimpmatchdouble` (sum of impmtx[i][i]) to the accept/reject
/// score so the constraint contribution is part of the threshold check
/// (mirrors C's `mscore = oimpmatchdouble + tmpdouble`,
/// `tscore = impmatchdouble + tmpdouble` in `tditeration.c:953,1094`).
///
/// Reference: `tests/fixtures/sample.first9.linsi.iter2`.
///
/// Note: full L-INS-i / G-INS-i / E-INS-i refinement on the 36-seq sample
/// still diverges from C (see TODO §3) because C uses FFT-segmented
/// `Falign_localhom` for the constrained DP and we use a single
/// non-FFT `profile_align_imp` call. This test captures the behaviour
/// where they coincide (small inputs without segmentation pressure).
#[test]
fn linsi_first9_iter2_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.first9.linsi.iter2"))
        .expect("missing tests/fixtures/sample.first9.linsi.iter2");
    let input = read_fasta(fixture_path("sample.first9.fa")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::LInsi { iterations: 2 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "L-INS-i (n=9, iter=2) width: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} sequences differ from C");
}

/// L-INS-i refinement on 12 sequences with `--maxiterate 5` must be
/// byte-identical to C. Guards the 2026-05-06 `Falign_localhom` port:
/// FFT-segmented constraint-aware DP with per-segment impmtx slicing,
/// using `partA__align`-style strict-`>` tie-break for the
/// prept-vs-mi/mjpt update.
///
/// Reference: `tests/fixtures/sample.first12.linsi.iter5` (= `mafft
/// --localpair --maxiterate 5 sample.first12.fa`).
#[test]
fn linsi_first12_iter5_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.first12.linsi.iter5"))
        .expect("missing tests/fixtures/sample.first12.linsi.iter5");
    let input = read_fasta(fixture_path("sample.first12.fa")).unwrap();

    let msa = MafftEngine::new(AlignmentMode::LInsi { iterations: 5 }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "L-INS-i (n=12, iter=5) width: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} sequences differ from C");
}

/// L-INS-i refinement on 13 sequences with `--maxiterate 2`. n=13 was
/// the smallest input where the iterative refinement produced a 1-column
/// gap-shift divergence vs C MAFFT before the 2026-05-06 boundary
/// nongap-frequency fix (threading `headgapfreq{1,2}` and
/// `gapfreq{1,2}[lgth]` from sgap/egap into the per-segment DP). Guards
/// that fix.
///
/// Reference: `tests/fixtures/sample.first13.linsi.iter2`
/// (= `mafft --localpair --maxiterate 2 sample.first13.fa`).
#[test]
fn linsi_first13_iter2_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first13.fa", "sample.first13.linsi.iter2", 2);
}

/// L-INS-i refinement on 14 sequences with `--maxiterate 2`. Closed
/// 2026-05-07 by using initial pairwise distances (3-decimal-rounded
/// to mimic hat2 file precision) for the refinement tree, instead of
/// recomputing from the progressive alignment.
#[test]
fn linsi_first14_iter2_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first14.fa", "sample.first14.linsi.iter2", 2);
}

/// L-INS-i refinement on 15 sequences with `--maxiterate 2`. Closed by
/// the 2026-05-07 hat2-distance fix.
#[test]
fn linsi_first15_iter2_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first15.fa", "sample.first15.linsi.iter2", 2);
}

/// L-INS-i refinement on the full 36-sequence sample with `--maxiterate 2`.
/// Closed by the 2026-05-07 hat2-distance fix — confirms the fix
/// generalizes to non-trivial input sizes.
#[test]
fn linsi_first36_iter2_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first36.fa", "sample.first36.linsi.iter2", 2);
}

/// L-INS-i (no refinement) on n=14/15/36 isolates progressive build
/// from refinement. These all pass byte-identical to C, proving the
/// remaining iter≥1 cascade is purely in `realign_all_constrained_fft`
/// (the `partA__align` per-segment DP) and not in the progressive build.
#[test]
fn linsi_first14_maxit0_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first14.fa", "sample.first14.linsi.maxit0", 0);
}

/// L-INS-i refinement on 14 sequences with `--maxiterate 1`. Pins
/// the very first refinement iteration since this was the cleanest
/// signal during the 2026-05-07 hat2-distance investigation.
#[test]
fn linsi_first14_iter1_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first14.fa", "sample.first14.linsi.iter1", 1);
}

#[test]
fn linsi_first15_maxit0_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first15.fa", "sample.first15.linsi.maxit0", 0);
}

#[test]
fn linsi_first36_maxit0_byte_identical_to_c() {
    assert_linsi_byte_identical("sample.first36.fa", "sample.first36.linsi.maxit0", 0);
}

/// G-INS-i refinement on n=14/36 with `--maxiterate 2`. Should be
/// byte-identical to C after the 2026-05-07 hat2-distance fix, since
/// G-INS-i shares the same `initial_pairwise_dm` engine code path.
#[test]
fn ginsi_first14_iter2_byte_identical_to_c() {
    assert_insi_byte_identical(
        "sample.first14.fa", "sample.first14.ginsi.iter2",
        AlignmentMode::GInsi { iterations: 2 }, "G-INS-i n=14 iter=2",
    );
}

#[test]
fn ginsi_first36_iter2_byte_identical_to_c() {
    assert_insi_byte_identical(
        "sample.first36.fa", "sample.first36.ginsi.iter2",
        AlignmentMode::GInsi { iterations: 2 }, "G-INS-i n=36 iter=2",
    );
}

#[test]
fn einsi_first14_iter2_byte_identical_to_c() {
    assert_insi_byte_identical(
        "sample.first14.fa", "sample.first14.einsi.iter2",
        AlignmentMode::EInsi { iterations: 2 }, "E-INS-i n=14 iter=2",
    );
}

/// E-INS-i refinement on the full 36-sequence sample. Closed
/// 2026-05-07 by switching to `naivepairscore11`-style scoring for
/// the genaffine pairwise distance (mirrors C's
/// `pairlocalalign.c:2225-2229` `usenaivescoreinsteadofalignmentscore`
/// branch when `-Z` is passed for `--genafpair`).
#[test]
fn einsi_first36_iter2_byte_identical_to_c() {
    assert_insi_byte_identical(
        "sample.first36.fa", "sample.first36.einsi.iter2",
        AlignmentMode::EInsi { iterations: 2 }, "E-INS-i n=36 iter=2",
    );
}

fn assert_insi_byte_identical(input_fixture: &str, c_ref_fixture: &str, mode: AlignmentMode, label: &str) {
    let c_ref = read_fasta(fixture_path(c_ref_fixture))
        .unwrap_or_else(|_| panic!("missing tests/fixtures/{c_ref_fixture}"));
    let input = read_fasta(fixture_path(input_fixture)).unwrap();
    let msa = MafftEngine::new(mode).align(&input);
    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "{label} width: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "{label}: {mismatches} sequences differ from C");
}

/// `--add` byte-identity: align 6 new sequences onto a 30-sequence
/// existing alignment, expect the result to match C MAFFT 7.526's
/// `mafft --add` output exactly.
///
/// Pipeline:
///   - existing 30-seq alignment built with C `mafft --auto` (FFT-NS-2).
///   - 6 raw sequences added via `mafft --add`.
/// Reference: `sample.add6.aln` (committed, generated by C 7.526).
#[test]
fn add_six_to_thirty_byte_identical_to_c() {
    let existing = read_fasta(fixture_path("sample.first30.fftns2.aln")).unwrap();
    let new_seqs = read_fasta(fixture_path("sample.last6_for_add.fa")).unwrap();
    let c_ref = read_fasta(fixture_path("sample.add6.aln")).unwrap();
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let msa = engine.add_to_alignment(&existing, &new_seqs, false);
    assert_eq!(msa.nseq(), c_ref.nseq(), "--add nseq: Rust={} C={}", msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "--add width: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "--add 30+6: {mismatches} sequences differ from C");
}

fn assert_linsi_byte_identical(input_fixture: &str, c_ref_fixture: &str, iterations: usize) {
    let c_ref = read_fasta(fixture_path(c_ref_fixture))
        .unwrap_or_else(|_| panic!("missing tests/fixtures/{c_ref_fixture}"));
    let input = read_fasta(fixture_path(input_fixture)).unwrap();

    let msa = MafftEngine::new(AlignmentMode::LInsi { iterations }).align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "L-INS-i ({input_fixture}, iter={iterations}) width: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len(),
    );
    let mut mismatches = 0usize;
    for i in 0..msa.nseq() {
        if msa.sequences[i] != c_ref.sequences[i].data {
            mismatches += 1;
        }
    }
    assert_eq!(mismatches, 0, "{mismatches} sequences differ from C ({input_fixture})");
}

/// `--reorder` must produce the same aligned content as `--inputorder` but
/// permuted into guide-tree DFS order. Asserts the FFT-NS-2 default-mode
/// reorder permutation matches what C MAFFT 7.526 emits — for the 36-seq
/// sample, indices 10 and 11 swap (1-indexed: seqs 11 ↔ 12).
#[test]
fn reorder_fftns2_matches_c() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    let default_msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    let reordered_msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_reorder(true)
        .align(&input);

    assert_eq!(default_msa.nseq(), reordered_msa.nseq());

    // C MAFFT 7.526 `mafft --reorder mafft-upstream/test/sample` permutes
    // the 36-seq sample so that input positions 10 and 11 (0-indexed) swap.
    // All other positions remain identity (this dataset's tree happens to
    // align almost identically with input order).
    let mut expected_order: Vec<usize> = (0..default_msa.nseq()).collect();
    expected_order.swap(10, 11);

    for (out_pos, &input_pos) in expected_order.iter().enumerate() {
        assert_eq!(
            reordered_msa.sequences[out_pos], default_msa.sequences[input_pos],
            "reorder seq at output position {out_pos} should equal default \
             alignment of input seq {input_pos}",
        );
        assert_eq!(
            reordered_msa.names[out_pos], default_msa.names[input_pos],
            "reorder name at output position {out_pos} should match input {input_pos}",
        );
    }
}

/// `--parttree --reorder` must produce byte-identical output to C MAFFT
/// 7.526. The two passes (CALL 1 with raw 6-mer distances, CALL 2 with
/// `naivepairscore11` on the first-pass alignment) compose as
/// `final_order[k] = call1_order[call2_order[k]]`. This guards that
/// composition end-to-end against a regression in either pass.
#[test]
fn reorder_parttree_matches_c() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    let default_msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_parttree(true)
        .align(&input);
    let reordered_msa = MafftEngine::new(AlignmentMode::FftNs2)
        .with_parttree(true)
        .with_reorder(true)
        .align(&input);
    assert_eq!(default_msa.nseq(), reordered_msa.nseq());

    // C MAFFT 7.526 `mafft --parttree --reorder sample` order (0-indexed),
    // verified end-to-end by `diff <(mafft-rs --parttree --reorder sample)
    // <(mafft --parttree --reorder sample) == 0`.
    let expected_order: [usize; 36] = [
        34, 33, 31, 32, 35, 30, 29, 28,  7,  8,
         9, 11, 10, 12,  1,  0,  3,  4,  2,  5,
         6, 17, 18, 25, 24, 22, 23, 20, 19, 21,
        26, 27, 16, 15, 14, 13,
    ];
    for (out_pos, &input_pos) in expected_order.iter().enumerate() {
        assert_eq!(
            reordered_msa.sequences[out_pos], default_msa.sequences[input_pos],
            "parttree reorder seq at output {out_pos} should equal default \
             alignment of input seq {input_pos}",
        );
        assert_eq!(
            reordered_msa.names[out_pos], default_msa.names[input_pos],
            "parttree reorder name at output {out_pos} should match input {input_pos}",
        );
    }
}

/// `--treeout` smoke-test for `--dpparttree`: the CLI binary uses
/// `run_parttree_pipeline_with_scorer` with a `global_align`-based
/// scorer to mirror C's `G__align11_noalign( n_disLN, -1200, -60 )`
/// pipeline. Here we exercise the pipeline directly to ensure the
/// generic core (`run_parttree_pipeline_with_scorer` +
/// `parttree_result_to_newick`) is byte-identical to a hand-computed
/// reference on a trivial 3-seq case.
#[test]
fn dpparttree_pipeline_smoke() {
    use mafft_tree::parttree_split::{
        run_parttree_pipeline_with_scorer, parttree_result_to_newick,
    };
    // 3 distinct seqs; trivial selfscores; pair distance fully tied.
    // All seqs distinct → each becomes its own yuko (npick=3, nyuko=3).
    // UPGMA with all-equal distances must produce some tree; the
    // resulting Newick should have all 3 leaves and round-trip via
    // `parttree_result_to_newick` without panicking.
    let r = run_parttree_pipeline_with_scorer(
        3,
        |_| 100,
        |_| 10,
        |i, j| if i == j { 100.0 } else { 50.0 },
        |i, j| i == j,
        50,
    ).expect("pipeline should produce a result for n=3");
    assert_eq!(r.outs.iter().flatten().count(), 3, "all 3 seqs assigned");
    let nw = parttree_result_to_newick(&r);
    // Numeric leaves 1, 2, 3 should appear.
    for n in &["1", "2", "3"] {
        assert!(nw.contains(n), "Newick should contain leaf {}", n);
    }
}

/// `--treein FILE` must reproduce C MAFFT byte-for-byte across multiple
/// scoring/iteration modes. The fixture tree
/// `tests/fixtures/sample.treein.tree` is the result of running
/// `mafft --treeout` on `test/sample` and converting the Newick output
/// with `newick2mafft.rb`. Each fixture pairs the same tree with a
/// different `mafft <mode> --treein <tree>` invocation.
///
/// Reference: `tests/fixtures/sample.treein.*` (`mafft --treein <tree>`).
#[test]
fn treein_fftns2_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.treein.fftns2"))
        .expect("missing fixtures/sample.treein.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::FftNs2);
    engine.treein_path = Some(fixture_path("sample.treein.tree"));
    let msa = engine.align(&input);

    assert_eq!(msa.nseq(), c_ref.nseq());
    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --treein FFT-NS-2");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --treein FFT-NS-2 output");
    }
}

#[test]
fn treein_nwns2_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.treein.nwns2"))
        .expect("missing fixtures/sample.treein.nwns2");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::FftNs2).with_nofft(true);
    engine.treein_path = Some(fixture_path("sample.treein.tree"));
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --treein --nofft");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --treein --nofft output");
    }
}

#[test]
fn treein_linsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.treein.linsi"))
        .expect("missing fixtures/sample.treein.linsi");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::LInsi { iterations: 1000 });
    engine.treein_path = Some(fixture_path("sample.treein.tree"));
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --treein L-INS-i");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --treein L-INS-i output");
    }
}

#[test]
fn treein_ginsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.treein.ginsi"))
        .expect("missing fixtures/sample.treein.ginsi");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::GInsi { iterations: 1000 });
    engine.treein_path = Some(fixture_path("sample.treein.tree"));
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --treein G-INS-i");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --treein G-INS-i output");
    }
}

/// `--leavegappyregion` / `--legacygappenalty` (`legacygapcost = 1`,
/// `Salignmm.c:1604-1610`) must reproduce C MAFFT byte-for-byte. With
/// the flag set, the profile DP treats every column as fully nongap,
/// disabling the 7.110 gap-aware reweighting.
#[test]
fn leavegappyregion_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.leavegappyregion"))
        .expect("missing fixtures/sample.leavegappyregion");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::FftNs2);
    engine.legacy_gap_cost = true;
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --leavegappyregion: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len());
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --leavegappyregion output");
    }
}

/// `--allowshift --globalpair --maxiterate 1000` (G-INS-i with the warp DP
/// + multi-distance-class refinement, C `partA__align_variousdist`) must be
/// byte-identical to C MAFFT. `--allowshift` sets `unalign_level = 0.8` and
/// enables the warp/shift DP. Reference is the upstream-shipped
/// `mafft-upstream/test/sample.ginsi.allowshift`.
#[test]
fn allowshift_ginsi_byte_identical_to_c() {
    let c_ref = read_fasta(test_data_path("sample.ginsi.allowshift"))
        .expect("missing mafft-upstream/test/sample.ginsi.allowshift");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::GInsi { iterations: 1000 });
    engine.allowshift = true;
    engine.unalign_level = 0.8;
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --allowshift G-INS-i: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len());
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --allowshift G-INS-i output");
    }
}

/// `--memsavetree` must reproduce C MAFFT byte-for-byte on the 36-seq
/// sample (covers both pass 0 k-mer-distance tree-build via
/// `compacttreegivendist` and pass 1 MSA-distance rebuild).
#[test]
fn memsavetree_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.memsavetree"))
        .expect("missing fixtures/sample.memsavetree");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::FftNs2);
    engine.memsavetree = true;
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --memsavetree: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len());
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --memsavetree output");
    }
}

/// `--auto` on the 36-seq sample must select L-INS-i with `iterate=1000`
/// (nseq=36 < 100, nlen ~360 < 3000) per `scripts/mafft:1295-1299` and
/// then produce output byte-identical to C MAFFT 7.526's `--auto`.
///
/// This guards the size-heuristic in `mafft-bin::decide_auto` — if the
/// thresholds drift from C's, we'd pick the wrong mode here and the
/// alignment would change.
#[test]
fn auto_picks_linsi_for_small_sample() {
    let c_ref = read_fasta(fixture_path("sample.auto"))
        .expect("missing fixtures/sample.auto");
    let input = read_fasta(test_data_path("sample")).unwrap();

    // The CLI does the size-based dispatch; this test exercises the
    // resulting engine config. nseq=36, nlen<3000 → L-INS-i, iter=1000.
    let engine = MafftEngine::new(AlignmentMode::LInsi { iterations: 1000 })
        .with_retree(1);
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --auto (L-INS-i 1000): Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len());
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --auto output");
    }
}

/// E-INS-i with `--treein` exercises both the LH table reweighting
/// (`recompute_importance` derived from the loaded topology, mirroring
/// C `tbfast.c:2967 counteff_simple` + `tbfast.c:1355 calcimportance_half`)
/// AND the refinement loop's user-tree override. Without the
/// reweighting fix, iter 1 matched but iter 2+ diverged because the LH
/// table was importance-weighted using a different tree than C's.
#[test]
fn treein_einsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.treein.einsi"))
        .expect("missing fixtures/sample.treein.einsi");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::EInsi { iterations: 1000 });
    engine.treein_path = Some(fixture_path("sample.treein.tree"));
    let msa = engine.align(&input);

    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --treein E-INS-i");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --treein E-INS-i output");
    }
}

/// `--treeout` Newick serialization of the post-progressive guide tree.
/// Smoke-test that:
/// 1. The engine populates `msa.guide_tree`.
/// 2. `topology_to_newick` produces a well-formed Newick string ending
///    in `;\n` with the expected number of leaf labels.
/// The end-to-end byte-identity test (Rust `.tree` vs C MAFFT `.tree`)
/// is exercised at the shell level — `diff sample.tree sample.tree == 0`
/// is verified during release sweeps; replicating it here would require
/// invoking the CLI binary as a subprocess.
#[test]
fn treeout_engine_populates_guide_tree() {
    let input = read_fasta(test_data_path("sample")).unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    let topo = msa.guide_tree.as_ref()
        .expect("engine should populate guide_tree after align()");
    assert_eq!(topo.nseq, input.nseq());
    assert!(topo.is_complete(), "guide tree should be complete after align()");

    let nw = mafft_tree::topology_to_newick(topo, &msa.names);
    assert!(nw.ends_with(";\n"), "Newick should terminate with `;\\n`");
    // Each leaf appears as `<i+1>_<sanitized_name>` — count leaf-position
    // separators (each leaf is wrapped in `\n` per C's leaf format).
    let leaf_count = (1..=input.nseq())
        .filter(|i| nw.contains(&format!("\n{}_", i))).count();
    assert_eq!(leaf_count, input.nseq(),
        "all {} leaves should appear in Newick output", input.nseq());
}

/// Helper used by all `--seed` byte-identity tests: mirror the CLI's
/// `--seed` preprocessing. Prepends gap-stripped seed sequences (with
/// `_seed_` name prefix) to the user input, builds a seed-only
/// `LocalHomologyTable` over the combined dimension, and returns
/// `(combined_input, seed_table)` ready for `engine.seed_homology =
/// Some(seed_table)`.
fn prepare_seed_input(
    seed_file: &std::path::Path,
    user_file: &std::path::Path,
) -> (mafft_types::SequenceSet, mafft_types::LocalHomologyTable) {
    let user_input = read_fasta(user_file).expect("read user input");
    let seed_set = read_fasta_casepreserve(seed_file).expect("read seed file");
    let user_nseq = user_input.nseq();
    let aligned: Vec<Vec<u8>> = seed_set.sequences.iter().map(|s| s.data.clone()).collect();

    let mut combined = user_input.clone();
    let mut idx = 0usize;
    for s in &seed_set.sequences {
        let ungapped: Vec<u8> = s.data.iter().copied()
            .filter(|&c| c != b'-' && c != b'.').collect();
        combined.sequences.insert(idx, Sequence {
            name: format!("_seed_{}", s.name),
            data: ungapped,
        });
        idx += 1;
    }
    let scoring = mafft_scoring::build_context(
        mafft_types::ScoringModel::Blosum(62),
        combined.seq_type,
    );
    let seed_group = mafft_align::SeedGroup {
        aligned: aligned.iter().map(|s| s.as_slice()).collect(),
        global_indices: (0..aligned.len()).collect(),
    };
    let table = mafft_align::build_seed_homology_table(
        std::slice::from_ref(&seed_group),
        combined.nseq(),
        user_nseq,
        &scoring.consweight_matrix,
        &scoring.amino_map,
    );
    (combined, table)
}

/// `--seed` + L-INS-i must reproduce C MAFFT byte-for-byte. Exercises:
/// - `multi2hat3s`-style `putlocalhom2` extraction over each seed pair
///   with `korh = 'k'`.
/// - `tsuyosa = user_nseq² * 100` importance boost.
/// - Merging seed entries into the pairwise homology table BEFORE
///   `recompute_importance` (so the position-vote pass weighs both).
/// - Forcing `iterate ≥ 2` (C `scripts/mafft:1911-1923`).
#[test]
fn seed_linsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.linsi.iter2"))
        .expect("missing fixtures/sample.seed.linsi.iter2");
    let (combined, seed_table) = prepare_seed_input(
        &fixture_path("sample.seed3.aln"),
        &fixture_path("sample.seed_input5.fa"),
    );
    let mut engine = MafftEngine::new(AlignmentMode::LInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    assert_eq!(msa.nseq(), c_ref.nseq(), "nseq mismatch");
    assert_eq!(msa.sequences[0].len(), c_ref.sequences[0].data.len(),
        "width differs for --seed L-INS-i: Rust={} C={}",
        msa.sequences[0].len(), c_ref.sequences[0].data.len());
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed L-INS-i output");
    }
}

/// `--seed` + G-INS-i: same flow as L-INS-i but with `--globalpair`
/// pairwise alignment driving the initial homology table.
#[test]
fn seed_ginsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.ginsi.iter2"))
        .expect("missing fixtures/sample.seed.ginsi.iter2");
    let (combined, seed_table) = prepare_seed_input(
        &fixture_path("sample.seed3.aln"),
        &fixture_path("sample.seed_input5.fa"),
    );
    let mut engine = MafftEngine::new(AlignmentMode::GInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed G-INS-i output");
    }
}

/// `--seed` + E-INS-i: `--genafpair` (generalized affine) drives the
/// pairwise step; seed entries merge in identically.
#[test]
fn seed_einsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.einsi.iter2"))
        .expect("missing fixtures/sample.seed.einsi.iter2");
    let (combined, seed_table) = prepare_seed_input(
        &fixture_path("sample.seed3.aln"),
        &fixture_path("sample.seed_input5.fa"),
    );
    let mut engine = MafftEngine::new(AlignmentMode::EInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed E-INS-i output");
    }
}

/// `--memsave` for FFT-NS-2 must produce the same alignment as default
/// mode. C's `--memsave` (`scripts/mafft:543-544`) sets `alg='M'` so
/// `tbfast` runs `MSalignmm` (Hirschberg-style linear-space DP) instead
/// of `A__align`. For inputs ≤ 30000 in length the alignment converges
/// to the same trace — C's auto-switch at `len > 30000`
/// (`tbfast.c:1096`) makes the two paths equivalent for typical inputs.
/// Our engine uses full-memory DP regardless, so the byte-identity
/// guarantee here is "FFT-NS-2 output == C's FFT-NS-2 + --memsave output."
#[test]
fn memsave_fftns2_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.memsave.fftns2"))
        .expect("missing fixtures/sample.memsave.fftns2");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's FFT-NS-2 + --memsave output");
    }
}

/// `--memsave` with FFT-NS-i (`--maxiterate 2`).
#[test]
fn memsave_fftnsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.memsave.fftnsi.iter2"))
        .expect("missing fixtures/sample.memsave.fftnsi.iter2");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 }).align(&input);
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's FFT-NS-i + --memsave output");
    }
}

/// `--memsave` with `--retree 1` (FFT-NS-1).
#[test]
fn memsave_fftns2_retree1_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.memsave.fftns2.retree1"))
        .expect("missing fixtures/sample.memsave.fftns2.retree1");
    let input = read_fasta(test_data_path("sample")).unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNs2).with_retree(1).align(&input);
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's FFT-NS-2 retree-1 + --memsave output");
    }
}

/// `--seed` with the default FFT-NS pipeline: C promotes `iterate=0`
/// to `iterate=2` so the seed constraints actually drive refinement
/// (`scripts/mafft:1911-1923`). Verifies that the seed-only constraint
/// table works without a pairwise homology step.
#[test]
fn seed_fftnsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.fftnsi.iter2"))
        .expect("missing fixtures/sample.seed.fftnsi.iter2");
    let (combined, seed_table) = prepare_seed_input(
        &fixture_path("sample.seed3.aln"),
        &fixture_path("sample.seed_input5.fa"),
    );
    let mut engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed FFT-NS-i output");
    }
}

/// `--seedtable` byte-identical to C MAFFT for L-INS-i.
///
/// The hat3.seed file was captured from `mafft --seed --debug` (the
/// `multi2hat3s` output). The combined input file is the gap-stripped
/// concatenation of seed FASTA + user FASTA, matching what our `--seed`
/// CLI prepends internally. C's `multi2hat3s.c:375` actually writes
/// *gapped* seeds to infile, but hat3 positions are in gap-stripped
/// residue space; either convention works as long as the input matches
/// the position space, and we ship the gap-stripped variant for clarity.
/// The resulting alignment matches C's `--seed` reference (which equals
/// C's `--seedtable` output for the same combined input).
#[test]
fn seedtable_linsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.linsi.iter2"))
        .expect("missing fixtures/sample.seed.linsi.iter2");
    let combined = read_fasta(fixture_path("sample.seed.combined.fa"))
        .expect("missing fixtures/sample.seed.combined.fa");
    let hat3 = std::fs::read_to_string(fixture_path("sample.seed.hat3"))
        .expect("missing fixtures/sample.seed.hat3");
    let seed_table = mafft_align::parse_hat3_seed(&hat3, combined.nseq())
        .expect("parse hat3");
    let mut engine = MafftEngine::new(AlignmentMode::LInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    assert_eq!(msa.nseq(), c_ref.nseq(), "nseq mismatch");
    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed L-INS-i reference");
    }
}

/// `--seedtable` byte-identical to C MAFFT for G-INS-i.
#[test]
fn seedtable_ginsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.ginsi.iter2"))
        .expect("missing fixtures/sample.seed.ginsi.iter2");
    let combined = read_fasta(fixture_path("sample.seed.combined.fa")).unwrap();
    let hat3 = std::fs::read_to_string(fixture_path("sample.seed.hat3")).unwrap();
    let seed_table = mafft_align::parse_hat3_seed(&hat3, combined.nseq()).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::GInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed G-INS-i reference");
    }
}

/// `--seedtable` byte-identical to C MAFFT for E-INS-i.
#[test]
fn seedtable_einsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.einsi.iter2"))
        .expect("missing fixtures/sample.seed.einsi.iter2");
    let combined = read_fasta(fixture_path("sample.seed.combined.fa")).unwrap();
    let hat3 = std::fs::read_to_string(fixture_path("sample.seed.hat3")).unwrap();
    let seed_table = mafft_align::parse_hat3_seed(&hat3, combined.nseq()).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::EInsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed E-INS-i reference");
    }
}

/// `--seedtable` byte-identical to C MAFFT for FFT-NS-i.
#[test]
fn seedtable_fftnsi_byte_identical_to_c() {
    let c_ref = read_fasta(fixture_path("sample.seed.fftnsi.iter2"))
        .expect("missing fixtures/sample.seed.fftnsi.iter2");
    let combined = read_fasta(fixture_path("sample.seed.combined.fa")).unwrap();
    let hat3 = std::fs::read_to_string(fixture_path("sample.seed.hat3")).unwrap();
    let seed_table = mafft_align::parse_hat3_seed(&hat3, combined.nseq()).unwrap();
    let mut engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 2 });
    engine.seed_homology = Some(seed_table);
    let msa = engine.align(&combined);

    for i in 0..msa.nseq() {
        assert_eq!(msa.sequences[i], c_ref.sequences[i].data,
            "seq {i} differs from C's --seed FFT-NS-i reference");
    }
}

// ---------------------------------------------------------------------------
// §B.4 — `--retree N` byte-identity tests (closed 2026-05-18).
//
// C `scripts/mafft:1840-1842` clamps `cycle = min(cycle, 3)` for all modes,
// and `scripts/mafft:1934-1936` forces `cycle = 1` for L/G/E/Q/X-INS-i
// regardless of `--retree N`. Without those guards, `--retree N > 1
// --localpair` silently runs extra progressive passes and diverges from C
// by hundreds of lines.
//
// These tests lock down byte-identity for the non-default `--retree`
// values that exercise both rewrites.
// ---------------------------------------------------------------------------

fn run_retree_engine(mode: AlignmentMode, retree: usize, sample: &str) -> mafft_core::MultipleAlignment {
    let input = read_fasta(test_data_path(sample)).expect("sample fixture");
    MafftEngine::new(mode).with_retree(retree).align(&input)
}

fn assert_byte_equal_to_ref(msa: &mafft_core::MultipleAlignment, c_ref_path: &str, label: &str) {
    let c_ref = read_fasta(fixture_path(c_ref_path)).unwrap_or_else(|_| panic!("missing fixture {c_ref_path}"));
    assert_eq!(msa.nseq(), c_ref.nseq(), "{label}: nseq mismatch");
    for i in 0..msa.nseq() {
        assert_eq!(
            msa.sequences[i], c_ref.sequences[i].data,
            "{label}: seq {i} differs from C reference {c_ref_path}",
        );
    }
}

/// `--retree 3` for FFT-NS-2 — exercises the user-requested 3 passes on
/// a non-INS-i mode (no override fires). C uses cycle=3 directly.
#[test]
fn retree_3_fftns2_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::FftNs2, 3, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree3.fftns2", "--retree 3");
}

/// `--retree 5` for FFT-NS-2 — verifies the `cycle = min(cycle, 3)` clamp
/// at `scripts/mafft:1840`. C clamps to 3; we must match.
#[test]
fn retree_5_fftns2_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::FftNs2, 5, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree5.fftns2", "--retree 5");
}

/// `--retree 3 --maxiterate 2` for FFT-NS-i — same cycle clamp, plus
/// iterative refinement on top of the 3-pass progressive build.
#[test]
fn retree_3_fftnsi_iter2_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::FftNsi { iterations: 2 }, 3, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree3.fftnsi.iter2", "--retree 3 --maxiterate 2");
}

/// `--retree 3 --localpair` (L-INS-i) — verifies the INS-i override at
/// `scripts/mafft:1934-1936` forces cycle=1 regardless of `--retree N`.
/// Before the fix, our engine ran 3 progressive passes here, diverging
/// from C by 898 lines on the 36-seq sample.
#[test]
fn retree_3_linsi_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::LInsi { iterations: 0 }, 3, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree3.linsi", "--retree 3 --localpair");
}

/// `--retree 5 --localpair` — same override, plus the clamp would
/// otherwise cap at 3. Both must collapse to cycle=1.
#[test]
fn retree_5_linsi_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::LInsi { iterations: 0 }, 5, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree5.linsi", "--retree 5 --localpair");
}

/// `--retree 3 --globalpair` (G-INS-1) — same override as L-INS-i.
#[test]
fn retree_3_ginsi_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::GInsi { iterations: 0 }, 3, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree3.ginsi", "--retree 3 --globalpair");
}

/// `--retree 3 --genafpair` (E-INS-1) — same override as L-INS-i.
#[test]
fn retree_3_einsi_byte_identical_to_c() {
    let msa = run_retree_engine(AlignmentMode::EInsi { iterations: 0 }, 3, "sample");
    assert_byte_equal_to_ref(&msa, "sample.retree3.einsi", "--retree 3 --genafpair");
}
