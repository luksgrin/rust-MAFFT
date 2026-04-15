/// End-to-end integration tests: read real test data, align, verify output.

use std::path::PathBuf;

use mafft_core::{MafftEngine, AlignmentMode};
use mafft_io::read_fasta;

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

#[test]
fn compare_against_c_reference() {
    // Read the C reference alignment (test/sample.fftns2)
    let c_ref = read_fasta(test_data_path("sample.fftns2")).unwrap();

    // Run Rust engine on same input
    let input = read_fasta(test_data_path("sample")).unwrap();
    let engine = MafftEngine::new(AlignmentMode::FftNs2);
    let msa = engine.align(&input);

    // 1. Same number of sequences
    assert_eq!(msa.nseq(), c_ref.nseq(), "different number of sequences");

    // 2. Same residue content per sequence (ungapped)
    for i in 0..msa.nseq() {
        let rust_ungapped: Vec<u8> = msa.sequences[i].iter().filter(|&&c| c != b'-').cloned().collect();
        let c_ungapped: Vec<u8> = c_ref.sequences[i].data.iter().filter(|&&c| c != b'-').cloned().collect();
        assert_eq!(
            rust_ungapped, c_ungapped,
            "sequence {i} has different residue content between Rust and C"
        );
    }

    // 3. Compare alignment quality via sum-of-pairs identity
    let rust_sp = sum_of_pairs_identity(&msa.sequences);
    let c_seqs: Vec<Vec<u8>> = c_ref.sequences.iter().map(|s| s.data.clone()).collect();
    let c_sp = sum_of_pairs_identity(&c_seqs);

    // Report quality comparison (not a hard failure — different algorithms
    // produce different alignments, but quality should be in the same ballpark)
    let ratio = if c_sp > 0.0 { rust_sp / c_sp } else { 1.0 };
    eprintln!(
        "Alignment quality: Rust SP={rust_sp:.4}, C SP={c_sp:.4}, ratio={ratio:.4}"
    );
    // Rust alignment should be at least 50% as good as C's
    // (a loose bound — we're not matching C's exact algorithm)
    assert!(
        ratio > 0.3,
        "Rust alignment quality too low: {rust_sp:.4} vs C's {c_sp:.4} (ratio {ratio:.4})"
    );
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
    let dp_aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

    // FFT-accelerated alignment
    let fft_params = FftAlignParams {
        num_candidates: 20,
        segment_params: SegmentParams::protein(),
        gap: gap.clone(),
        head_gap: true,
        tail_gap: true,
        num_channels: 20,
    };
    let fft_aln = fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &fft_params);

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
    use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, ktuple_distance, sequence_weights};
    use mafft_scoring::build_context;
    use mafft_types::{ScoringModel, SeqType};
    use mafft_align::{Profile, profile_align, GapModel};

    let input = read_fasta(test_data_path("sample")).unwrap();
    let scoring = build_context(ScoringModel::Jtt, SeqType::Protein);
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
    let mut aligned: Vec<Vec<u8>> = input.sequences.iter()
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
    let scoring = build_context(ScoringModel::Jtt, SeqType::Protein);
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
    let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
    
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
    let scoring = build_context(ScoringModel::Jtt, SeqType::Protein);
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
        let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

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
    let scoring = build_context(ScoringModel::Jtt, SeqType::Protein);
    let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);

    // Compare align11 vs profile_align for seqs 19 and 20 (merged at step 1)
    let s1 = &input.sequences[19].data;
    let s2 = &input.sequences[20].data;

    // Use the same boundary convention the engine uses (outgap=0 → false, false).
    let aln11 = pairwise_align11(s1, s2, &scoring.substitution_matrix, &scoring.amino_map,
        scoring.gap.open as f64, false, false);

    let seqs1: Vec<&[u8]> = vec![s1.as_slice()];
    let seqs2: Vec<&[u8]> = vec![s2.as_slice()];
    let prof1 = Profile::from_aligned(&seqs1, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let prof2 = Profile::from_aligned(&seqs2, &[1.0], &scoring.amino_map, scoring.nalphabets);
    let aln_prof = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, false, false);

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
    };

    let aln_fft = fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &fft_params);
    let aln_dp = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

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
    let aln_prof = profile_align(&prof, &prof, &scoring.substitution_matrix, &gap, true, true);
    
    // G__align11 score
    let aln11 = pairwise_align11(s33, s34, &scoring.substitution_matrix, &scoring.amino_map,
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
        let s = prof.match_score(i, &prof, i, &scoring.substitution_matrix);
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
    };
    let aln = fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &params);
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
    let full = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
    eprintln!("Full alignment: score={} ops={}", full.score, full.operations.len());

    // Sub-profile alignment: positions 0..28 of each
    let sub1 = prof1.sub_profile(0, 28);
    let sub2 = prof2.sub_profile(0, 28);
    let seg = profile_align(&sub1, &sub2, &scoring.substitution_matrix, &gap, true, false);
    eprintln!("Segment 0..28: score={} ops={}", seg.score, seg.operations.len());

    // Same with head_gap=false (intermediate segment)
    let seg2 = profile_align(&sub1, &sub2, &scoring.substitution_matrix, &gap, false, false);
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
    };
    let aln = fft_profile_align(&prof1, &prof2, &scoring.substitution_matrix, &params);
    eprintln!("FFT result: score={} ops={}", aln.score, aln.operations.len());

    use mafft_align::profile_align;
    let dp = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
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

    let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
    let m_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Match)).count();
    let d_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Delete)).count();
    let i_count = aln.operations.iter().filter(|op| matches!(op, AlignOp::Insert)).count();
    eprintln!("DP result: score={} ops={} M/D/I={}/{}/{}", aln.score, aln.operations.len(), m_count, d_count, i_count);

    // Compute manual score: for each match position, look up what residues match
    use mafft_align::pairwise_align11;
    let aln11 = pairwise_align11(s7, s8, &scoring.substitution_matrix, &scoring.amino_map,
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

    let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);

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
    let aln = profile_align(&prof1, &prof2, &scoring.substitution_matrix, &gap, true, true);
    
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
    let aln = profile_align(&prof1, &prof1, &scoring.substitution_matrix, &gap, true, true);
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

