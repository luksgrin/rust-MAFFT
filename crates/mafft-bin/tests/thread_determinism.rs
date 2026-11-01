//! `--thread N` must not change a single output byte.
//!
//! The pairwise all-vs-all phase, the progressive merges of independent
//! subtrees and FFT segments, the FFT-NS-i refinement segments and the
//! speculative evaluation of refinement branches all run on the per-call
//! `--thread` pool, but only over independent units whose results land in
//! fixed slots and are consumed in a fixed order. Each test aligns one
//! input in-process through `run_from_seqs` at `--thread 1, 2, 4, 8` and
//! requires identical rows. (`--thread 1` is itself pinned to C by the
//! fixture suites; this file pins every other thread count to it.)

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use mafft_rs::{run_from_seqs, SeqName, SequenceSet, SilentProgress};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures")
}

fn upstream_sample() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../mafft-upstream/test/sample")
        .canonicalize()
        .expect("mafft-upstream submodule checked out")
}

fn align(flags: &[&str], threads: usize, set: &SequenceSet) -> Vec<(SeqName, Vec<u8>)> {
    let mut argv = vec![OsString::from("mafft-rs"), OsString::from("--quiet")];
    argv.extend(flags.iter().map(OsString::from));
    argv.push("--thread".into());
    argv.push(threads.to_string().into());
    let msa = run_from_seqs(argv, set, &SilentProgress)
        .unwrap_or_else(|e| panic!("run_from_seqs {flags:?} --thread {threads}: {}", e.message()));
    msa.names.into_iter().zip(msa.sequences).collect()
}

fn assert_thread_invariant(input: &Path, flags: &[&str]) {
    let set = mafft_io::read_fasta(input).expect("read fixture");
    let reference = align(flags, 1, &set);
    assert!(!reference.is_empty());
    for threads in [2, 4, 8] {
        let rows = align(flags, threads, &set);
        assert_eq!(rows.len(), reference.len());
        for ((name, row), (ref_name, ref_row)) in rows.iter().zip(&reference) {
            assert_eq!(name, ref_name);
            assert!(
                row == ref_row,
                "{} {flags:?}: --thread {threads} row {name} differs from --thread 1",
                input.display(),
            );
        }
    }
}

#[test]
fn linsi_sample_is_thread_invariant() {
    assert_thread_invariant(&upstream_sample(), &["--localpair", "--maxiterate", "1000"]);
}

#[test]
fn ginsi_sample_is_thread_invariant() {
    assert_thread_invariant(&upstream_sample(), &["--globalpair", "--maxiterate", "1000"]);
}

#[test]
fn einsi_sample_is_thread_invariant() {
    assert_thread_invariant(&upstream_sample(), &["--genafpair", "--maxiterate", "1000"]);
}

#[test]
fn fftnsi_sample_is_thread_invariant() {
    assert_thread_invariant(&upstream_sample(), &["--maxiterate", "1000"]);
}

#[test]
fn fftns2_sample_is_thread_invariant() {
    assert_thread_invariant(&upstream_sample(), &[]);
}

/// Segmented FFT-NS-i refinement (85 segments) and a DNA progressive pass.
#[test]
fn fftnsi_dna_cds_is_thread_invariant() {
    assert_thread_invariant(&fixtures().join("mtb_cds_120x1400.fa"), &["--maxiterate", "2"]);
}

#[test]
fn linsi_dna_is_thread_invariant() {
    assert_thread_invariant(
        &fixtures().join("mtb_cds_first8.fa"),
        &["--localpair", "--maxiterate", "1000"],
    );
}
