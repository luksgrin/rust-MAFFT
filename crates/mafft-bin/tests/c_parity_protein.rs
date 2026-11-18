//! Byte-for-byte parity with C MAFFT 7.526 on protein input.
//!
//! Every `.expected` file here is the verbatim stdout of the in-tree C build
//! (`mafft-upstream` at the pinned commit, `scripts/mafft` with the named
//! flags, `--thread 1`), so these tests need no C installation.
//!
//! * `protein_dp_branch_n200` — 200 synthetic sequences of 20–30 residues
//!   from two ancestors with subfamily-wide block indels. Clusters outgrow
//!   the representative length `nlen` that `disttbfast` keeps per sequence
//!   (`disttbfast.c:2850`, `ffttry = nlen[m1] > clus1 && nlen[m2] >
//!   clus2`), so the late merges take the DP branch, where `treebase`
//!   hands the cached profile of the previous DP merge (`cpmxhist`,
//!   `:2889`/`:2913`) to the next. The port tried the FFT aligner on
//!   every merge and rebuilt every profile; it diverged from C on this
//!   set under FFT-NS-2 (and on the PF00069 corpus from 500 sequences
//!   under `--memsave` and from 4000 under the defaults).

use std::ffi::OsString;
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures")
}

/// Run `mafft-rs <flags> <input>` in-process and return the FASTA bytes.
fn run(flags: &[&str], input: &Path) -> Vec<u8> {
    let mut argv: Vec<OsString> = vec![OsString::from("mafft-rs"), OsString::from("--quiet")];
    argv.extend(flags.iter().map(OsString::from));
    argv.push(input.as_os_str().to_os_string());
    let mut out = Vec::new();
    mafft_rs::run_from(argv, &mut out)
        .unwrap_or_else(|e| panic!("mafft-rs {flags:?} {}: {}", input.display(), e.message()));
    out
}

fn expect_identical(flags: &[&str], input: &Path, expected: &Path) {
    let got = run(flags, input);
    let want = std::fs::read(expected).expect("read expected");
    if got != want {
        panic!(
            "mafft-rs {:?} {} differs from C MAFFT 7.526 ({})\n--- C ---\n{}--- rust ---\n{}",
            flags,
            input.display(),
            expected.display(),
            String::from_utf8_lossy(&want),
            String::from_utf8_lossy(&got),
        );
    }
}

/// Bug: the progressive stage tried the FFT aligner on every merge and
/// never reused a DP merge's cached profile; `disttbfast` gates the FFT
/// aligner on `ffttry` and blends the cached child profile into the next
/// DP merge (`disttbfast.c:2850`, `:2925`, `:2889`).
#[test]
fn fftns2_dp_branch_merges_match_c() {
    let f = fixtures();
    expect_identical(
        &["--thread", "1"],
        &f.join("protein_dp_branch_n200.fa"),
        &f.join("protein_dp_branch_n200.fftns2.expected"),
    );
}

/// The DP-branch schedule and the profile cache are decided from the
/// guide tree alone, so the output must not depend on the thread count.
#[test]
fn fftns2_dp_branch_merges_thread_independent() {
    let f = fixtures();
    let input = f.join("protein_dp_branch_n200.fa");
    let t1 = run(&["--thread", "1"], &input);
    let t4 = run(&["--thread", "4"], &input);
    assert!(t1 == t4, "--thread 4 output differs from --thread 1");
}
