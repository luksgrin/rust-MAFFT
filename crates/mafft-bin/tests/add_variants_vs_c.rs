//! Issue #50: `--add` variants against the in-tree C build.
//!
//! Each test compares `mafft-rs` stdout with the C reference in
//! `fixtures/add_variants/` on two inputs: the canonical 30 + 6 sample
//! and the issue #40 reproducer (30 aligned rows + 6 new). Tests that
//! do not match C yet are `#[ignore]`d so the gap stays visible.
//!
//! `--addfragments` uses the default 6-mer `addsingle` path. The
//! `--localpair` combination is not pinned: C reads off the end of a
//! buffer on some of those inputs (`docs/architecture/memory-safety.md`).

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

fn repo_fixture(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
}

fn expected(name: &str) -> PathBuf {
    repo_fixture("tests/fixtures/add_variants").join(name)
}

fn scratch(tag: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    let n = N.fetch_add(1, Ordering::Relaxed);
    let mut dir = std::env::temp_dir();
    dir.push(format!("mafft-rs-add50-{}-{}-{tag}", std::process::id(), n));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run `mafft-rs --quiet <flag> <new> [extra...] <existing>` and return stdout.
fn run_add(existing: &Path, new: &Path, flag: &str, extra: &[&str]) -> Vec<u8> {
    let dir = scratch(flag.trim_start_matches('-'));
    let ex = dir.join("existing.fa");
    let nw = dir.join("new.fa");
    std::fs::copy(existing, &ex).unwrap();
    std::fs::copy(new, &nw).unwrap();
    let mut argv: Vec<OsString> = vec!["mafft-rs".into(), "--quiet".into(), flag.into(), nw.into()];
    argv.extend(extra.iter().map(OsString::from));
    argv.push(ex.into());
    let mut out = Vec::new();
    mafft_rs::run_from(argv, &mut out).unwrap_or_else(|e| panic!("{flag} {extra:?}: {e}"));
    let _ = std::fs::remove_dir_all(&dir);
    out
}

fn expect_both(flag: &str, extra: &[&str], sample_name: &str, sentinel_name: &str) {
    let sample_ex = repo_fixture("../mafft-core/tests/fixtures/sample.first30.fftns2.aln");
    let sample_nw = repo_fixture("../mafft-core/tests/fixtures/sample.last6_for_add.fa");
    let sent_ex = repo_fixture("tests/fixtures/memsafety/add_sentinel_existing.fa");
    let sent_nw = repo_fixture("tests/fixtures/memsafety/add_sentinel_new.fa");
    for (label, existing, new, want_name) in [
        ("30+6", sample_ex, sample_nw, sample_name),
        ("#40", sent_ex, sent_nw, sentinel_name),
    ] {
        let got = run_add(&existing, &new, flag, extra);
        let want = std::fs::read(expected(want_name)).unwrap_or_else(|e| panic!("{want_name}: {e}"));
        assert_eq!(
            got, want,
            "{flag} {extra:?} on {label} differs from C ({want_name}): rust {} bytes, C {} bytes",
            got.len(),
            want.len(),
        );
    }
}

/// C `dvtditr -K nadd -B` after the progressive add (`scripts/mafft:2076`,
/// `:2691-2693`). No `--thread`, so `-C 0` (single-threaded rules).
#[test]
#[ignore = "differs from in-tree C --add --maxiterate 2 (#50)"]
fn add_maxiterate_2_matches_c() {
    expect_both(
        "--add",
        &["--maxiterate", "2"],
        "sample.add6.maxiterate2.aln",
        "sentinel.maxiterate2.aln",
    );
}

/// C writes `order` from `topolorderz` of disttbfast's second guide tree
/// (`-E 2`, `disttbfast.c:4081-4102`) and `f2cl -r` applies it. On these
/// inputs that order is input order with sequences 11 and 12 swapped.
#[test]
fn add_reorder_matches_c() {
    expect_both("--add", &["--reorder"], "sample.add6.reorder.aln", "sentinel.reorder.aln");
}

/// C `addsingle` (`scripts/mafft:2671`), not the progressive `--add` path.
#[test]
#[ignore = "--addfragments still follows the progressive --add path; C runs addsingle (addsingle.c, scripts/mafft:2671). #50"]
fn addfragments_matches_c() {
    expect_both(
        "--addfragments",
        &[],
        "sample.add6.addfragments.aln",
        "sentinel.addfragments.aln",
    );
}

/// Same `addsingle` path plus `-Y` (`--keeplength`).
#[test]
#[ignore = "--addfragments --keeplength still follows the progressive --add path; C runs addsingle -Y (addsingle.c, scripts/mafft:2671). #50"]
fn addfragments_keeplength_matches_c() {
    expect_both(
        "--addfragments",
        &["--keeplength"],
        "sample.add6.addfragments.keeplength.aln",
        "sentinel.addfragments.keeplength.aln",
    );
}

/// C's wrapper sets `-z -Y` for `--compactmapout` (`scripts/mafft:699-700`)
/// and `reconstructdeletemap_compact` rewrites added-sequence names to
/// `Nins:<first-run>,etc|<original>` (`addfunctions.c:2159`, `MODIFYNAME`).
#[test]
fn add_compactmapout_matches_c() {
    expect_both(
        "--add",
        &["--compactmapout"],
        "sample.add6.compactmapout.aln",
        "sentinel.compactmapout.aln",
    );
}
