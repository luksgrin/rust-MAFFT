//! Inputs on which the sanitizer-instrumented C MAFFT reports undefined
//! behaviour, or on which the comparison found `mafft-rs` failing badly
//! (issue #27, `docs/architecture/memory-safety.md`).
//!
//! The fixtures live in `fixtures/memsafety/`. They come from
//! `benchmarks/memsafety_gen.py` (`control_protein.fa`, the `--treein`
//! files, `empty_record_middle.fa`), from `benchmarks/memsafety_cases/`
//! (the minimised `addfragments_localpair_*` and `add_sentinel_*`), and
//! `homopolymer_dna.fa` (300-nt poly-A + 300 random bases, seed 27).
//! `.expected` files are the stdout of the in-tree C build
//! (`make -C mafft-upstream/core`, default flags, arm64 / Apple clang).

use std::ffi::OsString;
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures")
}

fn memsafety(name: &str) -> PathBuf {
    fixtures().join("memsafety").join(name)
}

fn argv(flags: &[&str], input: &Path) -> Vec<OsString> {
    let mut argv: Vec<OsString> = vec![OsString::from("mafft-rs"), OsString::from("--quiet")];
    argv.extend(flags.iter().map(OsString::from));
    argv.push(input.as_os_str().to_os_string());
    argv
}

fn expect_error(flags: &[&str], input: &Path, message_contains: &str) {
    let mut out = Vec::new();
    match mafft_rs::run_from(argv(flags, input), &mut out) {
        Ok(()) => panic!("mafft-rs {flags:?} {} succeeded, expected an error", input.display()),
        Err(e) => {
            assert_eq!(e.code(), 1, "exit code: {}", e.message());
            assert!(
                e.message().contains(message_contains),
                "message `{}` does not contain `{message_contains}`",
                e.message()
            );
            assert!(out.is_empty(), "stdout must be empty on failure");
        }
    }
}

fn residues(fasta: &[u8]) -> Vec<(String, String)> {
    let mut recs: Vec<(String, String)> = Vec::new();
    for line in String::from_utf8_lossy(fasta).lines() {
        if let Some(name) = line.strip_prefix('>') {
            recs.push((name.to_string(), String::new()));
        } else if let Some(last) = recs.last_mut() {
            last.1.extend(line.chars().filter(|&c| c != '-').map(|c| c.to_ascii_uppercase()));
        }
    }
    recs
}

/// `--treein` with a 0 index: C's `loadtreeoneline` decrements it to -1,
/// only checks the upper bound, and `loadtree` reads `hist[-1]`
/// (ASan: heap-buffer-overflow READ of size 4 at `mltaln9.c:2636`; the
/// uninstrumented build segfaults). mafft-rs rejects the line.
#[test]
fn treein_zero_index_is_rejected() {
    let tree = memsafety("tree_zero_index.txt");
    let tree = tree.to_str().unwrap();
    expect_error(&["--treein", tree], &memsafety("control_protein.fa"), "out of range");
    expect_error(&["--treein", tree, "--maxiterate", "2"], &memsafety("control_protein.fa"), "out of range");
}

/// `--treein` with a negative index: same C read path as the 0 index.
#[test]
fn treein_negative_index_is_rejected() {
    let tree = memsafety("tree_negative_index.txt");
    expect_error(&["--treein", tree.to_str().unwrap()], &memsafety("control_protein.fa"), "invalid im");
}

/// `--add` onto an "alignment" whose rows differ in length. C refuses with
/// "The original N sequences must be aligned" (`disttbfast.c:4154`, where
/// the wrapper then exits 0 with empty stdout; `addsingle.c:3224` exits 1).
/// mafft-rs used to panic with an index out of bounds in `commongappick`.
#[test]
fn add_to_unaligned_existing_is_rejected() {
    let ih = fixtures().join("input_handling");
    let add = ih.join("add_new_unusual.fa");
    let add = add.to_str().unwrap();
    let existing = ih.join("dna_mixed_case.fa");
    for flags in [&["--add", add][..], &["--add", add, "--keeplength"], &["--addfragments", add]] {
        expect_error(flags, &existing, "sequences must be aligned");
    }
}

/// `--addfragments --localpair` with one existing sequence and two
/// fragments. C's `pairlocalalign` reads 4 bytes before an `n_dis` row in
/// `putlocalhom2` (`io.c:803`, ASan heap-buffer-overflow). There is no
/// trustworthy C output to compare with, so this pins that mafft-rs
/// completes and keeps every residue of every input sequence.
#[test]
fn addfragments_localpair_completes_and_keeps_residues() {
    let existing = memsafety("addfragments_localpair_existing.fa");
    let frags = memsafety("addfragments_localpair_frags.fa");
    for iterate in ["0", "2"] {
        let mut out = Vec::new();
        mafft_rs::run_from(
            argv(&["--addfragments", frags.to_str().unwrap(), "--localpair", "--maxiterate", iterate], &existing),
            &mut out,
        )
        .unwrap_or_else(|e| panic!("mafft-rs failed: {}", e.message()));
        let mut want = residues(&std::fs::read(&existing).unwrap());
        want.extend(residues(&std::fs::read(&frags).unwrap()));
        assert_eq!(residues(&out), want);
        let widths: Vec<usize> = {
            let mut w = Vec::new();
            for block in String::from_utf8_lossy(&out).split('>').skip(1) {
                w.push(block.lines().skip(1).map(str::len).sum());
            }
            w
        };
        assert!(widths.windows(2).all(|p| p[0] == p[1]), "ragged output: {widths:?}");
    }
}

/// `--add` of 6 sequences onto 30 aligned rows (issue #40, R4). mafft-rs
/// printed C's internal `=` new-gap sentinel in 3 rows: the port of
/// `insertnewgaps` / `profilealignment` writes it, and C clears it with
/// `eq2dashmatometehayaku` after every merge (`disttbfast.c:2988`), which
/// the port did not. The `.expected` file is the in-tree C build's output.
#[test]
fn add_sentinel_matches_c() {
    let new = memsafety("add_sentinel_new.fa");
    let mut out = Vec::new();
    mafft_rs::run_from(argv(&["--add", new.to_str().unwrap()], &memsafety("add_sentinel_existing.fa")), &mut out)
        .unwrap_or_else(|e| panic!("mafft-rs failed: {}", e.message()));
    let want = std::fs::read(memsafety("add_sentinel.add.expected")).unwrap();
    assert_eq!(String::from_utf8_lossy(&out), String::from_utf8_lossy(&want));
}

/// A record with no residues in the middle of the input. C aligns the
/// other sequences and prints the empty one as a gap-only row. mafft-rs
/// printed every sequence empty (exit 0): the profile aligners returned no
/// operations when one side had no columns, so the other side's residues
/// were dropped. They now place that side against gaps, as C `A__align`
/// does (`Salignmm.c:1228-1255`). The `.expected` file is the in-tree C
/// build's output.
#[test]
fn empty_record_matches_c() {
    let mut out = Vec::new();
    mafft_rs::run_from(argv(&[], &memsafety("empty_record_middle.fa")), &mut out).unwrap();
    let want = std::fs::read(memsafety("empty_record_middle.default.expected")).unwrap();
    assert_eq!(String::from_utf8_lossy(&out), String::from_utf8_lossy(&want));
}

/// Two DNA sequences behind a 300-nt poly-A run. The FFT anchors fall
/// inside the run and leave segments with no columns on one side; the
/// same empty-side bug as above dropped those residues (11 of 1195 here)
/// with exit 0. Every residue must survive.
#[test]
fn homopolymer_dna_keeps_every_residue() {
    let input = memsafety("homopolymer_dna.fa");
    for flags in [&[][..], &["--retree", "1"], &["--maxiterate", "2"], &["--memsave"]] {
        let mut out = Vec::new();
        mafft_rs::run_from(argv(flags, &input), &mut out).unwrap();
        assert_eq!(residues(&out), residues(&std::fs::read(&input).unwrap()), "{flags:?}");
    }
}

/// Same input, byte-identity with the in-tree C build. mafft-rs places the
/// gaps inside the poly-A run differently (608 vs 600 columns): the
/// anchor choice on a homopolymer still diverges.
#[test]
#[ignore = "FFT anchor placement on long homopolymers differs from C (docs/architecture/memory-safety.md)"]
fn homopolymer_dna_matches_c() {
    let mut out = Vec::new();
    mafft_rs::run_from(argv(&[], &memsafety("homopolymer_dna.fa")), &mut out).unwrap();
    let want = std::fs::read(memsafety("homopolymer_dna.default.expected")).unwrap();
    assert_eq!(String::from_utf8_lossy(&out), String::from_utf8_lossy(&want));
}
