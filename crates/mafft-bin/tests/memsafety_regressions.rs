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

/// `--treein` branch lengths `nan`, `inf`, `-inf`, and `1e308`. C accepts
/// them (`loadtree` keeps `nan`/`+inf` and stores `-inf` as 0) and the
/// alignment is 132 columns. mafft-rs used to succeed with a different
/// alignment (288 columns).
#[test]
fn treein_nan_lengths_match_c() {
    let tree = memsafety("tree_nan_lengths.txt");
    let mut out = Vec::new();
    mafft_rs::run_from(
        argv(&["--thread", "1", "--treein", tree.to_str().unwrap()], &memsafety("control_protein.fa")),
        &mut out,
    )
    .unwrap_or_else(|e| panic!("mafft-rs failed: {}", e.message()));
    let want = std::fs::read(memsafety("tree_nan_lengths.fftns2.expected")).unwrap();
    assert_eq!(String::from_utf8_lossy(&out), String::from_utf8_lossy(&want));
}

/// `--treein` with a negative index: same C read path as the 0 index.
#[test]
fn treein_negative_index_is_rejected() {
    let tree = memsafety("tree_negative_index.txt");
    expect_error(&["--treein", tree.to_str().unwrap()], &memsafety("control_protein.fa"), "invalid im");
}

/// `--seed` whose rows differ in length. `multi2hat3s` exits 1
/// (`multi2hat3s.c:199-204`): `Not aligned,  =sd0 - =sd1`. The `=` is
/// what `readData_pointer` stores in place of `>`.
#[test]
fn seed_ragged_is_rejected() {
    let seed = memsafety("seed_ragged.fa");
    let mut out = Vec::new();
    match mafft_rs::run_from(argv(&["--seed", seed.to_str().unwrap()], &memsafety("control_protein.fa")), &mut out) {
        Ok(()) => panic!("ragged --seed was accepted"),
        Err(e) => {
            assert_eq!(e.code(), 1, "{}", e.message());
            assert!(
                e.message().contains("## ERROR  ###\nNot aligned,  =sd0 - =sd1\n## ERROR  ###"),
                "{}",
                e.message()
            );
            assert!(out.is_empty(), "stdout must be empty on failure");
        }
    }
}

/// `--seed` with one sequence. `multi2hat3s` exits 1
/// (`multi2hat3s.c:326-330`).
#[test]
fn seed_single_is_rejected() {
    let seed = memsafety("seed_single.fa");
    let mut out = Vec::new();
    match mafft_rs::run_from(argv(&["--seed", seed.to_str().unwrap()], &memsafety("control_protein.fa")), &mut out) {
        Ok(()) => panic!("single-sequence --seed was accepted"),
        Err(e) => {
            assert_eq!(e.code(), 1, "{}", e.message());
            assert!(
                e.message().contains("At least 2 sequences should be input!\nOnly 1 sequence found."),
                "{}",
                e.message()
            );
            assert!(out.is_empty(), "stdout must be empty on failure");
        }
    }
}

/// `--treein` that merges node 1 with node 2 on every line. C's `loadtree`
/// exits 0 and writes an alignment, with no diagnostic. mafft-rs rejects
/// the second merge. That rejection stays.
#[test]
fn treein_reused_node_is_rejected() {
    let tree = memsafety("tree_reused_node.txt");
    expect_error(
        &["--treein", tree.to_str().unwrap()],
        &memsafety("control_protein.fa"),
        "already merged",
    );
}

/// `--treein` whose first line is longer than C's 999-byte `fgets` buffer
/// (`mltaln9.c:1488`). C reports `Format error (2)` and exits 1. The line
/// is `1 2` plus 5000 zeros plus `.1 0.1`, which is the same merge as
/// `1 2 0.1 0.1` once the whole line is read. mafft-rs keeps it and aligns.
#[test]
fn treein_overlong_line_is_kept() {
    let tree = memsafety("tree_overlong_line.txt");
    let mut out = Vec::new();
    mafft_rs::run_from(
        argv(&["--treein", tree.to_str().unwrap()], &memsafety("control_protein.fa")),
        &mut out,
    )
    .unwrap_or_else(|e| panic!("overlong --treein was rejected: {}", e.message()));
    assert!(out.starts_with(b">"), "expected an alignment, got {:?}", String::from_utf8_lossy(&out));
    let widths: Vec<usize> = String::from_utf8_lossy(&out)
        .split('>')
        .skip(1)
        .map(|block| block.lines().skip(1).map(str::len).sum())
        .collect();
    assert!(widths.len() > 1, "too few rows: {widths:?}");
    assert!(widths.windows(2).all(|p| p[0] == p[1]), "ragged output: {widths:?}");
}

/// Empty file, a file with no header line, and a file of only NUL bytes.
/// The wrapper exits 0 and writes nothing (`disttbfast.c:3352-3365`).
#[test]
fn empty_headerless_and_nul_only_match_c() {
    let dir = std::env::temp_dir().join(format!("mafft-rs-empty-header-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let empty = dir.join("empty.fa");
    let nul = dir.join("nul_only.fa");
    std::fs::write(&empty, b"").unwrap();
    std::fs::write(&nul, vec![0u8; 64]).unwrap();

    for path in [empty, memsafety("no_header.fa"), nul] {
        let mut out = Vec::new();
        mafft_rs::run_from(argv(&[], &path), &mut out)
            .unwrap_or_else(|e| panic!("{}: {}", path.display(), e.message()));
        assert!(out.is_empty(), "{}", path.display());
    }
    let _ = std::fs::remove_dir_all(&dir);
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

fn protein_pair() -> mafft_rs::SequenceSet {
    mafft_rs::SequenceSet {
        sequences: vec![
            mafft_rs::Sequence { name: "a".into(), data: b"ACDEF".to_vec() },
            mafft_rs::Sequence { name: "b".into(), data: b"ACD-F".to_vec() },
        ],
        seq_type: mafft_rs::SeqType::Protein,
    }
}

/// `--treein` combined with `--add` skips the check that runs before
/// `align` (the combined sequence count is not known yet). A tree that
/// does not parse used to kill the process from inside the engine.
/// `run_from_seqs` must return the error and stay able to align.
#[test]
fn treein_with_add_returns_error() {
    let dir = std::env::temp_dir().join(format!("mafft-treein-add-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let add = dir.join("add.fa");
    let tree = dir.join("bad.tree");
    std::fs::write(&add, b">n\nACDE\n").unwrap();
    std::fs::write(&tree, b"0 2 0.1 0.1\n").unwrap();
    let existing = protein_pair();
    let err = mafft_rs::run_from_seqs(
        [
            "mafft-rs",
            "--quiet",
            "--treein",
            tree.to_str().unwrap(),
            "--add",
            add.to_str().unwrap(),
        ],
        &existing,
        &mafft_rs::SilentProgress,
    )
    .expect_err("failing --treein --add must return Err");
    assert_eq!(err.code(), 1, "{}", err.message());
    assert!(err.message().contains("--treein:"), "{}", err.message());
    let ok = mafft_rs::run_from_seqs(["mafft-rs", "--quiet"], &existing, &mafft_rs::SilentProgress)
        .expect("a later call must still align");
    assert_eq!(ok.nseq(), 2);
    let _ = std::fs::remove_dir_all(&dir);
}

/// Q-INS-i computes base-pair probabilities after the progressive merge
/// and used to `process::exit` when that step failed. Without
/// `mxscarnamod` on `PATH` the library must return that error. Two short
/// sequences, so the merge before the failure is small.
#[test]
fn qinsi_bpp_failure_returns_error() {
    // A real `mxscarnamod` would run the success path. The failure under
    // test is the missing tool, which is how this step errors in a normal
    // build.
    if std::process::Command::new("mxscarnamod").output().is_ok() {
        eprintln!("mxscarnamod is installed; RNA failure path not exercised");
        return;
    }
    let seqs = mafft_rs::SequenceSet {
        sequences: vec![
            mafft_rs::Sequence { name: "a".into(), data: b"acgtacgt".to_vec() },
            mafft_rs::Sequence { name: "b".into(), data: b"acgtacga".to_vec() },
        ],
        seq_type: mafft_rs::SeqType::Dna,
    };
    let err = mafft_rs::run_from_seqs(
        ["mafft-rs", "--quiet", "--nuc", "--qinsi", "--maxiterate", "0"],
        &seqs,
        &mafft_rs::SilentProgress,
    )
    .expect_err("RNA base-pair failure must return Err");
    assert_eq!(err.code(), 1, "{}", err.message());
    assert!(err.message().contains("mxscarnamod"), "{}", err.message());
}
