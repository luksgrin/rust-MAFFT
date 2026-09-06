//! `run_from_seqs` must return exactly the rows `run_from` prints.
//!
//! The in-memory entry point shares every stage after input parsing with
//! the FASTA path, so the two can only differ in how the input is brought
//! into the engine. Each test therefore runs the same flags through both
//! paths on the same sequences and compares (name, row) lists: the FASTA
//! path's bytes are parsed back with a deliberately naive parser (no
//! normalisation, no case fold) so nothing can mask a difference.
//!
//! The in-memory side is exercised in two forms of the same file:
//!
//! * `canonical` — the set `read_fasta` produces, which already carries
//!   the reader's normalisation and a detected type. This is the zero-copy
//!   route (`Cow::Borrowed`).
//! * `raw` — the file's residues verbatim (whatever case the file used)
//!   with `SeqType::Unknown`, so the entry point has to normalise, detect
//!   and fold exactly as the reader would.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use mafft_rs::{
    run_from_seqs, run_from_with_progress, Mafft, SeqType, Sequence, SequenceSet, SilentProgress,
};

type Rows = Vec<(String, Vec<u8>)>;

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures")
}

/// C MAFFT's 36-sequence protein sample, from the `mafft-upstream` submodule.
fn upstream_sample() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../mafft-upstream/test/sample")
        .canonicalize()
        .expect("mafft-upstream submodule checked out")
}

fn argv(flags: &[&str]) -> Vec<OsString> {
    let mut argv = vec![OsString::from("mafft-rs"), OsString::from("--quiet")];
    argv.extend(flags.iter().map(OsString::from));
    argv
}

/// Header text after `>` is the name (minus a CRLF's `\r`, which the
/// reader's `lines()` also drops); data lines are concatenated verbatim,
/// `\r` included.
fn parse_fasta_raw(bytes: &[u8]) -> Rows {
    let mut rows: Rows = Vec::new();
    for line in bytes.split(|&b| b == b'\n') {
        if let Some(name) = line.strip_prefix(b">") {
            let name = name.strip_suffix(b"\r").unwrap_or(name);
            rows.push((String::from_utf8(name.to_vec()).unwrap(), Vec::new()));
        } else if let Some((_, data)) = rows.last_mut() {
            data.extend_from_slice(line);
        }
    }
    rows
}

fn fasta_path_rows(flags: &[&str], input: &Path) -> Rows {
    let mut argv = argv(flags);
    argv.push(input.as_os_str().to_os_string());
    let mut out = Vec::new();
    run_from_with_progress(argv, &mut out, &SilentProgress)
        .unwrap_or_else(|e| panic!("run_from {flags:?} {}: {}", input.display(), e.message()));
    parse_fasta_raw(&out)
}

fn memory_path_rows(flags: &[&str], set: &SequenceSet) -> Rows {
    let msa = run_from_seqs(argv(flags), set, &SilentProgress)
        .unwrap_or_else(|e| panic!("run_from_seqs {flags:?}: {}", e.message()));
    assert_eq!(msa.names.len(), msa.sequences.len());
    msa.names.into_iter().zip(msa.sequences).collect()
}

/// The file's residues exactly as written, typed `Unknown`.
fn raw_set(input: &Path) -> SequenceSet {
    let bytes = std::fs::read(input).expect("read fixture");
    SequenceSet {
        sequences: parse_fasta_raw(&bytes)
            .into_iter()
            .map(|(name, data)| Sequence { name, data })
            .collect(),
        seq_type: SeqType::Unknown,
    }
}

fn variants(input: &Path) -> Vec<(&'static str, SequenceSet)> {
    vec![
        ("canonical", mafft_io::read_fasta(input).expect("read fixture")),
        ("raw", raw_set(input)),
    ]
}

fn assert_parity(flags: &[&str], input: &Path) {
    let want = fasta_path_rows(flags, input);
    assert!(!want.is_empty());
    for (label, set) in variants(input) {
        let got = memory_path_rows(flags, &set);
        assert_rows_equal(&got, &want, &format!("{flags:?} on {} ({label})", input.display()));
    }
}

fn assert_rows_equal(got: &Rows, want: &Rows, what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: row count");
    for (i, (g, w)) in got.iter().zip(want).enumerate() {
        assert_eq!(g.0, w.0, "{what}: name of row {i}");
        assert!(
            g.1 == w.1,
            "{what}: row {i} ({})\n memory: {}\n  fasta: {}",
            g.0,
            String::from_utf8_lossy(&g.1),
            String::from_utf8_lossy(&w.1)
        );
    }
}

// --- protein: the 36-sequence upstream sample --------------------------

#[test]
fn sample_default_fftns2() {
    assert_parity(&[], &upstream_sample());
}

#[test]
fn sample_reorder() {
    assert_parity(&["--reorder"], &upstream_sample());
}

#[test]
fn sample_amino_forced() {
    assert_parity(&["--amino"], &upstream_sample());
}

#[test]
fn sample_bl50_retree2_progressive() {
    assert_parity(&["--bl", "50", "--retree", "2", "--maxiterate", "0"], &upstream_sample());
}

#[test]
fn sample_localpair_progressive() {
    assert_parity(&["--localpair", "--maxiterate", "0"], &upstream_sample());
}

// --- nucleotide: the 120 x 1400 bp CDS fixture --------------------------

#[test]
fn mtb_auto_adjustdirection_thread1_nuc() {
    // The exact command line panaroo-rs issues per gene cluster.
    assert_parity(
        &["--auto", "--adjustdirection", "--thread", "1", "--nuc"],
        &fixtures().join("mtb_cds_120x1400.fa"),
    );
}

#[test]
fn mtb_retree2_maxiterate2() {
    assert_parity(&["--retree", "2", "--maxiterate", "2"], &fixtures().join("mtb_cds_120x1400.fa"));
}

#[test]
fn mtb_nuc_reorder() {
    assert_parity(&["--nuc", "--reorder"], &fixtures().join("mtb_cds_120x1400.fa"));
}

// --- small nucleotide fixtures ----------------------------------------

#[test]
fn tiny_pair_auto_and_localpair() {
    assert_parity(&["--auto"], &fixtures().join("refine_njob2_min.fa"));
    assert_parity(&["--localpair"], &fixtures().join("dna_pair_gapscale_min.fa"));
}

// --- adjustdirection actually flipping a strand ------------------------

const DNA: [&str; 4] = [
    "ATGGCTAGCTTGGACCATTGCAGGTACCCATGGAACTTGGGCCATTAGGCATTGACCTAGATTGACC",
    "ATGGCTAGCTTGGACCATTGCAGGTACCCTTGGAACTTGGCCATTAGGCATTGACCTAGGATTGACC",
    "ATGGCAAGCTTAGACCTTTGCAGGTACGCATGGAACTAGGGCCTTTAGGCATTGACCTAGATTGACC",
    "TTGGCTAGCTTGGACCATTGCAGCTACCCATGGAACTTGGGCCATTAGGCTTTGACGTAGATTGACC",
];

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' => 'A',
            b'G' => 'C',
            b'C' => 'G',
            other => other as char,
        })
        .collect()
}

fn write_tmp_fasta(tag: &str, records: &[(&str, &str)]) -> PathBuf {
    let mut path = std::env::temp_dir();
    path.push(format!("mafft-rs-inmemory-{}-{tag}.fa", std::process::id()));
    let mut text = String::new();
    for (name, seq) in records {
        text.push_str(&format!(">{name}\n{seq}\n"));
    }
    std::fs::write(&path, text).expect("write temp fasta");
    path
}

#[test]
fn reverse_complemented_input_gets_the_r_prefix_on_both_paths() {
    let flipped = revcomp(DNA[2]);
    let records = [("a", DNA[0]), ("b", DNA[1]), ("c flipped", flipped.as_str()), ("d", DNA[3])];
    let path = write_tmp_fasta("revcomp", &records);
    for flags in [
        &["--adjustdirection", "--nuc"][..],
        &["--adjustdirection", "--reorder"],
        &["--adjustdirectionaccurately", "--nuc"],
    ] {
        let want = fasta_path_rows(flags, &path);
        assert!(
            want.iter().any(|(name, _)| name == "_R_c flipped"),
            "{flags:?}: the FASTA path must flip sequence c: {:?}",
            want.iter().map(|(n, _)| n).collect::<Vec<_>>()
        );
        for (label, set) in variants(&path) {
            let got = memory_path_rows(flags, &set);
            assert_rows_equal(&got, &want, &format!("{flags:?} revcomp ({label})"));
        }
    }
    std::fs::remove_file(&path).ok();
}

// --- case handling ----------------------------------------------------

#[test]
fn preservecase_keeps_the_callers_case_on_both_paths() {
    let mixed: Vec<String> = DNA
        .iter()
        .enumerate()
        .map(|(i, s)| if i % 2 == 0 { s.to_ascii_lowercase() } else { s.to_string() })
        .collect();
    let records: Vec<(&str, &str)> =
        ["a", "b", "c", "d"].iter().zip(&mixed).map(|(n, s)| (*n, s.as_str())).collect();
    let path = write_tmp_fasta("preservecase", &records);
    // Only the raw variant is the same input here: `read_fasta` has already
    // folded the case away, so the "canonical" set is a different file.
    let raw = raw_set(&path);
    for flags in [&["--preservecase"][..], &["--anysymbol", "--nuc"]] {
        let want = fasta_path_rows(flags, &path);
        let got = memory_path_rows(flags, &raw);
        assert_rows_equal(&got, &want, &format!("{flags:?} preservecase (raw)"));
    }
    // And the case really was kept.
    let rows = memory_path_rows(&["--preservecase"], &raw);
    assert!(rows[0].1.iter().any(|c| c.is_ascii_lowercase()));
    assert!(rows[1].1.iter().any(|c| c.is_ascii_uppercase()));
    std::fs::remove_file(&path).ok();
}

#[test]
fn declared_type_is_honoured_like_the_matching_force_flag() {
    // A declared `seq_type` stands in for `--nuc` / `--amino`, so aligning
    // a set typed `Protein` must equal aligning the same residues with
    // `--amino`; both differ from the nucleotide alignment.
    let path = write_tmp_fasta("declared", &[("a", DNA[0]), ("b", DNA[1]), ("c", DNA[2])]);
    let as_protein = SequenceSet { seq_type: SeqType::Protein, ..raw_set(&path) };
    let got = memory_path_rows(&[], &as_protein);
    let want = fasta_path_rows(&["--amino"], &path);
    assert_rows_equal(&got, &want, "declared Protein vs --amino");
    assert_ne!(got, fasta_path_rows(&[], &path), "protein and nucleotide runs must differ");
    // ... and the flag still wins over the declaration.
    let forced = memory_path_rows(&["--nuc"], &as_protein);
    assert_rows_equal(&forced, &fasta_path_rows(&["--nuc"], &path), "--nuc over declared Protein");
    std::fs::remove_file(&path).ok();
}

// --- builder and output flag -----------------------------------------

#[test]
fn builder_run_seqs_matches_run_from_seqs() {
    let set = mafft_io::read_fasta(fixtures().join("mtb_cds_120x1400.fa")).unwrap();
    let via_fn = memory_path_rows(&["--auto", "--adjustdirection", "--thread", "1", "--nuc"], &set);
    let msa = Mafft::new()
        .quiet()
        .auto()
        .adjust_direction()
        .thread(1)
        .nuc()
        .progress(SilentProgress)
        .run_seqs(&set)
        .expect("builder run_seqs");
    let via_builder: Rows = msa.names.into_iter().zip(msa.sequences).collect();
    assert_rows_equal(&via_builder, &via_fn, "builder vs run_from_seqs");
}

#[test]
fn output_flag_still_writes_the_file_and_the_rows_come_back() {
    let set = mafft_io::read_fasta(fixtures().join("refine_njob2_min.fa")).unwrap();
    let mut out_path = std::env::temp_dir();
    out_path.push(format!("mafft-rs-inmemory-{}-output.fa", std::process::id()));
    let out_str = out_path.to_str().unwrap().to_string();
    let msa = run_from_seqs(argv(&["--output", &out_str]), &set, &SilentProgress).unwrap();
    let written = std::fs::read(&out_path).expect("--output file written");
    let mut expected = Vec::new();
    let mut argv_file = argv(&[]);
    argv_file.push(fixtures().join("refine_njob2_min.fa").into_os_string());
    run_from_with_progress(argv_file, &mut expected, &SilentProgress).unwrap();
    assert_eq!(written, expected, "--output must hold what the CLI prints");
    assert_eq!(msa.names.len(), 2);
    std::fs::remove_file(&out_path).ok();
}

// --- errors -----------------------------------------------------------

fn dna_set() -> SequenceSet {
    SequenceSet {
        sequences: DNA
            .iter()
            .enumerate()
            .map(|(i, s)| Sequence { name: format!("s{i}"), data: s.as_bytes().to_vec() })
            .collect(),
        seq_type: SeqType::Unknown,
    }
}

#[test]
fn an_input_path_in_argv_is_rejected() {
    let err = run_from_seqs(argv(&["some.fa"]), &dna_set(), &SilentProgress)
        .expect_err("INPUT alongside in-memory sequences must fail");
    assert_eq!(err.code(), 1);
    assert!(err.message().starts_with("run_from_seqs: "), "{}", err.message());
    assert!(err.message().contains("some.fa"), "{}", err.message());
}

#[test]
fn errors_match_run_from() {
    // clap failure: same code, same rendered text.
    let file = run_from_with_progress(argv(&["--no-such-flag"]), &mut Vec::new(), &SilentProgress)
        .unwrap_err();
    let mem = run_from_seqs(argv(&["--no-such-flag"]), &dna_set(), &SilentProgress).unwrap_err();
    assert_eq!(mem.code(), 2);
    assert_eq!(mem.code(), file.code());
    assert_eq!(mem.message(), file.message());

    // Pre-input validation.
    let mem = run_from_seqs(argv(&["--nodeout", "--maxiterate", "5"]), &dna_set(), &SilentProgress)
        .unwrap_err();
    assert_eq!(mem.code(), 1);
    assert_eq!(
        mem.message(),
        "The --nodeout option supports only progressive method (--maxiterate 0) for now."
    );

    // Post-input validation.
    let mem = run_from_seqs(argv(&["--memsave", "--localpair"]), &dna_set(), &SilentProgress)
        .unwrap_err();
    assert_eq!((mem.code(), mem.message()), (1, "Impossible"));

    // Successful early exits are reported the way `--pdbidlist` is.
    let mem = run_from_seqs(argv(&["--pdbidlist", "x"]), &dna_set(), &SilentProgress).unwrap_err();
    assert_eq!(mem.code(), 0);
    let cite = run_from_seqs(argv(&["--cite"]), &dna_set(), &SilentProgress).unwrap_err();
    assert_eq!(cite.code(), 0);
    assert!(cite.message().contains("Katoh"), "{}", cite.message());
}

#[test]
fn an_empty_set_is_reported_like_an_empty_file() {
    let empty = SequenceSet { sequences: Vec::new(), seq_type: SeqType::Dna };
    let err = run_from_seqs(argv(&[]), &empty, &SilentProgress).unwrap_err();
    assert_eq!(err.code(), 1);
    assert_eq!(err.message(), "Error: no sequences found in input");
}

// --- the input-handling corpus: both paths, success and failure alike --

/// What a run produces: the rows, or the `MafftError` (code, message).
type Outcome = Result<Rows, (i32, String)>;

fn corpus_files() -> Vec<PathBuf> {
    let dir = fixtures().join("input_handling");
    let mut files: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("input_handling corpus")
        .map(|e| e.expect("dir entry").path())
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .collect();
    files.sort();
    assert_eq!(files.len(), 22, "corpus size");
    files
}

fn fasta_path_outcome(mut argv: Vec<OsString>, input: &Path) -> Outcome {
    argv.push(input.as_os_str().to_os_string());
    let mut out = Vec::new();
    match run_from_with_progress(argv, &mut out, &SilentProgress) {
        Ok(()) => Ok(parse_fasta_raw(&out)),
        Err(e) => {
            assert!(out.is_empty(), "{}: stdout must be empty on failure", input.display());
            Err((e.code(), e.message().to_string()))
        }
    }
}

fn memory_path_outcome(argv: Vec<OsString>, set: &SequenceSet) -> Outcome {
    run_from_seqs(argv, set, &SilentProgress)
        .map(|msa| msa.names.into_iter().zip(msa.sequences).collect())
        .map_err(|e| (e.code(), e.message().to_string()))
}

fn assert_outcomes_equal(got: &Outcome, want: &Outcome, what: &str) {
    match (got, want) {
        (Ok(g), Ok(w)) => assert_rows_equal(g, w, what),
        (Err(g), Err(w)) => assert_eq!(g, w, "{what}: error"),
        (Ok(_), Err(w)) => panic!("{what}: memory path succeeded, fasta path failed with {w:?}"),
        (Err(g), Ok(_)) => panic!("{what}: memory path failed with {g:?}, fasta path succeeded"),
    }
}

/// Every corpus file, through both paths, with the reader the CLI would
/// use for `flags` (the `Cow::Borrowed` route) and as raw bytes typed
/// `Unknown`. Files C refuses must fail on both paths with the same code
/// and message: the alphabet check (`seqcheck`) and the case-preserving
/// reader's `= < >` rule are applied to in-memory input exactly as to a
/// file.
///
/// The one thing an in-memory caller cannot reproduce is a FASTA framing
/// error — a blank before `>` is not a residue — so when the reader
/// rejects a file, the raw variant is skipped and the canonical variant
/// asserts the CLI reported that reader error.
fn assert_corpus_parity(flags: &[&str]) {
    let casepreserve = flags.contains(&"--anysymbol") || flags.contains(&"--preservecase");
    for input in corpus_files() {
        let what = format!("{flags:?} on {}", input.file_name().unwrap().to_string_lossy());
        let want = fasta_path_outcome(argv(flags), &input);

        let read = if casepreserve {
            mafft_io::read_fasta_casepreserve(&input)
        } else {
            mafft_io::read_fasta(&input)
        };
        match read {
            Ok(set) => {
                let got = memory_path_outcome(argv(flags), &set);
                assert_outcomes_equal(&got, &want, &format!("{what} (canonical)"));
            }
            Err(e) => {
                let (code, msg) = want
                    .as_ref()
                    .expect_err(&format!("{what}: the reader refused the file ({e}) but the CLI aligned it"));
                assert_eq!(*code, 1, "{what}");
                assert!(
                    msg.ends_with(&e.to_string()),
                    "{what}: CLI error `{msg}` is not the reader's `{e}`"
                );
                // The raw variant would have to invent a framing for the
                // file; the failure is not about residues.
                if msg.starts_with("The first character of a description line") {
                    continue;
                }
            }
        }

        let got = memory_path_outcome(argv(flags), &raw_set(&input));
        assert_outcomes_equal(&got, &want, &format!("{what} (raw)"));
    }
}

#[test]
fn corpus_default() {
    assert_corpus_parity(&[]);
}

#[test]
fn corpus_anysymbol() {
    assert_corpus_parity(&["--anysymbol"]);
}

#[test]
fn corpus_localpair_progressive() {
    assert_corpus_parity(&["--localpair", "--maxiterate", "0"]);
}

#[test]
fn corpus_nuc_forced() {
    assert_corpus_parity(&["--nuc"]);
}

/// The exit-1 cases spelled out, quiet and with C's banner (which names
/// the site and sequence, so it must be computed on the same residues).
#[test]
fn corpus_illegal_input_errors_are_identical() {
    let cases: [(&str, &[&str], &str); 5] = [
        ("prot_unusual_UOJBZX", &[], "Illegal character U"),
        ("prot_unusual_UOJBZX", &["--localpair", "--maxiterate", "0"], "Illegal character U"),
        ("dna_gaps_inside", &[], "Illegal character ."),
        ("dna_gaps_inside", &["--nuc"], "Illegal character ."),
        ("prot_equals", &["--anysymbol"], "Characters '= < >' can be used only in the title lines"),
    ];
    for (name, flags, fragment) in cases {
        let path = fixtures().join("input_handling").join(format!("{name}.fa"));
        for quiet in [true, false] {
            let mut argv = vec![OsString::from("mafft-rs")];
            if quiet {
                argv.push(OsString::from("--quiet"));
            }
            argv.extend(flags.iter().map(OsString::from));
            let what = format!("{name} {flags:?} quiet={quiet}");
            let want = fasta_path_outcome(argv.clone(), &path);
            let (code, msg) = want.as_ref().expect_err(&format!("{what}: C exits 1"));
            assert_eq!(*code, 1, "{what}");
            assert!(msg.contains(fragment), "{what}: `{msg}`");
            if !quiet && fragment.starts_with("Illegal") {
                assert!(msg.contains("Please check site"), "{what}: banner expected: `{msg}`");
            }
            let got = memory_path_outcome(argv, &raw_set(&path));
            assert_eq!(got, want, "{what}: memory path must fail identically");
        }
    }
}

/// A caller that declares the type is held to that type's alphabet, like
/// `--nuc` / `--amino` on a file: protein `U` is fatal, nucleotide `.` is
/// fatal, and each is legal under the other type.
#[test]
fn declared_type_is_seqchecked_like_the_force_flag() {
    let corpus = fixtures().join("input_handling");
    for (name, declared, flag) in [
        ("prot_unusual_UOJBZX", SeqType::Protein, "--amino"),
        ("dna_gaps_inside", SeqType::Dna, "--nuc"),
    ] {
        let path = corpus.join(format!("{name}.fa"));
        let want = fasta_path_outcome(argv(&[flag]), &path);
        assert!(want.is_err(), "{name} {flag} must fail");
        let set = SequenceSet { seq_type: declared, ..raw_set(&path) };
        assert_eq!(memory_path_outcome(argv(&[]), &set), want, "{name} declared {declared:?}");
    }
    // `.` is a legal protein residue: the DNA file with dots aligns when
    // declared Protein, exactly as with `--amino`.
    let path = corpus.join("dna_gaps_inside.fa");
    let want = fasta_path_outcome(argv(&["--amino"]), &path);
    assert!(want.is_ok(), "dna_gaps_inside --amino: {want:?}");
    let set = SequenceSet { seq_type: SeqType::Protein, ..raw_set(&path) };
    assert_outcomes_equal(&memory_path_outcome(argv(&[]), &set), &want, "dots as protein");
}
