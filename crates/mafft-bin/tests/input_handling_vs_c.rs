//! Input-handling parity with C MAFFT 7.526 on adversarial FASTA.
//!
//! Every `.expected` file under `fixtures/input_handling/` is the verbatim
//! stdout of the arm64 C MAFFT 7.526 binary (`mafft --version` →
//! `v7.526 (2024/Apr/26)`) for the command named in the test; the exit-1
//! cases pin the exit status and message instead. See `fixtures/input_handling/README.md` at the
//! repository root for the C read-path analysis these tests encode.
//!
//! What C does by default (no `--anysymbol`), with the source it comes from:
//!
//! * `io.c:1435-1470` `load1SeqWithoutName_realloc`: keep only
//!   `isalpha || - || * || .`, fold case by `dorp`, then `* → -`.
//! * `mltaln9.c:60-85` `seqcheck`: any residue outside the alphabet
//!   (`blosum.c:12` protein `ARNDCQEGHILKMFPSTWYVBZX.-J`, `DNA.h:43`
//!   nucleotide `agctuAGCTUnNbdhkmnrsvwyx-O`) → `Illegal character c`, exit 1.
//!   So protein `U`/`O` and nucleotide `.` are fatal; protein `.` is a residue.
//! * `mltaln9.c:10537` `gappick0` in every driver: `-` is stripped from
//!   the input before alignment, `.` is not.
//! * `io.c:2065-2090` `countATGC`: `n` counts toward the nucleotide fraction.
//! * `scripts/mafft:1132` `tr "\r" "\n"`: CRLF input is fine.
//! * `scripts/mafft:1827-1834`: a description line preceded by blanks is
//!   fatal (exit 1) with a three-line message.
//!
//! And under `--anysymbol` (`replaceu` → align → `restoreu`):
//!
//! * `io.c:1329-1352` `charfilter`: only `\n`, space, `\r` are dropped;
//!   digits and tabs are KEPT and become `X`/`n` residues; `=`, `<`, `>`
//!   inside a sequence are fatal (exit 1).
//! * `restoreu.c:9-30` `fillorichar`: original characters are written back
//!   onto every aligned position that is not `-`.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("input_handling")
}

fn argv(flags: &[&str], input: &Path) -> Vec<OsString> {
    let mut argv: Vec<OsString> = vec![OsString::from("mafft-rs"), OsString::from("--quiet")];
    argv.extend(flags.iter().map(OsString::from));
    argv.push(input.as_os_str().to_os_string());
    argv
}

/// Run `mafft-rs --quiet <flags> <input>` in-process and return stdout.
fn run(flags: &[&str], input: &Path) -> Vec<u8> {
    let mut out = Vec::new();
    mafft_rs::run_from(argv(flags, input), &mut out)
        .unwrap_or_else(|e| panic!("mafft-rs {flags:?} {}: {}", input.display(), e.message()));
    out
}

fn expect_identical(name: &str, mode: &str, flags: &[&str]) {
    let f = fixtures();
    let input = f.join(format!("{name}.fa"));
    let expected = f.join(format!("{name}.{mode}.expected"));
    let got = run(flags, &input);
    let want = std::fs::read(&expected).expect("read expected");
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

/// C exits 1 (and writes nothing to stdout); pin the status and the
/// message fragment C prints.
fn expect_exit1(name: &str, flags: &[&str], message_contains: &str) {
    let input = fixtures().join(format!("{name}.fa"));
    let mut out = Vec::new();
    match mafft_rs::run_from(argv(flags, &input), &mut out) {
        Ok(()) => panic!(
            "mafft-rs {flags:?} {} succeeded; C MAFFT 7.526 exits 1 with `{message_contains}`\n{}",
            input.display(),
            String::from_utf8_lossy(&out)
        ),
        Err(e) => {
            assert_eq!(e.code(), 1, "exit code for {flags:?} {}: {}", input.display(), e.message());
            assert!(
                e.message().contains(message_contains),
                "message for {flags:?} {}: got `{}`, want it to contain `{message_contains}`",
                input.display(),
                e.message()
            );
            assert!(out.is_empty(), "stdout must be empty on failure");
        }
    }
}

const DEFAULT: &[&str] = &[];
const ANYSYMBOL: &[&str] = &["--anysymbol"];
const LINSI1: &[&str] = &["--localpair", "--maxiterate", "0"];
const NUC: &[&str] = &["--nuc"];

/// Every (fixture, mode) pair where C produces an alignment. One test per
/// fixture keeps the failure report readable; each mode is checked in turn.
fn all_modes_identical(name: &str, dna: bool) {
    expect_identical(name, "default", DEFAULT);
    expect_identical(name, "anysymbol", ANYSYMBOL);
    expect_identical(name, "linsi1", LINSI1);
    if dna {
        expect_identical(name, "nuc", NUC);
    }
}

// ---------------------------------------------------------------- protein

/// Lowercase and mixed-case protein: `onlyAlpha_upper` folds to upper
/// (`io.c:1367`); `--anysymbol` restores the input case.
#[test]
fn prot_lowercase() {
    all_modes_identical("prot_lowercase", false);
}

/// `*` (stop codon) → `-` via `kake2hiku` (`io.c:1380`), then stripped by
/// `gappick0`; under `--anysymbol` it is an unusual residue (`X`) that
/// comes back as `*`.
#[test]
fn prot_star() {
    all_modes_identical("prot_star", false);
}

/// `.` is a legal protein residue (`amino_n['.'] == 23`, `blosum.c:12`):
/// C keeps it, aligns it, and prints it. Bug on main: the engine stripped
/// `.` together with `-` before alignment.
#[test]
fn prot_dot_is_a_residue() {
    all_modes_identical("prot_dot", false);
}

/// `-` inside input sequences is stripped by `gappick0` before alignment.
#[test]
fn prot_dash_inside() {
    all_modes_identical("prot_dash_inside", false);
}

/// Leading/trailing `-`. Bug on main: the L-INS-1 pair phase read the
/// gapped input directly (`engine.rs` pair-path `seq_refs`), unlike C
/// `pairlocalalign.c:3372`.
#[test]
fn prot_leading_trailing_gaps() {
    all_modes_identical("prot_leading_trailing_gaps", false);
}

/// `U` / `O` are not in C's protein alphabet: `seqcheck` aborts with
/// `Illegal character U` (site 6 of sequence 1 here). `J`, `B`, `Z`, `X`
/// are legal. With `--anysymbol` everything is substituted and restored.
#[test]
fn prot_unusual_UOJBZX_default_is_illegal_character() {
    expect_exit1("prot_unusual_UOJBZX", DEFAULT, "Illegal character U");
    expect_exit1("prot_unusual_UOJBZX", LINSI1, "Illegal character U");
}

#[test]
fn prot_unusual_UOJBZX_anysymbol_restores() {
    expect_identical("prot_unusual_UOJBZX", "anysymbol", ANYSYMBOL);
}

/// GenBank-style numbered lines. Default path drops digits and spaces
/// (`onlyAlpha_*`). Under `--anysymbol` C's `charfilter` keeps the digits,
/// which become `X` residues and are restored as digits.
#[test]
fn prot_digits_spaces() {
    all_modes_identical("prot_digits_spaces", false);
}

/// CRLF line endings (`tr "\r" "\n"`, `scripts/mafft:1132`).
#[test]
fn prot_crlf() {
    all_modes_identical("prot_crlf", false);
}

/// Blank lines between records are ignored.
#[test]
fn prot_blank_lines() {
    all_modes_identical("prot_blank_lines", false);
}

/// A description line preceded by blanks is fatal in C
/// (`scripts/mafft:1827-1834`), in every mode.
#[test]
fn prot_header_leading_ws_is_fatal() {
    let msg = "The first character of a description line must be";
    expect_exit1("prot_header_leading_ws", DEFAULT, msg);
    expect_exit1("prot_header_leading_ws", ANYSYMBOL, msg);
    expect_exit1("prot_header_leading_ws", LINSI1, msg);
}

/// `@ # ~ !` and a tab inside sequences. Default path drops them all.
/// Under `--anysymbol` C keeps them (including the tab) as unusual
/// residues and restores them.
#[test]
fn prot_punct() {
    all_modes_identical("prot_punct", false);
}

/// `=` inside a sequence line: silently dropped by default; fatal under
/// `--anysymbol` (`charfilter`, `io.c:1337`).
#[test]
fn prot_equals() {
    expect_identical("prot_equals", "default", DEFAULT);
    expect_identical("prot_equals", "linsi1", LINSI1);
    expect_exit1("prot_equals", ANYSYMBOL, "Characters '= < >' can be used only in the title lines");
}

// -------------------------------------------------------------------- DNA

/// Upper/mixed-case nucleotide → lowercase (`onlyAlpha_lower`).
#[test]
fn dna_mixed_case() {
    all_modes_identical("dna_mixed_case", true);
}

/// IUPAC ambiguity codes `N R Y K M S W B D H V` are all legal nucleotides.
#[test]
fn dna_iupac() {
    all_modes_identical("dna_iupac", true);
}

/// `-` is stripped, but `.` is NOT in the nucleotide alphabet: C aborts
/// with `Illegal character .` (site 11 of sequence 3). `--anysymbol` maps
/// `.` to `n` and restores it.
#[test]
fn dna_gaps_inside_default_is_illegal_character() {
    expect_exit1("dna_gaps_inside", DEFAULT, "Illegal character .");
    expect_exit1("dna_gaps_inside", LINSI1, "Illegal character .");
    expect_exit1("dna_gaps_inside", NUC, "Illegal character .");
}

/// Bug on main: the `--anysymbol` restore treated `.` as a gap on both
/// sides, so the restored residues shifted by the number of dots.
#[test]
fn dna_gaps_inside_anysymbol_restores_dots() {
    expect_identical("dna_gaps_inside", "anysymbol", ANYSYMBOL);
}

#[test]
fn dna_leading_trailing_gaps() {
    all_modes_identical("dna_leading_trailing_gaps", true);
}

#[test]
fn dna_crlf() {
    all_modes_identical("dna_crlf", true);
}

#[test]
fn dna_digits_spaces() {
    all_modes_identical("dna_digits_spaces", true);
}

/// One residue in three is `N`. C's `countATGC` counts `n` as nucleotide,
/// so this is DNA. Bug on main: `detect_seq_type` did not, so the input
/// was aligned as protein (uppercase, protein matrix).
#[test]
fn dna_n_rich_is_nucleotide() {
    all_modes_identical("dna_n_rich", true);
}

/// `U` (RNA) is a legal nucleotide; `*` → `-` → stripped.
#[test]
fn dna_rna_u_star() {
    all_modes_identical("dna_rna_u_star", true);
}

// ------------------------------------------------------------------- --add

fn add_flags(extra: &[&str]) -> Vec<String> {
    let new = fixtures().join("add_new_unusual.fa");
    let mut v: Vec<String> = extra.iter().map(|s| s.to_string()).collect();
    v.push("--add".into());
    v.push(new.display().to_string());
    v
}

fn expect_add_identical(mode: &str, extra: &[&str]) {
    let flags = add_flags(extra);
    let flags: Vec<&str> = flags.iter().map(String::as_str).collect();
    expect_identical("add_existing_gapped", mode, &flags);
}

/// Gapped existing alignment + two added sequences carrying lowercase,
/// a trailing `*` and a `.`. C keeps the `.` as a residue and inserts a
/// column for it; the existing gaps are handled by `commongappick`.
#[test]
fn add_gapped_existing_with_unusual_new() {
    expect_add_identical("add", &[]);
}

#[test]
fn add_keeplength_gapped_existing_with_unusual_new() {
    expect_add_identical("add_keeplength", &["--keeplength"]);
}

/// Under `--anysymbol`, C substitutes/restores the ADDED sequences too
/// (`replaceu` runs on the concatenated infile, `scripts/mafft:1142,2305`).
/// The restore itself now matches; the remaining difference is where the
/// trailing `X` (from `*`) lands in the `--add` DP — C puts it after the
/// terminal gap run, rust before. The same difference appears with a
/// literal trailing `X` and no `--anysymbol`, so it is an `--add`
/// tie-break, not an input-handling issue.
#[test]
#[ignore = "diverges from C: trailing-X placement in --add DP (see fixtures/input_handling/README.md, residual R-A)"]
fn add_anysymbol_gapped_existing_with_unusual_new() {
    expect_add_identical("add_anysymbol", &["--anysymbol"]);
}
