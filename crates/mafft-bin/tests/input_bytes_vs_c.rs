//! Bytes the C wrapper accepts and `mafft-rs` used to refuse.
//!
//! The in-tree `scripts/mafft` runs `tr "\r" "\n"` before the engines
//! (`scripts/mafft:1132`) and the engines read bytes, not UTF-8. These
//! tests write the same file for both programs and require the same
//! stdout. A Latin-1 byte inside a sequence is not a letter, so C drops
//! it and aligns what remains; the alignment, not merely a successful
//! exit, is what is compared.
//!
//! The C command is run with `LANG=C` and `LC_ALL=C`. The wrapper exports
//! `LANG=C`, but a UTF-8 `LC_ALL` from the caller overrides that, and
//! BSD `tr` then stops at the first byte that is not valid UTF-8
//! (`tr: Illegal byte sequence`) and leaves a truncated file. Under the
//! C locale `tr` is the byte translation the wrapper is written as: a
//! non-UTF-8 name byte is kept, and a non-letter inside a sequence is
//! dropped.

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::process::Command;

fn worktree() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("worktree root")
}

fn c_mafft(input: &Path) -> Vec<u8> {
    let root = worktree();
    let script = root.join("mafft-upstream/scripts/mafft");
    let binaries = root.join("mafft-upstream/binaries");
    assert!(script.is_file(), "C wrapper missing at {}", script.display());
    assert!(
        binaries.join("disttbfast").is_file(),
        "in-tree C build missing at {} (make -C mafft-upstream/core)",
        binaries.display()
    );
    let out = Command::new(&script)
        .args(["--quiet", "--thread", "1"])
        .arg(input)
        .env("MAFFT_BINARIES", &binaries)
        // BSD `tr` in a UTF-8 locale stops at the first illegal byte and
        // truncates the file. `LC_ALL` overrides the wrapper's `LANG=C`.
        .env("LANG", "C")
        .env("LC_ALL", "C")
        .output()
        .unwrap_or_else(|e| panic!("run {}: {e}", script.display()));
    assert!(
        out.status.success(),
        "C mafft failed ({})\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn rust_mafft(input: &Path) -> Vec<u8> {
    let mut argv = vec![
        OsString::from("mafft-rs"),
        OsString::from("--quiet"),
        OsString::from("--thread"),
        OsString::from("1"),
    ];
    argv.push(input.as_os_str().to_os_string());
    let mut stdout = Vec::new();
    mafft_rs::run_from(argv, &mut stdout)
        .unwrap_or_else(|e| panic!("mafft-rs: {}", e.message()));
    stdout
}

fn assert_same_alignment(label: &str, bytes: &[u8]) {
    let dir = std::env::temp_dir().join(format!(
        "mafft-rs-input-bytes-{}-{}",
        label,
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("in.fa");
    std::fs::write(&path, bytes).unwrap();
    let c = c_mafft(&path);
    let rust = rust_mafft(&path);
    let _ = std::fs::remove_dir_all(&dir);
    assert!(
        !c.is_empty() && c.starts_with(b">"),
        "{label}: C produced no alignment\n{}",
        String::from_utf8_lossy(&c)
    );
    if rust != c {
        panic!(
            "{label}: mafft-rs stdout differs from the in-tree C build\n--- C ---\n{}--- rust ---\n{}",
            String::from_utf8_lossy(&c),
            String::from_utf8_lossy(&rust)
        );
    }
}

/// Latin-1 e-acute between protein residues. C drops the byte
/// (`onlyAlpha_*`) and aligns the rest. The header also carries a
/// non-UTF-8 byte, which both sides keep.
#[test]
fn non_utf8_sequence_matches_c() {
    let mut fasta = b">s1 \xff\nMNGTEGDNFYVPFSNKTGLA".to_vec();
    fasta.push(0xe9);
    fasta.extend_from_slice(b"RSPYEY\n>s2\nMNGTEGDNFYVPFSNKTGVARSPYEY\n");
    assert_same_alignment("nonutf8", &fasta);
}

/// Two records separated only by CR. The wrapper's `tr` makes that two
/// sequences; a reader that splits only on LF would keep one record.
#[test]
fn cr_only_fasta_matches_c() {
    let fasta = b">s1\rMNGTEGDNFYVPFSNKTGLARSPYEY\r>s2\rMNGTEGDNFYVPFSNKTGVARSPYEY\r";
    assert_same_alignment("cronly", fasta);
}
