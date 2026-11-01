//! `--add` / `--addfragments` option combinations against the in-tree C
//! build: `--anysymbol` restore under `--keeplength` and the map options,
//! the existing alignment's all-gap columns under `--keeplength`, the
//! `--addfragments --compactmapout` coordinates, and `--addfragments`
//! with `--reorder` and the linkage options.
//!
//! References and the commands that produced them are in
//! `fixtures/add_followups/`. Where C writes `<addfile>.map`, the map is
//! compared too, not only stdout.

use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/add_followups").join(name)
}

fn scratch(tag: &str) -> PathBuf {
    static N: AtomicU64 = AtomicU64::new(0);
    let n = N.fetch_add(1, Ordering::Relaxed);
    let mut dir = std::env::temp_dir();
    dir.push(format!("mafft-rs-addfu-{}-{}-{tag}", std::process::id(), n));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Run `mafft-rs --quiet <flag> <case>.new.fa [extra...] <case>.existing.fa`
/// on copies of the inputs and compare stdout with `<tag>.aln`, and the
/// written map with `<tag>.map` when that reference exists.
fn expect(tag: &str, case: &str, flag: &str, extra: &[&str]) {
    let dir = scratch(tag);
    let ex = dir.join("existing.fa");
    let nw = dir.join("new.fa");
    std::fs::copy(fixture(&format!("{case}.existing.fa")), &ex).unwrap();
    std::fs::copy(fixture(&format!("{case}.new.fa")), &nw).unwrap();
    let mut argv: Vec<OsString> =
        vec!["mafft-rs".into(), "--quiet".into(), flag.into(), nw.clone().into()];
    argv.extend(extra.iter().map(OsString::from));
    argv.push(ex.into());
    let mut out = Vec::new();
    mafft_rs::run_from(argv, &mut out).unwrap_or_else(|e| panic!("{tag}: {e}"));

    let want = std::fs::read(fixture(&format!("{tag}.aln"))).unwrap();
    assert!(
        out == want,
        "{tag}: stdout differs from C (rust {} bytes, C {} bytes)",
        out.len(),
        want.len(),
    );

    let want_map = fixture(&format!("{tag}.map"));
    if want_map.exists() {
        let mut map_path = nw.into_os_string();
        map_path.push(".map");
        let got = std::fs::read(&map_path).unwrap_or_else(|e| panic!("{tag}: no map written: {e}"));
        let want = std::fs::read(&want_map).unwrap();
        assert!(
            got == want,
            "{tag}: map differs from C\n--- rust\n{}\n--- C\n{}",
            String::from_utf8_lossy(&got),
            String::from_utf8_lossy(&want),
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// back after the insertions are deleted (`disttbfast.c:4862`).
#[test]
fn add_keeplength_keeps_all_gap_columns() {
    expect("sample33.add.keeplength", "sample33", "--add", &["--keeplength"]);
}

/// The full map's columns are counted on the restored-width rows.
#[test]
fn add_mapout_matches_c() {
    expect("sample33.add.mapout", "sample33", "--add", &["--mapout"]);
}
