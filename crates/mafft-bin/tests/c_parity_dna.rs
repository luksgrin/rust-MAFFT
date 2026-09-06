//! Byte-for-byte parity with C MAFFT 7.526 on nucleotide input.
//!
//! Every `.expected` file here is the verbatim stdout of C MAFFT 7.526 for
//! the command named in the test, so these tests need no C installation.
//! The files were first captured from an x86-64 build and then verified
//! byte-identical (2026-09-06) against the reference build defined in
//! `docs/architecture/byte-identity.md`: the pinned `mafft-upstream` source
//! compiled in-tree with the upstream Makefile's default flags, on arm64.
//! None of these inputs is sensitive to the floating-point contraction
//! policy, so a single variant serves both platforms; if that ever changes,
//! add the fixture to `crates/mafft-core/tests/fixtures/policy_fixtures.tsv`
//! and regenerate with `scripts/regen_policy_fixtures.sh`.
//!
//! They pin two fixes and the class of input that detects them:
//!
//! * `dna_pair_gapscale_min` — the smallest input on which the pair-phase
//!   gap penalties for nucleotide were a third of C's (missing
//!   `constants.c:316-322` `3 *`). Two 15 bp sequences; C leaves the pair
//!   ungapped, the old code opened a gap.
//! * `refine_njob2_min` — the smallest input showing that C refines a
//!   *pair* (`dvtditr.c:704-708`: `njob == 2 → weight = 0; niter = 1`),
//!   where the old code skipped refinement at `nseq <= 2`.
//! * `r2_dna_clusters/` — 30 clusters of 2–6 sequences, 300–900 bp,
//!   85–99 % identity, evolved from real biological ancestors with
//!   substitutions and indels. Uniform-random or fixture-only corpora did
//!   not reliably surface the gap-scale bug; this shape did (16/30 before
//!   the fix), so it is the regression guard for that whole class.
//! * `panaroo_tiny_dna_clusters/` — five real 4-sequence gene clusters from
//!   the panaroo-rs tiny/core parity corpus (issue #1), the exact workload
//!   panaroo-rs hands to MAFFT: `--auto --adjustdirection --thread 1 --nuc`.
//!   See the README in that directory for provenance and regeneration.

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

/// Bug: nucleotide pair-phase gap penalties lacked C's `3 *`
/// (`constants.c:316-322`), so L-INS-1 opened a gap C refuses.
#[test]
fn linsi1_dna_pair_gap_scale_matches_c() {
    let f = fixtures();
    expect_identical(
        &["--localpair", "--maxiterate", "0"],
        &f.join("dna_pair_gapscale_min.fa"),
        &f.join("dna_pair_gapscale_min.linsi1.expected"),
    );
}

/// Bug: refinement was skipped at `nseq <= 2`; C runs it once, unweighted
/// (`dvtditr.c:704-708`, `tditeration.c:772`). FFT-NS-i and `--auto`
/// (L-INS-i) both exercise it.
#[test]
fn fftnsi_refines_a_pair_like_c() {
    let f = fixtures();
    expect_identical(
        &["--maxiterate", "1000"],
        &f.join("refine_njob2_min.fa"),
        &f.join("refine_njob2_min.fftnsi.expected"),
    );
}

#[test]
fn auto_refines_a_pair_like_c() {
    let f = fixtures();
    expect_identical(
        &["--auto"],
        &f.join("refine_njob2_min.fa"),
        &f.join("refine_njob2_min.auto.expected"),
    );
}

/// The 120-sequence FFT-NS-i case that exposed the `dndpre` offset bug.
///
/// C's `dndpre` writes the refinement-tree distances with its DEFAULT
/// `poffset`, which differs by alphabet (`DEFAULTOFS_N = -369` vs
/// `DEFAULTOFS_B = -123`). Using the protein shift for DNA changed those
/// distances, which swapped two UPGMA merges, which changed the group on 131
/// of 237 refinement branches. This input is large enough for `--auto` to
/// select FFT-NS-i, which is the regime that matters when aligning a gene
/// family across a few hundred genomes.
///
/// Run without `--thread`: C takes its single-threaded `TreeDependentIteration`
/// path there. `fftnsi_120seq_dna_matches_c_under_thread_1` below pins the
/// `athread` path, which C selects for `--thread N >= 1`.
#[test]
fn fftnsi_120seq_dna_matches_c() {
    let f = fixtures();
    expect_identical(
        &["--adjustdirection", "--nuc", "--retree", "2", "--maxiterate", "2"],
        &f.join("mtb_cds_120x1400.fa"),
        &f.join("mtb_cds_120x1400.fftnsi.expected"),
    );
}

/// The same input under `--thread 1`, where C runs `athread`
/// (`tditeration.c:1433`) instead of `TreeDependentIteration`.
///
/// C's two implementations give different output on this input (6 lines
/// apart), so the two `.expected` files legitimately differ. `athread` with
/// one worker is deterministic (3/3 identical C runs) but walks the tree in
/// ascending order in EVERY cycle (`branchtable` is the identity unless
/// `randomseed != 0`, `:522`/`:728-730`), whereas the single-threaded loop
/// reverses odd cycles (`:1641-1648`). Reversing under `--thread 1` moved
/// one gap column in `s97` (`gggccc-gc-` vs `gggcccg-c-`); the fixed order
/// closes it. `athread` also converges once per cycle and stops on
/// `Converged2.` / `Oscillating?` rather than the single-threaded
/// per-branch rules — see `RefinementParams::per_cycle_convergence`.
///
/// The expected bytes come from the pinned `mafft-upstream` source built
/// in-tree on arm64 (`make -C mafft-upstream/core`), and were checked equal
/// to the reference 7.526 arm64 binary's output. `--maxiterate 2` is what
/// `--auto` picks here; C's output is the same for `--maxiterate 3` and
/// `1000` (the run converges in cycle 1).
#[test]
fn fftnsi_120seq_dna_matches_c_under_thread_1() {
    let f = fixtures();
    expect_identical(
        &["--thread", "1", "--retree", "2", "--maxiterate", "2"],
        &f.join("mtb_cds_120x1400.fa"),
        &f.join("mtb_cds_120x1400.thread1.expected"),
    );
}

/// The first 8 sequences of the same fixture, full FFT-NS-i (`--retree 2
/// --maxiterate 1000`), no `--thread`.
///
/// This pins the trailing-anchor flush of `searchAnchors` (`mltaln9.c`).
/// C's segment search runs `for( i=1; i<len-divWinSize; i++ )` and, when a
/// high-conservation region is still open at the end, closes it with
/// `seg->end = i` — the loop's EXIT value `len - divWinSize`, one past the
/// last iterated column. The port used the last iterated column, which
/// moves the final anchor one column left whenever `start + len` is even:
/// here C splits at 1412 and we split at 1411, so refinement segments 10
/// and 11 saw different columns and `s5` ended up with `aactca-cctgacc`
/// instead of C's `aactcacc-tgacc`. Segments 1–9 (and every branch score
/// in them) were already identical, and the result is the same with
/// `--nofft`, `--thread 1`, and any `--maxiterate >= 1`. Subsets of 5–7,
/// 9, 10 or 12 sequences do not reach the end of the alignment with an
/// open region, which is why the 120-sequence tests above did not catch
/// it.
///
/// Expected bytes: pinned `mafft-upstream` built in-tree on arm64, checked
/// equal to the reference 7.526 arm64 binary's output.
#[test]
fn fftnsi_first8_dna_trailing_anchor_matches_c() {
    let f = fixtures();
    expect_identical(
        &["--retree", "2", "--maxiterate", "1000"],
        &f.join("mtb_cds_first8.fa"),
        &f.join("mtb_cds_first8.fftnsi.expected"),
    );
}

/// Differential test on realistic clusters under the pipeline invocation
/// `--auto --adjustdirection --thread 1 --nuc`. Every cluster must be
/// byte-identical to C MAFFT 7.526.
#[test]
fn r2_real_ancestor_clusters_match_c_under_auto() {
    let dir = fixtures().join("r2_dna_clusters");
    let mut inputs: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("r2_dna_clusters")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .collect();
    inputs.sort();
    assert_eq!(inputs.len(), 30, "expected 30 clusters in {}", dir.display());
    let mut failures = Vec::new();
    for input in &inputs {
        let expected = input.with_extension("expected");
        let got = run(&["--auto", "--adjustdirection", "--thread", "1", "--nuc"], input);
        let want = std::fs::read(&expected).expect("read expected");
        if got != want {
            failures.push(input.file_name().unwrap().to_string_lossy().into_owned());
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} clusters differ from C MAFFT 7.526: {failures:?}",
        failures.len(),
        inputs.len()
    );
}

/// The five panaroo-rs clusters, sorted, with a sanity check on the count.
fn panaroo_cluster_inputs() -> Vec<PathBuf> {
    let dir = fixtures().join("panaroo_tiny_dna_clusters");
    let mut inputs: Vec<PathBuf> = std::fs::read_dir(&dir)
        .expect("panaroo_tiny_dna_clusters")
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|x| x == "fa"))
        .collect();
    inputs.sort();
    assert_eq!(inputs.len(), 5, "expected 5 clusters in {}", dir.display());
    inputs
}

/// Run every panaroo cluster through `run_from` with `flags` and compare
/// against the reference file named by `expected_ext`.
fn panaroo_clusters_match(flags: &[&str], expected_ext: &str) {
    let inputs = panaroo_cluster_inputs();
    let mut failures = Vec::new();
    for input in &inputs {
        let expected = input.with_extension(expected_ext);
        let got = run(flags, input);
        let want = std::fs::read(&expected).expect("read expected");
        if got != want {
            failures.push(input.file_name().unwrap().to_string_lossy().into_owned());
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} panaroo clusters differ from C MAFFT 7.526 \
         ({flags:?}, .{expected_ext}): {failures:?}",
        failures.len(),
        inputs.len()
    );
}

/// panaroo-rs invokes MAFFT once per gene family as
/// `mafft --auto --adjustdirection --thread 1 --nuc <cluster.fa>`. Every
/// cluster must be byte-identical to the in-tree C reference build.
#[test]
fn panaroo_tiny_clusters_match_c_under_panaroo_argv() {
    panaroo_clusters_match(
        &["--auto", "--adjustdirection", "--thread", "1", "--nuc"],
        "expected",
    );
}

/// Same clusters without `--thread 1`: C takes its single-threaded
/// `TreeDependentIteration` path instead of `athread`, whose convergence
/// semantics differ (see `fftnsi_120seq_dna_matches_c`). The two references
/// are byte-identical on these inputs today; pinning both keeps either
/// C path from drifting unnoticed.
#[test]
fn panaroo_tiny_clusters_match_c_without_thread() {
    panaroo_clusters_match(&["--auto", "--adjustdirection", "--nuc"], "nothread.expected");
}

/// The in-memory entry point on the panaroo workload: parse each cluster
/// with `read_fasta`, align it via `run_from_seqs` with the panaroo argv (no
/// input path), and require the (name, row) list to equal the C reference
/// parsed the same way. `read_fasta` keeps `-`, keeps the full header
/// (including the `_R_` prefix `--adjustdirection` adds) and folds DNA to
/// lowercase, which is already the case C writes, so the comparison is exact.
#[test]
fn panaroo_tiny_clusters_run_from_seqs_matches_c() {
    use mafft_rs::{run_from_seqs, SilentProgress};

    let argv: Vec<OsString> =
        ["mafft-rs", "--quiet", "--auto", "--adjustdirection", "--thread", "1", "--nuc"]
            .iter()
            .map(OsString::from)
            .collect();
    let mut failures = Vec::new();
    for input in &panaroo_cluster_inputs() {
        let name = input.file_name().unwrap().to_string_lossy().into_owned();
        let set = mafft_io::read_fasta(input).expect("read cluster");
        let msa = run_from_seqs(argv.clone(), &set, &SilentProgress)
            .unwrap_or_else(|e| panic!("run_from_seqs {name}: {}", e.message()));
        let got: Vec<(String, Vec<u8>)> = msa.names.into_iter().zip(msa.sequences).collect();
        let want: Vec<(String, Vec<u8>)> = mafft_io::read_fasta(input.with_extension("expected"))
            .expect("read expected")
            .sequences
            .into_iter()
            .map(|s| (s.name, s.data))
            .collect();
        if got != want {
            failures.push(name);
        }
    }
    assert!(
        failures.is_empty(),
        "{} of 5 panaroo clusters differ from C MAFFT 7.526 via run_from_seqs: {failures:?}",
        failures.len()
    );
}
