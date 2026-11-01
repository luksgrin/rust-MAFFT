# Upstream policy

rust-MAFFT is a port of a tool that is still maintained upstream
([GitLab `sysimm/mafft`](https://gitlab.com/sysimm/mafft)). This page
states how the port follows it: which upstream version is the
reference, what happens when upstream moves, which way bug fixes flow,
and what the port will not do.

## Scope

rust-MAFFT tracks MAFFT. It does not add alignment behaviour that MAFFT
does not have: every mode, flag, default, scoring matrix and gap
penalty it implements is the one the pinned upstream source implements,
with the same meaning.

The correctness criterion is byte-identity with the pinned upstream
version, built in-tree on the same platform (see
[Byte-identity](byte-identity.md#the-reference-build)). Because the
reference is a specific commit and the test suite compares against it
at the function level (FFI cross-validation, 82 test functions across
15 `cross_validate*.rs` files) and the output level (C-canonical
fixtures, the `--bl 50` sentinel, `regen_policy_fixtures.sh check`),
the port cannot drift in output without a test failing.

Engineering choices that do not change output bytes — no global state,
Rayon parallelism, flattened DP scratch buffers, the Rust and Python
APIs — are not behaviour and are outside this policy. So are the
places where C MAFFT itself has no single answer to match:
`--thread N ≥ 2` (nondeterministic in C; rust-MAFFT is deterministic and
reproduces C's `--thread 1` path) and the `qsort` tie order behind
`--parttree --reorder` on duplicate sequences (platform-dependent in C;
rust-MAFFT ships the BSD algorithm and matches macOS C MAFFT 7.526
everywhere).

## Version pinning

The reference is the `mafft-upstream/` git submodule
(`url = https://gitlab.com/sysimm/mafft.git` in `.gitmodules`). Its
gitlink is the single source of truth: `mafft-c-bindings/build.rs`
compiles the FFI shim from it, `make -C mafft-upstream/core` builds the
C binaries CI compares against from it, and every C-canonical fixture
is generated from it.

The current pin is commit `0a2319b` (`git submodule status` reports
`v7.526-2-g0a2319b`): the upstream HEAD of branch `main` as of
2026-09-06, two commits past the `v7.526` tag (`ee97999`). Both extra
commits are the Alpine 3.23 build fix and carry no algorithm change, so
the version is reported as "MAFFT 7.526" throughout.

The pin is restated, and must be kept consistent with the gitlink, in:

- `README.md`, section "Upstream MAFFT" (commit, tag, date checked);
- the "7.526" version string in `README.md`, `docs/index.md`,
  `mkdocs.yml` (`site_description`), the `abstract` in `CITATION.cff`
  and the test-fixture READMEs;
- the `CHANGELOG.md` section of every release, which names the upstream
  commit and tag the release was validated against. A release is never
  cut with a pin that differs from the one its notes name.

## Upgrade procedure

The [Check MAFFT upstream](https://github.com/luksgrin/rust-MAFFT/blob/main/.github/workflows/check-mafft-upstream.yml)
workflow runs weekly (Monday 06:00 UTC) and on manual dispatch. When
the pinned commit is behind upstream HEAD, or upstream has a version
tag the pin does not contain, it opens an issue titled
`upstream: mafft-upstream pin is behind …` (or updates the open one)
naming the pinned and upstream commits and tags. It never moves the
submodule. That issue is the trigger for the procedure below.

1. **Read the upstream diff.**
   `git -C mafft-upstream log --stat <pinned>..<new>`. Classify each
   commit as build/packaging only, or touching `core/` sources, the
   `scripts/mafft` wrapper, or the matrices.
2. **Bump the submodule** on a branch:
   ```bash
   git -C mafft-upstream fetch --tags origin
   git -C mafft-upstream checkout <new-commit-or-tag>
   git add mafft-upstream
   ```
3. **Rebuild the C reference in-tree** and record the contraction
   census, which must stay in the documented range (~1140 fused ops per
   engine binary on arm64 clang, 0 on x86-64 gcc):
   ```bash
   make -C mafft-upstream/core
   scripts/fma_census.sh mafft-upstream/binaries
   ```
4. **Run the FFI cross-validation and the full suite.** The FFI shim is
   recompiled from the new submodule, so these compare the Rust engine
   against the new C directly:
   ```bash
   cargo test --workspace --exclude pymafft --release
   ```
   If a changed C function is one of those copied into
   `crates/mafft-c-bindings/wrappers/` (copies of `static` or
   `main()`-adjacent code from `Salignmm.c`, `MSalignmm.c`, `mltaln9.c`
   and `addsingle.c`, so the FFI tests can link it),
   re-sync the copy first, and check that `build.rs`'s source list still
   matches `core/`.
5. **Re-derive the C-canonical fixtures.** The fixture-based tests
   compare against committed outputs, so they keep passing against the
   *old* C until the fixtures are regenerated. Run the policy-sensitive
   ones through the script, on an arm64 host for `.fma` and an x86-64
   host for `.nofma` (or let the CI legs run `check`):
   ```bash
   scripts/regen_policy_fixtures.sh mafft-upstream/scripts/mafft check
   ```
   and regenerate the rest with the commands in
   `crates/mafft-core/tests/fixtures/README.md` and the per-corpus
   READMEs under `crates/mafft-bin/tests/fixtures/`. Any byte difference
   from the committed fixture is an upstream output change.
6. **Run the BAliBASE sweep** with `scripts/balibase_parity.py` over
   RV11–RV50 (386 protein test sets) for the five modes in the README
   table (FFT-NS-i, L-INS-i, G-INS-i, E-INS-i, G-INS-i + `--allowshift`,
   all at `--maxiterate 1000`), with `--mafft` pointing at the new
   in-tree build. The bar is 1930/1930 again.
7. **Triage every divergence** found in steps 4–6:
   - *Upstream algorithm change or upstream bug fix* — port it. The new
     upstream behaviour becomes the reference, including when it changes
     output; the fixtures move with it.
   - *Upstream regression* (a change that looks unintended) — report it
     on the upstream tracker. Until upstream resolves it, either stay on
     the previous pin, or bump and mirror the regression; whichever is
     chosen, rust-MAFFT keeps matching its pinned C. It does not ship a
     pin it does not match.
8. **Update the restatements** listed under
   [Version pinning](#version-pinning), and add a CHANGELOG entry that
   names the old and new commits and tags and lists any output change.
   An upgrade that changes output for some command line is called out
   the way 0.2.0's "Output changes" section is.

Expected turnaround: an upstream change that is build-only, or that
leaves every test, fixture and the BAliBASE sweep unchanged, is bumped
within two weeks of the drift issue being opened. A change to the
alignment code takes as long as porting it takes; until it lands, the
previous pin remains the documented reference and releases keep naming
it.

## Bug-fix direction

Fixes flow upstream first. If work on the port finds a defect that is
also in the C code (an out-of-bounds read, an uninitialised value, a
wrong result on some input), it is reported on the
[upstream tracker](https://gitlab.com/sysimm/mafft/-/issues) with a
reproducer against the C build. The port keeps reproducing the C
behaviour until upstream changes it; the fix then arrives through the
[upgrade procedure](#upgrade-procedure) like any other upstream change.
The one exception is a defect that would make the Rust code unsound or
panic (an out-of-bounds index C gets away with): the Rust side handles
it without changing the output C produces on inputs where C is well
defined.

rust-MAFFT does not knowingly carry behaviour differences from its
pinned upstream. The known divergences are listed in the README's
[Known limitations](https://github.com/luksgrin/rust-MAFFT#known-limitations);
in summary:

- **Two tied-trace residuals**
  ([B](https://github.com/luksgrin/rust-MAFFT#b-intentional-residuals-tied-trace)).
  `--exp ≥ 4.30` on FFT-NS-2 without refinement (same width, 517, on
  the 36-seq sample, 513 of 517 columns identical, `--scoreout` 55049.75
  for C vs 55081.77 for rust-MAFFT; closes with `--nofft` or any
  `--maxiterate ≥ 1`), and `--add` with an added sequence whose last
  residue is zero-scoring (C writes `…KWRR-----X`, rust-MAFFT
  `…KWRRX-----`). Both are DP tie-breaks on zero-scoring cells at the
  same width; the shifted residue changes the sum-of-pairs score by a
  fraction of a percent, so they are not score-equivalent. They
  are carried because the first is below per-step FFI replay resolution
  (each `fft_profile_align` step is byte-identical to C's `Falign`; the
  shift comes from accumulated FP order the replay cannot reproduce) and
  neither is reachable from a realistic command line. They are
  residuals, not features: a fix that closes either without changing
  any other output is taken.
- **External-dependency-blocked modes**
  ([A](https://github.com/luksgrin/rust-MAFFT#a-external-dependency-blocked-modes)).
  `--xinsi`, `--qinsi`, `--scarnalike` are wired as in C and validated
  only where the third-party binary is present.
- **Out-of-scope modes**
  ([C](https://github.com/luksgrin/rust-MAFFT#c-out-of-scope-by-design)).
  RNA-structure pipelines, external pairwise aligners and `--mpi` are
  not implemented; `--pdbidlist` / `--pdbfilelist` reproduce upstream's
  own "temporarily unavailable, 2018/Dec." exit. Not implementing an
  upstream mode is a gap, not a different behaviour: it never makes a
  supported command line produce different output.
- **Other known differences**
  ([E](https://github.com/luksgrin/rust-MAFFT#e-other-known-differences)).
  `--add` runs one guide-tree cycle where C runs two, `--add --localpair`
  is not implemented, `--groupsize` is ignored, `--partsize` is not
  accepted, and `--memsave` does not change the FFT merges. A few inputs
  outside the parity corpora (3000 short proteins, a 120 kb sequence,
  low-complexity and homopolymer DNA, `nan` / `inf` `--treein` lengths)
  give different bytes. These are gaps in the port, carried until fixed.
  Keeping long names whole where C truncates them is intentional.

The floating-point contraction policy is not a residual. C MAFFT 7.526
gives different bytes on arm64 clang and x86-64 gcc; rust-MAFFT mirrors
the same-platform C build by default, and the `fp-contract-fma` /
`fp-contract-none` features reproduce the other platform's build. Both
outputs are upstream behaviour.

## What would make it a fork

rust-MAFFT would be a fork, rather than a port, if its output for some
command line stopped being defined by the pinned upstream source. Any of
the following would do that, and the project commits to none of them:

- **New default behaviour.** No change to what any existing command line
  produces, other than moving toward the pinned C output or following an
  upstream change through the [upgrade procedure](#upgrade-procedure).
- **Changed flag semantics.** Every flag accepted by both tools means
  what it means in the pinned `scripts/mafft` and the C engines,
  including defaults, `argv[0]` shortcuts (`linsi`, `ginsi`, …) and
  `--auto`'s size thresholds. Controls that exist only in rust-MAFFT are
  limited to ones that do not change which alignment is computed:
  output plumbing (`-o/--output`, `--format`, `--linewidth`,
  `--namelength`), `--c-compat` (replicates C's static-TLS
  `reuseprofiles` memoisation, i.e. moves *toward* C), the
  `fp-contract-fma` / `fp-contract-none` cargo features (select between
  the two upstream builds) and diagnostic environment variables such as
  `MAFFT_RS_REFINE_STATS`.
- **Scoring changes outside upstream.** No substitution matrix, gap
  penalty, distance or weighting change that upstream has not made.
- **Algorithm changes outside upstream.** No "improvement" to an
  algorithm in a mode that shares a name with an upstream mode.

[Byte-identity](byte-identity.md#what-this-doesnt-mean) notes that
algorithmic work is still possible behind its own modes or flags, as
long as the existing modes stay byte-identical. Under this policy such
work is proposed upstream first; it is not added to rust-MAFFT's
alignment engine as rust-only behaviour.

## Relationship with upstream and attribution

MAFFT is the work of Kazutaka Katoh and colleagues; every algorithmic
decision in rust-MAFFT traces back to it. The project states this in
the README ("Acknowledgments and Citation"), on the
[Citation](../citation.md) page and in `CITATION.cff`, and requires
anyone who uses rust-MAFFT in published work to cite the original MAFFT
paper (Katoh & Standley 2013, doi:10.1093/molbev/mst010), plus the
mode-specific references listed on the Citation page.

The ported algorithmic constructs remain under upstream's BSD-3-Clause
licence (copyright Kazutaka Katoh; `LICENSE-BSD`), alongside MIT for
the original Rust code, and the submodule carries upstream's own
licence unmodified.

Defects found in the C code and questions about upstream behaviour go
to the [upstream GitLab tracker](https://gitlab.com/sysimm/mafft/-/issues).
Issues that are specific to the port go to
[github.com/luksgrin/rust-MAFFT](https://github.com/luksgrin/rust-MAFFT/issues),
not upstream.
