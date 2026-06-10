# Changelog

All notable changes to rust-MAFFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
from 1.0.0 onwards. Pre-1.0 releases may contain breaking changes between
minor versions; the `mafft` and `mafft-rs` crates aim to keep their
public surfaces stable from 0.1.0 anyway.

## [Unreleased]

## [0.1.2] - 2026-06-10

Metadata fixes. No engine, library, or CLI behaviour changes from
0.1.1 / 0.1.0; bumped so that all distribution channels can publish
together cleanly.

### Fixed

- `CITATION.cff`: dropped the SPDX expression `MIT AND BSD-3-Clause`
  in favour of the single SPDX identifier `MIT` (Zenodo's `cffconvert`
  pipeline rejected the expression form with "Citation metadata load
  failed"). The BSD-3-Clause attribution for the algorithmic constructs
  ported from upstream MAFFT remains in `LICENSE-BSD` and the workspace
  `Cargo.toml` `license = "MIT AND BSD-3-Clause"` (where cargo accepts
  expressions fine).
- Maintainer email updated to `lucas.goiriz@csic.es` in workspace
  authors and `pymafft` pyproject metadata.

### Notes

- v0.1.1 published partially: docs + GH-release binaries succeeded;
  crates.io got 5 of 10 crates before hitting the new-crate rate
  limit; PyPI publish was gated by a transient quay.io docker-pull
  flake on one Linux wheel. v0.1.2 retries all channels with a valid
  CITATION.cff.

## [0.1.1] - 2026-06-10

Release-pipeline fixes. No engine, library, or CLI behaviour changes
from `0.1.0`; bumped only because `pymafft 0.1.0` was published to
PyPI before the rest of the release-day workflows could be fixed.

### Fixed

- `release.yml`: added `permissions: contents: write` so
  `softprops/action-gh-release@v2` can attach the cross-compiled
  binaries to the GitHub release (previously failed with
  "Resource not accessible by integration").
- `docs.yml`: the `github-pages` deployment environment is now scoped
  to allow tag-triggered deploys (`v*`) in addition to `main` pushes.
- `python.yml`: x86_64-apple-darwin wheel is now cross-compiled from
  the Apple Silicon `macos-latest` runner (the previously-targeted
  `macos-13` Intel runner pool was retired by GitHub in 2026).
- `python.yml`: Linux wheel build no longer changes the cwd inside the
  manylinux container (`cd ../..` had broken maturin's manifest
  resolution; replaced with `--manifest-path ../../Cargo.toml`).

### Notes

- `pymafft 0.1.0` was yanked on PyPI after this release. Use 0.1.1+
  for any new installs. Programmatic and CLI behaviour is unchanged.

## [0.1.0] - 2026-06-10 [YANKED]

Initial release attempt. `pymafft` published to PyPI successfully but
the crates.io publish (email verification not yet completed), GitHub
release binary upload (missing `contents: write` permission), and
docs deploy (environment-rule rejected the tag) all failed. Superseded
by 0.1.1.

### Engine

- **Byte-identical to C MAFFT 7.526** across the full BAliBASE 3
  fixture set (1930/1930) for every supported mode: FFT-NS-2, FFT-NS-i,
  G-INS-i, L-INS-i, E-INS-i, PartTree, DPPartTree, plus `--add`,
  `--adjustdirection` (k-mer mode), `--oneiteration`, `--bestfirst`,
  `--allowshift`, `--unalignlevel`, `--skipiterate`, `--pileup`,
  `--youngestlinkage` / `--averagelinkage` / `--minimumlinkage` /
  `--mixedlinkage`.
- **Stub-parity** for `--pdbidlist` and `--pdbfilelist`: print the
  upstream "temporarily unavailable, 2018/Dec." message and exit 0 to
  match C MAFFT 7.526 behaviour verbatim (these were disabled
  upstream).
- **Pure Rust runtime**: the release binary compiles zero C code.
  `mafft-c-bindings` (the FFI cross-validation shim) is a dev-only
  internal crate.

### Distribution

Four parallel release channels off a single GitHub release tag:

- **Rust library**: `cargo add mafft` (top-level re-export of
  `mafft-core`, `mafft-types`, `mafft-io`). Sub-crates also published
  individually for finer-grained dependence.
- **Rust CLI**: `cargo install mafft-rs`.
- **Python package**: `pip install pymafft` ships 25 wheels (Linux
  x86_64 + aarch64, macOS Intel + Apple Silicon, Windows x86_64 ×
  Python 3.9–3.13) plus an sdist.
- **Python executable**: the `pymafft` wheel bundles the `mafft-rs`
  binary inside it; `pip install pymafft` puts `mafft-rs` on `$PATH`
  via a console-script entry. No Rust toolchain required.
- **Pre-built CLI binary**: attached to every GitHub release for 5
  targets (Linux x86_64 + aarch64, macOS Intel + Apple Silicon,
  Windows x86_64).

### Testing

- ~430 Rust tests including ~140 FFI cross-validation tests that
  compile MAFFT C source in-tree and compare per-function outputs
  byte-for-byte.
- 76 Python tests across binding parity, biopython interop, and the
  bundled-CLI smoke path.
- Full CI matrix on every push / PR: Linux byte-identity gate, plus
  macOS + Windows build + lib-test cross-platform sanity.

### Documentation

- Site at https://luksgrin.github.io/rust-MAFFT (Material for MkDocs)
- Per-crate `cargo doc` published to https://docs.rs

[Unreleased]: https://github.com/luksgrin/rust-MAFFT/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.2
[0.1.1]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.1
[0.1.0]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.0
