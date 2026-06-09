# Changelog

All notable changes to rust-MAFFT will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
from 1.0.0 onwards. Pre-1.0 releases may contain breaking changes between
minor versions; the `mafft` and `mafft-rs` crates aim to keep their
public surfaces stable from 0.1.0 anyway.

## [Unreleased]

## [0.1.0] - TBD

Initial release.

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

[Unreleased]: https://github.com/luksgrin/rust-MAFFT/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/luksgrin/rust-MAFFT/releases/tag/v0.1.0
