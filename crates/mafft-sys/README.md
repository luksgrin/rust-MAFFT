# mafft-sys

**Reserved namespace.** This crate is intentionally empty at `0.0.1` —
it exists to claim the `mafft-sys` name on crates.io ahead of any real
FFI shim release. By `-sys` convention, that slot is for MAFFT's native
C bindings.

If you're looking for the alignment engine, use:

- [`mafft`](https://crates.io/crates/mafft) — ergonomic top-level crate (re-exports `mafft-core`)
- [`mafft-core`](https://crates.io/crates/mafft-core) — direct dependency for finer control
- [`mafft-rs`](https://crates.io/crates/mafft-rs) — install the standalone CLI: `cargo install mafft-rs`

A real FFI shim with vendored C source will land at `0.1.0`.
