//! Reserved namespace for future raw FFI bindings to MAFFT's C source.
//!
//! **This crate is intentionally empty** at version `0.0.1`. It exists
//! only to claim the `mafft-sys` name on crates.io ahead of any real
//! FFI shim release — by `-sys` convention, that slot belongs to MAFFT's
//! native C bindings. Squatters would otherwise be free to take it.
//!
//! For high-level alignment, use [`mafft`](https://docs.rs/mafft) or
//! [`mafft-core`](https://docs.rs/mafft-core). For the standalone CLI
//! binary, install [`mafft-rs`](https://crates.io/crates/mafft-rs).
//!
//! When a real FFI shim is ready (with vendored C source so the crate
//! builds standalone from a published tarball, no submodule required),
//! this crate will be bumped to `0.1.0` and the bindings exposed here.

#![no_std]
