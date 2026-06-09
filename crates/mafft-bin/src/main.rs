//! Thin wrapper around `mafft_rs::run()`.
//!
//! The CLI logic lives in `lib.rs` so that it can be reused by `pymafft`
//! (bundled into the wheel) without duplicating argument parsing.

fn main() {
    mafft_rs::run();
}
