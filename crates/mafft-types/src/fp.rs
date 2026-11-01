//! Floating-point contraction policy: mirror the reference C build of the
//! target you run on.
//!
//! # The problem
//!
//! Many MAFFT inner loops accumulate `acc = a * b + acc`. Whether that is one
//! rounding (a fused multiply-add) or two (separate multiply, then add)
//! changes the last ulp of the result, and in tie-break-sensitive DP that
//! last ulp changes the alignment. C MAFFT 7.526 leaves `-ffp-contract`
//! unset, so the compiler, not the source, chooses the rounding. The same
//! C source therefore produces both alignments.
//!
//! The **reference build** is the pinned `mafft-upstream` source, built
//! in-tree with the upstream Makefile's default flags
//! (`make -C mafft-upstream/core`), on the platform you run on. Nothing
//! else — in particular no downloaded binary — defines the reference. That
//! Makefile passes `-O3 -std=c99` and does not set `-ffp-contract`:
//!
//! - **Apple clang on arm64** contracts `a * b + c` into one rounding
//!   (`FP_CONTRACT` on by default). The reference `disttbfast`, `dvtditr` and `tbfast`
//!   binaries each contain roughly 1025 `fmadd`/`fmsub` instructions (~1140
//!   counting the vector `fmla`/`fmls` forms, as CI's census reports).
//! - **gcc in that ISO C mode** does not contract. `-std=c99` implies
//!   `-ffp-contract=off`, so gcc emits no fused ops even on aarch64, where
//!   the CPU has FMA hardware. A gcc 15 build checked previously emitted 0
//!   fused instructions and matched the x86-64 fixtures. The Linux x86-64
//!   reference is the same: 0 `vfmadd`/`vfmsub`. **MSVC** does not contract
//!   by default either. **Apple clang on x86-64** does not contract
//!   either: without `-mfma` the baseline x86-64 target has no FMA
//!   instruction, so `-O3 -std=c99` emits `mulsd` + `addsd`.
//!
//! The visible consequence on the 36-sequence protein sample
//! (`mafft-upstream/test/sample`) with `--bl 50 --retree 2 --maxiterate 0`:
//! the Apple-clang C binary produces an alignment of width **712**; a gcc
//! or MSVC C binary produces width **738**. Both are "C MAFFT 7.526".
//!
//! Rust's [`f64::mul_add`] is *always* fused. Plain `a * b + c` is never
//! contracted by rustc. The port therefore picks, per target, which
//! arithmetic to emit, and mirrors the C toolchain usual for that platform:
//! Apple clang on macOS (fused on Apple Silicon only), gcc on Linux (including aarch64; distro and
//! bioconda packages are gcc builds), MSVC on Windows.
//!
//! # The policy
//!
//! Every site that Apple clang on arm64 contracts is written as
//! [`fmadd(a, b, c)`](fmadd). Under [`CONTRACTS_FMA`] it lowers to
//! `a.mul_add(b, c)`; otherwise to `a * b + c`. Sites that Apple clang does
//! *not* contract (e.g. `rootnode[s] += len * eff[s]` in
//! `mafft-tree::weighting::sequence_weights`) stay as plain `*` and `+` on
//! every target and must not be routed through this module.
//!
//! [`CONTRACTS_FMA`] is a compile-time constant, so the inner loop does not
//! grow a per-cell check. There is no runtime flag. It is selected by cargo
//! features on `mafft-types`, which every workspace crate forwards under
//! the same names:
//!
//! - neither: fuse on Apple Silicon (`target_vendor = "apple"` and
//!   `target_arch = "aarch64"`), and not on any other target, Intel macOS
//!   included. That is what `cargo build` and the published
//!   wheels do.
//! - `fp-contract-fma`: always fuse, on any target (match an Apple-clang C
//!   build).
//! - `fp-contract-none`: never fuse, on any target (match a gcc or MSVC C
//!   build).
//! - both: compile error.
//!
//! Someone who installed a wheel and wants the other policy rebuilds from
//! source with the feature.
//!
//! Fixtures whose bytes depend on the policy are committed twice, as
//! `<name>.fma` (in-tree Apple-clang build) and `<name>.nofma` (in-tree gcc
//! build, including Linux aarch64), and the test helpers pick the variant
//! matching [`CONTRACTS_FMA`]. `scripts/regen_policy_fixtures.sh` regenerates
//! and checks them from the in-tree build.

#[cfg(all(feature = "fp-contract-fma", feature = "fp-contract-none"))]
compile_error!(
    "features `fp-contract-fma` and `fp-contract-none` are mutually exclusive; \
     enable at most one (default: fused on Apple Silicon, not fused elsewhere)"
);

/// `true` when this build fuses `a * b + c` into a single-rounding FMA at
/// every [`fmadd`] site, `false` when it rounds twice.
///
/// Resolution order: `fp-contract-fma` => `true`; `fp-contract-none` =>
/// `false`; neither => `cfg!(all(target_vendor = "apple", target_arch = "aarch64"))`. See the
/// [module docs](self) for why. This is a compile-time branch: there is
/// no runtime flag.
pub const CONTRACTS_FMA: bool = cfg!(feature = "fp-contract-fma")
    || (!cfg!(feature = "fp-contract-none")
        && cfg!(all(target_vendor = "apple", target_arch = "aarch64")));

/// Policy-dependent multiply-add: `a * b + c`.
///
/// One rounding (`a.mul_add(b, c)`) when [`CONTRACTS_FMA`] is set, two
/// roundings (`a * b + c`) otherwise. Operand order matters for the
/// two-rounding form only through the usual IEEE rules (multiplication and
/// addition commute), but callers preserve the C source's operand order and
/// nesting anyway so the fused form reproduces clang's exact `fmadd`
/// sequence, e.g. `fmadd(a, b, fmadd(b, c, a * c))` for
/// `s = b*c + c*a + a*b` in `calcW`.
#[inline(always)]
pub fn fmadd(a: f64, b: f64, c: f64) -> f64 {
    if CONTRACTS_FMA {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmadd_matches_selected_policy() {
        // 0.1 * 0.2 + 0.3 rounds differently fused vs unfused; whichever the
        // policy is, `fmadd` must agree with the corresponding primitive.
        let (a, b, c) = (0.1f64, 0.2f64, 0.3f64);
        let expected = if CONTRACTS_FMA {
            a.mul_add(b, c)
        } else {
            a * b + c
        };
        assert_eq!(fmadd(a, b, c).to_bits(), expected.to_bits());
    }

    #[test]
    fn fmadd_is_exact_when_no_rounding_needed() {
        assert_eq!(fmadd(2.0, 3.0, 4.0), 10.0);
        assert_eq!(fmadd(-1.5, 2.0, 1.0), -2.0);
    }
}
