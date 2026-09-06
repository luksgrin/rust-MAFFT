//! Floating-point contraction policy: mirror the reference C build of the
//! target you run on.
//!
//! # The problem
//!
//! Many MAFFT inner loops accumulate `acc = a * b + acc`. Whether that is one
//! rounding (a fused multiply-add) or two (separate multiply, then add)
//! changes the last ulp of the result, and in tie-break-sensitive DP that
//! last ulp changes the alignment. C MAFFT 7.526 leaves the choice to the
//! compiler, so **C MAFFT itself is not bit-reproducible across CPU
//! architectures**:
//!
//! The **reference build** is the pinned `mafft-upstream` source, built
//! in-tree with the upstream Makefile's default flags
//! (`make -C mafft-upstream/core`), on the platform you run on. Nothing
//! else — in particular no downloaded binary — defines the reference. That
//! build behaves differently per platform:
//!
//! - **arm64 macOS** (Apple clang, `-O3`, `FP_CONTRACT=on` by default): the
//!   reference `disttbfast`, `dvtditr` and `tbfast` binaries each contain
//!   roughly 1025 `fmadd`/`fmsub` instructions (~1140 counting the vector
//!   `fmla`/`fmls` forms, as CI's census reports). clang contracts
//!   `a * b + c` into FMA wherever the source expression allows it.
//! - **baseline x86-64** (gcc `-O3`, the same pinned source and default
//!   flags on Linux): zero `vfmadd`/`vfmsub` instructions, against ~1250
//!   `mulsd` and ~1550 `addsd`. Baseline x86-64 has no FMA unit, so gcc
//!   cannot contract.
//!
//! The visible consequence on the 36-sequence protein sample
//! (`mafft-upstream/test/sample`) with `--bl 50 --retree 2 --maxiterate 0`:
//! the arm64 C binary produces an alignment of width **712**, the x86-64 C
//! binary produces width **738**. Both are "C MAFFT 7.526".
//!
//! Rust's [`f64::mul_add`] is *always* fused: a native `fmadd` on aarch64,
//! but a slow software `fma()` on baseline x86-64 (no hardware FMA to lower
//! to). Plain `a * b + c` is never contracted by rustc. So the port has to
//! pick, per target, which arithmetic to emit — and the only sensible rule
//! is to mirror the reference C build of the target it runs on, because that
//! is the binary users compare against.
//!
//! # The policy
//!
//! Every site that C's arm64 clang build contracts is written as
//! [`fmadd(a, b, c)`](fmadd). Under [`CONTRACTS_FMA`] it lowers to
//! `a.mul_add(b, c)`; otherwise to `a * b + c`. Sites that arm64 clang does
//! *not* contract (e.g. `rootnode[s] += len * eff[s]` in
//! `mafft-tree::weighting::sequence_weights`) stay as plain `*` and `+` on
//! every target and must not be routed through this module.
//!
//! [`CONTRACTS_FMA`] is selected by cargo features on `mafft-types`, which
//! every workspace crate forwards under the same names:
//!
//! - `fp-contract-fma`: always fuse (matches the arm64 clang reference).
//! - `fp-contract-none`: never fuse (matches the baseline x86-64 gcc
//!   reference).
//! - neither: choose by target — fused on `aarch64`, not fused elsewhere.
//! - both: compile error.
//!
//! Fixtures whose bytes depend on the policy are committed twice, as
//! `<name>.fma` (in-tree C build on arm64) and `<name>.nofma` (in-tree C
//! build on x86-64), and the test helpers pick the variant matching
//! [`CONTRACTS_FMA`]. `scripts/regen_policy_fixtures.sh` regenerates and
//! checks them from the in-tree build; CI does so on both platforms.


#[cfg(all(feature = "fp-contract-fma", feature = "fp-contract-none"))]
compile_error!(
    "features `fp-contract-fma` and `fp-contract-none` are mutually exclusive; \
     enable at most one (default: fused on aarch64, not fused elsewhere)"
);

/// `true` when this build fuses `a * b + c` into a single-rounding FMA at
/// every [`fmadd`] site, `false` when it rounds twice.
///
/// Resolution order: `fp-contract-fma` => `true`; `fp-contract-none` =>
/// `false`; neither => `cfg!(target_arch = "aarch64")`. See the
/// [module docs](self) for why.
pub const CONTRACTS_FMA: bool = cfg!(feature = "fp-contract-fma")
    || (!cfg!(feature = "fp-contract-none") && cfg!(target_arch = "aarch64"));

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
        let expected = if CONTRACTS_FMA { a.mul_add(b, c) } else { a * b + c };
        assert_eq!(fmadd(a, b, c).to_bits(), expected.to_bits());
    }

    #[test]
    fn fmadd_is_exact_when_no_rounding_needed() {
        assert_eq!(fmadd(2.0, 3.0, 4.0), 10.0);
        assert_eq!(fmadd(-1.5, 2.0, 1.0), -2.0);
    }
}
