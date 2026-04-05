/// Complex number type for FFT computation.
///
/// Re-exports `num_complex::Complex64` as the replacement for the C
/// `Fukusosuu` struct. This gives us full arithmetic, `exp`, `from_polar`,
/// `norm`, and everything else we'll need for the FFT port in Phase 3.
pub use num_complex::Complex64;

// ---------------------------------------------------------------------------
// FFI conversion: Fukusosuu <-> Complex64
// ---------------------------------------------------------------------------

/// Convert a C `Fukusosuu` to a Rust `Complex64`.
pub fn from_fukusosuu(c: mafft_sys::Fukusosuu) -> Complex64 {
    Complex64::new(c.R, c.I)
}

/// Convert a Rust `Complex64` to a C `Fukusosuu`.
pub fn to_fukusosuu(c: Complex64) -> mafft_sys::Fukusosuu {
    mafft_sys::Fukusosuu { R: c.re, I: c.im }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_arithmetic() {
        let a = Complex64::new(3.0, 4.0);
        let b = Complex64::new(1.0, -2.0);

        assert_eq!(a + b, Complex64::new(4.0, 2.0));
        assert_eq!(a - b, Complex64::new(2.0, 6.0));
        // (3+4i)(1-2i) = 11 - 2i
        assert_eq!(a * b, Complex64::new(11.0, -2.0));
        assert!((a.norm_sqr() - 25.0).abs() < 1e-12);
    }

    #[test]
    fn ffi_roundtrip() {
        let rust = Complex64::new(1.5, -2.5);
        let c = to_fukusosuu(rust);
        let back = from_fukusosuu(c);
        assert_eq!(rust, back);
    }
}
