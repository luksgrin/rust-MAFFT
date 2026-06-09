/// Complex number type for FFT computation.
///
/// Re-exports `num_complex::Complex64` as the replacement for the C
/// `Fukusosuu` struct. This gives us full arithmetic, `exp`, `from_polar`,
/// `norm`, and everything else we'll need for the FFT port in Phase 3.
pub use num_complex::Complex64;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_arithmetic() {
        let a = Complex64::new(3.0, 4.0);
        let b = Complex64::new(1.0, -2.0);

        assert_eq!(a + b, Complex64::new(4.0, 2.0));
        assert_eq!(a - b, Complex64::new(2.0, 6.0));
        assert_eq!(a * b, Complex64::new(11.0, -2.0));
        assert!((a.norm_sqr() - 25.0).abs() < 1e-12);
    }
}
