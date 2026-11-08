//! The `hat2` distance-file round trip of C MAFFT.

/// The distance C MAFFT reads back after writing `d` to a `hat2` file.
///
/// C writes every distance with `fprintf( hat2p, DFORMAT, mtx[i][j] )`,
/// `DFORMAT = "%#6.3f"` (`mltaln.h:55`, `io.c::WriteHat2_pointer`), and
/// `io.c::readhat2_pointer` reads each 6-character field back with `atof`
/// (`io.c::input_new`). The value a later stage sees (`dvtditr`'s
/// refinement tree, `addsingle`'s placement) is therefore the 3-decimal
/// text, rounded from the *exact* binary value of `d` with ties to even,
/// then rounded back to the nearest double.
///
/// `(d * 1000.0).round() / 1000.0` is not equivalent, for two reasons:
///
/// * the product is rounded before `round` sees it, and can land exactly on
///   a `.5`: `0.5505` is stored as `0.550499999999999989…`, C prints
///   `0.550`, but `0.5505 * 1000.0` is exactly `550.5`, which rounds to 551;
/// * `f64::round` breaks exact ties away from zero (`0.0625` -> `0.063`),
///   printf breaks them to even (`0.062`).
///
/// Rust's `{:.3}` rounds the exact value with ties to even, like printf, and
/// `str::parse::<f64>` is correctly rounded, like `atof`, so formatting and
/// parsing gives C's value bit for bit.
pub(crate) fn hat2_round(d: f64) -> f64 {
    format!("{d:.3}").parse().unwrap_or(d)
}

#[cfg(test)]
mod tests {
    use super::hat2_round;

    #[test]
    fn product_rounding_onto_a_tie_follows_printf() {
        // The FFT-NS-i refinement distance that split C (612 columns) from
        // rust-MAFFT (611) on a 30-sequence BAliBASE batch.
        let d = f64::from_bits(0x3fe1_9db2_2d0e_5604);
        assert_eq!(d * 1000.0, 550.5);
        assert_eq!(format!("{d:.3}"), "0.550");
        assert_eq!(hat2_round(d), 0.55);
    }

    #[test]
    fn exact_ties_round_to_even_like_printf() {
        for (d, c) in [
            (0.0625, 0.062),
            (0.1875, 0.188),
            (0.3125, 0.312),
            (0.4375, 0.438),
            (0.5625, 0.562),
            (0.8125, 0.812),
            (1.0625, 1.062),
        ] {
            assert_eq!(hat2_round(d), c, "{d}");
        }
    }

    #[test]
    fn ordinary_values() {
        for (d, c) in [(0.0, 0.0), (0.0004, 0.0), (0.0006, 0.001), (0.27549, 0.275), (9.9, 9.9)] {
            assert_eq!(hat2_round(d), c, "{d}");
        }
    }
}
