use std::io::{BufRead, BufReader, Read, Write};

use crate::error::IoError;

/// MAFFT's hat2 distance matrix format.
///
/// Stores an upper-triangular pairwise distance matrix with sequence names.
/// The on-disk format is:
/// ```text
/// 1
///    <nseq>
///  <scaled_max>
///    1. name1
///    2. name2
///    ...
/// <d(0,1)> <d(0,2)> ... (12 per line, 6-char fixed-width fields)
/// <d(1,2)> <d(1,3)> ...
/// ```
#[derive(Debug, Clone)]
pub struct Hat2Matrix {
    /// Sequence names.
    pub names: Vec<String>,
    /// Upper-triangular distances: `distances[i][j]` is the distance
    /// between sequence `i` and sequence `i + j + 1`.
    /// So row `i` has `nseq - i - 1` elements.
    pub distances: Vec<Vec<f64>>,
}

impl Hat2Matrix {
    pub fn nseq(&self) -> usize {
        self.names.len()
    }

    /// Get distance between sequences i and j (i != j).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i < j {
            self.distances[i][j - i - 1]
        } else if i > j {
            self.distances[j][i - j - 1]
        } else {
            0.0
        }
    }
}

/// Values per line in hat2 output.
const VALS_PER_LINE: usize = 12;

/// Read a hat2 distance matrix.
pub fn read_hat2<R: Read>(reader: R) -> Result<Hat2Matrix, IoError> {
    let reader = BufReader::new(reader);
    let mut lines = reader.lines();

    // Line 1: format identifier (skip)
    lines.next().ok_or(IoError::Hat2Format("empty file".into()))??;

    // Line 2: nseq
    let nseq_line = lines
        .next()
        .ok_or(IoError::Hat2Format("missing nseq line".into()))??;
    let nseq: usize = nseq_line
        .trim()
        .parse()
        .map_err(|_| IoError::Hat2Format(format!("invalid nseq: '{nseq_line}'")))?;

    // Line 3: scaled max (informational, we don't need it)
    lines.next().ok_or(IoError::Hat2Format("missing max line".into()))??;

    // Lines 4..4+nseq: "   N. name"
    let mut names = Vec::with_capacity(nseq);
    for _ in 0..nseq {
        let line = lines
            .next()
            .ok_or(IoError::Hat2Format("truncated name section".into()))??;
        // Format: "   1. sequence_name" -- extract everything after ". "
        let name = if let Some(pos) = line.find(". ") {
            line[pos + 2..].to_string()
        } else {
            line.trim().to_string()
        };
        names.push(name);
    }

    // Distance values: upper triangle, 6-char fixed-width fields.
    // Collect all remaining non-empty content into a flat token stream.
    let mut all_values: Vec<f64> = Vec::new();
    for line_result in lines {
        let line = line_result?;
        // Parse 6-char fixed-width fields, or fall back to whitespace splitting
        if !line.trim().is_empty() {
            for token in line.split_whitespace() {
                let val: f64 = token
                    .parse()
                    .map_err(|_| IoError::Hat2Format(format!("invalid float: '{token}'")))?;
                all_values.push(val);
            }
        }
    }

    // Distribute into upper-triangular rows
    let expected = nseq * (nseq - 1) / 2;
    if all_values.len() != expected {
        return Err(IoError::Hat2Format(format!(
            "expected {} distance values, got {}",
            expected,
            all_values.len()
        )));
    }

    let mut distances = Vec::with_capacity(nseq);
    let mut offset = 0;
    for i in 0..nseq {
        let row_len = nseq - i - 1;
        distances.push(all_values[offset..offset + row_len].to_vec());
        offset += row_len;
    }

    Ok(Hat2Matrix { names, distances })
}

/// Write a hat2 distance matrix.
pub fn write_hat2<W: Write>(mat: &Hat2Matrix, writer: &mut W) -> Result<(), IoError> {
    let nseq = mat.nseq();

    // Find max distance for the header
    let max_dist = mat
        .distances
        .iter()
        .flat_map(|row| row.iter())
        .cloned()
        .fold(0.0_f64, f64::max);

    // Header
    writeln!(writer, "    1")?;
    writeln!(writer, "{nseq:5}")?;
    writeln!(writer, " {:.3}", max_dist * 2.5)?;

    // Names
    for (i, name) in mat.names.iter().enumerate() {
        writeln!(writer, "{:>4}. {name}", i + 1)?;
    }

    // Distance values: upper triangle, rows 0..nseq-1
    for i in 0..nseq.saturating_sub(1) {
        let row = &mat.distances[i];
        for (col, val) in row.iter().enumerate() {
            write!(writer, "{val:6.3}")?;
            if (col + 1) % VALS_PER_LINE == 0 || col == row.len() - 1 {
                writeln!(writer)?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn sample_matrix() -> Hat2Matrix {
        Hat2Matrix {
            names: vec!["seq1".into(), "seq2".into(), "seq3".into()],
            distances: vec![
                vec![0.123, 0.456], // (0,1), (0,2)
                vec![0.789],        // (1,2)
            ],
        }
    }

    #[test]
    fn roundtrip_hat2() {
        let mat = sample_matrix();

        let mut buf = Vec::new();
        write_hat2(&mat, &mut buf).unwrap();

        let parsed = read_hat2(Cursor::new(&buf)).unwrap();
        assert_eq!(parsed.nseq(), 3);
        assert_eq!(parsed.names, mat.names);
        assert!((parsed.get(0, 1) - 0.123).abs() < 0.001);
        assert!((parsed.get(0, 2) - 0.456).abs() < 0.001);
        assert!((parsed.get(1, 2) - 0.789).abs() < 0.001);
        // Symmetric access
        assert!((parsed.get(2, 0) - 0.456).abs() < 0.001);
    }

    #[test]
    fn self_distance_is_zero() {
        let mat = sample_matrix();
        assert_eq!(mat.get(0, 0), 0.0);
        assert_eq!(mat.get(1, 1), 0.0);
    }
}
