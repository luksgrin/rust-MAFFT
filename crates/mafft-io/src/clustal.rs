use std::io::Write;

use mafft_types::SequenceSet;

use crate::error::IoError;

/// Default block width for Clustal output (matches C code: 60).
const BLOCK_WIDTH: usize = 60;

/// Default name field width.
const DEFAULT_NAME_LEN: usize = 15;

/// Write an alignment in Clustal format.
///
/// Matches the output of the C `clustalout_pointer()` function.
///
/// - `order`: optional sequence ordering (indices into `seqs.sequences`).
///   If `None`, sequences are written in their natural order.
/// - `name_len`: field width for names. Pass `None` for the default (15).
/// - `conservation`: optional conservation mark string (same length as alignment).
pub fn write_clustal<W: Write>(
    seqs: &SequenceSet,
    writer: &mut W,
    order: Option<&[usize]>,
    name_len: Option<usize>,
    conservation: Option<&str>,
) -> Result<(), IoError> {
    let nseq = seqs.nseq();
    if nseq == 0 {
        return Err(IoError::EmptyInput);
    }

    let name_len = name_len.unwrap_or(DEFAULT_NAME_LEN);
    let max_len = seqs.max_len();

    // Header
    writeln!(writer, "CLUSTAL format alignment by MAFFT (v7.526)")?;
    writeln!(writer)?;

    let default_order: Vec<usize> = (0..nseq).collect();
    let order = order.unwrap_or(&default_order);

    let mut pos = 0;
    while pos < max_len {
        writeln!(writer)?;

        let end = (pos + BLOCK_WIDTH).min(max_len);

        for &idx in order {
            let seq = &seqs.sequences[idx];
            let first_word = first_word_of(&seq.name);
            write!(writer, "{:<width$} ", first_word, width = name_len)?;

            let chunk_end = end.min(seq.data.len());
            if pos < chunk_end {
                writer.write_all(&seq.data[pos..chunk_end])?;
            }
            writeln!(writer)?;
        }

        if let Some(marks) = conservation {
            let mark_bytes = marks.as_bytes();
            write!(writer, "{:<width$} ", "", width = name_len)?;
            let chunk_end = end.min(mark_bytes.len());
            if pos < chunk_end {
                writer.write_all(&mark_bytes[pos..chunk_end])?;
            }
            writeln!(writer)?;
        }

        pos += BLOCK_WIDTH;
    }

    Ok(())
}

/// Extract the first whitespace-delimited word from a name string.
fn first_word_of(name: &str) -> &str {
    name.split_whitespace().next().unwrap_or(name)
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_types::{SeqType, Sequence};

    #[test]
    fn clustal_output_format() {
        let seqs = SequenceSet {
            sequences: vec![
                Sequence {
                    name: "seq1 some description".into(),
                    data: b"ACGT-ACGT".to_vec(),
                },
                Sequence {
                    name: "seq2".into(),
                    data: b"ACGTAAC-T".to_vec(),
                },
            ],
            seq_type: SeqType::Dna,
        };

        let mut buf = Vec::new();
        write_clustal(&seqs, &mut buf, None, None, None).unwrap();
        let output = String::from_utf8(buf).unwrap();

        assert!(output.starts_with("CLUSTAL format alignment by MAFFT"));
        assert!(output.contains("seq1"));
        assert!(output.contains("seq2"));
        assert!(output.contains("ACGT-ACGT"));
    }
}
