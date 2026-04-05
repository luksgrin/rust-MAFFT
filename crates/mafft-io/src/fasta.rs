use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use mafft_types::{Sequence, SequenceSet};

use crate::detect::detect_seq_type;
use crate::error::IoError;

/// Default line width for FASTA output (matches C macro `C = 60`).
const DEFAULT_LINE_WIDTH: usize = 60;

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

/// Read a FASTA file from a path into a `SequenceSet`.
///
/// Performs the same normalization as the C code:
/// - Strips non-alphabetic characters (except '-', '.') from sequences.
/// - Converts '*' to '-'.
/// - Auto-detects DNA vs protein via ATGC frequency.
///
/// Uses a lenient parser that handles MAFFT's non-standard headers
/// (e.g. `>     1== name ...` with leading spaces).
pub fn read_fasta(path: impl AsRef<Path>) -> Result<SequenceSet, IoError> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    read_fasta_from_reader(reader)
}

/// Read FASTA from any buffered reader.
///
/// Handles non-standard headers with leading whitespace that strict parsers
/// (like `noodles-fasta`) reject. The full text after '>' is preserved as
/// the sequence name, matching MAFFT's C behavior.
pub fn read_fasta_from_reader<R: BufRead>(reader: R) -> Result<SequenceSet, IoError> {
    let mut sequences = Vec::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = Vec::new();

    for line_result in reader.lines() {
        let line = line_result?;

        if let Some(header) = line.strip_prefix('>') {
            // Flush previous sequence
            if let Some(name) = current_name.take() {
                sequences.push(Sequence {
                    name,
                    data: normalize_sequence(&current_seq),
                });
                current_seq.clear();
            }
            current_name = Some(header.to_string());
        } else if current_name.is_some() {
            // Sequence data line
            current_seq.extend_from_slice(line.as_bytes());
        }
        // Lines before the first '>' are ignored
    }

    // Flush last sequence
    if let Some(name) = current_name.take() {
        sequences.push(Sequence {
            name,
            data: normalize_sequence(&current_seq),
        });
    }

    if sequences.is_empty() {
        return Err(IoError::EmptyInput);
    }

    let seq_type = detect_seq_type(
        &sequences.iter().map(|s| s.data.clone()).collect::<Vec<_>>(),
    );

    Ok(SequenceSet {
        sequences,
        seq_type,
    })
}

/// Normalize a raw sequence: keep only alpha + gap chars, convert '*' to '-'.
///
/// Mirrors the C `onlyAlpha_lower()` + `kake2hiku()` pipeline.
fn normalize_sequence(raw: &[u8]) -> Vec<u8> {
    raw.iter()
        .filter_map(|&ch| {
            if ch.is_ascii_alphabetic() {
                Some(ch.to_ascii_uppercase())
            } else if ch == b'-' || ch == b'.' {
                Some(ch)
            } else if ch == b'*' {
                Some(b'-') // kake2hiku: * → -
            } else {
                None // strip digits, whitespace, etc.
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// Write a `SequenceSet` as FASTA to a file path.
pub fn write_fasta(seqs: &SequenceSet, path: impl AsRef<Path>) -> Result<(), IoError> {
    let file = std::fs::File::create(path)?;
    let writer = io::BufWriter::new(file);
    write_fasta_to_writer(seqs, writer)
}

/// Write a `SequenceSet` as FASTA to any writer.
///
/// Uses 60-character line width by default (matching the C output).
pub fn write_fasta_to_writer<W: Write>(
    seqs: &SequenceSet,
    mut writer: W,
) -> Result<(), IoError> {
    write_fasta_to_writer_with_width(seqs, &mut writer, DEFAULT_LINE_WIDTH)
}

/// Write FASTA with a custom line width. Pass `0` for unlimited (single line).
pub fn write_fasta_to_writer_with_width<W: Write>(
    seqs: &SequenceSet,
    writer: &mut W,
    line_width: usize,
) -> Result<(), IoError> {
    for seq in &seqs.sequences {
        writeln!(writer, ">{}", seq.name)?;

        if line_width == 0 {
            writer.write_all(&seq.data)?;
            writeln!(writer)?;
        } else {
            for chunk in seq.data.chunks(line_width) {
                writer.write_all(chunk)?;
                writeln!(writer)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_strips_and_converts() {
        let raw = b"MNG*T.E-G 123\n";
        let result = normalize_sequence(raw);
        assert_eq!(result, b"MNG-T.E-G");
    }

    #[test]
    fn roundtrip_fasta() {
        let original = SequenceSet {
            sequences: vec![
                Sequence {
                    name: "seq1 description".into(),
                    data: b"ACGTACGTACGT".to_vec(),
                },
                Sequence {
                    name: "seq2".into(),
                    data: b"MNGTEGDNFYVP".to_vec(),
                },
            ],
            seq_type: mafft_types::SeqType::Protein,
        };

        let mut buf = Vec::new();
        write_fasta_to_writer(&original, &mut buf).unwrap();

        let parsed = read_fasta_from_reader(io::Cursor::new(&buf)).unwrap();
        assert_eq!(parsed.sequences.len(), 2);
        assert_eq!(parsed.sequences[0].name, "seq1 description");
        assert_eq!(parsed.sequences[0].data, b"ACGTACGTACGT");
        assert_eq!(parsed.sequences[1].name, "seq2");
        assert_eq!(parsed.sequences[1].data, b"MNGTEGDNFYVP");
    }

    #[test]
    fn handles_mafft_style_headers() {
        let input = b">     1== M63632 rhodopsin\nMNGTEGDNFYVP\n>     2== U22180 rat opsin\nACGT\n";
        let seqs = read_fasta_from_reader(io::Cursor::new(&input[..])).unwrap();
        assert_eq!(seqs.nseq(), 2);
        assert!(seqs.sequences[0].name.contains("M63632"));
        assert!(seqs.sequences[1].name.contains("U22180"));
    }
}
