use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

use mafft_types::{SeqName, Sequence, SequenceSet};

use crate::detect::detect_seq_type;
use crate::error::IoError;

/// Default line width for FASTA output (matches C macro `C = 60`).
const DEFAULT_LINE_WIDTH: usize = 60;

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

/// Read a FASTA file from a path into a `SequenceSet`.
///
/// The file is bytes. A non-UTF-8 byte is not an error: names keep the
/// header bytes, and the residue filter below decides which sequence
/// bytes survive. Every CR is turned into LF first, matching the C
/// wrapper's `tr "\r" "\n"`, so a lone CR is a line break.
///
/// Performs the same normalization as the C code:
/// - Strips non-alphabetic characters (except '-', '.') from sequences.
/// - Converts '*' to '-'.
/// - Auto-detects DNA vs protein via ATGC frequency.
/// - Canonicalises residue case per the detected type — lowercase for
///   DNA/RNA, uppercase for protein (see [`apply_case_convention`]).
///
/// Uses a lenient parser that handles MAFFT's non-standard headers
/// (e.g. `>     1== name ...` with leading spaces).
pub fn read_fasta(path: impl AsRef<Path>) -> Result<SequenceSet, IoError> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    read_fasta_from_reader(reader)
}

/// Read a FASTA file preserving case and non-standard residues — used
/// by `--anysymbol`/`--preservecase`. Strips only newline, space and
/// carriage return (matching C MAFFT's `readData_pointer_casepreserve`
/// → `charfilter`, `io.c:1329-1352`); any other character — digits,
/// tabs, punctuation, lowercase — is kept verbatim so the
/// post-alignment restore pass can put the originals back.
pub fn read_fasta_casepreserve(path: impl AsRef<Path>) -> Result<SequenceSet, IoError> {
    let file = std::fs::File::open(path)?;
    let reader = BufReader::new(file);
    read_fasta_from_reader_casepreserve(reader)
}

/// Like `read_fasta_from_reader` but preserves case and non-standard
/// residues (see `read_fasta_casepreserve`).
pub fn read_fasta_from_reader_casepreserve<R: BufRead>(reader: R) -> Result<SequenceSet, IoError> {
    read_fasta_bytes(reader, true)
}

/// Case-preserving sequence normaliser — mirrors C `charfilter`
/// (`io.c:1329-1352`, reached via `load1SeqWithoutName_realloc_casepreserve`):
/// drops only `\n`, space and `\r`; rejects `=`, `<`, `>`; keeps every
/// other byte — `*`, `@`, digits, tabs, lowercase, IUPAC — so
/// `--anysymbol`/`--preservecase` can replace then restore them. C keeps
/// digits and tabs as (unusual) residues here, so a GenBank-style
/// numbered sequence gains `X`/`n` columns under `--anysymbol`; we
/// reproduce that rather than second-guess it.
fn normalize_sequence_casepreserve(raw: &[u8]) -> Result<Vec<u8>, IoError> {
    if raw.iter().any(|&c| c == b'=' || c == b'<' || c == b'>') {
        return Err(IoError::IllegalTitleCharInSequence);
    }
    Ok(raw
        .iter()
        .copied()
        .filter(|&c| keep_casepreserve(c))
        .collect())
}

/// Bytes C's `charfilter` (`io.c:1329-1352`) keeps: everything except
/// `\n`, space and `\r`. Digits, tabs and punctuation survive as
/// (unusual) residues; `=`, `<`, `>` pass this test but are rejected
/// by [`normalize_sequence_casepreserve`] / [`is_title_only_char`].
fn keep_casepreserve(c: u8) -> bool {
    c != b'\n' && c != b' ' && c != b'\r'
}

/// `=`, `<`, `>` — legal only in description lines on the
/// case-preserving path (`charfilter`, `io.c:1337`).
fn is_title_only_char(c: u8) -> bool {
    c == b'=' || c == b'<' || c == b'>'
}

/// C MAFFT refuses input where a description line is preceded by blanks
/// (`scripts/mafft:1827-1834`: `grep -c '^[[:blank:]]\+>'` → exit 1).
/// Without this check a lenient line parser would treat the line as
/// sequence data and silently glue two records together.
///
/// `text` in the error is lossy UTF-8: it is the diagnostic C prints,
/// not a sequence or a header that has to round-trip.
fn reject_blank_before_header(idx: usize, line: &[u8]) -> Result<(), IoError> {
    let trimmed = trim_ascii_space_tab(line);
    if trimmed.len() != line.len() && trimmed.first() == Some(&b'>') {
        return Err(IoError::BlankBeforeHeader {
            line: idx + 1,
            text: String::from_utf8_lossy(line).into_owned(),
        });
    }
    Ok(())
}

fn trim_ascii_space_tab(line: &[u8]) -> &[u8] {
    let mut i = 0;
    while i < line.len() && (line[i] == b' ' || line[i] == b'\t') {
        i += 1;
    }
    &line[i..]
}

/// Read FASTA from raw bytes.
///
/// The file is not text. C's wrapper runs `tr "\r" "\n"`
/// (`scripts/mafft:1132`) before any engine sees it, so every CR,
/// including the one in CRLF, becomes LF and a CR-only file splits into
/// one record per line. A NUL does not end a sequence. C's string
/// handling stops at it and drops the rest of the record (`io.c:1355-1377`,
/// `io.c:1751`); this reader keeps those bytes,
/// which is an improvement. The residue filter then drops a NUL inside
/// a sequence and keeps what follows (`AC\0DE` stays `ACDE`). A NUL in
/// a name is kept, with the bytes after it.
fn read_fasta_bytes<R: BufRead>(mut reader: R, casepreserve: bool) -> Result<SequenceSet, IoError> {
    let mut raw = Vec::new();
    reader.read_to_end(&mut raw)?;
    for b in &mut raw {
        if *b == b'\r' {
            *b = b'\n';
        }
    }

    let mut sequences = Vec::new();
    let mut current_name: Option<SeqName> = None;
    let mut current_seq = Vec::new();

    for (idx, line) in raw.split(|b| *b == b'\n').enumerate() {
        reject_blank_before_header(idx, line)?;
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(name) = current_name.take() {
                sequences.push(Sequence {
                    name,
                    data: finish_sequence(&current_seq, casepreserve)?,
                });
                current_seq.clear();
            }
            current_name = Some(SeqName::from(header));
        } else if current_name.is_some() {
            current_seq.extend_from_slice(line);
        }
    }
    if let Some(name) = current_name.take() {
        sequences.push(Sequence {
            name,
            data: finish_sequence(&current_seq, casepreserve)?,
        });
    }
    if sequences.is_empty() {
        return Err(IoError::EmptyInput);
    }

    let seq_type = detect_seq_type(sequences.iter().map(|s| &s.data));
    let mut set = SequenceSet { sequences, seq_type };
    if !casepreserve {
        apply_case_convention(&mut set);
    }
    Ok(set)
}

fn finish_sequence(raw: &[u8], casepreserve: bool) -> Result<Vec<u8>, IoError> {
    if casepreserve {
        normalize_sequence_casepreserve(raw)
    } else {
        Ok(normalize_sequence(raw))
    }
}

/// Apply the FASTA reader's residue filter to residues that are already in
/// memory, so an in-memory caller ends up with exactly the bytes a FASTA
/// round trip would have produced — and fails exactly where the reader
/// would.
///
/// `casepreserve = false` is the default reader's rule (keep letters, `-`
/// and `.`, turn `*` into `-`, drop everything else — see
/// [`read_fasta`]; it cannot fail); `casepreserve = true` is the
/// `--anysymbol` / `--preservecase` rule (drop only `\n`, space and `\r`,
/// reject `=`, `<`, `>` with [`IoError::IllegalTitleCharInSequence`] — see
/// [`read_fasta_casepreserve`]). Neither touches case; that is
/// [`apply_case_convention`]'s job and depends on the sequence type.
pub fn normalize_residues(raw: &[u8], casepreserve: bool) -> Result<Vec<u8>, IoError> {
    if casepreserve {
        normalize_sequence_casepreserve(raw)
    } else {
        Ok(normalize_sequence(raw))
    }
}

/// `true` when [`normalize_residues`] would return `raw` unchanged, i.e.
/// the residues already look like they came out of the FASTA reader.
/// Lets a caller holding borrowed data skip the copy when nothing needs
/// to change. Residues the reader would reject are not "unchanged".
pub fn residues_are_normalized(raw: &[u8], casepreserve: bool) -> bool {
    if casepreserve {
        raw.iter().all(|&c| keep_casepreserve(c) && !is_title_only_char(c))
    } else {
        // `*` is not dropped but it is rewritten, so it is not "unchanged".
        raw.iter().all(|&c| c.is_ascii_alphabetic() || c == b'-' || c == b'.')
    }
}

/// `true` when [`apply_case_convention`] would leave `data` unchanged for
/// a set of type `seq_type`: no uppercase letters for nucleotides, no
/// lowercase letters otherwise.
pub fn residues_follow_case_convention(data: &[u8], seq_type: mafft_types::SeqType) -> bool {
    if seq_type.is_nucleotide() {
        !data.iter().any(|c| c.is_ascii_uppercase())
    } else {
        !data.iter().any(|c| c.is_ascii_lowercase())
    }
}

/// Read FASTA from any buffered reader.
///
/// The stream is bytes, not UTF-8 text. Handles non-standard headers
/// with leading whitespace that strict parsers (like `noodles-fasta`)
/// reject. The bytes after `>` are preserved as the sequence name,
/// including a non-UTF-8 byte and the bytes after a NUL. A lone CR is a
/// line break, matching the C wrapper's `tr "\r" "\n"`.
pub fn read_fasta_from_reader<R: BufRead>(reader: R) -> Result<SequenceSet, IoError> {
    read_fasta_bytes(reader, false)
}

/// Apply C MAFFT's residue-case convention to an already-parsed set:
/// lowercase for DNA/RNA, uppercase for everything else.
///
/// C canonicalises case as it reads — `io.c:1462-1467`
/// (`load1SeqWithoutName_realloc`) calls `onlyAlpha_lower` when
/// `dorp == 'd'` and `onlyAlpha_upper` otherwise, and `readData_pointer`
/// repeats the nucleotide pass with `seqLower` (`io.c:1755`). The
/// `upperCase != -1` guard there is only reachable from the legacy
/// non-FASTA `FRead` header parser (`io.c:1174-1184`), so for FASTA input
/// it is always true. Net effect: C MAFFT's default output is lowercase
/// for DNA/RNA and uppercase for protein, whatever case the input used.
///
/// The fold is idempotent, so it is safe to re-apply after `--nuc` /
/// `--amino` override the detected type — which is what C does, since
/// `$seqtype` fixes `dorp` before any sequence is read
/// (`scripts/mafft:547-550`).
///
/// Deliberately NOT applied by [`read_fasta_casepreserve`]: on the
/// `--anysymbol` / `--preservecase` path C reads with
/// `readData_pointer_casepreserve` and restores the original characters
/// after alignment (`replaceu` + `restoreu`), so the input case survives.
pub fn apply_case_convention(set: &mut SequenceSet) {
    let nucleotide = set.seq_type.is_nucleotide();
    for seq in set.sequences.iter_mut() {
        for ch in seq.data.iter_mut() {
            *ch = if nucleotide {
                ch.to_ascii_lowercase()
            } else {
                ch.to_ascii_uppercase()
            };
        }
    }
}

/// Normalize a raw sequence: keep only alpha + gap chars, convert '*' to '-'.
///
/// Mirrors the character-filtering half of C's `onlyAlpha_lower()` /
/// `onlyAlpha_upper()` plus `kake2hiku()` (`io.c:1425-1470`). Case is left
/// alone here because C picks the case fold from `dorp`, which is only
/// known once the sequence type has been detected (or forced by
/// `--nuc` / `--amino`); [`apply_case_convention`] applies it afterwards.
fn normalize_sequence(raw: &[u8]) -> Vec<u8> {
    raw.iter()
        .filter_map(|&ch| {
            if ch.is_ascii_alphabetic() {
                Some(ch)
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
        writer.write_all(b">")?;
        writer.write_all(seq.name.as_bytes())?;
        writeln!(writer)?;

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

    // --- C MAFFT residue-case convention (io.c:1462-1467, io.c:1755) ---

    #[test]
    fn nucleotide_input_is_lowercased_whatever_the_input_case() {
        let input = b">a\nATGGCtagcTTGGACCATTGCAGG\n>b\nATGGCTAGCTTGGACCATTGCAGG\n";
        let seqs = read_fasta_from_reader(io::Cursor::new(&input[..])).unwrap();
        assert_eq!(seqs.seq_type, mafft_types::SeqType::Dna);
        assert_eq!(seqs.sequences[0].data, b"atggctagcttggaccattgcagg".to_vec());
        assert_eq!(seqs.sequences[1].data, b"atggctagcttggaccattgcagg".to_vec());
    }

    #[test]
    fn protein_input_is_uppercased_whatever_the_input_case() {
        let input = b">a\nMNGTegdnFYVPFSNKTGLARSPYEY\n>b\nMNGTEGDNFYVPFSNKTGLARSPYEY\n";
        let seqs = read_fasta_from_reader(io::Cursor::new(&input[..])).unwrap();
        assert_eq!(seqs.seq_type, mafft_types::SeqType::Protein);
        assert_eq!(seqs.sequences[0].data, b"MNGTEGDNFYVPFSNKTGLARSPYEY".to_vec());
    }

    #[test]
    fn casepreserve_reader_keeps_the_input_case() {
        // `--anysymbol` / `--preservecase` restore the originals after
        // alignment, so this reader must not fold anything.
        let input = b">a\nATGGCtagcTTGGACCATTGCAGG\n";
        let seqs = read_fasta_from_reader_casepreserve(io::Cursor::new(&input[..])).unwrap();
        assert_eq!(seqs.sequences[0].data, b"ATGGCtagcTTGGACCATTGCAGG".to_vec());
    }

    #[test]
    fn apply_case_convention_is_idempotent_and_follows_seq_type() {
        // Safe to re-apply after `--nuc` / `--amino` override the type.
        let mut set = SequenceSet {
            sequences: vec![Sequence { name: "a".into(), data: b"AtGc".to_vec() }],
            seq_type: mafft_types::SeqType::Dna,
        };
        apply_case_convention(&mut set);
        assert_eq!(set.sequences[0].data, b"atgc".to_vec());
        apply_case_convention(&mut set);
        assert_eq!(set.sequences[0].data, b"atgc".to_vec());

        set.seq_type = mafft_types::SeqType::Protein;
        apply_case_convention(&mut set);
        assert_eq!(set.sequences[0].data, b"ATGC".to_vec());
    }

    #[test]
    fn high_byte_in_the_sequence_is_dropped_like_any_other_non_letter() {
        // Latin-1 e-acute sits between residues. C's default filter
        // (`onlyAlpha_*`) drops it; the reader used to refuse the file
        // as invalid UTF-8 before it got that far.
        let with = read_fasta_from_reader(io::Cursor::new(&b">a\nACDE\xe9FGHIK\n>b\nACDEFGHIKL\n"[..])).unwrap();
        let without = read_fasta_from_reader(io::Cursor::new(&b">a\nACDEFGHIK\n>b\nACDEFGHIKL\n"[..])).unwrap();
        assert_eq!(with.nseq(), 2);
        assert_eq!(with.sequences[0].data, without.sequences[0].data);
        assert_eq!(with.sequences[1].data, without.sequences[1].data);
    }

    #[test]
    fn casepreserve_keeps_a_high_byte_in_the_sequence() {
        let seqs = read_fasta_from_reader_casepreserve(io::Cursor::new(&b">a\nAC\xe9DE\n"[..])).unwrap();
        assert_eq!(seqs.sequences[0].data, b"AC\xe9DE");
    }

    #[test]
    fn non_utf8_name_roundtrips_as_the_original_bytes() {
        let seqs = read_fasta_from_reader(io::Cursor::new(&b">a\xff\nACDEFGHIKLMNPQRSTVWY\n"[..])).unwrap();
        assert_eq!(seqs.sequences[0].name.as_bytes(), b"a\xff");
        let mut buf = Vec::new();
        write_fasta_to_writer(&seqs, &mut buf).unwrap();
        assert!(buf.starts_with(b">a\xff\n"), "{buf:?}");
        let again = read_fasta_from_reader(io::Cursor::new(&buf)).unwrap();
        assert_eq!(again.sequences[0].name.as_bytes(), b"a\xff");
    }

    #[test]
    fn nul_does_not_end_the_line() {
        // Improvement over C: C stops at a NUL and drops the rest of the
        // record. This reader drops a NUL inside a sequence (it is
        // not a letter) and keeps the residues after it, and it keeps a
        // NUL and the following bytes in a name.
        let seq = read_fasta_from_reader(io::Cursor::new(&b">ab\x00cd\nAC\x00DE\n"[..])).unwrap();
        let plain = read_fasta_from_reader(io::Cursor::new(&b">abcd\nACDE\n"[..])).unwrap();
        assert_eq!(seq.sequences[0].name.as_bytes(), b"ab\x00cd");
        assert_eq!(seq.sequences[0].data, plain.sequences[0].data);
    }

    #[test]
    fn lone_cr_is_a_line_break_and_crlf_stays_one_record_boundary() {
        let lf = b">a\nACDEFGHIK\n>b\nACDEFGHIKL\n";
        let cr = b">a\rACDEFGHIK\r>b\rACDEFGHIKL\r";
        let crlf = b">a\r\nACDEFGHIK\r\n>b\r\nACDEFGHIKL\r\n";
        let mixed = b">a\r\nACDEFGHIK\n>b\rACDEFGHIKL\n";
        let via_lf = read_fasta_from_reader(io::Cursor::new(&lf[..])).unwrap();
        for bytes in [&cr[..], &crlf[..], &mixed[..]] {
            let got = read_fasta_from_reader(io::Cursor::new(bytes)).unwrap();
            assert_eq!(got.nseq(), 2, "{bytes:?}");
            assert_eq!(got.sequences[0].name, via_lf.sequences[0].name);
            assert_eq!(got.sequences[1].name, via_lf.sequences[1].name);
            assert_eq!(got.sequences[0].data, via_lf.sequences[0].data);
            assert_eq!(got.sequences[1].data, via_lf.sequences[1].data);
        }
    }

    #[test]
    fn case_fold_does_not_disturb_type_detection() {
        // Detection runs on the pre-fold residues and is case-insensitive,
        // so a lowercase and an uppercase copy detect the same type.
        let upper = read_fasta_from_reader(io::Cursor::new(&b">a\nACGTACGTACGTACGT\n"[..])).unwrap();
        let lower = read_fasta_from_reader(io::Cursor::new(&b">a\nacgtacgtacgtacgt\n"[..])).unwrap();
        assert_eq!(upper.seq_type, lower.seq_type);
        assert_eq!(upper.sequences[0].data, lower.sequences[0].data);
    }
}

/// The in-memory helpers must agree with the readers byte for byte: an
/// in-memory caller relies on them to reproduce a FASTA round trip.
#[cfg(test)]
mod in_memory_helper_tests {
    use super::*;
    use mafft_types::SeqType;

    #[test]
    fn normalize_residues_matches_the_readers() {
        let raw = b"MNG*T.E-G 123\t\r\n@x";
        assert_eq!(normalize_residues(raw, false).unwrap(), normalize_sequence(raw));
        assert_eq!(
            normalize_residues(raw, true).unwrap(),
            normalize_sequence_casepreserve(raw).unwrap()
        );
        assert_eq!(normalize_residues(raw, false).unwrap(), b"MNG-T.E-Gx");
        // `charfilter`: digits and the tab are residues; only \n, space, \r go.
        assert_eq!(normalize_residues(raw, true).unwrap(), b"MNG*T.E-G123\t@x");
        // `= < >` inside a sequence are fatal on the case-preserving path
        // only, exactly as in the reader.
        assert!(matches!(
            normalize_residues(b"MN=G", true),
            Err(IoError::IllegalTitleCharInSequence)
        ));
        assert_eq!(normalize_residues(b"MN=G", false).unwrap(), b"MNG");
    }

    #[test]
    fn residues_are_normalized_iff_normalize_is_identity() {
        for casepreserve in [false, true] {
            for raw in [
                &b"ACGT-acgt."[..],
                b"MNG*T",
                b"AC GT",
                b"AC1GT",
                b"AC\tGT",
                b"@x",
                b"A=C",
                b"A<C>",
                b"",
            ] {
                let identity = normalize_residues(raw, casepreserve).is_ok_and(|v| v == raw);
                assert_eq!(
                    residues_are_normalized(raw, casepreserve),
                    identity,
                    "casepreserve={casepreserve} raw={raw:?}"
                );
            }
        }
    }

    #[test]
    fn case_convention_check_matches_apply() {
        for (data, seq_type) in [
            (&b"acgt-"[..], SeqType::Dna),
            (b"ACGT-", SeqType::Dna),
            (b"MKV-", SeqType::Protein),
            (b"mkv-", SeqType::Protein),
            (b"MkV", SeqType::Unknown),
        ] {
            let mut set = SequenceSet {
                sequences: vec![Sequence { name: "s".into(), data: data.to_vec() }],
                seq_type,
            };
            apply_case_convention(&mut set);
            let unchanged = set.sequences[0].data == data;
            assert_eq!(residues_follow_case_convention(data, seq_type), unchanged, "{data:?}");
        }
    }

    #[test]
    fn detect_accepts_borrowed_rows() {
        let set = SequenceSet {
            sequences: vec![Sequence { name: "s".into(), data: b"ATGCGATCGATCG".to_vec() }],
            seq_type: SeqType::Unknown,
        };
        let borrowed = detect_seq_type(set.sequences.iter().map(|s| &s.data));
        let owned: Vec<Vec<u8>> = set.sequences.iter().map(|s| s.data.clone()).collect();
        assert_eq!(borrowed, detect_seq_type(&owned));
        assert_eq!(borrowed, SeqType::Dna);
    }
}
