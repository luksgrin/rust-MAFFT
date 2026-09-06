use mafft_types::SeqType;

/// Detect whether sequences are DNA/RNA or protein by ATGC frequency.
///
/// Mirrors the C `getnumlen()` heuristic: if >75% of the first `limit`
/// residues are A, T, G, C, U or N (case-insensitive), the input is
/// nucleotide. `N` counts as nucleotide exactly as in C `countATGC`
/// (`io.c:2065-2090`), so N-rich DNA is not mistaken for protein.
///
/// Accepts anything that iterates over byte slices — `&[Vec<u8>]`,
/// `&Vec<Vec<u8>>`, or an iterator of borrowed rows such as
/// `set.sequences.iter().map(|s| &s.data)` — so callers need not copy the
/// residues into a fresh `Vec<Vec<u8>>` just to ask the question.
///
/// Only alphabetic bytes are counted (gap characters and everything else
/// are skipped), so the answer is the same whether the residues have been
/// through the FASTA reader's normalisation or not.
pub fn detect_seq_type<I>(sequences: I) -> SeqType
where
    I: IntoIterator,
    I::Item: AsRef<[u8]>,
{
    detect_seq_type_with_limit(sequences, 1_000_000)
}

pub fn detect_seq_type_with_limit<I>(sequences: I, limit: usize) -> SeqType
where
    I: IntoIterator,
    I::Item: AsRef<[u8]>,
{
    let mut atgc = 0u64;
    let mut total = 0u64;

    'outer: for seq in sequences {
        for &ch in seq.as_ref() {
            match ch.to_ascii_uppercase() {
                b'A' | b'T' | b'G' | b'C' | b'U' | b'N' => {
                    atgc += 1;
                    total += 1;
                }
                b'-' | b'.' | b'*' => {} // skip gaps
                _ if ch.is_ascii_alphabetic() => {
                    total += 1;
                }
                _ => {}
            }
            if total as usize >= limit {
                break 'outer;
            }
        }
    }

    if total == 0 {
        return SeqType::Unknown;
    }

    let freq = atgc as f64 / total as f64;
    if freq > 0.75 {
        SeqType::Dna
    } else {
        SeqType::Protein
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_protein() {
        let seqs = vec![b"MNGTEGDNFYVPFSNKTGLAR".to_vec()];
        assert_eq!(detect_seq_type(&seqs), SeqType::Protein);
    }

    #[test]
    fn detects_dna() {
        let seqs = vec![b"ATGCGATCGATCGATCGATCG".to_vec()];
        assert_eq!(detect_seq_type(&seqs), SeqType::Dna);
    }

    #[test]
    fn empty_is_unknown() {
        let seqs: Vec<Vec<u8>> = vec![];
        assert_eq!(detect_seq_type(&seqs), SeqType::Unknown);
    }
}
