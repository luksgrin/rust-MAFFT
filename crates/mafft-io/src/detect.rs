use mafft_types::SeqType;

/// Detect whether sequences are DNA/RNA or protein by ATGC frequency.
///
/// Mirrors the C `getnumlen()` heuristic: if >75% of the first `limit`
/// residues are A, T, G, C, U (case-insensitive), the input is nucleotide.
pub fn detect_seq_type(sequences: &[Vec<u8>]) -> SeqType {
    detect_seq_type_with_limit(sequences, 1_000_000)
}

pub fn detect_seq_type_with_limit(sequences: &[Vec<u8>], limit: usize) -> SeqType {
    let mut atgc = 0u64;
    let mut total = 0u64;

    'outer: for seq in sequences {
        for &ch in seq {
            match ch.to_ascii_uppercase() {
                b'A' | b'T' | b'G' | b'C' | b'U' => {
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
