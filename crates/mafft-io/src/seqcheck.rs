//! C `seqcheck` (`mltaln9.c:60-85`): the alphabet test every aligner
//! driver runs on the sequences it has just read (`disttbfast.c:3470`,
//! `tbfast.c:2483`, `pairlocalalign.c:3229`, `dvtditr.c:714`).
//!
//! This is the one place that decides which residues the program accepts
//! for a given sequence type; the file reader and the in-memory entry
//! point of `mafft-rs` both go through it, so an input can never be legal
//! on one path and fatal on the other. The alphabets are C's `locaminod`
//! (`blosum.c:12`, `ARNDCQEGHILKMFPSTWYVBZX.-J`) and `locaminon`
//! (`DNA.h:43`, `agctuAGCTUnNbdhkmnrsvwyx-O`), taken from `mafft_scoring`
//! rather than restated here.
//!
//! Consequences worth knowing: protein `U`/`O` are illegal (C tells the
//! user to try `--anysymbol`), protein `.` is a legal residue, and
//! nucleotide `.` is illegal. C runs the check after the reader has folded
//! case (`io.c:1462-1467`), so `dorp`-dependent legality (e.g. nucleotide
//! `O` → `o` → illegal) falls out of the same test; callers here are
//! expected to do the same and pass an already case-folded set.

use mafft_types::SequenceSet;

/// The first residue `seqcheck` would reject, in reading order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IllegalResidue {
    /// 0-based index of the sequence in the set.
    pub seq_index: usize,
    /// 0-based position of the residue within that sequence.
    pub position: usize,
    /// The offending byte.
    pub byte: u8,
}

/// The alphabet `seqcheck` tests `set` against — `mafft_scoring`'s DNA
/// alphabet for nucleotide sets, its protein alphabet otherwise.
pub fn seqcheck_alphabet(set: &SequenceSet) -> &'static [u8] {
    if set.seq_type.is_nucleotide() {
        &mafft_scoring::DNA_ALPHABET.chars
    } else {
        &mafft_scoring::PROTEIN_ALPHABET.chars
    }
}

/// Scan `set` the way C `seqcheck` does and return the first residue that
/// is not in the alphabet for `set.seq_type`, or `None` when every
/// residue is legal. Gap characters (`-`) are part of both alphabets, so
/// a gapped input passes.
pub fn find_illegal_residue(set: &SequenceSet) -> Option<IllegalResidue> {
    let alphabet = seqcheck_alphabet(set);
    for (seq_index, s) in set.sequences.iter().enumerate() {
        if let Some((position, &byte)) =
            s.data.iter().enumerate().find(|(_, c)| !alphabet.contains(c))
        {
            return Some(IllegalResidue { seq_index, position, byte });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_types::{SeqType, Sequence};

    fn set(rows: &[&[u8]], seq_type: SeqType) -> SequenceSet {
        SequenceSet {
            sequences: rows
                .iter()
                .enumerate()
                .map(|(i, d)| Sequence { name: format!("s{i}"), data: d.to_vec() })
                .collect(),
            seq_type,
        }
    }

    #[test]
    fn protein_alphabet_is_c_locaminod() {
        // Legal: the 20 amino acids, B Z X J, `.` and `-`.
        assert_eq!(find_illegal_residue(&set(&[b"ARNDCQEGHILKMFPSTWYVBZXJ.-"], SeqType::Protein)), None);
        // U and O are fatal (site 6 of sequence 1 in the corpus fixture).
        assert_eq!(
            find_illegal_residue(&set(&[b"MKALVUWQ", b"MKALVOWQ"], SeqType::Protein)),
            Some(IllegalResidue { seq_index: 0, position: 5, byte: b'U' })
        );
        assert_eq!(
            find_illegal_residue(&set(&[b"MKALVWQ", b"MKALVOWQ"], SeqType::Protein)),
            Some(IllegalResidue { seq_index: 1, position: 5, byte: b'O' })
        );
        // Lowercase is illegal too: C checks after the case fold.
        assert_eq!(
            find_illegal_residue(&set(&[b"MKalv"], SeqType::Protein)).map(|r| r.byte),
            Some(b'a')
        );
    }

    #[test]
    fn nucleotide_alphabet_is_c_locaminon() {
        assert_eq!(find_illegal_residue(&set(&[b"agctunbdhkmrsvwyx-"], SeqType::Dna)), None);
        assert_eq!(find_illegal_residue(&set(&[b"agcu"], SeqType::Rna)), None);
        // `.` is not a nucleotide.
        assert_eq!(
            find_illegal_residue(&set(&[b"acgt", b"acg.t"], SeqType::Dna)),
            Some(IllegalResidue { seq_index: 1, position: 3, byte: b'.' })
        );
        // Uppercase `O` is in `locaminon`, lowercase `o` is not — which is
        // what a case-folded nucleotide `O` becomes.
        assert_eq!(find_illegal_residue(&set(&[b"acgtO"], SeqType::Dna)), None);
        assert_eq!(find_illegal_residue(&set(&[b"acgto"], SeqType::Dna)).map(|r| r.byte), Some(b'o'));
    }

    #[test]
    fn unknown_type_uses_the_protein_alphabet() {
        assert_eq!(find_illegal_residue(&set(&[b"ACGT."], SeqType::Unknown)), None);
        assert_eq!(find_illegal_residue(&set(&[b"acgt"], SeqType::Unknown)).map(|r| r.byte), Some(b'a'));
    }
}
