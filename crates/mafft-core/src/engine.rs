/// High-level MAFFT alignment engine.

use rayon::prelude::*;

use mafft_io::read_fasta;
use mafft_scoring::build_context;
use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, pairwise_identity_distance};
use mafft_align::{build_local_homology_table, GapModel};
use mafft_types::{ScoringModel, SeqType, SequenceSet, LocalHomologyTable};

use crate::progressive::{progressive_align, MultipleAlignment};
use crate::refinement::{iterative_refine, RefinementParams};

/// Alignment mode (strategy).
#[derive(Debug, Clone)]
pub enum AlignmentMode {
    /// FFT-NS-2: fast progressive (default).
    FftNs2,
    /// FFT-NS-i: progressive + limited iterative refinement.
    FftNsi { iterations: usize },
    /// G-INS-i: global alignment + iterative refinement.
    GInsi { iterations: usize },
    /// L-INS-i: local alignment + iterative refinement.
    LInsi { iterations: usize },
    /// E-INS-i: generalized affine + iterative refinement.
    EInsi { iterations: usize },
}

impl Default for AlignmentMode {
    fn default() -> Self {
        Self::FftNs2
    }
}

/// The main MAFFT alignment engine.
#[derive(Debug, Clone)]
pub struct MafftEngine {
    pub mode: AlignmentMode,
    pub scoring_model: ScoringModel,
}

impl Default for MafftEngine {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::FftNs2,
            scoring_model: ScoringModel::Jtt,
        }
    }
}

impl MafftEngine {
    pub fn new(mode: AlignmentMode) -> Self {
        Self { mode, scoring_model: ScoringModel::Jtt }
    }

    /// Align a set of sequences.
    pub fn align(&self, input: &SequenceSet) -> MultipleAlignment {
        let seq_type = input.seq_type;
        let scoring_model = if seq_type.is_nucleotide() {
            ScoringModel::Dna
        } else {
            self.scoring_model
        };

        let scoring = build_context(scoring_model, seq_type);
        let nseq = input.nseq();

        // Step 1: Compute pairwise distances (parallel)
        let pairs: Vec<(usize, usize, f64)> = (0..nseq)
            .into_par_iter()
            .flat_map(|i| {
                let seqs = &input.sequences;
                ((i + 1)..nseq).into_par_iter().map(move |j| {
                    let d = pairwise_identity_distance(&seqs[i].data, &seqs[j].data);
                    (i, j, d)
                })
            })
            .collect();

        let mut dm = DistanceMatrix::new(nseq);
        for (i, j, d) in pairs {
            dm.set(i, j, d);
        }

        // Step 2: Build guide tree
        let topo = musclesupg(&dm, ClusterMethod::default());

        // Step 3: Build local homology table (for constrained modes)
        let uses_constraints = matches!(
            self.mode,
            AlignmentMode::LInsi { .. } | AlignmentMode::EInsi { .. }
        );
        let local_hom = if uses_constraints {
            let seq_refs: Vec<&[u8]> = input.sequences.iter().map(|s| s.data.as_slice()).collect();
            let gap = GapModel::new(scoring.gap.open as f64, scoring.gap.extend as f64);
            let (table, _dist) = build_local_homology_table(
                &seq_refs,
                &scoring.substitution_matrix,
                &scoring.amino_map,
                &gap,
                0.0,
            );
            Some(table)
        } else {
            None
        };

        // Step 4: Progressive alignment
        let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();
        let names: Vec<String> = input.sequences.iter().map(|s| s.name.clone()).collect();

        let use_fft = matches!(
            self.mode,
            AlignmentMode::FftNs2 | AlignmentMode::FftNsi { .. }
        );

        let mut msa = progressive_align(&sequences, &names, &topo, &scoring, use_fft);

        // Step 5: Iterative refinement (if mode requires it)
        match &self.mode {
            AlignmentMode::FftNs2 => {}
            AlignmentMode::FftNsi { iterations }
            | AlignmentMode::GInsi { iterations }
            | AlignmentMode::LInsi { iterations }
            | AlignmentMode::EInsi { iterations } => {
                let params = RefinementParams {
                    max_iterations: *iterations,
                    cut: 0.0001,
                    use_fft,
                };
                iterative_refine(
                    &mut msa,
                    &topo,
                    &scoring,
                    &params,
                    local_hom.as_ref(),
                );
            }
        }

        msa
    }

    /// Convenience: read FASTA file and align.
    pub fn align_file(&self, path: &std::path::Path) -> Result<MultipleAlignment, mafft_io::IoError> {
        let input = read_fasta(path)?;
        Ok(self.align(&input))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mafft_types::{Sequence, SeqType};

    fn make_test_input() -> SequenceSet {
        SequenceSet {
            sequences: vec![
                Sequence { name: "s1".into(), data: b"ACDEFGHIKLMNP".to_vec() },
                Sequence { name: "s2".into(), data: b"ACDEFHIKLMNP".to_vec() },
                Sequence { name: "s3".into(), data: b"ACDEHIKLMNP".to_vec() },
            ],
            seq_type: SeqType::Protein,
        }
    }

    #[test]
    fn engine_fftns2() {
        let engine = MafftEngine::new(AlignmentMode::FftNs2);
        let msa = engine.align(&make_test_input());
        assert_eq!(msa.nseq(), 3);
        let w = msa.width();
        assert!(w >= 13);
        for seq in &msa.sequences { assert_eq!(seq.len(), w); }
    }

    #[test]
    fn engine_with_refinement() {
        let engine = MafftEngine::new(AlignmentMode::FftNsi { iterations: 5 });
        let msa = engine.align(&make_test_input());
        assert_eq!(msa.nseq(), 3);
        let w = msa.width();
        for seq in &msa.sequences { assert_eq!(seq.len(), w); }
    }

    #[test]
    fn engine_two_sequences() {
        let input = SequenceSet {
            sequences: vec![
                Sequence { name: "a".into(), data: b"ACDEFGHIK".to_vec() },
                Sequence { name: "b".into(), data: b"ACDEFGHIK".to_vec() },
            ],
            seq_type: SeqType::Protein,
        };
        let engine = MafftEngine::default();
        let msa = engine.align(&input);
        assert_eq!(msa.nseq(), 2);
        assert_eq!(msa.sequences[0], msa.sequences[1]);
    }

    #[test]
    fn engine_linsi_mode() {
        let engine = MafftEngine::new(AlignmentMode::LInsi { iterations: 2 });
        let msa = engine.align(&make_test_input());
        assert_eq!(msa.nseq(), 3);
        let w = msa.width();
        for seq in &msa.sequences { assert_eq!(seq.len(), w); }
    }
}
