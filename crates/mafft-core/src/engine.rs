/// High-level MAFFT alignment engine.

use rayon::prelude::*;

use mafft_io::read_fasta;
use mafft_scoring::build_context;
use mafft_tree::{DistanceMatrix, musclesupg, ClusterMethod, pairwise_identity_distance, ktuple_distance};
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
    /// Number of guide tree rebuilds. C's FFT-NS-2 default is 2.
    pub retree: usize,
    /// Gap opening penalty override (positive float, e.g. 1.53 → internal -1530).
    /// None = use default for the scoring model.
    pub gap_open: Option<f64>,
    /// Offset/extension penalty override (positive float, e.g. 0.123 → internal -123).
    /// None = use default.
    pub gap_offset: Option<f64>,
    /// Disable FFT: force pure DP for all alignment steps.
    pub nofft: bool,
}

impl Default for MafftEngine {
    fn default() -> Self {
        Self {
            mode: AlignmentMode::FftNs2,
            scoring_model: ScoringModel::Jtt,
            retree: 2,
            gap_open: None,
            gap_offset: None,
            nofft: false,
        }
    }
}

impl MafftEngine {
    pub fn new(mode: AlignmentMode) -> Self {
        Self { mode, scoring_model: ScoringModel::Jtt, retree: 2, gap_open: None, gap_offset: None, nofft: false }
    }

    /// Set the number of guide tree rebuilds.
    pub fn with_retree(mut self, retree: usize) -> Self {
        self.retree = retree;
        self
    }

    /// Set gap opening penalty (positive float, e.g. 1.53).
    pub fn with_gap_open(mut self, op: f64) -> Self {
        self.gap_open = Some(op);
        self
    }

    /// Set offset/extension penalty (positive float, e.g. 0.123).
    pub fn with_gap_offset(mut self, ep: f64) -> Self {
        self.gap_offset = Some(ep);
        self
    }

    /// Disable FFT: force pure DP for all alignment steps.
    pub fn with_nofft(mut self, nofft: bool) -> Self {
        self.nofft = nofft;
        self
    }

    /// Set scoring model (e.g. BLOSUM with specific number).
    pub fn with_scoring_model(mut self, model: ScoringModel) -> Self {
        self.scoring_model = model;
        self
    }

    /// Align a set of sequences.
    pub fn align(&self, input: &SequenceSet) -> MultipleAlignment {
        let seq_type = input.seq_type;
        let scoring_model = if seq_type.is_nucleotide() {
            ScoringModel::Dna
        } else {
            self.scoring_model
        };

        let mut scoring = build_context(scoring_model, seq_type);

        // Apply gap penalty overrides if set.
        // C convention: --op 1.53 means ppenalty = -1530 (multiply by -1000).
        // After scaling: penalty = (int)(600/1000 * ppenalty + 0.5).
        if let Some(op) = self.gap_open {
            let ppenalty = -(op * 1000.0) as i32;
            let scale = if seq_type.is_nucleotide() { 3.0 * 600.0 / 1000.0 } else { 600.0 / 1000.0 };
            scoring.gap.open = (scale * ppenalty as f64 + 0.5) as i32;
        }
        if let Some(ep) = self.gap_offset {
            let poffset = -(ep * 1000.0) as i32;
            let scale = 600.0 / 1000.0;
            scoring.gap.offset = (scale * poffset as f64 + 0.5) as i32;
        }

        let nseq = input.nseq();
        let sequences: Vec<Vec<u8>> = input.sequences.iter().map(|s| s.data.clone()).collect();
        let names: Vec<String> = input.sequences.iter().map(|s| s.name.clone()).collect();

        let use_fft = !self.nofft && matches!(
            self.mode,
            AlignmentMode::FftNs2 | AlignmentMode::FftNsi { .. }
        );

        // Step 1: Initial pairwise distances (parallel)
        let mut dm = compute_distance_matrix_from_seqs(&sequences);

        // Step 2: Build guide tree and progressive align, repeating `retree` times.
        // Each iteration after the first computes distances from the ALIGNMENT
        // (not the raw sequences), producing a better tree.
        let retree = self.retree.max(1);
        let mut msa = MultipleAlignment {
            sequences: sequences.clone(),
            names: names.clone(),
            score: 0.0,
        };

        for pass in 0..retree {
            // Build guide tree from current distance matrix
            let topo = musclesupg(&dm, ClusterMethod::default());

            // Progressive alignment
            let input_seqs = if pass == 0 {
                // First pass: use raw sequences
                sequences.clone()
            } else {
                // Subsequent passes: strip gaps from previous alignment
                // to get updated unaligned sequences (same content, different
                // order may produce better tree-guided alignment)
                sequences.clone()
            };

            msa = progressive_align(&input_seqs, &names, &topo, &scoring, use_fft);

            // If there's another pass, compute new distances from the alignment
            if pass + 1 < retree {
                dm = compute_distance_matrix_from_alignment(&msa.sequences);
            }
        }

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

        // Step 4: Iterative refinement (if mode requires it)
        match &self.mode {
            AlignmentMode::FftNs2 => {}
            AlignmentMode::FftNsi { iterations }
            | AlignmentMode::GInsi { iterations }
            | AlignmentMode::LInsi { iterations }
            | AlignmentMode::EInsi { iterations } => {
                // Rebuild tree one more time for refinement
                let dm = compute_distance_matrix_from_alignment(&msa.sequences);
                let topo = musclesupg(&dm, ClusterMethod::default());
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

/// Compute pairwise 6-tuple distances from raw (unaligned) sequences.
/// Matches C's default distance computation using `commonsextet_p`.
fn compute_distance_matrix_from_seqs(sequences: &[Vec<u8>]) -> DistanceMatrix {
    let nseq = sequences.len();
    let pairs: Vec<(usize, usize, f64)> = (0..nseq)
        .into_par_iter()
        .flat_map(|i| {
            let seqs = sequences;
            ((i + 1)..nseq).into_par_iter().map(move |j| {
                let d = ktuple_distance(&seqs[i], &seqs[j], 6);
                (i, j, d)
            })
        })
        .collect();

    let mut dm = DistanceMatrix::new(nseq);
    for (i, j, d) in pairs {
        dm.set(i, j, d);
    }
    dm
}

/// Compute pairwise identity distances from aligned sequences (with gaps).
fn compute_distance_matrix_from_alignment(sequences: &[Vec<u8>]) -> DistanceMatrix {
    let nseq = sequences.len();
    let pairs: Vec<(usize, usize, f64)> = (0..nseq)
        .into_par_iter()
        .flat_map(|i| {
            let seqs = sequences;
            ((i + 1)..nseq).into_par_iter().map(move |j| {
                let d = pairwise_identity_distance(&seqs[i], &seqs[j]);
                (i, j, d)
            })
        })
        .collect();

    let mut dm = DistanceMatrix::new(nseq);
    for (i, j, d) in pairs {
        dm.set(i, j, d);
    }
    dm
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

    #[test]
    fn engine_retree_1_vs_2() {
        let input = make_test_input();
        let msa1 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(1).align(&input);
        let msa2 = MafftEngine::new(AlignmentMode::FftNs2).with_retree(2).align(&input);
        // Both should produce valid alignments
        assert_eq!(msa1.nseq(), 3);
        assert_eq!(msa2.nseq(), 3);
        let w1 = msa1.width();
        let w2 = msa2.width();
        for seq in &msa1.sequences { assert_eq!(seq.len(), w1); }
        for seq in &msa2.sequences { assert_eq!(seq.len(), w2); }
    }
}
