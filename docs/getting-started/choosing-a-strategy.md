# Choosing a strategy

MAFFT's strategy choice trades **speed** against **accuracy**. The right
pick depends on (i) how many sequences you have, (ii) how diverged they
are, and (iii) whether you suspect large insertions or alignable regions
embedded in unrelated flanks.

## Quick decision table

| Sequences | Divergence | Insertions / unalignable flanks | Strategy | CLI flag |
|--|--|--|--|--|
| any | any | no | **FFT-NS-2** (default) | (none — default) |
| many, fast iter | any | no | **FFT-NS-i** | `--maxiterate 1000` |
| < 200 | any | no | **L-INS-i** | `--localpair --maxiterate 1000` |
| < 200 | similar lengths | no | **G-INS-i** | `--globalpair --maxiterate 1000` |
| < 200 | mixed | yes | **E-INS-i** | `--genafpair --maxiterate 1000` |
| 10 000+ | any | no | **PartTree** | `--parttree` |

## Cheat sheet

=== "Rust"

    ```rust
    use mafft::AlignmentMode;

    // Default
    AlignmentMode::FftNs2

    // Iterative refinement
    AlignmentMode::FftNsi { iterations: 1000 }

    // Most accurate (< 200 seqs)
    AlignmentMode::LInsi  { iterations: 1000 }   // local pair
    AlignmentMode::GInsi  { iterations: 1000 }   // global pair
    AlignmentMode::EInsi  { iterations: 1000 }   // generalized affine
    ```

=== "Python"

    ```python
    pymafft.align(seqs)                              # FFT-NS-2
    pymafft.align(seqs, strategy="fftnsi", maxiterate=1000)
    pymafft.align(seqs, strategy="linsi", maxiterate=1000)
    pymafft.align(seqs, strategy="ginsi", maxiterate=1000)
    pymafft.align(seqs, strategy="einsi", maxiterate=1000)
    ```

=== "CLI"

    ```sh
    mafft-rs in.fa                                      # FFT-NS-2
    mafft-rs --maxiterate 1000 in.fa                    # FFT-NS-i
    mafft-rs --localpair  --maxiterate 1000 in.fa       # L-INS-i
    mafft-rs --globalpair --maxiterate 1000 in.fa       # G-INS-i
    mafft-rs --genafpair  --maxiterate 1000 in.fa       # E-INS-i
    ```

## Pragmatic guidance

- **You don't know which to pick** → start with **FFT-NS-2**. It's the
  default for a reason; on most inputs it's within 1–2 SP-score points
  of the iterative methods.
- **You have time and < 200 sequences** → run **L-INS-i**. Highest
  average accuracy on BAliBASE.
- **Your sequences clearly differ in length but the alignable region is
  obvious** → **L-INS-i** (local pairwise reduces noise from flanks).
- **You have 10 000+ sequences** → **PartTree** is the only practical
  option; `O(N log N)` distance estimation instead of `O(N²)`.

The MAFFT [official docs](https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)
go into more depth on the algorithmic differences. rust-MAFFT is
byte-identical to the C reference for every strategy listed above, so
guidance written for upstream MAFFT applies verbatim here.
