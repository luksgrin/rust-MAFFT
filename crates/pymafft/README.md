# pymafft

Python bindings for [MAFFT-rs](https://github.com/luksgrin/rust-MAFFT), a Rust reimplementation of the [MAFFT](https://mafft.cbrc.jp/alignment/software/) multiple sequence alignment tool.

## Installation

```bash
pip install pymafft
```

## Quick start

```python
import pymafft

# Align protein sequences
result = pymafft.align([
    "ACDEFGHIKLMNPQRSTVWY",
    "ACDEFHIKLMNPQRSTVWY",
    "ACDEHIKLMNPQRSTVWY",
])

for seq in result:
    print(f"{seq.name}: {seq.sequence}")

print(f"Alignment width: {result.width}")
```

## Usage

### Align from a list of sequences

```python
# Auto-named (seq_1, seq_2, ...)
result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])

# Named sequences
result = pymafft.align([
    ("human", "ACDEFGHIKLMNP"),
    ("mouse", "ACDEFHIKLMNP"),
    ("fish",  "ACDEHIKLMNP"),
])
```

### Align from a FASTA file

```python
result = pymafft.align_file("sequences.fasta")
```

### Align from a FASTA string

```python
fasta = ">s1\nACDEFGHIK\n>s2\nACDEFHIK\n"
result = pymafft.align_fasta_string(fasta)
```

### Choose an alignment strategy

```python
# FFT-NS-2: fast default
result = pymafft.align(seqs, strategy="fftns2")

# G-INS-i: global iterative (accurate, for globally alignable sequences)
result = pymafft.align(seqs, strategy="ginsi")

# L-INS-i: local iterative (most accurate for <200 sequences)
result = pymafft.align(seqs, strategy="linsi")

# E-INS-i: generalized affine (for sequences with large internal gaps)
result = pymafft.align(seqs, strategy="einsi")

# FFT-NS-i: progressive + iterative refinement
result = pymafft.align(seqs, strategy="fftnsi", maxiterate=100)

# --auto: let MAFFT pick from the number and length of the sequences
result = pymafft.align(seqs, strategy="auto")
```

### CLI-equivalent options

Every keyword is one `mafft-rs` command-line flag, applied by the CLI's own
flag layer, so the result is byte-identical to the command line's for the
same flags.

```python
# `mafft --auto --adjustdirection --thread 1 --nuc` over in-memory records
result = pymafft.align(genes, strategy="auto", adjust_direction=True,
                       threads=1, seq_type="nuc")

# Substitution matrix and gap penalties (`--bl 50 --op 1.0 --ep 0.2`)
result = pymafft.align(seqs, scoring="bl50", gap_open=1.0, gap_extend=0.2)

# Guide-tree order output, one more tree rebuild
result = pymafft.align(seqs, reorder=True, retree=3)

# Progress lines (the CLI's stderr) to a callable; quiet by default
result = pymafft.align(seqs, strategy="auto", progress=print)

# Any other flag, verbatim
result = pymafft.run(["--nofft", "--maxiterate", "2", "--quiet"], seqs)
```

Failures raise `pymafft.MafftError` (a `ValueError`) with the CLI's exit
code in `.code` and stderr text in `.message`, e.g. `Illegal character U`
for a residue outside the alphabet. Nucleotide output is lowercase and
protein output uppercase, as with C MAFFT.

### Work with results

```python
result = pymafft.align(["ACDEFGHIK", "ACDEFHIK"])

# Access sequences
result.nseq              # number of sequences
result.width             # alignment width (columns)
result[0].name           # sequence name
result[0].sequence       # aligned sequence (with gaps)
result[0].ungapped()     # sequence without gaps

# Export
result.to_fasta()        # FASTA-formatted string
result.to_tuples()       # list of (name, sequence) tuples
result.to_biopython()    # Bio.Align.MultipleSeqAlignment (Biopython imported lazily)

# Iterate
for seq in result:
    print(seq.name, len(seq))
```

## Supported strategies

| Strategy | Flag | Description |
|----------|------|-------------|
| `fftns2` | default | Fast progressive alignment |
| `fftnsi` | `--maxiterate N` | Progressive + iterative refinement |
| `ginsi` | `--globalpair` | Global iterative (accurate) |
| `linsi` | `--localpair` | Local iterative (most accurate for <200 seqs) |
| `einsi` | `--genafpair` | Generalized affine (large internal gaps) |
| `auto` | `--auto` | Chosen from sequence count and length |

## License

- Rust code: MIT
- Original MAFFT algorithm: BSD-3-Clause

Based on [MAFFT](https://mafft.cbrc.jp/alignment/software/) by Kazutaka Katoh.
