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
```

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

## License

- Rust code: MIT
- Original MAFFT algorithm: BSD-3-Clause

Based on [MAFFT](https://mafft.cbrc.jp/alignment/software/) by Kazutaka Katoh.
