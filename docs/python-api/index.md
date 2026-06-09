# Python API

The full surface of the `pymafft` package, auto-generated from the
PyO3 docstrings in `crates/pymafft/src/lib.rs`. Three entry-point
functions, two result types, and a stable iteration protocol.

## Functions

### `pymafft.align`

::: pymafft.align
    options:
      heading_level: 4

### `pymafft.align_file`

::: pymafft.align_file
    options:
      heading_level: 4

### `pymafft.align_fasta_string`

::: pymafft.align_fasta_string
    options:
      heading_level: 4

## Result types

### `AlignmentResult`

The object returned by every alignment function. It's iterable, length-
addressable, indexable, and convertible to FASTA / tuples.

::: pymafft.AlignmentResult
    options:
      heading_level: 4
      show_root_heading: false
      members:
        - sequences
        - width
        - score
        - nseq
        - to_fasta
        - to_tuples

### `AlignedSequence`

A single aligned row.

::: pymafft.AlignedSequence
    options:
      heading_level: 4
      show_root_heading: false
      members:
        - name
        - sequence
        - ungapped

## Input shapes

`pymafft.align` accepts three shapes interchangeably:

```python
# 1. Bare strings — auto-named seq_1, seq_2, ...
pymafft.align(["ACDEFGHIK", "ACDEFHIK"])

# 2. (name, seq) tuples
pymafft.align([("alpha", "ACDEFGHIK"), ("beta", "ACDEFHIK")])

# 3. Any object with `.id` and `.seq` attributes — duck-types
#    Biopython's SeqRecord, scikit-bio's Sequence, etc.
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
records = [
    SeqRecord(Seq("ACDEFGHIK"), id="alpha"),
    SeqRecord(Seq("ACDEFHIK"),  id="beta"),
]
pymafft.align(records)
```

For full Biopython interop including `MultipleSeqAlignment` round-trips,
see [Biopython interop](biopython.md).

## Output shapes

```python
result = pymafft.align(seqs)

# Iterate
for aligned in result:
    print(aligned.name, aligned.sequence)

# Index
result[0].sequence       # first aligned row
len(result)              # number of sequences

# Aggregate views
result.to_fasta()        # FASTA-formatted string
result.to_tuples()       # [(name, gapped_seq), ...]
result.width             # alignment width (uniform across rows)
result.score             # final SP score
```

## Strategy parameter

```python
pymafft.align(seqs, strategy="fftns2", maxiterate=0)   # default
pymafft.align(seqs, strategy="fftnsi", maxiterate=1000)
pymafft.align(seqs, strategy="linsi",  maxiterate=1000)
pymafft.align(seqs, strategy="ginsi",  maxiterate=1000)
pymafft.align(seqs, strategy="einsi",  maxiterate=1000)
```

See [Choosing a strategy](../getting-started/choosing-a-strategy.md)
for the trade-offs.

## CLI from Python

`pip install pymafft` also drops a `mafft-rs` console script on `$PATH`
(via the wheel-bundled native binary). Three equivalent invocations
after install:

```sh
mafft-rs --help
python -m pymafft --help
```

See the [CLI reference](../cli-reference.md) for the flag list.
