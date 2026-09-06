# Python API

The full surface of the `pymafft` package, generated from the docstrings in
`crates/pymafft/python/pymafft/__init__.py` and the PyO3 classes in
`crates/pymafft/src/lib.rs`. Four entry-point functions, two result types,
one exception, and a stable iteration protocol.

## One flag layer, two front ends

Every pymafft function builds a `mafft-rs` command line from its keyword
arguments and runs it in-process through the CLI's own flag layer
(`mafft_rs::run_from_seqs`). Nothing about a flag's meaning is re-derived in
Python — `--auto`'s size heuristic, `--adjustdirection`'s strand detection,
`--nuc`'s type forcing and case fold, `--thread`'s pool, `--reorder`, the
scoring matrices and gap penalties are all decided by the same code the
binary runs — so for the same input and flags the Python result is
byte-identical to the CLI's stdout:

```python
import pymafft

# == `mafft-rs --auto --adjustdirection --thread 1 --nuc --quiet genes.fa`
result = pymafft.align(
    [("gene_1", "ATGGCTAGCTTGGACC"), ("gene_2", "GGTCCAAGCTAGCCAT")],
    strategy="auto", adjust_direction=True, threads=1, seq_type="nuc",
)
```

That is the flag set [panaroo](https://github.com/gtonkinhill/panaroo)
passes for its gene-cluster alignments; with pymafft it runs over an
in-memory list, with no temporary FASTA files and no FASTA re-parsing.

| Keyword | CLI flag | Values |
|---------|----------|--------|
| `strategy` | — / `--maxiterate N` / `--globalpair` / `--localpair` / `--genafpair` / `--auto` | `"fftns2"` (default), `"fftnsi"`, `"ginsi"`, `"linsi"`, `"einsi"`, `"auto"` |
| `maxiterate` | `--maxiterate N` | `0` = strategy default (100 for `fftnsi`, 1000 for `ginsi` / `linsi` / `einsi`); ignored by `fftns2` |
| `seq_type` | `--nuc` / `--amino` | `None` (detect from the residues, as the CLI does), `"nuc"`, `"amino"` |
| `adjust_direction` | `--adjustdirection` / `--adjustdirectionaccurately` | `False`, `True`, `"accurately"` |
| `threads` | `--thread N` | `0` = all cores |
| `reorder` | `--reorder` | `bool` |
| `retree` | `--retree N` | `None` = CLI default (2) |
| `scoring` | `--bl N` / `--jtt N` / `--tm N` | `"bl30"` … `"bl80"`, `"jtt200"`, `"tm100"`, … |
| `gap_open` / `gap_extend` | `--op X` / `--ep X` | floats (CLI defaults 1.53 / 0.123) |
| `quiet` | `--quiet` | default `True`; ignored when `progress` is given |
| `progress` | (stderr) | callable receiving each progress line |

Flags without a keyword go through [`pymafft.run`](#pymafftrun), which takes
the flag list verbatim.

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

### `pymafft.run`

::: pymafft.run
    options:
      heading_level: 4

## Result types

### `AlignmentResult`

The object returned by every alignment function. It's iterable, length-
addressable, indexable, and convertible to FASTA / tuples / Biopython.

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
        - to_biopython

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

## Errors

### `MafftError`

::: pymafft.MafftError
    options:
      heading_level: 4
      show_root_heading: false

Everything the command line would exit on is raised as `pymafft.MafftError`,
with the CLI's exit code in `.code` and its stderr text in `.message`
(also `str(err)`). It subclasses `ValueError`, so existing
`except ValueError` handlers keep working.

```python
try:
    pymafft.align([("p1", "MKTAYUAKQR"), ("p2", "MKTAYIAKQR")])
except pymafft.MafftError as e:
    print(e.code, e.message)   # 1 Illegal character U
```

With `quiet=False` the message carries the same multi-line banner the CLI
prints (`=== Alphabet 'U' is unknown. …`). Residues outside the alphabet can
be kept with `--anysymbol` through `run`.

Errors that are not the CLI's — an empty sequence list, an unsupported input
shape, an unreadable file, an invalid keyword value — stay plain
`ValueError` / `TypeError`.

## Progress

The CLI reports what it is doing on stderr (`mafft-rs v0.2.0`,
`36 sequences (aa), strategy: FFT-NS-2`, `Alignment: 717 columns`).
pymafft is quiet by default; pass a callable to receive exactly those lines,
one call per line without the trailing newline:

```python
lines = []
pymafft.align(seqs, strategy="auto", progress=lines.append)
# or straight into a logger
pymafft.align(seqs, progress=logging.getLogger("mafft").info)
```

`quiet=False` without a callable writes them to the process's stderr, as
the command line does (file descriptor 2, not `sys.stderr`). An exception
raised inside the callable is re-raised after the run returns.

The alignment runs with the GIL released, so several Python threads can
align concurrently; `threads` controls the Rust-side pool of each call.

## Input shapes

`pymafft.align` and `pymafft.run` accept three shapes interchangeably:

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

Residues are handled exactly as the CLI's FASTA reader handles them: the
sequence type is detected from the ATGC frequency unless `seq_type` forces
it, nucleotide residues are lowercased and protein residues uppercased
(C MAFFT's convention, so the output case matches the CLI's), `*` becomes
`-`, and a residue outside the alphabet raises `MafftError`.

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
result.to_fasta()        # FASTA-formatted string (one line per row)
result.to_tuples()       # [(name, gapped_seq), ...]
result.to_biopython()    # Bio.Align.MultipleSeqAlignment (lazy import)
result.width             # alignment width (uniform across rows)
result.score             # final SP score
```

## Strategy parameter

```python
pymafft.align(seqs, strategy="fftns2", maxiterate=0)   # default, == `mafft-rs`
pymafft.align(seqs, strategy="fftnsi", maxiterate=1000) # == --maxiterate 1000
pymafft.align(seqs, strategy="linsi",  maxiterate=1000) # == --localpair --maxiterate 1000
pymafft.align(seqs, strategy="ginsi",  maxiterate=1000) # == --globalpair --maxiterate 1000
pymafft.align(seqs, strategy="einsi",  maxiterate=1000) # == --genafpair --maxiterate 1000
pymafft.align(seqs, strategy="auto")                    # == --auto
```

See [Choosing a strategy](../getting-started/choosing-a-strategy.md)
for the trade-offs.

## CLI from Python

`pip install pymafft` also drops a `mafft-rs` console script on `$PATH`
(via the wheel-bundled native binary). Two equivalent invocations
after install:

```sh
mafft-rs --help
python -m pymafft --help
```

See the [CLI reference](../cli-reference.md) for the flag list. Any flag
listed there can also be passed in-process through `pymafft.run`.
