# Biopython interop

`pymafft` duck-types `Bio.SeqRecord.SeqRecord` — anything with `.id` and
`.seq` attributes works as input. Biopython itself is not a runtime
dependency. The cooperation is one-way at the input boundary
(`pymafft.align` reads `SeqRecord`s) and two-way at the output boundary
(`AlignmentResult` produces FASTA that re-parses cleanly).

## SeqRecord → pymafft

```python
from Bio import SeqIO
import pymafft

records = list(SeqIO.parse("input.fasta", "fasta"))
result  = pymafft.align(records, strategy="linsi", maxiterate=1000)
```

**The id mapping**: `SeqRecord.id` becomes the alignment row name —
NOT `.name` (which is set to `.id` if not provided) and NOT
`.description` (which includes the post-id header text). This matches
how Biopython itself treats `id` as the stable identifier.

## pymafft → MultipleSeqAlignment

Two routes — pick by style preference, the result is identical.

=== "Via FASTA string"

    ```python
    from io import StringIO
    from Bio import AlignIO

    msa = AlignIO.read(StringIO(result.to_fasta()), "fasta")
    assert msa.get_alignment_length() == result.width
    ```

=== "Direct construction"

    ```python
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    msa = MultipleSeqAlignment([
        SeqRecord(Seq(s.sequence), id=s.name, description="")
        for s in result.sequences
    ])
    ```

## Full round-trip

```python
from Bio import SeqIO, AlignIO
from io import StringIO
import pymafft

# Read → align → write
records = list(SeqIO.parse("opsins.fasta", "fasta"))
result  = pymafft.align(records, strategy="linsi", maxiterate=1000)

msa = AlignIO.read(StringIO(result.to_fasta()), "fasta")
AlignIO.write([msa], "opsins.aligned.aln", "clustal")
```

## Header caveat

`SeqIO.write(records, ..., "fasta")` emits headers of the form
`id description`. If you pass a SeqRecord directly to `pymafft.align`,
the row name is just `id`. If you serialize the records to FASTA first
and then call `pymafft.align_fasta_string(...)`, the row name is the
whole header (`id description`).

When this matters, prefer the direct `pymafft.align(records)` path — the
test suite verifies both produce equivalent aligned sequences, but the
names differ.

## Tests

The interop guarantees above are exercised by 17 tests in
`crates/pymafft/tests/test_biopython.py`. The suite skips cleanly when
Biopython isn't installed (via `pytest.importorskip("Bio")`).
