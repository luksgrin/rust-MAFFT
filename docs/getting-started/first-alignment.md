# First alignment

The shortest viable program in each language.

## Rust

```rust
use mafft::{MafftEngine, AlignmentMode, read_fasta};

fn main() {
    let input = read_fasta("sequences.fasta").unwrap();
    let msa = MafftEngine::new(AlignmentMode::FftNs2).align(&input);
    for (name, row) in msa.names.iter().zip(msa.sequences.iter()) {
        println!(">{name}");
        println!("{}", std::str::from_utf8(row).unwrap());
    }
}
```

For finer control over the engine — custom scoring, guide tree
construction, iterative refinement — see the [`mafft` crate on
docs.rs](https://docs.rs/mafft).

## Python

```python
import pymafft

# 1) Pass raw strings (auto-named seq_1, seq_2, ...)
result = pymafft.align(["ACDEFGHIK", "ACDEFHIK", "ACDHIK"])

# 2) Or (name, seq) tuples
result = pymafft.align([
    ("alpha", "ACDEFGHIK"),
    ("beta",  "ACDEFHIK"),
    ("gamma", "ACDHIK"),
])

# 3) Or Biopython SeqRecord objects
from Bio import SeqIO
records = list(SeqIO.parse("sequences.fasta", "fasta"))
result = pymafft.align(records, strategy="linsi", maxiterate=1000)

# Print FASTA
print(result.to_fasta())

# Or iterate
for aligned in result:
    print(f">{aligned.name}")
    print(aligned.sequence)
```

See the [Python API reference](../python-api/index.md) for the full
surface area and the [Biopython interop guide](../python-api/biopython.md)
for round-tripping through `SeqRecord` / `MultipleSeqAlignment`.

## CLI

```sh
# FFT-NS-2 (default — fast progressive alignment)
mafft-rs sequences.fasta > aligned.fasta

# L-INS-i (slower; most accurate for < 200 sequences)
mafft-rs --localpair --maxiterate 1000 sequences.fasta > aligned.fasta

# Suppress progress output
mafft-rs --quiet sequences.fasta > aligned.fasta
```

The CLI is a drop-in replacement for `mafft` from the original C
toolchain — the only difference is the executable name. See the
[CLI reference](../cli-reference.md) for the full flag list.
