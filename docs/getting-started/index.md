# Getting Started

Three short pages that take you from zero to an aligned set of
sequences in your preferred environment.

The short version, for Python and the bundled CLI (install, then check
it works):

```sh
pip install pymafft
python -c "import pymafft; print(pymafft.align(['ACDEFGHIK', 'ACDEFHIK', 'ACDHIK']).to_fasta())"
mafft-rs --version
```

In more detail:

1. [**Install**](install.md) — choose your channel (Rust, Python, CLI
   binary), get it running, and
   [verify the install](install.md#verify-your-installation) in seconds.
2. [**First alignment**](first-alignment.md) — five-line examples for
   each channel.
3. [**Choosing a strategy**](choosing-a-strategy.md) — when to use
   FFT-NS-2 vs L-INS-i vs G-INS-i vs the others.
