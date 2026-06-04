# Architecture

Three pages on how rust-MAFFT is put together and the design decisions
that fall out of the byte-identity constraint.

1. [**Workspace overview**](workspace.md) — the 12 crates, who depends
   on whom, and what gets published.
2. [**Byte-identity**](byte-identity.md) — why we held to it, and the
   classes of bugs we hit along the way (FP order, tie-breaks,
   FMA-aware compilation).
3. [**Cross-validation**](cross-validation.md) — how ~140 FFI tests
   compare each Rust function to its C counterpart at the value level,
   not just the alignment level.

## What's still out of scope

A short, honest list of things rust-MAFFT does NOT do yet:

- **`--adjustdirectionaccurately`** — the DP-based variant (the
  k-mer-based `--adjustdirection` IS byte-identical).
- **Structure-aware alignment** (`--pdbidlist`, `--pdbfilelist`) —
  upstream C MAFFT disabled these in Dec 2018. rust-MAFFT matches the
  upstream "temporarily unavailable" exit verbatim. A future port via
  a Rust-native structural-aligner dependency (e.g. one of the TMalign
  ports) is plausible but not on the roadmap.
- **Real-time progress in pymafft** — the Python API runs to
  completion; there's no streaming progress callback yet.

The full status matrix lives in
[TODO.md](https://github.com/luksgrin/rust-MAFFT/blob/main/TODO.md) on
the repo.
