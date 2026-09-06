# Architecture

Three pages on how rust-MAFFT is put together and the design decisions
that fall out of the byte-identity constraint.

1. [**Workspace overview**](workspace.md) — the 12 crates, who depends
   on whom, and what gets published.
2. [**Byte-identity**](byte-identity.md) — why we held to it, and the
   classes of bugs we hit along the way (FP order, tie-breaks,
   FMA-aware compilation).
3. [**Cross-validation**](cross-validation.md) — how 99 FFI test functions
   compare each Rust function to its C counterpart at the value level,
   not just the alignment level.

## What's still out of scope

A short, honest list of things rust-MAFFT does NOT do, or does not
match C on (as of 2026-09-06; no other divergence is known):

- **Modes that need third-party binaries** — `--xinsi` (CONTRAfold),
  `--qinsi` (`mxscarnamod`), `--scarnalike` (`dash_client`). Wired
  exactly as in C; validated only where the binary is present.
- **Structure-aware alignment** (`--pdbidlist`, `--pdbfilelist`) —
  upstream C MAFFT disabled these in Dec 2018. rust-MAFFT matches the
  upstream "temporarily unavailable" exit verbatim. A future port via
  a Rust-native structural-aligner dependency (e.g. one of the TMalign
  ports) is plausible but not on the roadmap.
- **Out of scope by design** — RNA-structure pipelines (DAFS, FoldAlign,
  LARA, SCARNA, …), external pairwise aligners (`--blastpair`,
  `--lastpair`, …), `--mpi`.
- **Two tied-trace residuals** — `--exp ≥ 4.30` on FFT-NS-2 without
  refinement, and a trailing zero-scoring residue under `--add`; same
  width and score, one residue shifted. Neither is reachable from a
  realistic command line.
- **`--thread N ≥ 2`** — C MAFFT itself is nondeterministic there, so
  there is no single C output to match. rust-MAFFT is deterministic at
  every thread count and reproduces C's `--thread 1` (`athread`) path.

The full status matrix lives in
[TODO.md](https://github.com/luksgrin/rust-MAFFT/blob/main/TODO.md) on
the repo.
