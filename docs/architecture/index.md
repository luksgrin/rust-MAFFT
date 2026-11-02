# Architecture

Five pages on how rust-MAFFT is put together, the design decisions
that fall out of the byte-identity constraint, and how it follows
upstream MAFFT.

1. [**Workspace overview**](workspace.md) — the 12 crates, who depends
   on whom, and what gets published.
2. [**Byte-identity**](byte-identity.md) — why we held to it, and the
   classes of bugs we hit along the way (FP order, tie-breaks,
   FMA-aware compilation).
3. [**Cross-validation**](cross-validation.md) — how 99 FFI test functions
   compare each Rust function to its C counterpart at the value level,
   not just the alignment level.
4. [**Upstream policy**](upstream-policy.md) — how the `mafft-upstream`
   pin is recorded and bumped, which way bug fixes flow, and what the
   port commits not to change.
5. [**Memory safety**](memory-safety.md) — the C reference under
   AddressSanitizer / UndefinedBehaviorSanitizer and `mafft-rs` on the
   same corpora and on malformed inputs: what each side reported, what
   was fixed, what is open.

## What's still out of scope

A short list of things rust-MAFFT does NOT do, or does not match C on.
The full list, with the inputs, is the README's
[Known limitations](https://github.com/luksgrin/rust-MAFFT#known-limitations):

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
  width, one residue shifted, so the sum-of-pairs score differs slightly.
  Neither is reachable from a realistic command line.
- **Partly implemented options** — `--add --localpair` / `--globalpair`
  / `--genafpair` are not implemented; `--add --maxiterate N` (N ≥ 2)
  can differ in the refinement; `--groupsize` is accepted and ignored;
  `--partsize` and `--6merpair` are not accepted; `--maxiterate N
  --memsave` does not follow `dvtditr -M`; `--nofft --memsave` differs on
  one CDS set
  ([E](https://github.com/luksgrin/rust-MAFFT#e-other-known-differences)).
- **Inputs outside the parity corpora** — the inputs from the
  memory-safety run now match C; names longer than about 250 bytes are
  kept whole where C truncates them (intentional).
- **`--thread N ≥ 2`** — C MAFFT itself is nondeterministic there, so
  there is no single C output to match. rust-MAFFT is deterministic at
  every thread count and reproduces C's `--thread 1` (`athread`) path.

The full status matrix lives in
[TODO.md](https://github.com/luksgrin/rust-MAFFT/blob/main/TODO.md) on
the repo.
