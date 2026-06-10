# rust-MAFFT manuscript

Orphan branch for the rust-MAFFT *Bioinformatics* Application Note
manuscript. **No code lives here**, no shared history with `main` — it's
a deliberately separate timeline for the paper itself, its LaTeX
sources, figures, response-to-reviewers letters, and any auxiliary
artifacts (e.g. the dependency-graph figure, benchmark tables).

The companion code branch is [`main`](https://github.com/luksgrin/rust-MAFFT/tree/main).
Releases that this manuscript cites are pinned by Zenodo concept DOI
recorded in the `main` branch's [`CITATION.cff`](https://github.com/luksgrin/rust-MAFFT/blob/main/CITATION.cff).

## Why an orphan branch

A paper-as-code-in-the-software-repo pattern keeps the manuscript
near its software without polluting `main`'s history with LaTeX
files, build artifacts, or response-letter drafts. Reviewers and
readers who clone `main` get only what they need to use the
software; those interested in the paper itself can check out
`manuscript`.

## Layout (TBD as work progresses)

```
manuscript/
├── main.tex                  application-note source
├── refs.bib                  bibliography
├── figures/
│   └── workspace-deps.svg    dependency graph (from docs/architecture/workspace.md)
├── tables/
│   └── benchmarks.tex        rust-MAFFT vs C MAFFT timing comparison
├── cover-letter.tex
├── response-letter.tex       added at first round of revisions
└── README.md                 this file
```

## Working with this branch

```sh
# fetch the branch
git fetch origin manuscript

# check it out into a separate worktree (recommended — keeps main intact)
git worktree add ../rust-MAFFT-manuscript manuscript

# work, commit, push as normal
cd ../rust-MAFFT-manuscript
# ... edit, commit, push
```

## See also

- [`main`](https://github.com/luksgrin/rust-MAFFT) — the software itself
- [`archive/pre-release`](https://github.com/luksgrin/rust-MAFFT/tree/archive/pre-release) — 187-commit engineering history of the port
- [Project docs](https://luksgrin.github.io/rust-MAFFT/)
- [Citation guide](https://luksgrin.github.io/rust-MAFFT/citation/)
