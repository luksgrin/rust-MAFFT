# Citing rust-MAFFT

rust-MAFFT is a port. Every algorithmic decision in this codebase
traces back to work by **Kazutaka Katoh** and colleagues at CBRC.
The FFT-anchored alignment, the iterative-refinement strategies
(L-INS-i / G-INS-i / E-INS-i), the scoring matrices, the PartTree
heuristic for large datasets — all of it is theirs. This project is
engineering on top of their science.

If you use rust-MAFFT in published work, **you must cite the original
MAFFT paper**. The rust-MAFFT manuscript (in preparation for
*Bioinformatics*) should be cited additionally once it is published —
but the algorithmic foundation always remains the Katoh et al. work,
and that citation never goes away.

## Primary citation — MAFFT v7

**Always required**, regardless of which mode you use.

!!! quote ""

    Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence
    alignment software version 7: improvements in performance and
    usability. *Molecular Biology and Evolution*, 30(4), 772–780.
    [10.1093/molbev/mst010](https://doi.org/10.1093/molbev/mst010)

=== "BibTeX"

    ```bibtex
    @article{Katoh2013,
      author  = {Katoh, Kazutaka and Standley, Daron M.},
      title   = {{MAFFT} Multiple Sequence Alignment Software Version 7:
                 Improvements in Performance and Usability},
      journal = {Molecular Biology and Evolution},
      volume  = {30},
      number  = {4},
      pages   = {772--780},
      year    = {2013},
      doi     = {10.1093/molbev/mst010}
    }
    ```

=== "RIS"

    ```
    TY  - JOUR
    AU  - Katoh, Kazutaka
    AU  - Standley, Daron M.
    PY  - 2013
    TI  - MAFFT Multiple Sequence Alignment Software Version 7: Improvements in Performance and Usability
    JO  - Molecular Biology and Evolution
    VL  - 30
    IS  - 4
    SP  - 772
    EP  - 780
    DO  - 10.1093/molbev/mst010
    ER  -
    ```

## Foundational citation — original MAFFT (2002)

Cite this **additionally** if your work uses the FFT-NS-1 or FFT-NS-2
modes. These are the original FFT-based strategies that the 2002 paper
introduced.

!!! quote ""

    Katoh, K., Misawa, K., Kuma, K., & Miyata, T. (2002). MAFFT: a
    novel method for rapid multiple sequence alignment based on fast
    Fourier transform. *Nucleic Acids Research*, 30(14), 3059–3066.
    [10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436)

=== "BibTeX"

    ```bibtex
    @article{Katoh2002,
      author  = {Katoh, Kazutaka and Misawa, Kazuharu and
                 Kuma, Kei-ichi and Miyata, Takashi},
      title   = {{MAFFT}: a novel method for rapid multiple sequence
                 alignment based on fast {F}ourier transform},
      journal = {Nucleic Acids Research},
      volume  = {30},
      number  = {14},
      pages   = {3059--3066},
      year    = {2002},
      doi     = {10.1093/nar/gkf436}
    }
    ```

## Mode-specific citations

Cite these **in addition** to the primary citation when the
corresponding mode is part of your analysis.

### Iterative refinement (L-INS-i / G-INS-i / E-INS-i)

Introduced in MAFFT v5.

!!! quote ""

    Katoh, K., Kuma, K., Toh, H., & Miyata, T. (2005). MAFFT version 5:
    improvement in accuracy of multiple sequence alignment. *Nucleic
    Acids Research*, 33(2), 511–518.
    [10.1093/nar/gki198](https://doi.org/10.1093/nar/gki198)

```bibtex
@article{Katoh2005,
  author  = {Katoh, Kazutaka and Kuma, Kei-ichi and
             Toh, Hiroyuki and Miyata, Takashi},
  title   = {{MAFFT} version 5: improvement in accuracy of multiple
             sequence alignment},
  journal = {Nucleic Acids Research},
  volume  = {33},
  number  = {2},
  pages   = {511--518},
  year    = {2005},
  doi     = {10.1093/nar/gki198}
}
```

### PartTree (`--parttree` / `--dpparttree`)

The fast guide-tree algorithm used for 10K+ sequence datasets.

!!! quote ""

    Katoh, K., & Toh, H. (2007). PartTree: an algorithm to build an
    approximate tree from a large number of unaligned sequences.
    *Bioinformatics*, 23(3), 372–374.
    [10.1093/bioinformatics/btl592](https://doi.org/10.1093/bioinformatics/btl592)

```bibtex
@article{Katoh2007,
  author  = {Katoh, Kazutaka and Toh, Hiroyuki},
  title   = {{PartTree}: an algorithm to build an approximate tree
             from a large number of unaligned sequences},
  journal = {Bioinformatics},
  volume  = {23},
  number  = {3},
  pages   = {372--374},
  year    = {2007},
  doi     = {10.1093/bioinformatics/btl592}
}
```

## CITATION.cff

For citation managers that support the [Citation File Format](https://citation-file-format.github.io/),
the canonical machine-readable record is at
[`CITATION.cff`](https://github.com/luksgrin/rust-MAFFT/blob/main/CITATION.cff)
in the repository root. GitHub also exposes a "Cite this repository"
button on the [project page](https://github.com/luksgrin/rust-MAFFT)
that reads this file.

You can also generate plain-text or BibTeX from the CFF directly:

```sh
pip install cffconvert
cffconvert -f bibtex -i CITATION.cff
cffconvert -f apalike -i CITATION.cff
```

## Citing the rust-MAFFT port

The rust-MAFFT manuscript is **in preparation for *Bioinformatics***.
Once published, the citation will appear here and `CITATION.cff` will
be updated. Until then, please:

1. **Always cite Katoh & Standley 2013** — it is the algorithmic
   foundation.
2. Optionally cite the rust-MAFFT GitHub release if you want to pin a
   reproducible version of the port (e.g. *Goiriz, rust-MAFFT v0.2.0,
   https://github.com/luksgrin/rust-MAFFT, 2026*). DOI-minting via
   Zenodo is on the roadmap.

## Why both should be cited

A port that produces byte-identical output to its reference does not
create new science — it engineers new accessibility, integrability,
and language-ecosystem reach for existing science. The science
belongs to the upstream authors; the engineering belongs to the port
authors. Citing both reflects this division honestly and is the
practice we ask everyone using rust-MAFFT to follow.
