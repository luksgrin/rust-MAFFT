#!/usr/bin/env python3
"""Print the peak-RSS summary of a `memory.py` TSV as Markdown tables.

    python3 benchmarks/summarize_memory.py benchmarks/results/memory_<date>.tsv [more.tsv ...]

Rows from later files replace earlier ones for the same
(input, mode, threads, impl) cell, so a re-measured cell can be layered on
top of the main run.
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from collections import defaultdict

MiB = 2**20


def main(paths: list[str]) -> int:
    rss = defaultdict(list)
    meta = {}
    host = ""
    for path in paths:
        seen = set()
        for r in csv.DictReader(open(path), delimiter="\t"):
            k = (r["input"], r["mode"], int(r["threads"]), r["impl"])
            if k not in seen:
                rss[k] = []
                seen.add(k)
            rss[k].append(int(r["peak_rss_bytes"]))
            meta[r["input"]] = (r["n"], r["median_len"])
            host = r["host"]

    modes = sorted({k[1] for k in rss}, key=lambda m: ("+" in m, m))
    print(f"host: {host}\n")
    print("Median peak RSS in MiB over the runs (max in parentheses when it differs by > 5 %); "
          "C / rs = ratio of medians.\n")
    print("| input (n, median len) | mode | threads | C MAFFT | mafft-rs | C / rs |")
    print("|---|---|---:|---:|---:|---:|")
    for name in sorted(meta, key=lambda n: (int(meta[n][0]) * float(meta[n][1]))):
        n, ml = meta[name]
        for mode in modes:
            for t in sorted({k[2] for k in rss if k[0] == name and k[1] == mode}):
                c, r = rss.get((name, mode, t, "c")), rss.get((name, mode, t, "rust"))

                def fmt(v):
                    if not v:
                        return "–"
                    med, mx = st.median(v), max(v)
                    s = f"{med / MiB:,.0f}" if med >= 10 * MiB else f"{med / MiB:.1f}"
                    if mx > 1.05 * med:
                        s += f" ({mx / MiB:,.0f})"
                    return s

                ratio = f"{st.median(c) / st.median(r):.2f}" if c and r else "–"
                print(f"| {name} ({n}, {float(ml):g}) | {mode} | {t} | {fmt(c)} | {fmt(r)} | {ratio} |")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
