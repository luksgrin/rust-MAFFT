#!/usr/bin/env python3
"""`pymafft.align(..., threads=N)` from Python: determinism, timing, GIL release.

Needs `pymafft` importable (``maturin develop --release`` in crates/pymafft)
and the `mafft-rs` release binary (the CLI output is the reference).

1. For every input x strategy x thread count, `pymafft.align` is called
   `--runs` times; every result must equal, row for row, the
   `mafft-rs --thread 1` CLI output for the same flags (exit status 1 if
   not). Wall-clock medians are recorded.
2. GIL release: one call is timed alone, then two identical calls are
   started in two Python threads. If the GIL is released during the
   alignment, the pair takes about as long as one call, not two.

Output: benchmarks/results/pymafft_threads_<date>.tsv with columns
input, n, median_len, strategy, threads, run, wall_s, equals_cli_t1,
load1, host, fp_policy (GIL rows use strategy `gil:<strategy>` and
threads = number of concurrent calls).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as cm  # noqa: E402

try:
    import pymafft
except ImportError:  # pragma: no cover
    sys.exit("pymafft is not importable: run `maturin develop --release` in crates/pymafft first")

STRATEGY_MODE = {"fftnsi": "fftnsi", "linsi": "linsi", "ginsi": "ginsi"}
COLUMNS = ["input", "n", "median_len", "strategy", "threads", "run", "wall_s",
           "equals_cli_t1", "load1", "host", "fp_policy"]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            recs.append((line[1:].rstrip(), []))
        elif recs:
            recs[-1][1].append(line.strip())
    return [(n, "".join(s)) for n, s in recs]


def parse_fasta_bytes(data: bytes) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    for line in data.decode().splitlines():
        if line.startswith(">"):
            recs.append((line[1:], []))
        elif recs:
            recs[-1][1].append(line.strip())
    return [(n, "".join(s)) for n, s in recs]


def align(records, strategy: str, threads: int):
    return pymafft.align(records, strategy=strategy, maxiterate=1000, threads=threads)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="*", default=["sample36", "mtb_cds_120x1400.fa", "BB40049.tfa"])
    ap.add_argument("--strategies", nargs="*", default=["fftnsi", "linsi", "ginsi"], choices=list(STRATEGY_MODE))
    ap.add_argument("--threads", nargs="*", type=int, default=[1, 2, 4, 8, 10])
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--gil-input", default="sample36")
    ap.add_argument("--gil-strategy", default="linsi", choices=list(STRATEGY_MODE))
    ap.add_argument("--gil-runs", type=int, default=3)
    ap.add_argument("--skip-sweep", action="store_true", help="only run the GIL demonstration")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    cm.check_binaries()

    host, fp = cm.host(), cm.fp_policy()
    out = args.out or cm.RESULTS / f"pymafft_threads_{dt.date.today().isoformat()}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    bad = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)

        if not args.skip_sweep:
            print("| input | strategy | threads | median wall s | speedup vs threads=1 | all runs = CLI --thread 1 |")
            print("|---|---|---:|---:|---:|---|")
        for inp in ([] if args.skip_sweep else cm.resolve_inputs(args.inputs)):
            records = read_fasta(inp.path)
            for strat in args.strategies:
                flags = cm.MODES[STRATEGY_MODE[strat]] + ["--thread", "1"]
                ref = subprocess.run(cm.command("rust", flags, inp.path), capture_output=True, check=True).stdout
                ref_rows = parse_fasta_bytes(ref)
                base = None
                for t in args.threads:
                    walls, same = [], True
                    for r in range(1, args.runs + 1):
                        t0 = time.perf_counter()
                        res = align(records, strat, t)
                        wall = time.perf_counter() - t0
                        eq = res.to_tuples() == ref_rows
                        same &= eq
                        walls.append(wall)
                        w.writerow([inp.name, inp.n, inp.median_len, strat, t, r, f"{wall:.3f}",
                                    int(eq), f"{cm.loadavg():.2f}", host, fp])
                        fh.flush()
                    med = cm.median(walls)
                    base = med if t == args.threads[0] else base
                    bad += 0 if same else 1
                    print(f"| {inp.name} | {strat} | {t} | {med:.2f} | {base / med:.2f}× | "
                          f"{'yes' if same else 'NO'} |", flush=True)

        # GIL demonstration.
        inp = cm.resolve_inputs([args.gil_input])[0]
        records = read_fasta(inp.path)
        align(records, args.gil_strategy, 1)  # warm-up
        singles, pairs = [], []
        for r in range(1, args.gil_runs + 1):
            t0 = time.perf_counter()
            align(records, args.gil_strategy, 1)
            singles.append(time.perf_counter() - t0)
            w.writerow([inp.name, inp.n, inp.median_len, f"gil:{args.gil_strategy}", 1, r,
                        f"{singles[-1]:.3f}", "", f"{cm.loadavg():.2f}", host, fp])

            results = [None, None]

            def work(i):
                results[i] = align(records, args.gil_strategy, 1)

            ths = [threading.Thread(target=work, args=(i,)) for i in range(2)]
            t0 = time.perf_counter()
            for th in ths:
                th.start()
            for th in ths:
                th.join()
            pairs.append(time.perf_counter() - t0)
            same = results[0].to_tuples() == results[1].to_tuples()
            bad += 0 if same else 1
            w.writerow([inp.name, inp.n, inp.median_len, f"gil:{args.gil_strategy}", 2, r,
                        f"{pairs[-1]:.3f}", int(same), f"{cm.loadavg():.2f}", host, fp])
            fh.flush()
        s, p = cm.median(singles), cm.median(pairs)
        print(f"\nGIL: {inp.name} {args.gil_strategy}, threads=1 per call, {args.gil_runs} runs")
        print(f"  one call                      {s:.2f} s (median)")
        print(f"  two calls in two Python threads {p:.2f} s (median); sequential would be ~{2 * s:.2f} s")
        print(f"  pair / single = {p / s:.2f}  (1.0 = full overlap, 2.0 = serialised)")
        print(f"  load average now {cm.loadavg():.1f}")
    print(f"\nwrote {out}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
