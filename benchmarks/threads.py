#!/usr/bin/env python3
"""`--thread N`: determinism and wall-clock speedup, C MAFFT vs mafft-rs.

For every input x mode x thread count the driver runs C and mafft-rs
`--runs` times each, interleaved (C, rs, rs, C, C, rs, ...), and records the
wall-clock time and a hash of stdout. From those runs it checks:

* determinism of mafft-rs: every run at every thread count is
  byte-identical to `mafft-rs --thread 1` (exit status 1 if not);
* C at the same thread count: how many distinct outputs its runs gave, and
  whether they equal C `--thread 1` (C's `--thread N >= 2` is known to be
  nondeterministic; this is recorded, not asserted);
* mafft-rs `--thread 1` vs C `--thread 1` (the byte-identity reference).

`--impls rust` runs a determinism-only sweep (no C), e.g. over a whole
BAliBASE subset: `--inputs ~/data/balibase/bb3_release/RV11/*.tfa`.

Output: benchmarks/results/threads_<date>.tsv with columns
input, n, median_len, mode, flags, threads, impl, run, wall_s, out_sha,
load1, host, fp_policy.
`--summarize FILE` re-prints the summary tables of an existing TSV.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as cm  # noqa: E402

COLUMNS = [
    "input", "n", "median_len", "mode", "flags", "threads", "impl", "run",
    "wall_s", "out_sha", "load1", "host", "fp_policy",
]


def wait_for_load(max_load: float, max_wait: float) -> float:
    """Sleep until the 1-min load average is <= max_load (or max_wait s passed)."""
    waited = 0.0
    while cm.loadavg() > max_load and waited < max_wait:
        time.sleep(15)
        waited += 15
    return cm.loadavg()


def summarize(path: Path) -> int:
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    walls = defaultdict(list)
    shas = defaultdict(set)
    loads = []
    meta = {}
    for r in rows:
        key = (r["input"], r["mode"], int(r["threads"]), r["impl"])
        walls[key].append(float(r["wall_s"]))
        shas[key].add(r["out_sha"])
        loads.append(float(r["load1"]))
        meta[r["input"]] = (r["n"], r["median_len"])

    bad = 0
    groups = sorted({(k[0], k[1]) for k in walls})
    print("## Determinism\n")
    print("| input | mode | mafft-rs: distinct outputs over all thread counts | = rs t1 | "
          "C t1 = rs t1 | C: distinct outputs per thread count |")
    print("|---|---|---|---|---|---|")
    for name, mode in groups:
        ts = sorted({k[2] for k in walls if k[0] == name and k[1] == mode})
        rs_all = set().union(*(shas.get((name, mode, t, "rust"), set()) for t in ts))
        rs1 = shas.get((name, mode, 1, "rust"), set())
        c1 = shas.get((name, mode, 1, "c"), set())
        ok = len(rs_all) == 1
        bad += 0 if ok else 1
        c_cells = ", ".join(
            f"t{t}:{len(shas[(name, mode, t, 'c')])}{'' if shas[(name, mode, t, 'c')] <= c1 else '≠t1'}"
            for t in ts if (name, mode, t, "c") in shas
        ) or "–"
        c_eq = "–" if not c1 else ("yes" if c1 == rs1 else "NO")
        print(f"| {name} | {mode} | {len(rs_all)} | {'yes' if ok else 'NO'} | {c_eq} | {c_cells} |")

    print("\n## Wall clock (median of runs, seconds; speedup vs own --thread 1)\n")
    print("| input | mode | threads | C | C speedup | mafft-rs | rs speedup | C / rs |")
    print("|---|---|---:|---:|---:|---:|---:|---:|")
    for name, mode in groups:
        ts = sorted({k[2] for k in walls if k[0] == name and k[1] == mode})
        c1 = cm.median(walls.get((name, mode, 1, "c"), []))
        r1 = cm.median(walls.get((name, mode, 1, "rust"), []))
        for t in ts:
            c = cm.median(walls.get((name, mode, t, "c"), []))
            r = cm.median(walls.get((name, mode, t, "rust"), []))
            print(f"| {name} | {mode} | {t} | {c:.2f} | {c1 / c:.2f}× | {r:.2f} | {r1 / r:.2f}× | {c / r:.2f} |")
    if loads:
        print(f"\nload average (1 min) during runs: median {cm.median(loads):.1f}, "
              f"min {min(loads):.1f}, max {max(loads):.1f}")
    return 1 if bad else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="*", help="default-set names or paths")
    ap.add_argument("--modes", nargs="*", default=["fftnsi", "linsi", "ginsi"], choices=list(cm.MODES))
    ap.add_argument("--threads", nargs="*", type=int, default=[1, 2, 4, 8, 10])
    ap.add_argument("--impls", nargs="*", default=["c", "rust"], choices=["c", "rust"])
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--max-load", type=float, default=None,
                    help="before each input x mode batch, wait until the 1-min load average is at most this")
    ap.add_argument("--max-wait", type=float, default=1800, help="give up waiting for --max-load after this many s")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--append", action="store_true")
    ap.add_argument("--summarize", type=Path, default=None, help="print the summary of an existing TSV and exit")
    args = ap.parse_args()
    if args.summarize:
        return summarize(args.summarize)
    cm.check_binaries()

    inputs = cm.resolve_inputs(args.inputs)
    today = dt.date.today().isoformat()
    out = args.out or cm.RESULTS / f"threads_{today}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    host, fp = cm.host(), cm.fp_policy()
    mode_flag = "a" if args.append and out.exists() else "w"
    with open(out, mode_flag, newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if mode_flag == "w":
            w.writerow(COLUMNS)
        for inp in inputs:
            for mode in args.modes:
                if args.max_load is not None:
                    load = wait_for_load(args.max_load, args.max_wait)
                    print(f"# {inp.name} {mode}: load {load:.1f}", flush=True)
                for t in args.threads:
                    flags = cm.MODES[mode] + ["--thread", str(t)]
                    for r in range(1, args.runs + 1):
                        for impl in (args.impls if r % 2 else list(reversed(args.impls))):
                            res = cm.run(cm.command(impl, flags, inp.path))
                            if res.returncode != 0:
                                print(f"FAIL {impl} {inp.name} {flags}: {res.stderr[-500:]}", file=sys.stderr)
                                continue
                            load = cm.loadavg()
                            w.writerow([inp.name, inp.n, inp.median_len, mode, " ".join(flags), t, impl, r,
                                        f"{res.wall_s:.3f}", cm.sha(res.stdout), f"{load:.2f}", host, fp])
                            fh.flush()
                            print(f"{inp.name:<22} {mode:<7} t{t:<2} {impl:<4} run{r} {res.wall_s:8.2f} s "
                                  f"{cm.sha(res.stdout)} load {load:.1f}", flush=True)
    print(f"\nwrote {out}\n")
    return summarize(out)


if __name__ == "__main__":
    sys.exit(main())
