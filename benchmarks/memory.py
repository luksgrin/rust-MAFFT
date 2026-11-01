#!/usr/bin/env python3
"""Peak resident set size of C MAFFT (in-tree build) vs mafft-rs.

Each cell (input x mode x threads x impl) is run `--runs` times under
`/usr/bin/time -l` (macOS) or `/usr/bin/time -v` (Linux) and the reported
maximum RSS is recorded. For C the measured process is the `mafft` wrapper
script; the kernel reports the maximum RSS over the process and all of its
waited-for descendants, so the figure is the peak of the largest engine
binary the wrapper spawned (`disttbfast`, `tbfast`, `dvtditr`,
`pairlocalalign`, ...). `--verify-wrapper` checks that claim on one input
by shimming every engine binary with its own `/usr/bin/time`.

Output: benchmarks/results/memory_<date>.tsv with columns
input, n, median_len, mode, flags, threads, impl, run, peak_rss_bytes,
wall_s, host, fp_policy.

Examples:
    python3 benchmarks/memory.py                       # full default grid
    python3 benchmarks/memory.py --inputs sample36 --modes fftns2 --runs 1
    python3 benchmarks/memory.py --verify-wrapper
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import shutil
import stat
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as cm  # noqa: E402

COLUMNS = [
    "input", "n", "median_len", "mode", "flags", "threads", "impl", "run",
    "peak_rss_bytes", "wall_s", "host", "fp_policy",
]


def verify_wrapper(inp: cm.Input, flags: list[str]) -> int:
    """Run C once through the wrapper with every engine binary shimmed by
    `/usr/bin/time`, and compare the max per-binary RSS with the wrapper's."""
    tmp = Path(tempfile.mkdtemp(prefix="mafft-shim-"))
    logs = tmp / "logs"
    logs.mkdir()
    bins = tmp / "bin"
    bins.mkdir()
    flag = "-l" if sys.platform == "darwin" else "-v"
    for b in cm.C_BINARIES.iterdir():
        dst = bins / b.name
        if b.suffix == ".pl" or not os.access(b, os.X_OK):
            shutil.copy2(b, dst)
            continue
        dst.write_text(
            "#!/bin/sh\n"
            f'exec /usr/bin/time {flag} -o "{logs}/{b.name}.$$" "{b}" "$@"\n'
        )
        dst.chmod(dst.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    env = cm.env()
    env["MAFFT_BINARIES"] = str(bins)
    import subprocess
    import time

    cmd = cm.time_prefix() + cm.command("c", flags, inp.path)
    t0 = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, env=env)
    wall = time.perf_counter() - t0
    if p.returncode != 0:
        sys.exit(p.stderr.decode(errors="replace")[-2000:])
    wrapper_rss = cm.parse_rss(p.stderr.decode(errors="replace"))
    per_bin = {}
    for log in logs.iterdir():
        name = log.name.rsplit(".", 1)[0]
        rss = cm.parse_rss(log.read_text())
        per_bin[name] = max(per_bin.get(name, 0), rss)
    shutil.rmtree(tmp)
    print(f"{inp.name} {' '.join(flags)}  ({wall:.1f} s)")
    for name, rss in sorted(per_bin.items(), key=lambda kv: -kv[1]):
        print(f"  {name:<16} {rss / 2**20:9.1f} MiB")
    mx = max(per_bin.values())
    print(f"  max over engines  {mx / 2**20:9.1f} MiB")
    print(f"  wrapper (time)    {wrapper_rss / 2**20:9.1f} MiB")
    ok = wrapper_rss >= mx
    print("  wrapper figure covers the engines:", "yes" if ok else "NO")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="*", help="default-set names (sample36, mtb_cds_120x1400.fa, BB30003.tfa, ...) or paths")
    ap.add_argument("--modes", nargs="*", default=["fftns2", "fftnsi", "linsi", "parttree"], choices=list(cm.MODES))
    ap.add_argument("--threads", nargs="*", type=int, default=[1, 4, 8])
    ap.add_argument("--impls", nargs="*", default=["c", "rust"], choices=["c", "rust"])
    ap.add_argument("--runs", type=int, default=3)
    ap.add_argument("--memsave-inputs", nargs="*", default=None,
                    help="inputs that also get `--memsave` cells (default: the two largest BAliBASE sets)")
    ap.add_argument("--memsave-modes", nargs="*", default=["fftns2", "fftnsi"], choices=list(cm.MODES))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--append", action="store_true", help="append to --out instead of overwriting")
    ap.add_argument("--verify-wrapper", action="store_true",
                    help="check that the wrapper's RSS covers its engine binaries, then exit")
    args = ap.parse_args()
    cm.check_binaries()

    if args.verify_wrapper:
        inp = cm.resolve_inputs(args.inputs or ["BB30003.tfa"])[0]
        rc = 0
        for mode in args.modes:
            for t in args.threads:
                rc |= verify_wrapper(inp, cm.MODES[mode] + ["--thread", str(t)])
        return rc

    inputs = cm.resolve_inputs(args.inputs)
    memsave_names = args.memsave_inputs
    if memsave_names is None:
        memsave_names = [i.name for i in sorted(inputs, key=lambda i: -i.total_len)[:2]]
    memsave = [i for i in inputs if i.name in set(memsave_names)]

    cells = []
    for inp in inputs:
        for mode in args.modes:
            cells.append((inp, mode, cm.MODES[mode]))
    for inp in memsave:
        for mode in args.memsave_modes:
            cells.append((inp, mode + "+memsave", cm.MODES[mode] + ["--memsave"]))

    today = dt.date.today().isoformat()
    out = args.out or cm.RESULTS / f"memory_{today}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    host, fp = cm.host(), cm.fp_policy()
    mode_flag = "a" if args.append and out.exists() else "w"
    summary = defaultdict(list)
    with open(out, mode_flag, newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        if mode_flag == "w":
            w.writerow(COLUMNS)
        for inp, mode, flags in cells:
            for t in args.threads:
                tflags = flags + ["--thread", str(t)]
                for r in range(1, args.runs + 1):
                    # Alternate the order so neither side always runs on a
                    # cold or a warm file cache.
                    order = args.impls if r % 2 else list(reversed(args.impls))
                    for impl in order:
                        res = cm.run(cm.command(impl, tflags, inp.path), measure_rss=True)
                        if res.returncode != 0:
                            print(f"FAIL {impl} {inp.name} {tflags}: {res.stderr[-500:]}", file=sys.stderr)
                            continue
                        w.writerow([inp.name, inp.n, inp.median_len, mode, " ".join(tflags), t, impl, r,
                                    res.peak_rss, f"{res.wall_s:.3f}", host, fp])
                        fh.flush()
                        summary[(inp.name, mode, t, impl)].append(res.peak_rss)
                        print(f"{inp.name:<22} {mode:<16} t{t:<2} {impl:<4} run{r} "
                              f"{res.peak_rss / 2**20:8.1f} MiB {res.wall_s:7.2f} s  load {cm.loadavg():.1f}",
                              flush=True)

    print(f"\nwrote {out}\n")
    print(f"{'input':<22} {'mode':<16} {'thr':>3} {'C med MiB':>10} {'rs med MiB':>10} {'C/rs':>6}")
    for (name, mode, t, impl), vals in sorted(summary.items()):
        if impl != "c":
            continue
        rs = summary.get((name, mode, t, "rust"))
        if not rs:
            continue
        c_med, r_med = cm.median(vals), cm.median(rs)
        print(f"{name:<22} {mode:<16} {t:>3} {c_med / 2**20:10.1f} {r_med / 2**20:10.1f} {c_med / r_med:6.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
