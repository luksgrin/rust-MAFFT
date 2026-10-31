#!/usr/bin/env python3
"""End-to-end high-throughput workload: thousands of per-cluster alignments (issue #23).

Aligns every cluster of the workload corpus (`benchmarks/htworkload_corpus.py`:
2695 S. pneumoniae gene clusters, DNA; 1500 OrthoDB orthogroups, protein)
one alignment per cluster, with the same flags every way:

    dna      --auto --adjustdirection --thread 1 --nuc   (the panaroo command)
    protein  --auto

    pymafft.align(seqs, strategy="auto", adjust_direction=True, threads=1,
                  seq_type="nuc")      # dna
    pymafft.align(seqs, strategy="auto")   # protein

Implementations (`--impl`):

    snake-c     Snakemake (`benchmarks/htworkload/Snakefile`), one job per
                cluster calling the in-tree C wrapper
                (`mafft-upstream/scripts/mafft`, `MAFFT_BINARIES` set);
                `--cores <workers>`
    snake-rust  the same Snakefile with `impl=rust` (`target/release/mafft-rs`)
    loop-c      a plain Python loop (thread pool of <workers>) running the C
                wrapper as a subprocess per cluster: process spawn and FASTA
                round-trip without the Snakemake scheduler
    loop-rust   the same loop with `mafft-rs`
    py-serial   one Python process calling `pymafft.align()` per cluster
                (workers = 1)
    py-process  `pymafft.align()` in a `concurrent.futures.ProcessPoolExecutor`
    py-thread   `pymafft.align()` in a `ThreadPoolExecutor` (the GIL is
                released during alignment, so threads scale like processes
                if nothing else serialises them)

Every path reads the cluster FASTA from disk and writes one aligned FASTA
per cluster to a scratch directory under `target/htworkload/`. CLI paths
discard MAFFT's progress output (stderr to /dev/null); pymafft runs with
its default `quiet=True`.

Each (molecule, impl, workers, rep) cell is one whole-workload run, measured
as a single process tree: `/usr/bin/time -l` (`-v` on Linux) around the
driver (for Snakemake, around `snakemake` itself) gives wall-clock, user and
system CPU summed over every process of the run, and `max_rss_mb`, the peak
RSS of the largest single process. `tree_rss_mb` is the peak of the summed
RSS of the whole process tree, sampled with `ps` every `--sample-interval`
seconds (the memory footprint of N concurrent workers). After the run every
output is hashed and compared with the C outputs of the same clusters
(`identical_c`: byte-identical files for snake-rust / loop-rust; same names
and gapped sequences, case-folded, for pymafft, which does not wrap lines).

USAGE
-----

    make -C mafft-upstream/core
    cargo build --release -p mafft-rs
    python3.12 -m venv .venv && .venv/bin/pip install maturin biopython snakemake
    (cd crates/pymafft && ../../.venv/bin/maturin develop --release)
    .venv/bin/python benchmarks/htworkload_corpus.py build

    # Pipeline check: 50 clusters per molecule, every impl, workers 1 and 4
    .venv/bin/python benchmarks/htworkload.py run --smoke

    # Full workload on a quiet machine (cached per cell under
    # target/htworkload/<host>/; re-running fills in what is missing)
    .venv/bin/python benchmarks/htworkload.py run --reps 3

    # Subsets
    .venv/bin/python benchmarks/htworkload.py run --molecule dna \\
        --impl snake-c py-process --workers 4

    .venv/bin/python benchmarks/htworkload.py summarize \\
        benchmarks/results/htworkload_<host>_<date>.tsv

Default grid per molecule: snake-c, snake-rust, loop-c, loop-rust and
py-process at workers 1, 4, 8 (py-process at 1 is py-serial), py-thread at
4 and 8. `run` writes `benchmarks/results/htworkload_<host>_<date>.tsv`
(`htworkload_smoke_<date>.tsv` with `--smoke`), merging with an existing
file of the same name.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import os
import platform
import random
import re
import shutil
import statistics
import subprocess
import sys
import threading
import time
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "benchmarks"))
import htworkload_corpus as corpus  # noqa: E402

SNAKEFILE = REPO / "benchmarks" / "htworkload" / "Snakefile"
CLI_FLAGS = {"dna": ["--auto", "--adjustdirection", "--thread", "1", "--nuc"],
             "protein": ["--auto"]}
PY_KW = {"dna": {"strategy": "auto", "adjust_direction": True, "threads": 1,
                 "seq_type": "nuc"},
         "protein": {"strategy": "auto"}}
IMPLS = ("snake-c", "snake-rust", "loop-c", "loop-rust", "py-serial", "py-process",
         "py-thread")
C_IMPLS = ("snake-c", "loop-c")
SMOKE_N = 50
COLUMNS = ["host", "uname_m", "cpu", "molecule", "n_clusters", "impl", "workers", "rep",
           "status", "wall_s", "user_s", "sys_s", "max_rss_mb", "tree_rss_mb",
           "n_ok", "identical_c", "load1_start"]


# ---------------------------------------------------------------- host info

def cpu_model() -> str:
    try:
        if sys.platform == "darwin":
            return subprocess.run(["sysctl", "-n", "machdep.cpu.brand_string"],
                                  capture_output=True, text=True).stdout.strip()
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "NA"


# ------------------------------------------------------------ inner drivers

def read_list(path: Path) -> list[tuple[str, str, str]]:
    return [tuple(l.split("\t")) for l in path.read_text().splitlines() if l.strip()]


def cli_argv(impl: str, args, molecule: str, fasta: str) -> tuple[list[str], dict]:
    env = dict(os.environ)
    if impl == "c":
        env["MAFFT_BINARIES"] = str(args.mafft_binaries)
        return [str(args.mafft), *CLI_FLAGS[molecule], fasta], env
    return [str(args.mafft_rs), *CLI_FLAGS[molecule], fasta], env


def loop_child(args) -> int:
    """`htworkload.py _loop`: subprocess per cluster, `workers` at a time."""
    jobs = read_list(Path(args.list))
    out = Path(args.out)

    def one(job):
        cid, mol, fasta = job
        argv, env = cli_argv(args.impl, args, mol, fasta)
        with (out / f"{cid}.fa").open("wb") as fh:
            return cid, subprocess.run(argv, stdout=fh, stderr=subprocess.DEVNULL,
                                       env=env).returncode

    if args.workers == 1:
        results = [one(j) for j in jobs]
    else:
        with cf.ThreadPoolExecutor(args.workers) as ex:
            results = list(ex.map(one, jobs))
    failed = [c for c, rc in results if rc != 0]
    for c in failed:
        print(f"failed: {c}", file=sys.stderr)
    return 1 if failed else 0


def _py_align(job: tuple[str, str, str], out: str) -> tuple[str, str | None]:
    import pymafft
    cid, mol, fasta = job
    try:
        recs = [(h, s) for h, s in corpus.parse_fasta(Path(fasta).read_text())]
        res = pymafft.align(recs, **PY_KW[mol])
        Path(out, f"{cid}.fa").write_text(res.to_fasta())
        return cid, None
    except Exception as e:  # noqa: BLE001 - recorded, the run goes on
        return cid, f"{type(e).__name__}: {e}"


def py_child(args) -> int:
    """`htworkload.py _py`: pymafft.align() per cluster in one Python process."""
    import pymafft  # noqa: F401 - import cost is part of the run
    jobs = read_list(Path(args.list))
    if args.mode == "serial":
        results = [_py_align(j, args.out) for j in jobs]
    else:
        pool = cf.ProcessPoolExecutor if args.mode == "process" else cf.ThreadPoolExecutor
        with pool(args.workers) as ex:
            results = list(ex.map(_py_align, jobs, [args.out] * len(jobs)))
    failed = [(c, e) for c, e in results if e]
    for c, e in failed:
        print(f"failed: {c}: {e}", file=sys.stderr)
    return 1 if failed else 0


# ----------------------------------------------------------------- measuring

class TreeRss:
    """Peak of the summed RSS of a process and its descendants, via `ps`."""

    def __init__(self, pid: int, interval: float):
        self.pid, self.interval, self.peak_kb = pid, interval, 0
        self._stop = threading.Event()
        self._t = threading.Thread(target=self._run, daemon=True)
        if interval > 0:
            self._t.start()

    def _sample(self) -> int:
        out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,rss="], capture_output=True,
                             text=True).stdout
        kids, rss = defaultdict(list), {}
        for line in out.splitlines():
            p = line.split()
            if len(p) == 3:
                kids[int(p[1])].append(int(p[0]))
                rss[int(p[0])] = int(p[2])
        total, stack = 0, [self.pid]
        while stack:
            p = stack.pop()
            total += rss.get(p, 0)
            stack += kids.get(p, [])
        return total

    def _run(self):
        while not self._stop.wait(self.interval):
            try:
                self.peak_kb = max(self.peak_kb, self._sample())
            except (OSError, ValueError):
                pass

    def stop(self) -> float:
        self._stop.set()
        if self._t.is_alive():
            self._t.join()
        return round(self.peak_kb / 1024, 1)


def parse_time(text: str) -> dict:
    """User / sys seconds and max RSS (MB) from BSD `time -l` or GNU `time -v`."""
    out = {}
    m = re.search(r"([\d.]+) real\s+([\d.]+) user\s+([\d.]+) sys", text)
    if m:
        out["user_s"], out["sys_s"] = float(m.group(2)), float(m.group(3))
        r = re.search(r"(\d+)\s+maximum resident set size", text)
        if r:
            out["max_rss_mb"] = round(int(r.group(1)) / (1 << 20), 1)
        return out
    u = re.search(r"User time \(seconds\): ([\d.]+)", text)
    s = re.search(r"System time \(seconds\): ([\d.]+)", text)
    r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if u and s:
        out["user_s"], out["sys_s"] = float(u.group(1)), float(s.group(1))
    if r:
        out["max_rss_mb"] = round(int(r.group(1)) / 1024, 1)
    return out


def measure(argv: list[str], err_path: Path, interval: float, cwd: Path) -> dict:
    flag = "-l" if sys.platform == "darwin" else "-v"
    timer = ["/usr/bin/time", flag] if Path("/usr/bin/time").is_file() else []
    with err_path.open("wb") as err:
        t0 = time.perf_counter()
        p = subprocess.Popen(timer + argv, stdout=err, stderr=subprocess.STDOUT, cwd=cwd)
        tree = TreeRss(p.pid, interval)
        _, status, ru = os.wait4(p.pid, 0)
        wall = time.perf_counter() - t0
        tree_mb = tree.stop()
    rc = os.waitstatus_to_exitcode(status)
    text = err_path.read_text(errors="replace")
    meta = {"wall_s": round(wall, 3), "tree_rss_mb": tree_mb, **parse_time(text)}
    # Fallback (no /usr/bin/time): rusage of the reaped tree.
    meta.setdefault("user_s", round(ru.ru_utime, 3))
    meta.setdefault("sys_s", round(ru.ru_stime, 3))
    meta.setdefault("max_rss_mb", round(ru.ru_maxrss / ((1 << 20) if sys.platform == "darwin"
                                                       else 1024), 1))
    meta["status"] = "ok" if rc == 0 else f"exit_{rc}"
    return meta


def rows_of(text: str) -> list[tuple[str, str]]:
    return [(h.split()[0], s.upper()) for h, s in corpus.parse_fasta(text)]


def output_hashes(out: Path, ids: list[str]) -> dict[str, list[str]]:
    h = {}
    for cid in ids:
        p = out / f"{cid}.fa"
        if p.is_file() and p.stat().st_size > 0:
            b = p.read_bytes()
            r = rows_of(b.decode())
            h[cid] = [hashlib.sha256(b).hexdigest(), hashlib.sha256(
                "\n".join(f"{a}\t{s}" for a, s in r).encode()).hexdigest()]
    return h


# ------------------------------------------------------------------- running

def cells(args) -> list[tuple[str, str, int]]:
    out = []
    for mol in args.molecule:
        for impl in args.impl:
            for w in args.workers:
                if impl == "py-serial" and w != 1:
                    continue
                if impl == "py-process" and w == 1:
                    continue  # that is py-serial
                if impl == "py-thread" and w == 1 and not args.explicit_workers:
                    continue
                out.append((mol, impl, w))
    return out


def select_clusters(args, mol: str) -> list[dict]:
    rows = [r for r in corpus.read_manifest() if r["molecule"] == mol]
    if args.smoke:
        rows = sorted(random.Random(f"htworkload:smoke:{mol}").sample(rows, SMOKE_N),
                      key=lambda r: r["cluster"])
    elif args.limit:
        rows = rows[:args.limit]
    return rows


def run_cell(args, mol: str, impl: str, workers: int, rep: int,
             clusters: list[dict]) -> dict:
    d = args.cache / mol / f"{impl}_w{workers}_rep{rep}"
    meta_p = d / "meta.json"
    if meta_p.is_file():
        meta = json.loads(meta_p.read_text())
        if meta["status"] == "ok" or args.no_retry:
            return meta
    if d.is_dir():
        shutil.rmtree(d)
    out = d / "out"
    out.mkdir(parents=True)
    lst = d / "clusters.tsv"
    lst.write_text("".join(f"{r['cluster']}\t{mol}\t{corpus.fasta_path(r['cluster'], mol)}\n"
                           for r in clusters))
    me = [str(args.python), str(Path(__file__).resolve())]
    tools = ["--mafft", str(args.mafft), "--mafft-binaries", str(args.mafft_binaries),
             "--mafft-rs", str(args.mafft_rs)]
    if impl.startswith("snake-"):
        argv = [str(args.snakemake), "-s", str(SNAKEFILE), "--cores", str(workers),
                "--directory", str(d), "--keep-going",
                "--config", f"impl={impl.split('-')[1]}", f"clusters={lst}", f"outdir={out}",
                f"mafft={args.mafft}", f"mafft_binaries={args.mafft_binaries}",
                f"mafft_rs={args.mafft_rs}"]
    elif impl.startswith("loop-"):
        argv = me + ["_loop", "--impl", impl.split("-")[1], "--workers", str(workers),
                     "--list", str(lst), "--out", str(out), *tools]
    else:
        argv = me + ["_py", "--mode", impl.split("-")[1], "--workers", str(workers),
                     "--list", str(lst), "--out", str(out)]
    load1 = os.getloadavg()[0]
    meta = measure(argv, d / "driver.log", args.sample_interval, d)
    meta.update({"load1_start": round(load1, 2), "n_clusters": len(clusters),
                 "argv": argv})
    ids = [r["cluster"] for r in clusters]
    hashes = output_hashes(out, ids)
    meta["n_ok"] = len(hashes)
    (d / "hashes.json").write_text(json.dumps(hashes))
    if not args.keep:
        shutil.rmtree(out)
        shutil.rmtree(d / ".snakemake", ignore_errors=True)
    meta_p.write_text(json.dumps(meta))
    return meta


def cmd_run(args) -> int:
    args.explicit_workers = args.workers is not None
    if args.smoke:
        args.workers = args.workers or [1, 4]
        args.reps = 1
        args.cache = args.cache / args.host / "smoke"
    else:
        args.workers = args.workers or [1, 4, 8]
        args.cache = args.cache / args.host / ("full" if not args.limit
                                               else f"limit{args.limit}")
    for p, what in ((args.mafft, "C mafft wrapper (make -C mafft-upstream/core)"),
                    (args.mafft_rs, "mafft-rs (cargo build --release -p mafft-rs)"),
                    (args.python, "venv python with pymafft (see USAGE)"),
                    (args.snakemake, "snakemake (pip install snakemake in the venv)")):
        if not Path(p).is_file():
            print(f"error: {what} not found at {p}", file=sys.stderr)
            return 2
    manifest = corpus.read_manifest()
    args.molecule = [m for m in args.molecule if any(r["molecule"] == m for r in manifest)]
    sel = {m: select_clusters(args, m) for m in args.molecule}
    for m, rows in sel.items():
        missing = [r["cluster"] for r in rows
                   if not corpus.fasta_path(r["cluster"], m).is_file()]
        if missing:
            print(f"error: {len(missing)} {m} cluster files missing; run "
                  "benchmarks/htworkload_corpus.py build", file=sys.stderr)
            return 2
    grid = cells(args)
    if args.list:
        for mol, impl, w in grid:
            print(f"{mol}\t{len(sel[mol])}\t{impl}\t{w}")
        return 0
    info = {"host": args.host, "uname_m": platform.machine(), "cpu": cpu_model()}
    print(f"{info}; load {os.getloadavg()}; {len(grid)} cells x {args.reps} reps",
          file=sys.stderr)
    t_start = time.time()
    for rep in range(args.reps):
        # Rotate the order per replicate so drift does not favour one impl.
        order = grid[rep % len(grid):] + grid[:rep % len(grid)] if rep else grid
        # C first in replicate 0: the other impls are compared with it.
        if rep == 0:
            order = sorted(order, key=lambda c: c[1] not in C_IMPLS)
        for k, (mol, impl, w) in enumerate(order, 1):
            m = run_cell(args, mol, impl, w, rep, sel[mol])
            print(f"[rep{rep} {k}/{len(order)} {time.time() - t_start:7.0f}s] {mol} {impl:<10} "
                  f"w{w}: {m['status']} {m['wall_s']:.1f}s wall, "
                  f"{m['user_s'] + m['sys_s']:.1f}s CPU, {m['max_rss_mb']:.0f} MB max, "
                  f"{m['tree_rss_mb']:.0f} MB tree, {m['n_ok']}/{m['n_clusters']} outputs",
                  file=sys.stderr)
    rows = collect(args, grid, info)
    out = args.output or (REPO / "benchmarks" / "results" /
                          (f"htworkload_smoke_{dt.date.today().isoformat()}.tsv" if args.smoke
                           else f"htworkload_{args.host}_{dt.date.today().isoformat()}.tsv"))
    write_tsv(out, rows)
    print(f"wrote {out}", file=sys.stderr)
    print(summarize(rows))
    return 0


def collect(args, grid, info) -> list[dict]:
    rows = []
    for mol in args.molecule:
        cell_meta, cell_hash = {}, {}
        for m2, impl, w in grid:
            if m2 != mol:
                continue
            for rep in range(args.reps):
                d = args.cache / mol / f"{impl}_w{w}_rep{rep}"
                if (d / "meta.json").is_file():
                    cell_meta[(impl, w, rep)] = json.loads((d / "meta.json").read_text())
                    cell_hash[(impl, w, rep)] = json.loads((d / "hashes.json").read_text())
        # Reference: the first C cell, preferring Snakemake at 1 worker.
        refs = sorted((k for k in cell_hash if k[0] in C_IMPLS),
                      key=lambda k: (C_IMPLS.index(k[0]), k[1], k[2]))
        ref = cell_hash[refs[0]] if refs else None
        for (impl, w, rep), m in cell_meta.items():
            ident = "NA"
            if ref is not None:
                h = cell_hash[(impl, w, rep)]
                col = 1 if impl.startswith("py-") else 0
                same = sum(1 for c, v in h.items() if c in ref and ref[c][col] == v[col])
                ident = f"{same}/{m['n_clusters']}"
            rows.append({**info, "molecule": mol, "n_clusters": m["n_clusters"],
                         "impl": impl, "workers": w, "rep": rep, "status": m["status"],
                         "wall_s": m["wall_s"], "user_s": m["user_s"], "sys_s": m["sys_s"],
                         "max_rss_mb": m["max_rss_mb"], "tree_rss_mb": m["tree_rss_mb"],
                         "n_ok": m["n_ok"], "identical_c": ident,
                         "load1_start": m["load1_start"]})
    return rows


def key(r: dict) -> tuple:
    return (r["host"], r["molecule"], str(r["n_clusters"]), r["impl"], str(r["workers"]),
            str(r["rep"]))


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    merged = {key(r): r for r in (read_tsv(path) if path.is_file() else [])}
    merged.update({key(r): {k: str(v) for k, v in r.items()} for r in rows})
    with path.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for k in sorted(merged, key=lambda k: (k[0], corpus.MOLECULES.index(k[1]), int(k[2]),
                                               IMPLS.index(k[3]), int(k[4]), int(k[5]))):
            fh.write("\t".join(merged[k][c] for c in COLUMNS) + "\n")


# ----------------------------------------------------------------- reporting

def summarize(rows: list[dict]) -> str:
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        if r["status"] == "ok":
            groups[(r["molecule"], int(r["n_clusters"]), r["impl"], int(r["workers"]))].append(r)
    med = {k: {f: statistics.median(float(r[f]) for r in v)
               for f in ("wall_s", "user_s", "sys_s", "max_rss_mb", "tree_rss_mb")}
           for k, v in groups.items()}
    lines = ["| Molecule | Clusters | Impl | Workers | Wall s | CPU s | Max RSS MB "
             "| Tree RSS MB | ms / alignment | Speed-up vs snake-c (same workers) "
             "| vs snake-c w1 | Identical to C |",
             "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for k in sorted(med, key=lambda k: (corpus.MOLECULES.index(k[0]), k[1],
                                        IMPLS.index(k[2]), k[3])):
        mol, n, impl, w = k
        m = med[k]
        base = med.get((mol, n, "snake-c", w))
        base1 = med.get((mol, n, "snake-c", 1))
        sp = f"{base['wall_s'] / m['wall_s']:.2f}" if base else "–"
        sp1 = f"{base1['wall_s'] / m['wall_s']:.2f}" if base1 else "–"
        idv = sorted({r["identical_c"] for r in groups[k]})
        lines.append(f"| {mol} | {n} | {impl} | {w} | {m['wall_s']:.1f} | "
                     f"{m['user_s'] + m['sys_s']:.1f} | {m['max_rss_mb']:.0f} | "
                     f"{m['tree_rss_mb']:.0f} | {1000 * m['wall_s'] / n:.1f} | {sp} | {sp1} | "
                     f"{', '.join(idv)} |")
    return "\n".join(lines)


def cmd_summarize(args) -> int:
    print(summarize(read_tsv(Path(args.tsv))))
    return 0


# ------------------------------------------------------------------------ CLI

def add_tools(p) -> None:
    p.add_argument("--mafft", type=Path, default=REPO / "mafft-upstream/scripts/mafft")
    p.add_argument("--mafft-binaries", type=Path, default=REPO / "mafft-upstream/binaries")
    p.add_argument("--mafft-rs", type=Path, default=REPO / "target/release/mafft-rs")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split("Every path")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="run the workload grid (cached, resumable)")
    r.add_argument("--smoke", action="store_true",
                   help=f"{SMOKE_N} clusters per molecule, workers 1 and 4, one replicate")
    r.add_argument("--molecule", nargs="+", default=list(corpus.MOLECULES),
                   choices=corpus.MOLECULES)
    r.add_argument("--impl", nargs="+", default=list(IMPLS), choices=IMPLS)
    r.add_argument("--workers", nargs="+", type=int, default=None,
                   help="default 1 4 8 (1 4 with --smoke)")
    r.add_argument("--reps", type=int, default=1)
    r.add_argument("--limit", type=int, default=0,
                   help="only the first N clusters of each molecule")
    r.add_argument("--host", default=f"{platform.system().lower()}-{platform.machine()}")
    add_tools(r)
    r.add_argument("--python", type=Path, default=REPO / ".venv/bin/python",
                   help="interpreter with pymafft installed")
    r.add_argument("--snakemake", type=Path, default=REPO / ".venv/bin/snakemake")
    r.add_argument("--sample-interval", type=float, default=0.5,
                   help="seconds between process-tree RSS samples (0: off)")
    r.add_argument("--cache", type=Path, default=REPO / "target/htworkload")
    r.add_argument("--output", type=Path, default=None)
    r.add_argument("--keep", action="store_true", help="keep the aligned FASTA files")
    r.add_argument("--no-retry", action="store_true", help="do not retry failed cells")
    r.add_argument("--list", action="store_true", help="print the cells and exit")
    r.set_defaults(func=cmd_run)

    s = sub.add_parser("summarize", help="Markdown table of a results TSV")
    s.add_argument("tsv")
    s.set_defaults(func=cmd_summarize)

    lp = sub.add_parser("_loop", help=argparse.SUPPRESS)
    lp.add_argument("--impl", choices=("c", "rust"), required=True)
    lp.add_argument("--workers", type=int, required=True)
    lp.add_argument("--list", required=True)
    lp.add_argument("--out", required=True)
    add_tools(lp)
    lp.set_defaults(func=loop_child)

    y = sub.add_parser("_py", help=argparse.SUPPRESS)
    y.add_argument("--mode", choices=("serial", "process", "thread"), required=True)
    y.add_argument("--workers", type=int, required=True)
    y.add_argument("--list", required=True)
    y.add_argument("--out", required=True)
    y.set_defaults(func=py_child)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
