"""Shared helpers for the benchmark drivers in this directory.

Paths are resolved relative to the repository root, so the scripts can be
run from anywhere. The C reference is the in-tree build
(`make -C mafft-upstream/core`), driven through the upstream wrapper script
`mafft-upstream/scripts/mafft` with `MAFFT_BINARIES` pointing at
`mafft-upstream/binaries`, exactly as CI does
(`.github/workflows/build-test.yml`).
"""

from __future__ import annotations

import hashlib
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
C_WRAPPER = ROOT / "mafft-upstream" / "scripts" / "mafft"
C_BINARIES = ROOT / "mafft-upstream" / "binaries"
RUST_BIN = ROOT / "target" / "release" / "mafft-rs"
FIXTURES = ROOT / "crates" / "mafft-bin" / "tests" / "fixtures"
BALIBASE = Path(os.environ.get("BALIBASE_DIR", Path.home() / "data" / "balibase" / "bb3_release"))
RESULTS = Path(__file__).resolve().parent / "results"

# Mode name -> mafft flags (identical on both sides).
MODES = {
    "fftns2": [],
    "fftnsi": ["--maxiterate", "1000"],
    "linsi": ["--localpair", "--maxiterate", "1000"],
    "ginsi": ["--globalpair", "--maxiterate", "1000"],
    "parttree": ["--parttree"],
}

# The five largest BAliBASE 3 sets by total residues, plus the default
# small/DNA inputs. Missing files are skipped with a warning.
BALIBASE_LARGEST = [
    "RV30/BB30003.tfa",
    "RV30/BB30013.tfa",
    "RV30/BB30001.tfa",
    "RV30/BB30005.tfa",
    "RV40/BB40049.tfa",
]


@dataclass(frozen=True)
class Input:
    name: str
    path: Path
    n: int
    total_len: int
    median_len: float


def read_lengths(path: Path) -> list[int]:
    lens: list[int] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                lens.append(0)
            elif lens:
                lens[-1] += len(line.replace("-", ""))
    return lens


def make_input(path: Path, name: str | None = None) -> Input:
    lens = read_lengths(path)
    return Input(
        name=name or path.name,
        path=path,
        n=len(lens),
        total_len=sum(lens),
        median_len=statistics.median(lens),
    )


def default_inputs() -> list[Input]:
    out = [
        make_input(ROOT / "mafft-upstream" / "test" / "sample", "sample36"),
        make_input(FIXTURES / "mtb_cds_120x1400.fa"),
    ]
    for rel in BALIBASE_LARGEST:
        p = BALIBASE / rel
        if p.exists():
            out.append(make_input(p))
        else:
            print(f"warning: {p} not found, skipped", file=sys.stderr)
    return out


def resolve_inputs(names: list[str] | None) -> list[Input]:
    """`None` -> the default set; otherwise paths or names from the default set."""
    if not names:
        return default_inputs()
    by_name = {i.name: i for i in default_inputs()}
    out = []
    for n in names:
        if n in by_name:
            out.append(by_name[n])
        else:
            out.append(make_input(Path(n)))
    return out


def command(impl: str, flags: list[str], path: Path) -> list[str]:
    exe = str(C_WRAPPER) if impl == "c" else str(RUST_BIN)
    return [exe, "--quiet", *flags, str(path)]


def env() -> dict[str, str]:
    e = dict(os.environ)
    e["MAFFT_BINARIES"] = str(C_BINARIES)
    e["LC_ALL"] = "C"
    return e


def check_binaries() -> None:
    missing = [p for p in (C_WRAPPER, C_BINARIES / "disttbfast", RUST_BIN) if not p.exists()]
    if missing:
        sys.exit(
            "missing: " + ", ".join(map(str, missing)) + "\n"
            "build with `make -C mafft-upstream/core` and `cargo build --release -p mafft-rs`"
        )


def host() -> str:
    cpu = platform.processor() or platform.machine()
    if sys.platform == "darwin":
        try:
            cpu = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
            ).stdout.strip() or cpu
        except OSError:
            pass
    else:
        try:
            for line in open("/proc/cpuinfo"):
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
        except OSError:
            pass
    return f"{cpu}/{platform.system()}-{platform.machine()}/{os.cpu_count()}cpu".replace(" ", "_")


def fp_policy() -> str:
    """Floating-point contraction policy of the default builds on this host.

    The in-tree C build contracts a*b+c into FMA on arm64/clang and not on
    x86-64/gcc; `mafft-rs` mirrors the host default (see
    docs/architecture/byte-identity.md).
    """
    return "fma" if platform.machine() in ("arm64", "aarch64") else "none"


def loadavg() -> float:
    try:
        return os.getloadavg()[0]
    except OSError:
        return float("nan")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()[:16]


_RSS_RE_MAC = re.compile(r"^\s*(\d+)\s+maximum resident set size", re.M)
_RSS_RE_LINUX = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")


def time_prefix() -> list[str]:
    return ["/usr/bin/time", "-l"] if sys.platform == "darwin" else ["/usr/bin/time", "-v"]


def parse_rss(stderr: str) -> int:
    """Peak RSS in bytes from `/usr/bin/time -l` (macOS, bytes) or `-v` (Linux, KiB)."""
    m = _RSS_RE_MAC.search(stderr)
    if m:
        return int(m.group(1))
    m = _RSS_RE_LINUX.search(stderr)
    if m:
        return int(m.group(1)) * 1024
    raise RuntimeError("no maximum resident set size in time output:\n" + stderr[-2000:])


@dataclass
class RunResult:
    returncode: int
    stdout: bytes
    stderr: str
    wall_s: float
    peak_rss: int | None


def run(cmd: list[str], measure_rss: bool = False) -> RunResult:
    full = (time_prefix() + cmd) if measure_rss else cmd
    t0 = time.perf_counter()
    p = subprocess.run(full, capture_output=True, env=env())
    wall = time.perf_counter() - t0
    err = p.stderr.decode(errors="replace")
    rss = parse_rss(err) if measure_rss and p.returncode == 0 else None
    return RunResult(p.returncode, p.stdout, err, wall, rss)


def median(xs):
    return statistics.median(xs) if xs else float("nan")
