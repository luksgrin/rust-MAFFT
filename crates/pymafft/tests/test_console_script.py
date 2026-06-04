"""Smoke tests for the bundled `mafft-rs` CLI.

`pip install pymafft` exposes a `mafft-rs` command on `$PATH` via the
`[project.scripts]` entry in pyproject.toml; that launcher (defined in
`pymafft.__main__`) hands off to the native binary bundled inside the
wheel at `pymafft/_bin/mafft-rs`.

These tests verify:

  * The launcher is importable and resolves the bundled binary.
  * `python -m pymafft --help` works (Python entry-point path).
  * `mafft-rs --help` works when the wheel is installed (PATH path).

If the binary is not bundled (e.g. running directly against a `maturin
develop` build that skipped the pre-build step), the tests skip rather
than fail — the launcher is correct, the binary is just absent.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import pymafft


def _bin_dir() -> Path:
    return Path(pymafft.__file__).parent / "_bin"


def _bundled_present() -> bool:
    d = _bin_dir()
    return (d / "mafft-rs").exists() or (d / "mafft-rs.exe").exists()


pytestmark = pytest.mark.skipif(
    not _bundled_present(),
    reason=(
        "bundled mafft-rs binary not present under pymafft/_bin/ — "
        "this only ships in wheels built with the pre-build step. "
        "Run `cargo build --release -p mafft-rs && cp target/release/"
        "mafft-rs crates/pymafft/python/pymafft/_bin/` to enable."
    ),
)


def test_launcher_locates_binary():
    """`pymafft.__main__._bundled_binary()` returns an existing path."""
    from pymafft.__main__ import _bundled_binary

    binary = _bundled_binary()
    assert binary.exists()
    assert binary.is_file()


def test_python_m_pymafft_help():
    """`python -m pymafft --help` exits 0 and prints something."""
    result = subprocess.run(
        [sys.executable, "-m", "pymafft", "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, (
        f"returncode={result.returncode}\nstdout={result.stdout}\n"
        f"stderr={result.stderr}"
    )
    # Clap's help banner starts with "Usage:" — that's the cheapest stable
    # marker. We don't pin the rest because the help text is rich and may
    # legitimately change between releases.
    out = result.stdout + result.stderr
    assert "Usage:" in out


def test_python_m_pymafft_version():
    """`python -m pymafft --version` prints a version string."""
    result = subprocess.run(
        [sys.executable, "-m", "pymafft", "--version"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0
    out = (result.stdout + result.stderr).strip()
    # Either "mafft-rs 0.1.0" or just "0.1.0" — both are fine.
    assert "0.1" in out or "0.2" in out or "1." in out


def test_console_script_on_path():
    """If installed via pip, `mafft-rs` should be on $PATH (`[project.scripts]`)."""
    if not shutil.which("mafft-rs"):
        pytest.skip(
            "`mafft-rs` not on $PATH — this only happens after `pip install "
            "pymafft` registers the console script. Local `maturin develop` "
            "doesn't create $PATH entries."
        )
    result = subprocess.run(
        ["mafft-rs", "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0
    assert "Usage:" in result.stdout + result.stderr


def test_console_script_aligns_real_input(tmp_path):
    """End-to-end: the bundled `mafft-rs` actually aligns a FASTA file."""
    fasta = tmp_path / "in.fa"
    fasta.write_text(">a\nACDEFGHIK\n>b\nACDEFHIK\n>c\nACDHIK\n")

    # Use `python -m pymafft` so the test works pre-`pip install` (when
    # the console script isn't on $PATH yet). Add --quiet to suppress
    # progress output.
    result = subprocess.run(
        [sys.executable, "-m", "pymafft", "--quiet", str(fasta)],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, (
        f"returncode={result.returncode}\nstderr={result.stderr}"
    )
    # Output is FASTA: 3 headers, 3 same-length sequence lines.
    lines = [l for l in result.stdout.splitlines() if l]
    headers = [l for l in lines if l.startswith(">")]
    seqs = [l for l in lines if not l.startswith(">")]
    assert len(headers) == 3
    # Aligned rows must all share the same width.
    widths = {len(s) for s in seqs}
    assert len(widths) == 1, f"ragged alignment widths: {widths}"
