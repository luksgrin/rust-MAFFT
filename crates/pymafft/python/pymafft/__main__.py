"""Console-script entry point for the bundled `mafft-rs` binary.

The wheel ships the compiled `mafft-rs` binary inside the package at
`pymafft/_bin/mafft-rs` (or `mafft-rs.exe` on Windows). `pip install
pymafft` exposes a `mafft-rs` command that runs this module, which then
hands off (via `os.execvp` on POSIX, `subprocess.run` on Windows) to the
native binary. Startup overhead is one Python interpreter launch (~30 ms);
the alignment itself runs at native speed.

Equivalent invocations after install:

    mafft-rs --help
    python -m pymafft --help
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def _bundled_binary() -> Path:
    """Locate the `mafft-rs` binary bundled inside the wheel.

    Raises:
        RuntimeError: if the binary is missing — typically because the
            wheel was built without the pre-build step that copies
            `target/release/mafft-rs` into `python/pymafft/_bin/`. Source
            installs (`pip install --no-binary :all: pymafft`) currently
            do not produce the CLI launcher; use the published wheels or
            `cargo install mafft-rs` instead.
    """
    bin_dir = Path(__file__).parent / "_bin"
    for name in ("mafft-rs", "mafft-rs.exe"):
        candidate = bin_dir / name
        if candidate.exists():
            return candidate
    raise RuntimeError(
        f"mafft-rs binary not found under {bin_dir}. The CLI is only "
        "available from pre-built wheels (or `cargo install mafft-rs`); "
        "source installs do not bundle the binary."
    )


def main() -> int:
    binary = _bundled_binary()
    # Ensure +x. pip-installed files on POSIX should already be executable
    # (maturin marks the wheel data as such), but be defensive.
    if os.name != "nt":
        mode = binary.stat().st_mode
        if not (mode & 0o111):
            binary.chmod(mode | 0o755)
        # execvp replaces this process; signals/stdio pass through cleanly.
        os.execv(str(binary), [str(binary), *sys.argv[1:]])
        # Unreachable.
        return 0
    # Windows has no execv-by-name with stdio inheritance that behaves the
    # same as subprocess; spawn a child and propagate its exit code.
    completed = subprocess.run([str(binary), *sys.argv[1:]])
    return completed.returncode


if __name__ == "__main__":
    sys.exit(main())
