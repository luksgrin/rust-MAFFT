"""pymafft: Python bindings for MAFFT-rs multiple sequence alignment.

Every function here builds a ``mafft-rs`` command line from its keyword
arguments and runs it in-process through the same flag layer as the CLI
(``mafft_rs::run_from_seqs``). Nothing about a flag's meaning is decided in
Python: ``strategy="auto"`` is ``--auto``, ``adjust_direction=True`` is
``--adjustdirection``, ``seq_type="nuc"`` is ``--nuc``, and so on, so the
alignment you get back is byte-identical to what ``mafft-rs`` prints for
the same input and flags. Flags without a keyword are reachable through
:func:`run`.
"""

from __future__ import annotations

import re
from typing import Any, Callable, Iterable, Optional, Sequence, Union

from .pymafft import (
    AlignedSequence,
    AlignmentResult,
    MafftError,
    __version__,
    _run,
    _run_fasta_string,
    _run_file,
)

__all__ = [
    "AlignedSequence",
    "AlignmentResult",
    "MafftError",
    "align",
    "align_file",
    "align_fasta_string",
    "run",
    "__version__",
]

ProgressCallback = Callable[[str], None]

# strategy -> (mode flag, refinement cycles used when maxiterate == 0).
# `None` cycles means "do not pass --maxiterate": FFT-NS-2 is the CLI's
# default and ignores `maxiterate` (as pymafft always has); `auto` lets
# `--auto` pick the strategy from the input size.
_STRATEGIES = {
    "fftns2": (None, None),
    "default": (None, None),
    "": (None, None),
    "fftnsi": (None, 100),
    "ginsi": ("--globalpair", 1000),
    "globalpair": ("--globalpair", 1000),
    "linsi": ("--localpair", 1000),
    "localpair": ("--localpair", 1000),
    "einsi": ("--genafpair", 1000),
    "genafpair": ("--genafpair", 1000),
    "auto": ("--auto", None),
}

_VALID_STRATEGIES = "fftns2, fftnsi, ginsi, linsi, einsi, auto"

# scoring="bl62" -> --bl 62, "jtt200" -> --jtt 200, "tm100" -> --tm 100.
_SCORING_RE = re.compile(r"(blosum|bl|jtt|tm)\s*(\d+)")
_SCORING_FLAG = {"blosum": "--bl", "bl": "--bl", "jtt": "--jtt", "tm": "--tm"}


def _flags(
    strategy: str,
    maxiterate: int,
    *,
    seq_type: Optional[str],
    adjust_direction: Union[bool, str],
    threads: int,
    reorder: bool,
    retree: Optional[int],
    scoring: Optional[str],
    gap_open: Optional[float],
    gap_extend: Optional[float],
    quiet: bool,
    progress: Optional[ProgressCallback],
) -> list[str]:
    """Turn the keyword options into the equivalent ``mafft-rs`` flags."""
    if not isinstance(strategy, str):
        raise TypeError(f"strategy must be a str, got {type(strategy).__name__}")
    key = strategy.lower()
    if key not in _STRATEGIES:
        raise ValueError(f"unknown strategy '{strategy}'. Valid: {_VALID_STRATEGIES}")
    if isinstance(maxiterate, bool) or not isinstance(maxiterate, int):
        raise TypeError(f"maxiterate must be an int, got {type(maxiterate).__name__}")
    if maxiterate < 0:
        raise ValueError("maxiterate must be >= 0")

    mode_flag, default_cycles = _STRATEGIES[key]
    args: list[str] = []
    if mode_flag is not None:
        args.append(mode_flag)
    if key == "auto":
        cycles = maxiterate if maxiterate > 0 else None
    elif default_cycles is None:
        cycles = None  # FFT-NS-2: `maxiterate` is ignored, use "fftnsi" to refine.
    else:
        cycles = maxiterate if maxiterate > 0 else default_cycles
    if cycles is not None:
        args += ["--maxiterate", str(cycles)]

    if adjust_direction is True:
        args.append("--adjustdirection")
    elif adjust_direction == "accurately":
        args.append("--adjustdirectionaccurately")
    elif adjust_direction is not False:
        raise ValueError(
            f"adjust_direction must be True, False or 'accurately', got {adjust_direction!r}"
        )

    if isinstance(threads, bool) or not isinstance(threads, int):
        raise TypeError(f"threads must be an int, got {type(threads).__name__}")
    if threads < 0:
        raise ValueError("threads must be >= 0 (0 = all cores)")
    if threads > 0:
        args += ["--thread", str(threads)]

    if seq_type is not None:
        if seq_type not in ("nuc", "amino"):
            raise ValueError(f"seq_type must be None, 'nuc' or 'amino', got {seq_type!r}")
        args.append(f"--{seq_type}")

    if reorder:
        args.append("--reorder")

    if retree is not None:
        if isinstance(retree, bool) or not isinstance(retree, int) or retree < 0:
            raise ValueError(f"retree must be a non-negative int or None, got {retree!r}")
        args += ["--retree", str(retree)]

    if scoring is not None:
        m = _SCORING_RE.fullmatch(scoring.strip().lower())
        if m is None:
            raise ValueError(
                f"scoring must look like 'bl62', 'jtt200' or 'tm100', got {scoring!r}"
            )
        args += [_SCORING_FLAG[m.group(1)], m.group(2)]

    if gap_open is not None:
        args += ["--op", repr(float(gap_open))]
    if gap_extend is not None:
        args += ["--ep", repr(float(gap_extend))]

    # A progress callable wants the messages, so `--quiet` is only passed
    # when nobody is listening.
    if quiet and progress is None:
        args.append("--quiet")

    return args


def align(
    sequences: Union[Iterable[str], Iterable[tuple[str, str]], Iterable[Any]],
    strategy: str = "fftns2",
    maxiterate: int = 0,
    *,
    seq_type: Optional[str] = None,
    adjust_direction: Union[bool, str] = False,
    threads: int = 0,
    reorder: bool = False,
    retree: Optional[int] = None,
    scoring: Optional[str] = None,
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    quiet: bool = True,
    progress: Optional[ProgressCallback] = None,
) -> AlignmentResult:
    """Align sequences.

    Every option maps to one ``mafft-rs`` command-line flag and is applied
    by the CLI's own flag layer, so the result is byte-identical to running
    ``mafft-rs`` with those flags on a FASTA file of the same sequences.

    Args:
        sequences: The sequences to align. Accepted shapes:
            a list of strings (auto-named ``seq_1``, ``seq_2``, ...);
            a list of ``(name, sequence)`` tuples;
            an iterable of objects with ``.id`` and ``.seq`` attributes
            (e.g. Biopython ``SeqRecord``).
        strategy: ``"fftns2"`` (default; FFT-NS-2), ``"fftnsi"``
            (``--maxiterate N``), ``"ginsi"`` (``--globalpair``),
            ``"linsi"`` (``--localpair``), ``"einsi"`` (``--genafpair``) or
            ``"auto"`` (``--auto``: strategy chosen from the input size).
        maxiterate: Refinement cycles (``--maxiterate``). ``0`` means the
            strategy's default: 100 for ``fftnsi``, 1000 for ``ginsi`` /
            ``linsi`` / ``einsi``; ``fftns2`` never refines and ``auto``
            decides for itself.
        seq_type: ``None`` detects the type from the residues as the CLI
            does; ``"nuc"`` / ``"amino"`` force it (``--nuc`` / ``--amino``).
        adjust_direction: ``True`` for ``--adjustdirection`` (k-mer strand
            detection, DNA only), ``"accurately"`` for
            ``--adjustdirectionaccurately`` (DP-based).
        threads: ``--thread N``; ``0`` (default) uses all cores.
        reorder: ``--reorder``: output rows in guide-tree order.
        retree: ``--retree N`` guide-tree rebuilds (CLI default 2).
        scoring: Substitution matrix: ``"bl30"`` ... ``"bl80"`` (``--bl``),
            ``"jtt200"`` (``--jtt``), ``"tm100"`` (``--tm``).
        gap_open: ``--op`` gap opening penalty (CLI default 1.53).
        gap_extend: ``--ep`` offset / gap extension penalty (CLI default 0.123).
        quiet: Pass ``--quiet``, suppressing progress messages. Ignored
            when ``progress`` is given.
        progress: Callable receiving each progress line (the text the CLI
            prints on stderr without ``--quiet``, one call per line, no
            trailing newline). With ``progress=None`` and ``quiet=False``
            the lines go to the process's stderr, as on the command line.

    Returns:
        AlignmentResult with the aligned sequences. Nucleotide output is
        lowercase and protein output uppercase, as with C MAFFT.

    Raises:
        MafftError: When the alignment fails, e.g. on an illegal residue;
            ``.code`` and ``.message`` are the CLI's exit code and stderr
            text. A ``ValueError`` subclass.
        ValueError: If ``sequences`` is empty or has an unsupported shape,
            or an option value is invalid.
    """
    args = _flags(
        strategy, maxiterate, seq_type=seq_type, adjust_direction=adjust_direction,
        threads=threads, reorder=reorder, retree=retree, scoring=scoring,
        gap_open=gap_open, gap_extend=gap_extend, quiet=quiet, progress=progress,
    )
    return _run(args, sequences, progress)


def align_file(
    path: str,
    strategy: str = "fftns2",
    maxiterate: int = 0,
    *,
    seq_type: Optional[str] = None,
    adjust_direction: Union[bool, str] = False,
    threads: int = 0,
    reorder: bool = False,
    retree: Optional[int] = None,
    scoring: Optional[str] = None,
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    quiet: bool = True,
    progress: Optional[ProgressCallback] = None,
) -> AlignmentResult:
    """Align sequences from a FASTA file.

    Same options as :func:`align`; the file is read with the CLI's FASTA
    reader.

    Args:
        path: Path to the FASTA file.

    Raises:
        ValueError: If the file cannot be read or parsed.
        MafftError: If the alignment fails (see :func:`align`).
    """
    args = _flags(
        strategy, maxiterate, seq_type=seq_type, adjust_direction=adjust_direction,
        threads=threads, reorder=reorder, retree=retree, scoring=scoring,
        gap_open=gap_open, gap_extend=gap_extend, quiet=quiet, progress=progress,
    )
    return _run_file(args, str(path), progress)


def align_fasta_string(
    fasta_string: str,
    strategy: str = "fftns2",
    maxiterate: int = 0,
    *,
    seq_type: Optional[str] = None,
    adjust_direction: Union[bool, str] = False,
    threads: int = 0,
    reorder: bool = False,
    retree: Optional[int] = None,
    scoring: Optional[str] = None,
    gap_open: Optional[float] = None,
    gap_extend: Optional[float] = None,
    quiet: bool = True,
    progress: Optional[ProgressCallback] = None,
) -> AlignmentResult:
    """Align sequences from a FASTA-formatted string.

    Same options as :func:`align`; the text is parsed with the CLI's FASTA
    reader.

    Args:
        fasta_string: FASTA-formatted text.

    Raises:
        ValueError: If the string cannot be parsed as FASTA.
        MafftError: If the alignment fails (see :func:`align`).
    """
    args = _flags(
        strategy, maxiterate, seq_type=seq_type, adjust_direction=adjust_direction,
        threads=threads, reorder=reorder, retree=retree, scoring=scoring,
        gap_open=gap_open, gap_extend=gap_extend, quiet=quiet, progress=progress,
    )
    return _run_fasta_string(args, fasta_string, progress)


def run(
    args: Sequence[str],
    sequences: Union[Iterable[str], Iterable[tuple[str, str]], Iterable[Any]],
    *,
    progress: Optional[ProgressCallback] = None,
) -> AlignmentResult:
    """Align with an explicit ``mafft-rs`` command line.

    The escape hatch for flags that have no keyword in :func:`align`:
    ``args`` is passed verbatim (without the program name and without an
    input path; the sequences take its place) to the same entry point
    :func:`align` uses, so ``run(["--auto", "--nuc", "--quiet"], seqs)``
    equals ``align(seqs, strategy="auto", seq_type="nuc")``.

    Args:
        args: The flags, e.g. ``["--localpair", "--maxiterate", "1000",
            "--quiet"]``. Progress goes to stderr unless ``--quiet`` is
            present or ``progress`` is given.
        sequences: As for :func:`align`.
        progress: As for :func:`align`. The flags are not modified, so add
            ``--quiet`` yourself to silence stderr when ``progress`` is
            ``None``.

    Raises:
        MafftError: For any failure the CLI would exit on, including
            unknown flags (``.message`` is clap's usage text) and an INPUT
            path in ``args``.
    """
    if isinstance(args, (str, bytes)):
        raise TypeError("args must be a sequence of flag strings, not a single string")
    return _run([str(a) for a in args], sequences, progress)
