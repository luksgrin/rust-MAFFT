"""Type stubs for pymafft.

Every alignment function builds a ``mafft-rs`` command line from its keyword
arguments and runs it in-process through the CLI's flag layer
(``mafft_rs::run_from_seqs``), so the result is byte-identical to the
command line's for the same input and flags.
"""

from typing import Any, Callable, Iterator, Optional, Sequence, Union

__version__: str

ProgressCallback = Callable[[str], None]
SequenceInput = Union[Sequence[str], Sequence[tuple[str, str]], Sequence[Any]]

class MafftError(ValueError):
    """Raised when the alignment fails.

    Carries the exact stderr text and exit code the ``mafft-rs`` command line
    would have produced for the same input and flags. A ``ValueError``
    subclass, so ``except ValueError`` keeps working.
    """

    code: int
    """Exit code the CLI would terminate with (1 for an illegal residue,
    2 for a clap usage error, 0 for an informational early exit)."""
    message: str
    """The CLI's stderr text, without the trailing newline; also ``str(err)``."""

class AlignedSequence:
    """A single aligned sequence with name and gapped data."""

    @property
    def name(self) -> str: ...
    @property
    def sequence(self) -> str: ...
    def ungapped(self) -> str:
        """Return the sequence without gap characters."""
        ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...

class AlignmentResult:
    """Result of a multiple sequence alignment."""

    @property
    def sequences(self) -> list[AlignedSequence]: ...
    @property
    def width(self) -> int: ...
    @property
    def score(self) -> float: ...
    @property
    def nseq(self) -> int: ...
    def to_fasta(self) -> str:
        """Return the alignment as a FASTA-formatted string."""
        ...
    def to_tuples(self) -> list[tuple[str, str]]:
        """Return the alignment as a list of (name, sequence) tuples."""
        ...
    def to_biopython(self) -> Any:
        """Return the alignment as a ``Bio.Align.MultipleSeqAlignment``.

        Each row becomes a ``Bio.SeqRecord.SeqRecord`` with ``id`` set to the
        sequence name and an empty ``description``. Biopython is imported
        lazily here and only here; it is not a dependency of pymafft.

        Raises:
            ImportError: If Biopython is not installed.
        """
        ...
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __getitem__(self, idx: int) -> AlignedSequence: ...
    def __iter__(self) -> Iterator[AlignedSequence]: ...

def align(
    sequences: SequenceInput,
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

    Every option maps to one ``mafft-rs`` command-line flag and is applied by
    the CLI's own flag layer, so the result is byte-identical to running
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
    ...

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
    ...

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
    ...

def run(
    args: Sequence[str],
    sequences: SequenceInput,
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
    ...
