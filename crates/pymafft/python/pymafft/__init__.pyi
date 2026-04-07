"""Type stubs for pymafft."""

from typing import Sequence, Union, Iterator

__version__: str

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
    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...
    def __getitem__(self, idx: int) -> AlignedSequence: ...
    def __iter__(self) -> Iterator[AlignedSequence]: ...

def align(
    sequences: Union[Sequence[str], Sequence[tuple[str, str]]],
    strategy: str = "fftns2",
    maxiterate: int = 0,
) -> AlignmentResult:
    """Align sequences.

    Args:
        sequences: List of sequences. Can be:
            - List of strings (auto-named "seq_1", "seq_2", ...)
            - List of (name, sequence) tuples
        strategy: Alignment strategy. One of:
            "fftns2" (default, fast), "fftnsi", "ginsi", "linsi", "einsi"
        maxiterate: Maximum refinement iterations (0 = default for strategy).

    Returns:
        AlignmentResult with aligned sequences.

    Raises:
        ValueError: If sequences is empty or strategy is invalid.
    """
    ...

def align_file(
    path: str,
    strategy: str = "fftns2",
    maxiterate: int = 0,
) -> AlignmentResult:
    """Align sequences from a FASTA file.

    Args:
        path: Path to FASTA file.
        strategy: Alignment strategy (default: "fftns2").
        maxiterate: Maximum refinement iterations.

    Returns:
        AlignmentResult with aligned sequences.

    Raises:
        ValueError: If the file cannot be read or parsed.
    """
    ...

def align_fasta_string(
    fasta_string: str,
    strategy: str = "fftns2",
    maxiterate: int = 0,
) -> AlignmentResult:
    """Align sequences from a FASTA-formatted string.

    Args:
        fasta_string: FASTA-formatted string.
        strategy: Alignment strategy (default: "fftns2").
        maxiterate: Maximum refinement iterations.

    Returns:
        AlignmentResult with aligned sequences.

    Raises:
        ValueError: If the string cannot be parsed as FASTA.
    """
    ...
