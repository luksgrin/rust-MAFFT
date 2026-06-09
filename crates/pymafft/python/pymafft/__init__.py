"""pymafft: Python bindings for MAFFT-rs multiple sequence alignment."""

from .pymafft import (
    AlignedSequence,
    AlignmentResult,
    align,
    align_file,
    align_fasta_string,
    __version__,
)

__all__ = [
    "AlignedSequence",
    "AlignmentResult",
    "align",
    "align_file",
    "align_fasta_string",
    "__version__",
]
