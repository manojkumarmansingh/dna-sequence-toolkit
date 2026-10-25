"""DNA Sequence Analysis Toolkit.

A small, dependency-free toolkit for common DNA sequence analyses:
composition and GC content, reverse complement, transcription, translation,
open-reading-frame finding, restriction-site mapping and FASTA parsing.

Example:
    >>> from dna_toolkit import gc_content, reverse_complement
    >>> gc_content("ATGCGC")
    66.66666666666666
    >>> reverse_complement("ATGC")
    'GCAT'
"""

from .core import (
    DNA_ALPHABET,
    RESTRICTION_ENZYMES,
    START_CODON,
    STOP_CODONS,
    base_composition,
    clean_sequence,
    find_orfs,
    find_restriction_sites,
    gc_content,
    kmer_counts,
    molecular_weight,
    read_fasta,
    read_fasta_string,
    reverse_complement,
    summarize,
    transcribe,
    translate,
    validate_sequence,
)

__version__ = "0.1.0"

__all__ = [
    "DNA_ALPHABET",
    "RESTRICTION_ENZYMES",
    "START_CODON",
    "STOP_CODONS",
    "base_composition",
    "clean_sequence",
    "find_orfs",
    "find_restriction_sites",
    "gc_content",
    "kmer_counts",
    "molecular_weight",
    "read_fasta",
    "read_fasta_string",
    "reverse_complement",
    "summarize",
    "transcribe",
    "translate",
    "validate_sequence",
    "__version__",
]
