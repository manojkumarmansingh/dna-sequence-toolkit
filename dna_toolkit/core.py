"""Core sequence-analysis functions for the DNA Sequence Analysis Toolkit.

Implemented with the Python standard library only, so the package runs
anywhere without third-party dependencies.

All coordinate values returned by this module are 1-based and inclusive,
following the convention used in biology.

Public API
----------
This module is intentionally stable: the functions listed in ``__all__`` keep
their existing signatures. New behaviour is added through optional keyword
arguments with backwards-compatible defaults.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Set

__all__ = [
    "DNA_ALPHABET",
    "STOP_CODONS",
    "START_CODON",
    "START_CODONS",
    "RESTRICTION_ENZYMES",
    "clean_sequence",
    "validate_sequence",
    "gc_content",
    "base_composition",
    "reverse_complement",
    "transcribe",
    "translate",
    "find_orfs",
    "find_restriction_sites",
    "molecular_weight",
    "kmer_counts",
    "read_fasta",
    "read_fasta_string",
    "summarize",
]

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

DNA_ALPHABET = set("ACGTN")
RNA_ALPHABET = set("ACGUN")

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_TRANSCRIBE = str.maketrans("Tt", "Uu")

# Standard genetic code (translation table 1).
CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODON = "ATG"
START_CODONS: Set[str] = {START_CODON}
STOP_CODONS = {"TAA", "TAG", "TGA"}

# A small, commonly used set of restriction enzymes and their recognition sites.
RESTRICTION_ENZYMES: Dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "PstI": "CTGCAG",
    "SalI": "GTCGAC",
    "XhoI": "CTCGAG",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SmaI": "CCCGGG",
    "XbaI": "TCTAGA",
    "NcoI": "CCATGG",
}

# Average residue masses (Da) for single-stranded DNA, plus the mass of a
# water molecule that is lost per phosphodiester bond.
_DNA_RESIDUE_MASS = {"A": 313.21, "T": 304.2, "C": 289.18, "G": 329.21}
_WATER_MASS = 18.015


# --------------------------------------------------------------------------
# Input handling
# --------------------------------------------------------------------------

def clean_sequence(sequence: str) -> str:
    """Return the sequence uppercased with all whitespace removed."""
    return "".join(sequence.split()).upper()


def validate_sequence(sequence: str, alphabet: Set[str] = DNA_ALPHABET) -> str:
    """Validate and normalise a sequence.

    Args:
        sequence: the raw sequence; case and whitespace are ignored.
        alphabet: the allowed set of characters. Defaults to ACGTN.

    Returns:
        The cleaned (uppercase, whitespace-free) sequence.

    Raises:
        ValueError: if the sequence is empty or contains characters that are
            not part of the allowed alphabet.
    """
    seq = clean_sequence(sequence)
    if not seq:
        raise ValueError("Sequence is empty.")
    invalid = sorted(set(seq) - alphabet)
    if invalid:
        raise ValueError(
            "Invalid character(s) in sequence: "
            + ", ".join(repr(c) for c in invalid)
            + f". Allowed characters: {''.join(sorted(alphabet))}."
        )
    return seq


# --------------------------------------------------------------------------
# Composition and properties
# --------------------------------------------------------------------------

def base_composition(sequence: str) -> Dict[str, object]:
    """Return nucleotide counts and fractions for a DNA sequence.

    Ambiguous bases (N) are counted separately and included in the length and
    fractions, so the fractions always sum to 1.
    """
    seq = validate_sequence(sequence)
    counts = {base: seq.count(base) for base in sorted(DNA_ALPHABET)}
    total = len(seq)
    fractions = {base: counts[base] / total for base in counts}
    return {"length": total, "counts": counts, "fractions": fractions}


def gc_content(sequence: str) -> float:
    """Return the GC content as a percentage (0-100).

    Ambiguous bases (N) are excluded from the denominator, because counting
    them as non-GC would understate the GC fraction of the known bases. A
    sequence with no unambiguous bases returns 0.0.
    """
    seq = validate_sequence(sequence)
    unambiguous = [b for b in seq if b in "ACGT"]
    if not unambiguous:
        return 0.0
    gc = sum(1 for b in unambiguous if b in "GC")
    return 100.0 * gc / len(unambiguous)


def molecular_weight(sequence: str) -> float:
    """Return the approximate molecular weight of single-stranded DNA, in Da.

    This is an approximation. Ambiguous bases (N) have no defined mass and
    contribute zero, so sequences containing N will be underestimated.
    """
    seq = validate_sequence(sequence)
    weight = sum(_DNA_RESIDUE_MASS.get(b, 0.0) for b in seq)
    return weight - _WATER_MASS * max(len(seq) - 1, 0)


def kmer_counts(sequence: str, k: int = 3) -> Dict[str, int]:
    """Return counts of every k-mer in the sequence.

    If k is longer than the sequence, an empty mapping is returned.

    Raises:
        ValueError: if k is not a positive integer.
    """
    if not isinstance(k, int) or k < 1:
        raise ValueError("k must be a positive integer.")
    seq = validate_sequence(sequence)
    counts: Dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        counts[kmer] = counts.get(kmer, 0) + 1
    return dict(sorted(counts.items()))


# --------------------------------------------------------------------------
# Transformation
# --------------------------------------------------------------------------

def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence.

    Accepts mixed case and ambiguous bases (N); N complements to N.
    """
    return clean_sequence(sequence).translate(_COMPLEMENT)[::-1]


def transcribe(sequence: str) -> str:
    """Return the RNA transcript of a DNA sequence (T -> U)."""
    return clean_sequence(sequence).translate(_TRANSCRIBE)


def translate(sequence: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate a DNA sequence into a protein string.

    Codons that are not in the standard table (for example ones containing N)
    are rendered as 'X', which is the conventional symbol for an unknown
    residue.

    Args:
        sequence: DNA sequence (may contain U, which is treated as T).
        frame: reading frame offset, 0, 1 or 2.
        stop_at_stop: if True, stop at the first stop codon; otherwise the
            stop codon is rendered as '*'.

    Raises:
        ValueError: if frame is not 0, 1 or 2.
    """
    if frame not in (0, 1, 2):
        raise ValueError("frame must be 0, 1 or 2.")
    seq = clean_sequence(sequence).replace("U", "T")
    protein: List[str] = []
    for i in range(frame, len(seq) - 2, 3):
        amino_acid = CODON_TABLE.get(seq[i:i + 3], "X")
        if amino_acid == "*" and stop_at_stop:
            break
        protein.append(amino_acid)
    return "".join(protein)


# --------------------------------------------------------------------------
# Feature finding
# --------------------------------------------------------------------------

def _orfs_on_strand(
    seq: str, strand: str, start_codons: Set[str]
) -> List[Dict[str, object]]:
    """Find ATG-to-stop ORFs on one strand, reported in original coordinates."""
    n = len(seq)
    hits: List[Dict[str, object]] = []
    for frame in range(3):
        i = frame
        while i + 3 <= n:
            if seq[i:i + 3] not in start_codons:
                i += 3
                continue
            j = i
            while j + 3 <= n:
                if seq[j:j + 3] in STOP_CODONS:
                    break
                j += 3
            else:
                j = None  # no in-frame stop codon before the end of the sequence
            if j is None:
                i = n
                continue
            end = j + 3  # exclusive, on the searched strand
            length_aa = (j - i) // 3
            if strand == "+":
                start1, end1 = i + 1, end
            else:
                # Map reverse-complement coordinates back to the original.
                start1, end1 = n - end + 1, n - i
            hits.append({
                "strand": strand,
                "frame": frame + 1,
                "start": start1,
                "end": end1,
                "length_aa": length_aa,
                "protein": translate(seq[i:end]),
            })
            i = j + 3
    return hits


def find_orfs(
    sequence: str,
    min_aa: int = 1,
    start_codons: Optional[Iterable[str]] = None,
) -> List[Dict[str, object]]:
    """Find open reading frames on both strands.

    An ORF runs from a start codon to the first in-frame stop codon. ORFs are
    non-overlapping within a reading frame: after a stop codon is found,
    scanning resumes after it.

    Args:
        sequence: DNA sequence.
        min_aa: minimum ORF length in amino acids (excluding the stop).
        start_codons: codons that may start an ORF. Defaults to ATG only.
            Bacterial sequences often use GTG or TTG as alternative starts.

    Returns:
        A list of dicts sorted by descending length, each with strand, frame,
        start, end, length_aa and protein.
    """
    if min_aa < 1:
        raise ValueError("min_aa must be >= 1.")
    starts = {clean_sequence(c) for c in (start_codons or START_CODONS)}
    bad = [c for c in starts if len(c) != 3]
    if bad:
        raise ValueError(f"Start codons must be 3 bases long: {bad}")
    seq = validate_sequence(sequence)
    hits = (
        _orfs_on_strand(seq, "+", starts)
        + _orfs_on_strand(reverse_complement(seq), "-", starts)
    )
    hits = [h for h in hits if h["length_aa"] >= min_aa]
    hits.sort(key=lambda h: (-int(h["length_aa"]), int(h["start"])))
    return hits


def find_restriction_sites(
    sequence: str,
    enzymes: Optional[Dict[str, str]] = None,
    circular: bool = False,
) -> Dict[str, List[int]]:
    """Return the 1-based start positions of restriction sites on the
    forward strand.

    Args:
        sequence: DNA sequence.
        enzymes: mapping of enzyme name to recognition site; defaults to
            RESTRICTION_ENZYMES.
        circular: if True, also detect sites that span the origin of a
            circular molecule. Positions are still reported within 1..len.
    """
    seq = validate_sequence(sequence)
    enzymes = enzymes or RESTRICTION_ENZYMES
    n = len(seq)
    result: Dict[str, List[int]] = {}
    for name, site in enzymes.items():
        site = clean_sequence(site)
        if not site:
            continue
        # For circular molecules, search a doubled sequence but only keep
        # positions that begin inside the original molecule.
        search_space = seq + seq[: len(site) - 1] if circular else seq
        positions: List[int] = []
        start = search_space.find(site)
        while start != -1:
            if start < n:
                positions.append(start + 1)
            start = search_space.find(site, start + 1)
        if positions:
            result[name] = sorted(set(positions))
    return result


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def read_fasta_string(text: str) -> Dict[str, str]:
    """Parse FASTA text into an ordered mapping of header -> sequence.

    Handles CRLF line endings, blank lines and lowercase sequence. Every
    record is validated, so an invalid character is reported together with the
    record it came from.

    Raises:
        ValueError: if no FASTA record is found, if sequence data appears
            before the first header, or if a record fails validation.
    """
    records: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []

    def flush() -> None:
        if name is None:
            return
        raw = "".join(chunks)
        try:
            records[name] = validate_sequence(raw)
        except ValueError as exc:
            raise ValueError(f"In record {name!r}: {exc}") from exc

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            name = line[1:].strip() or f"seq{len(records) + 1}"
            chunks = []
        else:
            if name is None:
                raise ValueError(
                    "Sequence data found before any '>' header. "
                    "A FASTA file must start with a header line."
                )
            chunks.append(line)
    flush()

    if not records:
        raise ValueError("No FASTA records found.")
    return records


def read_fasta(path: str) -> Dict[str, str]:
    """Read a FASTA file into a mapping of header -> sequence."""
    with open(path, "r", encoding="utf-8") as handle:
        return read_fasta_string(handle.read())


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

def summarize(sequence: str, min_orf_aa: int = 30) -> Dict[str, object]:
    """Return a compact report of the main statistics for a sequence."""
    seq = validate_sequence(sequence)
    comp = base_composition(seq)
    return {
        "length": comp["length"],
        "gc_content_percent": round(gc_content(seq), 2),
        "molecular_weight_da": round(molecular_weight(seq), 1),
        "base_counts": comp["counts"],
        "orf_count": len(find_orfs(seq, min_aa=min_orf_aa)),
        "restriction_sites": find_restriction_sites(seq),
    }
