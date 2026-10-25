"""Core sequence-analysis functions for the DNA Sequence Analysis Toolkit.

Implemented with the Python standard library only, so the package runs
anywhere without third-party dependencies.

All coordinate values returned by this module are 1-based and inclusive,
following the convention used in biology.
"""

from __future__ import annotations

from typing import Dict, Iterator, List, Tuple

__all__ = [
    "DNA_ALPHABET",
    "STOP_CODONS",
    "START_CODON",
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


def validate_sequence(sequence: str, alphabet: set = DNA_ALPHABET) -> str:
    """Validate and normalise a sequence.

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
    """Return nucleotide counts and fractions for a DNA sequence."""
    seq = validate_sequence(sequence)
    counts = {base: seq.count(base) for base in sorted(DNA_ALPHABET)}
    total = len(seq)
    fractions = {base: counts[base] / total for base in counts}
    return {"length": total, "counts": counts, "fractions": fractions}


def gc_content(sequence: str) -> float:
    """Return the GC content as a percentage (0-100).

    Ambiguous bases (N) are excluded from the denominator.
    """
    seq = validate_sequence(sequence)
    unambiguous = [b for b in seq if b in "ACGT"]
    if not unambiguous:
        return 0.0
    gc = sum(1 for b in unambiguous if b in "GC")
    return 100.0 * gc / len(unambiguous)


def molecular_weight(sequence: str) -> float:
    """Return the approximate molecular weight of a single-stranded DNA
    sequence in Daltons."""
    seq = validate_sequence(sequence)
    weight = sum(_DNA_RESIDUE_MASS.get(b, 0.0) for b in seq)
    return weight - _WATER_MASS * max(len(seq) - 1, 0)


def kmer_counts(sequence: str, k: int = 3) -> Dict[str, int]:
    """Return counts of every k-mer in the sequence.

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
    """Return the reverse complement of a DNA sequence."""
    return clean_sequence(sequence).translate(_COMPLEMENT)[::-1]


def transcribe(sequence: str) -> str:
    """Return the RNA transcript of a DNA sequence (T -> U)."""
    return clean_sequence(sequence).translate(_TRANSCRIBE)


def translate(sequence: str, frame: int = 0, stop_at_stop: bool = False) -> str:
    """Translate a DNA sequence into a protein string.

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

def _orfs_on_strand(seq: str, strand: str) -> List[Dict[str, object]]:
    """Find all ATG-to-stop ORFs on one strand, in original coordinates."""
    n = len(seq)
    hits: List[Dict[str, object]] = []
    for frame in range(3):
        i = frame
        while i + 3 <= n:
            if seq[i:i + 3] != START_CODON:
                i += 3
                continue
            j = i
            while j + 3 <= n:
                if seq[j:j + 3] in STOP_CODONS:
                    break
                j += 3
            else:
                j = None  # no in-frame stop codon before the end
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


def find_orfs(sequence: str, min_aa: int = 1) -> List[Dict[str, object]]:
    """Find open reading frames (ATG ... stop) on both strands.

    Args:
        sequence: DNA sequence.
        min_aa: minimum ORF length in amino acids (excluding the stop).

    Returns:
        A list of dicts sorted by descending length, each with strand, frame,
        start, end, length_aa and protein.
    """
    if min_aa < 1:
        raise ValueError("min_aa must be >= 1.")
    seq = validate_sequence(sequence)
    hits = _orfs_on_strand(seq, "+") + _orfs_on_strand(reverse_complement(seq), "-")
    hits = [h for h in hits if h["length_aa"] >= min_aa]
    hits.sort(key=lambda h: (-int(h["length_aa"]), int(h["start"])))
    return hits


def find_restriction_sites(
    sequence: str, enzymes: Dict[str, str] | None = None
) -> Dict[str, List[int]]:
    """Return the 1-based start positions of restriction sites on the
    forward strand.

    Args:
        sequence: DNA sequence.
        enzymes: mapping of enzyme name to recognition site; defaults to
            RESTRICTION_ENZYMES.
    """
    seq = validate_sequence(sequence)
    enzymes = enzymes or RESTRICTION_ENZYMES
    result: Dict[str, List[int]] = {}
    for name, site in enzymes.items():
        site = clean_sequence(site)
        positions = []
        start = seq.find(site)
        while start != -1:
            positions.append(start + 1)
            start = seq.find(site, start + 1)
        if positions:
            result[name] = positions
    return result


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def read_fasta_string(text: str) -> Dict[str, str]:
    """Parse FASTA text into an ordered mapping of header -> sequence.

    Raises:
        ValueError: if no FASTA record is found.
    """
    records: Dict[str, str] = {}
    name: str | None = None
    chunks: List[str] = []

    def flush() -> None:
        if name is not None:
            records[name] = clean_sequence("".join(chunks))

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
                raise ValueError("Sequence data found before any '>' header.")
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
