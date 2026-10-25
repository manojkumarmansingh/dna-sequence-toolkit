"""Tests for the DNA Sequence Analysis Toolkit core module."""

import math

import pytest

from dna_toolkit.core import (
    base_composition,
    clean_sequence,
    find_orfs,
    find_restriction_sites,
    gc_content,
    kmer_counts,
    molecular_weight,
    read_fasta_string,
    reverse_complement,
    summarize,
    transcribe,
    translate,
    validate_sequence,
)


# --------------------------------------------------------------------------
# Input handling
# --------------------------------------------------------------------------

def test_clean_sequence_strips_whitespace_and_uppercases():
    assert clean_sequence(" at\ngc c ") == "ATGCC"


def test_validate_sequence_rejects_invalid_characters():
    with pytest.raises(ValueError):
        validate_sequence("ATGCX")


def test_validate_sequence_rejects_empty():
    with pytest.raises(ValueError):
        validate_sequence("   ")


def test_validate_sequence_accepts_n():
    assert validate_sequence("ATNGCN") == "ATNGCN"


# --------------------------------------------------------------------------
# Composition
# --------------------------------------------------------------------------

def test_gc_content_basic():
    assert gc_content("GCGC") == pytest.approx(100.0)
    assert gc_content("ATAT") == pytest.approx(0.0)
    assert gc_content("ATGC") == pytest.approx(50.0)


def test_gc_content_ignores_n():
    assert gc_content("GCNN") == pytest.approx(100.0)


def test_base_composition():
    comp = base_composition("AATTGGCC")
    assert comp["length"] == 8
    assert comp["counts"]["A"] == 2
    assert comp["counts"]["C"] == 2
    assert comp["fractions"]["G"] == pytest.approx(0.25)


def test_molecular_weight_is_positive():
    assert molecular_weight("ATGCATGC") > 0


# --------------------------------------------------------------------------
# Transformation
# --------------------------------------------------------------------------

def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAA") == "TTTT"


def test_reverse_complement_of_reverse_complement_is_identity():
    seq = "ATGCGTACGTTAGC"
    assert reverse_complement(reverse_complement(seq)) == seq


def test_transcribe():
    assert transcribe("ATGC") == "AUGC"


def test_translate_standard_codons():
    # ATG = M, TTT = F, GGG = G, then a stop codon.
    assert translate("ATGTTTGGGTAA") == "MFG*"
    assert translate("ATGTTTGGGTAA", stop_at_stop=True) == "MFG"


def test_translate_frame():
    assert translate("AATGTTTGGGTAA", frame=1, stop_at_stop=True) == "MFG"


def test_translate_rejects_bad_frame():
    with pytest.raises(ValueError):
        translate("ATGC", frame=3)


# --------------------------------------------------------------------------
# ORF finding
# --------------------------------------------------------------------------

def test_find_orfs_finds_forward_orf():
    seq = "CCC" + "ATG" + "AAA" * 4 + "TAA" + "GGG"
    orfs = find_orfs(seq, min_aa=1)
    forward = [o for o in orfs if o["strand"] == "+"]
    assert forward, "expected at least one forward ORF"
    best = forward[0]
    # ATG + 4 x AAA, then a stop: M K K K K.
    assert best["length_aa"] == 5
    assert best["start"] == 4
    assert best["protein"] == "M" + "K" * 4 + "*"


def test_find_orfs_respects_min_length():
    seq = "ATG" + "AAA" + "TAA"  # ATG + AAA = a 2 aa ORF
    assert find_orfs(seq, min_aa=3) == []
    assert len(find_orfs(seq, min_aa=2)) >= 1


def test_find_orfs_reverse_strand_coordinates():
    # An ORF on the reverse strand should map back to valid coordinates.
    seq = "CCC" + "TTA" + "TTT" * 4 + "CAT" + "GGG"  # revcomp contains ATG...(stop)
    orfs = find_orfs(seq, min_aa=1)
    for orf in orfs:
        assert 1 <= orf["start"] <= len(seq)
        assert orf["start"] <= orf["end"] <= len(seq)


# --------------------------------------------------------------------------
# Restriction sites
# --------------------------------------------------------------------------

def test_find_restriction_sites():
    seq = "AAAGAATTCAAAA"
    sites = find_restriction_sites(seq)
    assert sites["EcoRI"] == [4]


def test_find_restriction_sites_none():
    assert find_restriction_sites("AAAAAAAA") == {}


def test_find_restriction_sites_custom_enzymes():
    sites = find_restriction_sites("GGTACC", {"KpnI": "GGTACC"})
    assert sites == {"KpnI": [1]}


# --------------------------------------------------------------------------
# k-mers
# --------------------------------------------------------------------------

def test_kmer_counts():
    counts = kmer_counts("ATAT", k=2)
    assert counts == {"AT": 2, "TA": 1}


def test_kmer_counts_rejects_bad_k():
    with pytest.raises(ValueError):
        kmer_counts("ATGC", k=0)


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def test_read_fasta_string_multiple_records():
    text = ">seq1\nATGC\n>seq2\nGGGG\nTTTT\n"
    records = read_fasta_string(text)
    assert records["seq1"] == "ATGC"
    assert records["seq2"] == "GGGGTTTT"


def test_read_fasta_string_rejects_empty():
    with pytest.raises(ValueError):
        read_fasta_string("no headers here")


# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------

def test_summarize_keys():
    report = summarize("ATGAAATAA", min_orf_aa=1)
    assert report["length"] == 9
    assert math.isclose(report["gc_content_percent"], 11.11, abs_tol=0.01)
    assert set(report) >= {
        "length",
        "gc_content_percent",
        "molecular_weight_da",
        "base_counts",
        "orf_count",
        "restriction_sites",
    }
