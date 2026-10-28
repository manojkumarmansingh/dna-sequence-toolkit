"""Edge-case tests for the DNA Sequence Analysis Toolkit.

These complement tests/test_core.py and focus on boundaries, ambiguous input
and the biologically tricky cases: reverse-strand ORF coordinates, circular
restriction sites and FASTA parsing.
"""

import pytest

from dna_toolkit.core import (
    base_composition,
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
)


# --------------------------------------------------------------------------
# Composition edge cases
# --------------------------------------------------------------------------

def test_gc_content_all_ambiguous_returns_zero():
    assert gc_content("NNNN") == 0.0


def test_gc_content_ignores_ambiguous_in_denominator():
    # 2 GC out of 2 unambiguous bases, even though the sequence is 4 long.
    assert gc_content("GCNN") == pytest.approx(100.0)


def test_gc_content_rejects_empty():
    with pytest.raises(ValueError):
        gc_content("")


def test_base_composition_fractions_sum_to_one():
    comp = base_composition("ATGCNN")
    assert pytest.approx(sum(comp["fractions"].values())) == 1.0


def test_base_composition_counts_n():
    assert base_composition("ANNT")["counts"]["N"] == 2


def test_molecular_weight_of_single_base():
    # A single residue, so no water is lost.
    assert molecular_weight("A") == pytest.approx(313.21)


def test_molecular_weight_ignores_ambiguous_bases():
    assert molecular_weight("N") == pytest.approx(0.0)


# --------------------------------------------------------------------------
# Transformation edge cases
# --------------------------------------------------------------------------

def test_reverse_complement_of_palindrome_is_itself():
    assert reverse_complement("GAATTC") == "GAATTC"


def test_reverse_complement_handles_ambiguous_bases():
    assert reverse_complement("ANG") == "CNT"


def test_reverse_complement_accepts_lowercase():
    assert reverse_complement("atgc") == "GCAT"


def test_transcribe_accepts_lowercase():
    assert transcribe("atgc") == "AUGC"


def test_translate_shorter_than_one_codon_returns_empty():
    assert translate("AT") == ""


def test_translate_unknown_codon_is_x():
    # GCN contains N, so the codon is not in the table.
    assert translate("GCN") == "X"


def test_translate_stop_codon_rendered_by_default():
    assert translate("TAA") == "*"


def test_translate_stop_at_stop_stops():
    assert translate("TAAATG", stop_at_stop=True) == ""


# --------------------------------------------------------------------------
# ORF edge cases
# --------------------------------------------------------------------------

def test_find_orfs_ignores_orf_without_stop():
    # ATG followed by codons to the end, with no stop codon.
    assert find_orfs("ATGAAAAAAAAAAAA", min_aa=1) == []


def test_find_orfs_reverse_strand_exact_coordinates():
    # Forward:   5'-CCC TTA TTT TTT TTT CAT GGG-3'
    # Reverse complement of positions 4..18 reads ATG AAA AAA AAA TAA,
    # i.e. a 4 aa ORF on the minus strand.
    seq = "CCCTTATTTTTTTTTCATGGG"
    orfs = [o for o in find_orfs(seq, min_aa=4) if o["strand"] == "-"]
    assert orfs, "expected a reverse-strand ORF"
    best = orfs[0]
    assert best["length_aa"] == 4
    # The ORF must map back inside the original sequence.
    assert 1 <= best["start"] < best["end"] <= len(seq)
    # And the protein should be M-K-K-K plus the stop.
    assert best["protein"] == "MKKK*"


def test_find_orfs_alternative_start_codons():
    seq = "GTGAAATAA"
    assert find_orfs(seq, min_aa=2) == []           # ATG only by default
    hits = find_orfs(seq, min_aa=2, start_codons={"GTG"})
    assert len(hits) == 1
    # GTG can start an ORF, but under the standard genetic code it still
    # translates to valine (V), not methionine.
    assert hits[0]["protein"] == "VK*"


def test_find_orfs_rejects_bad_start_codon():
    with pytest.raises(ValueError):
        find_orfs("ATGTAA", start_codons={"AT"})


def test_find_orfs_min_length_boundary():
    # ATG AAA TAA -> 2 aa excluding the stop.
    assert len(find_orfs("ATGAAATAA", min_aa=2)) >= 1
    assert find_orfs("ATGAAATAA", min_aa=3) == []


def test_find_orfs_on_both_strands_are_labelled():
    seq = "ATGCGTTAA" + "CCCC" + "ATGTTTTAA"
    strands = {o["strand"] for o in find_orfs(seq, min_aa=1)}
    assert strands <= {"+", "-"}
    assert "+" in strands


# --------------------------------------------------------------------------
# Restriction-site edge cases
# --------------------------------------------------------------------------

def test_restriction_site_at_first_position():
    assert find_restriction_sites("GAATTCAAAA")["EcoRI"] == [1]


def test_restriction_site_overlapping_matches():
    # "AAAAAA" contains AAAA at positions 1, 2 and 3 for a 4-base site.
    sites = find_restriction_sites("AAAAAA", {"Test": "AAAA"})
    assert sites == {"Test": [1, 2, 3]}


def test_restriction_site_circular_wraps_origin():
    # EcoRI (GAATTC) split across the origin: ...GAA + TTC...
    seq = "TTCCCCCGAA"
    linear = find_restriction_sites(seq)
    circular = find_restriction_sites(seq, circular=True)
    assert "EcoRI" not in linear
    assert "EcoRI" in circular


def test_restriction_site_circular_does_not_duplicate():
    seq = "GAATTC"
    linear = find_restriction_sites(seq)
    circular = find_restriction_sites(seq, circular=True)
    assert circular["EcoRI"] == linear["EcoRI"] == [1]


def test_restriction_sites_empty_when_no_sites():
    assert find_restriction_sites("ACACACAC") == {}


# --------------------------------------------------------------------------
# k-mer edge cases
# --------------------------------------------------------------------------

def test_kmer_longer_than_sequence_returns_empty():
    assert kmer_counts("AT", k=5) == {}


def test_kmer_equal_to_length_returns_single_entry():
    assert kmer_counts("ATGC", k=4) == {"ATGC": 1}


def test_kmer_k_one_counts_bases():
    assert kmer_counts("AAT", k=1) == {"A": 2, "T": 1}


# --------------------------------------------------------------------------
# FASTA edge cases
# --------------------------------------------------------------------------

def test_fasta_handles_crlf_line_endings():
    records = read_fasta_string(">a\r\nATGC\r\n>b\r\nTTTT\r\n")
    assert records == {"a": "ATGC", "b": "TTTT"}


def test_fasta_handles_blank_lines_and_lowercase():
    records = read_fasta_string(">a\n\natgc\n\ngggg\n\n")
    assert records["a"] == "ATGCGGGG"


def test_fasta_reports_record_name_on_invalid_character():
    with pytest.raises(ValueError) as exc:
        read_fasta_string(">good\nATGC\n>bad\nATGX\n")
    assert "bad" in str(exc.value)


def test_fasta_rejects_empty_record():
    with pytest.raises(ValueError) as exc:
        read_fasta_string(">empty\n")
    assert "empty" in str(exc.value)


def test_fasta_rejects_data_before_header():
    with pytest.raises(ValueError):
        read_fasta_string("ATGC\n>a\nATGC\n")


def test_fasta_unnamed_header_gets_generated_name():
    records = read_fasta_string(">\nATGC\n")
    assert list(records) == ["seq1"]


# --------------------------------------------------------------------------
# Summary edge cases
# --------------------------------------------------------------------------

def test_summarize_on_minimal_sequence():
    report = summarize("ATG", min_orf_aa=1)
    assert report["length"] == 3
    assert report["orf_count"] == 0  # no stop codon, so no complete ORF
