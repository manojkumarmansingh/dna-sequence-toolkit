"""Sample sequences used by the Streamlit app and the tests.

These are short, synthetic sequences written to demonstrate the toolkit's
features. They are NOT real genes or real genomic records, and no biological
claim should be made about them.
"""

from typing import Dict

SAMPLES: Dict[str, Dict[str, str]] = {
    "short_demo": {
        "label": "Short demo (forward ORF + restriction sites)",
        "sequence": "ATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA",
        "note": "42 bp. Contains a forward ORF and EcoRI, BamHI and HindIII sites.",
    },
    "gc_rich": {
        "label": "GC-rich example",
        "sequence": "GCGCGCGGCCGCGCGCGCGCGGCCGCGCGC",
        "note": "28 bp. High GC content with two NotI sites.",
    },
    "gc_poor": {
        "label": "AT-rich example",
        "sequence": "ATATATAAATTTATATAAATTTATATAAA",
        "note": "28 bp. Low GC content, no sites from the default enzyme set.",
    },
    "ambiguous": {
        "label": "With ambiguous bases (N)",
        "sequence": "ATGGCNTAANNCGTACGNNN",
        "note": "20 bp. Shows how N is handled in GC content and translation.",
    },
    "multi_record": {
        "label": "Multi-record FASTA",
        "sequence": ">demo_orf\nATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA\n>gc_rich\nGCGCGCGGCCGCGCGCGCGCGGCCGCGCGC\n",
        "note": "Two records, for testing the FASTA parser and record selector.",
    },
}


def get_sample(key: str) -> str:
    """Return the raw sequence text for a named sample."""
    return SAMPLES[key]["sequence"]
