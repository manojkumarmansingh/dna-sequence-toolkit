# DNA Sequence Analysis Toolkit

A small, dependency-free Python toolkit for common DNA sequence analyses:
composition and GC content, reverse complement, transcription, translation,
open-reading-frame (ORF) finding, restriction-site mapping, k-mer counting and
FASTA parsing.

It ships with both a Python API and a command-line interface, and the whole
package runs on the Python standard library alone — nothing to install beyond
Python itself.

[![tests](https://github.com/manojkumarmansingh/dna-sequence-toolkit/actions/workflows/tests.yml/badge.svg)](https://github.com/manojkumarmansingh/dna-sequence-toolkit/actions/workflows/tests.yml)

## Features

- **Composition** — nucleotide counts and fractions, GC content (with `N`
  handled correctly), approximate molecular weight.
- **Transformation** — reverse complement, DNA→RNA transcription, and
  translation in any reading frame.
- **ORF finding** — scans all three reading frames on both strands and reports
  coordinates mapped back to the original sequence.
- **Restriction mapping** — a built-in set of common enzymes (EcoRI, BamHI,
  HindIII, NotI, PstI, SalI, XhoI, KpnI, SacI, SmaI, XbaI, NcoI) plus support
  for your own.
- **FASTA** — reads single- or multi-record FASTA files.
- **CLI + library** — usable from the shell or imported into a notebook.
- **Tested** — 25 unit tests, run in CI across Python 3.9–3.12.

## Installation

```bash
git clone https://github.com/manojkumarmansingh/dna-sequence-toolkit.git
cd dna-sequence-toolkit
pip install -e .
```

Or just run it in place — there are no third-party runtime dependencies:

```bash
python -m dna_toolkit --help
```

## Command-line usage

All commands accept either a sequence directly (`-s`) or a FASTA file
(`-i`), and `--record NAME` to pick a specific record from a FASTA file.

### Summary statistics

```bash
python -m dna_toolkit stats -i examples/example.fasta --min-orf-aa 1
```

```
Length (bp): 42
GC content (%): 38.1
Molecular weight (Da): 12252.9
Base counts: {'A': 15, 'C': 8, 'G': 8, 'N': 0, 'T': 11}
ORFs (>= 1 aa): 1
Restriction sites: {'EcoRI': [25], 'BamHI': [31], 'HindIII': [18]}
```

Add `--json` for machine-readable output.

### GC content

```bash
python -m dna_toolkit gc -s ATGCGCAT
```

```
GC content (%): 50.0
Base counts: {'A': 2, 'C': 2, 'G': 2, 'N': 0, 'T': 2}
```

### Reverse complement

```bash
python -m dna_toolkit revcomp -s ATGCGT
```

```
ACGCAT
```

### Translation

```bash
python -m dna_toolkit translate -s ATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA
```

```
Frame: 0
Protein: MAKRTLSFEFGSK
```

### Open reading frames

```bash
python -m dna_toolkit orfs -i examples/example.fasta --min-aa 1
```

```
1 ORF(s) of at least 1 aa:
  strand +  frame 1  1-42  13 aa
    MAKRTLSFEFGSK*
```

### Restriction sites

```bash
python -m dna_toolkit digest -s GAATTCGGATCCAAGCTT
```

```
Restriction sites (1-based positions on the forward strand):
  EcoRI    [1]
  BamHI    [7]
  HindIII  [13]
```

### k-mer counts

```bash
python -m dna_toolkit kmers -s ATAT -k 2
```

```
AT	2
TA	1
```

## Python usage

```python
from dna_toolkit import gc_content, reverse_complement, translate, find_orfs

seq = "ATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA"

print(gc_content(seq))          # 38.095238095238095
print(reverse_complement(seq))  # TTATTTGGATCCGAATTCGAAGCTTAACGTACGTTTAGCCAT
print(translate(seq, stop_at_stop=True))  # MAKRTLSFEFGSK

for orf in find_orfs(seq, min_aa=1):
    print(orf["strand"], orf["start"], orf["end"], orf["protein"])
```

## API overview

| Function | Purpose |
| --- | --- |
| `clean_sequence(seq)` | Uppercase and strip whitespace |
| `validate_sequence(seq)` | Check the sequence is non-empty and valid DNA |
| `base_composition(seq)` | Counts and fractions per base |
| `gc_content(seq)` | GC content as a percentage |
| `molecular_weight(seq)` | Approximate single-strand DNA mass (Da) |
| `kmer_counts(seq, k)` | Counts of every k-mer |
| `reverse_complement(seq)` | Reverse complement |
| `transcribe(seq)` | DNA → RNA |
| `translate(seq, frame, stop_at_stop)` | Translate to protein |
| `find_orfs(seq, min_aa)` | ORFs on both strands |
| `find_restriction_sites(seq, enzymes)` | Restriction-site positions |
| `read_fasta(path)` / `read_fasta_string(text)` | Parse FASTA |
| `summarize(seq, min_orf_aa)` | One-call summary report |

Coordinates are **1-based and inclusive** throughout.

## Project structure

```
dna-sequence-toolkit/
├── dna_toolkit/
│   ├── __init__.py      # public API
│   ├── core.py          # sequence-analysis functions
│   ├── cli.py           # command-line interface
│   └── __main__.py      # enables `python -m dna_toolkit`
├── tests/
│   └── test_core.py     # 25 unit tests
├── examples/
│   └── example.fasta
├── .github/workflows/tests.yml
├── pyproject.toml
└── README.md
```

## Running the tests

```bash
pip install -r requirements-dev.txt
python -m pytest -q
```

## Notes and limitations

- ORF finding uses the standard genetic code (translation table 1) with `ATG`
  as the only start codon and `TAA`/`TAG`/`TGA` as stops.
- Restriction-site mapping is done on the forward strand only; it does not
  account for methylation or star activity.
- Molecular weight is an approximation for single-stranded DNA.

## License

Released under the MIT License. See [LICENSE](LICENSE).

## Author

Manoj Kumar Mansingh — [GitHub](https://github.com/manojkumarmansingh)
