"""Command-line interface for the DNA Sequence Analysis Toolkit.

Run ``python -m dna_toolkit --help`` for usage.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Optional

from . import __version__
from .core import (
    base_composition,
    find_orfs,
    find_restriction_sites,
    gc_content,
    kmer_counts,
    molecular_weight,
    read_fasta,
    reverse_complement,
    summarize,
    transcribe,
    translate,
    validate_sequence,
)


def _load_sequence(args: argparse.Namespace) -> str:
    """Return a single sequence from --sequence, or the first FASTA record."""
    if args.sequence:
        return validate_sequence(args.sequence)
    records = read_fasta(args.input)
    if not records:
        raise ValueError("No sequences found in input file.")
    name = args.record
    if name is None:
        name = next(iter(records))
    if name not in records:
        raise ValueError(
            f"Record {name!r} not found. Available: {', '.join(records)}"
        )
    return validate_sequence(records[name])


def _print(label: str, value: object) -> None:
    print(f"{label}: {value}")


def cmd_stats(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    report = summarize(seq, min_orf_aa=args.min_orf_aa)
    if args.json:
        print(json.dumps(report, indent=2))
        return 0
    _print("Length (bp)", report["length"])
    _print("GC content (%)", report["gc_content_percent"])
    _print("Molecular weight (Da)", report["molecular_weight_da"])
    _print("Base counts", report["base_counts"])
    _print(f"ORFs (>= {args.min_orf_aa} aa)", report["orf_count"])
    sites = report["restriction_sites"]
    _print("Restriction sites", sites if sites else "none found")
    return 0


def cmd_gc(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    _print("GC content (%)", round(gc_content(seq), 2))
    comp = base_composition(seq)
    _print("Base counts", comp["counts"])
    return 0


def cmd_revcomp(args: argparse.Namespace) -> int:
    print(reverse_complement(_load_sequence(args)))
    return 0


def cmd_transcribe(args: argparse.Namespace) -> int:
    print(transcribe(_load_sequence(args)))
    return 0


def cmd_translate(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    _print("Frame", args.frame)
    _print("Protein", translate(seq, frame=args.frame, stop_at_stop=True))
    return 0


def cmd_orfs(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    orfs = find_orfs(seq, min_aa=args.min_aa)
    if args.json:
        print(json.dumps(orfs, indent=2))
        return 0
    if not orfs:
        print(f"No ORFs of at least {args.min_aa} aa found.")
        return 0
    print(f"{len(orfs)} ORF(s) of at least {args.min_aa} aa:")
    for orf in orfs:
        print(
            f"  strand {orf['strand']}  frame {orf['frame']}  "
            f"{orf['start']}-{orf['end']}  {orf['length_aa']} aa"
        )
        print(f"    {orf['protein']}")
    return 0


def cmd_digest(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    sites = find_restriction_sites(seq)
    if not sites:
        print("No restriction sites found for the default enzyme set.")
        return 0
    print("Restriction sites (1-based positions on the forward strand):")
    for enzyme, positions in sites.items():
        print(f"  {enzyme:8s} {positions}")
    return 0


def cmd_kmers(args: argparse.Namespace) -> int:
    seq = _load_sequence(args)
    counts = kmer_counts(seq, k=args.k)
    for kmer, count in counts.items():
        print(f"{kmer}\t{count}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dna_toolkit",
        description="DNA Sequence Analysis Toolkit - common sequence analyses.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_input(p: argparse.ArgumentParser) -> None:
        group = p.add_mutually_exclusive_group(required=True)
        group.add_argument("-s", "--sequence", help="DNA sequence given directly")
        group.add_argument("-i", "--input", help="path to a FASTA file")
        p.add_argument("--record", help="FASTA record name to use (default: first)")

    p_stats = sub.add_parser("stats", help="summary statistics for a sequence")
    add_input(p_stats)
    p_stats.add_argument("--min-orf-aa", type=int, default=30, help="minimum ORF length in amino acids")
    p_stats.add_argument("--json", action="store_true", help="output as JSON")
    p_stats.set_defaults(func=cmd_stats)

    p_gc = sub.add_parser("gc", help="GC content and base composition")
    add_input(p_gc)
    p_gc.set_defaults(func=cmd_gc)

    p_rc = sub.add_parser("revcomp", help="reverse complement")
    add_input(p_rc)
    p_rc.set_defaults(func=cmd_revcomp)

    p_tr = sub.add_parser("transcribe", help="DNA to RNA")
    add_input(p_tr)
    p_tr.set_defaults(func=cmd_transcribe)

    p_tl = sub.add_parser("translate", help="translate to protein")
    add_input(p_tl)
    p_tl.add_argument("--frame", type=int, choices=[0, 1, 2], default=0)
    p_tl.set_defaults(func=cmd_translate)

    p_orf = sub.add_parser("orfs", help="find open reading frames")
    add_input(p_orf)
    p_orf.add_argument("--min-aa", type=int, default=30, help="minimum ORF length in amino acids")
    p_orf.add_argument("--json", action="store_true", help="output as JSON")
    p_orf.set_defaults(func=cmd_orfs)

    p_dig = sub.add_parser("digest", help="map restriction sites")
    add_input(p_dig)
    p_dig.set_defaults(func=cmd_digest)

    p_k = sub.add_parser("kmers", help="k-mer counts")
    add_input(p_k)
    p_k.add_argument("-k", type=int, default=3)
    p_k.set_defaults(func=cmd_kmers)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
