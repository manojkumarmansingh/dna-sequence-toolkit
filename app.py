"""Streamlit web interface for the DNA Sequence Analysis Toolkit.

Run locally with:

    streamlit run app.py

The app is a thin layer over ``dna_toolkit`` - all the analysis is done by the
package itself, so the command line, the library and this interface always give
the same answers.
"""

from __future__ import annotations

import io
import json
from typing import Dict, List

import pandas as pd
import streamlit as st

from dna_toolkit import __version__
from dna_toolkit.core import (
    RESTRICTION_ENZYMES,
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
    validate_sequence,
)

try:  # examples/ is shipped with the repository, not the installed package
    from examples.samples import SAMPLES
except ImportError:  # pragma: no cover - fallback keeps the app usable
    SAMPLES = {}

MAX_INPUT_CHARS = 200_000


# --------------------------------------------------------------------------
# Pure helpers (kept separate from the UI so they can be tested)
# --------------------------------------------------------------------------

def parse_input(text: str) -> Dict[str, str]:
    """Turn raw user input into a mapping of record name -> sequence.

    Input that starts with '>' is treated as FASTA; anything else is treated
    as a single bare sequence.
    """
    if not text or not text.strip():
        raise ValueError("Please provide a sequence or upload a FASTA file.")
    if len(text) > MAX_INPUT_CHARS:
        raise ValueError(
            f"Input is too large ({len(text):,} characters). "
            f"The limit is {MAX_INPUT_CHARS:,}."
        )
    if text.lstrip().startswith(">"):
        return read_fasta_string(text)
    return {"input": validate_sequence(text)}


def orf_rows(orfs: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """Shape ORF results for display in a table."""
    return [
        {
            "Strand": o["strand"],
            "Frame": o["frame"],
            "Start": o["start"],
            "End": o["end"],
            "Length (aa)": o["length_aa"],
            "Protein": o["protein"],
        }
        for o in orfs
    ]


def digest_rows(sites: Dict[str, List[int]]) -> List[Dict[str, object]]:
    """Shape restriction-site results for display in a table."""
    rows = []
    for enzyme, positions in sorted(sites.items()):
        rows.append({
            "Enzyme": enzyme,
            "Recognition site": RESTRICTION_ENZYMES.get(enzyme, ""),
            "Positions (1-based)": ", ".join(str(p) for p in positions),
            "Count": len(positions),
        })
    return rows


def to_csv_bytes(rows: List[Dict[str, object]]) -> bytes:
    """Serialise a list of row dicts as CSV bytes."""
    if not rows:
        return b""
    return pd.DataFrame(rows).to_csv(index=False).encode("utf-8")


def as_fasta(name: str, sequence: str, width: int = 60) -> str:
    """Format a sequence as a single FASTA record."""
    lines = [f">{name}"]
    for i in range(0, len(sequence), width):
        lines.append(sequence[i:i + width])
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# UI sections
# --------------------------------------------------------------------------

def render_overview(seq: str, name: str) -> None:
    st.subheader("Overview")
    report = summarize(seq, min_orf_aa=30)
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Length", f"{report['length']} bp")
    c2.metric("GC content", f"{report['gc_content_percent']}%")
    c3.metric("Molecular weight", f"{report['molecular_weight_da']:,.0f} Da")
    c4.metric("ORFs ≥ 30 aa", report["orf_count"])

    st.caption(
        "Molecular weight is an approximation for single-stranded DNA. "
        "GC content excludes ambiguous bases (N) from the denominator."
    )

    st.markdown("**Base composition**")
    comp = base_composition(seq)
    counts = comp["counts"]
    fractions = comp["fractions"]
    df = pd.DataFrame({
        "Base": list(counts),
        "Count": list(counts.values()),
        "Percent": [round(fractions[b] * 100, 2) for b in counts],
    })
    st.dataframe(df, hide_index=True, use_container_width=True)

    st.download_button(
        "Download summary (JSON)",
        data=json.dumps(report, indent=2).encode("utf-8"),
        file_name=f"{name}_summary.json",
        mime="application/json",
    )


def render_composition(seq: str, name: str) -> None:
    st.subheader("GC content and composition")
    st.metric("GC content", f"{gc_content(seq):.2f}%")
    comp = base_composition(seq)
    st.write(
        f"Length {comp['length']} bp; counts "
        + ", ".join(f"{b}={n}" for b, n in comp["counts"].items())
    )
    k = st.slider("k-mer size", min_value=1, max_value=6, value=3)
    counts = kmer_counts(seq, k=k)
    if counts:
        rows = [{"k-mer": kmer, "Count": c} for kmer, c in counts.items()]
        st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
        st.download_button(
            f"Download {k}-mer counts (CSV)",
            data=to_csv_bytes(rows),
            file_name=f"{name}_{k}mers.csv",
            mime="text/csv",
        )
    else:
        st.info(f"No {k}-mers: the sequence is shorter than k.")


def render_transform(seq: str, name: str) -> None:
    st.subheader("Reverse complement and transcription")
    rc = reverse_complement(seq)
    rna = transcribe(seq)
    st.markdown("**Reverse complement**")
    st.code(rc, language=None)
    st.markdown("**RNA transcript**")
    st.code(rna, language=None)
    c1, c2 = st.columns(2)
    c1.download_button(
        "Download reverse complement (FASTA)",
        data=as_fasta(f"{name}_revcomp", rc).encode("utf-8"),
        file_name=f"{name}_revcomp.fasta",
        mime="text/plain",
    )
    c2.download_button(
        "Download transcript (FASTA)",
        data=as_fasta(f"{name}_rna", rna).encode("utf-8"),
        file_name=f"{name}_rna.fasta",
        mime="text/plain",
    )


def render_translate(seq: str, name: str) -> None:
    st.subheader("Translation")
    frame = st.radio("Reading frame", options=[0, 1, 2], horizontal=True)
    stop_at_stop = st.checkbox("Stop at the first stop codon", value=True)
    protein = translate(seq, frame=frame, stop_at_stop=stop_at_stop)
    st.markdown(f"**Protein (frame {frame})**")
    st.code(protein or "(no complete codon in this frame)", language=None)
    st.caption("Codons containing ambiguous bases are shown as X.")
    st.download_button(
        "Download protein (text)",
        data=(protein + "\n").encode("utf-8"),
        file_name=f"{name}_frame{frame}.protein.txt",
        mime="text/plain",
    )


def render_orfs(seq: str, name: str) -> None:
    st.subheader("Open reading frames")
    min_aa = st.number_input("Minimum length (amino acids)", min_value=1, value=30, step=5)
    alternative = st.checkbox(
        "Include alternative start codons (GTG, TTG)", value=False,
        help="Bacterial sequences often use GTG or TTG as start codons.",
    )
    starts = {"ATG", "GTG", "TTG"} if alternative else None
    orfs = find_orfs(seq, min_aa=int(min_aa), start_codons=starts)

    if not orfs:
        st.info(f"No ORFs of at least {int(min_aa)} aa were found.")
        return

    rows = orf_rows(orfs)
    st.write(f"{len(orfs)} ORF(s) found, longest first.")
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.download_button(
        "Download ORFs (CSV)",
        data=to_csv_bytes(rows),
        file_name=f"{name}_orfs.csv",
        mime="text/csv",
    )


def render_digest(seq: str, name: str) -> None:
    st.subheader("Restriction sites")
    circular = st.checkbox(
        "Treat as a circular molecule", value=False,
        help="Also detects sites that span the origin.",
    )
    sites = find_restriction_sites(seq, circular=circular)
    if not sites:
        st.info("No sites from the default enzyme set were found.")
        return
    rows = digest_rows(sites)
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.download_button(
        "Download restriction sites (CSV)",
        data=to_csv_bytes(rows),
        file_name=f"{name}_restriction_sites.csv",
        mime="text/csv",
    )
    st.caption(
        "Sites are reported on the forward strand only and do not account for "
        "DNA methylation or star activity."
    )


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(
        page_title="DNA Sequence Analysis Toolkit",
        page_icon="🧬",
        layout="wide",
    )

    st.title("DNA Sequence Analysis Toolkit")
    st.caption(
        f"Version {__version__} · A small, dependency-free toolkit for common "
        "DNA sequence analyses."
    )

    with st.sidebar:
        st.header("Input")
        source = st.radio(
            "Where is your sequence?",
            options=["Sample", "Paste", "Upload FASTA"],
            index=0,
        )

        text = ""
        if source == "Sample":
            if SAMPLES:
                key = st.selectbox(
                    "Sample sequence",
                    options=list(SAMPLES),
                    format_func=lambda k: SAMPLES[k]["label"],
                )
                st.caption(SAMPLES[key]["note"])
                text = SAMPLES[key]["sequence"]
            else:
                st.warning("Sample sequences are not available in this install.")
        elif source == "Paste":
            text = st.text_area(
                "Paste a DNA sequence or FASTA",
                height=200,
                placeholder="ATGGCTAAACGT...\n\nor\n\n>my_sequence\nATGGCTAAACGT...",
            )
        else:
            uploaded = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "txt"])
            if uploaded is not None:
                try:
                    text = uploaded.getvalue().decode("utf-8", errors="replace")
                except Exception as exc:  # pragma: no cover
                    st.error(f"Could not read the uploaded file: {exc}")

        st.divider()
        with st.expander("Scientific limitations"):
            st.markdown(
                "- Translation uses the standard genetic code (table 1) only.\n"
                "- ORFs default to ATG as the only start codon.\n"
                "- Restriction sites are forward-strand only, with no methylation "
                "or star-activity modelling.\n"
                "- Molecular weight is an approximation for single-stranded DNA.\n"
                "- Sequences are validated as ACGTN only."
            )

    # ---- parse and validate ----
    try:
        records = parse_input(text)
    except (ValueError, FileNotFoundError) as exc:
        st.error(f"Could not read that sequence: {exc}")
        st.stop()

    if len(records) > 1:
        record_name = st.selectbox("Record", options=list(records))
    else:
        record_name = next(iter(records))

    seq = records[record_name]
    st.caption(f"Loaded **{record_name}** — {len(seq):,} bp")

    tabs = st.tabs(
        ["Overview", "Composition", "Transform", "Translate", "ORFs", "Restriction sites"]
    )
    with tabs[0]:
        render_overview(seq, record_name)
    with tabs[1]:
        render_composition(seq, record_name)
    with tabs[2]:
        render_transform(seq, record_name)
    with tabs[3]:
        render_translate(seq, record_name)
    with tabs[4]:
        render_orfs(seq, record_name)
    with tabs[5]:
        render_digest(seq, record_name)


if __name__ == "__main__":
    main()
