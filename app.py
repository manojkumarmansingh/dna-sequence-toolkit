"""Streamlit web interface for the DNA Sequence Analysis Toolkit.

Run locally with:

    streamlit run app.py

The app is a thin presentation layer over ``dna_toolkit``: every number on
screen comes from the package itself, so the command line, the library and this
interface always give the same answers. No analysis is reimplemented here.

Layout
------
    Sidebar ............ sequence source (paste, FASTA upload, example), limits
    Workspace .......... monospaced input, Analyse / Load example / Clear
    Navigation ......... Overview, Composition, Transform, Translation,
                         ORFs, Restriction sites, k-mer analysis

The pure helpers at the top of this module contain no Streamlit calls, so they
can be imported and tested without a running app.
"""

from __future__ import annotations

import json
from typing import Dict, List, Sequence

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

# Scientific UI palette: deep navy, cyan accents, soft off-white panels.
NAVY = "#0B1F3A"
CYAN = "#00B8D4"
TEAL = "#0E9AA7"
INK = "#0F172A"
MUTED = "#5B6B7F"
BORDER = "#DCE6F0"
CANVAS = "#F6F9FC"
PANEL = "#FFFFFF"

BASE_COLORS = {
    "A": CYAN,
    "T": NAVY,
    "G": TEAL,
    "C": "#5B8DEF",
    "N": "#94A3B8",
}

SECTIONS = [
    "Overview",
    "Composition",
    "Transform",
    "Translation",
    "ORFs",
    "Restriction sites",
    "k-mer analysis",
]


# --------------------------------------------------------------------------
# Pure helpers (no Streamlit calls, safe to unit-test)
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


def inspect_input(text: str) -> Dict[str, object]:
    """Describe the current input for the live status line under the editor."""
    raw = (text or "").strip()
    if not raw:
        return {
            "state": "empty",
            "message": "No sequence yet. Paste DNA, upload a FASTA file or load an example.",
            "records": 0,
            "length": 0,
        }
    try:
        records = parse_input(text)
    except (ValueError, FileNotFoundError) as exc:
        message = str(exc)
        if not raw.startswith(">"):
            message = f"{message} Use only the letters A, C, G, T and N."
        return {"state": "error", "message": message, "records": 0, "length": 0}
    total = sum(len(seq) for seq in records.values())
    if len(records) > 1:
        message = f"FASTA detected: {len(records)} records, {total:,} bp in total."
    else:
        message = f"{total:,} bp, ready to analyse."
    return {"state": "ok", "message": message, "records": len(records), "length": total}


def fmt_pct(value: float) -> str:
    """Format a percentage with one decimal place."""
    return f"{value:.1f}%"


def fmt_count(value: int) -> str:
    """Format an integer with thousands separators."""
    return f"{int(value):,}"


def fmt_weight(value: float) -> str:
    """Format a molecular weight in daltons."""
    return f"{value:,.0f} Da"


def composition_rows(sequence: str) -> List[Dict[str, object]]:
    """Return one display row per base, with counts and percentages."""
    comp = base_composition(sequence)
    counts = comp["counts"]
    total = comp["length"]
    rows = []
    for base, count in counts.items():
        rows.append({
            "Base": base,
            "Count": count,
            "Percent": round(100.0 * count / total, 2) if total else 0.0,
            "Share": (100.0 * count / total) if total else 0.0,
        })
    return rows


def kmer_rows(sequence: str, k: int, top: int | None = None) -> List[Dict[str, object]]:
    """Return k-mer counts as display rows, most frequent first."""
    counts = kmer_counts(sequence, k=k)
    if not counts:
        return []
    total = sum(counts.values())
    rows = [
        {
            "k-mer": kmer,
            "Count": count,
            "Percent": round(100.0 * count / total, 2) if total else 0.0,
            "Share": (100.0 * count / total) if total else 0.0,
        }
        for kmer, count in counts.items()
    ]
    rows.sort(key=lambda r: (-int(r["Count"]), str(r["k-mer"])))
    return rows[:top] if top else rows


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
# Styling
# --------------------------------------------------------------------------

CSS = f"""
<style>
:root {{
  --navy: {NAVY};
  --cyan: {CYAN};
  --teal: {TEAL};
  --ink: {INK};
  --muted: {MUTED};
  --border: {BORDER};
  --panel: {PANEL};
}}

.stApp {{ background: {CANVAS}; }}
.block-container {{ padding-top: 1.4rem; padding-bottom: 3rem; max-width: 1180px; }}

h1, h2, h3, h4 {{ color: var(--navy); letter-spacing: -0.01em; }}

/* ---------- hero ---------- */
.dna-hero {{
  display: flex; flex-wrap: wrap; align-items: center; justify-content: space-between;
  gap: 14px; padding: 18px 22px; border-radius: 14px;
  background: linear-gradient(135deg, {NAVY} 0%, #143055 62%, #1B3F6B 100%);
  border: 1px solid #1B3F6B;
}}
.dna-hero-left {{ display: flex; align-items: center; gap: 14px; min-width: 0; }}
.dna-hero-title {{ color: #FFFFFF; font-size: 1.32rem; font-weight: 700; line-height: 1.2; }}
.dna-hero-sub {{ color: #A9C4DC; font-size: 0.86rem; margin-top: 3px; }}
.dna-hero-right {{ display: flex; align-items: center; gap: 8px; flex-wrap: wrap; }}

.chip {{
  display: inline-flex; align-items: center; gap: 6px;
  padding: 4px 10px; border-radius: 999px; font-size: 0.74rem;
  font-weight: 600; letter-spacing: 0.02em; white-space: nowrap;
}}
.chip-ver {{ background: rgba(255,255,255,0.10); color: #CFE6F5; border: 1px solid rgba(255,255,255,0.18); }}
.chip-on {{ background: rgba(0,184,212,0.16); color: #9FE9F7; border: 1px solid rgba(0,184,212,0.45); }}
.chip-off {{ background: rgba(255,255,255,0.08); color: #A9C4DC; border: 1px solid rgba(255,255,255,0.16); }}
.dot {{ width: 7px; height: 7px; border-radius: 50%; display: inline-block; }}
.dot-on {{ background: #35D0A5; box-shadow: 0 0 0 3px rgba(53,208,165,0.18); }}
.dot-off {{ background: #7E93A8; }}

/* ---------- cards ---------- */
.card {{
  background: var(--panel); border: 1px solid var(--border); border-radius: 12px;
  padding: 14px 16px; box-shadow: 0 1px 2px rgba(11,31,58,0.04);
}}
.card-grid {{
  display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
  gap: 12px; margin: 4px 0 6px 0;
}}
.metric-label {{
  font-size: 0.71rem; text-transform: uppercase; letter-spacing: 0.07em;
  color: var(--muted); font-weight: 600;
}}
.metric-value {{ font-size: 1.42rem; font-weight: 700; color: var(--navy); line-height: 1.25; margin-top: 4px; }}
.metric-note {{ font-size: 0.72rem; color: var(--muted); margin-top: 2px; }}

/* ---------- bars ---------- */
.bars {{ margin-top: 6px; }}
.bar-row {{
  display: grid; grid-template-columns: 22px 1fr 104px;
  align-items: center; gap: 10px; margin-bottom: 7px;
}}
.bar-label {{ font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 0.82rem; font-weight: 700; color: var(--navy); }}
.bar-track {{ background: #EDF3F9; border-radius: 6px; height: 9px; overflow: hidden; }}
.bar-fill {{ display: block; height: 100%; border-radius: 6px; }}
.bar-val {{ font-size: 0.76rem; color: var(--muted); text-align: right; white-space: nowrap; }}

/* ---------- misc ---------- */
.panel-note {{
  background: #F1F7FB; border-left: 3px solid {CYAN}; border-radius: 6px;
  padding: 9px 12px; font-size: 0.82rem; color: #33506B; margin: 6px 0 2px 0;
}}
.empty-state {{
  background: var(--panel); border: 1px dashed #C9D8E6; border-radius: 12px;
  padding: 26px 22px; text-align: center; color: var(--muted); font-size: 0.9rem;
}}
.empty-state strong {{ color: var(--navy); }}

/* ---------- widgets ---------- */
.stTextArea textarea {{
  font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 0.86rem; line-height: 1.55; background: #FFFFFF;
}}
.stButton > button, .stDownloadButton > button {{
  border-radius: 8px; border: 1px solid #D2E0EC; font-weight: 600; font-size: 0.86rem;
}}
.stButton > button[kind="primary"] {{ background: {NAVY}; border-color: {NAVY}; }}
.stButton > button[kind="primary"]:hover {{ background: #143055; border-color: #143055; }}
div[data-testid="stSidebar"] {{ background: #EFF4FA; border-right: 1px solid #DCE6F0; }}
div[data-testid="stSidebar"] h2, div[data-testid="stSidebar"] h3 {{ color: var(--navy); }}
div[role="radiogroup"] {{ flex-wrap: wrap; gap: 6px; }}
div[data-testid="stDataFrame"] {{ border-radius: 10px; }}
code {{ color: #143055; }}
</style>
"""


def inject_styles() -> None:
    """Load the stylesheet once per run."""
    st.markdown(CSS, unsafe_allow_html=True)


def hero(status_on: bool, status_text: str) -> None:
    """Render the navy header with wordmark, version and status chip."""
    dot = "dot-on" if status_on else "dot-off"
    chip = "chip-on" if status_on else "chip-off"
    st.markdown(
        f"""
<div class="dna-hero">
  <div class="dna-hero-left">
    <svg width="36" height="36" viewBox="0 0 40 40" fill="none" aria-hidden="true">
      <path d="M12 4 C24 12 24 28 12 36" stroke="{CYAN}" stroke-width="2.4" stroke-linecap="round"/>
      <path d="M28 4 C16 12 16 28 28 36" stroke="#7FE3F0" stroke-width="2.4" stroke-linecap="round"/>
      <path d="M14.6 10.5 H25.4 M17 16 H23 M17 20 H23 M14.6 25.5 H25.4"
            stroke="#8FB6D6" stroke-width="1.5" stroke-linecap="round" opacity="0.75"/>
    </svg>
    <div>
      <div class="dna-hero-title">DNA Sequence Toolkit</div>
      <div class="dna-hero-sub">Explore DNA sequences with interactive bioinformatics tools.</div>
    </div>
  </div>
  <div class="dna-hero-right">
    <span class="chip chip-ver">v{__version__}</span>
    <span class="chip {chip}"><span class="dot {dot}"></span>{status_text}</span>
  </div>
</div>
""",
        unsafe_allow_html=True,
    )


def metric_cards(cards: Sequence[Dict[str, str]]) -> None:
    """Render a responsive grid of metric cards."""
    cells = "".join(
        f"""
    <div class="card">
      <div class="metric-label">{c['label']}</div>
      <div class="metric-value">{c['value']}</div>
      <div class="metric-note">{c.get('note', '')}</div>
    </div>"""
        for c in cards
    )
    st.markdown(f'<div class="card-grid">{cells}</div>', unsafe_allow_html=True)


def bar_chart(rows: Sequence[Dict[str, object]], label_key: str, value_key: str,
              share_key: str, colors: Dict[str, str] | None = None,
              max_share: float | None = None) -> None:
    """Render a lightweight horizontal bar chart in pure HTML/CSS."""
    scale = max_share if max_share is not None else 100.0
    parts = []
    for row in rows:
        share = float(row[share_key])
        width = max(1.0, min(100.0, 100.0 * share / scale)) if scale else 1.0
        color = colors.get(str(row[label_key]), CYAN) if colors else CYAN
        parts.append(
            f"""
    <div class="bar-row">
      <span class="bar-label">{row[label_key]}</span>
      <span class="bar-track"><span class="bar-fill" style="width:{width:.1f}%;background:{color}"></span></span>
      <span class="bar-val">{fmt_count(row[value_key])} &middot; {fmt_pct(share)}</span>
    </div>"""
        )
    st.markdown(f'<div class="bars">{"".join(parts)}</div>', unsafe_allow_html=True)


def panel_note(text: str) -> None:
    """Render a soft note strip used for assumptions and limitations."""
    st.markdown(f'<div class="panel-note">{text}</div>', unsafe_allow_html=True)


def section_header(title: str, blurb: str) -> None:
    """Section title plus a one-line explanation of what the tool does."""
    st.markdown(f"#### {title}")
    st.caption(blurb)


def empty_state() -> None:
    """Friendly placeholder shown when no sequence is loaded."""
    st.markdown(
        """
<div class="empty-state">
  <strong>No sequence loaded yet.</strong><br/>
  Paste a DNA sequence or FASTA text in the workspace above, upload a FASTA file
  from the sidebar, or press <em>Load example</em> to try the tools.
</div>
""",
        unsafe_allow_html=True,
    )


# --------------------------------------------------------------------------
# Analysis sections
# --------------------------------------------------------------------------

def render_overview(seq: str, name: str) -> None:
    """Headline statistics and nucleotide composition."""
    section_header(
        "Overview",
        "A one-screen summary of the loaded sequence: length, GC content, "
        "approximate mass and open reading frames of 30 amino acids or longer.",
    )
    report = summarize(seq, min_orf_aa=30)
    metric_cards([
        {"label": "Length", "value": f"{fmt_count(report['length'])} bp", "note": f"Record: {name}"},
        {"label": "GC content", "value": fmt_pct(report["gc_content_percent"]),
         "note": "N excluded from denominator"},
        {"label": "Molecular weight", "value": fmt_weight(report["molecular_weight_da"]),
         "note": "Single-stranded approximation"},
        {"label": "ORFs", "value": fmt_count(report["orf_count"]),
         "note": "At least 30 aa, ATG start"},
    ])

    st.markdown("**Nucleotide composition**")
    rows = composition_rows(seq)
    bar_chart(rows, label_key="Base", value_key="Count", share_key="Share", colors=BASE_COLORS)

    panel_note(
        "Molecular weight is an approximation for single-stranded DNA. GC content "
        "excludes ambiguous bases (N) from the denominator, so the percentages "
        "describe known bases only."
    )

    st.download_button(
        "Download summary (JSON)",
        data=json.dumps(report, indent=2).encode("utf-8"),
        file_name=f"{name}_summary.json",
        mime="application/json",
    )


def render_composition(seq: str, name: str) -> None:
    """Full base composition table with GC content."""
    section_header(
        "Composition",
        "Base counts and fractions for A, C, G, T and any ambiguous N bases, "
        "with GC content calculated over unambiguous bases.",
    )
    gc = gc_content(seq)
    comp = base_composition(seq)
    metric_cards([
        {"label": "GC content", "value": fmt_pct(gc), "note": "Known bases only"},
        {"label": "Length", "value": f"{fmt_count(comp['length'])} bp", "note": "Including N"},
        {"label": "Ambiguous bases", "value": fmt_count(comp["counts"].get("N", 0)),
         "note": "Counted as N"},
    ])

    rows = composition_rows(seq)
    st.markdown("**Base composition**")
    bar_chart(rows, label_key="Base", value_key="Count", share_key="Share", colors=BASE_COLORS)
    st.dataframe(
        pd.DataFrame([{k: r[k] for k in ("Base", "Count", "Percent")} for r in rows]),
        hide_index=True,
        use_container_width=True,
    )
    st.download_button(
        "Download composition (CSV)",
        data=to_csv_bytes([{k: r[k] for k in ("Base", "Count", "Percent")} for r in rows]),
        file_name=f"{name}_composition.csv",
        mime="text/csv",
    )
    panel_note(
        "Fractions are calculated over the full length of the sequence, so they "
        "always sum to 100%. GC content uses only A, C, G and T in the denominator."
    )


def render_transform(seq: str, name: str) -> None:
    """Reverse complement and RNA transcript."""
    section_header(
        "Transform",
        "Reverse complement (the opposite strand, written 5' to 3') and the RNA "
        "transcript of the forward strand, with T replaced by U.",
    )
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
        use_container_width=True,
    )
    c2.download_button(
        "Download transcript (FASTA)",
        data=as_fasta(f"{name}_rna", rna).encode("utf-8"),
        file_name=f"{name}_rna.fasta",
        mime="text/plain",
        use_container_width=True,
    )
    panel_note(
        "Both strands are written 5' to 3'. Ambiguous bases complement to N and "
        "are preserved as N in the reverse complement."
    )


def render_translate(seq: str, name: str) -> None:
    """Six-frame capable translation in a chosen reading frame."""
    section_header(
        "Translation",
        "Translate the forward strand in one reading frame using the standard "
        "genetic code (translation table 1).",
    )
    c1, c2 = st.columns([1, 2])
    frame = c1.radio("Reading frame", options=[0, 1, 2], horizontal=True, index=0)
    stop_at_stop = c2.checkbox("Stop at the first stop codon", value=True)

    protein = translate(seq, frame=frame, stop_at_stop=stop_at_stop)
    stops = protein.count("*")
    metric_cards([
        {"label": "Residues", "value": fmt_count(len(protein)), "note": f"Frame {frame}"},
        {"label": "Stop codons shown", "value": fmt_count(stops),
         "note": "Disabled when stopping early"},
        {"label": "Unknown residues", "value": fmt_count(protein.count("X")),
         "note": "Codons with N or gaps"},
    ])

    st.markdown(f"**Protein, frame {frame}**")
    st.code(protein or "(no complete codon in this frame)", language=None)
    st.caption("Codons containing ambiguous bases are shown as X. `*` marks a stop codon.")
    st.download_button(
        "Download protein (text)",
        data=(protein + "\n").encode("utf-8"),
        file_name=f"{name}_frame{frame}.protein.txt",
        mime="text/plain",
    )
    panel_note(
        "Only the forward strand is translated here. Use the ORF finder for "
        "both strands, alternative start codons and coordinate mapping."
    )


def render_orfs(seq: str, name: str) -> None:
    """Open reading frames on both strands."""
    section_header(
        "Open reading frames",
        "Find start-to-stop ORFs on all three reading frames of both strands, "
        "with coordinates mapped back to the original sequence.",
    )
    c1, c2 = st.columns([1, 2])
    min_aa = c1.number_input(
        "Minimum length (amino acids)", min_value=1, value=30, step=5,
        help="Length excludes the stop codon.",
    )
    alternative = c2.checkbox(
        "Include alternative start codons (GTG, TTG)", value=False,
        help="Bacterial sequences often use GTG or TTG as start codons.",
    )
    starts = {"ATG", "GTG", "TTG"} if alternative else None
    orfs = find_orfs(seq, min_aa=int(min_aa), start_codons=starts)

    if not orfs:
        st.info(f"No ORFs of at least {int(min_aa)} aa were found in this sequence.")
        return

    longest = max(int(o["length_aa"]) for o in orfs)
    metric_cards([
        {"label": "ORFs found", "value": fmt_count(len(orfs)), "note": f"At least {int(min_aa)} aa"},
        {"label": "Longest", "value": f"{fmt_count(longest)} aa", "note": "Excluding stop"},
        {"label": "Start codons", "value": "ATG, GTG, TTG" if alternative else "ATG",
         "note": "Standard genetic code"},
    ])

    rows = orf_rows(orfs)
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.download_button(
        "Download ORFs (CSV)",
        data=to_csv_bytes(rows),
        file_name=f"{name}_orfs.csv",
        mime="text/csv",
    )
    panel_note(
        "Coordinates are 1-based and inclusive on the original sequence. ORFs are "
        "non-overlapping within a reading frame: scanning resumes after each stop."
    )


def render_digest(seq: str, name: str) -> None:
    """Restriction-site mapping."""
    section_header(
        "Restriction sites",
        "Positions of recognition sites for 12 common enzymes on the forward "
        "strand, with an option for circular molecules.",
    )
    circular = st.checkbox(
        "Treat as a circular molecule", value=False,
        help="Also detects sites that span the origin.",
    )
    sites = find_restriction_sites(seq, circular=circular)
    if not sites:
        st.info("No sites from the default enzyme set were found in this sequence.")
        return

    total = sum(len(v) for v in sites.values())
    metric_cards([
        {"label": "Enzymes with sites", "value": fmt_count(len(sites)), "note": "Of 12 in the set"},
        {"label": "Total sites", "value": fmt_count(total), "note": "Forward strand"},
        {"label": "Topology", "value": "Circular" if circular else "Linear", "note": "Affects origin-spanning sites"},
    ])

    rows = digest_rows(sites)
    st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)
    st.download_button(
        "Download restriction sites (CSV)",
        data=to_csv_bytes(rows),
        file_name=f"{name}_restriction_sites.csv",
        mime="text/csv",
    )
    panel_note(
        "Sites are reported on the forward strand only. DNA methylation and star "
        "activity are not modelled."
    )


def render_kmers(seq: str, name: str) -> None:
    """k-mer counting with a frequency chart."""
    section_header(
        "k-mer analysis",
        "Count every k-mer of a chosen length. Useful for quick composition "
        "checks and for spotting repeats.",
    )
    k = st.slider("k-mer size", min_value=1, max_value=6, value=3)
    rows = kmer_rows(seq, k)
    if not rows:
        st.info(f"No {k}-mers: the sequence is shorter than k.")
        return

    top = rows[:12]
    metric_cards([
        {"label": "Distinct k-mers", "value": fmt_count(len(rows)), "note": f"k = {k}"},
        {"label": "Total positions", "value": fmt_count(sum(int(r["Count"]) for r in rows)),
         "note": "Overlapping windows"},
        {"label": "Most frequent", "value": str(top[0]["k-mer"]),
         "note": f"{fmt_count(top[0]['Count'])} occurrences"},
    ])

    st.markdown(f"**Most frequent {k}-mers**")
    bar_chart(top, label_key="k-mer", value_key="Count", share_key="Share",
              max_share=float(top[0]["Share"]))
    st.dataframe(
        pd.DataFrame([{kk: r[kk] for kk in ("k-mer", "Count", "Percent")} for r in rows]),
        hide_index=True,
        use_container_width=True,
    )
    st.download_button(
        f"Download {k}-mer counts (CSV)",
        data=to_csv_bytes([{kk: r[kk] for kk in ("k-mer", "Count", "Percent")} for r in rows]),
        file_name=f"{name}_{k}mers.csv",
        mime="text/csv",
    )
    panel_note(
        "Counts include overlapping windows, so a sequence of length L contains "
        "L - k + 1 k-mers. Percentages are relative to that total."
    )


RENDERERS = {
    "Overview": render_overview,
    "Composition": render_composition,
    "Transform": render_transform,
    "Translation": render_translate,
    "ORFs": render_orfs,
    "Restriction sites": render_digest,
    "k-mer analysis": render_kmers,
}


# --------------------------------------------------------------------------
# Input workspace
# --------------------------------------------------------------------------

def activate_text(text: str) -> None:
    """Load a sequence into the editor and analyse it in one step."""
    st.session_state["seq_text"] = text
    try:
        st.session_state["records"] = parse_input(text)
        st.session_state["load_error"] = None
    except (ValueError, FileNotFoundError) as exc:
        st.session_state["records"] = None
        st.session_state["load_error"] = str(exc)


def _load_example() -> None:
    """Callback: put the short demo sequence in the editor."""
    if SAMPLES:
        st.session_state["seq_text"] = SAMPLES["short_demo"]["sequence"]
    st.session_state["records"] = None


def _clear_input() -> None:
    """Callback: empty the editor and forget the loaded sequence."""
    st.session_state["seq_text"] = ""
    st.session_state["records"] = None


def render_workspace() -> None:
    """Sequence editor with Analyse, Load example and Clear controls."""
    st.markdown("#### Sequence workspace")
    st.caption("Paste a raw sequence or FASTA text, or upload a file from the sidebar.")

    if "seq_text" not in st.session_state:
        st.session_state["seq_text"] = ""
    st.session_state.setdefault("records", None)

    st.text_area(
        "DNA sequence or FASTA",
        key="seq_text",
        height=150,
        placeholder="ATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA\n\n"
                    "or FASTA:\n>my_sequence\nATGGCTAAACGT...",
        label_visibility="collapsed",
    )

    c1, c2, c3 = st.columns([1, 1, 1])
    analyse = c1.button("Analyse sequence", type="primary", use_container_width=True)
    c2.button("Load example", on_click=_load_example, use_container_width=True)
    c3.button("Clear", on_click=_clear_input, use_container_width=True)

    if st.session_state.get("load_error"):
        st.error(f"Could not read that sequence. {st.session_state['load_error']}")

    info = inspect_input(st.session_state.get("seq_text", ""))
    if info["state"] == "error":
        st.error(f"Could not read that sequence. {info['message']}")
    elif info["state"] == "ok":
        st.caption(info["message"])
    else:
        st.caption(info["message"])

    if analyse:
        try:
            st.session_state["records"] = parse_input(st.session_state.get("seq_text", ""))
        except (ValueError, FileNotFoundError) as exc:
            st.session_state["records"] = None
            if info["state"] != "error":  # avoid repeating the inline message
                st.error(
                    f"Could not read that sequence. {exc} "
                    "Use only the letters A, C, G, T and N."
                )
        else:
            # Re-run so the header, record selector and navigation pick up the
            # freshly loaded sequence in the same interaction.
            st.rerun()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main() -> None:
    st.set_page_config(
        page_title="DNA Sequence Toolkit",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    inject_styles()

    with st.sidebar:
        st.markdown("### Load sequence")
        source = st.radio(
            "Where is your sequence?",
            options=["Sample", "Paste", "Upload FASTA"],
            index=0,
        )
        if source == "Sample":
            if SAMPLES:
                key = st.selectbox(
                    "Sample sequence",
                    options=list(SAMPLES),
                    format_func=lambda k: SAMPLES[k]["label"],
                )
                st.caption(SAMPLES[key]["note"])
                if st.button("Use this sample", use_container_width=True):
                    activate_text(SAMPLES[key]["sequence"])
                    st.rerun()
            else:
                st.warning("Sample sequences are not available in this install.")
        elif source == "Upload FASTA":
            uploaded = st.file_uploader("Upload a FASTA file", type=["fasta", "fa", "txt"])
            if uploaded is not None:
                try:
                    text = uploaded.getvalue().decode("utf-8", errors="replace")
                except Exception as exc:  # pragma: no cover
                    st.error(f"Could not read the uploaded file. {exc}")
                else:
                    if text != st.session_state.get("seq_text"):
                        activate_text(text)
                        st.rerun()
                    st.caption(f"Loaded {uploaded.name}.")
        else:
            st.caption("Paste your sequence in the workspace below.")

        st.divider()
        with st.expander("Scientific limitations", expanded=False):
            st.markdown(
                "- Translation uses the standard genetic code (table 1) only.\n"
                "- ORFs default to ATG as the only start codon.\n"
                "- Restriction sites are forward-strand only, with no methylation "
                "or star-activity modelling.\n"
                "- Molecular weight is an approximation for single-stranded DNA.\n"
                "- Sequences are validated as ACGTN only."
            )
        st.divider()
        st.markdown("### About")
        st.caption(
            "A small, dependency-free toolkit for common DNA sequence analyses. "
            "The same functions power the Python API, the command line and this app."
        )
        st.markdown(
            "[Repository](https://github.com/manojkumarmansingh/dna-sequence-toolkit) · "
            "[Deployment notes](https://github.com/manojkumarmansingh/dna-sequence-toolkit/blob/main/docs/DEPLOYMENT.md)"
        )

    records = st.session_state.get("records")
    if records:
        name = next(iter(records)) if len(records) == 1 else st.selectbox(
            "Record", options=list(records)
        )
        seq = records[name]
        hero(True, f"Loaded · {name} · {fmt_count(len(seq))} bp")
    else:
        name, seq = "", ""
        hero(False, "No sequence loaded")

    render_workspace()

    if not records:
        empty_state()
        return

    st.divider()
    choice = st.radio(
        "Analysis tool",
        options=SECTIONS,
        horizontal=True,
        label_visibility="collapsed",
    )
    RENDERERS[choice](seq, name)


if __name__ == "__main__":
    main()
