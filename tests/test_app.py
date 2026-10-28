"""Smoke tests for the Streamlit interface.

These use Streamlit's own ``AppTest`` harness, so they exercise the real script
the way the browser would: build the page, click the buttons, switch tabs and
check that nothing raises. Skipped automatically when Streamlit is not
installed (the core test job only installs pytest).
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("streamlit")

from streamlit.testing.v1 import AppTest  # noqa: E402

import app as app_module  # noqa: E402

APP_PATH = Path(__file__).resolve().parents[1] / "app.py"
DEMO = "ATGGCTAAACGTACGTTAAGCTTCGAATTCGGATCCAAATAA"
SECTIONS = app_module.SECTIONS


def _button(at: AppTest, label: str):
    for button in at.button:
        if button.label == label:
            return button
    raise AssertionError(f"button {label!r} not found")


def _radio(at: AppTest, label: str):
    for radio in at.radio:
        if radio.label == label:
            return radio
    raise AssertionError(f"radio {label!r} not found")


def _fresh() -> AppTest:
    at = AppTest.from_file(str(APP_PATH), default_timeout=60)
    at.run()
    return at


# --------------------------------------------------------------------------
# Pure helpers
# --------------------------------------------------------------------------

def test_parse_input_bare_sequence():
    assert app_module.parse_input("atg cgt\n") == {"input": "ATGCGT"}


def test_parse_input_fasta():
    records = app_module.parse_input(">a\nATGC\n>b\nGGCC\n")
    assert records == {"a": "ATGC", "b": "GGCC"}


def test_parse_input_rejects_empty_and_invalid():
    with pytest.raises(ValueError):
        app_module.parse_input("   ")
    with pytest.raises(ValueError):
        app_module.parse_input("ATGCZ")


def test_inspect_input_states():
    assert app_module.inspect_input("")["state"] == "empty"
    assert app_module.inspect_input("ATGCZ")["state"] == "error"
    ok = app_module.inspect_input("ATGCATGC")
    assert ok["state"] == "ok" and ok["length"] == 8


def test_formatters():
    assert app_module.fmt_pct(38.095) == "38.1%"
    assert app_module.fmt_count(1234) == "1,234"
    assert app_module.fmt_weight(12252.9) == "12,253 Da"


def test_composition_and_kmer_rows():
    rows = app_module.composition_rows(DEMO)
    assert [r["Base"] for r in rows] == ["A", "C", "G", "N", "T"]
    assert sum(r["Count"] for r in rows) == len(DEMO)
    kmers = app_module.kmer_rows(DEMO, 3, top=5)
    assert len(kmers) == 5
    assert kmers == sorted(kmers, key=lambda r: -r["Count"])


def test_as_fasta_wraps_at_sixty():
    text = app_module.as_fasta("demo", "A" * 130)
    lines = text.strip().splitlines()
    assert lines[0] == ">demo"
    assert [len(line) for line in lines[1:]] == [60, 60, 10]


# --------------------------------------------------------------------------
# App behaviour
# --------------------------------------------------------------------------

def test_app_starts_without_a_sequence():
    at = _fresh()
    assert not at.exception
    assert "No sequence loaded" in " ".join(m.value for m in at.markdown)


def test_load_example_then_analyse():
    at = _fresh()
    _button(at, "Load example").click().run()
    assert not at.exception
    _button(at, "Analyse sequence").click().run()
    assert not at.exception
    assert any("Loaded" in m.value for m in at.markdown)
    assert at.text_area[0].value == DEMO


def test_every_section_renders():
    at = _fresh()
    at.text_area[0].set_value(DEMO).run()
    _button(at, "Analyse sequence").click().run()
    for section in SECTIONS:
        _radio(at, "Analysis tool").set_value(section).run()
        assert not at.exception, f"{section} raised {at.exception}"


def test_invalid_input_shows_friendly_error():
    at = _fresh()
    at.text_area[0].set_value("ATGCZZZ").run()
    _button(at, "Analyse sequence").click().run()
    assert not at.exception
    assert at.error
    assert "A, C, G, T and N" in at.error[0].value


def test_clear_resets_the_workspace():
    at = _fresh()
    _button(at, "Load example").click().run()
    _button(at, "Clear").click().run()
    assert not at.exception
    assert at.text_area[0].value == ""
    assert "No sequence loaded" in " ".join(m.value for m in at.markdown)


def test_multi_record_fasta_offers_record_selector():
    at = _fresh()
    at.text_area[0].set_value(">one\nATGCATGC\n>two\nGGCCGGCC\n").run()
    _button(at, "Analyse sequence").click().run()
    assert not at.exception
    labels = [box.label for box in at.selectbox]
    assert "Record" in labels


def test_sample_button_loads_and_analyses():
    at = _fresh()
    _button(at, "Use this sample").click().run()
    assert not at.exception
    assert at.text_area[0].value == DEMO
    assert any("Loaded" in m.value for m in at.markdown)
