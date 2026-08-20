"""The capability table replaces two figure panels, so it carries their claims.

Boundary snapping is the claim to watch. fastder snaps only the internal edges
of a multi-exon chain, never a monoexonic region and never the outer ends of a
chain, and the reviewer objected that the paper credited more to the mechanism
than it does. A bare "yes" in that cell would put the overclaim back.
"""

import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow" / "scripts" / "figures" / "make_capability_table.py"


def load_module():
    spec = importlib.util.spec_from_file_location("make_capability_table", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_every_tool_has_a_value_for_every_capability():
    m = load_module()
    for tool in m.TOOLS:
        assert len(m.CELLS[tool]) == len(m.COLUMNS), tool


def test_snapping_claim_stays_qualified():
    m = load_module()
    row = m.COLUMNS.index("Snaps boundaries to splice junctions")
    assert m.CELLS["fastder"][row] != "yes"
    assert "chain" in m.CELLS["fastder"][row]


def test_only_fastder_claims_multi_exon_output_and_strand():
    m = load_module()
    for capability in ("Emits multi-exon structure", "Assigns strand"):
        row = m.COLUMNS.index(capability)
        assert m.CELLS["fastder"][row].startswith("yes"), capability
        for tool in m.TOOLS:
            if tool != "fastder":
                assert not m.CELLS[tool][row].startswith("yes"), (tool, capability)


def test_written_files_match_the_declared_table(tmp_path):
    m = load_module()
    m.write_csv(tmp_path / "tool_capabilities.csv")
    m.write_tex(tmp_path / "tool_capabilities.tex")

    with open(tmp_path / "tool_capabilities.csv") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == len(m.TOOLS) * len(m.COLUMNS)

    tex = (tmp_path / "tool_capabilities.tex").read_text()
    for column in m.COLUMNS:
        assert column in tex
    for tool in m.TOOLS:
        assert tool in tex
    # One data row per capability, plus the header row.
    assert tex.count(r" \\") == len(m.COLUMNS) + 1
