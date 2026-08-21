"""The report palettes can fall behind the event classes a config simulates.

A sample missing from as_event_levels becomes NA in every plot faceted by
event class, and ComplexHeatmap stops on a colour it cannot map, which is how
job 11377514 lost its summary report after two hours of tool runs.
"""

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml", reason="PyYAML not installed in this env")

ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "config"
R_FILES = (
    ROOT / "workflow" / "reports" / "summary.Rmd",
    ROOT / "workflow" / "scripts" / "figures" / "figure_sim_event_jaccard.R",
)

VECTOR = re.compile(r"as_event_(levels|labels|palette) <- c\(")


def configured_samples():
    names = set()
    for path in sorted(CONFIG_DIR.glob("*.yaml")):
        config = yaml.safe_load(path.read_text()) or {}
        names |= set((config.get("asimulator") or {}).get("samples") or {})
    return names


def call_body(text, start):
    """The text between the c( at start and its matching paren.

    A label carries its own parentheses, as in "exon skip (ES)", so the first
    close paren is not the end of the call.
    """
    depth, quoted, i = 0, False, start
    while i < len(text):
        char = text[i]
        if char == '"':
            quoted = not quoted
        elif not quoted and char == "(":
            depth += 1
        elif not quoted and char == ")":
            depth -= 1
            if depth == 0:
                return text[start:i]
        i += 1
    raise AssertionError("unclosed c( at %d" % start)


def declared_keys(path):
    """Every as_event_* vector in one file, as name -> set of keys."""
    text = path.read_text()
    found = {}
    for match in VECTOR.finditer(text):
        body = call_body(text, match.end() - 1)
        if match.group(1) == "levels":
            found[match.group(1)] = set(re.findall(r'"([^"]+)"', body))
        else:
            found[match.group(1)] = set(re.findall(r"(\w+)\s*=", body))
    return found


def test_configs_declare_samples():
    assert configured_samples()


@pytest.mark.parametrize("path", R_FILES, ids=lambda p: p.name)
def test_every_configured_sample_has_a_level(path):
    for name, keys in declared_keys(path).items():
        assert configured_samples() <= keys, (name, sorted(configured_samples() - keys))


@pytest.mark.parametrize("path", R_FILES, ids=lambda p: p.name)
def test_the_vectors_in_one_file_agree(path):
    found = declared_keys(path)
    assert found, path
    assert len(set(map(frozenset, found.values()))) == 1, {k: sorted(v) for k, v in found.items()}
