"""fastder parameter grid: identifiers, CLI flags, and redundancy collapsing.

Lives here rather than in the Snakefile so it can be unit tested. The Snakefile
imports it and keeps the config handling.

A combination is a dict of parameter name to value. Two kinds of parameter
exist: values, passed as `--flag value`, and switches, passed as `--flag` when
true and omitted when false.
"""

import itertools
import re

# Parameter name to (identifier abbreviation, CLI flag).
PARAM_SPEC = {
    "min_coverage":       ("mc", "--min-coverage"),
    "min_length":         ("ml", "--min-length"),
    "position_tolerance": ("pt", "--position-tolerance"),
    "coverage_tolerance": ("ct", "--coverage-tolerance"),
    "min_junction_reads": ("mjr", "--min-junction-reads"),
    "no_stitch":          ("ns", "--no-stitch"),
}

# Switches rather than values. The identifier still records both states, so a
# switched-on and a switched-off run land in separate directories.
FLAG_PARAMS = {"no_stitch"}

# Parameters each switch makes inert, keyed by the switch. With --no-stitch
# nothing is joined across junctions, so neither tolerance can act.
IGNORED_WHEN_SET = {
    "no_stitch": {"position_tolerance", "coverage_tolerance"},
}


def param_id(combo):
    """Short, directory-safe identifier for a parameter combination."""
    parts = [f"{PARAM_SPEC[k][0]}{int(bool(v)) if k in FLAG_PARAMS else v}"
             for k, v in combo.items()]
    return "_".join(parts) if parts else "default"


def fastder_cli_args(combo):
    """CLI flags for this combination, covering only explicitly-set parameters."""
    args = []
    for k, v in combo.items():
        flag = PARAM_SPEC[k][1]
        if k in FLAG_PARAMS:
            if v:
                args.append(flag)
        else:
            args.append(f"{flag} {v}")
    return " ".join(args)


def drop_ignored(combo):
    """Remove parameters this combination cannot act on.

    Leaving them in would produce several identical runs per parameter point,
    wasting the downstream evaluation jobs and over-weighting the affected
    condition in any figure that aggregates over the grid.
    """
    ignored = set()
    for switch, names in IGNORED_WHEN_SET.items():
        if combo.get(switch):
            ignored |= names
    if not ignored:
        return combo
    return {k: v for k, v in combo.items() if k not in ignored}


def parse_param_id(param_id):
    """Recover a combination from its identifier.

    The inverse of param_id(), so a table carrying only param_id can be grouped
    by the axis that moved. A parameter absent from the identifier is absent
    here, which is not the same as zero. figures/helpers.R matches the same
    identifiers with regular expressions; keep the two in step.
    """
    combo = {}
    for name, (abbrev, _) in PARAM_SPEC.items():
        match = re.search(rf"(?:^|_){abbrev}([0-9]*\.?[0-9]+)(?:_|$)", param_id)
        if not match:
            continue
        raw = match.group(1)
        if name in FLAG_PARAMS:
            combo[name] = bool(int(float(raw)))
        elif "." in raw:
            combo[name] = float(raw)
        else:
            combo[name] = int(raw)
    return combo


def build_combos(axes_config):
    """Cross-product of the configured value lists, with redundancy collapsed.

    axes_config maps parameter names to lists of values; names absent from
    PARAM_SPEC are ignored. An empty mapping yields a single default combination.
    """
    axes = [[(key, val) for val in axes_config[key]]
            for key in PARAM_SPEC
            if key in axes_config]
    raw = [dict(items) for items in itertools.product(*axes)] if axes else [{}]

    combos = []
    seen = set()
    for combo in (drop_ignored(c) for c in raw):
        key = tuple(sorted(combo.items()))
        if key not in seen:
            seen.add(key)
            combos.append(combo)
    return combos
