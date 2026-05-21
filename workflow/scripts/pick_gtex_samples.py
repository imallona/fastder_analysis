"""Pick GTEx samples per tissue and emit YAML stanzas (or splice them
in-place) for the fastder-evaluation gtex configs.

Reads recount3 GTEx project metadata (the .recount_project.tsv files
already cached under workflow/data/recount3/ by the
recount3_fetch_metadata rule) and selects N samples for each requested
tissue. Selection is seeded so re-running yields the same set.

Two modes:

1. Print to stdout (default).  Prints YAML sub-group stanzas to standard
   output; useful for previewing.

2. --apply <config> [<config> ...].  Rewrites the recount3.groups block
   of each given config file in place.  Tissues already present in the
   config keep their existing sample IDs exactly (so already-downloaded
   BigWigs are not invalidated).  New tissues get a fresh deterministic
   pick from metadata.  Existing tissues stay in their original order in
   the file; new tissues are appended in --tissues order.

If a metadata TSV is missing for a tissue that needs picking, the script
prints the curl command (matching the recount3_fetch_metadata rule) and
exits.

Usage:
    python workflow/scripts/pick_gtex_samples.py \\
        --tissues BLOOD BRAIN HEART MUSCLE LIVER LUNG TESTIS ADIPOSE_TISSUE \\
        --metadata-dir workflow/data/recount3 \\
        --seed 10 --n-per-tissue 40 --subgroups 8 \\
        --apply config/config_gtex_comparison.yaml \\
                config/config_gtex_concordance.yaml
"""
import argparse
import csv
import io
import os.path as op
import random
import sys

import yaml


METADATA_URL_TEMPLATE = (
    "https://duffel.rail.bio/recount3/human/data_sources/gtex/metadata/"
    "{shard}/{study}/gtex.recount_project.{study}.MD.gz"
)


def _shard(study):
    return study[-2:]


def load_external_ids(tsv_path):
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [row["external_id"] for row in reader]


def pick_samples(samples, n, seed):
    if len(samples) < n:
        raise ValueError(f"only {len(samples)} samples available, need {n}")
    rng = random.Random(seed)
    return sorted(rng.sample(samples, n))


def emit_groups_block(groups_dict, indent_groups=2):
    """Render the recount3.groups block from a dict {sub_group_name:
    {'study': ..., 'samples': [...]}}, starting with '  groups:'."""
    pad_g = " " * indent_groups
    pad_i = " " * (indent_groups + 2)
    pad_f = " " * (indent_groups + 4)
    pad_s = " " * (indent_groups + 6)
    lines = [f"{pad_g}groups:"]
    for name, body in groups_dict.items():
        lines.append(f"{pad_i}{name}:")
        lines.append(f"{pad_f}study: {body['study']}")
        lines.append(f"{pad_f}samples:")
        for s in body["samples"]:
            lines.append(f"{pad_s}- {s}")
    return "\n".join(lines)


def emit_stanzas_only(groups_dict, indent_items=4):
    """Render only the sub-group stanzas (no 'groups:' header), for the
    print-to-stdout preview mode."""
    pad_i = " " * indent_items
    pad_f = " " * (indent_items + 2)
    pad_s = " " * (indent_items + 4)
    lines = []
    for name, body in groups_dict.items():
        lines.append(f"{pad_i}{name}:")
        lines.append(f"{pad_f}study: {body['study']}")
        lines.append(f"{pad_f}samples:")
        for s in body["samples"]:
            lines.append(f"{pad_s}- {s}")
    return "\n".join(lines)


def find_groups_block(text):
    """Return (start_char, end_char) of the recount3.groups block in
    `text`. start_char points to the first character of 'groups:' (with
    its leading indent); end_char points one past the last character of
    the block, i.e. to the next top-level (or 2-indent) key, or EOF."""
    lines = text.split("\n")
    in_recount3 = False
    groups_line = None
    for i, line in enumerate(lines):
        stripped = line.rstrip()
        if not stripped or stripped.lstrip().startswith("#"):
            continue
        leading = len(line) - len(line.lstrip())
        if leading == 0 and stripped.endswith(":"):
            in_recount3 = stripped.startswith("recount3:")
            continue
        if in_recount3 and leading == 2 and stripped.strip() == "groups:":
            groups_line = i
            break
    if groups_line is None:
        raise ValueError("could not find 'recount3.groups:' block")

    end_line = len(lines)
    for k in range(groups_line + 1, len(lines)):
        s = lines[k]
        if not s.strip() or s.lstrip().startswith("#"):
            continue
        leading = len(s) - len(s.lstrip())
        if leading <= 2:
            end_line = k
            break

    char_offsets = [0]
    for ln in lines:
        char_offsets.append(char_offsets[-1] + len(ln) + 1)
    return char_offsets[groups_line], char_offsets[end_line]


def parse_existing_groups(text):
    """Return {sub_group_name: {'study': ..., 'samples': [...]}} parsed
    from the recount3.groups block of `text`. Preserves the order in
    which sub-groups appear in the file."""
    start, end = find_groups_block(text)
    block_text = text[start:end]
    parsed = yaml.safe_load(block_text) or {}
    return parsed.get("groups") or {}


def generate_tissue_groups(tissue, n_per_tissue, subgroups, seed, metadata_dir):
    """Pick samples deterministically and split into sub-groups for one
    tissue. Returns an ordered dict {sub_group_name: body}."""
    tsv_path = op.join(metadata_dir, f"{tissue}.recount_project.tsv")
    samples = load_external_ids(tsv_path)
    picked = pick_samples(samples, n_per_tissue, seed)
    subgroup_size = n_per_tissue // subgroups
    tissue_lc = tissue.lower()
    out = {}
    for i in range(subgroups):
        block = picked[i * subgroup_size : (i + 1) * subgroup_size]
        out[f"{tissue_lc}_{i + 1}"] = {"study": tissue, "samples": block}
    return out


def build_merged_groups(existing, requested_tissues, n_per_tissue, subgroups,
                       seed, metadata_dir):
    """Existing tissues keep their original sample IDs (and their order
    in the file). New tissues get a fresh deterministic pick. Returns
    (merged_dict, added_tissues, preserved_tissues, missing_metadata)."""
    existing_by_tissue = {}
    for name in existing:
        if "_" not in name:
            continue
        tissue_lc = name.rsplit("_", 1)[0]
        existing_by_tissue.setdefault(tissue_lc, []).append(name)

    requested_lcs = {t.lower() for t in requested_tissues}
    new_tissues = [t for t in requested_tissues
                   if t.lower() not in existing_by_tissue]
    preserved_tissues = [t for t in requested_tissues
                         if t.lower() in existing_by_tissue]

    missing = [(t, op.join(metadata_dir, f"{t}.recount_project.tsv"))
               for t in new_tissues
               if not op.exists(op.join(metadata_dir,
                                        f"{t}.recount_project.tsv"))]
    if missing:
        return None, [], [], missing

    merged = {}
    for name, body in existing.items():
        merged[name] = body
    for tissue in new_tissues:
        new_groups = generate_tissue_groups(tissue, n_per_tissue, subgroups,
                                            seed, metadata_dir)
        merged.update(new_groups)
    return merged, new_tissues, preserved_tissues, []


def apply_to_config(config_path, requested_tissues, n_per_tissue, subgroups,
                    seed, metadata_dir):
    with open(config_path) as fh:
        text = fh.read()
    existing = parse_existing_groups(text)
    merged, added, preserved, missing = build_merged_groups(
        existing, requested_tissues, n_per_tissue, subgroups, seed,
        metadata_dir)
    if missing:
        return None, missing
    start, end = find_groups_block(text)
    new_block = emit_groups_block(merged)
    suffix = text[end:]
    if suffix and not suffix.startswith("\n"):
        new_block += "\n"
    new_text = text[:start] + new_block + ("\n" if not suffix.startswith("\n")
                                            else "") + suffix
    # Simpler: ensure exactly one newline between new_block and suffix.
    while new_text.endswith("\n\n\n"):
        new_text = new_text[:-1]
    with open(config_path, "w") as fh:
        fh.write(new_text)
    return (added, preserved), []


def report_missing_metadata(missing):
    print("missing metadata files:", file=sys.stderr)
    for tissue, path in missing:
        url = METADATA_URL_TEMPLATE.format(shard=_shard(tissue), study=tissue)
        print(f"  {path}", file=sys.stderr)
        print(f"    curl -fSL '{url}' | gunzip -c > '{path}'",
              file=sys.stderr)
    print("\nfetch them (one curl per tissue above, or trigger the "
          "recount3_fetch_metadata snakemake rule), then re-run.",
          file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tissues", nargs="+", required=True,
                    help="recount3 GTEx project names, uppercase "
                         "(e.g. BLOOD BRAIN HEART MUSCLE LIVER LUNG "
                         "TESTIS ADIPOSE_TISSUE)")
    ap.add_argument("--metadata-dir",
                    default=op.join("workflow", "data", "recount3"))
    ap.add_argument("--seed", type=int, default=10)
    ap.add_argument("--n-per-tissue", type=int, default=40)
    ap.add_argument("--subgroups", type=int, default=8)
    ap.add_argument("--apply", nargs="+", default=None,
                    help="config YAML file(s) to rewrite in place; "
                         "without this flag, the script prints stanzas "
                         "to stdout instead")
    args = ap.parse_args()

    if args.n_per_tissue % args.subgroups != 0:
        sys.exit(f"--n-per-tissue ({args.n_per_tissue}) must be divisible "
                 f"by --subgroups ({args.subgroups})")

    if args.apply:
        for cfg in args.apply:
            result, missing = apply_to_config(
                cfg, args.tissues, args.n_per_tissue, args.subgroups,
                args.seed, args.metadata_dir)
            if missing:
                report_missing_metadata(missing)
                sys.exit(1)
            added, preserved = result
            print(f"{cfg}: preserved {len(preserved)} tissues "
                  f"({', '.join(preserved) or '-'}); added "
                  f"{len(added)} tissues ({', '.join(added) or '-'})",
                  file=sys.stderr)
        return

    # No --apply: preview mode.  Show what *would* be merged into a fresh
    # config (no existing tissues).
    merged, added, preserved, missing = build_merged_groups(
        {}, args.tissues, args.n_per_tissue, args.subgroups, args.seed,
        args.metadata_dir)
    if missing:
        report_missing_metadata(missing)
        sys.exit(1)
    print(emit_stanzas_only(merged))


if __name__ == "__main__":
    main()
