"""Directory keys for intermediate files that configs may share.

Lives here rather than in the Snakefile so it can be unit tested.
"""
import hashlib
import json


def reads_key(design, seed):
    """Key of a set of simulated reads.

    `design` is the config's asimulator block. Depth and seed are spelled out
    for the reader; the digest covers every setting, so two configs share
    reads only when all of them agree. int() maps 1e7 and 10000000 together.
    """
    design = dict(design or {})
    depth = int(design.get("seq_depth", 0))
    design["seq_depth"] = depth
    digest = hashlib.sha1(json.dumps(design, sort_keys=True).encode()).hexdigest()[:8]
    return f"depth_{depth}_seed_{seed}_{digest}"


def sample_seed(run_seed, samples, sample):
    """Seed of one simulated sample.

    Under one shared seed the simulator draws the same genes and expression
    levels for every sample, so samples differ only by event class. Each
    sample gets its own seed, spaced so that runs under consecutive run seeds
    never share one.
    """
    return int(run_seed) * 1000 + list(samples).index(sample)


def index_key(ref_scope, annotated):
    """Key of a STAR index: its chromosomes, and whether it holds the annotation."""
    return ref_scope if annotated else f"{ref_scope}_unannotated"


def alignment_key(config_name, pump_source, design, seed, index):
    """Key of a set of alignments. `index` is an index_key().

    Simulated reads are aligned once per read set and index, so configs that
    share both share the alignments. Other input has no such key and stays
    with its config.
    """
    if pump_source == "asimulator":
        return f"{reads_key(design, seed)}_{index}"
    return config_name
