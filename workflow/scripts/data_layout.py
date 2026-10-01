"""Directory keys for intermediate files that configs may share.

Lives here rather than in the Snakefile so it can be unit tested.
"""


def reads_key(seq_depth, seed):
    """Key of a set of simulated reads. int() maps 1e7 and 10000000 together."""
    return f"depth_{int(seq_depth)}_seed_{seed}"


def alignment_key(config_name, pump_source, seq_depth, seed, ref_scope):
    """Key of a set of alignments.

    Simulated reads are aligned once per read set and index, so configs that
    share both share the alignments. Other input has no such key and stays
    with its config.
    """
    if pump_source == "asimulator":
        return f"{reads_key(seq_depth, seed)}_{ref_scope}"
    return config_name
