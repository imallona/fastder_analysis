"""Add the chr prefix to the sequence names of a GTF or GFF3 file."""


def prefix_line(line):
    if line.startswith("#") or not line.strip():
        return line
    return line if line.startswith("chr") else "chr" + line


def prefix_chromosomes(src, dst):
    """Write src to dst line by line, with chr-prefixed sequence names."""
    with open(src) as source, open(dst, "w") as target:
        for line in source:
            target.write(prefix_line(line))
