"""Machine a run's timings came from, as a field and value TSV.

A value the machine does not expose is written as NA.

Usage:
    python record_host_info.py --out <tsv>
"""
import argparse
import csv
import os
import os.path as op
import socket

MISSING = "NA"

# Intel's driver reports turbo inverted (1 means off); the generic cpufreq
# file reports it directly.
TURBO_FILES = (
    ("intel_pstate/no_turbo", {"0": "on", "1": "off"}),
    ("cpufreq/boost", {"1": "on", "0": "off"}),
)


def read_text(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ""


def cpu_topology(cpuinfo):
    """Model, logical cores and physical cores from /proc/cpuinfo text."""
    model = MISSING
    logical = 0
    physical = set()
    package = None
    for line in cpuinfo.splitlines():
        key, _, value = line.partition(":")
        key, value = key.strip(), value.strip()
        if key == "processor":
            logical += 1
        elif key == "model name" and model == MISSING:
            model = value
        elif key == "physical id":
            package = value
        elif key == "core id":
            physical.add((package, value))
    return model, logical, len(physical)


def threads_per_core(logical, physical):
    if not physical or logical % physical:
        return MISSING
    return logical // physical


def turbo_state(cpu_dir):
    for name, states in TURBO_FILES:
        state = states.get(read_text(op.join(cpu_dir, name)).strip())
        if state:
            return state
    return MISSING


def load_1min(loadavg):
    fields = loadavg.split()
    return fields[0] if fields else MISSING


def mem_total_kb(meminfo):
    for line in meminfo.splitlines():
        if line.startswith("MemTotal:"):
            return line.split()[1]
    return MISSING


def host_info(proc_dir="/proc", cpu_dir="/sys/devices/system/cpu", environ=os.environ):
    model, logical, physical = cpu_topology(read_text(op.join(proc_dir, "cpuinfo")))
    return [
        ("hostname", socket.gethostname()),
        ("cpu_model", model),
        ("cpu_cores_total", logical or MISSING),
        ("cpu_cores_physical", physical or MISSING),
        ("threads_per_core", threads_per_core(logical, physical)),
        # The scheduler affinity, which OMP_NUM_THREADS does not change.
        ("cpu_cores_available", len(os.sched_getaffinity(0))),
        ("turbo", turbo_state(cpu_dir)),
        ("load_1min", load_1min(read_text(op.join(proc_dir, "loadavg")))),
        ("mem_total_kb", mem_total_kb(read_text(op.join(proc_dir, "meminfo")))),
        ("slurm_job_id", environ.get("SLURM_JOB_ID", "none")),
        ("slurm_node", environ.get("SLURMD_NODENAME", "none")),
    ]


def write_tsv(rows, path):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    write_tsv(host_info(), args.out)


if __name__ == "__main__":
    main()
