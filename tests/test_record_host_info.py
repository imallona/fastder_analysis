"""Tests for record_host_info.py.

Covers a hyperthreaded two-core machine, both turbo files, and a machine that
exposes neither topology nor turbo.
"""

import csv

from record_host_info import (MISSING, cpu_topology, host_info, load_1min,
                              mem_total_kb, threads_per_core, turbo_state, write_tsv)


def cpuinfo(logical, cores, model="Test CPU @ 3.00GHz"):
    """One package with `cores` physical cores and `logical` processors."""
    blocks = []
    for processor in range(logical):
        blocks.append(f"processor\t: {processor}\n"
                      f"model name\t: {model}\n"
                      "physical id\t: 0\n"
                      f"core id\t\t: {processor % cores}\n")
    return "\n".join(blocks)


def test_hyperthreads_are_not_physical_cores():
    assert cpu_topology(cpuinfo(logical=4, cores=2)) == ("Test CPU @ 3.00GHz", 4, 2)
    assert threads_per_core(4, 2) == 2


def test_two_packages_count_their_cores_separately():
    text = ("processor: 0\nphysical id: 0\ncore id: 0\n\n"
            "processor: 1\nphysical id: 1\ncore id: 0\n")
    assert cpu_topology(text)[1:] == (2, 2)


def test_missing_topology_is_reported_as_missing():
    assert cpu_topology("processor\t: 0\nBogoMIPS\t: 50.00\n") == (MISSING, 1, 0)
    assert threads_per_core(1, 0) == MISSING


def test_intel_turbo_file_is_inverted(tmp_path):
    (tmp_path / "intel_pstate").mkdir()
    (tmp_path / "intel_pstate" / "no_turbo").write_text("0\n")
    assert turbo_state(str(tmp_path)) == "on"
    (tmp_path / "intel_pstate" / "no_turbo").write_text("1\n")
    assert turbo_state(str(tmp_path)) == "off"


def test_generic_boost_file(tmp_path):
    (tmp_path / "cpufreq").mkdir()
    (tmp_path / "cpufreq" / "boost").write_text("1\n")
    assert turbo_state(str(tmp_path)) == "on"


def test_turbo_is_missing_without_either_file(tmp_path):
    assert turbo_state(str(tmp_path)) == MISSING


def test_load_and_memory():
    assert load_1min("0.21 0.14 0.10 1/900 123\n") == "0.21"
    assert load_1min("") == MISSING
    assert mem_total_kb("MemTotal:       129234544 kB\nMemFree: 1 kB\n") == "129234544"
    assert mem_total_kb("") == MISSING


def test_table_round_trip(tmp_path):
    proc = tmp_path / "proc"
    proc.mkdir()
    (proc / "cpuinfo").write_text(cpuinfo(logical=4, cores=2))
    (proc / "loadavg").write_text("1.50 1.00 0.50 2/100 1\n")
    (proc / "meminfo").write_text("MemTotal:       1000 kB\n")
    rows = host_info(str(proc), str(tmp_path / "cpu"), environ={"SLURM_JOB_ID": "42"})
    out = tmp_path / "host_info.tsv"
    write_tsv(rows, out)

    with open(out, newline="") as fh:
        table = dict(csv.reader(fh, delimiter="\t"))
    assert table["field"] == "value"
    assert table["cpu_cores_physical"] == "2"
    assert table["threads_per_core"] == "2"
    assert table["turbo"] == MISSING
    assert table["load_1min"] == "1.50"
    assert table["slurm_job_id"] == "42"
    assert table["slurm_node"] == "none"
