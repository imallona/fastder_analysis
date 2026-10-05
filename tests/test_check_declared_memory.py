from check_declared_memory import declared_memory, excesses, measured_memory

RULES = """rule small_job:
    resources:
        mem_mb=2000,
        runtime=30,

rule small_job_scaling:
    resources:
        mem_mb=8000,
        runtime=30,

rule from_config:
    resources:
        mem_mb=SOME_CONSTANT,
"""

HEADER = "s\th:m:s\tmax_rss\tmax_vms\n"


def write_tree(tmp_path, peaks):
    rules = tmp_path / "rules"
    rules.mkdir()
    (rules / "jobs.smk").write_text(RULES)
    benchmarks = tmp_path / "benchmarks"
    for name, rows in peaks.items():
        path = benchmarks / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(HEADER + "".join(f"1.0\t0:00:01\t{rss}\t0\n" for rss in rows))
    return rules, benchmarks


def test_only_numeric_declarations_are_read(tmp_path):
    rules, _ = write_tree(tmp_path, {})
    assert declared_memory(rules) == {"small_job": 2000, "small_job_scaling": 8000}


def test_one_rule_over_and_one_under(tmp_path):
    rules, benchmarks = write_tree(tmp_path, {
        "small_job_brain_1.tsv": [150.0, 10500.0, 300.0],
        "small_job_scaling/cores4.tsv": [7000.0],
    })
    declared = declared_memory(rules)
    rows = excesses(declared, measured_memory(benchmarks, declared))
    assert [(rule, declared_mb, measured_mb) for rule, declared_mb, measured_mb, _ in rows] == [
        ("small_job", 2000, 10500.0)
    ]
    assert rows[0][3].endswith("small_job_brain_1.tsv")


def test_a_file_goes_to_the_longest_matching_rule(tmp_path):
    rules, benchmarks = write_tree(tmp_path, {"small_job_scaling_x.tsv": [5000.0]})
    declared = declared_memory(rules)
    assert set(measured_memory(benchmarks, declared)) == {"small_job_scaling"}
