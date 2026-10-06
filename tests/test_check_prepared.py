from check_prepared import (commit_problem, forbidden_problem, manifest_problem, plan_problem,
                            planned_rules, write_manifest)

PLAN = """Building DAG of jobs...
Job stats:
job              count
-------------  -------
build_fastder        1
ml_star_align        4
total                5

[Mon Oct  5 10:00:00 2026]
rule ml_star_align:
    jobid: 3
"""


def test_the_job_table_is_read():
    assert planned_rules(PLAN) == {"build_fastder": 1, "ml_star_align": 4}


def test_a_finished_preparation_has_no_problem():
    assert plan_problem("Nothing to be done (all requested files are present).", {"build_fastder"}) is None
    only_build = PLAN.replace("ml_star_align        4\n", "")
    assert plan_problem(only_build, {"build_fastder"}) is None


def test_a_rule_that_would_run_again_is_named():
    assert plan_problem(PLAN, {"build_fastder"}) == "preparation would be redone: ml_star_align (4)"


def test_a_listed_rule_left_to_run_is_named():
    assert forbidden_problem(PLAN, {"ml_star_align"}) == "left to run: ml_star_align (4)"
    assert forbidden_problem(PLAN, {"run_fastder"}) is None
    assert forbidden_problem("", {"run_fastder"}) == "the dry run printed no plan"


def test_a_failed_dry_run_is_a_problem():
    assert plan_problem("MissingInputException in rule x", {"build_fastder"}) == "the dry run printed no plan"


def test_the_prepared_commit_must_be_the_checked_out_one(tmp_path):
    commit_file = tmp_path / "prepared_commit.txt"
    commit_file.write_text("abc123\n")
    assert commit_problem(commit_file, "abc123") is None
    assert commit_problem(commit_file, "def456") == "prepared at abc123, checked out at def456"
    assert "is missing" in commit_problem(tmp_path / "absent.txt", "abc123")


def prepared_tree(tmp_path):
    root = tmp_path / "fastder"
    scenario = root / "config" / "brain_1"
    scenario.mkdir(parents=True)
    (scenario / "a.all.bw").write_bytes(b"x" * 100)
    (scenario / "b.all.bw").write_bytes(b"x" * 200)
    (root / "config" / "reference_label.gtf").write_text("chr1\n")
    (scenario / "reference_label.gtf").symlink_to("../reference_label.gtf")
    manifest = tmp_path / "prepared_manifest.tsv"
    write_manifest(root, manifest)
    return root, scenario, manifest


def test_a_complete_copy_matches_its_manifest(tmp_path):
    root, scenario, manifest = prepared_tree(tmp_path)
    (scenario / "runs").mkdir()
    (scenario / "runs" / "later_output.gtf").write_text("made after the copy\n")
    assert manifest_problem(root, manifest) is None


def test_a_missing_and_a_truncated_file_are_counted(tmp_path):
    root, scenario, manifest = prepared_tree(tmp_path)
    (scenario / "a.all.bw").unlink()
    (scenario / "b.all.bw").write_bytes(b"x" * 50)
    problem = manifest_problem(root, manifest)
    assert problem.startswith("1 files missing and 1 of another size")
    assert "config/brain_1/a.all.bw" in problem


def test_a_tree_without_manifest_is_a_problem(tmp_path):
    root, _, manifest = prepared_tree(tmp_path)
    manifest.unlink()
    assert "is missing" in manifest_problem(root, manifest)
