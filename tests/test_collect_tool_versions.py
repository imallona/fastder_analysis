"""Tests for collect_tool_versions.py.

Covers a built environment, a package whose name is a prefix of another, an
environment that was never built, and the two written files.
"""

import csv
import json
import subprocess

from collect_tool_versions import collect, environments, package_version, write_csv, write_tex


def make_env(root, digest, stem, packages):
    """An environment file, and its build the way Snakemake leaves it: a copy of
    the file as <hash>.yaml beside <hash>/."""
    declared = "dependencies:\n" + "".join(f"  - {package}\n" for package, _ in packages)
    for directory, name in (("envs", stem), ("conda", digest)):
        (root / directory).mkdir(exist_ok=True)
        (root / directory / f"{name}.yaml").write_text(declared)
    conda_dir = root / "conda"
    meta = conda_dir / digest / "conda-meta"
    meta.mkdir(parents=True)
    for package, version in packages:
        (meta / f"{package}-{version}-h0_0.json").write_text(
            json.dumps({"name": package, "version": version}))
    return conda_dir / digest


def make_checkout(path, tag):
    path.mkdir()
    def git(*args):
        subprocess.run(["git", "-C", str(path), *args], check=True, capture_output=True)
    git("init", "-q")
    git("-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "--allow-empty", "-m", "c")
    git("tag", tag)
    return path


def test_a_build_is_matched_to_its_environment_file_by_content(tmp_path):
    env_dir = make_env(tmp_path, "abc_", "star", [("star", "2.7.11b")])
    (tmp_path / "envs" / "never_built.yaml").write_text("dependencies:\n  - gffcompare\n")
    (tmp_path / "conda" / "def_.yaml").write_text("dependencies:\n  - something_else\n")
    found = environments(str(tmp_path / "conda"), str(tmp_path / "envs"))
    assert found == {"star": str(env_dir)}


def test_a_package_is_matched_by_name_not_by_prefix(tmp_path):
    env_dir = make_env(tmp_path, "abc_", "star",
                       [("star-fusion", "1.0"), ("star", "2.7.11b")])
    assert package_version(str(env_dir), "star") == "2.7.11b"
    assert package_version(str(env_dir), "samtools") is None


def test_only_what_is_installed_is_reported(tmp_path):
    make_env(tmp_path, "abc_", "star", [("star", "2.7.11b"), ("samtools", "1.24")])
    src = make_checkout(tmp_path / "fastder", "v0.1.0")
    version_file = tmp_path / "asimulator.txt"
    version_file.write_text("1.0.0\n")

    rows = collect(str(tmp_path / "conda"), str(tmp_path / "envs"), str(src), str(version_file))
    assert {r["tool"]: r["version"] for r in rows} == {
        "fastder": "v0.1.0", "ASimulatoR": "1.0.0", "STAR": "2.7.11b", "samtools": "1.24"}
    assert {r["source"] for r in rows} == {"git", "container", "conda:star"}


def test_missing_checkout_and_version_file_leave_rows_out(tmp_path):
    assert collect(str(tmp_path), str(tmp_path), str(tmp_path / "nothing"), None) == []


def test_written_files(tmp_path):
    rows = [{"tool": "fastder", "version": "v0.1.0-1-g8da02f5", "source": "git"},
            {"tool": "megadepth_baseline", "version": "1", "source": "conda:x"}]
    write_csv(rows, tmp_path / "v.csv")
    write_tex(rows, tmp_path / "v.tex")
    with open(tmp_path / "v.csv", newline="") as fh:
        assert list(csv.DictReader(fh)) == rows
    tex = (tmp_path / "v.tex").read_text()
    assert r"megadepth\_baseline & 1 \\" in tex
    assert tex.count(r"\\") == 3
