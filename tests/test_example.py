"""The bundled example: the dataset it writes, and the checks it runs on the results."""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

from primer_finder.seqio import find_genomes, iter_records

ROOT = Path(__file__).resolve().parents[1]


def load(name: str):
    """Import one of the example scripts, which are not part of the package."""
    spec = importlib.util.spec_from_file_location(name, ROOT / "example" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    folder = tmp_path_factory.mktemp("example")
    assert load("make_example").main([str(folder)]) == 0
    return folder


def test_the_groups_hold_four_genomes_each(dataset):
    assert [path.name for path in find_genomes(dataset / "inclusion")] == [
        "inclusion_1.fasta", "inclusion_2.fasta", "inclusion_3.fasta", "inclusion_4.fasta.gz"
    ]
    assert len(find_genomes(dataset / "exclusion")) == 4


def test_the_planted_variants_are_where_the_truth_says(dataset):
    make_example = load("make_example")
    truth = json.loads((dataset / "truth.json").read_text())
    inclusion = "".join(r.seq for r in iter_records(dataset / "inclusion" / "inclusion_1.fasta"))
    exclusion = "".join(r.seq for r in iter_records(dataset / "exclusion" / "exclusion_1.fasta"))
    assert len(exclusion) == make_example.LENGTH
    assert len(inclusion) == make_example.LENGTH + len(truth["insertion"])
    assert inclusion[truth["insertion_at"]:truth["insertion_at"] + len(truth["insertion"])] == truth["insertion"]
    for position in truth["mismatches"]:
        offset = len(truth["insertion"]) if position > truth["insertion_at"] else 0
        assert inclusion[position + offset] != exclusion[position]


def test_every_inclusion_genome_carries_the_same_specific_region(dataset):
    truth = json.loads((dataset / "truth.json").read_text())
    window = slice(truth["target"] - 50, truth["target"] + 150)
    regions = {
        "".join(r.seq for r in iter_records(path))[window]
        for path in find_genomes(dataset / "inclusion")
    }
    assert len(regions) == 1
    assert truth["insertion"] in regions.pop()


def test_the_dataset_is_the_same_every_time(dataset, tmp_path):
    assert load("make_example").main([str(tmp_path / "again")]) == 0
    first = (dataset / "inclusion" / "inclusion_1.fasta").read_text()
    assert (tmp_path / "again" / "inclusion" / "inclusion_1.fasta").read_text() == first


def test_the_checks_fail_when_there_is_no_result(dataset, tmp_path, capsys):
    """check_example.py must not pass when primer-finder found nothing."""
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text("")
    (results / "run_info.json").write_text(json.dumps({"counts": {"final": 0, "kmers": 0}, "reference": "x"}))
    assert load("check_example").main([str(results), str(dataset)]) == 1
    assert "FAIL" in capsys.readouterr().out


PROGRAMS = ("kmc", "kmc_tools", "skesa", "minimap2", "blastn", "makeblastdb")


@pytest.mark.tools
@pytest.mark.skipif(
    not all(shutil.which(program) for program in PROGRAMS),
    reason=f"needs {', '.join(PROGRAMS)} on PATH",
)
def test_the_whole_example_with_the_real_programs(tmp_path):
    """The real pipeline on the example dataset: `conda env create -f environment.yml` first."""
    import subprocess

    result = subprocess.run(
        ["bash", str(ROOT / "example" / "run_example.sh"), str(tmp_path / "example"), "2", "4"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "All checks passed" in result.stdout
