"""The command line."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

from primer_finder import __version__
from primer_finder.cli import main, with_default_command
from primer_finder.system import default_memory_gb, usable_cpus

ROOT = Path(__file__).resolve().parents[1]


def test_version_is_the_same_everywhere():
    """pyproject.toml, the package, CITATION.cff and the conda recipe must agree."""
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert re.search(rf'^version = "{re.escape(__version__)}"$', pyproject, re.M)
    citation = (ROOT / "CITATION.cff").read_text()
    assert f"\nversion: {__version__}\n" in citation
    recipe = (ROOT / "recipe" / "meta.yaml").read_text()
    assert f'set version = "{__version__}"' in recipe
    assert f"## {__version__}" in (ROOT / "CHANGELOG.md").read_text()


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exit_code:
        main(["--version"])
    assert exit_code.value.code == 0
    assert capsys.readouterr().out.strip() == f"primer-finder {__version__}"


def test_no_argument_prints_the_help(capsys):
    assert main([]) == 1
    assert "usage: primer-finder" in capsys.readouterr().out


def test_help_of_the_find_command(capsys):
    with pytest.raises(SystemExit):
        main(["find", "--help"])
    out = capsys.readouterr().out
    assert "--inclusion" in out and "--kmer_size" in out


@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (["-i", "a", "-e", "b", "-o", "c"], ["find", "-i", "a", "-e", "b", "-o", "c"]),
        (["find", "-i", "a"], ["find", "-i", "a"]),
        (["idt", "x", "y"], ["idt", "x", "y"]),
        (["--version"], ["--version"]),
        (["-h"], ["-h"]),
        (["nonsense"], ["nonsense"]),
        ([], []),
    ],
)
def test_the_find_command_is_the_default(argv, expected):
    assert with_default_command(argv) == expected


def test_an_unknown_command_is_refused(capsys):
    with pytest.raises(SystemExit):
        main(["nonsense"])
    assert "invalid choice" in capsys.readouterr().err


def test_a_full_run_from_the_command_line(stubs, genomes, tmp_path):
    inclusion, exclusion = genomes
    out = tmp_path / "out"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out), "-t", "1", "-m", "2"]) == 0
    info = json.loads((out / "run_info.json").read_text())
    assert info["command_line"][0] == "primer-finder"
    assert "--version" not in info["command_line"]
    assert info["parameters"]["kmer_size"] == 99


def test_errors_are_reported_without_a_traceback(stubs, genomes, tmp_path, caplog):
    inclusion, exclusion = genomes
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(tmp_path / "out"), "-k", "0"]) == 1
    assert "between 1 and 256" in caplog.text


def test_too_many_threads_and_too_much_memory_are_capped(stubs, genomes, tmp_path, caplog):
    inclusion, exclusion = genomes
    out = tmp_path / "out"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out),
                 "-t", "100000", "-m", "100000"]) == 0
    assert "only" in caplog.text
    info = json.loads((out / "run_info.json").read_text())
    assert info["parameters"]["threads"] == usable_cpus()
    assert info["parameters"]["memory_gb"] == default_memory_gb()


@pytest.mark.parametrize("argv", [["-t", "0"], ["-m", "0"]])
def test_threads_and_memory_must_be_positive(stubs, genomes, tmp_path, argv, capsys):
    inclusion, exclusion = genomes
    with pytest.raises(SystemExit):
        main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(tmp_path / "out"), *argv])
    assert "must be 1" in capsys.readouterr().err


def test_an_unknown_assembler_is_refused(stubs, genomes, tmp_path, capsys):
    inclusion, exclusion = genomes
    with pytest.raises(SystemExit):
        main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(tmp_path / "out"), "-a", "velvet"])
    assert "invalid choice" in capsys.readouterr().err


def test_the_idt_command(tmp_path):
    table = tmp_path / "order.tsv"
    table.write_text("Type\tSequence\tAmplicon\nForward Primer\tAAAA\t\nReverse Primer\tTTTT\t95\n")
    output = tmp_path / "assays.fasta"
    assert main(["idt", str(table), str(output), "pf"]) == 0
    assert output.read_text() == ">pf_0_95bp-F\nAAAA\n>pf_0_95bp-R\nTTTT\n"


def test_the_idt_command_reports_a_bad_file(tmp_path, caplog):
    assert main(["idt", str(tmp_path / "missing.xlsx"), str(tmp_path / "out.fasta")]) == 1
    assert "No such file" in caplog.text


def test_the_old_scripts_still_work(stubs, genomes, tmp_path):
    """primer_finder.py and IDT_results_converter.py are kept as entry points."""
    import subprocess

    inclusion, exclusion = genomes
    result = subprocess.run(
        [sys.executable, str(ROOT / "primer_finder.py"), "-i", str(inclusion), "-e", str(exclusion),
         "-o", str(tmp_path / "out"), "-t", "1", "-m", "2"],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out" / "final_kmers.fasta").exists()

    table = tmp_path / "order.tsv"
    table.write_text("Type\tSequence\tAmplicon\nForward Primer\tAAAA\t95\nReverse Primer\tTTTT\t\n")
    result = subprocess.run(
        [sys.executable, str(ROOT / "IDT_results_converter.py"), str(table), str(tmp_path / "a.fasta")],
        capture_output=True, text=True, cwd=ROOT,
    )
    assert result.returncode == 0, result.stderr
    assert "Deprecated" in result.stderr
    assert (tmp_path / "a.fasta").read_text().startswith(">0_95bp-F")


@pytest.mark.parametrize("given,expected", [("0.75", 0.75), ("1", 1.0), ("1.0", 1.0), (".9", 0.9),
                                            ("1e-3", 0.001)])
def test_the_threshold_is_read_as_a_fraction(stubs, genomes, tmp_path, given, expected):
    inclusion, exclusion = genomes
    out = tmp_path / f"out{given}"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out), "-p", given,
                 "-t", "1", "-m", "2"]) == 0
    parameters = json.loads((out / "run_info.json").read_text())["parameters"]
    assert parameters["min_inclusion"] == expected


@pytest.mark.parametrize("given", ["0", "-0.5", "1.5", "2", "abc", "", "nan", "inf"])
def test_a_threshold_the_option_cannot_take(stubs, genomes, tmp_path, given, capsys):
    inclusion, exclusion = genomes
    with pytest.raises(SystemExit):
        main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(tmp_path / "out"), "-p", given])
    error = capsys.readouterr().err
    assert "--min-inclusion" in error or "min_inclusion" in error


def test_the_long_and_short_threshold_options_agree(stubs, genomes, tmp_path):
    inclusion, exclusion = genomes
    for number, option in enumerate(("-p", "--min-inclusion")):
        out = tmp_path / f"out{number}"
        assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out), option, "0.5",
                     "-t", "1", "-m", "2"]) == 0
        assert json.loads((out / "run_info.json").read_text())["parameters"]["min_inclusion"] == 0.5


def test_the_design_command(stubs, genomes, tmp_path):
    """`primer-finder design` reads a finished run and writes the assays."""
    import json

    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text(">ctg1 [100, 101]\n" + "A" * 100 + "cg" + "A" * 100 + "\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1"]) == 0
    assert (out / "assays.tsv").is_file()
    info = json.loads((out / "design_info.json").read_text())
    assert info["command_line"][:2] == ["primer-finder", "design"]
    assert info["parameters"]["product_size"] == "70-150"


def iter_rows(path: Path):
    import csv

    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def a_finished_run(tmp_path, genomes) -> Path:
    """A results folder the design command will accept, with one region holding two differences."""
    import json

    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text(">ctg1 [100, 101]\n" + "A" * 100 + "cg" + "A" * 100 + "\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    return results


def primer3_inputs(output: Path) -> list[str]:
    """Exactly what was sent to Primer3, which the command keeps."""
    sent = sorted((output / "primer3").glob("*_input.txt"))
    assert sent, "the design command wrote no Primer3 input"
    return [path.read_text() for path in sent]


def test_the_reaction_conditions_reach_the_design_run(stubs, genomes, tmp_path):
    """They change both the melting temperatures Primer3 predicts and which oligos it will return, so what
    matters is that they arrive in the request -- not only that the run recorded them."""
    import json

    results = a_finished_run(tmp_path, genomes)
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1", "--max-hairpin-tm", "40",
                 "--max-dimer-tm", "35", "--monovalent", "60", "--divalent", "5", "--dntp", "1.2",
                 "--primer-conc", "500", "--probe-conc", "100"]) == 0
    conditions = json.loads((out / "design_info.json").read_text())["parameters"]["conditions"]
    assert conditions == {"max_hairpin_tm": 40.0, "max_dimer_tm": 35.0, "monovalent_mm": 60.0,
                          "divalent_mm": 5.0, "dntp_mm": 1.2, "primer_nm": 500.0, "probe_nm": 100.0}
    asked = {"PRIMER_MAX_HAIRPIN_TH=40.0", "PRIMER_MAX_SELF_ANY_TH=35.0", "PRIMER_SALT_MONOVALENT=60.0",
             "PRIMER_SALT_DIVALENT=5.0", "PRIMER_DNTP_CONC=1.2", "PRIMER_DNA_CONC=500.0",
             "PRIMER_INTERNAL_DNA_CONC=100.0"}
    sent = primer3_inputs(out)
    for request in sent:   # every request, the ones pinned on a difference included
        assert asked <= set(request.splitlines()), request[:80]
    assert any("SEQUENCE_FORCE_LEFT_END" in request for request in sent), "no forced request was made"


def test_the_gc_clamp_reaches_the_design_run_except_where_an_end_is_pinned(stubs, genomes, tmp_path):
    """PRIMER_GC_CLAMP is one tag for the whole request, so a request that pins a primer's 3' end on a
    differing base cannot ask for a G or a C there as well: the clamp is dropped for it."""
    results = a_finished_run(tmp_path, genomes)
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1", "--gc-clamp", "2"]) == 0
    pinned, free = [], []
    for request in primer3_inputs(out):
        lines = request.splitlines()
        clamp = next(line for line in lines if line.startswith("PRIMER_GC_CLAMP="))
        (pinned if any(line.startswith("SEQUENCE_FORCE_") for line in lines) else free).append(clamp)
    assert pinned and free
    assert set(free) == {"PRIMER_GC_CLAMP=2"}
    assert set(pinned) == {"PRIMER_GC_CLAMP=0"}


def test_the_minimum_differences_under_one_primer_reaches_the_design_run(stubs, genomes, tmp_path):
    import json

    results = a_finished_run(tmp_path, genomes)
    for asked, flag in ((1, "--min-primer-differences"), (3, "--min-oligo-differences")):
        out = tmp_path / f"assays{asked}"
        assert main(["design", str(results), "-o", str(out), "-t", "1", flag, str(asked)]) == 0
        recorded = json.loads((out / "design_info.json").read_text())["parameters"]
        assert recorded["min_primer_differences"] == asked   # the old name still works
    # These assays are specific by absence, which the filter does not apply to, so asking for three
    # differences under one oligo takes none of them away
    assert len(list(iter_rows(tmp_path / "assays1" / "assays.tsv"))) == \
        len(list(iter_rows(tmp_path / "assays3" / "assays.tsv")))


def test_another_mismatch_tolerance_reaches_insilico_pcr(stubs, genomes, tmp_path, insilico_pcr):
    import json

    calls = stubs()
    results = a_finished_run(tmp_path, genomes)
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1", "--mismatches", "2",
                 "--insilico-pcr", str(insilico_pcr)]) == 0
    assert json.loads((out / "design_info.json").read_text())["parameters"]["mismatches"] == 2
    asked = [line for line in calls.read_text().splitlines() if line.startswith("insilico_pcr")]
    assert asked and all(" -m 2" in line for line in asked), asked


def test_the_design_command_reports_a_folder_that_is_not_a_run(stubs, tmp_path, caplog):
    assert main(["design", str(tmp_path / "nothing"), "-o", str(tmp_path / "out")]) == 1
    assert "finished primer-finder run" in caplog.text


def test_design_is_a_command_not_the_default(capsys):
    """`primer-finder design ...` must not be read as the find command with a stray argument."""
    from primer_finder.cli import with_default_command

    assert with_default_command(["design", "results/", "-o", "assays/"])[0] == "design"
