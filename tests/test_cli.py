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
