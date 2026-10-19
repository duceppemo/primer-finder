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
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out), "-t", "1", "-m", "4"]) == 0
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


def test_the_last_resort_pass_is_on_unless_asked_otherwise(stubs, genomes, tmp_path, caplog):
    """A flag that does not reach the pipeline would leave the strict error unreachable, so this runs a
    case the strict rule turns down and checks both answers: the region, marked, by default; the error
    with `--no-last-resort`."""
    inclusion, exclusion = genomes
    stubs(sam=[["ctg1", 0, "50=1X69="]])  # one mismatch, which the strict rule does not keep

    on = tmp_path / "on"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(on)]) == 0
    assert "last_resort" in (on / "final_kmers.fasta").read_text()

    off = tmp_path / "off"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(off), "--no-last-resort"]) == 1
    assert "No contig differs enough" in caplog.text
    assert not (off / "final_kmers.fasta").exists()


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
         "-o", str(tmp_path / "out"), "-t", "1", "-m", "4"],
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
                 "-t", "1", "-m", "4"]) == 0
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
                     "-t", "1", "-m", "4"]) == 0
        assert json.loads((out / "run_info.json").read_text())["parameters"]["min_inclusion"] == 0.5


# Every one of these could be changed in `build_parser` -- a different default, or wired to nothing at all
# -- without a single test noticing, which is how the `--min-region` help text came to describe a
# restriction that no longer held. Each is passed on the command line and read back out of the run's own
# record of what it was asked for.
FIND_OPTIONS = [
    (["--assembler", "spades"], "assembler", "spades"),
    (["--min-region", "50"], "min_region", 50),
    (["-d", "3"], "duplication", 3),
    (["-k", "31"], "kmer_size", 31),
]


@pytest.mark.parametrize(("options", "name", "expected"), FIND_OPTIONS,
                         ids=[option[0].lstrip("-") for option, _, _ in FIND_OPTIONS])
def test_a_find_option_reaches_the_run(stubs, genomes, tmp_path, options, name, expected):
    inclusion, exclusion = genomes
    out = tmp_path / "out"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out),
                 "-t", "1", "-m", "4", *options]) == 0
    assert json.loads((out / "run_info.json").read_text())["parameters"][name] == expected


def test_the_no_last_resort_option_reaches_the_run(stubs, genomes, tmp_path):
    """Recorded at the top level rather than among the parameters, because it says what the run did and
    not only what it was asked for.

    The scenario matters: with a contig the strict rule accepts, `last_resort` is false whether or not
    the option is wired to anything. So this gives the contig a single difference, which only the
    last-resort pass accepts, and checks the option decides the outcome.
    """
    inclusion, exclusion = genomes
    # The whole contig, differing at one position only: too few for the strict rule, enough for the
    # last-resort one. A partial alignment would not do, since the unaligned rest counts as differing.
    whole = "ACGT" * 25
    one_difference = [["ctg1", 1, 100, 1e-30, whole, "ACGTA" + "A" + whole[6:]]]
    stubs(exclusion_hits={"exclusion_1": one_difference, "exclusion_2": one_difference})
    rescued, refused = tmp_path / "rescued", tmp_path / "refused"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(rescued),
                 "-t", "1", "-m", "4"]) == 0
    info = json.loads((rescued / "run_info.json").read_text())
    assert info["last_resort"] is True and info["counts"]["final"] == 1
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(refused),
                 "-t", "1", "-m", "4", "--no-last-resort"]) == 0
    info = json.loads((refused / "run_info.json").read_text())
    assert info["last_resort"] is False and info["counts"]["final"] == 0


def test_the_keep_intermediate_option_reaches_the_run(stubs, genomes, tmp_path):
    """Nothing in the record says whether it was given, so the files it keeps are what pins it."""
    inclusion, exclusion = genomes
    kept, tidied = tmp_path / "kept", tmp_path / "tidied"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(kept), "-t", "1", "-m", "4",
                 "--keep-intermediate"]) == 0
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(tidied), "-t", "1",
                 "-m", "4"]) == 0
    # `clean_up` empties `1_kmers` rather than removing it, so the databases are what to look for
    assert (kept / "1_kmers" / "inclusion.kmc_pre").is_file()
    assert not (tidied / "1_kmers" / "inclusion.kmc_pre").exists()


def test_the_reference_option_reaches_the_run(stubs, genomes, tmp_path):
    """Left alone the first exclusion genome is used; `-r` names another."""
    inclusion, exclusion = genomes
    chosen = sorted(exclusion.glob("*.fasta"))[-1]
    out = tmp_path / "out"
    assert main(["-i", str(inclusion), "-e", str(exclusion), "-o", str(out), "-t", "1", "-m", "4",
                 "-r", str(chosen)]) == 0
    assert Path(json.loads((out / "run_info.json").read_text())["reference"]).name == chosen.name


DESIGN_OPTIONS = [
    (["--max-regions", "7"], "max_regions", 7),
    (["--assays-per-region", "5"], "assays_per_region", 5),
    (["--gc-clamp", "0"], "gc_clamp", 0),
    (["--min-primer-differences", "1"], "min_primer_differences", 1),
]


@pytest.mark.parametrize(("options", "name", "expected"), DESIGN_OPTIONS,
                         ids=[option[0].lstrip("-") for option, _, _ in DESIGN_OPTIONS])
def test_a_design_option_reaches_the_run(stubs, genomes, tmp_path, options, name, expected):
    results = a_finished_run(tmp_path, genomes)
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1", *options]) == 0
    assert json.loads((out / "design_info.json").read_text())["parameters"][name] == expected


def test_the_design_groups_can_be_named_on_the_command_line(stubs, genomes, tmp_path):
    """`-i` and `-e` override what the run recorded, which is what lets a design be re-run against a
    different panel."""
    inclusion, exclusion = genomes
    results = a_finished_run(tmp_path, genomes)
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "1",
                 "-i", str(exclusion), "-e", str(inclusion)]) == 0  # Deliberately the other way round
    info = json.loads((out / "design_info.json").read_text())
    assert Path(info["inclusion"]).name == exclusion.name
    assert Path(info["exclusion"]).name == inclusion.name


BAD_DESIGN_OPTIONS = [
    (["--product-size", "20-40"], "starts below 27"),   # Primer3 refuses it outright
    (["--product-size", "30-53"], "ends below 54"),     # Primer3 takes it and returns nothing
    (["--assays-per-region", "0"], "--assays-per-region must be 1 or more"),
    (["--max-regions", "-1"], "--max-regions must be 0 or more"),  # 0 means every region
    (["--gc-clamp", "-1"], "--gc-clamp must be 0 or more"),
    (["--mismatches", "-1"], "--mismatches must be 0 or more"),
    (["--min-primer-differences", "-1"], "--min-primer-differences must be 0 or more"),
    (["--monovalent", "0"], "--monovalent must be above 0"),
    (["--dntp", "-1"], "--dntp must be above 0"),
]


@pytest.mark.parametrize(("options", "says"), BAD_DESIGN_OPTIONS,
                         ids=[" ".join(option) for option, _ in BAD_DESIGN_OPTIONS])
def test_a_design_setting_primer3_would_refuse_is_refused_here(stubs, genomes, tmp_path, capsys,
                                                               options, says):
    """Each of these used to reach Primer3, which refuses them in its own words and its own exit code --
    "PRIMER_NUM_RETURN < 1", "Illegal value for primer salt or dna concentration" -- or, for a product
    range under what three oligos can span, returns nothing region after region and leaves the run
    advising a wider range."""
    results = a_finished_run(tmp_path, genomes)
    with pytest.raises(SystemExit):
        main(["design", str(results), "-o", str(tmp_path / "assays"), "-t", "1", *options])
    assert says in capsys.readouterr().err


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


def test_a_product_size_range_primer3_cannot_read_stops_the_design_command(stubs, tmp_path, capsys):
    """Checking the string is not enough on its own -- the check has to be wired to the option."""
    with pytest.raises(SystemExit) as exit_code:
        main(["design", str(tmp_path / "anything"), "-o", str(tmp_path / "out"),
              "--product-size", "150-70"])
    assert exit_code.value.code == 2
    assert "runs backwards" in capsys.readouterr().err


@pytest.mark.parametrize("asked", ["0", "-1"])
def test_the_design_command_refuses_a_thread_count_below_one(stubs, tmp_path, capsys, asked):
    """The find command checks this; the design command runs blast and Primer3 too, so it checks it too."""
    with pytest.raises(SystemExit) as exit_code:
        main(["design", str(tmp_path / "anything"), "-o", str(tmp_path / "out"), "-t", asked])
    assert exit_code.value.code == 2
    assert "must be 1" in capsys.readouterr().err


def test_the_design_command_caps_the_thread_count(stubs, genomes, tmp_path, caplog):
    import json

    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text(">ctg1 [100, 101]\n" + "A" * 100 + "cg" + "A" * 100 + "\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    out = tmp_path / "assays"
    assert main(["design", str(results), "-o", str(out), "-t", "100000"]) == 0
    assert "only" in caplog.text
    assert json.loads((out / "design_info.json").read_text())["parameters"]["threads"] == usable_cpus()


def test_the_design_command_reports_a_folder_that_is_not_a_run(stubs, tmp_path, caplog):
    assert main(["design", str(tmp_path / "nothing"), "-o", str(tmp_path / "out")]) == 1
    assert "finished primer-finder run" in caplog.text


def test_design_is_a_command_not_the_default(capsys):
    """`primer-finder design ...` must not be read as the find command with a stray argument."""
    from primer_finder.cli import with_default_command

    assert with_default_command(["design", "results/", "-o", "assays/"])[0] == "design"
