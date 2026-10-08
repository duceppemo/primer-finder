"""Handing the designed assays to insilicoPCR, and reading its verdict."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.design import Assay, Oligo
from primer_finder.insilico import (
    REPORT,
    Verdict,
    command,
    count_genomes,
    read_report,
    read_terminal_only,
    resolve,
    run,
    verdicts,
    write_primer_fasta,
    write_script,
)


def oligo(name: str, sequence: str) -> Oligo:
    return Oligo(name, sequence, start=0, length=len(sequence), tm=60.0, gc=50.0)


def assay(name: str, probe: bool = True) -> Assay:
    one = Assay(region=name, number=0, forward=oligo(f"{name}-F", "AAAACCCC"),
                reverse=oligo(f"{name}-R", "TTTTGGGG"),
                probe=oligo(f"{name}-P", "CCCCAAAA") if probe else None,
                product_size=100, penalty=0.5)
    return one


def write_report(folder: Path, rows: list[tuple[str, str]]) -> Path:
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("Sample\tGene\tGenomeLocation\tAmpliconSize\n")
        for sample, gene in rows:
            fh.write(f"{sample}\t{gene}\t1-100\t100\n")
    return report


def test_the_two_kinds_of_primer_file(tmp_path):
    """insilicoPCR runs one mode at a time: a probe in the file puts the whole report in qPCR mode."""
    assays = [assay("withprobe"), assay("noprobe", probe=False)]
    qpcr, pcr = tmp_path / "q.fasta", tmp_path / "p.fasta"
    assert write_primer_fasta(assays, qpcr, "qpcr") == 1
    assert write_primer_fasta(assays, pcr, "pcr") == 2

    names = [line[1:] for line in qpcr.read_text().splitlines() if line.startswith(">")]
    assert names == ["withprobe_assay0-F", "withprobe_assay0-R", "withprobe_assay0-P"]
    names = [line[1:] for line in pcr.read_text().splitlines() if line.startswith(">")]
    assert names == ["withprobe_assay0-F", "withprobe_assay0-R", "noprobe_assay0-F", "noprobe_assay0-R"]
    assert "-P" not in pcr.read_text()


def test_nothing_to_write(tmp_path):
    path = tmp_path / "q.fasta"
    assert write_primer_fasta([assay("a", probe=False)], path, "qpcr") == 0
    assert not path.exists()


def test_an_unknown_kind(tmp_path):
    with pytest.raises(PrimerFinderError, match="Unknown assay type"):
        write_primer_fasta([assay("a")], tmp_path / "x.fasta", "lamp")


def test_resolve_a_portable_release(tmp_path):
    """The folder of an extracted release brings its own Java runtime."""
    (tmp_path / "insilicoPCR.jar").touch()
    java = tmp_path / "runtime" / "linux" / "jdk" / "bin" / "java"
    java.parent.mkdir(parents=True)
    java.touch()
    assert resolve(tmp_path) == [str(java), "-jar", str(tmp_path / "insilicoPCR.jar")]


def test_resolve_a_folder_without_the_jar(tmp_path):
    with pytest.raises(PrimerFinderError, match="No insilicoPCR.jar"):
        resolve(tmp_path)


def test_resolve_a_launcher_script(tmp_path):
    script = tmp_path / "run.sh"
    script.write_text("#!/bin/sh\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    assert resolve(script) == [str(script)]


def test_resolve_a_script_that_cannot_be_run(tmp_path):
    script = tmp_path / "run.sh"
    script.write_text("#!/bin/sh\n")
    script.chmod(0o644)
    with pytest.raises(PrimerFinderError, match="Not executable"):
        resolve(script)


def test_resolve_something_that_is_not_there(tmp_path):
    with pytest.raises(PrimerFinderError, match="No such file"):
        resolve(tmp_path / "gone.jar")


def test_resolve_a_jar_without_java(tmp_path, monkeypatch):
    jar = tmp_path / "insilicoPCR.jar"
    jar.touch()
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(PrimerFinderError, match="needs Java"):
        resolve(jar)


def test_the_command_line(tmp_path):
    assert command(["java", "-jar", "x.jar"], Path("genomes"), Path("p.fasta"), Path("out"), 8, 1) == [
        "java", "-jar", "x.jar", "-i", "genomes", "-p", "p.fasta", "-o", "out", "-t", "8", "-m", "1"
    ]


def test_read_report(tmp_path):
    write_report(tmp_path, [("sample1", "assayA"), ("sample2", "assayA"), ("sample1", "assayB")])
    assert read_report(tmp_path) == {"assayA": {"sample1", "sample2"}, "assayB": {"sample1"}}


def test_read_a_report_with_no_positive(tmp_path):
    write_report(tmp_path, [])
    assert read_report(tmp_path) == {}


def test_run_insists_on_a_report(stubs, tmp_path, insilico_pcr):
    """A run that writes no consolidated report is a failed run, however it exited."""
    launcher = resolve(insilico_pcr)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    assert run(launcher, genomes, primers, tmp_path / "out", threads=1) == tmp_path / "out"
    assert (tmp_path / "out" / REPORT).is_file()


def test_a_failing_run_is_reported(stubs, tmp_path, insilico_pcr):
    stubs(insilico_fail=True)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    with pytest.raises(PrimerFinderError, match="failed"):
        run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", threads=1)


def write_hits(folder: Path, rows: list[tuple[str, str, int, int]]) -> Path:
    """A report with the mismatch columns: (sample, gene, counted mismatches, bases trimmed off a 3' end)."""
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("Sample\tGene\tForwardMismatches\tReverseMismatches"
                 "\tForwardEndMismatch\tReverseEndMismatch\n")
        for sample, gene, mismatches, trimmed in rows:
            fh.write(f"{sample}\t{gene}\t{mismatches}\t0\t{-trimmed if trimmed else 0}\t0\n")
    return report


def test_an_amplicon_that_needed_a_trimmed_primer_end_is_recognised(tmp_path):
    """insilicoPCR does not count the last two bases of a primer: blast trims an unmatched one off and the
    primer is called bound. An amplification that needed that is one it could not have refused."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1),      # only there because a base was trimmed
                          ("g1", "a_assay0", 0, 0),      # a clean match: a real amplification
                          ("g2", "b_assay0", 1, 2)])     # trimmed as well, and a counted mismatch too
    found = read_terminal_only(tmp_path)
    assert found == {"a_assay0": {"g0"}, "b_assay0": {"g2"}}


def test_a_genome_with_one_clean_amplicon_does_not_count_as_terminal(tmp_path):
    """One assay can place several amplicons in the same genome. It only amplifies it through an ignored
    difference if every one of them needed the trim."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1), ("g0", "a_assay0", 0, 0)])
    assert read_terminal_only(tmp_path) == {}


def test_a_verdict_whose_only_exclusion_hits_are_terminal_is_the_models_blind_spot():
    """Not selective, but every exclusion genome it amplified rests on bases insilicoPCR cannot judge."""
    blind = Verdict(8, 8, 3, 17, exclusion_terminal_only=3)
    assert not blind.selective and blind.only_terminal
    partly = Verdict(8, 8, 3, 17, exclusion_terminal_only=2)
    assert not partly.only_terminal           # one of the three is a real amplification
    assert not Verdict(8, 8, 0, 17).only_terminal          # selective: nothing to explain
    assert not Verdict(7, 8, 1, 17, exclusion_terminal_only=1).only_terminal  # misses an inclusion genome


def test_verdicts_carry_the_terminal_only_count(tmp_path):
    inclusion, exclusion = tmp_path / "inclusion", tmp_path / "exclusion"
    for folder, how_many in ((inclusion, 3), (exclusion, 2)):
        folder.mkdir()
        for number in range(how_many):
            (folder / f"g{number}.fasta").write_text(">chr\nACGT\n")
    found = verdicts([assay("a")], inclusion, exclusion,
                     {"a_assay0": {"g0", "g1", "g2"}}, {"a_assay0": {"g0", "g1"}},
                     exclusion_terminal={"a_assay0": {"g0", "g1"}})
    assert found["a_assay0"].exclusion_terminal_only == 2
    assert found["a_assay0"].only_terminal


def test_verdicts_put_the_two_reports_together(tmp_path):
    inclusion, exclusion = tmp_path / "inclusion", tmp_path / "exclusion"
    for folder, how_many in ((inclusion, 3), (exclusion, 2)):
        folder.mkdir()
        for number in range(how_many):
            (folder / f"g{number}.fasta").write_text(">chr\nACGT\n")
    assays = [assay("a"), assay("b")]
    found = verdicts(
        assays, inclusion, exclusion,
        {"a_assay0": {"g0", "g1", "g2"}, "b_assay0": {"g0"}},
        {"b_assay0": {"g1"}},
    )
    assert found["a_assay0"].selective
    assert str(found["a_assay0"]) == "inclusion 3/3 (100%), exclusion 0/2"
    assert not found["b_assay0"].selective  # misses an inclusion genome and hits an exclusion one
    assert found["b_assay0"].inclusion_amplified == 1
    assert found["b_assay0"].exclusion_amplified == 1
    assert count_genomes(inclusion) == 3


def test_a_verdict_is_selective_only_when_everything_lines_up():
    assert Verdict(8, 8, 0, 17).selective
    assert not Verdict(7, 8, 0, 17).selective
    assert not Verdict(8, 8, 1, 17).selective


def test_the_script_runs_both_groups_and_both_files(tmp_path):
    path = write_script(tmp_path / "run.sh", ["java", "-jar", "x.jar"],
                        [tmp_path / "q.fasta", tmp_path / "p.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=4, mismatches=1)
    text = path.read_text()
    assert text.count("java -jar x.jar") == 4  # two files, two groups
    assert "-t 4 -m 1" in text
    assert os.access(path, os.X_OK)


def test_the_script_without_a_launcher_asks_for_one(tmp_path):
    path = write_script(tmp_path / "run.sh", None, [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=0)
    assert "INSILICO_PCR" in path.read_text()


def test_a_threshold_below_one_accepts_an_assay_that_misses_a_genome():
    """A region found with -p 0.75 may be absent from some inclusion genomes, so an assay on it is judged
    against the same share rather than against all of them."""
    almost = Verdict(inclusion_amplified=7, inclusion_total=8, exclusion_amplified=0, exclusion_total=17,
                     threshold=0.75)
    assert almost.selective and not almost.complete
    assert almost.label == "partial (87%)"
    assert almost.percent == 87

    strict = Verdict(7, 8, 0, 17)  # the default threshold is every genome
    assert not strict.selective and strict.label == "no"


def test_a_threshold_never_excuses_amplifying_the_exclusion_group():
    assert not Verdict(8, 8, 1, 17, threshold=0.5).selective


def test_below_the_threshold_is_not_selective():
    assert not Verdict(5, 8, 0, 17, threshold=0.75).selective


def test_the_percentage_is_rounded_down_so_it_never_overstates():
    assert Verdict(2, 3, 0, 1, threshold=0.6).percent == 66
    assert Verdict(0, 3, 0, 1, threshold=0.0001).percent == 0


def test_verdicts_pass_the_threshold_on(tmp_path):
    inclusion, exclusion = tmp_path / "inclusion", tmp_path / "exclusion"
    for folder, how_many in ((inclusion, 4), (exclusion, 1)):
        folder.mkdir()
        for number in range(how_many):
            (folder / f"g{number}.fasta").write_text(">chr\nACGT\n")
    found = verdicts([assay("a")], inclusion, exclusion, {"a_assay0": {"g0", "g1", "g2"}}, {},
                     threshold=0.75)
    assert found["a_assay0"].label == "partial (75%)"
