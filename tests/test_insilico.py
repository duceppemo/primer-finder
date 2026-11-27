"""Handing the designed assays to insilicoPCR, and reading its verdict."""

from __future__ import annotations

import dataclasses
import gzip
import os
import stat
import subprocess
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.design import Assay, Oligo
from primer_finder.insilico import (
    LENIENT_PENALTY,
    LONGEST_PRODUCT,
    REPORT,
    Reaction,
    Verdict,
    beside,
    check_samples,
    check_version,
    command,
    count_genomes,
    no_product,
    read_report,
    read_three_prime_held,
    resolve,
    run,
    stage_genomes,
    verdicts,
    weakest_tm,
    write_primer_fasta,
    write_script,
)

# The reaction primer-finder designs in, judged at insilicoPCR's defaults
REACTION = Reaction(anneal_temp=60.0, tm_margin=5.0, monovalent_mm=50.0, divalent_mm=3.0, dntp_mm=0.8,
                    oligo_nm=250.0)
REACTION_OPTIONS = ("--anneal-temp 60 --tm-margin 5 --na 50 --mg 3 --dntp 0.8 --oligo-conc 250")


def oligo(name: str, sequence: str) -> Oligo:
    return Oligo(name, sequence, start=0, length=len(sequence), tm=60.0, gc=50.0)


def assay(name: str, probe: bool = True) -> Assay:
    one = Assay(region=name, number=0, forward=oligo(f"{name}-F", "AAAACCCC"),
                reverse=oligo(f"{name}-R", "TTTTGGGG"),
                probe=oligo(f"{name}-P", "CCCCAAAA") if probe else None,
                product_size=100, penalty=0.5)
    return one


# What insilicoPCR 2.1.0 writes, in every report whatever the primer file holds
COLUMNS = ("Sample", "Assay", "GenomeLocation", "AmpliconSize", "Contig", "ContigDescription",
           "ForwardPrimer", "ReversePrimer", "ForwardMismatches", "ReverseMismatches",
           "Forward3PrimeMismatches", "Reverse3PrimeMismatches", "Strand", "Probe", "ProbeLocation",
           "ProbeSize", "ProbeMismatches", "Forward5PrimeUnpaired", "Reverse5PrimeUnpaired", "Pairing",
           "Source", "SizeCheck", "ForwardTm", "ReverseTm", "ProbeTm")


def report_row(sample: str, assay: str, description: str = "", windows: tuple[int, int] = (0, 0),
               pairing: str = "designed", tms: tuple[str, str, str] = ("60.0", "60.0", "")) -> str:
    """One row of a 2.1.0 report: a product of `assay` in `sample`."""
    w_f, w_r = windows
    return "\t".join([sample, assay, "1-100", "100", "chr", description, f"{assay}-F", f"{assay}-R",
                      str(w_f), str(w_r), str(w_f), str(w_r), "+", "", "", "", "", "0", "0", pairing,
                      "assembly", "target", *tms]) + "\n"


def write_report(folder: Path, rows: list[tuple[str, str]], description: str = "") -> Path:
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample, name in rows:
            fh.write(report_row(sample, name, description))
    return report


def test_the_two_kinds_of_primer_file(tmp_path):
    """One file per question: with the probe in the file a probe assay is only positive when the probe
    binds, so whether its primer pair discriminates alone needs a file with the probe stripped."""
    assays = [assay("withprobe"), assay("noprobe", probe=False)]
    qpcr, pcr = tmp_path / "q.fasta", tmp_path / "p.fasta"
    assert write_primer_fasta(assays, qpcr, "qpcr") == 1
    assert write_primer_fasta(assays, pcr, "pcr") == 2

    names = [line[1:].split()[0] for line in qpcr.read_text().splitlines() if line.startswith(">")]
    assert names == ["withprobe_assay0-F", "withprobe_assay0-R", "withprobe_assay0-P"]
    names = [line[1:].split()[0] for line in pcr.read_text().splitlines() if line.startswith(">")]
    assert names == ["withprobe_assay0-F", "withprobe_assay0-R", "noprobe_assay0-F", "noprobe_assay0-R"]
    assert "-P" not in pcr.read_text()


def test_the_forward_primer_declares_the_product_size(tmp_path):
    """insilicoPCR 0.11.0 reads `size=` after a primer's name as the product the assay is designed to
    give, and labels every product it reports against it. One declaration per assay is enough, and it
    sits in one place -- the forward primer -- in both files."""
    assays = [assay("withprobe"), assay("noprobe", probe=False)]
    for kind in ("qpcr", "pcr"):
        path = tmp_path / f"{kind}.fasta"
        write_primer_fasta(assays, path, kind)
        headers = [line[1:] for line in path.read_text().splitlines() if line.startswith(">")]
        assert "withprobe_assay0-F size=100" in headers
        # Nothing but the name on the others: a second declaration would have to agree with the first
        assert all(" " not in header for header in headers if not header.endswith("size=100"))
        assert sum(header.endswith("size=100") for header in headers) == len(headers) // (3 if kind == "qpcr"
                                                                                         else 2)


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
    """No mismatch budget -- insilicoPCR 2.0 refuses `-m` -- and every setting a site is judged by named
    explicitly: the annealing temperature and margin, the reaction the oligos were designed in, and `-l`.
    insilicoPCR's defaults are its own to move, and the verdicts should not move with them."""
    assert command(["java", "-jar", "x.jar"], Path("genomes"), Path("p.fasta"), Path("out"), 8,
                   REACTION) == [
        "java", "-jar", "x.jar", "-i", "genomes", "-p", "p.fasta", "-o", "out", "-t", "8",
        "-l", str(LONGEST_PRODUCT), *REACTION_OPTIONS.split(),
    ]


def test_the_reaction_is_the_one_asked_for(tmp_path):
    other = Reaction(anneal_temp=58.5, tm_margin=3.0, monovalent_mm=40.0, divalent_mm=2.5, dntp_mm=1.2,
                     oligo_nm=400.0)
    assert other.needed_tm == 55.5
    assert other.options() == ["--anneal-temp", "58.5", "--tm-margin", "3", "--na", "40", "--mg", "2.5",
                               "--dntp", "1.2", "--oligo-conc", "400"]


def test_the_command_line_with_a_three_prime_penalty(tmp_path):
    """The lenient rerun names the penalty and changes nothing else; the strict runs leave insilicoPCR's
    default 0 alone."""
    strict = command(["x"], Path("g"), Path("p"), Path("o"), 1, REACTION)
    lenient = command(["x"], Path("g"), Path("p"), Path("o"), 1, REACTION, three_prime_penalty=LENIENT_PENALTY)
    assert lenient == [*strict, "--3p-penalty", str(LENIENT_PENALTY)]
    assert "--3p-penalty" not in strict


def write_lenient_report(folder: Path, rows: list[tuple[str, str, int, int]]) -> Path:
    """A consolidated report as the lenient rerun writes it: (sample, assay, forward 3'-window
    imperfections, reverse 3'-window imperfections)."""
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample, name, w_f, w_r in rows:
            fh.write(report_row(sample, name, windows=(w_f, w_r)))
    return report


def test_a_genome_held_off_only_by_the_terminal_base_is_counted(tmp_path):
    """The lenient run differs from the strict one in its 3' penalty alone, the Tm rule being the same in
    both. A sample whose every amplicon carries a window imperfection is held off by that base alone --
    whatever else the primer carries, since the strict run accepted the rest of the site on the same Tm."""
    report = write_lenient_report(tmp_path, [("g0", "a", 1, 0),    # the readmitted terminal base
                                             ("g2", "a", 0, 1)])   # the reverse primer's terminal base
    assert read_three_prime_held(report) == {"a": {"g0", "g2"}}


def test_a_clean_lenient_amplicon_is_amplification_for_another_reason(tmp_path):
    """A row without a window imperfection did not need the 3'-end tolerance, so the sample is not held
    by the 3' end -- one such row beside a rescued one decides it."""
    report = write_lenient_report(tmp_path, [("g0", "a", 1, 0),
                                             ("g0", "a", 0, 0)])   # a second site, clean
    assert read_three_prime_held(report) == {}


def test_two_imperfections_in_one_window_are_not_the_3prime_end_alone(tmp_path):
    """The lenient run's penalty of 4 also lets through two mismatches in one window that weigh no more
    than a terminal base (the fourth and fifth bases in, 1 each). That is not a base a polymerase
    occasionally reads through, and on the Xylella panel under 0.13.0 it was a fifth of the held pairs.
    One on each primer is still one per window."""
    report = write_lenient_report(tmp_path, [("g0", "a", 2, 0),    # two in the forward window
                                             ("g1", "a", 1, 1)])   # one in each primer's window
    assert read_three_prime_held(report) == {"a": {"g1"}}


def test_an_unreadable_count_is_not_claimed_as_the_terminal_base(tmp_path):
    report = write_lenient_report(tmp_path, [("g0", "a", 1, 0)])
    report.write_text(report.read_text().replace("\t1\t0\t+", "\tNA\t0\t+"))
    assert read_three_prime_held(report) == {}


def test_verdicts_leave_real_cross_reactions_out_of_the_held_count(tmp_path):
    """A genome the strict run already counts as amplified is a cross-reaction on the model's own terms;
    only genomes the lenient run ADDS are held by the 3' end."""
    inclusion, exclusion = tmp_path / "inclusion", tmp_path / "exclusion"
    for folder, how_many in ((inclusion, 2), (exclusion, 3)):
        folder.mkdir()
        for number in range(how_many):
            (folder / f"g{number}.fasta").write_text(">chr\nACGT\n")
    found = verdicts([assay("a")], inclusion, exclusion,
                     {"a_assay0": {"g0", "g1"}}, {"a_assay0": {"g0"}},
                     exclusion_held={"a_assay0": {"g0", "g1"}})
    assert found["a_assay0"].exclusion_amplified == 1
    assert found["a_assay0"].exclusion_held_by_three_prime == 1  # g1 alone: g0 is a real cross-reaction
    assert found["a_assay0"].label == "no"  # the held count never moves the verdict


def test_read_report(tmp_path):
    report = write_report(tmp_path, [("sample1", "assayA"), ("sample2", "assayA"), ("sample1", "assayB")])
    assert read_report(report) == {"assayA": {"sample1", "sample2"}, "assayB": {"sample1"}}


def test_read_a_report_with_no_positive(tmp_path):
    assert read_report(write_report(tmp_path, [])) == {}


def test_a_product_of_one_primer_alone_is_amplification_but_not_the_target(tmp_path):
    """insilicoPCR 1.0+ reports a primer that binds both strands, facing itself at an inverted repeat, as
    a product of one primer (`Pairing` `one-primer`). It is a real product of the tube -- in an exclusion
    genome a cross-reaction -- but not the product the assay is designed to give, so it does not count as
    amplifying an inclusion genome."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True)
    report.write_text("\t".join(COLUMNS) + "\n" + report_row("g0", "a") + report_row("g1", "a", pairing="one-primer"))
    assert read_report(report) == {"a": {"g0", "g1"}}
    assert read_report(report, designed_only=True) == {"a": {"g0"}}


def test_a_quote_in_a_genomes_header_does_not_swallow_the_rows_after_it(tmp_path):
    """insilicoPCR writes the fasta description as it is, unquoted. Read with the default quote handling,
    a description beginning with a quote character swallows the rest of its row and every row up to the
    next one -- so amplifications vanish and a verdict moves."""
    report = write_report(tmp_path, [("g0", "a_assay0"), ("g1", "a_assay0")], description='"strain x')
    assert read_report(report) == {"a_assay0": {"g0", "g1"}}


def test_a_cell_written_for_spreadsheets_is_read_back(tmp_path):
    """insilicoPCR writes a cell a spreadsheet would take for a formula after an apostrophe (a genome
    named `-ctrl.fasta` is the sample `'-ctrl`), and quotes a value that starts and ends with a double
    quote. Read as written, the sample names would not match the genomes staged, nor the samples file."""
    report = write_report(tmp_path, [("'-ctrl", "a"), ('"""isolate 7"""', "a"), ("'-5x", "a"), ("-5", "a"),
                                     ("''-a", "a"), ("-1.5e3", "a"), ("'-inf", "a"), ("'-.5", "a"),
                                     ("'-1_0", "a")])
    # insilicoPCR's test of a number is its own (ReportColumn.formulaLike), not Python's float(): `-inf`,
    # `-.5` and `-1_0` are written after an apostrophe, `-1.5e3` is not
    assert read_report(report) == {"a": {"-ctrl", '"isolate 7"', "-5x", "-5", "'-a", "-1.5e3", "-inf", "-.5",
                                         "-1_0"}}


def test_a_report_without_the_columns_the_counts_come_from_is_an_error(tmp_path):
    """Read as an empty report, every assay would amplify nothing: a clean-looking answer out of a report
    nothing could be read from."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("sample\tassay\tForwardMismatches\n")  # lower case, as a later version might
    with pytest.raises(PrimerFinderError, match="no Sample or Assay column"):
        read_report(report)


def test_a_report_from_before_1_0_is_refused(tmp_path):
    """Before 1.0 the assay column was `Gene`, and every release before 2.0 accepted a site by mismatch
    counts. Counting such a report would mix two rules without showing anywhere."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    old = ["Gene" if column == "Assay" else column for column in COLUMNS[:22]]
    report.write_text("\t".join(old) + "\n")
    with pytest.raises(PrimerFinderError, match="older than 1.0.*Gene"):
        read_report(report)
    with pytest.raises(PrimerFinderError, match="older than 1.0"):
        read_three_prime_held(report)


def test_a_report_from_1_x_is_refused(tmp_path):
    """1.x writes the 1.0 columns, `Assay` included, but no Tm: its sites were accepted by mismatch
    counts, which the same assays answer otherwise."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\t".join(COLUMNS[:22]) + "\n" + "\t".join(report_row("g0", "a").split("\t")[:22]) + "\n")
    with pytest.raises(PrimerFinderError, match="no ForwardTm, ReverseTm, ProbeTm columns.*older than 2.0"):
        read_report(report)


def test_a_report_missing_another_column_is_refused(tmp_path):
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\t".join(column for column in COLUMNS if column != "Pairing") + "\n")
    with pytest.raises(PrimerFinderError, match="no Pairing column"):
        read_report(report)


def test_a_report_without_the_3prime_window_columns_is_refused(tmp_path):
    """The 3'-end check reads them; read as blank, every lenient row would count as clean and the check
    would say nothing is held, without a word."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\t".join(column for column in COLUMNS if column != "Forward3PrimeMismatches") + "\n")
    with pytest.raises(PrimerFinderError, match="no Forward3PrimeMismatches column"):
        read_three_prime_held(report)


def test_an_unexpected_size_is_still_an_amplification(tmp_path):
    """The size check labels a product and removes none, and neither does the reading: a product of
    another size than the one designed is still a product, and on the exclusion group still a
    cross-reaction."""
    report = write_report(tmp_path, [("g0", "a_assay0"), ("g1", "a_assay0")])
    lines = report.read_text().splitlines()
    lines[2] = lines[2].replace("\t100\t", "\t583\t").replace("\ttarget", "\tunexpected")
    report.write_text("\n".join(lines) + "\n")
    assert read_report(report, designed_only=True) == {"a_assay0": {"g0", "g1"}}


def test_the_weakest_site_an_assay_rests_on(tmp_path):
    """A product rests on its weakest oligo, a genome on its strongest designed product, and the assay on
    the genome where that is lowest: how close it comes to failing. A product of one primer alone is not
    the assay's, and a probe's Tm counts when there is one."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True)
    report.write_text("\t".join(COLUMNS) + "\n"
                      + report_row("g0", "a", tms=("61.0", "59.5", "68.0"))
                      + report_row("g1", "a", tms=("57.2", "60.0", "68.0"))   # weakest genome ...
                      + report_row("g1", "a", tms=("58.0", "60.0", "66.0"))   # ... at its best product
                      + report_row("g2", "a", tms=("45.0", "45.0", ""), pairing="one-primer")
                      + report_row("g0", "b", tms=("60.0", "60.0", "56.1")))
    assert weakest_tm(report) == {"a": 58.0, "b": 56.1}


def write_no_product(report: Path, rows: list[tuple[str, str, str, str]]) -> Path:
    """(sample, assay, oligo, best site Tm) rows of the file of the assays without a product, filled as
    2.1.0 fills them: a Tm with its site; "3p", a site that fails the 3' rule, which has no Tm; "", no
    site at all (or, with no oligo, no pair)."""
    path = beside(report, "no_product")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("Sample\tAssay\tOligo\tReason\tBestSite\tBestSiteTm\tNeededTm\n" + "".join(
        f"{sample}\t{name}\t{oligo_id}\twhy\t{'chr:1-20(+)' if tm else ''}\t{'' if tm == '3p' else tm}\t55.0\n"
        for sample, name, oligo_id, tm in rows))
    return path


def test_why_an_assay_gave_no_product(tmp_path):
    """insilicoPCR 2.0+ says, for each sample an assay gives no product in, which oligo has no accepted
    site and why: its best site melts too low (with that Tm), fails the 3' rule (no Tm), or there is no
    site at all. Summed per assay: the oligos, in how many genomes, for which reasons, and the best of the
    best sites -- with `no pair` where every oligo has a site but none face each other. A genome whose
    best site fails the 3' rule is not put under another genome's Tm."""
    report = tmp_path / REPORT
    write_no_product(report, [("g0", "a", "a-R", "52.1"), ("g1", "a", "a-R", "50.3"),
                              ("g1", "a", "a-F", "54.0"), ("g2", "a", "", ""),
                              ("g0", "b", "b-P_1", ""), ("g0", "c", "c-R", "3p"),
                              ("g0", "d", "d-R", "53.8"), ("g1", "d", "d-R", "3p"), ("g2", "d", "d-R", ""),
                              ("g0", "e", "e-F_1", "3p"), ("g0", "e", "e-F_2", "51.0")])
    assert no_product(report) == {
        "a": "forward: 1 genome, best site 54.0 °C; reverse: 2 genomes, best site 52.1 °C; no pair: 1 genome",
        "b": "probe: 1 genome, no site",
        "c": "reverse: 1 genome, 3' rule",
        "d": "reverse: 3 genomes (best site 53.8 °C in 1; 3' rule in 1; no site in 1)",
        # two variants of a degenerate primer in one genome: counted once, by the site with a Tm
        "e": "forward: 1 genome, best site 51.0 °C",
    }
    assert beside(tmp_path / "report(2).tsv", "no_product").name == "report(2)_no_product.tsv"


def write_samples(report: Path, rows: list[tuple[str, str, str]]) -> None:
    path = beside(report, "samples")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("Sample\tStatus\tInput\tFiles\tProducts\tReason\n" + "".join(
        f"{sample}\t{status}\tassembly\t{sample}.fasta\t\t{reason}\n" for sample, status, reason in rows))


def test_every_staged_genome_must_have_been_analysed(tmp_path):
    """insilicoPCR leaves out a genome it cannot read and still exits 0. Without a row such a genome
    reads as "does not amplify": on the exclusion group a selective assay, on the inclusion group a miss.
    So the check fails, naming it, rather than count verdicts on a smaller group."""
    report = tmp_path / REPORT
    write_samples(report, [("g0", "analysed", ""), ("g1", "not-analysed", "the file holds no sequence")])
    with pytest.raises(PrimerFinderError, match=r"left 1 genome\(s\) out.*g1 \(the file holds no sequence\)"):
        check_samples(report, staged=2)
    write_samples(report, [("g0", "analysed", ""), ("'-ctrl", "analysed", "")])
    assert check_samples(report, staged=2) == {"g0", "-ctrl"}
    with pytest.raises(PrimerFinderError, match="analysed 2 sample.*3 genome"):
        check_samples(report, staged=3)


def test_a_run_without_its_samples_file_is_refused(tmp_path):
    report = write_report(tmp_path, [("g0", "a")])
    with pytest.raises(PrimerFinderError, match="no list of the samples"):
        check_samples(report, staged=1)


def script_that(tmp_path, body: str) -> list[str]:
    launcher = tmp_path / "insilico.sh"
    launcher.write_text(f"#!/bin/sh\n{body}\n")
    launcher.chmod(0o755)
    return [str(launcher)]


def test_check_version_accepts_a_current_release(tmp_path):
    assert check_version(script_that(tmp_path, 'echo "insilicoPCR 2.1.0"')) == "insilicoPCR 2.1.0"
    assert check_version(script_that(tmp_path, 'echo "insilicoPCR 3.0"')) == "insilicoPCR 3.0"


def test_check_version_refuses_a_release_that_does_not_know_itself(tmp_path):
    """0.6.x has no -v at all: it prints an "Unrecognized option" usage dump and exits non-zero. That is
    the answer, and it has to arrive here -- before Primer3 and every blast -- not at the last step."""
    old = script_that(tmp_path, "echo 'Unrecognized option: -v' >&2\nexit 2")
    with pytest.raises(PrimerFinderError, match="needs 2.1.0 or newer"):
        check_version(old)
    with pytest.raises(PrimerFinderError, match="no version"):
        check_version(script_that(tmp_path, 'echo "ready"'))


@pytest.mark.parametrize("release", ["0.8.0", "0.12.0", "0.14.0", "1.0.3"])
def test_check_version_refuses_a_release_that_counts_mismatches(tmp_path, release):
    """Before 2.0 a site was accepted by a mismatch count, and a run with the options passed here stops
    at `--anneal-temp` -- at the last step of a design. The version says so first, and names the
    primer-finder that runs those releases."""
    with pytest.raises(PrimerFinderError, match=f"insilicoPCR {release} is older than the 2.1.0.*mismatches.*"
                                                "primer-finder 1.13.0"):
        check_version(script_that(tmp_path, f'echo "insilicoPCR {release}"'))


def test_check_version_refuses_2_0_0(tmp_path):
    """2.0.0 judges sites by their Tm as 2.1.0 does, but it is not the release this was checked against,
    and it counts a mismatch among a primer's first five bases twice where it says why an assay gave no
    product."""
    with pytest.raises(PrimerFinderError, match="insilicoPCR 2.0.0 is older than the 2.1.0.*once"):
        check_version(script_that(tmp_path, 'echo "insilicoPCR 2.0.0"'))


def staged_genomes(tmp_path: Path, *names: str) -> tuple[Path, Path]:
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    for name in names or ("g1",):
        (genomes / f"{name}.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F size=100\nAAAA\n>a-R\nTTTT\n")
    return genomes, primers


def test_run_reads_the_report_its_output_names(stubs, tmp_path, insilico_pcr):
    """insilicoPCR says where it wrote its report, and that is the report read -- not a guess at
    report.tsv or report(2).tsv. A run that says nothing, or names a file that is not there, failed."""
    launcher = resolve(insilico_pcr)
    genomes, primers = staged_genomes(tmp_path)
    assert run(launcher, genomes, primers, tmp_path / "out", 1, REACTION) == tmp_path / "out" / REPORT
    stubs(insilico_silent=True)
    with pytest.raises(PrimerFinderError, match="did not say where it wrote"):
        run(launcher, genomes, primers, tmp_path / "silent", 1, REACTION)
    stubs(insilico_no_report=True)
    with pytest.raises(PrimerFinderError, match="did not say where it wrote"):
        run(launcher, genomes, primers, tmp_path / "quiet", 1, REACTION)


def test_a_stale_report_is_set_aside_before_the_run(stubs, tmp_path, insilico_pcr):
    """insilicoPCR never overwrites a report: beside an existing report.tsv it writes report(1).tsv. The
    folder keeps the latest run's report at report.tsv, which the rerun script and anyone looking rely
    on. The old report may be hours of insilicoPCR a failed rerun would otherwise leave nothing of, so it
    is renamed rather than deleted -- and a file a user saved beside it is not insilicoPCR's to lose, so
    only the names it writes itself are touched."""
    launcher = resolve(insilico_pcr)
    genomes, primers = staged_genomes(tmp_path)
    out = tmp_path / "out"
    write_report(out, [("g1", "stale_assay")])
    parent = (out / REPORT).parent
    (parent / "report(1).tsv").write_text((out / REPORT).read_text())
    (parent / "report(1)_samples.tsv").write_text("stale\n")
    (parent / "baseline_report.tsv").write_text("saved by hand: not a name insilicoPCR writes\n")
    (parent / "report(annotated).xlsx").write_text("saved by hand, and numbered otherwise\n")
    report = run(launcher, genomes, primers, out, 1, REACTION)
    assert report == out / REPORT
    assert read_report(report) == {"a": {"g1"}}
    assert sorted(path.name for path in parent.iterdir()) == [
        "baseline_report.tsv", "report(annotated).xlsx", "report.tsv", "report.tsv.previous",
        "report_no_product.tsv", "report_samples.tsv"]
    assert not (parent / "report(1).tsv").exists() and not (parent / "report(1)_samples.tsv").exists()
    assert "stale_assay" in (parent / "report.tsv.previous").read_text()
    assert "saved by hand" in (parent / "baseline_report.tsv").read_text()


def test_a_refused_command_line_points_at_insilico_pcrs_message(stubs, tmp_path, insilico_pcr):
    """Exit 2 is insilicoPCR's "the command line cannot be used as given". The version check has passed
    by then, so it is a value it cannot use, and its own message names it."""
    genomes, primers = staged_genomes(tmp_path)
    reaction = dataclasses.replace(REACTION, tm_margin=60.0)
    with pytest.raises(PrimerFinderError, match=r"(?s)--tm-margin must be below.*refused the command line "
                                                r"\(exit 2\): its message above"):
        run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", 1, reaction)


def test_a_value_reaches_insilico_pcr_as_it_was_checked():
    """A margin just below the anneal temperature passes primer-finder's check; rounded to 60 it would be
    one insilicoPCR refuses, after the whole design."""
    options = dataclasses.replace(REACTION, tm_margin=59.9999999).options()
    assert options[options.index("--tm-margin") + 1] == "59.9999999"
    assert " ".join(REACTION.options()) == REACTION_OPTIONS


def test_a_mismatch_budget_is_never_passed(stubs, tmp_path, insilico_pcr):
    """The stub, as 2.1.0, exits 2 on `-m` or `--probe-mismatches`: a run that gets through passed none."""
    genomes, primers = staged_genomes(tmp_path)
    report = run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", 1, REACTION,
                 three_prime_penalty=LENIENT_PENALTY)
    assert report.is_file()
    with pytest.raises(PrimerFinderError, match=r"(?s)-m was removed in 2.0.*exit 2"):
        run([*resolve(insilico_pcr), "-m", "1"], genomes, primers, tmp_path / "budget", 1, REACTION)


def test_a_failing_run_is_reported(stubs, tmp_path, insilico_pcr):
    stubs(insilico_fail=True)
    genomes, primers = staged_genomes(tmp_path)
    with pytest.raises(PrimerFinderError, match="failed"):
        run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", 1, REACTION)


def test_a_genome_that_cannot_be_read_is_reported_before_the_design_is_lost(tmp_path):
    """Staging happens after Primer3 and every blast. A gzip that stops early passes the fasta check,
    which reads only the first line, and used to raise a bare EOFError from there."""
    folder, into = tmp_path / "genomes", tmp_path / "staged"
    folder.mkdir()
    whole = gzip.compress(b">chr\nACGTACGTACGT\n")
    (folder / "g1.fasta.gz").write_bytes(whole[: len(whole) // 2])
    with pytest.raises(PrimerFinderError, match="Cannot read .*g1.fasta.gz"):
        stage_genomes(folder, into)


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


def test_a_verdict_on_no_inclusion_genome_at_all_is_not_a_pass():
    """0 of 0 is not every one of them. `fraction` guards the division; `complete` has to guard the
    comparison, or an empty group would read as `yes`."""
    nothing = Verdict(0, 0, 0, 17)
    assert not nothing.complete and not nothing.selective
    assert nothing.label == "no" and nothing.support == 0


def test_the_labels_of_a_verdict():
    assert Verdict(8, 8, 0, 17).label == "yes"
    assert Verdict(7, 8, 0, 17, threshold=0.75).label == "partial (87%)"
    assert Verdict(8, 8, 3, 17).label == "no"


def test_what_the_ranking_takes_from_a_verdict():
    """Amplifying the whole inclusion group and nothing else outranks merely meeting the run's threshold,
    which outranks a missed group or any exclusion amplification."""
    assert Verdict(8, 8, 0, 17).support == 2
    assert Verdict(7, 8, 0, 17, threshold=0.75).support == 1
    assert Verdict(8, 8, 1, 17).support == 0
    assert Verdict(4, 8, 0, 17).support == 0  # does not cover the inclusion group at all



def test_verdicts_carry_the_weakest_site_and_why_a_genome_was_missed(tmp_path):
    inclusion, exclusion = tmp_path / "inclusion", tmp_path / "exclusion"
    for folder, how_many in ((inclusion, 2), (exclusion, 1)):
        folder.mkdir()
        for number in range(how_many):
            (folder / f"g{number}.fasta").write_text(">chr\nACGT\n")
    found = verdicts([assay("a"), assay("b")], inclusion, exclusion,
                     {"a_assay0": {"g0", "g1"}, "b_assay0": {"g0"}}, {},
                     inclusion_tm={"a_assay0": 57.25, "b_assay0": 59.0},
                     inclusion_missing={"b_assay0": "reverse: 1 genome, best site 52.1 °C"})
    assert found["a_assay0"].inclusion_weakest_tm == 57.25 and found["a_assay0"].inclusion_no_product == ""
    assert found["b_assay0"].inclusion_no_product == "reverse: 1 genome, best site 52.1 °C"


def test_the_script_runs_both_groups_and_both_files(tmp_path):
    path = write_script(tmp_path / "run.sh", ["java", "-jar", "x.jar"],
                        [tmp_path / "q.fasta", tmp_path / "p.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", 4, REACTION)
    text = path.read_text()
    # two files, two groups and the lenient exclusion rerun -- after asking the release its version once
    assert text.count("java -jar x.jar") == 7 and text.count("java -jar x.jar -v") == 1
    assert text.count(f"-t 4 -l {LONGEST_PRODUCT} {REACTION_OPTIONS}\n") == 4
    # The lenient rerun: the far 3' end readmitted, and nothing else changed
    assert text.count(f"-t 4 -l {LONGEST_PRODUCT} {REACTION_OPTIONS} --3p-penalty {LENIENT_PENALTY}") == 2
    assert " -m " not in text and "--probe-mismatches" not in text
    assert "reaches 55 C" in text
    assert text.count("_exclusion_3prime") == 5  # the header's word, then two files x (dir= and -o)
    # insilicoPCR writes report(1).tsv beside an existing report, which the find below never names -- the
    # old report is set aside, not deleted, and nothing a user saved beside it is touched
    assert text.count('mv -f "${dir}/report.tsv" "${dir}/report.tsv.previous"') == 6
    assert text.count('if [[ "${stale##*/}" =~ ${numbered} ]]; then rm -f "${stale}"; fi') == 6
    assert os.access(path, os.X_OK)


def test_the_script_runs_from_a_launcher_path_with_a_space(stubs, tmp_path, insilico_pcr):
    """The launcher is quoted, and a rerun removes the numbered reports insilicoPCR wrote and nothing a
    user saved beside them."""
    folder = tmp_path / "my tools"
    folder.mkdir()
    (folder / "insilicoPCR").symlink_to(insilico_pcr)
    genomes, primers = staged_genomes(tmp_path)
    path = write_script(tmp_path / "run.sh", [str(folder / "insilicoPCR")], [primers], genomes, genomes, 1,
                        REACTION, three_prime_check=False)
    out = tmp_path / "out"
    kept = out / "p_inclusion" / "consolidated_report"
    kept.mkdir(parents=True)
    (kept / "report(3).tsv").write_text("stale\n")
    (kept / "report(3)_samples.tsv").write_text("stale\n")
    (kept / "report(annotated).xlsx").write_text("saved by hand\n")
    ran = subprocess.run(["bash", str(path), str(out)], capture_output=True, text=True)
    assert ran.returncode == 0, ran.stderr
    assert sorted(child.name for child in kept.iterdir()) == [
        "report(annotated).xlsx", "report.tsv", "report_no_product.tsv", "report_samples.tsv"]


def test_the_script_without_the_three_prime_check_runs_two_groups_only(tmp_path):
    path = write_script(tmp_path / "run.sh", ["x"], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", 2, REACTION,
                        three_prime_check=False)
    text = path.read_text()
    assert "--3p-penalty" not in text and "_exclusion_3prime" not in text


def test_the_script_without_a_launcher_asks_for_one(tmp_path):
    path = write_script(tmp_path / "run.sh", None, [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", 2, REACTION)
    assert "INSILICO_PCR" in path.read_text()
    ran = subprocess.run(["bash", str(path), str(tmp_path / "out")], capture_output=True, text=True,
                         env={key: value for key, value in os.environ.items() if key != "INSILICO_PCR"})
    assert ran.returncode != 0 and "point INSILICO_PCR at" in ran.stderr


@pytest.mark.parametrize("says, runs", [("insilicoPCR 2.1.0", True), ("insilicoPCR 3.0", True),
                                        ("insilicoPCR 2.0.0", False), ("insilicoPCR 1.0.3", False),
                                        ("insilicoPCR 0.12.0", False), ("", False)])
def test_the_script_refuses_a_release_older_than_the_one_it_was_written_for(tmp_path, says, runs):
    """A design checks the release before it starts; the script is run by hand, later, against whatever
    INSILICO_PCR or the recorded launcher is by then. A release older than 2.0 stops on the options
    passed; 2.0.0 is not the release the settings were checked against. So the script asks first -- and
    runs nothing when the answer is too old or no answer."""
    launcher = tmp_path / "insilico.sh"
    launcher.write_text("#!/bin/sh\n"
                        f'if [ "$1" = "-v" ]; then echo "{says}"; exit 0; fi\n'
                        f'echo ran >> "{tmp_path}/calls"\n')
    launcher.chmod(0o755)
    path = write_script(tmp_path / "run.sh", [str(launcher)], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", 2, REACTION)
    (tmp_path / "out").mkdir()  # The real program makes it; the script ends by listing what is in it
    ran = subprocess.run(["bash", str(path), str(tmp_path / "out")], capture_output=True, text=True)
    if runs:
        assert ran.returncode == 0, ran.stderr
        assert (tmp_path / "calls").read_text().count("ran") == 3
    else:
        assert ran.returncode == 1 and "2.1.0 or newer is needed" in ran.stderr
        assert not (tmp_path / "calls").exists()


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


def test_the_percentage_is_exact_for_a_fraction_a_float_cannot_hold():
    """`int(29 / 100 * 100)` is 28, so an assay that amplified 29 of 100 genomes used to read
    `partial (28%)` -- under a `-p 0.29` the user had asked for. Twelve such pairs occur under 300
    genomes; these are five of them."""
    for amplified, total in ((29, 100), (57, 100), (58, 100), (29, 50), (87, 150)):
        exact = amplified * 100 // total
        assert Verdict(amplified, total, 0, 1).percent == exact
        assert int(amplified / total * 100) < exact  # what it used to give


def test_staging_forgets_a_genome_that_has_left_the_input_folder(tmp_path):
    """The staged folder is what the verdicts are counted against, so a genome left behind by an earlier
    run would be counted as one more inclusion genome and make a selective assay read `partial`."""
    source = tmp_path / "inclusion"
    source.mkdir()
    for name in ("a", "b"):
        (source / f"{name}.fasta").write_text(f">chr\n{name.upper() * 4}\n")
    staged = tmp_path / "staged"
    assert stage_genomes(source, staged) == 2
    (source / "b.fasta").unlink()
    assert stage_genomes(source, staged) == 1
    assert sorted(path.name for path in staged.glob("*.fasta")) == ["a.fasta"]
    assert count_genomes(staged) == 1  # the number the verdicts divide by


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


def test_staging_flattens_subfolders_and_decompresses(tmp_path):
    """insilicoPCR reads only the top level of the folder it is given, so a genome in a subfolder would be
    skipped. Gzipped genomes are written out plain too, so a staged sample is one plain fasta whatever it
    came from."""
    import gzip

    source = tmp_path / "inclusion"
    (source / "batch1").mkdir(parents=True)
    (source / "plain.fasta").write_text(">chr\nACGT\n")
    (source / "batch1" / "nested.fasta").write_text(">chr\nTTTT\n")
    with gzip.open(source / "packed.fasta.gz", "wt") as fh:
        fh.write(">chr\nGGGG\n")

    staged = tmp_path / "staged"
    assert stage_genomes(source, staged) == 3
    assert sorted(path.name for path in staged.iterdir()) == ["nested.fasta", "packed.fasta", "plain.fasta"]
    assert (staged / "packed.fasta").read_text() == ">chr\nGGGG\n"  # written out plain, not a link
    assert (staged / "plain.fasta").is_symlink()


def test_two_genomes_with_the_same_name_stay_two_samples(tmp_path):
    """Flattening brings together names that were unique only within their own subfolder. Letting one
    overwrite the other would quietly shrink the panel the verdicts are counted against."""
    source = tmp_path / "exclusion"
    for sample in ("a", "b"):
        (source / sample).mkdir(parents=True)
        (source / sample / "contigs.fasta").write_text(f">chr\n{sample.upper() * 4}\n")

    staged = tmp_path / "staged"
    assert stage_genomes(source, staged) == 2
    assert sorted(path.name for path in staged.iterdir()) == ["contigs.fasta", "contigs_2.fasta"]


def test_staging_can_be_run_again_over_the_same_folder(tmp_path):
    source = tmp_path / "inclusion"
    source.mkdir()
    (source / "g.fasta").write_text(">chr\nACGT\n")
    staged = tmp_path / "staged"
    stage_genomes(source, staged)
    assert stage_genomes(source, staged) == 1


def test_staging_into_a_folder_that_holds_the_genomes_is_refused(tmp_path):
    """`stage_genomes` empties what it writes to, so it must refuse to write into the folder it reads
    from -- or into a parent of it. The design command refuses the arrangement too; this is the layer
    that protects any other caller."""
    folder = tmp_path / "genomes"
    folder.mkdir()
    (folder / "g.fasta").write_text(">chr\nACGT\n")
    for into in (folder, folder / "staged", tmp_path):
        with pytest.raises(PrimerFinderError, match="one is inside the other"):
            stage_genomes(folder, into)
    assert (folder / "g.fasta").is_file()
    # A folder that is neither is fine
    assert stage_genomes(folder, tmp_path / "elsewhere" / "staged") == 1
