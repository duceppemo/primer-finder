"""Handing the designed assays to insilicoPCR, and reading its verdict."""

from __future__ import annotations

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
    Verdict,
    check_version,
    command,
    count_genomes,
    read_report,
    read_three_prime_held,
    resolve,
    run,
    stage_genomes,
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


# What insilicoPCR 0.12.0 writes for a primer file that declares product sizes, as every file written
# here does; primer-finder reads Sample and Gene and insists the rest is not an older release's report.
COLUMNS = ("Sample", "Gene", "GenomeLocation", "AmpliconSize", "Contig", "Contig Description",
           "ForwardPrimers", "ReversePrimers", "ForwardMismatches", "ReverseMismatches",
           "Forward3PrimeMismatches", "Reverse3PrimeMismatches", "Strand",
           "Forward5PrimeUnaligned", "Reverse5PrimeUnaligned", "SizeCheck")


def write_report(folder: Path, rows: list[tuple[str, str]], description: str = "") -> Path:
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample, gene in rows:
            fh.write(f"{sample}\t{gene}\t1-100\t100\tchr\t{description}"
                     f"\t{gene}-F\t{gene}-R\t0\t0\t0\t0\t+\t0\t0\ttarget\n")
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
    """`-l` is always passed: insilicoPCR 0.8.0 applies the amplicon size cap itself, and naming it keeps
    the verdicts counted the same however its default moves."""
    assert command(["java", "-jar", "x.jar"], Path("genomes"), Path("p.fasta"), Path("out"), 8, 1) == [
        "java", "-jar", "x.jar", "-i", "genomes", "-p", "p.fasta", "-o", "out", "-t", "8", "-m", "1",
        "-l", str(LONGEST_PRODUCT),
    ]


def test_the_command_line_with_a_probe_tolerance(tmp_path):
    """insilicoPCR's own default is the `-m` value, so the flag is only passed when the user chose one --
    a None must not become the string \"None\" on the command line."""
    with_probe = command(["x"], Path("g"), Path("p"), Path("o"), 1, 1, probe_mismatches=2)
    assert with_probe[-2:] == ["--probe-mismatches", "2"]
    assert "--probe-mismatches" not in command(["x"], Path("g"), Path("p"), Path("o"), 1, 1)


def test_the_command_line_with_a_three_prime_penalty(tmp_path):
    """The lenient rerun names the penalty; the strict runs leave insilicoPCR's default 0 alone."""
    lenient = command(["x"], Path("g"), Path("p"), Path("o"), 1, 2, three_prime_penalty=LENIENT_PENALTY)
    assert lenient[-2:] == ["--3p-penalty", str(LENIENT_PENALTY)]
    assert "--3p-penalty" not in command(["x"], Path("g"), Path("p"), Path("o"), 1, 2)


def write_lenient_report(folder: Path, rows: list[tuple[str, str, int, int, int, int]]) -> Path:
    """A consolidated report as the lenient rerun writes it: (sample, gene, forward mismatches, forward
    3'-window imperfections, reverse mismatches, reverse 3'-window imperfections)."""
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample, gene, m_f, w_f, m_r, w_r in rows:
            fh.write(f"{sample}\t{gene}\t1-100\t100\tchr\t\t{gene}-F\t{gene}-R"
                     f"\t{m_f}\t{m_r}\t{w_f}\t{w_r}\t+\t0\t0\ttarget\n")
    return report


def test_a_genome_held_off_only_by_the_terminal_base_is_counted(tmp_path):
    """The lenient run readmits the 3'-terminal base and the one unit of -m the trim consumes. A sample
    whose every amplicon carries a window imperfection, with everything outside the window inside the
    user's own budget, is held off by that base alone."""
    write_lenient_report(tmp_path, [("g0", "a", 1, 1, 0, 0),    # the readmitted terminal base
                                    ("g1", "a", 2, 1, 0, 0),    # window plus one internal: budget 1 holds
                                    ("g2", "a", 0, 0, 1, 1)])   # the reverse primer's terminal base
    assert read_three_prime_held(tmp_path, budget=1) == {"a": {"g0", "g1", "g2"}}


def test_a_clean_lenient_amplicon_is_amplification_for_another_reason(tmp_path):
    """A row without a window imperfection did not need the 3'-end tolerance, so the sample is not held
    by the 3' end -- one such row beside a rescued one decides it."""
    write_lenient_report(tmp_path, [("g0", "a", 1, 1, 0, 0),
                                    ("g0", "a", 0, 0, 0, 0)])   # a second site, clean
    assert read_three_prime_held(tmp_path, budget=1) == {}


def test_an_internal_mismatch_over_the_budget_is_not_the_terminal_base(tmp_path):
    """The lenient run's extra -m unit exists to pay for the trimmed terminal base. A row that spends it
    on a second internal mismatch instead would not amplify at the user's own budget however the 3' end
    read, so it must not be attributed to the 3' end."""
    write_lenient_report(tmp_path, [("g0", "a", 3, 1, 0, 0)])   # window 1, internals 2, budget 1
    assert read_three_prime_held(tmp_path, budget=1) == {}
    assert read_three_prime_held(tmp_path, budget=2) == {"a": {"g0"}}


def test_two_imperfections_in_one_window_are_not_the_3prime_end_alone(tmp_path):
    """insilicoPCR 0.12.0 counts a mismatch near a primer's 3' end once, so the lenient run's penalty of
    4 also lets through two mismatches in one window that weigh no more than a terminal base (the fourth
    and fifth bases in, 1 each). That is not a base a polymerase occasionally reads through, and 0.11.0
    never reported it; on the Xylella panel under 0.13.0 it was a fifth of the held pairs. One on each
    primer is still one per window, as it always was."""
    write_lenient_report(tmp_path, [("g0", "a", 2, 2, 0, 0),    # two in the forward window
                                    ("g1", "a", 1, 1, 1, 1)])   # one in each primer's window
    assert read_three_prime_held(tmp_path, budget=1) == {"a": {"g1"}}
    assert read_three_prime_held(tmp_path, budget=2) == {"a": {"g1"}}


def test_an_unreadable_count_is_not_claimed_as_the_terminal_base(tmp_path):
    report = write_lenient_report(tmp_path, [("g0", "a", 1, 1, 0, 0)])
    report.write_text(report.read_text().replace("\t1\t0\t1", "\tNA\t0\t1"))
    assert read_three_prime_held(tmp_path, budget=1) == {}


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
    write_report(tmp_path, [("sample1", "assayA"), ("sample2", "assayA"), ("sample1", "assayB")])
    assert read_report(tmp_path) == {"assayA": {"sample1", "sample2"}, "assayB": {"sample1"}}


def test_read_a_report_with_no_positive(tmp_path):
    write_report(tmp_path, [])
    assert read_report(tmp_path) == {}


def test_a_quote_in_a_genomes_header_does_not_swallow_the_rows_after_it(tmp_path):
    """insilicoPCR writes the fasta description straight through, unquoted. Read with the default quote
    handling, a description beginning with a quote character swallows the rest of its row and every row up
    to the next one -- so amplifications vanish and a verdict moves."""
    write_report(tmp_path, [("g0", "a_assay0"), ("g1", "a_assay0")], description='"strain x')
    assert read_report(tmp_path) == {"a_assay0": {"g0", "g1"}}


def test_a_report_without_the_columns_the_counts_come_from_is_an_error(tmp_path):
    """Read as an empty report, every assay would amplify nothing: a clean-looking answer out of a report
    nothing could be read from."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("sample\tgene\tForwardMismatches\n")  # lower case, as a later version might
    with pytest.raises(PrimerFinderError, match="no Sample or Gene column"):
        read_report(tmp_path)


def test_an_old_releases_report_is_refused(tmp_path):
    """Before 0.8.0 a row was every pairing of a forward and a reverse alignment on a contig -- on one
    real 3.3 Mb genome, 21,014 rows for 28 assays of which 71% were pairings PCR cannot make. primer-finder
    1.9.x re-checked them against the detailed reports; this version counts rows as amplicons, so counting
    that report would call phantom pairings cross-reactions. The EndMismatch columns only those releases
    write are the tell."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    columns = ("Sample", "Gene", "GenomeLocation", "AmpliconSize", "Contig", "Contig Description",
               "ForwardPrimers", "ReversePrimers", "ForwardMismatches", "ReverseMismatches",
               "ForwardEndMismatch", "ReverseEndMismatch")
    report.write_text("\t".join(columns) + "\n" + "g0\ta\t1-100\t100\tchr\t\ta-F\ta-R\t0\t0\t0\t0\n")
    with pytest.raises(PrimerFinderError, match="older than 0.8.0"):
        read_report(tmp_path)


def test_a_report_without_any_column_of_0_8_0_is_refused(tmp_path):
    """The other way a report can be too old to count: the EndMismatch columns are a 0.6.x tell, but a
    report carrying neither those nor anything 0.8.0 writes is just as unknowable, and counting it anyway
    would fail open -- raw pairings read as cross-reactions."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("Sample\tGene\tGenomeLocation\tAmpliconSize\n" + "g0\ta\t1-100\t100\n")
    with pytest.raises(PrimerFinderError, match="none of the columns insilicoPCR 0.8.0 writes"):
        read_report(tmp_path)


def test_a_report_without_the_size_column_is_refused(tmp_path):
    """Every primer file written here declares product sizes, and 0.11.0 answers with a SizeCheck column
    -- on an empty report too. A release from 0.8.0 to 0.10.1 reads `<name>-F size=<bp>` as the primer's
    name, matches no alignment to it and writes an empty table, exit 0: measured on the example with
    0.10.1, 0 rows where 0.11.0 writes 68. Read as a report, that is every exclusion genome silent."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\t".join(COLUMNS[:-1]) + "\n")
    with pytest.raises(PrimerFinderError, match="no SizeCheck column.*older than 0.11.0"):
        read_report(tmp_path)
    with pytest.raises(PrimerFinderError, match="no SizeCheck column"):
        read_three_prime_held(tmp_path, budget=1)


def test_a_report_of_0_11_0_is_refused(tmp_path):
    """0.11.0 reads the primer file and its rows are amplicons, but it lets other bindings through than
    0.12.0 does -- a probe mismatched in its first or last two bases for nothing, none mismatched three
    or four bases in, and in the lenient pass no primer mismatched two to four bases from its 3' end --
    so the same assays get other verdicts and other held counts. Measured on the *Xylella* panel in
    validation/results/2026-11-13_v1.12.0_insilico0130/. The two columns 0.12.0 added are the tell, and
    they are there on an empty report too."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    old = [column for column in COLUMNS if not column.endswith("5PrimeUnaligned")]
    report.write_text("\t".join(old) + "\n" + "g0\ta\t1-100\t100\tchr\t\ta-F\ta-R\t0\t0\t0\t0\t+\ttarget\n")
    with pytest.raises(PrimerFinderError, match="no Forward5PrimeUnaligned.*older than 0.12.0"):
        read_report(tmp_path)
    with pytest.raises(PrimerFinderError, match="older than 0.12.0"):
        read_three_prime_held(tmp_path, budget=1)


def test_an_unexpected_size_is_still_an_amplification(tmp_path):
    """The size check labels a product and removes none, and neither does the reading: a product of
    another size than the one designed is still a product, and on the exclusion group still a
    cross-reaction."""
    report = write_report(tmp_path, [("g0", "a_assay0"), ("g1", "a_assay0")])
    lines = report.read_text().splitlines()
    lines[2] = lines[2].replace("\t100\t", "\t583\t").replace("\ttarget", "\tunexpected")
    report.write_text("\n".join(lines) + "\n")
    assert read_report(tmp_path) == {"a_assay0": {"g0", "g1"}}


def script_that(tmp_path, body: str) -> list[str]:
    launcher = tmp_path / "insilico.sh"
    launcher.write_text(f"#!/bin/sh\n{body}\n")
    launcher.chmod(0o755)
    return [str(launcher)]


def test_check_version_accepts_a_current_release(tmp_path):
    assert check_version(script_that(tmp_path, 'echo "insilicoPCR 0.12.0"')) == "insilicoPCR 0.12.0"
    assert check_version(script_that(tmp_path, 'echo "insilicoPCR 1.2"')) == "insilicoPCR 1.2"


def test_check_version_refuses_a_release_that_does_not_know_itself(tmp_path):
    """0.6.x has no -v at all: it prints an "Unrecognized option" usage dump and exits non-zero. That is
    the answer, and it has to arrive here -- before Primer3 and every blast -- not at the last step."""
    old = script_that(tmp_path, "echo 'Unrecognized option: -v' >&2\nexit 2")
    with pytest.raises(PrimerFinderError, match="needs 0.12.0 or newer"):
        check_version(old)


def test_check_version_refuses_a_release_below_0_8_0(tmp_path):
    with pytest.raises(PrimerFinderError, match="older than the 0.12.0.*raw primer pairings"):
        check_version(script_that(tmp_path, 'echo "insilicoPCR 0.7.0"'))
    with pytest.raises(PrimerFinderError, match="no version"):
        check_version(script_that(tmp_path, 'echo "ready"'))


@pytest.mark.parametrize("release", ["0.8.0", "0.9.0", "0.10.1"])
def test_check_version_refuses_a_release_that_cannot_read_the_declared_size(tmp_path, release):
    """These releases write amplicon reports, know every option a run passes and exit 0 -- having read
    `<name>-F size=<bp>` as a primer's name and so reported nothing. Nothing downstream would stop a
    design from calling every assay silent on the exclusion group, so the version has to."""
    with pytest.raises(PrimerFinderError, match=f"insilicoPCR {release} is older than the 0.12.0.*size="):
        check_version(script_that(tmp_path, f'echo "insilicoPCR {release}"'))


def test_check_version_refuses_a_release_that_counts_what_blast_trimmed(tmp_path):
    """0.11.0 reads the declared size, knows every option and writes amplicons -- and counts a mismatch
    near the end of an oligo by what blast trimmed off the alignment with it. Nothing in its report says
    so before the two columns 0.12.0 added are missed, at the last step; the version says so first, and
    names the primer-finder that was written for it."""
    with pytest.raises(PrimerFinderError, match="insilicoPCR 0.11.0 is older than the 0.12.0.*blast "
                                                "trimmed.*primer-finder 1.12.0"):
        check_version(script_that(tmp_path, 'echo "insilicoPCR 0.11.0"'))


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
    # And the part that matters: a run that exits 0 having written nothing is still a failed run
    stubs(insilico_no_report=True)
    with pytest.raises(PrimerFinderError, match="wrote no consolidated report"):
        run(launcher, genomes, primers, tmp_path / "quiet", threads=1)


def test_a_stale_report_is_set_aside_before_the_run(stubs, tmp_path, insilico_pcr):
    """insilicoPCR 0.8.0 never overwrites a report: beside an existing report.tsv it writes
    report(1).tsv, which the reader never looks at. Without the clearing, a rerun into the same folder
    read the previous run's verdicts under the new run's name. The old report may be hours of insilicoPCR
    a failed rerun would otherwise leave nothing of, so it is renamed rather than deleted -- and a file a
    user saved beside it is not insilicoPCR's to lose, so only the names it writes itself are touched."""
    launcher = resolve(insilico_pcr)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F size=100\nAAAA\n>a-R\nTTTT\n")
    out = tmp_path / "out"
    write_report(out, [("g1", "stale_assay")])
    parent = (out / REPORT).parent
    (parent / "report(1).tsv").write_text((out / REPORT).read_text())
    (parent / "report_baseline.tsv").write_text("saved by hand: not a name insilicoPCR writes\n")
    run(launcher, genomes, primers, out, threads=1)
    assert read_report(out) == {"a": {"g1"}}
    assert sorted(path.name for path in parent.glob("report*")) == [
        "report.tsv", "report.tsv.previous", "report_baseline.tsv"]
    assert "stale_assay" in (parent / "report.tsv.previous").read_text()
    assert "saved by hand" in (parent / "report_baseline.tsv").read_text()


def test_an_old_release_that_rejects_the_command_names_the_upgrade(tmp_path):
    """An insilicoPCR older than 0.8.0 stops at the `-l` this run passes, before writing anything. The
    usage dump it leaves never names a version, so the error has to say what is actually wrong. Measured
    against a real 0.6.1 portable release, which prints exactly this line and exits 2."""
    launcher = tmp_path / "old.sh"
    launcher.write_text("#!/bin/sh\necho 'Unrecognized option: -l' >&2\nexit 2\n")
    launcher.chmod(0o755)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    with pytest.raises(PrimerFinderError, match="older than 0.8.0"):
        run([str(launcher)], genomes, primers, tmp_path / "out", threads=1)


def test_a_failing_run_is_reported(stubs, tmp_path, insilico_pcr):
    stubs(insilico_fail=True)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    with pytest.raises(PrimerFinderError, match="failed"):
        run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", threads=1)


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


def test_the_script_runs_both_groups_and_both_files(tmp_path):
    path = write_script(tmp_path / "run.sh", ["java", "-jar", "x.jar"],
                        [tmp_path / "q.fasta", tmp_path / "p.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=4, mismatches=1)
    text = path.read_text()
    # two files, two groups and the lenient exclusion rerun -- after asking the release its version once
    assert text.count("java -jar x.jar") == 7 and text.count("java -jar x.jar -v") == 1
    assert f"-t 4 -m 1 -l {LONGEST_PRODUCT}" in text
    # The lenient rerun: the terminal base readmitted, the budget unit it consumes, probes pinned strict
    assert text.count(f"-t 4 -m 2 -l {LONGEST_PRODUCT} --probe-mismatches 1 "
                      f"--3p-penalty {LENIENT_PENALTY}") == 2
    assert text.count("_exclusion_3prime") == 5  # the header's word, then two files x (dir= and -o)
    # insilicoPCR writes report(1).tsv beside an existing report, which the find below never names -- the
    # old report is set aside, not deleted, and nothing a user saved beside it is touched
    assert text.count('mv -f "${dir}/report.tsv" "${dir}/report.tsv.previous"') == 6
    assert text.count('rm -f "${dir}/report("*").tsv"') == 6
    assert os.access(path, os.X_OK)


def test_the_script_without_the_three_prime_check_runs_two_groups_only(tmp_path):
    path = write_script(tmp_path / "run.sh", ["x"], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=1,
                        three_prime_check=False)
    text = path.read_text()
    assert "--3p-penalty" not in text and "_exclusion_3prime" not in text


def test_the_script_carries_the_probe_tolerance(tmp_path):
    text = write_script(tmp_path / "run.sh", ["x"], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=1,
                        probe_mismatches=2).read_text()
    assert "--probe-mismatches 2" in text
    # and the lenient rerun keeps the chosen probe budget rather than loosening it with -m
    assert f"-m 2 -l {LONGEST_PRODUCT} --probe-mismatches 2 --3p-penalty {LENIENT_PENALTY}" in text


def test_the_script_without_a_launcher_asks_for_one(tmp_path):
    path = write_script(tmp_path / "run.sh", None, [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=0)
    assert "INSILICO_PCR" in path.read_text()
    ran = subprocess.run(["bash", str(path), str(tmp_path / "out")], capture_output=True, text=True,
                         env={key: value for key, value in os.environ.items() if key != "INSILICO_PCR"})
    assert ran.returncode != 0 and "point INSILICO_PCR at" in ran.stderr


@pytest.mark.parametrize("says, runs", [("insilicoPCR 0.12.0", True), ("insilicoPCR 1.0", True),
                                        ("insilicoPCR 0.11.0", False), ("insilicoPCR 0.10.1", False),
                                        ("insilicoPCR 0.8.0", False), ("", False)])
def test_the_script_refuses_a_release_older_than_the_one_it_was_written_for(tmp_path, says, runs):
    """A design checks the release before it starts; the script is run by hand, later, against whatever
    INSILICO_PCR or the recorded launcher is by then. An older release would exit 0, with empty reports
    or with reports counted on other terms, so the script asks first -- and runs nothing when the answer
    is too old or no answer."""
    launcher = tmp_path / "insilico.sh"
    launcher.write_text("#!/bin/sh\n"
                        f'if [ "$1" = "-v" ]; then echo "{says}"; exit 0; fi\n'
                        f'echo ran >> "{tmp_path}/calls"\n')
    launcher.chmod(0o755)
    path = write_script(tmp_path / "run.sh", [str(launcher)], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=1)
    (tmp_path / "out").mkdir()  # The real program makes it; the script ends by listing what is in it
    ran = subprocess.run(["bash", str(path), str(tmp_path / "out")], capture_output=True, text=True)
    if runs:
        assert ran.returncode == 0, ran.stderr
        assert (tmp_path / "calls").read_text().count("ran") == 3
    else:
        assert ran.returncode == 1 and "0.12.0 or newer is needed" in ran.stderr
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
