"""Handing the designed assays to insilicoPCR, and reading its verdict."""

from __future__ import annotations

import gzip
import os
import stat
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.design import Assay, Oligo
from primer_finder.insilico import (
    LONGEST_PRODUCT,
    REPORT,
    Verdict,
    command,
    count_genomes,
    read_report,
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


# What insilicoPCR 0.8.0 writes; primer-finder reads Sample and Gene and insists the rest is not an older
# release's report.
COLUMNS = ("Sample", "Gene", "GenomeLocation", "AmpliconSize", "Contig", "Contig Description",
           "ForwardPrimers", "ReversePrimers", "ForwardMismatches", "ReverseMismatches",
           "Forward3PrimeMismatches", "Reverse3PrimeMismatches", "Strand")


def write_report(folder: Path, rows: list[tuple[str, str]], description: str = "") -> Path:
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for sample, gene in rows:
            fh.write(f"{sample}\t{gene}\t1-100\t100\tchr\t{description}"
                     f"\t{gene}-F\t{gene}-R\t0\t0\t0\t0\t+\n")
    return report


def test_the_two_kinds_of_primer_file(tmp_path):
    """One file per question: with the probe in the file a probe assay is only positive when the probe
    binds, so whether its primer pair discriminates alone needs a file with the probe stripped."""
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


def test_a_stale_report_is_cleared_before_the_run(stubs, tmp_path, insilico_pcr):
    """insilicoPCR 0.8.0 never overwrites a report: beside an existing report.tsv it writes
    report(1).tsv, which the reader never looks at. Without the clearing, a rerun into the same folder
    read the previous run's verdicts under the new run's name."""
    launcher = resolve(insilico_pcr)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    out = tmp_path / "out"
    write_report(out, [("g1", "stale_assay")])
    (out / REPORT).with_name("report(1).tsv").write_text((out / REPORT).read_text())
    run(launcher, genomes, primers, out, threads=1)
    assert read_report(out) == {"a": {"g1"}}
    assert [path.name for path in (out / REPORT).parent.glob("report*.tsv")] == ["report.tsv"]


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
    assert text.count("java -jar x.jar") == 4  # two files, two groups
    assert f"-t 4 -m 1 -l {LONGEST_PRODUCT}" in text
    assert "--probe-mismatches" not in text  # not chosen, so insilicoPCR's own default stands
    # insilicoPCR writes report(1).tsv beside an existing report, which the find below never names
    assert text.count("rm -f") == 4
    assert os.access(path, os.X_OK)


def test_the_script_carries_the_probe_tolerance(tmp_path):
    path = write_script(tmp_path / "run.sh", ["x"], [tmp_path / "q.fasta"],
                        tmp_path / "inclusion", tmp_path / "exclusion", threads=2, mismatches=1,
                        probe_mismatches=2)
    assert "--probe-mismatches 2" in path.read_text()


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
