"""Handing the designed assays to insilicoPCR, and reading its verdict."""

from __future__ import annotations

import gzip
import os
import shutil
import stat
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.design import Assay, Oligo
from primer_finder.insilico import (
    REPORT,
    Verdict,
    _side_in,
    _trimmed_primer_ends,
    command,
    count_genomes,
    read_report,
    read_terminal_only,
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
    # And the part that matters: a run that exits 0 having written nothing is still a failed run
    stubs(insilico_no_report=True)
    with pytest.raises(PrimerFinderError, match="wrote no consolidated report"):
        run(launcher, genomes, primers, tmp_path / "quiet", threads=1)


def test_a_failing_run_is_reported(stubs, tmp_path, insilico_pcr):
    stubs(insilico_fail=True)
    genomes = tmp_path / "inclusion"
    genomes.mkdir()
    (genomes / "g1.fasta").write_text(">chr\nACGT\n")
    primers = tmp_path / "p.fasta"
    primers.write_text(">a-F\nAAAA\n>a-R\nTTTT\n")
    with pytest.raises(PrimerFinderError, match="failed"):
        run(resolve(insilico_pcr), genomes, primers, tmp_path / "out", threads=1)


PRIMER_LENGTH = 20


def write_hits(folder: Path, rows: list[tuple[str, str, int, int, int]], probe: bool = False,
               end: str = "3") -> Path:
    """A report with the mismatch columns, one row per amplicon:
    (sample, gene, counted mismatches, bases trimmed off the forward primer, off the reverse).

    Also writes the detailed report, which is the only place insilicoPCR says WHICH end of a primer it
    trimmed -- `end` chooses, and the consolidated report reads the same either way. Each amplicon is
    given its own place in the genome so that the two reports can be matched up, as they are for a sample
    holding several amplicons of one assay.
    """
    report = folder / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    columns = ["Sample", "Gene", "Contig", "GenomeLocation", "ForwardPrimers", "ReversePrimers",
               "ForwardMismatches", "ReverseMismatches", "ForwardEndMismatch", "ReverseEndMismatch"]
    detailed: dict[str, list[list[str]]] = {}
    with report.open("w") as fh:
        fh.write("\t".join(columns + (["ProbeMismatches"] if probe else [])) + "\n")
        for number, (sample, gene, mismatches, forward, reverse) in enumerate(rows):
            start = 1000 * number + 1          # Somewhere of its own in the genome
            finish = start + 199
            # A 5' trim takes a base off the outer end of the amplicon, so the reported extent moves in;
            # a 3' trim leaves it where it is. Real insilicoPCR 0.6.1 reports exactly that.
            low_end, high_end = (start + forward, finish - reverse) if end == "5" else (start, finish)
            fields = [sample, gene, "chr", f"{low_end}-{high_end}", f"{gene}-F", f"{gene}-R",
                      str(mismatches), "0", str(-forward), str(-reverse)]
            fh.write("\t".join(fields + (["0"] if probe else [])) + "\n")
            for primer, trimmed in ((f"{gene}-F", forward), (f"{gene}-R", reverse)):
                # Trimmed off the 5' end the alignment starts at the primer's base 1 + n; off the 3' end
                # it stops at its base length - n. Each primer's alignment reaches its end of the
                # amplicon, which is how the binder is told from a second copy of the same site.
                low, high = (1 + trimmed, PRIMER_LENGTH) if end == "5" else (1, PRIMER_LENGTH - trimmed)
                span = high - low
                reverse_primer = primer.endswith("-R")
                where = high_end - span if reverse_primer else low_end
                # The reverse primer aligns reversed, so its positions along the primer descend -- which
                # is what the reader's min/max is for
                on_primer = (high, low) if reverse_primer else (low, high)
                detailed.setdefault(sample, []).append(
                    ["chr", primer, str(PRIMER_LENGTH), str(span + 1),
                     str(where), str(where + span), str(on_primer[0]), str(on_primer[1])])
    for sample, hits in detailed.items():
        path = folder / "detailed_report" / sample / f"{sample}.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as fh:
            fh.write("\t".join(["qseqid", "sseqid", "slen", "length",
                                "qstart", "qend", "sstart", "send"]) + "\n")
            for hit in hits:
                fh.write("\t".join(hit) + "\n")
    return report


def test_an_amplicon_that_needed_a_trimmed_primer_end_is_recognised(tmp_path):
    """insilicoPCR does not count the last two bases of a primer: blast trims an unmatched one off and the
    primer is called bound. An amplification that needed that is one it could not have refused."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0),   # the forward primer's last base
                          ("g1", "a_assay0", 0, 0, 2),   # the reverse primer's last two
                          ("g2", "a_assay0", 0, 0, 0)])  # a clean match: a real amplification
    assert read_terminal_only(tmp_path) == {"a_assay0": {"g0", "g1"}}


REAL_DETAILED_ROWS = """qseqid\tsseqid\tpositive\tmismatch\tgaps\tevalue\tbitscore\tslen\tlength\tqstart\tqend\tqseq\tsstart\tsend\tsseq
chr\tassay-R\t20\t0\t0\t4.28e-08\t40.1\t20\t20\t681\t700\tTTGCCC\t20\t1\tTTGCCC
chr\tassay-F\t19\t0\t0\t1.0e-07\t38.2\t20\t19\t501\t519\tCTGGAT\t1\t19\tCTGGAT
chr\tassay-R\t8\t0\t0\t0.62\t16.4\t20\t8\t944\t951\tAACTAG\t12\t5\tAACTAG
"""


def test_the_trimmed_end_is_read_off_rows_the_real_program_writes(tmp_path):
    """Copied from insilicoPCR 0.6.1's own output for a 3'-end mismatch under the forward primer. Note
    the reverse primer's positions along the primer run *down* (20 to 1) because it aligns reversed,
    while the forward primer's run up -- so the side cannot be read off the first and last values as
    written, and a short spurious hit elsewhere sits among the real ones."""
    sample = tmp_path / "detailed_report" / "g0"
    sample.mkdir(parents=True)
    (sample / "g0.tsv").write_text(REAL_DETAILED_ROWS)
    hits = _trimmed_primer_ends(tmp_path, "g0")
    assert _side_in(hits[("chr", "assay-F")], (501, 700)) == "3"   # 1-19 of 20: its last base is missing
    assert _side_in(hits[("chr", "assay-R")], (501, 700)) == "none"  # 20-1 of 20: all of it aligned
    # The spurious 8-base hit sits outside the amplicon and is not mistaken for the primer that bound
    assert len(hits[("chr", "assay-R")]) == 2


def test_the_primer_that_bound_is_the_one_reaching_the_end_of_the_amplicon(tmp_path):
    """A second copy of a primer's site inside the amplicon is not the primer that bound. Taking the
    longest alignment instead picked it whenever the two were the same length, which read a real 5'-end
    amplification as the 3'-end blind spot -- the misranking this check exists to prevent."""
    binder = (502, 520, 19, "5")        # At the amplicon's start, 5' base trimmed
    inner_copy = (601, 619, 19, "3")    # Further in, same length, 3' base trimmed
    assert _side_in([inner_copy, binder], (502, 700)) == "5"
    assert _side_in([binder, inner_copy], (502, 700)) == "5"


def test_a_quote_in_a_genomes_header_does_not_swallow_the_rows_after_it(tmp_path):
    """insilicoPCR writes the fasta description straight through, unquoted. Read with the default quote
    handling, a description beginning with a quote character swallows the rest of its row and every row up
    to the next one -- so amplifications vanish and a verdict moves."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 0, 0), ("g1", "a_assay0", 0, 0, 0)])
    report = tmp_path / REPORT
    lines = report.read_text().splitlines()
    columns = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    rows[0][columns.index("Contig Description") if "Contig Description" in columns else 2] = '"strain x'
    report.write_text("\n".join(["\t".join(columns)] + ["\t".join(row) for row in rows]) + "\n")
    assert read_report(tmp_path) == {"a_assay0": {"g0", "g1"}}


def test_a_quote_in_a_genomes_header_does_not_lose_the_trimmed_end_tally(tmp_path):
    """The same hazard in the other reader, which decides whether an amplification is excused."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0), ("g1", "a_assay0", 0, 1, 0)])
    report = tmp_path / REPORT
    lines = report.read_text().splitlines()
    columns = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    rows[0][columns.index("Contig")] = '"chr'
    report.write_text("\n".join(["\t".join(columns)] + ["\t".join(row) for row in rows]) + "\n")
    # g0's own row is unusable, but g1's must survive it
    assert read_terminal_only(tmp_path) == {"a_assay0": {"g1"}}


def test_a_report_without_the_columns_the_counts_come_from_is_an_error(tmp_path):
    """Read as an empty report, every assay would amplify nothing: a clean-looking answer out of a report
    nothing could be read from."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("sample\tgene\tForwardMismatches\n")  # lower case, as a later version might
    with pytest.raises(PrimerFinderError, match="no Sample or Gene column"):
        read_report(tmp_path)


def test_a_genome_that_cannot_be_read_is_reported_before_the_design_is_lost(tmp_path):
    """Staging happens after Primer3 and every blast. A gzip that stops early passes the fasta check,
    which reads only the first line, and used to raise a bare EOFError from there."""
    folder, into = tmp_path / "genomes", tmp_path / "staged"
    folder.mkdir()
    whole = gzip.compress(b">chr\nACGTACGTACGT\n")
    (folder / "g1.fasta.gz").write_bytes(whole[: len(whole) // 2])
    with pytest.raises(PrimerFinderError, match="Cannot read .*g1.fasta.gz"):
        stage_genomes(folder, into)


def test_a_difference_at_a_primers_five_prime_end_is_a_real_amplification(tmp_path):
    """`ForwardEndMismatch` says how many bases blast trimmed off a primer, not from which end, and the
    two ends mean opposite things. A difference at the 3' end is the base polymerase extends from, which
    is the blind spot this check exists for. One at the 5' end is a primer that binds and extends anyway,
    so the amplification is real and must not be excused.

    Measured on insilicoPCR 0.6.1, one mismatch under a 20-base forward primer: at its 5' base the
    consolidated report reads `ForwardEndMismatch -1` with the amplicon one base shorter, at its 3' base
    it reads `-1` with the amplicon unchanged. Only the detailed report tells them apart -- positions
    2-20 of the primer against 1-19."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0)], end="5")
    assert read_terminal_only(tmp_path) == {}
    # The same report, with the trim at the other end of the same primer
    other = tmp_path / "three_prime"
    write_hits(other, [("g0", "a_assay0", 0, 1, 0)], end="3")
    assert read_terminal_only(other) == {"a_assay0": {"g0"}}


def test_a_trimmed_end_that_cannot_be_placed_counts_as_a_real_amplification(tmp_path, caplog):
    """Without the detailed report there is no telling which end was trimmed. Excusing the amplification
    is the costly way to be wrong -- it ranks the assay above an honest cross-reaction -- so the benefit
    of the doubt goes the other way, and the run says so."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0)])
    shutil.rmtree(tmp_path / "detailed_report")
    with caplog.at_level("WARNING"):
        assert read_terminal_only(tmp_path) == {}
    assert "could not place on the primer" in caplog.text


def test_an_amplicon_that_needed_a_counted_mismatch_is_not_the_blind_spot(tmp_path):
    """It bound through a mismatch `-m` allowed, so at `-m 0` it would not amplify at all: that is the
    tolerance the user chose talking, not the bases this check cannot see. The answer has to be the same
    at every tolerance, which is only true if those amplicons are left out."""
    write_hits(tmp_path, [("g0", "a_assay0", 1, 1, 0),   # trimmed, but it also needed a mismatch
                          ("g1", "b_assay0", 2, 0, 0)])  # and this one needed two
    assert read_terminal_only(tmp_path) == {}


def test_a_genome_with_one_clean_amplicon_does_not_count_as_terminal(tmp_path):
    """One assay can place several amplicons in the same genome. It only amplifies it through an ignored
    difference if none of the amplicons that need no mismatch is an exact match."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0), ("g0", "a_assay0", 0, 0, 0)])
    assert read_terminal_only(tmp_path) == {}


def test_a_trimmed_amplicon_still_counts_beside_one_that_needed_a_mismatch(tmp_path):
    """The mismatch-bound amplicon is ignored rather than treated as a clean one, so it does not cancel
    the trimmed amplicon sitting beside it."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 2, 0), ("g0", "a_assay0", 1, 0, 0)])
    assert read_terminal_only(tmp_path) == {"a_assay0": {"g0"}}


def test_a_probe_bound_through_a_mismatch_is_not_the_blind_spot_either(tmp_path):
    """In qPCR mode a positive needs the probe too, so a probe that only matched because of the tolerance
    makes the amplification refusable in just the same way."""
    report = write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0)], probe=True)
    report.write_text(report.read_text().replace("\t0\n", "\t2\n"))  # ProbeMismatches=2
    assert read_terminal_only(tmp_path) == {}


def test_a_report_without_the_end_columns_says_so_rather_than_finding_nothing(tmp_path, caplog):
    """An insilicoPCR that stopped reporting the trim would otherwise turn every blind spot into a
    confident `no`, with nothing in the output to say the check had degraded."""
    report = tmp_path / REPORT
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("Sample\tGene\tForwardMismatches\n" + "g0\ta_assay0\t0\n")
    with caplog.at_level("WARNING"):
        assert read_terminal_only(tmp_path) == {}
    assert "ForwardEndMismatch" in caplog.text


def test_a_mismatch_count_that_is_not_a_number_is_left_out_and_reported(tmp_path, caplog):
    """One malformed cell must not lose a whole design: both in silico PCR runs are already done by the
    time this is read."""
    write_hits(tmp_path, [("g0", "a_assay0", 0, 1, 0), ("g1", "a_assay0", 0, 1, 0)])
    report = tmp_path / REPORT
    lines = report.read_text().splitlines()
    columns = lines[0].split("\t")
    rows = [line.split("\t") for line in lines[1:]]
    for row in rows:  # g1's forward mismatch count, whatever column that is
        if row[0] == "g1":
            row[columns.index("ForwardMismatches")] = "NA"
    report.write_text("\n".join(["\t".join(columns)] + ["\t".join(row) for row in rows]) + "\n")
    with caplog.at_level("WARNING"):
        assert read_terminal_only(tmp_path) == {"a_assay0": {"g0"}}
    assert "not a whole number" in caplog.text


def test_the_blind_spot_is_the_same_at_every_tolerance(tmp_path):
    """The set of amplifications insilicoPCR could not have refused cannot depend on what the user allowed
    it to bind through -- that is what makes it the blind spot and not the tolerance."""
    for tolerance in (0, 1, 2, 3):
        folder = tmp_path / f"m{tolerance}"
        # What a higher tolerance adds is rows with counted mismatches; the trimmed ones are always there
        rows = [("g0", "a_assay0", 0, 1, 0), ("g1", "a_assay0", 0, 0, 0)]
        rows += [(f"x{n}", "a_assay0", n, 1, 0) for n in range(1, tolerance + 1)]
        write_hits(folder, rows)
        assert read_terminal_only(folder) == {"a_assay0": {"g0"}}, tolerance


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


def test_a_verdict_on_no_inclusion_genome_at_all_is_not_a_pass():
    """0 of 0 is not every one of them. `fraction` guards the division; `complete` has to guard the
    comparison, or an empty group would read as `yes`."""
    nothing = Verdict(0, 0, 0, 17)
    assert not nothing.complete and not nothing.selective
    assert nothing.label == "no" and nothing.support == 0


def test_the_label_of_a_verdict_this_check_cannot_make():
    """`no` would say the assay amplifies an exclusion genome on evidence worth believing. It does not:
    every amplification rests on bases insilicoPCR does not count, so the check has no verdict to give."""
    assert Verdict(8, 8, 3, 17, exclusion_terminal_only=3).label == "undecided (3' end)"
    assert Verdict(8, 8, 3, 17, exclusion_terminal_only=2).label == "no"   # one is a real amplification
    assert Verdict(7, 8, 2, 17, threshold=0.75, exclusion_terminal_only=2).label == (
        "undecided (3' end, 87%)")                                        # and it is short of the group
    assert Verdict(8, 8, 0, 17).label == "yes"
    assert Verdict(7, 8, 0, 17, threshold=0.75).label == "partial (87%)"


def test_what_the_ranking_takes_from_a_verdict():
    """An assay this check could not refuse is a better bet than one it refused on evidence it can
    defend, and a worse one than an assay it cleared."""
    assert Verdict(8, 8, 0, 17).support == 3
    assert Verdict(7, 8, 0, 17, threshold=0.75).support == 2
    assert Verdict(8, 8, 1, 17, exclusion_terminal_only=1).support == 1
    assert Verdict(8, 8, 1, 17).support == 0
    assert Verdict(4, 8, 0, 17).support == 0  # does not cover the inclusion group at all


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
    """insilicoPCR reads only the top level of the folder it is given and hands each file to blast as it
    stands, so a genome in a subfolder would be skipped and a gzipped one would fail."""
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
