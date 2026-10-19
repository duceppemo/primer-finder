"""Handing the designed assays to insilicoPCR, and reading its verdict.

insilicoPCR (https://github.com/duceppemo/insilicoPCR) is a separate program: a portable Java application
that amplifies in silico against assemblies or reads with BBMap and BLAST+. primer-finder does not reimplement
it. It writes the primer fasta files insilicoPCR expects, and, when told where to find it, runs it over the
inclusion and the exclusion genomes and reads its consolidated report.

Its two rules that shape what is written here:

- an assay is `<name>-F`, `<name>-R` and optionally `<name>-P`;
- probe assays and plain primer-pair assays must not share a primer file, because one probe in the file puts
  the whole report in qPCR mode, where an assay without a probe is never positive.
"""

from __future__ import annotations

import csv
import io
import logging
import os
import shutil
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from primer_finder import PrimerFinderError, seqio, tools

if TYPE_CHECKING:  # Only for the annotations: design.py uses this module, so it cannot be imported here
    from primer_finder.design import Assay

log = logging.getLogger(__name__)

# Where insilicoPCR puts the table that says which assay amplified which sample.
REPORT = Path("consolidated_report") / "report.tsv"
SAMPLE_COLUMN, ASSAY_COLUMN = "Sample", "Gene"
CONTIG_COLUMN, LOCATION_COLUMN = "Contig", "GenomeLocation"
# Where insilicoPCR writes one row per primer alignment, which is the only place it says WHICH end of a
# primer blast trimmed. `<output>/detailed_report/<sample>/<sample>.tsv`.
DETAILED_FOLDER = "detailed_report"
# In it the genome is the query and the primer the subject, so `slen` is the primer's length and
# `sstart`/`send` are positions along the primer -- descending where the primer aligned reversed.
PRIMER_NAME_COLUMNS = {"ForwardEndMismatch": "ForwardPrimers", "ReverseEndMismatch": "ReversePrimers"}
# insilicoPCR reports, per oligo, how many mismatches it bound through and -- for the primers -- how many
# bases blast had to trim off the 3' end to place it. The trim is the interesting one: those bases are not
# counted as mismatches at any tolerance, so an amplicon that needed only the trim is one insilicoPCR could
# not have refused. An amplicon that also needed a counted mismatch is a different thing: at `-m 0` it is
# refused, so it is the tolerance talking and not the blind spot. See docs/wiki/Designing-assays.md.
MISMATCH_COLUMNS = ("ForwardMismatches", "ReverseMismatches", "ProbeMismatches")
END_COLUMNS = ("ForwardEndMismatch", "ReverseEndMismatch")
SIZE_COLUMN = "AmpliconSize"
# A row of the consolidated report is not an amplicon: insilicoPCR pairs **every** alignment of the
# forward primer with every alignment of the reverse primer on a contig, with no check that the two face
# each other or that the distance between them is a product. Measured on one 3.3 Mb Brucella genome with
# 28 IS711-targeted assays: 21,014 rows, of which 9,623 pair two alignments on the same strand and 5,256
# point away from each other -- 71% that PCR could not make -- with sizes from 61 bases to 1.9 Mb. One
# assay alone gave 820 rows from its 38 forward and 38 reverse alignments.
#
# Every such row used to count as "this genome amplified": a false cross-reaction in an exclusion genome,
# and in an inclusion genome a false amplification that can carry an assay over `-p`.
#
# The longest span that can be one product. Well above anything this tool designs for -- Primer3 is asked
# for 70-150 bases by default -- and far below the phantom pairings, whose median was 427,962 bases. A row
# longer than this is a pairing of two unrelated sites rather than an amplicon; one between stays counted,
# which is the cautious direction for an exclusion genome.
LONGEST_PRODUCT = 10_000
# The most bases blast trims off a primer's end rather than counting as a mismatch. Measured on
# insilicoPCR 0.6.1: one unmatched base at either end is trimmed and the primer is called bound, two are,
# three are not reported at all. So an alignment shorter than the primer by more than this did not bind.
LONGEST_TRIM = 2

QPCR_NAME = "assays_qpcr.fasta"
PCR_NAME = "assays_pcr.fasta"
SCRIPT_NAME = "run_insilico_pcr.sh"


@dataclass
class Verdict:
    """What in silico PCR says about one assay.

    `threshold` is the fraction of the inclusion genomes the assay has to amplify to count, which is the
    `-p/--min-inclusion` of the run the regions came from. With the default 1.0 an assay has to amplify
    every one of them; with a lower one, a region that some inclusion genomes lack was acceptable when it
    was found, so an assay on it is acceptable here on the same terms.
    """

    inclusion_amplified: int
    inclusion_total: int
    exclusion_amplified: int
    exclusion_total: int
    threshold: float = 1.0
    exclusion_terminal_only: int = 0  # Of those, the ones that rest on an uncounted 3'-end difference

    @property
    def fraction(self) -> float:
        """The share of the inclusion genomes it amplified."""
        return self.inclusion_amplified / self.inclusion_total if self.inclusion_total else 0.0

    @property
    def percent(self) -> int:
        """That share as a whole number, rounded down so that it never overstates the coverage.

        Taken from the two counts rather than from `fraction`: the division to a float loses the exact
        value, and `int(29 / 100 * 100)` is 28, so an assay that amplified 29 of 100 genomes was labelled
        `partial (28%)` -- below a `-p 0.29` the user had asked for. Twelve such pairs occur under 300
        genomes.
        """
        if not self.inclusion_total:
            return 0
        return self.inclusion_amplified * 100 // self.inclusion_total

    @property
    def complete(self) -> bool:
        """Amplifies every inclusion genome and no exclusion genome."""
        return (self.inclusion_total > 0 and self.inclusion_amplified == self.inclusion_total
                and self.exclusion_amplified == 0)

    @property
    def selective(self) -> bool:
        """Good enough for the threshold the regions were found with, and silent on the exclusion group."""
        return self.exclusion_amplified == 0 and self.fraction >= self.threshold

    @property
    def only_terminal(self) -> bool:
        """It covers the inclusion group, and every exclusion genome it amplified rests on a difference in
        the last two bases of a primer -- which insilicoPCR does not count. The `no` is then the model's
        blind spot: the assay may well discriminate, and in silico PCR cannot say either way."""
        return (self.fraction >= self.threshold and self.exclusion_amplified > 0
                and self.exclusion_terminal_only == self.exclusion_amplified)

    @property
    def label(self) -> str:
        """What the table says: `yes`, `partial (88%)` when it meets a threshold below 1,
        `undecided (3' end)` when the only exclusion amplifications are ones this check cannot refuse,
        or `no`."""
        if self.complete:
            return "yes"
        if self.selective:
            return f"partial ({self.percent}%)"
        if self.only_terminal:
            return "undecided (3' end)" if self.percent == 100 else f"undecided (3' end, {self.percent}%)"
        return "no"

    @property
    def support(self) -> int:
        """How far this check backs taking the assay to a bench, for the ranking: it amplified the whole
        inclusion group and nothing else (3), it met the run's threshold and nothing else (2), it amplified
        exclusion genomes but only through differences this check cannot see (1), or it amplified them on
        evidence the check can defend (0)."""
        if self.complete:
            return 3
        if self.selective:
            return 2
        return 1 if self.only_terminal else 0

    def __str__(self) -> str:
        return (f"inclusion {self.inclusion_amplified}/{self.inclusion_total} ({self.percent}%), "
                f"exclusion {self.exclusion_amplified}/{self.exclusion_total}")


def write_primer_fasta(assays: Sequence[Assay], path: Path, kind: str) -> int:
    """Write the assays as insilicoPCR reads them, in one of its two modes.

    insilicoPCR cannot run both modes at once: one probe in a file puts the whole report in qPCR mode,
    where an assay without a probe is never positive. So the same assays are written twice, and run twice:

    - `qpcr`: the assays that have a probe, with it. A positive needs the probe to bind as well.
    - `pcr`: every assay, primers only. This is what says whether the primer pair discriminates on its own,
      which is what matters when a mismatch under a primer is all the specificity there is.
    """
    if kind not in ("qpcr", "pcr"):  # pragma: no cover - the caller is this package
        raise PrimerFinderError(f"Unknown assay type {kind!r}")
    chosen = [assay for assay in assays if kind == "pcr" or assay.probe is not None]
    if not chosen:
        return 0
    records: list[seqio.Record] = []
    for assay in chosen:
        records.append(seqio.Record(f"{assay.name}-F", "", assay.forward.sequence))
        records.append(seqio.Record(f"{assay.name}-R", "", assay.reverse.sequence))
        if kind == "qpcr" and assay.probe is not None:
            records.append(seqio.Record(f"{assay.name}-P", "", assay.probe.sequence))
    seqio.write_fasta(path, records)
    return len(chosen)


STAGED = "insilico_genomes"


def stage_genomes(folder: Path, into: Path) -> int:
    """Lay a genome folder out the way insilicoPCR needs it: one flat folder of plain fasta files.

    It reads only the top level of what it is given, and hands each file to its own blast as it stands, so
    a genome in a subfolder is silently skipped and a gzipped one makes that blast call fail -- which takes
    the whole run with it. primer-finder accepts both, and the bundled example holds a `.fasta.gz`, so the
    genomes are staged rather than passed through: subfolders are flattened, gzipped genomes are written
    out plain, and a name used twice is numbered so that two files in different subfolders stay two
    samples.

    Returns how many genomes were staged, which is what the verdicts are counted against.
    """
    # This folder is emptied below, so it must not be, or contain, the folder being read from
    here, there = into.resolve(), folder.resolve()
    if here == there or there in here.parents or here in there.parents:
        raise PrimerFinderError(
            f"Refusing to stage {folder} into {into}: one is inside the other, and staging empties the "
            "folder it writes to, which would delete the genomes it is about to read."
        )
    into.mkdir(parents=True, exist_ok=True)
    # Clear what an earlier run staged. Only the names about to be written were replaced before, so a
    # genome that has since left the input folder stayed behind -- and `count_genomes` on this folder is
    # the denominator of every in-silico verdict, so a stale file made a selective assay read `partial`.
    for stale in into.glob("*.fasta"):
        stale.unlink()
    used: set[str] = set()
    for path in seqio.find_genomes(folder):
        stem = seqio.base_name(path)
        name, number = stem, 1
        while name in used:
            number += 1
            name = f"{stem}_{number}"
        used.add(name)
        target = into / f"{name}.fasta"
        if target.is_symlink() or target.exists():
            target.unlink()
        if path.suffix == ".gz":
            # A truncated or trailing-byte gzip passes the fasta check, which reads only the first line,
            # and then raises here. Left alone that is an EOFError out of a command that has already run
            # Primer3 and every blast, and the assays would never be written.
            try:
                with seqio.open_text(path) as source, target.open("w") as plain:
                    shutil.copyfileobj(source, plain)
            except (OSError, EOFError) as error:
                raise PrimerFinderError(f"Cannot read {path}: {error}") from error
        else:
            target.symlink_to(path.resolve())
    return len(used)


def resolve(given: Path) -> list[str]:
    """The command that starts insilicoPCR, from what the user pointed at: the folder of an extracted
    portable release, its jar, or a launcher script."""
    given = Path(given)
    if given.is_dir():
        jar = given / "insilicoPCR.jar"
        java = given / "runtime" / "linux" / "jdk" / "bin" / "java"
        if not jar.is_file():
            raise PrimerFinderError(f"No insilicoPCR.jar in {given}")
        if java.is_file():  # The bundled runtime of the portable release
            return [str(java), "-jar", str(jar)]
        return [_java(), "-jar", str(jar)]
    if given.suffix == ".jar":
        if not given.is_file():
            raise PrimerFinderError(f"No such file: {given}")
        return [_java(), "-jar", str(given)]
    if not given.is_file():
        raise PrimerFinderError(f"No such file: {given}")
    if not os.access(given, os.X_OK):
        raise PrimerFinderError(f"Not executable: {given}")
    return [str(given)]


def _java() -> str:
    java = shutil.which("java")
    if java is None:
        raise PrimerFinderError(
            "insilicoPCR needs Java, which is not on PATH. Point --insilico-pcr at the folder of an "
            "extracted portable release instead: it brings its own Java runtime."
        )
    return java


def command(launcher: Sequence[str], genomes: Path, primers: Path, output: Path,
            threads: int, mismatches: int) -> list[str]:
    return [*launcher, "-i", str(genomes), "-p", str(primers), "-o", str(output),
            "-t", str(threads), "-m", str(mismatches)]


def run(launcher: Sequence[str], genomes: Path, primers: Path, output: Path,
        threads: int, mismatches: int = 0) -> Path:
    """Run insilicoPCR over a folder of genomes. Returns the folder it wrote."""
    output.mkdir(parents=True, exist_ok=True)
    tools.run(command(launcher, genomes, primers, output, threads, mismatches))
    report = output / REPORT
    if not report.is_file():
        raise PrimerFinderError(f"insilicoPCR wrote no consolidated report ({report} missing)")
    return output


def read_report(output: Path, check_shape: bool = True) -> dict[str, set[str]]:
    """Which samples each assay amplified, out of insilicoPCR's consolidated report.

    A row is only counted when it could be a product. insilicoPCR pairs every forward-primer alignment
    with every reverse-primer alignment on a contig (see LONGEST_PRODUCT), so most rows of a repeat-rich
    genome are two unrelated sites rather than an amplicon. Two things rule those out, both from what the
    reports already say:

    - the span has to be short enough to be one product, which the consolidated report gives directly;
    - the two primers have to face each other, which the detailed report gives. That one is read per
      sample and only for the samples a row survives the first test for, so a specific assay pays almost
      nothing and the cost falls where the answer is in doubt.

    `check_shape=False` counts every row, which is what the readers did before 1.9.2.
    """
    report = output / REPORT
    amplified: dict[str, set[str]] = {}
    alignments: dict[str, dict[tuple[str, str], list]] = {}
    too_long = bent = 0
    with report.open(newline="") as fh:
        # QUOTE_NONE: insilicoPCR writes the fasta description straight through, so a genome whose header
        # holds a quote character would otherwise swallow the rest of its row and every row up to the next
        # one -- losing amplifications silently, which moves a verdict.
        rows = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        _insist_on_columns(rows.fieldnames, report)
        for row in rows:
            assay, sample = (row.get(ASSAY_COLUMN) or "").strip(), (row.get(SAMPLE_COLUMN) or "").strip()
            if not assay or not sample:
                continue
            if check_shape:
                size = _count(row, SIZE_COLUMN)
                contig = (row.get(CONTIG_COLUMN) or "").strip()
                span = _span((row.get(LOCATION_COLUMN) or "").strip())
                if size is not None and size > LONGEST_PRODUCT:
                    too_long += 1
                    continue
                if not _is_a_product(output, sample, contig, span, size, alignments):
                    bent += 1
                    continue
            amplified.setdefault(assay, set()).add(sample)
    if too_long or bent:
        log.info("%d row(s) of %s were two unrelated primer sites rather than an amplicon and are not "
                 "counted as an amplification: %d spanned more than %d bases, %d had the two primers "
                 "pointing the wrong way", too_long + bent, report, too_long, LONGEST_PRODUCT, bent)
    return amplified


def _insist_on_columns(fieldnames: Sequence[str] | None, report: Path) -> None:
    """The two columns every verdict is counted from. Without them every assay would read as amplifying
    nothing at all, which looks like a clean answer rather than an unreadable report."""
    missing = [name for name in (SAMPLE_COLUMN, ASSAY_COLUMN) if name not in (fieldnames or ())]
    if missing:
        raise PrimerFinderError(
            f"insilicoPCR's report has no {' or '.join(missing)} column, so there is no telling what "
            f"amplified: {report}"
        )


def _count(row: dict[str, str], column: str) -> int | None:
    """One of insilicoPCR's mismatch counts. A column it does not write, or leaves blank, is nothing
    reported and so a zero; anything that is not a whole number is None, which the caller must not read as
    one more than it says."""
    value = (row.get(column) or "").strip()
    if not value:
        return 0
    try:
        return int(value)
    except ValueError:
        return None


def read_terminal_only(output: Path) -> dict[str, set[str]]:
    """Which samples each assay amplified *only* through a difference in the last two bases of a primer.

    insilicoPCR does not count those bases: blast trims an unmatched base off the end of its alignment and
    the primer is called bound. That happens at every `-m`, so such an amplicon is one the tool could not
    have refused, and a verdict that rests entirely on them is the model's blind spot rather than a
    cross-reaction.

    Only the amplicons with no counted mismatch anywhere decide this, which is what makes the answer the
    same at every `-m`:

    - one of those with a trimmed primer end, and none without, is the blind spot;
    - one of those with nothing trimmed is an exact match, so the amplification is real and the sample does
      not count however many trimmed amplicons sit beside it;
    - an amplicon that needed a counted mismatch is the `-m` the user chose talking, not these bases, and is
      ignored here. A sample with nothing but those does not count: at `-m 0` it would not amplify at all.
    """
    terminal: set[tuple[str, str]] = set()
    exact: set[tuple[str, str]] = set()
    unreadable = 0
    untellable = 0
    shapeless = 0
    # One sample's primer alignments: the shape test and the trimmed-end test both read them, and the
    # report is grouped by sample, so one entry serves both.
    alignments: dict[str, dict[tuple[str, str], list]] = {}
    with (output / REPORT).open(newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        if rows.fieldnames is not None and not any(name in rows.fieldnames for name in END_COLUMNS):
            log.warning("insilicoPCR's report has no %s column, so which exclusion amplifications rest on "
                        "an ignored 3'-end difference cannot be told: %s",
                        " or ".join(END_COLUMNS), output / REPORT)
            return {}
        for row in rows:
            assay, sample = (row.get(ASSAY_COLUMN) or "").strip(), (row.get(SAMPLE_COLUMN) or "").strip()
            if not assay or not sample:
                continue
            contig = (row.get(CONTIG_COLUMN) or "").strip()
            span = _span((row.get(LOCATION_COLUMN) or "").strip())
            # The same test `read_report` applies: a row that cannot be a product is not an amplicon
            # here either. Without it a phantom pairing that happens to carry a 3' trim flipped a real
            # exact cross-reaction into "undecided (3' end)" -- excusing an assay as the model's blind
            # spot on the strength of two unrelated primer sites -- and a phantom *exact* row still
            # masked a real trimmed amplicon, which is what this check was added for.
            if not _is_a_product(output, sample, contig, span, _count(row, SIZE_COLUMN), alignments):
                shapeless += 1
                continue
            counted = [_count(row, name) for name in MISMATCH_COLUMNS]
            ends = [_count(row, name) for name in END_COLUMNS]
            if None in counted or None in ends:
                unreadable += 1
                continue
            if any(counted):  # it bound through a mismatch this check does count
                continue
            if not any(end < 0 for end in ends):
                exact.add((assay, sample))
                continue
            # Something was trimmed, but the report does not say which end of the primer: a difference at
            # the 5' end reads exactly the same here, and that one does not hinder extension, so the
            # amplification is real. Which end it was has to come from the detailed report.
            ends_here = alignments.setdefault(sample, _trimmed_primer_ends(output, sample))
            sides = set()
            for name, end in zip(END_COLUMNS, ends, strict=True):
                if end >= 0:
                    continue
                primer = (row.get(PRIMER_NAME_COLUMNS[name]) or "").strip()
                sides.add(_side_in(ends_here.get((contig, primer), []), span))
            if sides == {"3"}:  # Only the blind spot, so this amplicon is one insilicoPCR could not refuse
                terminal.add((assay, sample))
            else:
                # A 5' trim, or an end we could not place: either way not something to excuse. Counting it
                # as a real amplification is the cautious way to be wrong.
                untellable += "unknown" in sides
                exact.add((assay, sample))
    if shapeless:
        log.info("%d row(s) of %s could not be one product and were left out of the 3'-end difference "
                 "tally", shapeless, output / REPORT)
    if unreadable:
        log.warning("%d row(s) of %s give a mismatch count that is not a whole number and were left out of "
                    "the 3'-end difference tally", unreadable, output / REPORT)
    if untellable:
        log.warning("%d amplicon(s) had a primer end trimmed that %s could not place on the primer, and "
                    "are counted as real amplifications rather than as the 3'-end blind spot",
                    untellable, output / DETAILED_FOLDER)
    found: dict[str, set[str]] = {}
    for assay, sample in terminal - exact:
        found.setdefault(assay, set()).add(sample)
    return found


def _trimmed_primer_ends(output: Path, sample: str) -> dict[tuple[str, str], list[tuple[int, int, int, str]]]:
    """Which end of each primer blast left out of its alignment, per contig and primer, for one sample.

    "3" where the primer's last base did not align, "5" where its first did not, "both", or "none". The
    consolidated report's `ForwardEndMismatch` only says how many bases were trimmed, not from which end,
    and the difference decides whether an amplification is the model's blind spot or a real one: a primer
    whose 3' base does not match is the base polymerase extends from, while one whose 5' base does not is
    a primer that binds and extends anyway.

    Measured on insilicoPCR 0.6.1 with one mismatch under a 20-base primer: at its 5' base the alignment
    covers positions 2-20 of the primer, at its 3' base positions 1-19, and the consolidated report reads
    `-1` for both.
    """
    detailed = output / DETAILED_FOLDER / sample / f"{sample}.tsv"
    if not detailed.is_file():
        return {}
    found: dict[tuple[str, str], list[tuple[int, int, int, str]]] = {}
    try:
        with detailed.open(newline="") as fh:
            text = fh.read()
    except (OSError, UnicodeDecodeError) as problem:
        # A file written by another program, in a folder the user controls: unreadable or not text is
        # evidence we do not have, not a reason to lose a whole design run after insilicoPCR has already
        # run. Every caller treats "no alignments" as "cannot tell", which counts the amplicon.
        log.warning("%s could not be read (%s), so which end of a primer was trimmed in %s is not known",
                    detailed, problem, sample)
        return {}
    with io.StringIO(text, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            try:
                length, primer_length = int(row["length"]), int(row["slen"])
                on_primer = (int(row["sstart"]), int(row["send"]))
                in_genome = (int(row["qstart"]), int(row["qend"]))
            except (KeyError, TypeError, ValueError):
                continue
            low, high = min(on_primer), max(on_primer)
            side = {(True, True): "none", (False, True): "5", (True, False): "3",
                    (False, False): "both"}[(low <= 1, high >= primer_length)]
            # Only an alignment that covers nearly the whole primer can be the one that bound: blast trims
            # at most the last base or two, and a mismatch the run counts is aligned rather than trimmed.
            # Dropping the rest is what makes reading these reports affordable -- insilicoPCR's default
            # `-e 1e5` fills them with 7- and 8-base matches, 569,373 rows for one real genome, of which
            # 336 keys were bindings.
            if length < primer_length - LONGEST_TRIM:
                continue
            # Ascending on the primer means the primer matched the genome as written, descending means
            # reversed. That is the only orientation insilicoPCR gives, and a product needs one of each.
            ascending = on_primer[0] <= on_primer[1]
            key = ((row.get("qseqid") or "").strip(), (row.get("sseqid") or "").strip())
            found.setdefault(key, []).append((min(in_genome), max(in_genome), length, side, ascending))
    return found


def _is_a_product(output: Path, sample: str, contig: str, span: tuple[int, int] | None,
                  size: int | None, cache: dict[str, dict[tuple[str, str], list]]) -> bool:
    """Whether a consolidated-report row could be one product, by the two tests the reports allow.

    Both readers ask this. `read_terminal_only` had no such test until 1.9.2, so a phantom pairing that
    happened to carry a 3' trim excused an assay as the model's blind spot, and a phantom *exact* row
    masked a real trimmed amplicon.

    `cache` holds one sample's primer alignments. The consolidated report is written grouped by sample, so
    one entry is all that is ever wanted: keeping every sample's reached tens of GB on the 951- and
    1,290-genome panels of the Brucella records.
    """
    if size is not None and size > LONGEST_PRODUCT:
        return False
    if sample not in cache:
        cache.clear()
        cache[sample] = _trimmed_primer_ends(output, sample)
    here = [hit for (hit_contig, _), hits in cache[sample].items()
            if hit_contig == contig for hit in hits]
    return _convergent(span, here)


def _convergent(span: tuple[int, int] | None, alignments: Sequence[Sequence[object]]) -> bool:
    """Whether two of these primer alignments could have made the amplicon covering `span`.

    A product needs a primer at each end of the span pointing inwards: the one at the low end matched the
    genome as written and the one at the high end matched it reversed. Both orientations of the genome
    record give that same shape -- when the record is the reverse complement the two primers swap ends,
    and each swaps direction with it -- which is what makes this the test rather than "forward on the
    left". Measured on insilicoPCR 0.6.x output for an amplicon and for the same amplicon in a
    reverse-complemented genome.

    It takes positive evidence to say no. Without a span, or with no alignment reaching either end of it
    -- no detailed report, or a binder this dropped as too short -- the answer is True: an amplicon that
    cannot be placed is counted as real, which is the cautious direction for an exclusion genome. The
    answer is False only where an alignment does reach a boundary and faces the wrong way.
    """
    if span is None:
        return True
    low, high = span
    at_low = [hit for hit in alignments if hit[0] == low]
    at_high = [hit for hit in alignments if hit[1] == high]
    # The primer at the low end must have matched the genome as written and the one at the high end
    # reversed; an alignment there facing the other way is a site that points away from the span.
    outward_low = bool(at_low) and not any(hit[4] for hit in at_low)
    outward_high = bool(at_high) and not any(not hit[4] for hit in at_high)
    return not (outward_low or outward_high)


def _side_in(hits: list[tuple[int, int, int, str]], span: tuple[int, int] | None) -> str:
    """Which end was trimmed off the alignment that produced the amplicon covering `span`.

    One assay can place several amplicons in one genome, so the primer's alignments are matched to the
    amplicon by where they sit in the genome; the longest alignment inside the amplicon is the one that
    bound. Without a span to go on, or with no alignment inside it, the answer is not known.
    """
    if not hits:
        return "unknown"
    if span is None:
        return "unknown" if len({side for *_, side in hits}) > 1 else hits[0][3]
    low, high = span
    inside = [hit for hit in hits if low <= hit[0] and hit[1] <= high]
    if not inside:
        return "unknown"
    # The amplicon is reported as the extent of the two alignments that made it, so the primer that bound
    # reaches one of its ends. Any other alignment of the same primer inside the amplicon is a second copy
    # of that site, and taking the longest without regard to where it sits picks it whenever the two are
    # the same length: measured on real output, a genome with a 5'-trimmed F site at the amplicon's start
    # and an untouched-length 3'-trimmed copy of F further in was read as the 3' blind spot and excused,
    # which is the misranking this whole check exists to prevent.
    touching = [hit for hit in inside if hit[0] == low or hit[1] == high]
    return max(touching or inside, key=lambda hit: hit[2])[3]


def _span(text: str) -> tuple[int, int] | None:
    """insilicoPCR writes an amplicon's place in the genome as "501-700"."""
    first, _, second = text.partition("-")
    try:
        low, high = int(first), int(second)
    except ValueError:
        return None
    return (low, high) if low <= high else (high, low)


def count_genomes(folder: Path) -> int:
    return len(seqio.find_genomes(folder))


def verdicts(assays: Iterable[Assay], inclusion: Path, exclusion: Path,
             inclusion_report: dict[str, set[str]], exclusion_report: dict[str, set[str]],
             threshold: float = 1.0,
             exclusion_terminal: dict[str, set[str]] | None = None) -> dict[str, Verdict]:
    """Put the two reports together, one verdict per assay."""
    inclusion_total, exclusion_total = count_genomes(inclusion), count_genomes(exclusion)
    terminal = exclusion_terminal or {}
    return {
        assay.name: Verdict(
            inclusion_amplified=len(inclusion_report.get(assay.name, set())),
            inclusion_total=inclusion_total,
            exclusion_amplified=len(exclusion_report.get(assay.name, set())),
            exclusion_total=exclusion_total,
            threshold=threshold,
            exclusion_terminal_only=len(terminal.get(assay.name, set())),
        )
        for assay in assays
    }


def write_script(path: Path, launcher: Sequence[str] | None, primers: Sequence[Path],
                 inclusion: Path, exclusion: Path, threads: int, mismatches: int) -> Path:
    """A script that runs insilicoPCR over both groups, for running the check by hand later."""
    start = " ".join(launcher) if launcher else (
        '"${INSILICO_PCR:?point INSILICO_PCR at the folder of an extracted insilicoPCR release}"'
        '/runtime/linux/jdk/bin/java -jar "${INSILICO_PCR}"/insilicoPCR.jar'
    )
    lines = [
        "#!/usr/bin/env bash",
        "# In silico PCR of the designed assays against both groups of genomes, with insilicoPCR:",
        "#   https://github.com/duceppemo/insilicoPCR",
        "# An assay is selective when it amplifies every inclusion genome and no exclusion genome.",
        "set -euo pipefail",
        "",
        f'out="${{1:-{path.parent}/insilico_pcr}}"',
        "",
    ]
    for primer_file in primers:
        for group, folder in (("inclusion", inclusion), ("exclusion", exclusion)):
            lines += [
                f"{start} \\",
                f'    -i "{folder}" \\',
                f'    -p "{primer_file}" \\',
                f'    -o "${{out}}/{primer_file.stem}_{group}" \\',
                f"    -t {threads} -m {mismatches}",
                "",
            ]
    lines += [
        'echo "Consolidated reports:"',
        'find "${out}" -name report.tsv',
        "",
    ]
    path.write_text("\n".join(lines))
    path.chmod(path.stat().st_mode | 0o111)
    return path
