"""Handing the designed assays to insilicoPCR, and reading its verdict.

insilicoPCR (https://github.com/duceppemo/insilicoPCR) is a separate program: a portable Java application
that amplifies in silico against assemblies or reads with BBMap and BLAST+. primer-finder does not reimplement
it, and copies none of its code: it writes the primer fasta files insilicoPCR expects, and, when told where
to find it, runs it over the inclusion and the exclusion genomes and reads its consolidated report.

This module speaks insilicoPCR 2.1.0 or newer. A row of its consolidated report is a called amplicon --
two accepted primer sites on one contig, facing each other, within a size limit, and for a probe assay with
an accepted probe site inside the product -- so counting rows is counting amplifications.

Since 2.0 a primer or probe site is accepted when its melting temperature on the template, Primer3's
thermodynamic alignment worked out on the bases under the site, reaches the annealing temperature less a
margin, and, for a primer, when its 3' end passes the weighted 3' rule. That one rule replaced the
mismatch budgets of 1.x (`-m`, `--probe-mismatches`), which 2.0 refuses on the command line. The reaction
the Tm is worked out for is the one Primer3 designed the oligos in, passed through, so both programs judge
the oligos in one reaction. Older releases are refused up front: their command line is not this one, and
their reports, which name the assay `Gene` and carry no Tm, count other bindings.

Five rules of insilicoPCR shape what is written and read here:

- an assay is `<name>-F`, `<name>-R` and optionally `<name>-P`;
- what follows the name on a header is a description, of which `size=<bp>` declares the product size;
- an assay whose probe is in the primer file is only positive when the probe binds inside the product.
  Whether its primer pair discriminates *alone* is a different question, so the assays are also written
  with the probes stripped, and that file is run separately;
- a product can be the designed pair, or one primer alone binding both strands (`Pairing` `one-primer`):
  a real product in the tube, so it counts on the exclusion group, but not the target the inclusion group
  is asked about;
- a genome insilicoPCR cannot read is left out of the run, with exit 0, and is listed as `not-analysed`
  in `report_samples.tsv`. Counted as not amplifying, it would fake selectivity or a miss, so the check
  fails instead.
"""

from __future__ import annotations

import csv
import logging
import os
import re
import shlex
import shutil
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from primer_finder import PrimerFinderError, seqio, tools

if TYPE_CHECKING:  # Only for the annotations: design.py uses this module, so it cannot be imported here
    from primer_finder.design import Assay

log = logging.getLogger(__name__)

RELEASES = "https://github.com/duceppemo/insilicoPCR/releases"
# Where insilicoPCR puts the table that says which assay amplified which sample. The run's own report is
# the one its standard output names (REPORT_WRITTEN); this is where it lands in a folder of its own.
REPORT = Path("consolidated_report") / "report.tsv"
REPORT_WRITTEN = re.compile(r"^Report written to (.+?)\s*$", re.MULTILINE)
SAMPLE_COLUMN, ASSAY_COLUMN = "Sample", "Assay"
# What insilicoPCR named the assay column before 1.0: a report that has it is from a release that
# accepted sites by mismatch counts.
OLD_ASSAY_COLUMN = "Gene"
# How a product's primers came together (since 1.0): `designed` is the assay's forward primer with its
# reverse primer, the product the inclusion group is asked about. `one-primer` (one primer binding both
# strands, at an inverted repeat) and `two-assays` (assays sharing a reaction, which primer-finder never
# declares) are real products of the tube.
PAIRING_COLUMN, DESIGNED = "Pairing", "designed"
# Each oligo's Tm on its site (since 2.0). The primers' are always filled; the probe's on probe assays.
TM_COLUMNS = ("ForwardTm", "ReverseTm", "ProbeTm")
# The share of each primer's imperfections inside its 3' window.
WINDOW_COLUMNS = ("Forward3PrimeMismatches", "Reverse3PrimeMismatches")
# The columns this module reads, which a report must have rather than be read as empty.
REQUIRED_COLUMNS = (SAMPLE_COLUMN, ASSAY_COLUMN, PAIRING_COLUMN, *TM_COLUMNS, *WINDOW_COLUMNS)
# The files insilicoPCR writes for a numbered report (`report(1).tsv`, `report(1).html`,
# `report(1)_samples.tsv`, ...), and nothing else a user may have saved beside them.
NUMBERED_REPORT = re.compile(r"report\(\d+\)(\.tsv|\.html|_[a-z_]+\.tsv)")
# A cell insilicoPCR leaves without an apostrophe although it starts with `+` or `-`: a number, by the
# test its ReportColumn.formulaLike makes.
NUMBER = re.compile(r"[+-]\d+(\.\d+)?([eE][+-]?\d+)?")
# The oldest insilicoPCR this module runs. 2.0 is the first to accept a site by its Tm; 2.1.0 is the
# release this module was checked against, and the first whose report of the assays without a product
# counts a mismatch among a primer's first five bases once (2.0.0 counted it twice).
REQUIRED_VERSION = (2, 1, 0)
TM_VERSION = (2, 0, 0)
# The longest span that can be one product, passed to insilicoPCR as `-l`. Its default (5000) is its own
# to move, and the verdicts should not move with it. Well above anything this tool designs for --
# Primer3 is asked for 70-150 bases by default.
LONGEST_PRODUCT = 10_000
# What the 3'-end dependence check hands insilicoPCR as --3p-penalty: 4, the price of the worst single
# mismatch its weighted window knows, a G/C pair broken on the 3'-TERMINAL base (an A/T pair scores 3
# there). The lenient run differs from the strict one in this alone -- the Tm rule is the same in both --
# so the sites it accepts are the strict run's and those whose 3' window carries an imperfection that
# weighs no more than that; and a readmitted site must still melt at the minimum, the mismatch lowering its
# Tm as any other does. A genome only the lenient run amplifies is held off by the 3' rule and nothing
# else. The penalty also readmits two mismatches in one window when they come to no more than a terminal
# one -- measured, at the second and fourth bases from the 3' end, and at the third and fourth, whose
# sites still melt at the minimum -- which is not what the check asks about: those rows
# are left out when it is read (MOST_IN_WINDOW). The report does not say where in the window an
# imperfection sat, so the count reads "held off only by the last five bases of a primer", bases a real
# polymerase occasionally reads through. Measured on 2.1.0 in validation/results/2026-11-27_v1.13.0_insilico210/.
LENIENT_PENALTY = 4
# The imperfections one primer's 3' window may carry for a lenient row to count as held by the 3' end.
# One is the question the check asks; two in the last five bases of one primer is no longer a base a
# polymerase occasionally reads through. On the Xylella panel under insilicoPCR 0.13.0, 818 of 4,189 held
# (assay, genome) pairs in PCR mode rested on two in one window; insilicoPCR 0.11.0 gave none.
MOST_IN_WINDOW = 1

QPCR_NAME = "assays_qpcr.fasta"
PCR_NAME = "assays_pcr.fasta"
SCRIPT_NAME = "run_insilico_pcr.sh"


@dataclass(frozen=True)
class Reaction:
    """What insilicoPCR judges a site's Tm against, and the reaction it works the Tm out for.

    A site is accepted when its Tm reaches `anneal_temp` less `tm_margin`. The reaction is the one Primer3
    designed the oligos in. insilicoPCR has one concentration for every oligo, so the probe is judged at the
    primers' (250 nM by default, against the 200 nM Primer3 designs it at: a few tenths of a degree
    hotter).
    """

    anneal_temp: float
    tm_margin: float
    monovalent_mm: float
    divalent_mm: float
    dntp_mm: float
    oligo_nm: float

    @property
    def needed_tm(self) -> float:
        """The Tm a site needs to be accepted."""
        return self.anneal_temp - self.tm_margin

    def options(self) -> list[str]:
        """The insilicoPCR options that say it. Always given, as `-l` is: insilicoPCR's defaults are its
        own to move, and the verdicts should not move with them. Written to every digit a float holds, so
        that what insilicoPCR is given is what primer-finder checked (59.9999999 is not rounded to 60)."""
        return ["--anneal-temp", _number(self.anneal_temp), "--tm-margin", _number(self.tm_margin),
                "--na", _number(self.monovalent_mm), "--mg", _number(self.divalent_mm),
                "--dntp", _number(self.dntp_mm), "--oligo-conc", _number(self.oligo_nm)]


def _number(value: float) -> str:
    """A number as a command line gives it: `60`, `0.8`, `59.9999999`."""
    return format(value, ".15g")


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
    # Exclusion genomes this assay amplifies only once an imperfection at the far 3' end is tolerated:
    # genomes held off by bases a real polymerase occasionally reads through. An annotation, not a
    # verdict -- the label and support do not move with it, but the ranking breaks ties away from it.
    exclusion_held_by_three_prime: int = 0
    # Of the inclusion genomes it amplifies, the Tm of the weakest site the product rests on in the genome
    # where that is lowest: how close the assay is to failing. None when it amplifies none.
    inclusion_weakest_tm: float | None = None
    # Why it gave no product in the inclusion genomes it misses, from insilicoPCR's report of the assays
    # without a product: which oligo, in how many genomes, and why (its best site's Tm, or the 3' rule).
    # Empty when none.
    inclusion_no_product: str = ""

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
    def label(self) -> str:
        """What the table says: `yes`, `partial (88%)` when it meets a threshold below 1, or `no`."""
        if self.complete:
            return "yes"
        if self.selective:
            return f"partial ({self.percent}%)"
        return "no"

    @property
    def support(self) -> int:
        """How far this check backs taking the assay to a bench, for the ranking: it amplified the whole
        inclusion group and nothing else (2), it met the run's threshold and nothing else (1), or it
        missed the group or amplified an exclusion genome (0)."""
        if self.complete:
            return 2
        return 1 if self.selective else 0

    def __str__(self) -> str:
        return (f"inclusion {self.inclusion_amplified}/{self.inclusion_total} ({self.percent}%), "
                f"exclusion {self.exclusion_amplified}/{self.exclusion_total}")


def write_primer_fasta(assays: Sequence[Assay], path: Path, kind: str) -> int:
    """Write the assays as insilicoPCR reads them, one file per question asked of it:

    - `qpcr`: the assays that have a probe, with it. A positive needs the probe to bind inside the product.
    - `pcr`: every assay, primers only. This is what says whether the primer pair discriminates on its own,
      which is what matters when a mismatch under a primer is all the specificity there is.

    insilicoPCR can run both kinds of assay from one file -- the probe requirement follows each assay --
    but one file cannot ask both questions about the same assay, so the probe assays are written twice and
    the two files are run separately.

    The forward primer's header declares the product size the assay was designed to give, `size=<bp>`:
    insilicoPCR's report then says of every product whether it is that one (`target`, within 10 bases of
    the size) or another (`unexpected`), and its HTML report counts a sample as positive only on the
    first. Nothing is filtered by it, there or here. The declaration sits on the forward primer alone --
    one is enough, and one place is one convention. No `ta=` is written: every assay is judged at the
    annealing temperature of the run, which is the one the design was made for.
    """
    if kind not in ("qpcr", "pcr"):  # pragma: no cover - the caller is this package
        raise PrimerFinderError(f"Unknown assay type {kind!r}")
    chosen = [assay for assay in assays if kind == "pcr" or assay.probe is not None]
    if not chosen:
        return 0
    records: list[seqio.Record] = []
    for assay in chosen:
        records.append(seqio.Record(f"{assay.name}-F", f"size={assay.product_size}", assay.forward.sequence))
        records.append(seqio.Record(f"{assay.name}-R", "", assay.reverse.sequence))
        if kind == "qpcr" and assay.probe is not None:
            records.append(seqio.Record(f"{assay.name}-P", "", assay.probe.sequence))
    seqio.write_fasta(path, records)
    return len(chosen)


STAGED = "insilico_genomes"


def stage_genomes(folder: Path, into: Path) -> int:
    """Lay a genome folder out the way insilicoPCR needs it: one flat folder of plain fasta files.

    It reads only the top level of what it is given, so a genome in a subfolder is silently skipped.
    primer-finder accepts subfolders, so the genomes are staged rather than passed through: subfolders are
    flattened, and a name used twice is numbered so that two files in different subfolders stay two
    samples. Gzipped genomes are written out plain -- insilicoPCR reads them as they are, but a staged
    sample is then one plain fasta whatever it came from, and the report's sample names stay the file
    names.

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


def _dotted(version: tuple[int, ...]) -> str:
    return ".".join(map(str, version))


def check_version(launcher: Sequence[str]) -> str:
    """Ask the launcher its version, and insist on REQUIRED_VERSION or newer. Returns the version line.

    Called before any of the expensive work, so that a release primer-finder cannot run stops the design
    before Primer3 rather than at its last step. A release before 2.0 accepts sites by mismatch counts,
    which the options this module passes no longer ask for: it would stop on `--anneal-temp`, and its
    reports name the assay column otherwise. `-v` fails when a release is older than 0.8.0, which does not
    know it, and when the jar is run with a Java older than 25, which insilicoPCR needs since 1.0.1 (a
    portable release brings its own).
    """
    needed = _dotted(REQUIRED_VERSION)
    try:
        output = tools.run([*launcher, "-v"])
    except tools.ToolError as error:
        raise PrimerFinderError(
            "insilicoPCR could not be run with -v to say its version: either the release is older than "
            "0.8.0, which does not know -v, or its jar was run with a Java older than 25, which it needs "
            "since 1.0.1 (an UnsupportedClassVersionError below says so: give --insilico-pcr the folder of "
            "an extracted portable release instead, which brings its own Java). primer-finder needs "
            f"{needed} or newer: {RELEASES}\n{error}"
        ) from error
    match = re.search(r"\b(\d+)\.(\d+)(?:\.(\d+))?\b", output)
    found = tuple(int(part or 0) for part in match.groups()) if match else None
    said = match.group(0) if match else f"no version in {output.strip()!r}"
    if found is None or found < TM_VERSION:
        raise PrimerFinderError(
            f"insilicoPCR {said} is older than the {needed} primer-finder needs. Before 2.0 it accepted a "
            "primer or probe site by a count of its mismatches (-m, --probe-mismatches), which "
            "primer-finder no longer passes: since 2.0 a site is accepted when its Tm on the template "
            f"reaches the annealing temperature less a margin. Upgrade it from {RELEASES}, or design with "
            "primer-finder 1.13.0, which runs insilicoPCR 0.12.0 to 0.14.0."
        )
    if found < REQUIRED_VERSION:
        raise PrimerFinderError(
            f"insilicoPCR {said} is older than the {needed} primer-finder needs. {needed} is the release "
            "its verdicts were checked against, and the first to count a mismatch among a primer's first "
            "five bases once in its report of the assays without a product, which says why an assay "
            f"misses an inclusion genome. Upgrade it from {RELEASES}."
        )
    return output.strip()


def _java() -> str:
    java = shutil.which("java")
    if java is None:
        raise PrimerFinderError(
            "insilicoPCR needs Java, which is not on PATH. Point --insilico-pcr at the folder of an "
            "extracted portable release instead: it brings its own Java runtime."
        )
    return java


def command(launcher: Sequence[str], genomes: Path, primers: Path, output: Path, threads: int,
            reaction: Reaction, three_prime_penalty: int | None = None) -> list[str]:
    invocation = [*launcher, "-i", str(genomes), "-p", str(primers), "-o", str(output),
                  "-t", str(threads), "-l", str(LONGEST_PRODUCT), *reaction.options()]
    if three_prime_penalty is not None:
        invocation += ["--3p-penalty", str(three_prime_penalty)]
    return invocation


def run(launcher: Sequence[str], genomes: Path, primers: Path, output: Path, threads: int,
        reaction: Reaction, three_prime_penalty: int | None = None) -> Path:
    """Run insilicoPCR over a folder of genomes. Returns the consolidated report it wrote."""
    output.mkdir(parents=True, exist_ok=True)
    report = output / REPORT
    if report.parent.is_dir():
        # insilicoPCR never overwrites a report: beside an existing report.tsv it writes report(1).tsv,
        # with its own report(1)_samples.tsv and the rest. The run's report is read where its output says,
        # but the folder keeps the latest run's at report.tsv, which the rerun script and anyone looking
        # rely on. The old report is set aside rather than deleted -- it may be hours of insilicoPCR a
        # failed rerun would otherwise leave nothing of -- the files beside it are written over by the
        # run, and only the names insilicoPCR itself writes are touched.
        try:
            for stale in report.parent.glob("report(*"):
                if NUMBERED_REPORT.fullmatch(stale.name):
                    stale.unlink()
            if report.exists():
                report.replace(report.with_name(report.name + ".previous"))
        except OSError as error:
            # Not PrimerFinderError's subclass of its own accord, and anything else would escape the
            # design command's "a failed check must not lose the design" net.
            raise PrimerFinderError(
                f"Cannot set aside the previous consolidated report in {report.parent}: {error}"
            ) from error
    try:
        said = tools.run(command(launcher, genomes, primers, output, threads, reaction, three_prime_penalty))
    except tools.ToolError as error:
        hint = ""
        if error.returncode == 2:
            # insilicoPCR's "the command line cannot be used as given". The version check has passed, so
            # it is an option's value: its own message, above, names it.
            hint = (" -- insilicoPCR refused the command line (exit 2): its message above says which "
                    "option or value it cannot use")
        raise PrimerFinderError(f"{error}{hint}") from error
    written = REPORT_WRITTEN.findall(said)
    if not written:
        raise PrimerFinderError(
            f"insilicoPCR did not say where it wrote its consolidated report (no 'Report written to' "
            f"line), so there is no telling which report in {report.parent} is this run's"
        )
    found = Path(written[-1])
    if not found.is_file():
        raise PrimerFinderError(f"insilicoPCR wrote no consolidated report ({found} missing)")
    return found


def beside(report: Path, kind: str) -> Path:
    """A file insilicoPCR writes next to a report and names after it: `report_samples.tsv` beside
    `report.tsv`, `report(1)_no_product.tsv` beside `report(1).tsv`."""
    return report.with_name(f"{report.stem}_{kind}.tsv")


def _needs_escape(text: str) -> bool:
    """Whether insilicoPCR writes this cell after an apostrophe (Reports.md, "File format"): it starts
    with `=` or `@`, or with `+` or `-` and is not a number (NUMBER), or with an apostrophe followed by
    such text."""
    if not text:
        return False
    if text[0] in "=@":
        return True
    if text[0] in "+-" and len(text) > 1:
        return not NUMBER.fullmatch(text)
    return text[0] == "'" and _needs_escape(text[1:])


def _cell(value: str | None) -> str:
    """A cell as it was before insilicoPCR wrote it for spreadsheets: a value that starts and ends with a
    double quote is unquoted, its doubled quotes made one, and an apostrophe put before formula-like text
    (a genome named `-ctrl.fasta` is the sample `'-ctrl`) is taken off."""
    text = (value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1].replace('""', '"')
    if text.startswith("'") and _needs_escape(text[1:]):
        text = text[1:]
    return text


def _rows(path: Path, is_report: bool = False) -> Iterable[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        # QUOTE_NONE: insilicoPCR writes its cells as they are, quoting none, so a genome whose header
        # holds a quote character would otherwise swallow the rest of its row and every row up to the next
        # one -- losing amplifications silently, which moves a verdict.
        rows = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        if is_report:
            _insist_on_columns(rows.fieldnames, path)
        for row in rows:
            yield {key: _cell(value) for key, value in row.items() if key is not None}


def _report_rows(report: Path) -> Iterable[tuple[str, str, dict[str, str]]]:
    """The report's rows that name both an assay and a sample, with the two."""
    for row in _rows(report, is_report=True):
        assay, sample = row.get(ASSAY_COLUMN, ""), row.get(SAMPLE_COLUMN, "")
        if assay and sample:
            yield assay, sample, row


def read_report(report: Path, designed_only: bool = False) -> dict[str, set[str]]:
    """Which samples each assay amplified, out of insilicoPCR's consolidated report.

    Every row is a called amplicon, so a row is an amplification and nothing here second-guesses it. With
    `designed_only`, only the products of the assay's forward primer with its reverse primer count: what
    the inclusion group is asked about. Without, every product of the tube does, one primer alone too:
    what the exclusion group must not give.
    """
    amplified: dict[str, set[str]] = {}
    for assay, sample, row in _report_rows(report):
        if designed_only and row.get(PAIRING_COLUMN) != DESIGNED:
            continue
        amplified.setdefault(assay, set()).add(sample)
    return amplified


def _insist_on_columns(fieldnames: Sequence[str] | None, report: Path) -> None:
    """The columns the verdicts are read from. Without them every assay would read as amplifying nothing
    at all, which looks like a clean answer rather than an unreadable report. A report of an older
    insilicoPCR is refused for what it is: its sites were accepted by mismatch counts, not by their Tm,
    and counting its rows would mix two rules without showing anywhere."""
    names = set(fieldnames or ())
    if OLD_ASSAY_COLUMN in names and ASSAY_COLUMN not in names:
        raise PrimerFinderError(
            f"{report} was written by an insilicoPCR older than 1.0 (its {OLD_ASSAY_COLUMN} column says "
            "so), which accepted primer and probe sites by mismatch counts. primer-finder needs "
            f"{_dotted(REQUIRED_VERSION)} or newer: {RELEASES}"
        )
    missing = [name for name in (SAMPLE_COLUMN, ASSAY_COLUMN) if name not in names]
    if missing:
        raise PrimerFinderError(
            f"insilicoPCR's report has no {' or '.join(missing)} column, so there is no telling what "
            f"amplified: {report}"
        )
    if not all(name in names for name in TM_COLUMNS):
        raise PrimerFinderError(
            f"{report} has no {', '.join(TM_COLUMNS)} columns, so it was written by an insilicoPCR older "
            "than 2.0, which accepted primer and probe sites by mismatch counts rather than by their Tm. "
            f"primer-finder needs {_dotted(REQUIRED_VERSION)} or newer: {RELEASES}"
        )
    missing = [name for name in REQUIRED_COLUMNS if name not in names]
    if missing:
        raise PrimerFinderError(
            f"{report} has no {', '.join(missing)} column(s), which every report of insilicoPCR "
            f"{_dotted(REQUIRED_VERSION)} has, so it cannot be read as one."
        )


def check_samples(report: Path, staged: int) -> set[str]:
    """Insist that insilicoPCR analysed every genome it was given. Returns the samples' names.

    insilicoPCR leaves out a sample whose file it cannot read, says so in a warning, and still exits 0.
    Such a sample has no row in the report, which here would read as "does not amplify": on the
    exclusion group a selective assay, on the inclusion group a miss. Neither is what the genome says, so
    the check fails instead -- the verdicts are never counted on a smaller group than the one asked about.
    """
    path = beside(report, "samples")
    if not path.is_file():
        raise PrimerFinderError(
            f"insilicoPCR wrote no list of the samples it analysed ({path} missing), so there is no "
            "telling a genome without a product from one it left out"
        )
    rows = list(_rows(path))
    left_out = [row for row in rows if row.get("Status") != "analysed"]
    if left_out:
        raise PrimerFinderError(
            f"insilicoPCR left {len(left_out)} genome(s) out of the run: "
            + "; ".join(f"{row.get(SAMPLE_COLUMN)} ({row.get('Reason') or row.get('Status')})"
                        for row in left_out)
            + f". Counted as not amplifying they would move the verdicts, so none are given; fix or "
            f"remove the genome(s) and design again. ({path})"
        )
    names = {row.get(SAMPLE_COLUMN, "") for row in rows}
    if len(names) != staged:
        raise PrimerFinderError(
            f"insilicoPCR analysed {len(names)} sample(s) but {staged} genome(s) were staged, so the "
            f"verdicts would be counted on another group than the one asked about ({path})"
        )
    return names


def _float(value: str | None) -> float | None:
    """A cell already read by `_cell` as a number, or None when it is empty or not one."""
    try:
        return float(value or "")
    except ValueError:
        return None


def weakest_tm(report: Path) -> dict[str, float]:
    """For each assay, the Tm of the weakest site its products rest on, in the sample where that is
    lowest. A product rests on its weakest oligo; a sample on its strongest designed product, which is
    the one the reaction gives first. Only designed products count, as for the inclusion verdict."""
    best: dict[tuple[str, str], float] = {}
    for assay, sample, row in _report_rows(report):
        if row.get(PAIRING_COLUMN) != DESIGNED:
            continue
        tms = [tm for tm in (_float(row.get(column)) for column in TM_COLUMNS) if tm is not None]
        if not tms:
            continue
        key = (assay, sample)
        best[key] = max(best.get(key, min(tms)), min(tms))
    weakest: dict[str, float] = {}
    for (assay, _), tm in best.items():
        weakest[assay] = min(weakest.get(assay, tm), tm)
    return weakest


ROLES = {"F": "forward", "R": "reverse", "P": "probe"}


def _role(oligo: str) -> str:
    """`forward`, `reverse` or `probe`, from an oligo ID (`x-R`, a degenerate variant `x-F_1`)."""
    _, _, suffix = oligo.rpartition("-")
    return ROLES.get(suffix[:1], oligo)


def _genomes(count: int) -> str:
    return f"{count} genome{'s' if count != 1 else ''}"


def no_product(report: Path) -> dict[str, str]:
    """Why each assay gave no product where it gave none, from insilicoPCR's report of the assays without
    a product. Per oligo without an accepted site, in how many samples, and for what reason, told apart by
    the cells insilicoPCR fills (Reports.md, "Assays without a product"):

    - its best site melts too low: `BestSiteTm` is filled, and the highest of them is given;
    - its best site fails the 3' rule: `BestSite` is filled and `BestSiteTm` empty;
    - it has no site at all: both are empty.

    e.g. `reverse: 3 genomes, best site 52.1 °C`, or `reverse: 3 genomes (best site 53.8 °C in 2; 3' rule
    in 1)` when the reasons differ. `no pair` is a sample where every oligo has a site but none face each
    other within the size limit. A sample with several rows for one oligo (variants of a degenerate one)
    is counted once, by the reason that comes first above."""
    path = beside(report, "no_product")
    if not path.is_file():
        return {}
    # assay -> oligo role -> sample -> (rank of its reason, Tm of its best site)
    found: dict[str, dict[str, dict[str, tuple[int, float | None]]]] = {}
    for row in _rows(path):
        assay, sample = row.get(ASSAY_COLUMN, ""), row.get(SAMPLE_COLUMN, "")
        if not assay or not sample:
            continue
        oligo = row.get("Oligo", "")
        what = _role(oligo) if oligo else "no pair"
        tm = _float(row.get("BestSiteTm"))
        reason = (0, tm) if tm is not None else ((1, None) if row.get("BestSite") else (2, None))
        samples = found.setdefault(assay, {}).setdefault(what, {})
        had = samples.get(sample)
        if had is None or reason[0] < had[0] or (reason[0] == had[0] == 0 and tm > had[1]):
            samples[sample] = reason
    order = ["forward", "reverse", "probe", "no pair"]
    summary: dict[str, str] = {}
    for assay, whats in found.items():
        parts = []
        for what in sorted(whats, key=lambda w: order.index(w) if w in order else len(order)):
            reasons = list(whats[what].values())
            text = f"{what}: {_genomes(len(reasons))}"
            if what != "no pair":
                tms = [tm for rank, tm in reasons if rank == 0]
                said = [(len(tms), f"best site {max(tms):.1f} °C" if tms else ""),
                        (sum(1 for rank, _ in reasons if rank == 1), "3' rule"),
                        (sum(1 for rank, _ in reasons if rank == 2), "no site")]
                said = [(count, words) for count, words in said if count]
                if len(said) == 1:
                    text += f", {said[0][1]}"
                else:
                    text += " (" + "; ".join(f"{words} in {count}" for count, words in said) + ")"
            parts.append(text)
        summary[assay] = "; ".join(parts)
    return summary


def read_three_prime_held(report: Path) -> dict[str, set[str]]:
    """Which samples each assay amplifies *only* once an imperfection at the far 3' end is tolerated.

    Read from the lenient exclusion rerun, which differs from the strict run in its `--3p-penalty` alone
    (LENIENT_PENALTY, which readmits a mismatch in the last five bases of a primer whose site still melts
    at the minimum). A sample counts for an assay when every amplicon row of the pair is one that rescue
    explains: a primer carries a 3'-window imperfection, and no primer carries more than MOST_IN_WINDOW.
    A row that binds clean is one the strict run gives too -- amplification for another reason -- and one
    such row means the sample is not held off by the 3' end alone. Nor is a row with two imperfections in
    one window: it amplifies only at a tolerance the check does not mean to extend.

    A row whose counts cannot be read is not claimed.
    """
    rescued: dict[tuple[str, str], bool] = {}
    for assay, sample, row in _report_rows(report):
        key = (assay, sample)
        rescued[key] = rescued.get(key, True) and _three_prime_explains(row)
    found: dict[str, set[str]] = {}
    for (assay, sample), only in rescued.items():
        if only:
            found.setdefault(assay, set()).add(sample)
    return found


def _three_prime_explains(row: dict[str, str]) -> bool:
    """Whether the lenient run's 3'-end tolerance is what let this row through."""
    windows = [_int(row.get(column)) for column in WINDOW_COLUMNS]
    if any(w is None for w in windows):
        return False
    return any(windows) and all(w <= MOST_IN_WINDOW for w in windows)


def _int(value: str | None) -> int | None:
    """A count column: blank is nothing reported and so zero, anything not a whole number is None."""
    text = (value or "").strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        return None


def count_genomes(folder: Path) -> int:
    return len(seqio.find_genomes(folder))


def verdicts(assays: Iterable[Assay], inclusion: Path, exclusion: Path,
             inclusion_report: dict[str, set[str]], exclusion_report: dict[str, set[str]],
             threshold: float = 1.0,
             exclusion_held: dict[str, set[str]] | None = None,
             inclusion_tm: dict[str, float] | None = None,
             inclusion_missing: dict[str, str] | None = None) -> dict[str, Verdict]:
    """Put the two reports together, one verdict per assay.

    `exclusion_held` is `read_three_prime_held` of the lenient exclusion rerun; a genome the strict run
    already counts as amplified is a real cross-reaction, not a held one, so it is left out of the count
    however the lenient rows read. `inclusion_tm` and `inclusion_missing` are `weakest_tm` and
    `no_product` of the inclusion run.
    """
    inclusion_total, exclusion_total = count_genomes(inclusion), count_genomes(exclusion)
    held, tm, missing = exclusion_held or {}, inclusion_tm or {}, inclusion_missing or {}
    return {
        assay.name: Verdict(
            inclusion_amplified=len(inclusion_report.get(assay.name, set())),
            inclusion_total=inclusion_total,
            exclusion_amplified=len(exclusion_report.get(assay.name, set())),
            exclusion_total=exclusion_total,
            threshold=threshold,
            exclusion_held_by_three_prime=len(
                held.get(assay.name, set()) - exclusion_report.get(assay.name, set())),
            inclusion_weakest_tm=tm.get(assay.name),
            inclusion_no_product=missing.get(assay.name, ""),
        )
        for assay in assays
    }


def write_script(path: Path, launcher: Sequence[str] | None, primers: Sequence[Path],
                 inclusion: Path, exclusion: Path, threads: int, reaction: Reaction,
                 three_prime_check: bool = True) -> Path:
    """A script that runs insilicoPCR over both groups, for running the check by hand later."""
    start = shlex.join(launcher) if launcher else (
        '"${INSILICO_PCR:?point INSILICO_PCR at the folder of an extracted insilicoPCR release}"'
        '/runtime/linux/jdk/bin/java -jar "${INSILICO_PCR}"/insilicoPCR.jar'
    )
    settings = f"-t {threads} -l {LONGEST_PRODUCT} {' '.join(reaction.options())}"
    # The 3'-end dependence rerun: the same run with the far 3' end readmitted, and nothing else changed.
    lenient = f"{settings} --3p-penalty {LENIENT_PENALTY}"
    needed = _dotted(REQUIRED_VERSION)
    lines = [
        "#!/usr/bin/env bash",
        f"# In silico PCR of the designed assays against both groups of genomes, with insilicoPCR {needed}+:",
        "#   https://github.com/duceppemo/insilicoPCR",
        "# An assay is selective when it amplifies every inclusion genome and no exclusion genome. A site is",
        f"# accepted when its Tm on the genome reaches {reaction.needed_tm:g} C (--anneal-temp less --tm-margin), "
        "in the",
        "# reaction the assays were designed for.",
    ]
    if three_prime_check:
        lines += [
            "# The *_exclusion_3prime runs redo the exclusion group tolerating an imperfection at the "
            "far 3' end:",
            "# a genome that amplifies only there is held off by bases a polymerase can read "
            "through.",
        ]
    lines += [
        "set -euo pipefail",
        "",
        f'out="${{1:-{path.parent}/insilico_pcr}}"',
        "# The files insilicoPCR writes for a numbered report, report(1).tsv and the ones named after it.",
        "numbered='^report\\([0-9]+\\)(\\.tsv|\\.html|_[a-z_]+\\.tsv)$'",
        "",
        *([] if launcher else [
            ': "${INSILICO_PCR:?point INSILICO_PCR at the folder of an extracted insilicoPCR release}"', ""]),
        # A design checks the release before it does anything; this script runs without one.
        "# A release older than 2.0 accepts sites by mismatch counts and refuses the options below; one",
        f"# older than {needed} is not the release these settings were checked against. Both are refused here.",
        f"version=\"$({start} -v 2> /dev/null | grep -Eo '[0-9]+\\.[0-9]+(\\.[0-9]+)?' | head -n 1 || true)\"",
        f'if [ -z "${{version}}" ] || '
        f'[ "$(printf \'%s\\n\' {needed} "${{version}}" | sort -V | head -n 1)" != "{needed}" ]; then',
        f'    echo "insilicoPCR {needed} or newer is needed; this is '
        '${version:-one that could not be run with -v (older than 0.8.0, or its jar run with a Java '
        'older than 25)}:" >&2',
        '    echo "  https://github.com/duceppemo/insilicoPCR/releases" >&2',
        "    exit 1",
        "fi",
        "",
    ]
    passes: list[tuple[str, Path, str]] = []
    for group, folder in (("inclusion", inclusion), ("exclusion", exclusion)):
        passes.append((group, folder, settings))
    if three_prime_check:
        passes.append(("exclusion_3prime", exclusion, lenient))
    for primer_file in primers:
        for group, folder, chosen in passes:
            lines += [
                # insilicoPCR writes report(1).tsv beside an existing report.tsv rather than over it, so
                # a rerun would leave the old verdicts where the find below looks. The old report is set
                # aside, not deleted, and of the numbered reports only the files insilicoPCR writes are
                # removed: nothing a user saved beside them is touched.
                f'dir="${{out}}/{primer_file.stem}_{group}/consolidated_report"',
                'mv -f "${dir}/report.tsv" "${dir}/report.tsv.previous" 2>/dev/null || true',
                'for stale in "${dir}"/report\\(*; do',
                '    if [[ "${stale##*/}" =~ ${numbered} ]]; then rm -f "${stale}"; fi',
                "done",
                f"{start} \\",
                f'    -i "{folder}" \\',
                f'    -p "{primer_file}" \\',
                f'    -o "${{out}}/{primer_file.stem}_{group}" \\',
                f"    {chosen}",
                "",
            ]
    lines += [
        'echo "Consolidated reports (beside each: report_samples.tsv, every genome analysed or not, and"',
        'echo "report_no_product.tsv, why an assay gave no product):"',
        'find "${out}" -name report.tsv',
        "",
    ]
    path.write_text("\n".join(lines))
    path.chmod(path.stat().st_mode | 0o111)
    return path
