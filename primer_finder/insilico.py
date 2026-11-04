"""Handing the designed assays to insilicoPCR, and reading its verdict.

insilicoPCR (https://github.com/duceppemo/insilicoPCR) is a separate program: a portable Java application
that amplifies in silico against assemblies or reads with BBMap and BLAST+. primer-finder does not reimplement
it. It writes the primer fasta files insilicoPCR expects, and, when told where to find it, runs it over the
inclusion and the exclusion genomes and reads its consolidated report.

This module speaks insilicoPCR 0.8.0 or newer. Since that release a row of the consolidated report is a
called amplicon -- the two primers on the same contig, facing each other, within a size limit, and for a
probe assay with an accepted probe hit inside the product -- so counting rows is counting amplifications.
A report from an older release, whose rows were every pairing of two primer alignments and needed
re-checking here, is refused rather than counted; primer-finder 1.9.x carries that reading.

Two rules of insilicoPCR still shape what is written here:

- an assay is `<name>-F`, `<name>-R` and optionally `<name>-P`;
- an assay whose probe is in the primer file is only positive when the probe binds inside the product.
  Whether its primer pair discriminates *alone* is a different question, so the assays are also written
  with the probes stripped, and that file is run separately.
"""

from __future__ import annotations

import csv
import logging
import os
import re
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
# Columns only the reports of insilicoPCR 0.6.x and 0.7.x carry. Their rows were every pairing of a
# forward and a reverse primer alignment on a contig -- on one 3.3 Mb Brucella genome with 28 assays,
# 21,014 rows of which 71% were pairings PCR could not make -- and needed re-checking against the detailed
# reports here. 0.8.0 calls amplicons itself, dropped these columns, and made the re-checking obsolete, so
# finding one of them is finding a report this reader must not count rows from.
OLD_COLUMNS = ("ForwardEndMismatch", "ReverseEndMismatch")
# Columns only 0.8.0 and newer write. A report carrying neither kind is not known to be amplicons either,
# and counting it anyway would fail open -- phantom pairings read as cross-reactions -- so it is refused.
NEW_COLUMNS = ("Forward3PrimeMismatches", "Strand")
# The oldest insilicoPCR whose reports this module reads.
REQUIRED_VERSION = (0, 8, 0)
# The longest span that can be one product, passed to insilicoPCR as `-l`. It applies the limit itself
# since 0.8.0, but its default (5000) is its own to move, and the verdicts should not move with it.
# Well above anything this tool designs for -- Primer3 is asked for 70-150 bases by default.
LONGEST_PRODUCT = 10_000
# What the 3'-end dependence check hands insilicoPCR as --3p-penalty. Measured on 0.8.0
# (validation/results/2026-11-03_v1.10.0_insilico080/): 4 readmits the 3'-TERMINAL base, whatever it
# replaces -- an A/T pair scores 3 there, a G/C pair 4 -- and a mismatch five bases in (scored 1), but
# nothing 2-4 bases in, which blast trims as a tail whose stacked penalties pass 4. The report does not
# say where in the window the imperfection sat, so the count reads "held off only by the last five bases
# of a primer" -- usually the terminal base, the zone the ranking aims differences at, and in every
# admissible case a base a real polymerase occasionally reads through.
LENIENT_PENALTY = 4
# The per-primer mismatch count and the share of it inside the 3' window, per direction.
MISMATCH_COLUMNS = ("ForwardMismatches", "ReverseMismatches")
WINDOW_COLUMNS = ("Forward3PrimeMismatches", "Reverse3PrimeMismatches")

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
    # Exclusion genomes this assay amplifies only once an imperfection at the far 3' end is tolerated:
    # genomes held off by bases a real polymerase occasionally reads through. An annotation, not a
    # verdict -- the label and support do not move with it, but the ranking breaks ties away from it.
    exclusion_held_by_three_prime: int = 0

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

    insilicoPCR 0.8.0 can run both kinds of assay from one file -- the probe requirement follows each assay
    -- but one file cannot ask both questions about the same assay, so the probe assays are written twice
    and the two files are run separately.
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

    It reads only the top level of what it is given, so a genome in a subfolder is silently skipped.
    primer-finder accepts subfolders, so the genomes are staged rather than passed through: subfolders are
    flattened, and a name used twice is numbered so that two files in different subfolders stay two
    samples. Gzipped genomes are written out plain -- insilicoPCR 0.8.0 reads them as they are, but a
    staged sample is then one plain fasta whatever it came from, and the report's sample names stay the
    file names.

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


def check_version(launcher: Sequence[str]) -> str:
    """Ask the launcher its version, and insist on 0.8.0 or newer. Returns the version line.

    Called before any of the expensive work: an older release only announced itself at the in-silico step,
    after Primer3 and every blast had already run. 0.6.x does not know `-v` at all -- it prints an
    "Unrecognized option" usage dump and exits non-zero -- which is itself the answer.
    """
    needed = ".".join(map(str, REQUIRED_VERSION))
    try:
        output = tools.run([*launcher, "-v"])
    except tools.ToolError as error:
        raise PrimerFinderError(
            f"insilicoPCR could not say its version, which a release older than {needed} cannot "
            f"(it does not know -v). primer-finder needs {needed} or newer: "
            f"https://github.com/duceppemo/insilicoPCR/releases\n{error}"
        ) from error
    match = re.search(r"\b(\d+)\.(\d+)(?:\.(\d+))?\b", output)
    if match is None or tuple(int(part or 0) for part in match.groups()) < REQUIRED_VERSION:
        found = match.group(0) if match else f"no version in {output.strip()!r}"
        raise PrimerFinderError(
            f"insilicoPCR {found} is older than the {needed} primer-finder needs, whose report rows are "
            "amplicons rather than raw primer pairings. Upgrade it from "
            "https://github.com/duceppemo/insilicoPCR/releases, or design with primer-finder 1.9.x."
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


def command(launcher: Sequence[str], genomes: Path, primers: Path, output: Path,
            threads: int, mismatches: int, probe_mismatches: int | None = None,
            three_prime_penalty: int | None = None) -> list[str]:
    invocation = [*launcher, "-i", str(genomes), "-p", str(primers), "-o", str(output),
                  "-t", str(threads), "-m", str(mismatches), "-l", str(LONGEST_PRODUCT)]
    if probe_mismatches is not None:
        invocation += ["--probe-mismatches", str(probe_mismatches)]
    if three_prime_penalty is not None:
        invocation += ["--3p-penalty", str(three_prime_penalty)]
    return invocation


def run(launcher: Sequence[str], genomes: Path, primers: Path, output: Path,
        threads: int, mismatches: int = 0, probe_mismatches: int | None = None,
        three_prime_penalty: int | None = None) -> Path:
    """Run insilicoPCR over a folder of genomes. Returns the folder it wrote."""
    output.mkdir(parents=True, exist_ok=True)
    report = output / REPORT
    if report.parent.is_dir():
        # insilicoPCR 0.8.0 never overwrites a report: beside an existing report.tsv it writes
        # report(1).tsv, which the reader below would never look at -- so a rerun into the same folder
        # would hand this run the previous run's verdicts. The old report is set aside rather than
        # deleted -- it may be hours of insilicoPCR a failed rerun would otherwise leave nothing of --
        # and only the names insilicoPCR itself writes are touched, not anything a user saved beside them.
        try:
            if report.exists():
                report.replace(report.with_name(report.name + ".previous"))
            for stale in report.parent.glob("report(*).tsv"):
                stale.unlink()
        except OSError as error:
            # Not PrimerFinderError's subclass of its own accord, and anything else would escape the
            # design command's "a failed check must not lose the design" net.
            raise PrimerFinderError(
                f"Cannot set aside the previous consolidated report in {report.parent}: {error}"
            ) from error
    try:
        tools.run(command(launcher, genomes, primers, output, threads, mismatches, probe_mismatches,
                      three_prime_penalty))
    except tools.ToolError as error:
        hint = ""
        if "Unrecognized option" in str(error):
            # The one way an insilicoPCR older than 0.8.0 announces itself before writing anything: it
            # does not know the -l this run passes. Without the hint the user is left staring at a usage
            # message that never names the version.
            hint = (" -- the release looks older than 0.8.0, which primer-finder no longer speaks; "
                    "upgrade it from https://github.com/duceppemo/insilicoPCR/releases")
        raise PrimerFinderError(f"{error}{hint}") from error
    if not report.is_file():
        raise PrimerFinderError(f"insilicoPCR wrote no consolidated report ({report} missing)")
    return output


def read_report(output: Path) -> dict[str, set[str]]:
    """Which samples each assay amplified, out of insilicoPCR's consolidated report.

    Since insilicoPCR 0.8.0 every row is a called amplicon -- the two primers on the same contig, facing
    each other, no further apart than `-l`, each binding within the `-m` budget and with its last bases
    matching, and for a probe assay with the probe inside the product -- so a row is an amplification and
    nothing here second-guesses it. A report from an older release is refused: its rows were every pairing
    of two primer alignments on a contig, most of them nothing PCR could make.
    """
    report = output / REPORT
    amplified: dict[str, set[str]] = {}
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
            amplified.setdefault(assay, set()).add(sample)
    return amplified


def _insist_on_columns(fieldnames: Sequence[str] | None, report: Path) -> None:
    """The two columns every verdict is counted from. Without them every assay would read as amplifying
    nothing at all, which looks like a clean answer rather than an unreadable report. A report of an old
    insilicoPCR has both and is still refused: its rows are not amplicons, and counting them moves
    verdicts."""
    missing = [name for name in (SAMPLE_COLUMN, ASSAY_COLUMN) if name not in (fieldnames or ())]
    if missing:
        raise PrimerFinderError(
            f"insilicoPCR's report has no {' or '.join(missing)} column, so there is no telling what "
            f"amplified: {report}"
        )
    old = [name for name in OLD_COLUMNS if name in (fieldnames or ())]
    if old:
        raise PrimerFinderError(
            f"{report} was written by an insilicoPCR older than 0.8.0 (its {old[0]} column says so). "
            "Those releases reported every pairing of two primer alignments rather than amplicons, which "
            "this version of primer-finder no longer re-checks. Upgrade insilicoPCR "
            "(https://github.com/duceppemo/insilicoPCR/releases), or design with primer-finder 1.9.x."
        )
    if not any(name in (fieldnames or ()) for name in NEW_COLUMNS):
        # Failing open here would count rows that may be raw pairings; refusing is the honest answer.
        raise PrimerFinderError(
            f"{report} carries none of the columns insilicoPCR 0.8.0 writes "
            f"({' or '.join(NEW_COLUMNS)}), so there is no telling that its rows are amplicons rather "
            "than raw primer pairings. Upgrade insilicoPCR "
            "(https://github.com/duceppemo/insilicoPCR/releases)."
        )


def read_three_prime_held(output: Path, budget: int) -> dict[str, set[str]]:
    """Which samples each assay amplifies *only* once an imperfection at the far 3' end is tolerated.

    Read from the lenient exclusion rerun: `--3p-penalty` LENIENT_PENALTY, which readmits the terminal
    base or a mismatch five bases in and nothing between (measured, see LENIENT_PENALTY), and one more
    `-m` than the user's budget, which is what the readmitted base consumes. A sample counts for an assay when every amplicon row of the pair is
    one that rescue explains: a primer carries a 3'-window imperfection, and what each primer carries
    *outside* the window still fits the user's own budget. A row that binds clean at the lenient settings,
    or that needed the extra `-m` unit for an internal mismatch, is amplification for another reason, and
    one such row means the sample is not held off by the 3' end alone.

    The outside-the-window count is `ForwardMismatches - Forward3PrimeMismatches` (and the reverse pair),
    which counts an indel as one like everything else; a row whose counts cannot be read is not claimed.
    """
    report = output / REPORT
    rescued: dict[tuple[str, str], bool] = {}
    with report.open(newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        _insist_on_columns(rows.fieldnames, report)
        for row in rows:
            assay, sample = (row.get(ASSAY_COLUMN) or "").strip(), (row.get(SAMPLE_COLUMN) or "").strip()
            if not assay or not sample:
                continue
            key = (assay, sample)
            rescued[key] = rescued.get(key, True) and _three_prime_explains(row, budget)
    found: dict[str, set[str]] = {}
    for (assay, sample), only in rescued.items():
        if only:
            found.setdefault(assay, set()).add(sample)
    return found


def _three_prime_explains(row: dict[str, str], budget: int) -> bool:
    """Whether the lenient run's 3'-end tolerance is what let this row through."""
    counted = []
    for mismatches, window in zip(MISMATCH_COLUMNS, WINDOW_COLUMNS, strict=True):
        m, w = _int(row.get(mismatches)), _int(row.get(window))
        if m is None or w is None:
            return False
        counted.append((m, w))
    if not any(w for _, w in counted):
        return False
    return all(m - w <= budget for m, w in counted)


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
             exclusion_held: dict[str, set[str]] | None = None) -> dict[str, Verdict]:
    """Put the two reports together, one verdict per assay.

    `exclusion_held` is `read_three_prime_held` of the lenient exclusion rerun; a genome the strict run
    already counts as amplified is a real cross-reaction, not a held one, so it is left out of the count
    however the lenient rows read.
    """
    inclusion_total, exclusion_total = count_genomes(inclusion), count_genomes(exclusion)
    held = exclusion_held or {}
    return {
        assay.name: Verdict(
            inclusion_amplified=len(inclusion_report.get(assay.name, set())),
            inclusion_total=inclusion_total,
            exclusion_amplified=len(exclusion_report.get(assay.name, set())),
            exclusion_total=exclusion_total,
            threshold=threshold,
            exclusion_held_by_three_prime=len(
                held.get(assay.name, set()) - exclusion_report.get(assay.name, set())),
        )
        for assay in assays
    }


def write_script(path: Path, launcher: Sequence[str] | None, primers: Sequence[Path],
                 inclusion: Path, exclusion: Path, threads: int, mismatches: int,
                 probe_mismatches: int | None = None, three_prime_check: bool = True) -> Path:
    """A script that runs insilicoPCR over both groups, for running the check by hand later."""
    start = " ".join(launcher) if launcher else (
        '"${INSILICO_PCR:?point INSILICO_PCR at the folder of an extracted insilicoPCR release}"'
        '/runtime/linux/jdk/bin/java -jar "${INSILICO_PCR}"/insilicoPCR.jar'
    )
    settings = f"-t {threads} -m {mismatches} -l {LONGEST_PRODUCT}"
    if probe_mismatches is not None:
        settings += f" --probe-mismatches {probe_mismatches}"
    # The 3'-end dependence rerun: the far 3' end readmitted, the extra -m unit it consumes, and the
    # probes pinned at the strict budget so nothing else loosens with it.
    lenient = (f"-t {threads} -m {mismatches + 1} -l {LONGEST_PRODUCT} "
               f"--probe-mismatches {probe_mismatches if probe_mismatches is not None else mismatches} "
               f"--3p-penalty {LENIENT_PENALTY}")
    lines = [
        "#!/usr/bin/env bash",
        "# In silico PCR of the designed assays against both groups of genomes, with insilicoPCR 0.8.0+:",
        "#   https://github.com/duceppemo/insilicoPCR",
        "# An assay is selective when it amplifies every inclusion genome and no exclusion genome.",
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
                # aside, not deleted, and nothing a user saved beside it is touched.
                f'dir="${{out}}/{primer_file.stem}_{group}/consolidated_report"',
                'mv -f "${dir}/report.tsv" "${dir}/report.tsv.previous" 2>/dev/null || true',
                'rm -f "${dir}/report("*").tsv"',
                f"{start} \\",
                f'    -i "{folder}" \\',
                f'    -p "{primer_file}" \\',
                f'    -o "${{out}}/{primer_file.stem}_{group}" \\',
                f"    {chosen}",
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
