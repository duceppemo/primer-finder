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
# insilicoPCR reports, per primer, how many mismatches it bound through and how many bases blast had to
# trim off its 3' end to place it. The trim is the interesting one: those bases are not counted as
# mismatches at any tolerance, so an amplification that rests on one is an amplification insilicoPCR could
# not have refused. See docs/wiki/Designing-assays.md.
MISMATCH_COLUMNS = ("ForwardMismatches", "ReverseMismatches")
END_COLUMNS = ("ForwardEndMismatch", "ReverseEndMismatch")

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
        """That share as a whole number, rounded down so that it never overstates the coverage."""
        return int(self.fraction * 100)

    @property
    def complete(self) -> bool:
        """Amplifies every inclusion genome and no exclusion genome."""
        return self.inclusion_amplified == self.inclusion_total and self.exclusion_amplified == 0

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
        """What the table says: `yes`, `partial (88%)` when it meets a threshold below 1, or `no`."""
        if self.complete:
            return "yes"
        if self.selective:
            return f"partial ({self.percent}%)"
        return "no"

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


def read_report(output: Path) -> dict[str, set[str]]:
    """Which samples each assay amplified, out of insilicoPCR's consolidated report."""
    report = output / REPORT
    amplified: dict[str, set[str]] = {}
    with report.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            assay, sample = (row.get(ASSAY_COLUMN) or "").strip(), (row.get(SAMPLE_COLUMN) or "").strip()
            if assay and sample:
                amplified.setdefault(assay, set()).add(sample)
    return amplified


def _ends(row: dict[str, str]) -> bool:
    """Whether blast had to trim a base off the 3' end of either primer to place this amplicon."""
    return any(int((row.get(column) or "0").strip() or 0) != 0 for column in END_COLUMNS)


def read_terminal_only(output: Path) -> dict[str, set[str]]:
    """Which samples each assay amplified *only* through a difference in the last two bases of a primer.

    insilicoPCR does not count those bases: blast trims an unmatched base off the end of its alignment and
    the primer is called bound, at every `-m`. So an amplification whose every reported amplicon needed
    that trim is one the tool could not have refused, and a `no` verdict that rests entirely on them is the
    model's blind spot rather than a cross-reaction. A sample only counts when *all* of its amplicons needed
    it: one clean amplicon is a real one.
    """
    hits: dict[tuple[str, str], list[bool]] = {}
    with (output / REPORT).open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            assay, sample = (row.get(ASSAY_COLUMN) or "").strip(), (row.get(SAMPLE_COLUMN) or "").strip()
            if assay and sample:
                hits.setdefault((assay, sample), []).append(_ends(row))
    terminal: dict[str, set[str]] = {}
    for (assay, sample), ends in hits.items():
        if all(ends):
            terminal.setdefault(assay, set()).add(sample)
    return terminal


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
