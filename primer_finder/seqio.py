"""Reading and writing fasta files, gzipped or not."""

from __future__ import annotations

import gzip
import logging
import os
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path

from primer_finder import PrimerFinderError

log = logging.getLogger(__name__)

# How many offending file names to put in an error message or a warning.
MOST_LISTED = 5

# The file extensions accepted in the inclusion and exclusion folders.
EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


# Every IUPAC code and both cases, so that a masked or ambiguous base survives the round trip.
COMPLEMENT = str.maketrans("ACGTNacgtnRYKMBDHVrykmbdhv", "TGCANtgcanYRMKVHDByrmkvhdb")


def reverse_complement(seq: str) -> str:
    """The reverse complement of a sequence, keeping the case of each base."""
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class Record:
    """One fasta entry: its name (the first word of the header), the rest of the header, and its sequence."""

    name: str
    desc: str
    seq: str

    @property
    def header(self) -> str:
        return f"{self.name} {self.desc}".strip()


def open_text(path: Path | str):
    """Open a fasta file for reading, transparently decompressing a gzipped one.

    A byte that is not UTF-8 is replaced instead of raising: fasta files from older pipelines hold accented
    strain names in their headers, and a header is only a label here.
    """
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", errors="replace")
    return path.open("r", errors="replace")


def iter_records(path: Path | str) -> Iterator[Record]:
    """Yield the records of a fasta file, one at a time.

    A gzipped file that stops in the middle raises `EOFError` from the decompressor, which says nothing
    about which file it was. The genome folders of a real run hold hundreds, so the name is the whole of
    the useful information.
    """
    try:
        yield from _records(path)
    except EOFError as error:
        raise PrimerFinderError(
            f"{path} is a gzipped file that ends before its data does ({error}). It is truncated: "
            "download or copy it again."
        ) from error


def _records(path: Path | str) -> Iterator[Record]:
    name, desc, seq = None, "", []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield Record(name, desc, "".join(seq))
                header = line[1:].split(None, 1)
                name = header[0] if header else ""
                desc = header[1] if len(header) > 1 else ""
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield Record(name, desc, "".join(seq))


def read_fasta(path: Path | str) -> dict[str, Record]:
    """Read a fasta file into a dictionary keyed by record name.

    Duplicate names are reported: the records of such a file cannot be told apart, and the last one would
    silently stand for all of them.
    """
    records: dict[str, Record] = {}
    duplicates: list[str] = []
    for record in iter_records(path):
        if record.name in records:
            duplicates.append(record.name)
        records[record.name] = record
    if duplicates:
        log.warning("%s holds %d record(s) whose name is not unique (%s); only the last of each is used",
                    path, len(duplicates), ", ".join(sorted(set(duplicates))[:MOST_LISTED]))
    return records


def write_fasta(path: Path | str, records: Iterable[Record]) -> int:
    """Write records to a fasta file, one line per sequence. Returns how many were written."""
    written = 0
    with Path(path).open("w") as fh:
        for record in records:
            fh.write(f">{record.header}\n{record.seq}\n")
            written += 1
    return written


def count_records(path: Path | str) -> int:
    """Count the records of a fasta file without holding it in memory."""
    count = 0
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1
    except EOFError as error:
        raise PrimerFinderError(
            f"{path} is a gzipped file that ends before its data does ({error}). It is truncated: "
            "download or copy it again."
        ) from error
    return count


def is_fasta(path: Path | str) -> bool:
    """True if the file starts with a fasta header (an empty or unreadable file is not a fasta)."""
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.strip():
                    return line.startswith(">")
    except (OSError, EOFError, gzip.BadGzipFile):
        return False
    return False


def base_name(path: Path) -> str:
    """The name of a genome file without its extension: "genome.fasta.gz" gives "genome"."""
    name = path.name
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    return Path(name).stem


def find_genomes(folder: Path | str, extensions: tuple[str, ...] = EXTENSIONS) -> list[Path]:
    """Every fasta file in a folder and its subfolders, sorted, following symbolic links."""
    found: list[Path] = []
    for root, _, filenames in os.walk(folder, followlinks=True):
        for filename in filenames:
            if filename.endswith(extensions):
                found.append(Path(root) / filename)
    return sorted(found)


def require_genomes(folder: Path, group: str) -> list[Path]:
    """The fasta files of an input folder, or an error naming what was expected."""
    if not folder.is_dir():
        raise PrimerFinderError(f"The {group} folder does not exist: {folder}")
    genomes = find_genomes(folder)
    if not genomes:
        raise PrimerFinderError(
            f"The {group} folder has no file with an accepted extension "
            f"({', '.join(EXTENSIONS)}): {folder}"
        )
    bad = [str(genome) for genome in genomes if not is_fasta(genome)]
    if bad:
        listed = ", ".join(sorted(bad)[:MOST_LISTED])
        extra = f", and {len(bad) - MOST_LISTED} more" if len(bad) > MOST_LISTED else ""
        raise PrimerFinderError(
            f"{len(bad)} file(s) in the {group} folder are not fasta (no header on the first line): "
            f"{listed}{extra}"
        )
    # A fasta whose extension is not one of the accepted ones used to be left out in silence, and a
    # genome missing from the exclusion group is the dangerous direction: its kmers are then reported as
    # specific to the inclusion group. Say so rather than quietly working from a smaller panel.
    skipped = sorted(str(path) for path in _other_fasta_files(folder, genomes))
    if skipped:
        listed = ", ".join(skipped[:MOST_LISTED])
        extra = f", and {len(skipped) - MOST_LISTED} more" if len(skipped) > MOST_LISTED else ""
        log.warning("%d file(s) in the %s folder look like fasta but do not end in %s, so they are NOT "
                    "part of the %s group: %s%s", len(skipped), group, " or ".join(EXTENSIONS), group,
                    listed, extra)
    # Two paths to one genome -- a symbolic link beside the file, or a link to the folder itself, which
    # `os.walk(followlinks=True)` descends -- would count that genome twice, and once per level of a
    # loop. One genome, once.
    seen: dict[Path, Path] = {}
    for genome in genomes:
        seen.setdefault(genome.resolve(), genome)
    if len(seen) < len(genomes):
        log.warning("%d of the %d file(s) in the %s folder are further paths to a genome already counted "
                    "and were left out", len(genomes) - len(seen), len(genomes), group)
    return sorted(seen.values())


def _other_fasta_files(folder: Path | str, found: Sequence[Path]) -> list[Path]:
    """Files in the folder that read as fasta but whose extension is not an accepted one."""
    accepted = {path.resolve() for path in found}
    others: list[Path] = []
    for root, _, filenames in os.walk(folder, followlinks=True):
        for filename in filenames:
            path = Path(root) / filename
            if path.resolve() in accepted or filename.startswith("."):
                continue
            try:
                if is_fasta(path):
                    others.append(path)
            except OSError:  # Unreadable, or a dangling link: not this function's business
                continue
    return others

