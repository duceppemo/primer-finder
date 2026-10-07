"""Reading and writing fasta files, gzipped or not."""

from __future__ import annotations

import gzip
import os
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from primer_finder import PrimerFinderError

# The file extensions accepted in the inclusion and exclusion folders.
EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


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
    """Open a fasta file for reading, transparently decompressing a gzipped one."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("r")


def iter_records(path: Path | str) -> Iterator[Record]:
    """Yield the records of a fasta file, one at a time."""
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
    """Read a fasta file into a dictionary keyed by record name. The last entry wins on duplicate names."""
    return {record.name: record for record in iter_records(path)}


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
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
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
        raise PrimerFinderError("Not a fasta file (no header on the first line): " + ", ".join(bad))
    return genomes

