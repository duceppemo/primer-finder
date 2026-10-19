"""Reading and writing fasta files, gzipped or not."""

from __future__ import annotations

import gzip
import logging
import os
import stat
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
    if path.suffix.lower() == ".gz":
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


# The first bytes of the compressed formats a genome arrives in. A file in one of these is a genome by
# any reasonable reading, and it used to leave the group in silence: `is_fasta` opened it as text, found
# no `>`, and the "looks like fasta, wrong extension" warning therefore never named it.
COMPRESSION_MAGIC = {
    b"\x1f\x8b": "gzip", b"BZh": "bzip2", b"\xfd7zXZ\x00": "xz", b"PK\x03\x04": "zip",
    b"\x04\x22M\x18": "lz4", b"\x28\xb5/\xfd": "zstandard",
}


def _is_ordinary(path: Path | str) -> bool:
    """Whether a path is a regular file, following links. Anything else cannot hold a genome."""
    try:
        return stat.S_ISREG(os.stat(path).st_mode)
    except OSError:
        return False


def compressed_as(path: Path | str) -> str | None:
    """The compression format a file's first bytes say it is in, if any."""
    try:
        if not stat.S_ISREG(os.stat(path).st_mode):
            return None
        with open(path, "rb") as fh:
            head = fh.read(8)
    except OSError:
        return None
    return next((name for magic, name in COMPRESSION_MAGIC.items() if head.startswith(magic)), None)


def is_fasta(path: Path | str) -> bool:
    """True if the file starts with a fasta header (an empty or unreadable file is not a fasta)."""
    try:
        # Only an ordinary file is opened. Opening a named pipe with nothing on the other end never
        # returns, and this is called on every path with an accepted extension -- a FIFO called
        # `pipe.fasta` in a genome folder used to hang the run at discovery with no output at all.
        if not stat.S_ISREG(os.stat(path).st_mode):
            return False
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
    """Every fasta file in a folder and its subfolders, sorted, following symbolic links, each genome once.

    De-duplicated here rather than in `require_genomes` alone, because the staging and the counting that
    the in-silico verdicts are measured against call this directly: a link beside a genome used to make
    the denominator of those verdicts one larger than the panel the kmers were counted from.
    """
    return unique_paths(walk_genomes(folder, extensions))


def walk_genomes(folder: Path | str, extensions: tuple[str, ...] = EXTENSIONS) -> list[Path]:
    """Every path in a folder and its subfolders that ends in one of `extensions`, sorted, duplicates and
    all. `require_genomes` wants the duplicates so that it can say it dropped them."""
    found: list[Path] = []
    for root, _, filenames in os.walk(folder, followlinks=True):
        for filename in filenames:
            if filename.endswith(extensions):
                found.append(Path(root) / filename)
    return sorted(found)


def unique_paths(paths: Iterable[Path]) -> list[Path]:
    """One path per file, keeping the first name each was found under. Two paths to one genome -- a link
    beside the file, or a link to the folder itself, which the walk descends -- are one genome."""
    seen: dict[Path, Path] = {}
    for path in paths:
        try:
            key = path.resolve()
        except (OSError, RuntimeError):  # A symlink loop, or a name that cannot be resolved
            key = path.absolute()
        seen.setdefault(key, path)
    return sorted(seen.values())


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
    # A path that is not an ordinary file is refused on its own terms: "no header on the first line" is
    # not true of a named pipe, and a pipe is the one that used to hang the run rather than fail it.
    odd = [str(genome) for genome in genomes if not _is_ordinary(genome)]
    if odd:
        listed = ", ".join(sorted(odd)[:MOST_LISTED])
        extra = f", and {len(odd) - MOST_LISTED} more" if len(odd) > MOST_LISTED else ""
        raise PrimerFinderError(
            f"{len(odd)} path(s) in the {group} folder are not ordinary files, so they cannot be read as "
            f"a genome (a named pipe, a socket, a device or a broken link): {listed}{extra}"
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
        log.warning("%d file(s) in the %s folder look like a genome -- fasta, or compressed data -- but "
                    "do not end in %s, so they are NOT part of the %s group: %s%s", len(skipped), group,
                    " or ".join(EXTENSIONS), group, listed, extra)
    # `find_genomes` already dropped the second path to a genome; this only says so.
    duplicates = len(walk_genomes(folder)) - len(genomes)
    if duplicates > 0:
        log.warning("%d file(s) in the %s folder are further paths to a genome already counted -- a link "
                    "beside it, or a link to the folder itself -- and were left out", duplicates, group)
    return genomes


def _other_fasta_files(folder: Path | str, found: Sequence[Path]) -> list[Path]:
    """Files in the folder that read as fasta but whose extension is not an accepted one."""
    accepted = set()
    for path in found:
        try:
            accepted.add(path.resolve())
        except (OSError, RuntimeError):
            accepted.add(path.absolute())
    others: list[Path] = []
    for root, _, filenames in os.walk(folder, followlinks=True):
        for filename in filenames:
            path = Path(root) / filename
            if filename.startswith("."):
                continue
            try:
                if path.resolve() in accepted:
                    continue
                # Only ordinary files are opened below. Opening a named pipe with nothing on the other
                # end never returns, and this walks whatever happens to be in the folder.
                if not stat.S_ISREG(os.stat(path).st_mode):
                    continue
                if is_fasta(path) or compressed_as(path):
                    others.append(path)
            except (OSError, RuntimeError):  # Unreadable, a dangling link, a symlink loop
                continue
    return others
