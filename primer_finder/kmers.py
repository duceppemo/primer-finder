"""Counting kmers with KMC, and subtracting the exclusion kmers from the inclusion ones."""

from __future__ import annotations

import logging
from collections.abc import Iterable
from pathlib import Path

from primer_finder import PrimerFinderError, tools

log = logging.getLogger(__name__)

MIN_KMER_SIZE, MAX_KMER_SIZE = 1, 256  # What KMC accepts
MIN_MEMORY_GB, MAX_MEMORY_GB = 1, 1024  # What KMC accepts (-m)
# KMC counts up to 255 unless -cs says otherwise, so with more than 255 inclusion genomes -ci<N> would
# match nothing at all.
DEFAULT_COUNTER_MAX = 255
# KMC's largest accepted kmer count: enough to mean "no upper limit" for the exclusion group.
NO_MAX_COUNT = 1_000_000_000


def kmc_memory(memory_gb: int) -> int:
    """The memory to give KMC (-m), within what it accepts."""
    return max(MIN_MEMORY_GB, min(MAX_MEMORY_GB, memory_gb))


def write_file_list(paths: Iterable[Path], list_file: Path) -> Path:
    """Write the input file paths one per line, the file list KMC reads with @."""
    list_file.write_text("".join(f"{path}\n" for path in paths))
    return list_file


def count(
    list_file: Path,
    db_prefix: Path,
    work_dir: Path,
    kmer_size: int,
    threads: int,
    memory_gb: int,
    min_count: int,
    max_count: int,
) -> Path:
    """Count the kmers of the fasta files listed in `list_file`, keeping those seen `min_count`..`max_count`
    times. Returns the prefix of the KMC database (its `.kmc_pre` and `.kmc_suf` files)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "kmc",
        f"-k{kmer_size}",
        f"-t{threads}",
        f"-m{kmc_memory(memory_gb)}",
        "-fm",  # multi-fasta input
    ]
    if DEFAULT_COUNTER_MAX < max_count < NO_MAX_COUNT:
        # Count beyond 255, or no kmer would reach min_count and repeated kmers would slip past max_count.
        # When nothing is filtered out by max_count (the exclusion group) the default ceiling is harmless,
        # and leaving it alone keeps KMC's counters one byte wide.
        command.append(f"-cs{max_count}")
    command += [f"-ci{min_count}", f"-cx{max_count}", f"@{list_file}", str(db_prefix), str(work_dir)]
    tools.run(command)
    return require_db(db_prefix, "KMC counted no kmers")


def subtract(minuend: Path, subtrahend: Path, db_prefix: Path, threads: int) -> Path:
    """Keep the kmers of `minuend` that are absent from `subtrahend` (KMC's kmers_subtract)."""
    tools.run([
        "kmc_tools", f"-t{threads}", "simple",
        minuend, subtrahend,
        "kmers_subtract", db_prefix,
    ])
    return require_db(db_prefix, "KMC could not subtract the exclusion kmers")


def dump(db_prefix: Path, dump_file: Path, threads: int) -> Path:
    """Write a KMC database as text: one kmer and its count per line."""
    tools.run(["kmc_tools", f"-t{threads}", "transform", db_prefix, "dump", dump_file])
    if not dump_file.exists() or dump_file.stat().st_size == 0:
        raise PrimerFinderError(
            "No inclusion-specific kmer was found: every kmer shared by the inclusion genomes is also "
            "present in an exclusion genome. A smaller kmer size (-k) or a higher -d may help."
        )
    return dump_file


def require_db(db_prefix: Path, message: str) -> Path:
    """Check that KMC wrote both files of a database."""
    missing = [db_prefix.name + suffix for suffix in (".kmc_pre", ".kmc_suf")
               if not (db_prefix.parent / (db_prefix.name + suffix)).exists()]
    if missing:
        raise PrimerFinderError(f"{message} ({', '.join(missing)} missing)")
    return db_prefix


def dump_to_fasta(dump_file: Path, fasta_file: Path) -> int:
    """Convert a KMC dump to a fasta file of numbered kmers. Returns how many were written."""
    count = 0
    with dump_file.open() as src, fasta_file.open("w") as dst:
        for line in src:
            line = line.strip()
            if not line:
                continue
            dst.write(f">kmer_{count}\n{line.split()[0]}\n")
            count += 1
    if count == 0:  # pragma: no cover - dump() already refuses an empty dump
        raise PrimerFinderError("The KMC dump holds no kmer")
    return count
