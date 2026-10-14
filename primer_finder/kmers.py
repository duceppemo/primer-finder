"""Counting kmers with KMC, and subtracting the exclusion kmers from the inclusion ones."""

from __future__ import annotations

import logging
from collections.abc import Callable, Iterable
from pathlib import Path

from primer_finder import PrimerFinderError, tools

log = logging.getLogger(__name__)

MIN_KMER_SIZE, MAX_KMER_SIZE = 1, 256  # What KMC accepts
# What KMC accepts for -m. Its own help says "from 1 to 1024", but 1 is refused at run time with
# "min memory must be at least 2GB", so the floor here is 2 -- measured, not read off the help.
MIN_MEMORY_GB, MAX_MEMORY_GB = 2, 1024
# KMC stores a count in one byte unless -cs says otherwise, so a kmer seen more than 255 times is reported
# as 255. It applies -ci/-cx to the true count before saturating, so a high -ci still matches; what the
# ceiling would cost is the count itself, which the max_count check below needs to be exact.
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


def dump(db_prefix: Path, dump_file: Path, threads: int,
         count_shared: Callable[[], int] | None = None) -> Path:
    """Write a KMC database as text: one kmer and its count per line.

    An empty result has two quite different causes and they want different advice, so `count_shared` --
    how many kmers the inclusion group had before the subtraction -- decides which is given. Blaming the
    exclusion group when the inclusion genomes never shared a kmer in the first place sends the reader
    after the wrong problem. It is called only when the result is empty: counting means dumping the whole
    inclusion database, which is not worth doing on a run that found something.
    """
    tools.run(["kmc_tools", f"-t{threads}", "transform", db_prefix, "dump", dump_file])
    if not dump_file.exists() or dump_file.stat().st_size == 0:
        if count_shared is not None and count_shared() == 0:
            raise PrimerFinderError(
                "No kmer at all is shared by the inclusion genomes, so there was nothing to subtract the "
                "exclusion group from. They may be too diverse for this kmer size: a smaller -k, or "
                "-p below 1 so that a kmer need not be in every one of them, would give it room."
            )
        raise PrimerFinderError(
            "No inclusion-specific kmer was found: every kmer shared by the inclusion genomes is also "
            "present in an exclusion genome. A smaller kmer size (-k) may help, or -p below 1 to let a "
            "kmer through that most of the inclusion genomes share."
        )
    return dump_file


def count_in_db(db_prefix: Path, dump_file: Path, threads: int) -> int:
    """How many kmers a KMC database holds, by dumping it and counting the lines.

    Only used to tell the two causes of an empty result apart, so it is worth the dump: the alternative is
    parsing KMC's own statistics, which it prints rather than returns.
    """
    tools.run(["kmc_tools", f"-t{threads}", "transform", db_prefix, "dump", dump_file])
    if not dump_file.exists():
        return 0
    with dump_file.open() as fh:
        return sum(1 for line in fh if line.strip())


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
