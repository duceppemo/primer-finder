"""Checking the candidate contigs against every inclusion and exclusion genome with blast."""

from __future__ import annotations

import gzip
import logging
import os
import shutil
from collections.abc import Iterable, Sequence
from concurrent import futures
from dataclasses import dataclass, field
from pathlib import Path

from primer_finder import tools
from primer_finder.mapping import PRIMER_LENGTH
from primer_finder.seqio import base_name, iter_records

log = logging.getLogger(__name__)

# A hit has to be at least this good to count as "the contig is there".
MAX_EVALUE = 1e-10
# A variant position is kept when it differs in at least this fraction of the exclusion genomes that hold
# the region: a few exclusion genomes may carry the inclusion allele without making the assay useless.
MIN_EXCLUSION_FRACTION = 0.90

PRESENCE_FIELDS = ("qseqid", "evalue")
ALIGNMENT_FIELDS = ("qseqid", "qstart", "qend", "evalue", "qseq", "sseq")

# BLAST splits its path arguments on whitespace, so it is never given a path that could hold a space: each
# genome is linked into its own folder under these fixed names, and blast runs in that folder.
LOCAL_FASTA = "genome.fasta"
DB_NAME = "db"
HITS_NAME = "hits.tsv"


@dataclass
class Hit:
    """One blast high-scoring pair, with the aligned query and subject sequences."""

    query: str
    qstart: int
    qend: int
    evalue: float
    qseq: str = ""
    sseq: str = ""


@dataclass
class Alignments:
    """What one exclusion genome holds of one contig: every part of the contig it aligns, and where each
    of those alignments differs from it."""

    spans: list[tuple[int, int, set[int]]] = field(default_factory=list)  # start, end (0-based), variants

    def add(self, hit: Hit) -> None:
        self.spans.append((hit.qstart - 1, hit.qend, set(variant_positions(hit))))

    def covers(self, position: int) -> bool:
        return any(start <= position < end for start, end, _ in self.spans)

    def differs_at(self, position: int) -> bool:
        """True when every copy of this region in the genome differs from the contig at this position. One
        matching copy is enough for an assay to amplify the genome, so one is enough to say no."""
        covering = [(start, end, variants) for start, end, variants in self.spans if start <= position < end]
        return bool(covering) and all(position in variants for _, _, variants in covering)

    @property
    def variants(self) -> set[int]:
        return {position for _, _, variants in self.spans for position in variants}


@dataclass
class ExclusionResult:
    """What the exclusion genomes say about one contig: one entry per genome that holds any part of it."""

    genomes: list[Alignments] = field(default_factory=list)

    @property
    def genomes_hit(self) -> int:
        return len(self.genomes)


def make_db(genome: Path, work_dir: Path) -> Path:
    """Make a blast database for one genome, inside `work_dir` so that the input folders are left alone.

    The genome is linked into `work_dir` under a fixed name (a gzipped one is decompressed: makeblastdb does
    not read gzip), and makeblastdb runs there with relative paths, so neither a space in the path nor two
    genomes with the same file name can reach BLAST.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    local = work_dir / LOCAL_FASTA
    if local.exists() or local.is_symlink():
        local.unlink()
    if genome.name.endswith(".gz"):
        with gzip.open(genome, "rb") as src, local.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        try:
            local.symlink_to(genome.resolve())
        except OSError:  # pragma: no cover - a filesystem without symbolic links
            shutil.copyfile(genome, local)
    tools.run(["makeblastdb", "-in", LOCAL_FASTA, "-dbtype", "nucl", "-out", DB_NAME], cwd=work_dir)
    return work_dir / DB_NAME


def blastn(db: Path, query: Path, out_file: Path, fields: Sequence[str], max_targets: int = 10) -> Path:
    """Run blastn in the database's folder, writing the chosen tabular fields to `out_file`."""
    work_dir = db.parent
    tools.run([
        "blastn",
        "-db", db.name,
        "-query", os.path.relpath(query, work_dir),
        "-out", os.path.relpath(out_file, work_dir),
        "-evalue", str(MAX_EVALUE),
        "-max_target_seqs", str(max_targets),
        "-num_threads", "1",
        "-outfmt", "6 " + " ".join(fields),
    ], cwd=work_dir)
    return out_file


def parse_hits(out_file: Path, fields: Sequence[str]) -> list[Hit]:
    """Read a tabular blast output written with `fields`."""
    hits: list[Hit] = []
    index = {name: position for position, name in enumerate(fields)}
    with out_file.open() as fh:
        for line in fh:
            values = line.rstrip("\n").split("\t")
            if len(values) != len(fields):
                continue
            hits.append(Hit(
                query=values[index["qseqid"]],
                qstart=int(values[index["qstart"]]) if "qstart" in index else 0,
                qend=int(values[index["qend"]]) if "qend" in index else 0,
                evalue=float(values[index["evalue"]]),
                qseq=values[index["qseq"]] if "qseq" in index else "",
                sseq=values[index["sseq"]] if "sseq" in index else "",
            ))
    return hits


def variant_positions(hit: Hit) -> list[int]:
    """The positions of the contig (0-based) that differ from this exclusion sequence: mismatches, inserted
    bases, and the base next to a deletion."""
    positions: list[int] = []
    position = hit.qstart - 1
    for query_base, subject_base in zip(hit.qseq, hit.sseq, strict=False):
        if query_base == "-":  # The exclusion genome has bases the contig does not: mark where they fit in
            positions.append(position)
            continue
        if subject_base == "-" or query_base.upper() != subject_base.upper():
            positions.append(position)
        position += 1
    return sorted(set(positions))


def has_close_variants(positions: Sequence[int], window: int = PRIMER_LENGTH) -> bool:
    """True if at least two of these positions are less than `window` bases apart."""
    return any(second - first < window for first, second in zip(positions, positions[1:], strict=False))


def _parallel(work: Iterable[tuple], function, threads: int) -> list:
    """Run `function(*arguments)` for every item, at most `threads` at a time, keeping the input order."""
    work = list(work)
    if not work:
        return []
    with futures.ThreadPoolExecutor(max_workers=max(1, min(threads, len(work)))) as executor:
        return list(executor.map(lambda arguments: function(*arguments), work))


def genome_folder(work_dir: Path, index: int, genome: Path) -> Path:
    """A folder of its own for each genome: two genomes in different subfolders can have the same file
    name, so the position in the list is what makes it unique."""
    return work_dir / f"{index:04d}_{base_name(genome)}"


def presence_in_genomes(
    query: Path, genomes: Sequence[Path], work_dir: Path, threads: int
) -> dict[str, dict[str, bool]]:
    """For every contig of `query`, whether it is present in each genome. Keyed by contig, then by the
    genome's path, which is what tells two genomes with the same file name apart."""
    def one(index: int, genome: Path) -> set[str]:
        folder = genome_folder(work_dir, index, genome)
        db = make_db(genome, folder)
        hits = parse_hits(blastn(db, query, folder / HITS_NAME, PRESENCE_FIELDS, max_targets=1),
                          PRESENCE_FIELDS)
        return {hit.query for hit in hits if hit.evalue <= MAX_EVALUE}

    found = _parallel(enumerate(genomes), one, threads)
    contigs = _query_names(query)
    presence: dict[str, dict[str, bool]] = {contig: {} for contig in contigs}
    for genome, hits in zip(genomes, found, strict=True):
        for contig in contigs:
            presence[contig][str(genome)] = contig in hits
    return presence


def _query_names(query: Path) -> list[str]:
    """The names of the contigs in a fasta file, in file order."""
    return [record.name for record in iter_records(query)]


def exclusion_variants(
    query: Path, genomes: Sequence[Path], work_dir: Path, threads: int
) -> dict[str, ExclusionResult]:
    """For every contig of `query`, what each exclusion genome that holds part of it looks like there."""
    def one(index: int, genome: Path) -> list[Hit]:
        folder = genome_folder(work_dir, index, genome)
        db = make_db(genome, folder)
        return parse_hits(blastn(db, query, folder / HITS_NAME, ALIGNMENT_FIELDS), ALIGNMENT_FIELDS)

    results: dict[str, ExclusionResult] = {name: ExclusionResult() for name in _query_names(query)}
    for hits in _parallel(enumerate(genomes), one, threads):
        per_contig: dict[str, Alignments] = {}
        for hit in hits:
            if hit.evalue > MAX_EVALUE:
                continue
            per_contig.setdefault(hit.query, Alignments()).add(hit)
        for contig, alignments in per_contig.items():
            results.setdefault(contig, ExclusionResult()).genomes.append(alignments)
    return results


def shared_variants(result: ExclusionResult, fraction: float = MIN_EXCLUSION_FRACTION) -> list[int]:
    """The contig positions worth designing a primer on: those that differ from at least `fraction` of the
    exclusion genomes that hold that part of the contig.

    A genome that does not align a position says nothing about it — it does not hold that region, which only
    makes the assay more selective — so it is left out of the count rather than counted as identical.
    """
    kept: list[int] = []
    candidates = {position for alignments in result.genomes for position in alignments.variants}
    for position in sorted(candidates):
        covering = [alignments for alignments in result.genomes if alignments.covers(position)]
        differing = [alignments for alignments in covering if alignments.differs_at(position)]
        if covering and len(differing) >= fraction * len(covering):
            kept.append(position)
    return kept
