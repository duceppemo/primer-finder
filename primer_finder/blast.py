"""Checking the candidate contigs against every inclusion and exclusion genome with blast."""

from __future__ import annotations

import gzip
import logging
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
# A variant position is kept when it differs in at least this fraction of the exclusion genomes that the
# contig hits: a few exclusion genomes may carry the inclusion allele without making the assay useless.
MIN_EXCLUSION_FRACTION = 0.90

PRESENCE_FIELDS = ("qseqid", "evalue")
ALIGNMENT_FIELDS = ("qseqid", "qstart", "qend", "evalue", "qseq", "sseq")


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
class ExclusionResult:
    """What the exclusion genomes say about one contig."""

    genomes_hit: int = 0
    positions: list[list[int]] = field(default_factory=list)


def make_db(genome: Path, work_dir: Path) -> Path:
    """Make a blast database for one genome, inside `work_dir` so that the input folders are left alone.
    A gzipped genome is decompressed first: makeblastdb does not read gzip."""
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = genome
    if genome.name.endswith(".gz"):
        fasta = work_dir / genome.name[: -len(".gz")]
        with gzip.open(genome, "rb") as src, fasta.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    db = work_dir / base_name(genome)
    tools.run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", db])
    return db


def blastn(db: Path, query: Path, out_file: Path, fields: Sequence[str], max_targets: int = 10) -> Path:
    """Run blastn, writing the chosen tabular fields to `out_file`."""
    tools.run([
        "blastn",
        "-db", db,
        "-query", query,
        "-out", out_file,
        "-evalue", str(MAX_EVALUE),
        "-max_target_seqs", str(max_targets),
        "-num_threads", "1",
        "-outfmt", "6 " + " ".join(fields),
    ])
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
            if position < hit.qend:
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


def presence_in_genomes(
    query: Path, genomes: Sequence[Path], work_dir: Path, threads: int
) -> dict[str, dict[str, bool]]:
    """For every contig of `query`, whether it is present in each genome. Keyed by contig, then by genome
    file name."""
    def one(genome: Path) -> set[str]:
        db = make_db(genome, work_dir / base_name(genome))
        out_file = work_dir / f"{base_name(genome)}.tsv"
        blastn(db, query, out_file, PRESENCE_FIELDS, max_targets=1)
        return {hit.query for hit in parse_hits(out_file, PRESENCE_FIELDS) if hit.evalue <= MAX_EVALUE}

    found = _parallel(((genome,) for genome in genomes), one, threads)
    contigs = _query_names(query)
    presence: dict[str, dict[str, bool]] = {contig: {} for contig in contigs}
    for genome, hits in zip(genomes, found, strict=True):
        for contig in contigs:
            presence[contig][genome.name] = contig in hits
    return presence


def _query_names(query: Path) -> list[str]:
    """The names of the contigs in a fasta file, in file order."""
    return [record.name for record in iter_records(query)]


def exclusion_variants(
    query: Path, genomes: Sequence[Path], work_dir: Path, threads: int
) -> dict[str, ExclusionResult]:
    """For every contig of `query`, the variant positions found in each exclusion genome it hits."""
    def one(genome: Path) -> list[Hit]:
        db = make_db(genome, work_dir / base_name(genome))
        out_file = work_dir / f"{base_name(genome)}.tsv"
        blastn(db, query, out_file, ALIGNMENT_FIELDS)
        return parse_hits(out_file, ALIGNMENT_FIELDS)

    results: dict[str, ExclusionResult] = {name: ExclusionResult() for name in _query_names(query)}
    for hits in _parallel(((genome,) for genome in genomes), one, threads):
        best: dict[str, Hit] = {}
        for hit in hits:
            if hit.evalue > MAX_EVALUE:
                continue
            kept = best.get(hit.query)
            # The best alignment of this contig in this genome: the most significant, then the longest.
            if kept is None or (hit.evalue, -(hit.qend - hit.qstart)) < (kept.evalue, -(kept.qend - kept.qstart)):
                best[hit.query] = hit
        for contig, hit in best.items():
            result = results.setdefault(contig, ExclusionResult())
            result.genomes_hit += 1
            result.positions.append(variant_positions(hit))
    return results


def shared_variants(result: ExclusionResult, fraction: float = MIN_EXCLUSION_FRACTION) -> list[int]:
    """The contig positions that differ from at least `fraction` of the exclusion genomes it hits."""
    if not result.positions:
        return []
    counts: dict[int, int] = {}
    for positions in result.positions:
        for position in positions:
            counts[position] = counts.get(position, 0) + 1
    needed = fraction * len(result.positions)
    return sorted(position for position, count in counts.items() if count >= needed)
