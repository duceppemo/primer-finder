"""Mapping the assembled kmers to one exclusion genome, and keeping the contigs whose differences could
carry a selective assay."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from primer_finder import tools
from primer_finder.seqio import Record, reverse_complement

log = logging.getLogger(__name__)

# The widest primer Primer3 may return (`design.MAX_PRIMER_SIZE`), and so the widest stretch two
# differences can share and still both sit under one primer. The design step aims its pinned requests at
# windows of this width, and a contig is worth keeping on the same terms: being stricter here throws away
# pairs that a primer could carry. Counted as a span, so two differences qualify when the second is fewer
# than this many bases past the first.
MAX_PRIMER_SPAN = 25

# How much of a contig end the aligner must fail to place for the end alone to be worth an assay. A
# different quantity from the span above -- it is how much unplaceable sequence is enough, not how far
# apart two differences may be -- so it keeps the value it has always had rather than following it.
MIN_CLIPPED = 21

# The cigar operations that are a difference with the reference, those that are a match, those that are
# bases of the contig the aligner left out, and those that consume query bases.
DIFFERENCE_OPS = "XID"
MATCH_OPS = "=M"
CLIP_OPS = "SH"
QUERY_OPS = "MIS=X"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800



@dataclass
class Alignment:
    """One primary SAM record: the query name, whether it mapped, its strand and its cigar."""

    name: str
    mapped: bool
    reverse: bool
    cigar: list[tuple[int, str]]


@dataclass
class Candidate:
    """A contig kept for the blast steps: its name, what makes it different (the cigar string, or the
    length of the contig when it does not map at all), and its sequence with the differences in lowercase."""

    name: str
    desc: str
    seq: str

    def as_record(self) -> Record:
        return Record(self.name, self.desc, self.seq)


def parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """A cigar string as a list of (length, operation) pairs. "*" (no cigar) gives an empty list."""
    if cigar == "*":
        return []
    return [(int(length), op) for length, op in CIGAR.findall(cigar)]


def cigar_string(cigar: Iterable[tuple[int, str]]) -> str:
    return "".join(f"{length}{op}" for length, op in cigar)


def map_contigs(assembly: Path, reference: Path, sam_file: Path, threads: int) -> Path:
    """Map every contig to the reference exclusion genome, writing the alignments to a SAM file.
    `--eqx` is what makes minimap2 write = and X instead of M, so that mismatches can be told apart."""
    tools.run(
        ["minimap2", "-t", str(threads), "-a", "--eqx", reference, assembly],
        stdout_path=sam_file,
    )
    return sam_file


def parse_sam(sam_file: Path) -> Iterator[Alignment]:
    """Yield the primary alignment of every query in a SAM file (secondary and supplementary ones are
    skipped: a contig is judged on its best alignment)."""
    with sam_file.open() as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            flag = int(fields[1])
            if flag & (FLAG_SECONDARY | FLAG_SUPPLEMENTARY):
                continue
            yield Alignment(
                name=fields[0],
                mapped=not flag & FLAG_UNMAPPED,
                reverse=bool(flag & FLAG_REVERSE),
                cigar=parse_cigar(fields[5]),
            )


def is_candidate(cigar: list[tuple[int, str]]) -> bool:
    """True if the differences with the exclusion genome could fit in one primer or probe:

    - a mismatch run, an insertion or a deletion longer than one base;
    - two differences close enough to share one primer, that is with at most `MAX_PRIMER_SPAN - 2`
      matching bases between them, including two that touch;
    - a clipped end of at least `MIN_CLIPPED` bases, which the exclusion genome does not have at all.

    A difference is a mismatch, an insertion or a deletion; a deletion consumes no base of the contig, so
    two mismatches on either side of one are adjacent as far as a primer is concerned.
    """
    matches_since_difference = 0
    after_difference = False
    for length, op in cigar:
        if op in DIFFERENCE_OPS:
            if length > 1:
                return True
            # The pair spans the two differences and the matches between them, so it fits under one
            # primer when `matches + 2 <= MAX_PRIMER_SPAN`
            if after_difference and matches_since_difference <= MAX_PRIMER_SPAN - 2:
                return True
            after_difference, matches_since_difference = True, 0
        elif op in MATCH_OPS:
            matches_since_difference += length
        elif op in CLIP_OPS:
            # The aligner could not place these bases of the contig anywhere in the exclusion genome
            if length >= MIN_CLIPPED:
                return True
            after_difference, matches_since_difference = False, 0
        else:  # N and P: reference-only operations that leave the contig untouched
            after_difference, matches_since_difference = False, 0
    return False


def mark_differences(seq: str, cigar: list[tuple[int, str]]) -> str:
    """The sequence in upper case, with the bases that are not a match to the reference in lower case:
    mismatches, inserted bases, and the clipped ends the aligner could not place. Deleted bases are not in
    the sequence at all: they only show in the cigar string."""
    seq = seq.upper()
    marked: list[str] = []
    position = 0
    for length, op in cigar:
        if op not in QUERY_OPS:  # D, N, H and P consume no query base
            continue
        chunk = seq[position:position + length]
        marked.append(chunk if op in MATCH_OPS else chunk.lower())
        position += length
    marked.append(seq[position:])  # Whatever the cigar did not cover
    return "".join(marked)


def select_candidates(assembly: dict[str, Record], alignments: Iterable[Alignment]) -> dict[str, Candidate]:
    """The contigs worth taking to the blast steps, in the order the alignments came in."""
    candidates: dict[str, Candidate] = {}
    for alignment in alignments:
        contig = assembly.get(alignment.name)
        if contig is None:
            log.warning("The SAM file mentions %s, which is not in the assembly", alignment.name)
            continue
        if not alignment.mapped or not alignment.cigar:
            # Nothing in the exclusion genome looks like this contig: all of it is inclusion-specific.
            candidates[contig.name] = Candidate(contig.name, f"{len(contig.seq)}I", contig.seq.lower())
            continue
        if not is_candidate(alignment.cigar):
            continue
        # The cigar of a reverse alignment describes the reverse complement of the contig: mark the
        # differences on it, then put the contig back in its original orientation.
        seq = reverse_complement(contig.seq) if alignment.reverse else contig.seq
        marked = mark_differences(seq, alignment.cigar)
        if alignment.reverse:
            marked = reverse_complement(marked)
        candidates[contig.name] = Candidate(contig.name, cigar_string(alignment.cigar), marked)
    return candidates


def sort_candidates(candidates: dict[str, Candidate]) -> list[Candidate]:
    """The most promising contigs first: those with the most differing bases, which puts a contig that the
    exclusion genome does not hold at all at the top. Contigs that differ as much keep a stable order."""
    def key(candidate: Candidate) -> tuple[int, int, str]:
        differing = sum(1 for base in candidate.seq if base.islower())
        return -differing, -len(candidate.desc), candidate.name

    return sorted(candidates.values(), key=key)
