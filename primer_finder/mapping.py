"""Mapping the assembled kmers to one exclusion genome, and keeping the contigs whose differences could
carry a selective assay."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from primer_finder import tools
from primer_finder.seqio import Record

log = logging.getLogger(__name__)

# A primer is about 21 bases long: two differences closer than this can sit in the same primer or probe.
PRIMER_LENGTH = 21

# The cigar operations that are a difference with the reference, and those that consume query bases.
DIFFERENCE_OPS = "XID"
QUERY_OPS = "MIS=X"
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

FLAG_UNMAPPED = 0x4
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_SUPPLEMENTARY = 0x800

COMPLEMENT = str.maketrans("ACGTNacgtnRYKMBDHVrykmbdhv", "TGCANtgcanYRMKVHDByrmkvhdb")


@dataclass
class Alignment:
    """One primary SAM record: the query name, whether it mapped, its strand and its cigar."""

    name: str
    mapped: bool
    reverse: bool
    cigar: list[tuple[int, str]]
    length: int


@dataclass
class Candidate:
    """A contig kept for the blast steps: its name, what makes it different (the cigar string, or the
    length of the contig when it does not map at all), and its sequence with the differences in lowercase."""

    name: str
    desc: str
    seq: str

    def as_record(self) -> Record:
        return Record(self.name, self.desc, self.seq)


def reverse_complement(seq: str) -> str:
    """The reverse complement of a sequence, keeping the case of each base."""
    return seq.translate(COMPLEMENT)[::-1]


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
            seq = fields[9]
            cigar = parse_cigar(fields[5])
            length = len(seq) if seq != "*" else sum(n for n, op in cigar if op in QUERY_OPS)
            yield Alignment(
                name=fields[0],
                mapped=not flag & FLAG_UNMAPPED,
                reverse=bool(flag & FLAG_REVERSE),
                cigar=cigar,
                length=length,
            )


def is_candidate(cigar: list[tuple[int, str]]) -> bool:
    """True if the differences with the exclusion genome could fit in one primer or probe:
    an insertion, a deletion or a run of mismatches longer than one base, or two differences less than
    PRIMER_LENGTH bases apart."""
    for index, (length, op) in enumerate(cigar):
        if op in DIFFERENCE_OPS and length > 1:
            return True
        if op in DIFFERENCE_OPS and index + 2 < len(cigar):
            gap_length, gap_op = cigar[index + 1]
            next_op = cigar[index + 2][1]
            if gap_op in "=M" and gap_length < PRIMER_LENGTH and next_op in DIFFERENCE_OPS:
                return True
    return False


def mark_differences(seq: str, cigar: list[tuple[int, str]]) -> str:
    """The sequence in upper case, with the bases that differ from the reference (mismatches and inserted
    bases) in lower case. Deleted bases are not in the sequence: they only show in the cigar string."""
    seq = seq.upper()
    marked: list[str] = []
    position = 0
    for length, op in cigar:
        if op not in QUERY_OPS:  # D, N, H and P consume no query base
            continue
        chunk = seq[position:position + length]
        marked.append(chunk.lower() if op in DIFFERENCE_OPS else chunk)
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
    """The most promising contigs first: those with the most differences, which is to say the longest cigar
    string. Contigs with cigar strings of the same length keep a stable order (by name)."""
    return sorted(candidates.values(), key=lambda candidate: (-len(candidate.desc), candidate.name))
