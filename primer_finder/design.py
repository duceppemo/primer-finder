"""Designing assays on the candidate regions with Primer3, and judging where their specificity comes from.

A candidate region is inclusion-specific as a whole: every kmer of it is absent from the exclusion genomes.
An assay designed inside it is not specific by itself, though. Two things can make it so:

- the amplicon is absent from every exclusion genome, so there is nothing to amplify. This is how the
  published assays this tool is validated against work;
- the amplicon is there, but an oligo sits on a base that differs, which is what `final_kmers.fasta` marks
  in lower case. A mismatch under a primer, and especially at its 3' end, is what discriminates.

So Primer3 is asked for several assays per region, and each one is then checked against the exclusion
genomes. An assay whose amplicon exists in an exclusion genome and whose oligos cover no difference is
dropped: it would amplify both groups.
"""

from __future__ import annotations

import json
import logging
import re
import time
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from primer_finder import PrimerFinderError, __version__, blast, insilico, seqio, tools
from primer_finder.pipeline import add_log_file

log = logging.getLogger(__name__)

PROGRAM = "primer3_core"
# What Primer3 is asked for, unless the command line says otherwise.
DEFAULT_PRODUCT_SIZE = "70-150"
DEFAULT_ASSAYS_PER_REGION = 3
# A hit covering at least this much of an amplicon, at this identity, counts as "the exclusion genome has it".
MIN_AMPLICON_COVERAGE = 0.9
MIN_AMPLICON_IDENTITY = 90.0
# Blast hits to ask for per amplicon and genome, enough to count the copies of a repeated target.
MAX_COPIES = 100

BOULDER_LINE = re.compile(r"^([A-Z0-9_]+)=(.*)$")

# Primer3 is also asked for primers whose 3' end is forced onto a run of differences, which is what makes
# an assay allele-specific. The request id carries the region and what was forced, separated by this.
FORCE_SEP = "~force_"
MAX_FORCED_RUNS = 4  # Per region, longest runs first
MIN_PRIMER_SIZE = 18  # As the settings ask for: a forced end needs this much room behind it
MIN_PROBE_SIZE = 18
# A forced request only cares about the neighbourhood of the run it aims at. Letting Primer3 consider the
# whole of a several-kilobase region makes it slower by two orders of magnitude, since it weighs every
# pair of positions; this is how far either side of the forced position it is allowed to look.
FORCED_WINDOW_MARGIN = 60
# Each request is given its own run of Primer3 and its own time limit, so that a slow one costs only
# itself. Pinning the probe over a position (SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST) is the expensive
# constraint by far: on one real region, eleven of twelve requests took a tenth of a second each and two
# probe ones took 27 and 114 seconds. The rest are fast enough that their limit is only a safety net.
PRIMER3_TIMEOUT = 60.0
PROBE_TIMEOUT = 15.0

# Primer3's pair penalty is how far an assay is from the ideal chemistry. Assays are compared by band so
# that a probe covering more differences cannot win over one with clearly better chemistry; within a band
# the differences decide, and the exact penalty breaks what is left.
PENALTY_BANDS = (1.0, 2.0, 4.0)


# How close to the 3' end a difference still counts as "near the end", where it hinders extension most.
THREE_PRIME_WINDOW = 5

# A mismatch is worth more where the base it replaces is a G or a C. Those pair with three hydrogen bonds
# against two for A and T, so breaking one costs more, and which bases are involved is known to change what
# a mismatch does (Sharma et al. 2022, doi:10.1016/j.jmoldx.2022.08.005). This is reasoning from the
# chemistry, not a measurement: it only ever breaks a tie between assays with the same number of
# differences. A position whose base is not known counts in between.
STRONG_BASE_WEIGHT = 1.0   # the exclusion genome has a G or a C there
WEAK_BASE_WEIGHT = 0.5     # it has an A or a T
UNKNOWN_BASE_WEIGHT = 0.75


@dataclass
class Oligo:
    """One designed oligo, placed on the region it was designed in.

    `reverse` is set for the right primer, which anneals to the other strand: its 3' end is at its lowest
    coordinate in the region, not its highest. Where the differences sit relative to that end is what
    decides whether an exclusion template still amplifies, so the two cases must not be mixed up.
    """

    name: str
    sequence: str
    start: int  # 0-based in the region, of the leftmost base
    length: int
    tm: float
    gc: float
    variants: list[int] = field(default_factory=list)  # Positions in the region it covers that differ
    reverse: bool = False
    # What the exclusion genomes have at those positions, when it is known
    variant_bases: dict[int, str] = field(default_factory=dict)

    @property
    def end(self) -> int:
        return self.start + self.length

    @property
    def three_prime(self) -> int:
        """The position of the 3'-most base, in region coordinates."""
        return self.start if self.reverse else self.end - 1

    def distance_from_three_prime(self, position: int) -> int:
        return position - self.three_prime if self.reverse else self.three_prime - position

    @property
    def terminal_run(self) -> int:
        """How many bases differ in an unbroken run ending at the 3' end.

        One mismatch at the 3' end does not always stop amplification, so a run of two or three is worth
        much more than a single one. Zero means the 3'-most base matches the exclusion template.
        """
        covered = set(self.variants)
        run = 0
        step = 1 if self.reverse else -1
        position = self.three_prime
        while self.start <= position < self.end and position in covered:
            run += 1
            position += step
        return run

    @property
    def near_three_prime(self) -> int:
        """Differences within THREE_PRIME_WINDOW bases of the 3' end, run or not."""
        return sum(1 for position in self.variants
                   if self.distance_from_three_prime(position) < THREE_PRIME_WINDOW)

    @property
    def strong_variants(self) -> int:
        """Differences where the exclusion genome has a G or a C."""
        return sum(1 for position in self.variants if self.variant_bases.get(position) in ("G", "C"))

    @property
    def variant_weight(self) -> float:
        """The differences this oligo covers, weighted by the base each one replaces."""
        total = 0.0
        for position in self.variants:
            base = self.variant_bases.get(position)
            if base in ("G", "C"):
                total += STRONG_BASE_WEIGHT
            elif base in ("A", "T"):
                total += WEAK_BASE_WEIGHT
            else:
                total += UNKNOWN_BASE_WEIGHT
        return total


@dataclass
class Assay:
    """A designed assay: two primers, an optional probe, and what makes it specific."""

    region: str
    number: int
    forward: Oligo
    reverse: Oligo
    probe: Oligo | None
    product_size: int
    penalty: float
    exclusion_genomes_with_amplicon: int | None = None  # None until it has been checked
    inclusion_copies: list[int] | None = None  # Copies of the amplicon in each inclusion genome
    amplicon: str = ""

    @property
    def name(self) -> str:
        return f"{self.region}_assay{self.number}"

    @property
    def oligos(self) -> list[Oligo]:
        return [self.forward, self.reverse] + ([self.probe] if self.probe else [])

    @property
    def variants_covered(self) -> int:
        return sum(len(oligo.variants) for oligo in self.oligos)

    @property
    def primer_variants_covered(self) -> int:
        """Differences under the primers, which discriminate better than one under the probe."""
        return len(self.forward.variants) + len(self.reverse.variants)

    @property
    def best_terminal_run(self) -> int:
        """The longest run of differences ending at the 3' end of either primer. This is the number that
        matters most: a single 3'-terminal mismatch is often not enough to stop amplification."""
        return max(self.forward.terminal_run, self.reverse.terminal_run)

    @property
    def near_three_prime(self) -> int:
        """Differences close to the 3' end of either primer."""
        return self.forward.near_three_prime + self.reverse.near_three_prime

    @property
    def probe_variants_covered(self) -> int:
        return len(self.probe.variants) if self.probe else 0

    @property
    def primer_variant_weight(self) -> float:
        """The differences under the primers, weighted by which base each one replaces."""
        return self.forward.variant_weight + self.reverse.variant_weight

    @property
    def primer_strong_variants(self) -> int:
        return self.forward.strong_variants + self.reverse.strong_variants

    @property
    def penalty_band(self) -> int:
        """Which band of Primer3's pair penalty this assay falls in; lower is better chemistry."""
        for band, limit in enumerate(PENALTY_BANDS):
            if self.penalty < limit:
                return band
        return len(PENALTY_BANDS)

    @property
    def min_inclusion_copies(self) -> int:
        """How many copies of the amplicon the inclusion genome with the fewest of them holds.

        Two or more in every genome means a repeated target, which usually improves the limit of detection:
        more template per cell. Note that the default -d 1 of the `find` command throws such regions away
        before this step ever sees them, so this is only above 1 when the run used -d 2 or more.
        """
        return min(self.inclusion_copies) if self.inclusion_copies else 0

    @property
    def specific_by_absence(self) -> bool:
        return self.exclusion_genomes_with_amplicon == 0

    @property
    def usable(self) -> bool:
        """An assay that could tell the two groups apart: nothing to amplify in the exclusion genomes, or
        a difference under one of its oligos."""
        if self.exclusion_genomes_with_amplicon is None:
            return True  # Not checked: the caller decides
        return self.specific_by_absence or self.variants_covered > 0


def build_request(record: seqio.Record, product_size: str, how_many: int,
                  extra: dict[str, object] | None = None) -> str:
    """One Primer3 boulder-IO record for a candidate region, asking for a probe as well as primers."""
    settings: dict[str, object] = {
        "SEQUENCE_ID": record.name,
        "SEQUENCE_TEMPLATE": record.seq.upper(),
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 1,  # The probe of a qPCR assay
        "PRIMER_PRODUCT_SIZE_RANGE": product_size,
        "PRIMER_NUM_RETURN": how_many,
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_MIN_TM": 58.0,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_MIN_GC": 30.0,
        "PRIMER_MAX_GC": 70.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        # What an acceptable probe is. Primer3 returns nothing outside these, so a probe pushed onto a run
        # of differences is still a usable probe or there is no assay from that request at all.
        "PRIMER_INTERNAL_MIN_SIZE": 18,
        "PRIMER_INTERNAL_OPT_SIZE": 22,
        "PRIMER_INTERNAL_MAX_SIZE": 27,
        "PRIMER_INTERNAL_MIN_TM": 62.0,
        "PRIMER_INTERNAL_OPT_TM": 68.0,
        "PRIMER_INTERNAL_MAX_TM": 72.0,
        "PRIMER_INTERNAL_MIN_GC": 30.0,
        "PRIMER_INTERNAL_MAX_GC": 80.0,
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    settings.update(extra or {})
    return "".join(f"{key}={value}\n" for key, value in settings.items()) + "=\n"


def parse_records(text: str) -> list[dict[str, str]]:
    """Primer3's boulder-IO output, one dictionary per record."""
    records, current = [], {}
    for line in text.splitlines():
        if line == "=":
            if current:
                records.append(current)
            current = {}
            continue
        match = BOULDER_LINE.match(line)
        if match:
            current[match.group(1)] = match.group(2)
    if current:
        records.append(current)
    return records


def _oligo(fields: dict[str, str], kind: str, number: int, name: str, variants: Sequence[int]) -> Oligo | None:
    """One oligo out of a Primer3 record. `kind` is LEFT, RIGHT or INTERNAL."""
    position = fields.get(f"PRIMER_{kind}_{number}")
    sequence = fields.get(f"PRIMER_{kind}_{number}_SEQUENCE")
    if not position or not sequence:
        return None
    start, length = (int(part) for part in position.split(","))
    if kind == "RIGHT":  # Primer3 gives the rightmost base of a right primer
        start = start - length + 1
    covered = sorted(p for p in variants if start <= p < start + length)
    return Oligo(
        name=name,
        sequence=sequence,
        start=start,
        length=length,
        tm=float(fields.get(f"PRIMER_{kind}_{number}_TM", "nan")),
        gc=float(fields.get(f"PRIMER_{kind}_{number}_GC_PERCENT", "nan")),
        variants=covered,
        reverse=kind == "RIGHT",
    )


def assays_of(record: seqio.Record, fields: dict[str, str]) -> list[Assay]:
    """Every assay Primer3 returned for one region."""
    variants = [position for position, base in enumerate(record.seq) if base.islower()]
    assays = []
    for number in range(int(fields.get("PRIMER_PAIR_NUM_RETURNED", "0"))):
        forward = _oligo(fields, "LEFT", number, f"{record.name}_assay{number}-F", variants)
        reverse = _oligo(fields, "RIGHT", number, f"{record.name}_assay{number}-R", variants)
        if forward is None or reverse is None:  # pragma: no cover - Primer3 returns pairs
            continue
        probe = _oligo(fields, "INTERNAL", number, f"{record.name}_assay{number}-P", variants)
        assays.append(Assay(
            region=record.name,
            number=number,
            forward=forward,
            reverse=reverse,
            probe=probe,
            product_size=int(fields[f"PRIMER_PAIR_{number}_PRODUCT_SIZE"]),
            penalty=float(fields.get(f"PRIMER_PAIR_{number}_PENALTY", "nan")),
            amplicon=record.seq[forward.start:reverse.end].upper(),
        ))
    return assays


def variant_runs(seq: str) -> list[tuple[int, int]]:
    """The runs of consecutive differing (lower-case) bases, longest first, then leftmost.

    A run is what a primer's 3' end should sit on: three differences in a row at the 3' end stop extension
    on the exclusion template far more reliably than one.
    """
    runs: list[tuple[int, int]] = []
    start = None
    for position, base in enumerate(seq):
        if base.islower():
            start = position if start is None else start
        elif start is not None:
            runs.append((start, position - start))
            start = None
    if start is not None:
        runs.append((start, len(seq) - start))
    return sorted(runs, key=lambda run: (-run[1], run[0]))


def forced_requests(record: seqio.Record, product_size: str, how_many: int) -> list[str]:
    """Primer3 records aimed at a run of differences: up to three per run.

    `SEQUENCE_FORCE_LEFT_END` and `SEQUENCE_FORCE_RIGHT_END` name the position of a primer's 3'-most base.
    For the left primer that is the last base of the run; for the right primer, which reads the other way,
    it is the first. `SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST` makes the probe straddle a position, which
    is used to put it over the middle of the run.

    A request is only made when there is room for the rest of the assay on either side. Asking Primer3 for
    something that cannot exist -- a probe at base 15 of a region, with no room for a primer before it --
    does not return quickly: it searches the whole template first.

    Any of them may still come back empty, which is not an error: the run may sit where no oligo of the
    required size and melting temperature can be placed.
    """
    requests = []
    length = len(record.seq)
    smallest_product = minimum_product(product_size)
    for start, run in variant_runs(record.seq)[:MAX_FORCED_RUNS]:
        left_end = start + run - 1
        # A left primer ending here needs its own length behind it and room for an amplicon ahead
        if left_end >= MIN_PRIMER_SIZE - 1 and left_end + smallest_product <= length:
            requests.append(build_request(record, product_size, how_many, {
                "SEQUENCE_ID": f"{record.name}{FORCE_SEP}left{left_end}",
                "SEQUENCE_FORCE_LEFT_END": left_end,
                "SEQUENCE_INCLUDED_REGION": window(record, left_end, product_size),
            }))
        # A right primer ending here needs room for an amplicon behind it and its own length ahead
        if start >= smallest_product - 1 and start + MIN_PRIMER_SIZE <= length:
            requests.append(build_request(record, product_size, how_many, {
                "SEQUENCE_ID": f"{record.name}{FORCE_SEP}right{start}",
                "SEQUENCE_FORCE_RIGHT_END": start,
                "SEQUENCE_INCLUDED_REGION": window(record, start, product_size),
            }))
        # And one with the probe pinned over the middle of the run, so that it covers as many differences
        # as it can: a probe that cannot bind the exclusion template is what a qPCR assay relies on. It
        # needs a primer and half a probe on either side of it.
        middle = start + run // 2
        margin = MIN_PRIMER_SIZE + MIN_PROBE_SIZE
        if margin <= middle <= length - margin:
            requests.append(build_request(record, product_size, how_many, {
                "SEQUENCE_ID": f"{record.name}{FORCE_SEP}probe{middle}",
                "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST": middle,
                "SEQUENCE_INCLUDED_REGION": window(record, middle, product_size),
                "PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 3,
                "PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 3,
            }))
    return requests


def window(record: seqio.Record, position: int, product_size: str) -> str:
    """The stretch of the region Primer3 may place oligos in, around a forced position, as `start,length`.

    Wide enough for the largest amplicon asked for on either side of the position, and no wider: the work
    Primer3 does grows with the square of what it is given to look at.
    """
    reach = maximum_product(product_size) + FORCED_WINDOW_MARGIN
    start = max(0, position - reach)
    end = min(len(record.seq), position + reach)
    return f"{start},{end - start}"


def maximum_product(product_size: str) -> int:
    """The largest amplicon the product size range allows."""
    try:
        return max(int(part.split("-")[-1]) for part in product_size.split())
    except (ValueError, IndexError):  # pragma: no cover - the command line checks the format
        return 150


def minimum_product(product_size: str) -> int:
    """The smallest amplicon the product size range allows."""
    try:
        return min(int(part.split("-")[0]) for part in product_size.split())
    except (ValueError, IndexError):  # pragma: no cover - the command line checks the format
        return 70


def requests_for(record: seqio.Record, product_size: str, how_many: int,
                 force_ends: bool) -> list[tuple[str, float]]:
    """Every Primer3 request for one region, with how long each may take."""
    requests = [(build_request(record, product_size, how_many), PRIMER3_TIMEOUT)]
    if force_ends:
        for request in forced_requests(record, product_size, how_many):
            slow = "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST" in request
            requests.append((request, PROBE_TIMEOUT if slow else PRIMER3_TIMEOUT))
    return requests


def design(regions: Iterable[seqio.Record], work_dir: Path, product_size: str = DEFAULT_PRODUCT_SIZE,
           how_many: int = DEFAULT_ASSAYS_PER_REGION, force_ends: bool = True,
           threads: int = 1) -> list[Assay]:
    """Run Primer3 over every region and return the assays it proposes.

    Two kinds of request go in: a plain one, which lets Primer3 pick the best chemistry anywhere in the
    region, and one per run of differences with a primer's 3' end or the probe forced onto it. The second
    kind is what produces allele-specific assays; it often returns nothing, which is not an error.

    Each request is run on its own, several at a time, with a time limit. Primer3 can spend minutes on a
    request whose answer is hard to find, and one that runs out of time is dropped with a note rather than
    holding up the rest.
    """
    regions = list(regions)
    work_dir.mkdir(parents=True, exist_ok=True)
    work: list[tuple[int, seqio.Record, str, float]] = []
    for index, record in enumerate(regions):
        for number, (request, timeout) in enumerate(
                requests_for(record, product_size, how_many, force_ends)):
            work.append((f"{index:04d}_{number:02d}", record, request, timeout))

    def one(name: str, record: seqio.Record, request: str, timeout: float) -> list[Assay] | None:
        request_file = work_dir / f"{name}_input.txt"
        output = work_dir / f"{name}_output.txt"
        request_file.write_text(request)
        try:
            tools.run([PROGRAM, request_file], stdout_path=output, timeout=timeout)
        except tools.ToolTimeout:
            log.debug("Primer3 did not answer within %g s for %s; that request is dropped (%s)",
                      timeout, record.name, request_file)
            return None  # The caller counts these
        proposed: list[Assay] = []
        for fields in parse_records(output.read_text()):
            if fields.get("PRIMER_ERROR"):
                raise PrimerFinderError(
                    f"Primer3 failed on {fields.get('SEQUENCE_ID', record.name)}: {fields['PRIMER_ERROR']}"
                )
            proposed.extend(assays_of(record, fields))
        return proposed

    found: dict[str, list[Assay]] = {record.name: [] for record in regions}
    timed_out = 0
    for (_, record, _, _), proposed in zip(work, blast.parallel(work, one, threads), strict=True):
        if proposed is None:
            timed_out += 1
        else:
            found[record.name].extend(proposed)
    if timed_out:
        log.info("%d Primer3 request(s) ran out of time and were dropped; the regions they belong to are "
                 "still designed on from their other requests", timed_out)

    assays: list[Assay] = []
    for name, proposed in found.items():
        # The kinds of request overlap, so the same assay can come back more than once
        seen = set()
        number = 0
        for assay in proposed:
            key = (assay.forward.sequence, assay.reverse.sequence,
                   assay.probe.sequence if assay.probe else "")
            if key in seen:
                continue
            seen.add(key)
            assay.number = number
            for oligo, suffix in ((assay.forward, "F"), (assay.reverse, "R"), (assay.probe, "P")):
                if oligo is not None:
                    oligo.name = f"{name}_assay{number}-{suffix}"
            number += 1
            assays.append(assay)
    return assays


def amplicon_copies(assays: Sequence[Assay], genomes: Sequence[Path], work_dir: Path, threads: int,
                    with_bases: bool = True) -> tuple[dict[str, list[int]], dict[str, dict[int, str]]]:
    """How many copies of each amplicon each genome holds, and what those genomes have where they differ.

    A copy is a blast hit covering nearly the whole amplicon at high identity, so a genome that holds the
    target twice gives two. In the exclusion group that answers "is there anything to amplify"; in the
    inclusion group it is the copy number of the target, which bears on the limit of detection.

    The second value maps, for each assay, a position of the region to the base the genomes have there,
    taking the commonest when they disagree. It is only meaningful for the exclusion group, whose bases are
    what an oligo has to mismatch, so `with_bases` is off for the inclusion group: the aligned sequences are
    then not asked for and not read.
    """
    counts = {assay.name: [0] * len(genomes) for assay in assays}
    bases: dict[str, dict[int, str]] = {assay.name: {} for assay in assays}
    if not assays or not genomes:
        return counts, bases
    work_dir.mkdir(parents=True, exist_ok=True)
    query = work_dir / "amplicons.fasta"
    seqio.write_fasta(query, [seqio.Record(assay.name, "", assay.amplicon) for assay in assays])
    sizes = {assay.name: len(assay.amplicon) for assay in assays}
    starts = {assay.name: assay.forward.start for assay in assays}  # The amplicon starts here in the region
    fields = ("qseqid", "qstart", "qend", "evalue", "length", "pident")
    if with_bases:
        fields += ("qseq", "sseq")

    def one(index: int, genome: Path) -> tuple[dict[str, int], dict[str, dict[int, str]]]:
        folder = blast.genome_folder(work_dir, index, genome)
        db = blast.make_db(genome, folder)
        hits = blast.parse_hits(blast.blastn(db, query, folder / blast.HITS_NAME, fields,
                                             max_targets=MAX_COPIES), fields)
        found: dict[str, int] = {}
        seen: dict[str, dict[int, str]] = {}
        for hit in hits:
            if hit.evalue > blast.MAX_EVALUE:
                continue
            # A copy is the amplicon in one piece: that is what a primer pair can amplify
            if hit.length >= MIN_AMPLICON_COVERAGE * sizes[hit.query] and hit.identity >= MIN_AMPLICON_IDENTITY:
                found[hit.query] = found.get(hit.query, 0) + 1
                if with_bases:
                    offset = starts[hit.query]
                    seen.setdefault(hit.query, {}).update(
                        {position + offset: base for position, base in blast.variant_bases(hit).items()}
                    )
        return found, seen

    votes: dict[str, dict[int, Counter]] = {assay.name: {} for assay in assays}
    for index, (found, seen) in enumerate(blast.parallel(enumerate(genomes), one, threads)):
        for name, number in found.items():
            counts[name][index] = number
        for name, positions in seen.items():
            for position, base in positions.items():
                votes[name].setdefault(position, Counter())[base] += 1
    for name, positions in votes.items():
        bases[name] = {position: counter.most_common(1)[0][0] for position, counter in positions.items()}
    return counts, bases


def check_against_exclusion(assays: Sequence[Assay], exclusion: Sequence[Path], work_dir: Path,
                            threads: int) -> None:
    """Count how many exclusion genomes hold each amplicon, and note what they have where an oligo differs
    from them. Sets both in place."""
    counts, bases = amplicon_copies(assays, exclusion, work_dir, threads)
    for assay in assays:
        assay.exclusion_genomes_with_amplicon = sum(1 for number in counts[assay.name] if number)
        for oligo in assay.oligos:
            oligo.variant_bases = {position: base for position, base in bases[assay.name].items()
                                   if oligo.start <= position < oligo.end}


def count_inclusion_copies(assays: Sequence[Assay], inclusion: Sequence[Path], work_dir: Path,
                           threads: int) -> None:
    """Count the copies of each amplicon in every inclusion genome. Sets them in place."""
    counts, _ = amplicon_copies(assays, inclusion, work_dir, threads, with_bases=False)
    for assay in assays:
        assay.inclusion_copies = counts[assay.name]


def rank(assays: Iterable[Assay]) -> list[Assay]:
    """The assays worth looking at first. A heuristic order, not a prediction of what works at the bench.

    An assay whose amplicon no exclusion genome holds comes first: there is nothing for it to amplify, and
    that does not depend on how a reaction behaves. The rest rely on differences under their primers, and
    are ordered the way allele-specific design is usually approached:

    - how many differences the primers cover in total. A single mismatch, even at the 3' end, is often not
      enough on its own; adding a second one is the basis of double-mismatch allele-specific qPCR
      (Lefever et al. 2019, doi:10.1038/s41598-019-38581-z);
    - then the longest run of them ending at a 3' end, and how many sit near one, since position changes
      how much a mismatch costs (Sharma et al. 2022, doi:10.1016/j.jmoldx.2022.08.005);
    - then Primer3's penalty by band, so a probe covering more differences never beats an assay that is
      clearly better made;
    - then the copies of the amplicon per inclusion genome, which bear on the limit of detection;
    - then the probe's differences, the exact penalty, and the names, so that a run is reproducible.

    How much any mismatch actually costs depends on the assay -- annealing temperature, polymerase,
    magnesium, cycling -- not on the sequence alone. See docs/wiki/Designing-assays.md.
    """
    return sorted(assays, key=rank_key)





# ---------------------------------------------------------------------------------------------------
# The `design` command: Primer3 over a finished run, then the checks


ASSAYS_NAME = "assays.tsv"
DESIGN_INFO_NAME = "design_info.json"
LOG_NAME = "primer_finder_design.log"
DEFAULT_MAX_REGIONS = 50

COLUMNS = (
    "assay", "region", "specific_by", "exclusion_genomes_with_amplicon", "best_terminal_run",
    "inclusion_copies_min", "inclusion_copies_max", "product_size", "penalty", "penalty_band",
    "primer_variant_weight",
    "forward", "forward_start", "forward_tm", "forward_gc", "forward_variants",
    "forward_strong_variants", "forward_terminal_run", "forward_near_3prime",
    "reverse", "reverse_start", "reverse_tm", "reverse_gc", "reverse_variants",
    "reverse_strong_variants", "reverse_terminal_run", "reverse_near_3prime",
    "probe", "probe_start", "probe_tm", "probe_gc", "probe_variants",
)
VERDICT_COLUMNS = ("inclusion_total", "exclusion_total",
                   "qpcr_inclusion_amplified", "qpcr_inclusion_percent", "qpcr_exclusion_amplified",
                   "qpcr_selective",
                   "pcr_inclusion_amplified", "pcr_inclusion_percent", "pcr_exclusion_amplified",
                   "pcr_selective")


@dataclass
class DesignSettings:
    """What the `design` command needs."""

    results: Path  # A finished primer-finder output folder
    output: Path
    threads: int
    inclusion: Path | None = None  # Default: what the run recorded
    exclusion: Path | None = None
    product_size: str = DEFAULT_PRODUCT_SIZE
    assays_per_region: int = DEFAULT_ASSAYS_PER_REGION
    max_regions: int = DEFAULT_MAX_REGIONS  # 0 for every region
    insilico_pcr: Path | None = None
    mismatches: int = 0
    command_line: list[str] = field(default_factory=list)


def specific_by(assay: Assay) -> str:
    if assay.specific_by_absence:
        return "absence"
    if assay.variants_covered:
        return "difference"
    return "nothing"


def assay_row(assay: Assay, verdict: dict[str, object] | None = None) -> dict[str, object]:
    row: dict[str, object] = {
        "assay": assay.name, "region": assay.region, "specific_by": specific_by(assay),
        "exclusion_genomes_with_amplicon": assay.exclusion_genomes_with_amplicon,
        "best_terminal_run": assay.best_terminal_run,
        "inclusion_copies_min": assay.min_inclusion_copies,
        "inclusion_copies_max": max(assay.inclusion_copies) if assay.inclusion_copies else 0,
        "product_size": assay.product_size, "penalty": f"{assay.penalty:.4f}",
        "penalty_band": assay.penalty_band,
        "primer_variant_weight": f"{assay.primer_variant_weight:g}",
    }
    for part, oligo in (("forward", assay.forward), ("reverse", assay.reverse), ("probe", assay.probe)):
        row[part] = oligo.sequence if oligo else ""
        row[f"{part}_start"] = oligo.start if oligo else ""
        row[f"{part}_tm"] = f"{oligo.tm:.1f}" if oligo else ""
        row[f"{part}_gc"] = f"{oligo.gc:.1f}" if oligo else ""
        row[f"{part}_variants"] = len(oligo.variants) if oligo else ""
        if part != "probe":  # Where the differences sit matters only for the primers
            row[f"{part}_strong_variants"] = oligo.strong_variants if oligo else ""
            row[f"{part}_terminal_run"] = oligo.terminal_run if oligo else ""
            row[f"{part}_near_3prime"] = oligo.near_three_prime if oligo else ""
    for kind, one in (verdict or {}).items():
        row.update({
            "inclusion_total": one.inclusion_total,
            "exclusion_total": one.exclusion_total,
            f"{kind}_inclusion_amplified": one.inclusion_amplified,
            f"{kind}_inclusion_percent": one.percent,
            f"{kind}_exclusion_amplified": one.exclusion_amplified,
            f"{kind}_selective": one.label,
        })
    return row


def write_assays(path: Path, assays: Sequence[Assay], verdicts: dict[str, object] | None = None) -> None:
    import csv

    columns = list(COLUMNS) + (list(VERDICT_COLUMNS) if verdicts else [])
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for assay in assays:
            writer.writerow(assay_row(assay, (verdicts or {}).get(assay.name)))


def threshold_of(settings: DesignSettings) -> float:
    """The -p/--min-inclusion the run used, which is what an assay on its regions has to live up to.

    A region found with -p 0.9 is one that a tenth of the inclusion genomes may lack, so an assay on it
    cannot be expected to amplify all of them: in silico PCR judges it against the same fraction.
    """
    info_file = settings.results / "run_info.json"
    if not info_file.is_file():
        return 1.0
    try:
        recorded = json.loads(info_file.read_text()).get("parameters", {}).get("min_inclusion", 1.0)
        return float(recorded)
    except (ValueError, TypeError, json.JSONDecodeError):  # pragma: no cover - a hand-edited file
        return 1.0


def genome_folders(settings: DesignSettings) -> tuple[Path, Path]:
    """The two input folders: the ones given, or the ones the run recorded in run_info.json."""
    inclusion, exclusion = settings.inclusion, settings.exclusion
    if inclusion is None or exclusion is None:
        info_file = settings.results / "run_info.json"
        if not info_file.is_file():
            raise PrimerFinderError(
                f"{info_file} is missing, so the genome folders are not known. Give them with -i and -e."
            )
        info = json.loads(info_file.read_text())
        recorded = info.get("parameters", {})
        inclusion = inclusion or Path(recorded.get("inclusion", ""))
        exclusion = exclusion or Path(recorded.get("exclusion", ""))
    for group, folder in (("inclusion", inclusion), ("exclusion", exclusion)):
        if not folder or not Path(folder).is_dir():
            raise PrimerFinderError(
                f"The {group} folder of that run is gone ({folder}). Give it with "
                f"{'-i' if group == 'inclusion' else '-e'}."
            )
    return Path(inclusion), Path(exclusion)


def run(settings: DesignSettings) -> int:
    """Design assays on the regions of a finished run, and check them. Returns 0, or raises."""
    started = time.time()
    final = settings.results / "final_kmers.fasta"
    if not final.is_file():
        raise PrimerFinderError(
            f"{final} is missing: give the output folder of a finished primer-finder run."
        )
    inclusion, exclusion = genome_folders(settings)
    tools.require([PROGRAM, "makeblastdb", "blastn"])
    settings.output.mkdir(parents=True, exist_ok=True)
    add_log_file(settings.output / LOG_NAME)
    log.info("primer-finder %s, designing assays", __version__)

    regions = list(seqio.iter_records(final))
    if not regions:
        raise PrimerFinderError(f"{final} holds no candidate region: there is nothing to design on.")
    if settings.max_regions and len(regions) > settings.max_regions:
        log.info("Designing on the first %d of %d regions (--max-regions 0 for all of them)",
                 settings.max_regions, len(regions))
        regions = regions[:settings.max_regions]

    log.info("Running Primer3 on %d region(s), up to %d assay(s) each...", len(regions),
             settings.assays_per_region)
    assays = design(regions, settings.output / "primer3", settings.product_size,
                    settings.assays_per_region, threads=settings.threads)
    log.info("Primer3 proposed %d assay(s) on %d region(s)", len(assays),
             len({assay.region for assay in assays}))
    if not assays:
        raise PrimerFinderError(
            "Primer3 found no assay in any region. A wider --product-size may help; "
            f"its own explanation is in {settings.output / 'primer3' / 'primer3_output.txt'}."
        )

    exclusion_genomes = seqio.require_genomes(exclusion, "exclusion")
    inclusion_genomes = seqio.require_genomes(inclusion, "inclusion")
    log.info("Checking whether the %d amplicon(s) exist in the %d exclusion genomes...",
             len(assays), len(exclusion_genomes))
    check_against_exclusion(assays, exclusion_genomes, settings.output / "amplicons" / "exclusion",
                            settings.threads)
    log.info("Counting the copies of each amplicon in the %d inclusion genomes...", len(inclusion_genomes))
    count_inclusion_copies(assays, inclusion_genomes, settings.output / "amplicons" / "inclusion",
                           settings.threads)
    repeated = sum(1 for assay in assays if assay.min_inclusion_copies > 1)
    if repeated:
        log.info("%d assay(s) amplify a target present more than once in every inclusion genome, which "
                 "usually improves the limit of detection", repeated)
    ranked = rank(assays)
    usable = [assay for assay in ranked if assay.usable]
    by_absence = sum(1 for assay in usable if assay.specific_by_absence)
    log.info("%d assay(s) could tell the groups apart: %d because the exclusion genomes do not hold the "
             "amplicon, %d because an oligo sits on a difference", len(usable), by_absence,
             len(usable) - by_absence)
    dropped = len(ranked) - len(usable)
    if dropped:
        log.info("%d assay(s) would amplify both groups and are listed but not carried further", dropped)
    if not usable:
        raise PrimerFinderError(
            "Every assay Primer3 proposed would amplify the exclusion genomes too. The regions may differ "
            "from them by too little; --assays-per-region above the default gives Primer3 more tries."
        )

    # Both kinds of primer file, because insilicoPCR runs one mode at a time: the qPCR file says whether
    # the whole assay reports positive, the PCR file whether the primer pair alone discriminates.
    primer_files: dict[str, Path] = {}
    for kind, name in (("qpcr", insilico.QPCR_NAME), ("pcr", insilico.PCR_NAME)):
        path = settings.output / name
        written = insilico.write_primer_fasta(usable, path, kind)
        if written:
            primer_files[kind] = path
            log.info("%d assay(s) written for %s: %s", written, kind, path)

    verdicts = None
    if settings.insilico_pcr is not None:
        threshold = threshold_of(settings)
        if threshold < 1:
            log.info("The regions were found with -p %g, so an assay counts when it amplifies that share "
                     "of the inclusion genomes", threshold)
        verdicts = run_insilico_pcr(settings, usable, primer_files, inclusion, exclusion, threshold)
        # Amplifying every inclusion genome comes before merely meeting the threshold
        usable = sorted(usable, key=lambda assay: (
            -sum(2 if verdict.complete else 1 if verdict.selective else 0
                 for verdict in verdicts[assay.name].values()),
            rank_key(assay)))
    else:
        insilico.write_script(settings.output / insilico.SCRIPT_NAME, None, list(primer_files.values()),
                              inclusion, exclusion, settings.threads, settings.mismatches)
        log.info("In silico PCR was not run. To check the assays against both groups: %s",
                 settings.output / insilico.SCRIPT_NAME)

    write_assays(settings.output / ASSAYS_NAME, usable + [a for a in ranked if not a.usable], verdicts)
    (settings.output / DESIGN_INFO_NAME).write_text(json.dumps({
        "version": __version__,
        "command_line": settings.command_line,
        "results": str(settings.results),
        "inclusion": str(inclusion),
        "exclusion": str(exclusion),
        "parameters": {
            "product_size": settings.product_size,
            "assays_per_region": settings.assays_per_region,
            "max_regions": settings.max_regions,
            "mismatches": settings.mismatches,
            "threads": settings.threads,
            "min_inclusion": threshold_of(settings),
        },
        "programs": tools.versions([PROGRAM, "blastn"]),
        "counts": {
            "regions": len(regions),
            "assays": len(ranked),
            "usable": len(ranked) - dropped,
            "specific_by_absence": by_absence,
            "selective": {
                kind: {
                    "complete": sum(1 for modes in verdicts.values()
                                    if kind in modes and modes[kind].complete),
                    "at_threshold": sum(1 for modes in verdicts.values()
                                        if kind in modes and modes[kind].selective),
                }
                for kind in ("qpcr", "pcr")
            } if verdicts else None,
        },
        "seconds": round(time.time() - started, 1),
    }, indent=2) + "\n")
    log.info("Assays: %s", settings.output / ASSAYS_NAME)
    return 0


def rank_key(assay: Assay) -> tuple:
    """The order rank() uses; see its docstring for why it is this way round."""
    return (
        0 if assay.specific_by_absence else 1,
        -assay.primer_variants_covered,
        -assay.primer_variant_weight,  # At equal numbers, a difference over a G or a C is worth more
        -assay.best_terminal_run,
        -assay.near_three_prime,
        assay.penalty_band,            # Chemistry first, so a worse probe cannot win on mismatches alone
        -assay.min_inclusion_copies,   # A repeated target gives a better limit of detection
        -assay.probe_variants_covered,
        assay.penalty,
        assay.region,
        assay.number,
    )


def run_insilico_pcr(settings: DesignSettings, assays: Sequence[Assay], primer_files: dict[str, Path],
                     inclusion: Path, exclusion: Path,
                     threshold: float = 1.0) -> dict[str, dict[str, object]]:
    """Amplify the assays against both groups with insilicoPCR, once per mode, and put the reports
    together. Returns, per assay, a verdict per mode."""
    launcher = insilico.resolve(settings.insilico_pcr)
    insilico.write_script(settings.output / insilico.SCRIPT_NAME, launcher, list(primer_files.values()),
                          inclusion, exclusion, settings.threads, settings.mismatches)
    verdicts: dict[str, dict[str, object]] = {assay.name: {} for assay in assays}
    for kind, primer_file in primer_files.items():
        reports: dict[str, dict[str, set[str]]] = {}
        for group, folder in (("inclusion", inclusion), ("exclusion", exclusion)):
            log.info("In silico PCR (%s) against the %s genomes...", kind, group)
            output = settings.output / "insilico_pcr" / f"{kind}_{group}"
            insilico.run(launcher, folder, primer_file, output, settings.threads, settings.mismatches)
            reports[group] = insilico.read_report(output)
        mode = insilico.verdicts(assays, inclusion, exclusion, reports["inclusion"], reports["exclusion"],
                                 threshold)
        for name, verdict in mode.items():
            verdicts[name][kind] = verdict
        complete = sum(1 for verdict in mode.values() if verdict.complete)
        partial = sum(1 for verdict in mode.values() if verdict.selective and not verdict.complete)
        counted = sum(1 for assay in assays if kind == "pcr" or assay.probe is not None)
        log.info("In silico PCR (%s): %d of %d assay(s) amplify every inclusion genome and no exclusion "
                 "genome%s", kind, complete, counted,
                 f", and {partial} more reach the -p threshold" if partial else "")
    return verdicts
