#!/usr/bin/env python3
"""Check the example results against what make_example.py planted."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from primer_finder.mapping import MAX_PRIMER_SPAN, reverse_complement  # noqa: E402
from primer_finder.seqio import find_genomes, iter_records  # noqa: E402

FAILURES: list[str] = []


def check(condition: bool, message: str) -> None:
    print(f"{'ok  ' if condition else 'FAIL'}  {message}")
    if not condition:
        FAILURES.append(message)


def whole_sequence(path: Path) -> str:
    return "".join(record.seq for record in iter_records(path)).upper()


def locate(contig: str, genome: str) -> tuple[int, bool] | None:
    """Where the contig sits in the genome, and whether it is in the same orientation."""
    start = genome.find(contig)
    if start >= 0:
        return start, True
    start = genome.find(reverse_complement(contig))
    if start >= 0:
        return start, False
    return None


def shorten(positions: set[int], most: int = 8) -> str:
    """A set of positions, short enough to read."""
    listed = sorted(positions)
    if len(listed) <= most:
        return str(listed)
    return f"{listed[:most]}... and {len(listed) - most} more".replace("]...", ", ...]")


def genome_positions(contig_seq: str, start: int, forward: bool) -> set[int]:
    """The genome positions of the lower-case bases of the contig."""
    length = len(contig_seq)
    return {
        start + (index if forward else length - 1 - index)
        for index, base in enumerate(contig_seq)
        if base.islower()
    }


def equivalent_variant_sets(truth: dict, genome: str) -> list[set[int]]:
    """Every marked-position set that describes the same alignment as the planted one.

    An inserted block can sometimes be placed more than one way for the same resulting sequence: if the
    block's first base equals the base that follows it, sliding the whole block one to the right spells
    the genome unchanged. This walks that ambiguity to its end and returns one set per valid placement,
    planted placement first, so a checker can accept any of them and still reject a set that is shifted,
    mis-oriented or incomplete.
    """
    insertion_at = truth["insertion_at"]
    insertion = truth["insertion"]
    mismatches = set(truth["mismatches"])

    # Both ways from the planted placement, because an aligner may put the gap on either side of an
    # ambiguous run: sliding right spells the same genome when the base leaving the block at the front
    # equals the one joining at the back, and sliding left when the base joining at the front equals the
    # one leaving at the back. Bounded by the block length, since a longer run means periodic data this
    # example does not have.
    starts = [insertion_at]
    for direction in (1, -1):
        at = insertion_at
        for _ in range(len(insertion)):
            if direction == 1:
                arriving = genome[at + len(insertion):at + len(insertion) + 1]
                leaving = genome[at:at + 1]
            else:
                arriving = genome[at - 1:at]
                leaving = genome[at + len(insertion) - 1:at + len(insertion)]
            if not arriving or arriving != leaving:
                break
            at += direction
            if at < 0:  # pragma: no cover - not reachable from this planted position
                break
            starts.append(at)
    return [mismatches | set(range(start, start + len(insertion))) for start in starts]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="The folder primer-finder wrote.")
    parser.add_argument("data", type=Path, help="The folder make_example.py wrote.")
    args = parser.parse_args(argv)

    truth = json.loads((args.data / "truth.json").read_text())
    info = json.loads((args.output / "run_info.json").read_text())
    final = list(iter_records(args.output / "final_kmers.fasta"))

    check(len(final) == 1, f"one contig in final_kmers.fasta (got {len(final)})")
    if not final:
        return 1
    contig = final[0]
    sequence = contig.seq.upper()

    insertion = truth["insertion"]
    check(insertion in sequence or insertion in reverse_complement(sequence),
          f"the contig carries the planted insertion {insertion}")
    check(contig.desc.startswith("[") and contig.desc != "[]",
          f"the header lists the variant positions: {contig.desc[:40]}")

    # The lower-case bases must be exactly the planted variants, no others and none missing. Mapping them
    # back to the genome is what catches a shifted or mis-oriented set of positions.
    #
    # With one allowance, which is a property of the planted sequence rather than a weakening of the
    # check: the insertion is `GGCATT...AACC` and the base that follows it in the genome is also `G`, so
    # inserting those 15 bases *before* that G and inserting a 15-base rotation of them *after* it spell
    # the same genome. Both are correct alignments of the same pair of sequences, and which one an aligner
    # returns depends on its gap scoring -- megablast and blastn disagree here. So the insertion block is
    # accepted anywhere inside that ambiguity, computed from the sequence itself, while every mismatch
    # position must still be exact.
    found = locate(sequence, whole_sequence(args.data / "inclusion" / "inclusion_1.fasta"))
    check(found is not None, "the contig is present in inclusion_1.fasta")
    if found is not None:
        start, forward = found
        marked = genome_positions(contig.seq, start, forward)
        expected = set(truth["variants"])
        accepted = equivalent_variant_sets(truth, whole_sequence(
            args.data / "inclusion" / "inclusion_1.fasta"))
        if marked in accepted and marked != expected:
            detail = (f" (the insertion sits {accepted.index(marked)} base(s) along from where it was "
                      "planted, which spells the same genome)")
        elif marked == expected:
            detail = ""
        else:
            detail = f" (missing {shorten(expected - marked)}, extra {shorten(marked - expected)})"
        check(marked in accepted,
              "the lower-case bases are exactly the planted variants" + detail)
        clustered = sorted(set(truth["mismatches"]) & marked - set(range(
            truth["insertion_at"], truth["insertion_at"] + len(insertion))))
        close = [second - first for first, second in zip(clustered, clustered[1:], strict=False)
                 if second - first < MAX_PRIMER_SPAN]
        check(len(close) >= 2,
              f"three of the planted mismatches are marked within {MAX_PRIMER_SPAN} bases ({clustered})")

    for group, expected_presence in (("inclusion", True), ("exclusion", False)):
        for genome in find_genomes(args.data / group):
            present = locate(sequence, whole_sequence(genome)) is not None
            check(present is expected_presence,
                  f"the contig is {'present in' if expected_presence else 'absent from'} {genome.name}")

    # The run has two rescues built in -- the last-resort pass at one difference, and the unitig walk when
    # the assembler makes nothing of the kmers -- and either would carry this example on its own. So a
    # regression in the strict pass, or a SKESA collapse, would leave every check above passing. The
    # expected run uses neither.
    # Read with `get`, so a record missing a field fails the check it belongs to instead of raising
    parameters = info.get("parameters", {})
    asked, used = parameters.get("assembler"), parameters.get("assembler_used")
    check(info.get("last_resort") is False,
          f"the strict pass found the region, without the last-resort pass ({info.get('last_resort')})")
    check(parameters.get("min_differences") == 2,
          f"the reported region needs two differences (got {parameters.get('min_differences')})")
    check(asked is not None and used == asked,
          f"the chosen assembler produced the contigs (asked for {asked}, {used} did it)")
    check("last_resort" not in contig.desc.split(),
          f"the contig is not tagged as a last-resort candidate: {contig.desc[:60]}")

    counts = info["counts"]
    check(counts["final"] == 1, f"run_info.json counts one final contig (got {counts['final']})")
    check(counts["kmers"] > 0, f"inclusion-specific kmers were found ({counts['kmers']})")
    check(Path(info["reference"]).name == "exclusion_1.fasta",
          f"the first exclusion genome was used as the reference ({Path(info['reference']).name})")

    if FAILURES:
        print(f"\n{len(FAILURES)} check(s) failed", file=sys.stderr)
        return 1
    print("\nAll checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
