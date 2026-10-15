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
    sets: list[set[int]] = []
    shift = 0
    while True:
        block = set(range(insertion_at + shift, insertion_at + shift + len(insertion)))
        sets.append(mismatches | block)
        # Sliding one further is the same sequence only when the base leaving the block at the front
        # equals the base joining it at the back
        rotated = insertion[shift % len(insertion):] + insertion[:shift % len(insertion)]
        following = genome[insertion_at + shift + len(insertion):insertion_at + shift + len(insertion) + 1]
        if not following or rotated[0] != following:
            return sets
        shift += 1
        if shift > len(insertion):  # pragma: no cover - a run that long is not in this data
            return sets


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
