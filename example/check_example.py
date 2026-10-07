#!/usr/bin/env python3
"""Check the example results against what make_example.py planted."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from primer_finder.mapping import PRIMER_LENGTH, reverse_complement  # noqa: E402
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
    found = locate(sequence, whole_sequence(args.data / "inclusion" / "inclusion_1.fasta"))
    check(found is not None, "the contig is present in inclusion_1.fasta")
    if found is not None:
        start, forward = found
        marked = genome_positions(contig.seq, start, forward)
        expected = set(truth["variants"])
        check(marked == expected,
              "the lower-case bases are exactly the planted variants"
              + ("" if marked == expected else
                 f" (missing {shorten(expected - marked)}, extra {shorten(marked - expected)})"))
        clustered = sorted(set(truth["mismatches"]) & marked - set(range(
            truth["insertion_at"], truth["insertion_at"] + len(insertion))))
        close = [second - first for first, second in zip(clustered, clustered[1:], strict=False)
                 if second - first < PRIMER_LENGTH]
        check(len(close) >= 2,
              f"three of the planted mismatches are marked within {PRIMER_LENGTH} bases ({clustered})")

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
