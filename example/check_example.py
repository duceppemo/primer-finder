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
    both = (sequence, reverse_complement(sequence))

    insertion = truth["insertion"]
    check(insertion in both[0] or insertion in both[1], f"the contig carries the planted insertion {insertion}")

    lowercase = [position for position, base in enumerate(contig.seq) if base.islower()]
    close = [1 for first, second in zip(lowercase, lowercase[1:], strict=False) if second - first < PRIMER_LENGTH]
    check(len(lowercase) >= len(truth["mismatches"]) + len(insertion),
          f"the planted variants are in lower case ({len(lowercase)} bases)")
    check(len(close) >= 2, f"at least three variants within {PRIMER_LENGTH} bases")
    check(contig.desc.startswith("["), f"the header lists the variant positions: {contig.desc[:40]}")

    for group, expected in (("inclusion", True), ("exclusion", False)):
        for genome in find_genomes(args.data / group):
            genome_sequence = "".join(record.seq for record in iter_records(genome)).upper()
            found = any(part in genome_sequence for part in both)
            check(found is expected,
                  f"the contig is {'present in' if expected else 'absent from'} {genome.name}")

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
