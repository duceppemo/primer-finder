#!/usr/bin/env python3
"""What insilicoPCR's -m/--mismatches tolerance actually governs, measured.

A 22-mer primer pair is placed in a synthetic template, and the *template* is mutated: one mismatch at a
known distance from the left primer's 3' end, or a run of them ending at it. Each case is its own genome,
so every row of insilicoPCR's report names the case it came from.

    python mismatch_zones.py /path/to/insilicoPCR-linux-x64/run-insilicoPCR.sh [work_folder]

Prints which cases amplified at each tolerance. Section 1 of SUMMARY.md is this output.
"""

from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

TOLERANCES = (0, 1, 2, 3, 5, 10)
OFFSETS = range(1, 9)   # one mismatch, this far from the left primer's 3' end
RUNS = (2, 3, 4)        # a run of this many mismatches, ending at the 3' end
FLIP = {"A": "G", "C": "T", "G": "A", "T": "C"}


def build(folder: Path) -> list[str]:
    """Write the primer file and one genome per case. Returns the case names, in order."""
    genomes = folder / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    rng = random.Random(3)
    random_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    reverse_complement = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))  # noqa: E731

    left, right = random_seq(22), random_seq(22)
    amplicon = left + random_seq(80) + reverse_complement(right)
    (folder / "primers.fasta").write_text(f">test-F\n{left}\n>test-R\n{right}\n")

    def mutated(offsets):
        bases = list(amplicon)
        for offset in offsets:
            at = len(left) - offset
            bases[at] = FLIP[bases[at]]
        return "".join(bases)

    cases = {"perfect": []}
    cases.update({f"off{offset}": [offset] for offset in OFFSETS})
    cases.update({f"run{n}": list(range(1, n + 1)) for n in RUNS})
    for name, offsets in cases.items():
        (genomes / f"{name}.fasta").write_text(
            f">{name}\n{random_seq(600)}{mutated(offsets)}{random_seq(600)}\n")
    return list(cases)


def amplified(launcher: str, folder: Path, mismatches: int) -> dict[str, str]:
    """Run insilicoPCR at one tolerance: case -> what it reported for the left primer."""
    out = folder / f"m{mismatches}"
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run([launcher, "-i", str(folder / "genomes"), "-p", str(folder / "primers.fasta"),
                    "-o", str(out), "-t", "4", "-m", str(mismatches)],
                   check=True, capture_output=True)
    report = out / "consolidated_report" / "report.tsv"
    found = {}
    with report.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            found[row["Sample"].strip()] = (f"mismatches={row['ForwardMismatches'].strip()}"
                                            f" end={row['ForwardEndMismatch'].strip()}")
    return found


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    launcher = argv[0]
    folder = Path(argv[1] if len(argv) > 1 else "mismatch_zones_work").resolve()
    cases = build(folder)
    results = {m: amplified(launcher, folder, m) for m in TOLERANCES}

    width = max(len(name) for name in cases)
    print(f"{'case':<{width}} | " + " | ".join(f"m={m:<2}" for m in TOLERANCES) + " | reported")
    for name in cases:
        marks = " | ".join(f"{'yes' if name in results[m] else 'no':<4}" for m in TOLERANCES)
        reported = next((results[m][name] for m in TOLERANCES if name in results[m]), "")
        print(f"{name:<{width}} | {marks} | {reported}")
    print("\nThe last two bases of a primer are free at every tolerance; 3 and 4 bases in never amplify at\n"
          "any tolerance; 5 or more bases in is what -m decides. The tolerance is per primer.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
