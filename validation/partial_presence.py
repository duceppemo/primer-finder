#!/usr/bin/env python3
"""How much of each reported region every inclusion genome actually holds.

    python validation/partial_presence.py /path/to/results -t 32
    python validation/partial_presence.py /path/to/results --oligos ACGT...,TTGC...,GGAC...

A region is reported when enough inclusion genomes hold it, and "holds it" is decided on any blast hit at
all: no fraction of the region is required. This measures what that lets through on a finished run, which
is the evidence for deciding whether the rule should ask for a fraction -- and, with `--oligos`, whether
such a rule would have cost the run a published assay it did recover.

**A coverage rule is applied the way the run applies presence**, which is the only honest way to cost it:
a genome holds a region when it aligns at least the cut-off fraction of it, and the region survives when
at least `inclusion_genomes_required` genomes do. "Some genome holds less than the cut-off" is a different
and much larger set, and the two coincide only at `-p 1.0`. An earlier version of this script reported the
first as the cost of the second; `validation/results/2026-10-16_v1.9.1_presence/SUMMARY.md` carries the
erratum.

Writes `partial_presence.tsv` beside the results: one row per region, with the fraction the worst
inclusion genome holds and how many genomes hold each cut-off. Prints both readings, and if oligos are
given, the regions that carry them.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from primer_finder import blast, seqio  # noqa: E402
from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.pipeline import MOSTLY_PRESENT  # noqa: E402

CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def pattern(oligo: str) -> re.Pattern[str]:
    return re.compile("".join(f"[{CODES.get(base, base)}]" for base in oligo.upper()))


def holds(sequence: str, oligo: str) -> bool:
    """Whether a region carries an oligo, on either strand, IUPAC codes expanded.

    Weaker than `check_assay.py`, deliberately, since this script only needs to mark which rows to read:
    three unordered matches count, with no check that the primers face each other or that the amplicon is
    the published length, and an exact match misses a copy of a primer site that has diverged. Where it
    matters, the regions it marks were checked against `check_assay.py` and agree.
    """
    return bool(pattern(oligo).search(sequence)) or bool(pattern(reverse_complement(oligo)).search(sequence))


def inclusion_of(results: Path, given: Path | None) -> list[Path]:
    """The genomes the run was given, from `run_info.json` unless told otherwise.

    A relative path is resolved against the run folder's parent, which is what the command line that
    produced it meant. Testing `is_dir()` first, as this did, would silently accept a folder of that name
    under whatever directory the script happens to be run from.
    """
    if given is not None:
        return seqio.require_genomes(given, "inclusion")
    info = json.loads((results / "run_info.json").read_text())
    folder = Path(info["parameters"]["inclusion"])
    if not folder.is_absolute():
        folder = (results.parent / folder).resolve()
    return seqio.require_genomes(folder, "inclusion")


def required_of(results: Path, genomes: int) -> int:
    """How many inclusion genomes a region had to be in, as the run recorded it."""
    info = json.loads((results / "run_info.json").read_text())
    recorded = info["parameters"].get("inclusion_genomes_required")
    if recorded:
        return int(recorded)
    return math.ceil(float(info["parameters"].get("min_inclusion", 1.0)) * genomes)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="A finished run's output folder")
    parser.add_argument("-i", "--inclusion", type=Path, help="Default: what the run recorded")
    parser.add_argument("--oligos", help="Comma-separated published oligos, to mark the regions holding them")
    parser.add_argument("-t", "--threads", type=int, default=16)
    args = parser.parse_args(argv)

    final = args.results / "final_kmers.fasta"
    regions = {record.name: record.seq.upper() for record in seqio.iter_records(final)}
    genomes = inclusion_of(args.results, args.inclusion)
    print(f"{args.results}: {len(regions)} region(s), {len(genomes)} inclusion genome(s)")

    with tempfile.TemporaryDirectory(dir=args.results) as work:
        presence = blast.presence_in_genomes(final, genomes, Path(work) / "db", args.threads)

    required = required_of(args.results, len(genomes))
    print(f"  the run required {required} of them to hold a region")
    oligos = [piece.strip() for piece in (args.oligos or "").split(",") if piece.strip()]
    cuts = (0.9, 0.75, 0.5, 0.25)
    rows = []
    for name, sequence in regions.items():
        covered = presence[name]
        held = [value for value in covered.values() if value > 0]
        row = {
            "region": name,
            "length": len(sequence),
            "genomes_holding_any_of_it": len(held),
            # The worst holder, over the genomes that hold any of it; blank where none does. The mean
            # below is over *every* genome, non-holders included, so the two answer different questions
            # and a mean under the worst means some genomes hold none of it.
            "worst_fraction_held": f"{min(held):.3f}" if held else "",
            "mean_fraction_held": f"{statistics.mean(covered.values()):.3f}",
            "carries_the_oligos": "yes" if oligos and all(holds(sequence, o) for o in oligos) else "",
        }
        for cut in cuts:
            row[f"genomes_holding_{cut:g}"] = sum(1 for value in covered.values() if value >= cut)
        rows.append(row)

    out = args.results / "partial_presence.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    worsts = sorted(float(row["worst_fraction_held"]) for row in rows if row["worst_fraction_held"])
    nobody = sum(1 for row in rows if not row["worst_fraction_held"])
    partly = [value for value in worsts if value < MOSTLY_PRESENT]
    print(f"  worst fraction held, over the regions: median {statistics.median(worsts):.3f}, "
          f"lowest {worsts[0]:.3f}" + (f"; {nobody} region(s) held by no genome at all" if nobody else ""))
    print(f"  described: {len(partly)} of {len(rows)} ({100 * len(partly) / len(rows):.1f}%) are held "
          f"only in part (< {MOSTLY_PRESENT:g}) by a genome that holds them at all")
    print("  costed: what a coverage rule would actually drop, applied the way presence is applied "
          f"(fewer than {required} genomes holding that much):")
    for cut in cuts:
        key = f"genomes_holding_{cut:g}"
        would_go = sum(1 for row in rows if row[key] < required)
        print(f"    at {cut:g} of the region: {would_go} of {len(rows)} "
              f"({100 * would_go / len(rows):.1f}%)")
    carrying = [row for row in rows if row["carries_the_oligos"]]
    for row in carrying:
        reach = ", ".join(f"{cut:g}: {row[f'genomes_holding_{cut:g}']}" for cut in cuts)
        print(f"  the published assay is in {row['region']} ({row['length']} bp): worst holder "
              f"{row['worst_fraction_held']}, genomes reaching each cut-off -- {reach}"
              f" (needs {required})")
    if oligos and not carrying:
        print("  no region carries all the given oligos")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
