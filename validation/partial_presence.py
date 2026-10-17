#!/usr/bin/env python3
"""How much of each reported region every inclusion genome actually holds.

    python validation/partial_presence.py /path/to/results -t 32
    python validation/partial_presence.py /path/to/results --oligos ACGT...,TTGC...,GGAC...

A region is reported when every inclusion genome holds it, and "holds it" is decided on any blast hit at
all: no fraction of the region is required. This measures what that lets through on a finished run, which
is the evidence for deciding whether the rule should ask for a fraction -- and, with `--oligos`, whether
such a rule would have cost the run a published assay it did recover.

Writes `partial_presence.tsv` beside the results: one row per region, with the fraction the worst
inclusion genome holds. Prints the distribution and, if oligos are given, the regions that carry them.
"""

from __future__ import annotations

import argparse
import csv
import json
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
    return bool(pattern(oligo).search(sequence)) or bool(pattern(reverse_complement(oligo)).search(sequence))


def inclusion_of(results: Path, given: Path | None) -> list[Path]:
    """The genomes the run was given, from `run_info.json` unless told otherwise."""
    if given is not None:
        return seqio.require_genomes(given, "inclusion")
    info = json.loads((results / "run_info.json").read_text())
    folder = Path(info["parameters"]["inclusion"])
    if not folder.is_dir():  # A run started from its own folder recorded a relative path
        folder = (results.parent / folder).resolve()
    return seqio.require_genomes(folder, "inclusion")


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

    oligos = [piece.strip() for piece in (args.oligos or "").split(",") if piece.strip()]
    rows = []
    for name, sequence in regions.items():
        held = [covered for covered in presence[name].values() if covered > 0]
        worst = min(held) if held else 0.0
        rows.append({
            "region": name,
            "length": len(sequence),
            "genomes_holding_any_of_it": len(held),
            "worst_fraction_held": f"{worst:.3f}",
            "mean_fraction_held": f"{statistics.mean(presence[name].values()):.3f}",
            "carries_the_oligos": "yes" if oligos and all(holds(sequence, o) for o in oligos) else "",
        })

    out = args.results / "partial_presence.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    worsts = sorted(float(row["worst_fraction_held"]) for row in rows)
    partly = [value for value in worsts if value < MOSTLY_PRESENT]
    print(f"  worst fraction held, over the regions: median {statistics.median(worsts):.3f}, "
          f"lowest {worsts[0]:.3f}")
    print(f"  only partly in at least one genome (< {MOSTLY_PRESENT:g}): {len(partly)} of {len(rows)} "
          f"({100 * len(partly) / len(rows):.1f}%)")
    for cut in (0.9, 0.75, 0.5, 0.25):
        would_go = sum(1 for value in worsts if value < cut)
        print(f"    a rule asking for {cut:g} of the region would drop {would_go} of {len(rows)} "
              f"({100 * would_go / len(rows):.1f}%)")
    carrying = [row for row in rows if row["carries_the_oligos"]]
    for row in carrying:
        print(f"  the published assay is in {row['region']} ({row['length']} bp), worst genome holds "
              f"{row['worst_fraction_held']} of it")
    if oligos and not carrying:
        print("  no region carries all the given oligos")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
