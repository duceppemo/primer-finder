#!/usr/bin/env python3
"""Check that primer-finder reported a region holding a published clonal-complex assay.

For the target given, the primers and the probe of assays.tsv are looked for in the candidate regions of a
run (final_kmers.fasta). An assay is recovered when all three oligos are found in the same region, with the
primers facing each other, the probe between them, and the amplicon the published length.

Two differences from the Xylella check:

- it matches IUPAC codes, because one of the published primers carries an M;
- the amplicon is measured between the outer 5' ends of the two primers, inclusive, whichever of them the
  region happens to hold first. Felix et al. report a length exactly one base shorter than that for every
  assay checked here, so a span of `published` or `published + 1` counts as the published layout; the
  one-base gap is a convention in the paper, not a difference in the sequence.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.seqio import iter_records  # noqa: E402

HERE = Path(__file__).resolve().parent

# What each IUPAC code stands for, so that a published oligo written with one is still found
CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def pattern(oligo: str) -> re.Pattern[str]:
    out = []
    for base in oligo.upper():
        if base not in CODES:
            raise SystemExit(f"{oligo}: {base} is not an IUPAC code")
        bases = CODES[base]
        out.append(bases if len(bases) == 1 else f"[{bases}]")
    return re.compile("".join(out))


@dataclass
class Hit:
    region: str
    start: int
    strand: str
    region_length: int


def find(oligo: str, records: list[tuple[str, str]]) -> list[Hit]:
    """Every occurrence of an oligo in the regions, on either strand, IUPAC codes included."""
    hits = []
    for strand, query in (("+", oligo), ("-", reverse_complement(oligo))):
        search = pattern(query)
        for name, seq in records:
            for match in search.finditer(seq):
                hits.append(Hit(name, match.start(), strand, len(seq)))
    return hits


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="The output folder of a primer-finder run.")
    parser.add_argument("--target", required=True, help="Which assay of assays.tsv to look for.")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("--json", type=Path, help="Also write the findings as JSON.")
    args = parser.parse_args(argv)

    assays = {row["target"]: row for row in read_tsv(args.assays)}
    if args.target not in assays:
        parser.error(f"No assay for {args.target}; have: {', '.join(assays)}")
    assay = assays[args.target]

    final = args.results / "final_kmers.fasta"
    records = [(record.name, record.seq.upper()) for record in iter_records(final)]
    ranks = {name: position for position, (name, _) in enumerate(records, start=1)}
    print(f"{assay['assay']}: {len(records)} candidate region(s) in {final}")

    hits = {part: find(assay[part], records) for part in ("forward", "reverse", "probe")}
    for part, found in hits.items():
        where = ", ".join(f"{hit.region} at {hit.start} ({hit.strand})" for hit in found) or "not found"
        print(f"  {assay['assay']}-{part[0].upper()} {assay[part]:>32s}  {where}")

    shared = set.intersection(*({hit.region for hit in found} for found in hits.values())) \
        if all(hits.values()) else set()
    outcome: dict[str, object] = {
        "assay": assay["assay"], "target": args.target, "regions": len(records),
        "recovered": bool(shared), "region": None, "amplicon_bp": None,
        "expected_amplicon_bp": int(assay["amplicon_bp"]), "rank": None, "region_length": None,
    }
    if not shared:
        print("\nNOT RECOVERED: the three oligos are not together in one region")
    for region in sorted(shared):
        forward = next(hit for hit in hits["forward"] if hit.region == region)
        reverse = next(hit for hit in hits["reverse"] if hit.region == region)
        # The amplicon runs from the outer 5' end of one primer to the outer 5' end of the other, whichever
        # way round the region was assembled
        left = min(forward.start, reverse.start)
        right = max(forward.start + len(assay["forward"]), reverse.start + len(assay["reverse"])) - 1
        amplicon = right - left + 1
        published = int(assay["amplicon_bp"])
        facing = forward.strand != reverse.strand
        probe_inside = any(left <= hit.start <= right for hit in hits["probe"] if hit.region == region)
        as_published = amplicon in (published, published + 1)
        outcome.update({"region": region, "amplicon_bp": amplicon, "rank": ranks[region],
                        "region_length": forward.region_length})
        print(f"\nRECOVERED in {region} ({forward.region_length} bp, rank {ranks[region]} of {len(records)})")
        print(f"  amplicon {amplicon} bp (the paper reports {published} bp), "
              f"primers facing each other: {facing}, probe between them: {probe_inside}")
        if not (as_published and facing and probe_inside):
            outcome["recovered"] = False
            print("  but the layout does not match the published assay")
    if args.json:
        args.json.write_text(json.dumps(outcome, indent=2) + "\n")
    return 0 if outcome["recovered"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
