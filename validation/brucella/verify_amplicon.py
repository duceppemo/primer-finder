#!/usr/bin/env python3
"""Check the published amplicon lengths against a real genome, before trusting them as an answer key.

    python verify_amplicon.py /path/to/work

For each assay, the two primers are located in a genome of the assay's own group and the span between
their outer 5' ends is measured. If that span equals the length the paper reports, the paper counts the
amplicon the same way `check_assay.py` does and the +/-1 tolerance there is unused. If it is off by one,
the tolerance is accounted for. Anything else means the oligos or the lengths in assays.tsv are wrong, and
the run that follows would be checking the wrong thing.

It also reports which contig each primer is on, which is what the paper's two over-fragmented B. suis
assemblies are about.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.seqio import iter_records  # noqa: E402

HERE = Path(__file__).resolve().parent
CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def locate(oligo: str, records: list[tuple[str, str]]) -> list[tuple[str, int, str]]:
    found = []
    for strand, query in (("+", oligo.upper()), ("-", reverse_complement(oligo.upper()))):
        for name, seq in records:
            start = seq.find(query)
            while start != -1:
                found.append((name, start, strand))
                start = seq.find(query, start + 1)
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    args = parser.parse_args(argv)

    panel = read_tsv(args.work / "panel.tsv")
    assays = read_tsv(args.assays)
    status = 0
    for assay in assays:
        wanted = CLASSICAL if assay["target"] == "genus" else (assay["target"],)
        # A complete genome of the assay's own group, so that the amplicon is not split across contigs
        candidates = [row for row in panel if row["species"] in wanted
                      and row["level"] == "Complete Genome"] or \
                     [row for row in panel if row["species"] in wanted]
        measured = None
        for row in sorted(candidates, key=lambda r: r["accession"]):
            records = [(record.name, record.seq.upper()) for record in iter_records(Path(row["path"]))]
            forward = locate(assay["forward"], records)
            reverse = locate(assay["reverse"], records)
            probe = locate(assay["probe"], records)
            if not (forward and reverse and probe):
                continue
            same = {name for name, _, _ in forward} & {name for name, _, _ in reverse}
            if not same:
                print(f"{assay['assay']}: primers on different contigs in {row['accession']}")
                continue
            contig = sorted(same)[0]
            one = next(hit for hit in forward if hit[0] == contig)
            two = next(hit for hit in reverse if hit[0] == contig)
            left = min(one[1], two[1])
            right = max(one[1] + len(assay["forward"]), two[1] + len(assay["reverse"])) - 1
            measured = right - left + 1
            published = int(assay["amplicon_bp"])
            verdict = "as published" if measured == published else \
                      "one base out (counting convention)" if abs(measured - published) == 1 else "WRONG"
            print(f"{assay['assay']:>10s}  {row['accession']} ({row['species']}, {contig}): "
                  f"measured {measured} bp, published {published} bp -- {verdict}")
            if verdict == "WRONG":
                status = 1
            break
        if measured is None:
            print(f"{assay['assay']}: all three oligos found in no genome of its own group", file=sys.stderr)
            status = 1
    return status


if __name__ == "__main__":
    raise SystemExit(main())
