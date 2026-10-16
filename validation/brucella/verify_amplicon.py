#!/usr/bin/env python3
"""Check the published amplicon lengths against a real genome, before trusting them as an answer key.

    python verify_amplicon.py /path/to/work

For each assay, the two primers are located in a genome of the assay's own group and the span between
the outer ends of the closest inward-facing pair of them is measured. If that span equals the length the
paper reports, the paper counts the amplicon the same way `check_assay.py` does and the +/-1 tolerance
there is unused. If it is off by one, the tolerance is accounted for. Anything else means the oligos or
the lengths in assays.tsv are wrong, and the run that follows would be checking the wrong thing.

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
sys.path.insert(0, str(HERE))

from groups import inclusion_species  # noqa: E402


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


def facing_pair(
    forward: list[tuple[str, int, str]],
    reverse: list[tuple[str, int, str]],
    contigs: set[str],
) -> tuple[str, tuple[str, int, str], tuple[str, int, str]] | None:
    """The closest placement of the two primers that would amplify: one on each strand, with the
    plus-strand one to the left of the minus-strand one.

    Both primers being on one contig is not enough. A back-to-back pair -- the minus-strand primer to the
    left -- points outwards and amplifies nothing, so measuring the span between it would invent an
    amplicon and report it as the published length. Each primer can also occur more than once, and taking
    the first of each gives the span between two unrelated copies. This is the defect `check_assay.py` had
    and now checks for.
    """
    best = None
    for contig in sorted(contigs):
        for one in (hit for hit in forward if hit[0] == contig):
            for two in (hit for hit in reverse if hit[0] == contig):
                if one[2] == two[2]:  # both on the same strand: not a pair
                    continue
                plus, minus = (one, two) if one[2] == "+" else (two, one)
                if plus[1] >= minus[1]:  # back to back
                    continue
                span = minus[1] - plus[1]
                if best is None or span < best[0]:
                    best = (span, contig, one, two)
    return None if best is None else (best[1], best[2], best[3])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    args = parser.parse_args(argv)

    panel = read_tsv(args.work / "panel.tsv")
    assays = read_tsv(args.assays)
    status = 0
    for assay in assays:
        wanted = inclusion_species(assay["target"])
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
            pair = facing_pair(forward, reverse, same)
            if pair is None:
                print(f"{assay['assay']}: the primers share {sorted(same)[0]} in {row['accession']} but "
                      f"no placement of them faces inwards")
                continue
            contig, one, two = pair
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
