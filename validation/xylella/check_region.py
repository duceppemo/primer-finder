#!/usr/bin/env python3
"""Where a recovered region sits in the reference genome, next to the coordinates the paper reports.

This is context, not a verdict. The region primer-finder reports is not the same object as the specific
long-mer the paper reports: the exclusion genomes are not the same, and a long-mer is bounded by the
exclusion set rather than by an assembly. Whether the assay is inside the region is what check_assay.py
decides; this script only says how far the two regions overlap.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.seqio import read_fasta, write_fasta  # noqa: E402

HERE = Path(__file__).resolve().parent


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="The output folder of a primer-finder run.")
    parser.add_argument("genomes", type=Path, help="The per-subspecies genome folder.")
    parser.add_argument("--target", required=True)
    parser.add_argument("--region", type=Path, default=HERE / "regions.tsv")
    parser.add_argument("--check", type=Path, help="The JSON check_assay.py wrote, naming the region.")
    args = parser.parse_args(argv)

    expected = {row["target"]: row for row in read_tsv(args.region)}[args.target]
    name = json.loads(args.check.read_text())["region"] if args.check else None
    regions = read_fasta(args.results / "final_kmers.fasta")
    if name not in regions:
        print(f"{name} is not in the results", file=sys.stderr)
        return 1

    reference = next(args.genomes.rglob(f"{expected['reference_accession']}*.fasta"), None)
    if reference is None:
        print(f"No genome file for {expected['reference_accession']}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory() as work:
        work = Path(work)
        write_fasta(work / "region.fasta", [regions[name]])
        subprocess.run(["makeblastdb", "-in", reference.resolve(), "-dbtype", "nucl",
                        "-out", work / "ref"], check=True, capture_output=True)
        out = subprocess.run(
            ["blastn", "-db", work / "ref", "-query", work / "region.fasta", "-outfmt",
             "6 sseqid sstart send length pident", "-max_target_seqs", "5"],
            check=True, capture_output=True, text=True).stdout

    rows = [line.split("\t") for line in out.splitlines()]
    if not rows:
        print("The region does not align to the reference", file=sys.stderr)
        return 1
    contig, sstart, send, length, pident = rows[0]
    start, end = sorted((int(sstart), int(send)))
    want_start, want_end = int(expected["region_start"]), int(expected["region_end"])
    overlap = max(0, min(end, want_end) - max(start, want_start) + 1)
    covers = start <= want_start and end >= want_end
    print(f"{expected['assay']}: {name} ({len(regions[name].seq)} bp) aligns to {reference.name} "
          f"{contig}:{start}-{end} ({length} bp, {pident}% identity)")
    print(f"  the paper reports {want_start}-{want_end} ({expected['region_size']} bp) in "
          f"{expected['reference']}")
    print(f"  overlap with the published region: {overlap} bp of {expected['region_size']} "
          f"({100 * overlap // int(expected['region_size'])}%), "
          f"fully covered: {covers}")
    # Whether the two regions coincide is context, not a verdict: the exclusion sets differ, and a region
    # of specific kmers assembled here is not the same object as the long-mer the paper reports. What the
    # validation asks is whether the assay itself is inside the region, which check_assay.py decides.
    return 0 if overlap > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
