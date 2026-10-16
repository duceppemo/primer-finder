#!/usr/bin/env python3
"""Design on the regions of a run that hold the published IS711 assay, and print the copy columns.

    python validation/repeats/is711_design.py /path/to/run -o /path/to/design -t 32

The copy columns are the reason this is run at all. `inclusion_copies_min` and `inclusion_copies_max` say
how many times the designed amplicon occurs in the inclusion genome that holds it least and most often,
which for a repeat-targeted assay is the number the assay is chosen for -- and a minimum of 0 says an
inclusion genome does not hold the amplicon at all, which the find stage does not refuse: it counts a
region as present in a genome on any blast hit, so a genome holding a quarter of a region counts.

Only the regions holding all three published oligos are designed on, so the figures are about the
published target and not about whatever else the run reported.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.seqio import iter_records, write_fasta  # noqa: E402

HERE = Path(__file__).resolve().parent
COLUMNS = ("region", "product_size", "inclusion_copies_min", "inclusion_copies_max",
           "exclusion_genomes_with_amplicon", "specific_by")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def regions_holding(assay: dict[str, str], final: Path) -> list:
    """The records of `final` that hold all three oligos, on either strand."""
    oligos = [assay[part].upper() for part in ("forward", "reverse", "probe")]
    keep = []
    for record in iter_records(final):
        sequence = record.seq.upper()
        if all(oligo in sequence or reverse_complement(oligo) in sequence for oligo in oligos):
            keep.append(record)
    return keep


def groups_of(run: Path) -> tuple[Path, Path]:
    """The run's own inclusion and exclusion folders, as absolute paths.

    `run_info.json` records them as they were given on the command line, which for a run started from
    its own folder is relative -- and `design` then looks for them beside wherever it is being run from
    and stops. Resolving them against the run folder's parent is what the original command line meant.
    """
    recorded = json.loads((run / "run_info.json").read_text())["parameters"]
    folders = []
    for name in ("inclusion", "exclusion"):
        folder = Path(recorded[name])
        folders.append(folder if folder.is_absolute() else (run.parent / folder).resolve())
    return folders[0], folders[1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run", type=Path, help="A finished primer-finder run")
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--assays", type=Path, default=HERE / "is711_assay.tsv")
    parser.add_argument("-t", "--threads", type=int, default=16)
    parser.add_argument("--product-size", default="50-150",
                        help="The published amplicon is 63 bp, below the design default of 70-150")
    args = parser.parse_args(argv)

    assay = read_tsv(args.assays)[0]
    holding = regions_holding(assay, args.run / "final_kmers.fasta")
    if not holding:
        print(f"No region of {args.run} holds the three published oligos: nothing to design on.")
        return 1
    print(f"{len(holding)} region(s) hold the published assay: "
          + ", ".join(f"{record.name} ({len(record.seq)} bp)" for record in holding))

    staged = args.output.with_name(args.output.name + "_regions")
    staged.mkdir(parents=True, exist_ok=True)
    (staged / "run_info.json").write_bytes((args.run / "run_info.json").read_bytes())
    write_fasta(staged / "final_kmers.fasta", holding)

    inclusion, exclusion = groups_of(args.run)
    log = args.output.parent / f"{args.output.name}.log"
    with log.open("w") as fh:
        subprocess.run(["primer-finder", "design", str(staged), "-o", str(args.output),
                        "-i", str(inclusion), "-e", str(exclusion),
                        "-t", str(args.threads), "--product-size", args.product_size,
                        "--max-regions", "0"], check=True, stdout=fh, stderr=subprocess.STDOUT)
    rows = read_tsv(args.output / "assays.tsv")
    print(f"{len(rows)} assay(s) designed on them")
    for row in rows:
        print("  " + "  ".join(f"{column}={row[column]}" for column in COLUMNS if column in row))
    least = {row["inclusion_copies_min"] for row in rows}
    most = {row["inclusion_copies_max"] for row in rows}
    print(f"inclusion_copies_min over all of them: {sorted(least)}; max: {sorted(most)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
