#!/usr/bin/env python3
"""Two small inclusion groups for the IS711 question, cut out of the Brucella panel.

    python validation/repeats/is711_panels.py /path/to/work -o /path/to/is711

`complete/` is 30 complete classical genomes, up to six per species, each holding whole copies of the
element. On those the element is intact, so the only thing in the way of a repeat-targeted marker is
primer-finder's own handling of repeats -- which is what the 2 x 2 of `run_is711.sh` measures.

`draft/` is 15 of those complete genomes and 15 contig-level assemblies, which is what a real panel looks
like: `is711_census.py` finds no whole copy of the element in 557 of the 873 contig-level assemblies of
this panel. It asks the next question -- what the same run does when the repeat is broken in half the
group.

Both get the same exclusion group, 30 complete non-classical genomes: *Ochrobactrum* and its relatives,
which carry no IS711 at all. Choices are seeded, so the folders are the same on every run, and the
genomes are symbolic links into the panel.
"""

from __future__ import annotations

import argparse
import collections
import csv
import random
from pathlib import Path

SEED = 3
PER_SPECIES = 6
GROUP = 30
CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def link_group(rows: list[dict[str, str]], folder: Path) -> None:
    """One folder of symbolic links, named so that a log line says which genome it is."""
    folder.mkdir(parents=True, exist_ok=True)
    for old in folder.iterdir():
        old.unlink()
    for row in rows:
        path = Path(row["path"])
        suffixes = "".join(path.suffixes[-2:]) if path.suffixes[-1] == ".gz" else path.suffix
        (folder / f"{row['accession']}_{row['species']}{suffixes}").symlink_to(path)


def describe(name: str, rows: list[dict[str, str]], copies: dict[str, int]) -> None:
    counts = sorted(copies.get(row["accession"], 0) for row in rows)
    levels = collections.Counter(row["level"] for row in rows)
    print(f"{name}: {len(rows)} genomes, {dict(levels)}")
    print(f"  whole IS711 copies: pooled {sum(counts)}, mean {sum(counts)/len(counts):.2f}, "
          f"median {counts[len(counts)//2]}, range {counts[0]}-{counts[-1]}, "
          f"none at all in {sum(1 for c in counts if c == 0)}")
    print(f"  -d has to reach the mean, so -d {int(sum(counts)/len(counts)) + 1} for this group")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv and is711_census.tsv")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Folder to hold the groups")
    args = parser.parse_args(argv)

    panel = read_tsv(args.work / "panel.tsv")
    census = read_tsv(args.work / "is711_census.tsv")
    copies = {row["accession"]: int(row["element_copies"]) for row in census}
    if not copies:
        raise SystemExit(f"{args.work / 'is711_census.tsv'} is empty: run is711_census.py first.")

    rng = random.Random(SEED)
    complete = sorted((row for row in panel if row["species"] in CLASSICAL
                       and row["level"] == "Complete Genome" and copies.get(row["accession"], 0) > 0),
                      key=lambda row: row["accession"])
    drafts = sorted((row for row in panel if row["species"] in CLASSICAL and row["level"] == "Contig"),
                    key=lambda row: row["accession"])
    by_species = collections.defaultdict(list)
    for row in complete:
        by_species[row["species"]].append(row)

    chosen: list[dict[str, str]] = []
    for _species, rows in sorted(by_species.items()):
        take = min(len(rows), PER_SPECIES)
        chosen += rng.sample(rows, take) if len(rows) > take else rows
    mixed = rng.sample(complete, GROUP // 2) + rng.sample(drafts, GROUP // 2)

    others = sorted((row for row in panel if row["group"] == "non_classical"
                     and row["level"] == "Complete Genome"), key=lambda row: row["accession"])
    exclusion = rng.sample(others, min(GROUP, len(others)))

    for name, rows in (("complete", chosen), ("draft", mixed)):
        link_group(rows, args.output / name / "inclusion")
        link_group(exclusion, args.output / name / "exclusion")
        describe(name, rows, copies)
    print(f"exclusion group (the same for both): {len(exclusion)} non-classical genomes, "
          f"{sum(copies.get(row['accession'], 0) for row in exclusion)} IS711 copies between them")
    print(f"\nfolders under {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
