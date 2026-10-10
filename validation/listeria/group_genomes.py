#!/usr/bin/env python3
"""Build the inclusion and exclusion folders for one clonal complex, as symbolic links.

    python group_genomes.py /path/to/work/genomes.tsv /path/to/work/CC9 --target CC9

The inclusion group is every genome of that clonal complex; the exclusion group is every genome of a
*different* one. A genome the scheme cannot place in any clonal complex goes in neither: it might be the
target, and putting it in the exclusion group would throw away real regions.

`--max-exclusion` keeps the first N exclusion genomes, spread evenly over the clonal complexes present, for
a quicker run. The default is all of them.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def spread(rows: list[dict[str, str]], how_many: int) -> list[dict[str, str]]:
    """Up to `how_many` rows, taking them round-robin over the clonal complexes so that no single one
    fills the group."""
    by_cc: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_cc[row["cc"]].append(row)
    wheels = [iter(sorted(group, key=lambda r: r["accession"])) for _, group in sorted(by_cc.items())]
    taken: list[dict[str, str]] = []
    for wheel in itertools.cycle(wheels):
        if len(taken) >= how_many or not wheels:
            break
        row = next(wheel, None)
        if row is None:
            wheels = [w for w in wheels if w is not wheel]
            continue
        taken.append(row)
    return taken


def link(rows: list[dict[str, str]], folder: Path) -> int:
    folder.mkdir(parents=True, exist_ok=True)
    for row in rows:
        target = Path(row["path"]).resolve()
        name = folder / f"{row['accession']}_{row['cc']}.fasta"
        if name.is_symlink() or name.exists():
            name.unlink()
        name.symlink_to(target)
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("genomes", type=Path, help="genomes.tsv from type_genomes.py")
    parser.add_argument("output", type=Path, help="Folder to hold inclusion/ and exclusion/")
    parser.add_argument("--target", required=True, help="The clonal complex to put in the inclusion group")
    parser.add_argument("--max-exclusion", type=int, default=0, help="0 for every exclusion genome")
    args = parser.parse_args(argv)

    rows = read(args.genomes)
    placed = [row for row in rows if row["cc"]]
    inclusion = [row for row in placed if row["cc"] == args.target]
    exclusion = [row for row in placed if row["cc"] != args.target]
    if not inclusion:
        have = sorted({row["cc"] for row in placed})
        raise SystemExit(f"no genome of {args.target}; have: {', '.join(have)}")
    if args.max_exclusion:
        exclusion = spread(exclusion, args.max_exclusion)

    link(inclusion, args.output / "inclusion")
    link(exclusion, args.output / "exclusion")
    others = len({row["cc"] for row in exclusion})
    unplaced = len(rows) - len(placed)
    print(f"{args.target}: {len(inclusion)} inclusion, {len(exclusion)} exclusion from {others} other "
          f"clonal complex(es); {unplaced} genome(s) left out for having none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
