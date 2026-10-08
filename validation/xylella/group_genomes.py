#!/usr/bin/env python3
"""Build the inclusion and exclusion folders of one target out of the per-subspecies folders.

The folders hold symbolic links, which primer-finder follows, so the genomes are downloaded only once.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("genomes", type=Path, help="The folder get_genomes.py filled (one per subspecies).")
    parser.add_argument("output", type=Path, help="Folder to create inclusion/ and exclusion/ in.")
    parser.add_argument("--target", required=True, help="The subspecies that makes up the inclusion group.")
    args = parser.parse_args(argv)

    subspecies = sorted(folder.name for folder in args.genomes.iterdir() if folder.is_dir())
    if args.target not in subspecies:
        parser.error(f"No genomes for {args.target}; have: {', '.join(subspecies)}")
    counts = {}
    for group, names in (("inclusion", [args.target]),
                         ("exclusion", [name for name in subspecies if name != args.target])):
        folder = args.output / group
        shutil.rmtree(folder, ignore_errors=True)
        folder.mkdir(parents=True)
        for name in names:
            for genome in sorted((args.genomes / name).glob("*.fasta")):
                (folder / genome.name).symlink_to(genome.resolve())
        counts[group] = len(list(folder.glob("*.fasta")))
    print(f"{args.target}: {counts['inclusion']} inclusion, {counts['exclusion']} exclusion "
          f"({', '.join(name for name in subspecies if name != args.target)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
