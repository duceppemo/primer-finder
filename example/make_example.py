#!/usr/bin/env python3
"""Write the example dataset: four inclusion and four exclusion genomes with a known specific region.

Every genome is the same 20 kb backbone, with its own scattered single-base changes. The inclusion genomes
all carry the same planted variants around position 9,000: three mismatches within 16 bases, a 15-base
insertion, and one isolated mismatch. Those are the only differences the inclusion group shares and the
exclusion group lacks, so primer-finder must report that region and nothing else.
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
from pathlib import Path

LENGTH = 20_000  # The backbone, in bases
SEED = 42
CHROMOSOME = "chr"
TARGET = 9_000  # Where the planted variants start
# The planted variants, as offsets from TARGET: three mismatches within 16 bases, an insertion, and a
# mismatch far enough away to be useless on its own.
MISMATCHES = (0, 8, 15, 100)
INSERTION_AT = 50
INSERTION = "GGCATTACGTTAACC"  # 15 bases
PRIVATE_CHANGES = 20  # Per genome, away from the planted region
KEEP_CLEAR = 500  # No private change this close to the planted region
GROUP_SIZE = 4


def random_sequence(rng: random.Random, length: int) -> list[str]:
    return [rng.choice("ACGT") for _ in range(length)]


def other_base(rng: random.Random, base: str) -> str:
    return rng.choice([other for other in "ACGT" if other != base])


def add_private_changes(rng: random.Random, bases: list[str], how_many: int) -> list[int]:
    """Change single bases at random, away from the planted region, and return where."""
    positions: list[int] = []
    while len(positions) < how_many:
        position = rng.randrange(100, LENGTH - 100)
        if abs(position - TARGET) < KEEP_CLEAR or position in positions:
            continue
        bases[position] = other_base(rng, bases[position])
        positions.append(position)
    return sorted(positions)


def plant_variants(rng: random.Random, bases: list[str]) -> dict[str, object]:
    """Apply the inclusion-specific variants to a copy of the backbone.

    `mismatches` and `insertion_at` are positions in this sequence (the inclusion genomes), and `variants`
    is every base that differs from the exclusion genomes: the mismatches and the inserted bases. That is
    exactly what primer-finder must report in lower case, so check_example.py compares against it.
    """
    insertion_at = TARGET + INSERTION_AT
    mismatches = []
    for offset in MISMATCHES:
        position = TARGET + offset
        bases[position] = other_base(rng, bases[position])
        # A mismatch after the insertion point moves along once the insertion is spliced in
        mismatches.append(position + len(INSERTION) if position >= insertion_at else position)
    bases[insertion_at:insertion_at] = list(INSERTION)
    variants = sorted(set(mismatches) | set(range(insertion_at, insertion_at + len(INSERTION))))
    return {"mismatches": sorted(mismatches), "insertion_at": insertion_at, "insertion": INSERTION,
            "variants": variants}


def write_genome(path: Path, bases: list[str], line_length: int = 70) -> None:
    sequence = "".join(bases)
    lines = [sequence[start:start + line_length] for start in range(0, len(sequence), line_length)]
    text = f">{CHROMOSOME}\n" + "\n".join(lines) + "\n"
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", nargs="?", type=Path, default=Path(__file__).parent / "data",
                        help="Folder to write the inclusion/ and exclusion/ folders to.")
    args = parser.parse_args(argv)

    rng = random.Random(SEED)
    backbone = random_sequence(rng, LENGTH)
    inclusion_backbone = list(backbone)
    planted = plant_variants(rng, inclusion_backbone)

    truth: dict[str, object] = {"length": LENGTH, "chromosome": CHROMOSOME, "target": TARGET, **planted}
    for group, template in (("inclusion", inclusion_backbone), ("exclusion", backbone)):
        folder = args.output / group
        folder.mkdir(parents=True, exist_ok=True)
        for number in range(1, GROUP_SIZE + 1):
            bases = list(template)
            changes = add_private_changes(rng, bases, PRIVATE_CHANGES)
            # One genome of each group is gzipped: primer-finder has to handle both.
            name = f"{group}_{number}.fasta" + (".gz" if number == GROUP_SIZE else "")
            write_genome(folder / name, bases)
            truth.setdefault(f"{group}_private_changes", {})[name] = changes  # type: ignore[union-attr]
    (args.output / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(f"Wrote {GROUP_SIZE} inclusion and {GROUP_SIZE} exclusion genomes to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
