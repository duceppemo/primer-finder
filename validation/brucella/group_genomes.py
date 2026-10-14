#!/usr/bin/env python3
"""Build the inclusion and exclusion folders for one Brucella target, as symbolic links.

    python group_genomes.py /path/to/work/panel.tsv /path/to/work/suis --target suis

Four targets, the paper's four assays:

- `genus`: the inclusion group is the eight classical species the paper's panel holds, the exclusion group
  every non-classical genome (the former Ochrobactrum, now filed under Brucella, plus Pseudochrobactrum,
  Paenochrobactrum and Falsochrobactrum);
- `genus_wide`: the same, plus the other four species of the classical clade -- B. pinnipedialis,
  B. amazoniensis, B. vulpis and B. inopinata. A genus assay ought to detect those too, and the paper's
  panel had none of them, so this asks what spanning the whole clade costs;
- `abortus`, `melitensis`, `suis`: the inclusion group is that species, and the exclusion group is **the
  other classical species as well as** the non-classical genomes. A species assay has to tell its species
  from its close relatives first of all, which is the harder half of the job.

A genome labelled `held_out` is in neither group, whichever the target: "Brucella sp." is unidentified, and
about half of those assemblies are former Ochrobactrum. `census.py` reports them separately. The four
`classical_extra` species are in neither group except for `genus_wide`, where they are the point.

`--max-exclusion` keeps that many exclusion genomes, spread evenly over the species present, for a quicker
run. The default is all of them.

`--hold-out` takes a file of accessions to put in neither group, one per line. It is there for the
assemblies whose species label the census shows cannot be right -- a genome that carries exactly one
species assay, and not the one its label names. Leaving such a genome in the exclusion group removes a
region that is specific to the group it actually belongs to, so a run on the panel as labelled and a run
with those genomes held out answer two different questions, and the record reports both.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path

CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")
CLASSICAL_EXTRA = ("pinnipedialis", "amazoniensis", "vulpis", "inopinata")
TARGETS = ("genus", "genus_wide", "abortus", "melitensis", "suis")


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def split(rows: list[dict[str, str]], target: str) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    classical = [row for row in rows if row["species"] in CLASSICAL]
    extra = [row for row in rows if row["species"] in CLASSICAL_EXTRA]
    negative = [row for row in rows if row["group"] == "non_classical"]
    if target == "genus":
        return classical, negative
    if target == "genus_wide":
        return classical + extra, negative
    inclusion = [row for row in classical if row["species"] == target]
    exclusion = [row for row in classical if row["species"] != target] + negative
    return inclusion, exclusion


def spread(rows: list[dict[str, str]], how_many: int) -> list[dict[str, str]]:
    """Up to `how_many` rows, round-robin over the species, so no single one fills the group."""
    by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)
    wheels = [iter(sorted(group, key=lambda r: r["accession"])) for _, group in sorted(by_species.items())]
    taken: list[dict[str, str]] = []
    for wheel in itertools.cycle(wheels):
        if len(taken) >= how_many or not wheels:
            break
        row = next(wheel, None)
        if row is None:
            wheels = [one for one in wheels if one is not wheel]
            continue
        taken.append(row)
    return taken


def link(rows: list[dict[str, str]], folder: Path) -> int:
    folder.mkdir(parents=True, exist_ok=True)
    for stale in folder.glob("*.fasta"):
        stale.unlink()
    for row in rows:
        name = folder / f"{row['accession']}_{row['species']}.fasta"
        name.symlink_to(Path(row["path"]).resolve())
    return len(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("panel", type=Path, help="panel.tsv from panel.py")
    parser.add_argument("output", type=Path, help="Folder to hold inclusion/ and exclusion/")
    parser.add_argument("--target", required=True, choices=TARGETS)
    parser.add_argument("--max-exclusion", type=int, default=0, help="0 for every exclusion genome")
    parser.add_argument("--hold-out", type=Path,
                        help="File of accessions to put in neither group, one per line")
    args = parser.parse_args(argv)

    rows = read(args.panel)
    held_back = set()
    if args.hold_out:
        held_back = {line.split("#")[0].strip() for line in args.hold_out.read_text().splitlines()}
        held_back.discard("")
        known = {row["accession"] for row in rows}
        unknown = sorted(held_back - known)
        if unknown:
            raise SystemExit(f"--hold-out names accessions that are not in the panel: {', '.join(unknown)}")
        rows = [row for row in rows if row["accession"] not in held_back]
    inclusion, exclusion = split(rows, args.target)
    if not inclusion:
        raise SystemExit(f"no inclusion genome for {args.target}")
    if args.max_exclusion:
        exclusion = spread(exclusion, args.max_exclusion)

    link(inclusion, args.output / "inclusion")
    link(exclusion, args.output / "exclusion")
    held = sum(1 for row in rows if row["group"] == "held_out")
    species = len({row["species"] for row in exclusion})
    extra = f" (+{len(held_back)} held back by --hold-out)" if held_back else ""
    print(f"{args.target}: {len(inclusion)} inclusion, {len(exclusion)} exclusion from {species} other "
          f"species/group(s); {held} genome(s) held out{extra}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
