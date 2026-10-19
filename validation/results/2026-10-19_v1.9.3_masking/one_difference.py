#!/usr/bin/env python3
"""A low-complexity region that every genome holds, and what masking did to it. End to end.

    python one_difference.py [work_folder] [threads] [memory_GB]

Three inclusion and three exclusion genomes, 6 kb of random sequence either side of a 200-base AT element
that both groups carry. One base of the element differs between the groups, which makes the kmers over it
inclusion-specific: `find` assembles a 101-base region of pure low-complexity sequence, and the exclusion
genomes hold 100 of its 101 bases.

The run is done three ways -- with no unmasked second look, with one in the presence search only, and as
the code now stands -- and each is printed with what it said about the region.
"""

from __future__ import annotations

import json
import random
import shutil
import sys
from pathlib import Path

from primer_finder import blast, pipeline

ELEMENT = ("AAT" * 70)[:200]


def build(folder: Path) -> Path:
    rng = random.Random(23)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    left, right = dna(6000), dna(6000)
    changed = ELEMENT[:100] + ("C" if ELEMENT[100] != "C" else "G") + ELEMENT[101:]
    for group, element in (("inclusion", changed), ("exclusion", ELEMENT)):
        (folder / group).mkdir(parents=True, exist_ok=True)
        for number in range(3):
            (folder / group / f"{group[:3]}_{number}.fasta").write_text(f">chr\n{left}{element}{right}\n")
    return folder


def describe(output: Path) -> str:
    info_file = output / "run_info.json"
    if not info_file.is_file():
        return "nothing reported"
    info = json.loads(info_file.read_text())
    regions = {}
    for block in (output / "final_kmers.fasta").read_text().split(">")[1:]:
        head, *rest = block.splitlines()
        regions[head.split()[0]] = len("".join(rest))
    absent = info.get("absent_stretches") or {}
    return (f"{len(regions)} region(s) {sorted(regions.values())}, "
            f"absent from every exclusion genome: {absent or 'nothing'}, "
            f"last_resort={info.get('last_resort')}")


def main(argv: list[str]) -> int:
    work = Path(argv[0] if argv else "one_difference_work").resolve()
    threads = int(argv[1]) if len(argv) > 1 else 4
    memory = int(argv[2]) if len(argv) > 2 else 8
    build(work)
    real = blast.unmasked_hits

    def presence_only(db, folder, hits, sizes, sequences, masked, fields=blast.ALIGNMENT_FIELDS):
        if tuple(fields) == tuple(blast.ALIGNMENT_FIELDS):
            return []  # What the exclusion search did before
        return real(db, folder, hits, sizes, sequences, masked, fields)

    ways = (("no second look at all (every version up to 1.9.2)", lambda *_a, **_k: []),
            ("a second look in the presence search only", presence_only),
            ("as it stands", real))
    for label, patched in ways:
        blast.unmasked_hits = patched
        output = work / f"out_{ways.index((label, patched))}"
        shutil.rmtree(output, ignore_errors=True)
        try:
            pipeline.run(pipeline.Settings(inclusion=work / "inclusion", exclusion=work / "exclusion",
                                           output=output, threads=threads, memory_gb=memory,
                                           kmer_size=51))
        except Exception as error:  # noqa: BLE001 - "nothing to report" is one of the answers
            print(f"{label}: stopped -- {str(error).splitlines()[0][:120]}")
            continue
        print(f"{label}: {describe(output)}")
    blast.unmasked_hits = real
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
