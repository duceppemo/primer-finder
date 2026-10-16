#!/usr/bin/env python3
"""A dataset whose inclusion-specific marker is present in several copies per genome.

The shape a repeat-targeted assay has: a marker carried two to four times by each inclusion genome and
not at all by the exclusion genomes. Copy number varies between genomes, as a transposon's would, because
that is what decides whether `-d` keeps the marker's kmers.
"""
from __future__ import annotations

import gzip
import json
import random
import sys
from pathlib import Path

LENGTH = 20_000
MARKER = 300          # The repeated, inclusion-specific stretch
COPIES = (2, 3, 4, 3)  # Per inclusion genome
PRIVATE = 20
SEED = 42
GROUP = 4


def main(out: Path) -> int:
    rng = random.Random(SEED)
    pick = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    backbone = pick(LENGTH)
    marker = pick(MARKER)
    truth = {"length": LENGTH, "marker_length": MARKER, "marker": marker, "copies": list(COPIES)}

    for group in ("inclusion", "exclusion"):
        folder = out / group
        folder.mkdir(parents=True, exist_ok=True)
        for number in range(1, GROUP + 1):
            bases = list(backbone)
            if group == "inclusion":
                # Spread the copies out so no two sit within a kmer of each other
                for copy in range(COPIES[number - 1]):
                    at = 2_000 + copy * 4_000
                    bases[at:at] = list(marker)
            for _ in range(PRIVATE):
                where = rng.randrange(len(bases))
                bases[where] = rng.choice([b for b in "ACGT" if b != bases[where]])
            name = f"{group}_{number}.fasta" + (".gz" if number == GROUP else "")
            text = f">chr\n" + "\n".join("".join(bases)[i:i + 70]
                                         for i in range(0, len(bases), 70)) + "\n"
            path = folder / name
            if name.endswith(".gz"):
                path.write_bytes(gzip.compress(text.encode()))
            else:
                path.write_text(text)
    (out / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(f"marker {MARKER} bp, copies per inclusion genome {COPIES}, written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1])))
