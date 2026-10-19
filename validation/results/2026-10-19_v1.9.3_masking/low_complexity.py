#!/usr/bin/env python3
"""What blast's low-complexity masking hides, measured against blast itself.

    python low_complexity.py [work_folder]

For each of nine 200-base shapes, planted in the middle of a 6 kb random genome: whether `blastn` finds
the shape in that genome with its default masking, whether it finds it with `-dust no`, and what
`dustmasker` -- the same SDUST implementation blastn applies to its query -- says is masked. Two cheap
proxies are printed beside them to show why primer-finder asks the masker rather than guessing: the number
of distinct 11-base seeds (blastn's seed length) and a whole-sequence triplet score.

Needs blastn, makeblastdb and dustmasker on PATH. Section 1 of SUMMARY.md is this output.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

PERIOD = 200  # The length of every shape, so that only its composition differs


def shapes(rng: random.Random) -> dict[str, str]:
    dna = lambda n, alphabet="ACGT": "".join(rng.choice(alphabet) for _ in range(n))  # noqa: E731
    return {
        "(AAT) x 67": ("AAT" * 70)[:PERIOD],
        "random AT": dna(PERIOD, "AT"),
        "random AC": dna(PERIOD, "AC"),
        "A x 200": "A" * PERIOD,
        "random ACGT": dna(PERIOD),
        "80% AT, random": "".join(rng.choice("AATT" if rng.random() < 0.8 else "GC")
                                  for _ in range(PERIOD)),
        "(ACGTT) x 40": ("ACGTT" * 40)[:PERIOD],
        "(ACGTTGGCAT) x 20": ("ACGTTGGCAT" * 20)[:PERIOD],
        "half ACGT, half (AAT)": dna(PERIOD // 2) + ("AAT" * 40)[:PERIOD // 2],
    }


def triplet_score(sequence: str) -> float:
    """A whole-sequence version of the DUST score: how far the triplet composition is from even."""
    counts: dict[str, int] = {}
    for index in range(len(sequence) - 2):
        triplet = sequence[index:index + 3]
        counts[triplet] = counts.get(triplet, 0) + 1
    return sum(count * (count - 1) / 2 for count in counts.values()) / (len(sequence) - 3)


def main(argv: list[str]) -> int:
    folder = Path(argv[0] if argv else "low_complexity_work").resolve()
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(31)
    elements = shapes(rng)
    left, right = ("".join(rng.choice("ACGT") for _ in range(3000)) for _ in range(2))

    query = folder / "shapes.fasta"
    query.write_text("".join(f">s{n}\n{element}\n" for n, element in enumerate(elements.values())))
    masked = folder / "masked.txt"
    subprocess.run(["dustmasker", "-in", str(query), "-outfmt", "acclist", "-out", str(masked)],
                   check=True, capture_output=True)
    takes_out: dict[str, list[tuple[int, int]]] = {}
    for line in masked.read_text().splitlines():
        name, start, end = line.split("\t")
        takes_out.setdefault(name[1:], []).append((int(start), int(end)))

    print(f"{'shape':24} {'hits, masking on':>16} {'hits, -dust no':>15} {'dustmasker':>22} "
          f"{'distinct 11-mers':>17} {'triplet score':>14}")
    for n, (name, element) in enumerate(elements.items()):
        genome = folder / f"g{n}.fasta"
        genome.write_text(f">chr\n{left}{element}{right}\n")
        subprocess.run(["makeblastdb", "-in", str(genome), "-dbtype", "nucl",
                        "-out", str(folder / f"db{n}")], check=True, capture_output=True)
        one = folder / f"q{n}.fasta"
        one.write_text(f">q\n{element}\n")

        def hits(dust: bool, one: Path = one, n: int = n) -> int:
            done = subprocess.run(["blastn", "-task", "blastn", "-query", str(one),
                                   "-db", str(folder / f"db{n}"), "-evalue", "1e-10",
                                   "-outfmt", "6 qseqid length pident",
                                   *([] if dust else ["-dust", "no", "-soft_masking", "false"])],
                                  check=True, capture_output=True, text=True)
            return len(done.stdout.strip().splitlines()) if done.stdout.strip() else 0

        says = takes_out.get(f"s{n}", [])
        seeds = len({element[at:at + 11] for at in range(len(element) - 10)})
        print(f"{name:24} {hits(True):>16} {hits(False):>15} "
              f"{str(says) if says else 'nothing':>22} {seeds:>17} {triplet_score(element):>14.1f}")
    print(f"\nEvery shape is {PERIOD} bases, planted in the middle of the same 6 kb random genome, and is "
          f"held by it exactly.\nOf {len(elements)} shapes, dustmasker and blastn's own masking agree on "
          f"every one.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
