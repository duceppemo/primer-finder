#!/usr/bin/env python3
"""What insilicoPCR 0.8.0's hit-filtering knobs actually govern, measured.

The successor of ../2026-10-08_v1.3.0/mismatch_zones.py, which measured 0.6.1. 0.8.0 replaced the
alignment's accidental zones with a deliberate model -- unaligned 3' bases count against `-m`, a weighted
check guards the last `--3p-window` bases, and up to `--5p-truncation` unaligned 5' bases are free -- and
this measures each claim rather than quoting it. A 22-mer primer pair is placed in a synthetic template and
the *template* is mutated: one mismatch at a known distance from a primer end, or a run of them ending at
the 3' end. Each case is its own genome, so every row of the report names the case it came from.

    python mismatch_zones.py /path/to/insilicoPCR-linux-x64 [work_folder]

Prints one table per question. Section 1 of SUMMARY.md is this output.
"""

from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

TOLERANCES = (0, 1, 2, 3, 10)
OFFSETS = range(1, 9)   # one mismatch, this far from the left primer's 3' end (1 = the terminal base)
RUNS = (2, 3, 4)        # a run of this many mismatches, ending at the 3' end
PENALTIES = (0, 1, 2, 3, 4)
FLIP = {"A": "G", "C": "T", "G": "A", "T": "C"}

rng = random.Random(3)
random_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
reverse_complement = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))  # noqa: E731


def assay(folder: Path, name: str, left: str, cases: dict[str, list[int]], at_five: bool = False) -> str:
    """One primer pair and one genome per case, the template mutated at the given offsets from the left
    primer's 3' end -- or from its 5' end with `at_five`. Returns the primer fasta records."""
    right = random_seq(22)
    amplicon = left + random_seq(80) + reverse_complement(right)
    genomes = folder / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    for case, offsets in cases.items():
        bases = list(amplicon)
        for offset in offsets:
            at = (offset - 1) if at_five else (len(left) - offset)
            bases[at] = FLIP[bases[at]]
        (genomes / f"{case}.fasta").write_text(
            f">{case}\n{random_seq(600)}{''.join(bases)}{random_seq(600)}\n")
    return f">{name}-F\n{left}\n>{name}-R\n{right}\n"


def amplified(release: Path, folder: Path, out_name: str, extra: list[str]) -> dict[str, str]:
    """Run insilicoPCR once: case -> what it reported for the left primer."""
    out = folder / out_name
    out.mkdir(parents=True, exist_ok=True)
    java = release / "runtime" / "linux" / "jdk" / "bin" / "java"
    subprocess.run([str(java), "-jar", str(release / "insilicoPCR.jar"),
                    "-i", str(folder / "genomes"), "-p", str(folder / "primers.fasta"),
                    "-o", str(out), "-t", "4", *extra],
                   check=True, capture_output=True)
    found = {}
    with (out / "consolidated_report" / "report.tsv").open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            found[row["Sample"].strip()] = (f"m={row['ForwardMismatches'].strip()}"
                                            f" 3'window={row['Forward3PrimeMismatches'].strip()}")
    return found


def table(title: str, cases: list[str], columns: list[str], results: dict[str, dict[str, str]]) -> None:
    width = max(len(name) for name in cases + ["case"])
    print(f"\n== {title}")
    print(f"{'case':<{width}} | " + " | ".join(f"{c:<6}" for c in columns) + " | reported")
    for name in cases:
        marks = " | ".join(f"{'yes' if name in results[c] else 'no':<6}" for c in columns)
        reported = next((results[c][name] for c in columns if name in results[c]), "")
        print(f"{name:<{width}} | {marks} | {reported}")


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    release = Path(argv[0]).resolve()
    work = Path(argv[1] if len(argv) > 1 else "mismatch_zones_work").resolve()

    # -- 1. one 3'-side mismatch per genome, across -m, at the default --3p-penalty 0 ----------------
    folder = work / "three_prime"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    cases.update({f"run{n}": list(range(1, n + 1)) for n in RUNS})
    left = random_seq(22)
    (folder / "genomes").mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", left, cases))
    results = {f"m={m}": amplified(release, folder, f"m{m}", ["-m", str(m)]) for m in TOLERANCES}
    table("one 3'-side mismatch (offset from the 3' end), by -m, --3p-penalty 0 (default)",
          list(cases), list(results), results)

    # -- 2. the weighted 3' window: the same mismatches, across --3p-penalty, with -m headroom -------
    results = {f"p={p}": amplified(release, folder, f"p{p}", ["-m", "4", "--3p-penalty", str(p)])
               for p in PENALTIES}
    table("the same cases at -m 4, by --3p-penalty", list(cases), list(results), results)

    # -- 3. the G/C rule: a terminal mismatch under a primer ending in A versus one ending in G ------
    folder = work / "terminal_base"
    fasta = assay(folder, "endsat", random_seq(21) + "A", {"at_term": [1]})
    fasta += assay(folder, "endsgc", random_seq(21) + "G", {"gc_term": [1]})
    (folder / "primers.fasta").write_text(fasta)
    results = {f"p={p}": amplified(release, folder, f"p{p}", ["-m", "4", "--3p-penalty", str(p)])
               for p in PENALTIES}
    table("a terminal mismatch, primer ending A vs G, at -m 4, by --3p-penalty",
          ["at_term", "gc_term"], list(results), results)

    # -- 4. the harmless 5' end: one 5'-side mismatch per genome, at -m 0 ----------------------------
    folder = work / "five_prime"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    (folder / "primers.fasta").write_text(assay(folder, "test", random_seq(22), cases, at_five=True))
    results = {"m=0": amplified(release, folder, "m0", ["-m", "0"]),
               "m=1": amplified(release, folder, "m1", ["-m", "1"]),
               "trunc=2": amplified(release, folder, "t2", ["-m", "0", "--5p-truncation", "2"])}
    table("one 5'-side mismatch (offset from the 5' end), -m 0 and 1, and --5p-truncation 2",
          list(cases), list(results), results)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
