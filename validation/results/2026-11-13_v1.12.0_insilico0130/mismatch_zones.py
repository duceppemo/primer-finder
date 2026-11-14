#!/usr/bin/env python3
"""What insilicoPCR's hit filters govern, measured again for 0.12.0.

The successor of ../2026-11-03_v1.10.0_insilico080/mismatch_zones.py, which measured 0.8.0. 0.12.0 compares
the bases blast trims off the end of an alignment with the genome instead of treating every one of them as
unpaired, at the 3' end of a primer and at both ends of a probe, so the zones that script found moved. The
method is the same: an assay is placed in a synthetic template and the *template* is mutated at known
distances from an oligo's end, one genome per case, so every row of a report names the case it came from.

Sections 1-4 are the old script's, so the two outputs can be laid side by side. Sections 5-7 are new:

5. the exclusion rerun `primer-finder design` makes for its 3'-end dependence check, at the default
   `--mismatches 1`, read through `primer_finder.insilico` the way a design reads it;
6. one mismatch under a probe, by its distance from either end, across `--probe-mismatches`;
7. two mismatches under a probe, one of them at an end.

    python mismatch_zones.py /path/to/insilicoPCR-linux-x64 [work_folder]

Needs insilicoPCR 0.11.0 or newer (the primer files declare a product size, as a design's do) and the
primer-finder this script sits in. Under 0.11.0 that has to be the commit this record was added in,
which is primer-finder 1.12.0: a later one refuses the reports of 0.11.0. Prints one table per question;
run it under two releases and diff.
"""

from __future__ import annotations

import csv
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from primer_finder import insilico  # noqa: E402

TOLERANCES = (0, 1, 2, 3, 10)
OFFSETS = range(1, 9)   # one mismatch, this far from an oligo's end (1 = the terminal base)
RUNS = (2, 3, 4)        # a run of this many mismatches, ending at the 3' end
PENALTIES = (0, 1, 2, 3, 4)
PROBE_BUDGETS = (0, 1, 2)
BUDGET = 1              # primer-finder's default --mismatches, which section 5 is measured at
PRIMER, PROBE, SPACER = 22, 24, 28
FLIP = {"A": "G", "C": "T", "G": "A", "T": "C"}

rng = random.Random(3)
random_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
reverse_complement = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))  # noqa: E731


def assay(folder: Path, name: str, left: str, cases: dict[str, list[int]], at_five: bool = False,
          probe: bool = False) -> str:
    """One assay and one genome per case, the template mutated at the given offsets from the left
    primer's 3' end -- or from its 5' end with `at_five`. With `probe` the assay has one, and an offset
    is counted along the probe instead: positive from its 5' end, negative from its 3' end. Returns the
    primer fasta records."""
    right = random_seq(PRIMER)
    middle = random_seq(PROBE)
    amplicon = left + random_seq(SPACER) + middle + random_seq(SPACER) + reverse_complement(right)
    start = len(left) + SPACER  # Where the probe begins
    genomes = folder / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    for case, offsets in cases.items():
        bases = list(amplicon)
        for offset in offsets:
            if probe:
                at = start + offset - 1 if offset > 0 else start + PROBE + offset
            else:
                at = (offset - 1) if at_five else (len(left) - offset)
            bases[at] = FLIP[bases[at]]
        (genomes / f"{case}.fasta").write_text(
            f">{case}\n{random_seq(600)}{''.join(bases)}{random_seq(600)}\n")
    records = f">{name}-F size={len(amplicon)}\n{left}\n>{name}-R\n{right}\n"
    return records + (f">{name}-P\n{middle}\n" if probe else "")


def run(release: Path, folder: Path, out_name: str, extra: list[str]) -> Path:
    """Run insilicoPCR once over the folder's genomes. Returns the folder it wrote."""
    out = folder / out_name
    out.mkdir(parents=True, exist_ok=True)
    java = release / "runtime" / "linux" / "jdk" / "bin" / "java"
    subprocess.run([str(java), "-jar", str(release / "insilicoPCR.jar"),
                    "-i", str(folder / "genomes"), "-p", str(folder / "primers.fasta"),
                    "-o", str(out), "-t", "4", *extra],
                   check=True, capture_output=True)
    return out


def amplified(release: Path, folder: Path, out_name: str, extra: list[str]) -> dict[str, str]:
    """Run insilicoPCR once: case -> what it reported for the left primer, and for the probe if any."""
    found = {}
    report = run(release, folder, out_name, extra) / insilico.REPORT
    with report.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            said = f"m={row['ForwardMismatches'].strip()} 3'window={row['Forward3PrimeMismatches'].strip()}"
            if "Forward5PrimeUnaligned" in row:  # 0.12.0 and newer
                said += f" 5'unpaired={row['Forward5PrimeUnaligned'].strip()}"
            if (row.get("Probe") or "").strip():
                said += f" probe m={row['ProbeMismatches'].strip()} at {row['ProbeLocation'].strip()}"
            found[row["Sample"].strip()] = said
    return found


def table(title: str, cases: list[str], columns: list[str], results: dict[str, dict[str, str]]) -> None:
    width = max(len(name) for name in cases + ["case"])
    print(f"\n== {title}")
    print(f"{'case':<{width}} | " + " | ".join(f"{c:<6}" for c in columns) + " | reported")
    for name in cases:
        marks = " | ".join(f"{'yes' if name in results[c] else 'no':<{max(6, len(c))}}" for c in columns)
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
    left = random_seq(PRIMER)
    (folder / "genomes").mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", left, cases))
    results = {f"m={m}": amplified(release, folder, f"m{m}", ["-m", str(m)]) for m in TOLERANCES}
    table("1. one 3'-side mismatch (offset from the 3' end), by -m, --3p-penalty 0 (default)",
          list(cases), list(results), results)

    # -- 2. the weighted 3' window: the same mismatches, across --3p-penalty, with -m headroom -------
    results = {f"p={p}": amplified(release, folder, f"p{p}", ["-m", "4", "--3p-penalty", str(p)])
               for p in PENALTIES}
    table("2. the same cases at -m 4, by --3p-penalty", list(cases), list(results), results)

    # -- 3. the G/C rule: a terminal mismatch under a primer ending in A versus one ending in G ------
    folder = work / "terminal_base"
    fasta = assay(folder, "endsat", random_seq(PRIMER - 1) + "A", {"at_term": [1]})
    fasta += assay(folder, "endsgc", random_seq(PRIMER - 1) + "G", {"gc_term": [1]})
    (folder / "primers.fasta").write_text(fasta)
    results = {f"p={p}": amplified(release, folder, f"p{p}", ["-m", "4", "--3p-penalty", str(p)])
               for p in PENALTIES}
    table("3. a terminal mismatch, primer ending A vs G, at -m 4, by --3p-penalty",
          ["at_term", "gc_term"], list(results), results)

    # -- 4. the harmless 5' end: one 5'-side mismatch per genome, at -m 0 ----------------------------
    folder = work / "five_prime"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    (folder / "primers.fasta").write_text(assay(folder, "test", random_seq(PRIMER), cases, at_five=True))
    results = {"m=0": amplified(release, folder, "m0", ["-m", "0"]),
               "m=1": amplified(release, folder, "m1", ["-m", "1"]),
               "trunc=2": amplified(release, folder, "t2", ["-m", "0", "--5p-truncation", "2"])}
    table("4. one 5'-side mismatch (offset from the 5' end), -m 0 and 1, and --5p-truncation 2",
          list(cases), list(results), results)

    # -- 5. the 3'-end dependence check of a design, at --mismatches 1 -------------------------------
    # The strict pass and the lenient one, with the settings design.run_insilico_pcr gives them, and the
    # reading insilico.read_three_prime_held makes of the lenient report: a genome is "held" when it
    # amplifies only in the lenient pass and a 3'-window imperfection explains every row of it.
    folder = work / "held"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    cases.update({f"off{a}+off{b}": [a, b] for a, b in
                  ((1, 2), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),   # both in the window
                   (1, 10), (3, 10), (5, 10),                               # one in it, one internal
                   (8, 12))})                                               # both internal
    cases["off3+off4+off5"] = [3, 4, 5]
    (folder / "primers.fasta").write_text(assay(folder, "test", random_seq(PRIMER), cases))
    strict_args = ["-m", str(BUDGET), "-l", str(insilico.LONGEST_PRODUCT)]
    lenient_args = ["-m", str(BUDGET + 1), "-l", str(insilico.LONGEST_PRODUCT),
                    "--probe-mismatches", str(BUDGET), "--3p-penalty", str(insilico.LENIENT_PENALTY)]
    strict = amplified(release, folder, "strict", strict_args)
    lenient = amplified(release, folder, "lenient", lenient_args)
    counted = insilico.read_report(folder / "strict").get("test", set())
    held = insilico.read_three_prime_held(folder / "lenient", BUDGET).get("test", set()) - counted
    results = {"strict": strict, "lenient": lenient, "held": {name: lenient[name] for name in held}}
    table(f"5. a design's exclusion passes at --mismatches {BUDGET}: strict ({' '.join(strict_args)}), "
          f"lenient ({' '.join(lenient_args)}), and what is counted as held by the 3' end",
          list(cases), list(results), results)

    # -- 6. one mismatch under the probe, by its distance from either end ----------------------------
    folder = work / "probe"
    cases = {"perfect": []}
    cases.update({f"5p_off{o}": [o] for o in OFFSETS})
    cases.update({f"3p_off{o}": [-o] for o in OFFSETS})
    cases["middle"] = [PROBE // 2]
    (folder / "primers.fasta").write_text(assay(folder, "test", random_seq(PRIMER), cases, probe=True))
    results = {f"probe={p}": amplified(release, folder, f"probe{p}", ["-m", "1", "--probe-mismatches", str(p)])
               for p in PROBE_BUDGETS}
    table("6. one mismatch under the probe (offset from its 5' or its 3' end), at -m 1, "
          "by --probe-mismatches", list(cases), list(results), results)

    # -- 7. two mismatches under the probe, one of them at an end ------------------------------------
    folder = work / "probe_two"
    cases = {"middle+5p_off1": [PROBE // 2, 1], "middle+5p_off2": [PROBE // 2, 2],
             "middle+5p_off3": [PROBE // 2, 3], "middle+3p_off1": [PROBE // 2, -1],
             "middle+3p_off2": [PROBE // 2, -2], "middle+3p_off3": [PROBE // 2, -3],
             "5p_off1+3p_off1": [1, -1], "middle+middle": [PROBE // 2, PROBE // 2 + 4]}
    (folder / "primers.fasta").write_text(assay(folder, "test", random_seq(PRIMER), cases, probe=True))
    results = {f"probe={p}": amplified(release, folder, f"probe{p}", ["-m", "1", "--probe-mismatches", str(p)])
               for p in PROBE_BUDGETS}
    table("7. two mismatches under the probe, at -m 1, by --probe-mismatches",
          list(cases), list(results), results)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
