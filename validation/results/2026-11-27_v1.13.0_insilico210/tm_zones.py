#!/usr/bin/env python3
"""What insilicoPCR 2.1.0's site rule accepts, measured one base at a time.

The successor of ../2026-11-13_v1.12.0_insilico0130/mismatch_zones.py. insilicoPCR 2.0 replaced the
mismatch budgets (`-m`, `--probe-mismatches`) with one rule: a site is accepted when its Tm on the template
reaches `--anneal-temp` less `--tm-margin` (55 C by default), and a primer's 3' window passes the weighted
3' rule (`--3p-penalty`, default 0). The method is the old script's: an assay is placed in a synthetic
template and the *template* is mutated at known distances from an oligo's end, one genome per case, so every
row of a report names the case it came from. Where a case gives no product, the report of the assays
without a product (`report_no_product.tsv`) says which oligo failed and the Tm of its best site.

1. one mismatch under the forward primer, by its distance from the 3' end (and runs of 2-4 ending
   there), at the strict --3p-penalty 0 and at the lenient 4 a design's 3'-end check passes;
2. the G/C rule: a terminal mismatch under a primer ending in A against one ending in G, by --3p-penalty;
3. one mismatch by its distance from the 5' end, at the defaults;
4. the 3'-end dependence check of a design: the strict and the lenient exclusion passes with the settings
   `design.run_insilico_pcr` gives them, read through `primer_finder.insilico` as a design reads them;
5. one mismatch under a probe, by its distance from either end;
6. two mismatches under a probe.

    python tm_zones.py /path/to/insilicoPCR-linux-x64 [work_folder]

Needs insilicoPCR 2.1.0 or newer and the primer-finder this script sits in.
"""

from __future__ import annotations

import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from primer_finder import insilico  # noqa: E402

OFFSETS = (1, 2, 3, 4, 5, 6, 8, 10, 12)  # one mismatch, this far from an oligo's end (1 = the terminal base)
RUNS = (2, 3, 4)                         # a run of this many mismatches, ending at the 3' end
PENALTIES = (0, 1, 2, 3, 4)
SPACER = 28
# The first assay primer-finder designs on its example (Contig_1_67.5167_assay4), as Primer3 gave it: the
# primers at 59.8 and 61.0 C, the probe at 63.5 C, in the reaction below. Synthetic oligos of 50-55% G+C
# melt 5 C hotter than that and make a mismatch look cheaper than it is on a designed assay.
FORWARD, REVERSE, PROBE_SEQ = "ACTACTGTCCCATGTGAACT", "CACAGGCATTACGTTAACCG", "TCGGAGTACTAGACCTCAGAATGT"
PROBE = len(PROBE_SEQ)
FLIP = {"A": "G", "C": "T", "G": "A", "T": "C"}
# The reaction and temperature a design passes by default
REACTION = insilico.Reaction(anneal_temp=60.0, tm_margin=5.0, monovalent_mm=50.0, divalent_mm=3.0,
                             dntp_mm=0.8, oligo_nm=250.0)

rng = random.Random(3)
reverse_complement = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))  # noqa: E731


def random_seq(n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def assay(folder: Path, name: str, left: str, cases: dict[str, list[int]], at_five: bool = False,
          probe: bool = False) -> str:
    """One assay and one genome per case, the template mutated at the given offsets from the left
    primer's 3' end -- or from its 5' end with `at_five`. With `probe` the assay has one, and an offset
    is counted along the probe instead: positive from its 5' end, negative from its 3' end. Returns the
    primer fasta records."""
    right, middle = REVERSE, PROBE_SEQ
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


def run(release: Path, folder: Path, out_name: str, penalty: int | None = None) -> Path:
    """Run insilicoPCR once over the folder's genomes, as a design does. Returns its report."""
    java = release / "runtime" / "linux" / "jdk" / "bin" / "java"
    command = insilico.command([str(java), "-jar", str(release / "insilicoPCR.jar")], folder / "genomes",
                               folder / "primers.fasta", folder / out_name, 4, REACTION, penalty)
    said = subprocess.run(command, check=True, capture_output=True, text=True).stdout
    report = Path(insilico.REPORT_WRITTEN.findall(said)[-1])
    insilico.check_samples(report, insilico.count_genomes(folder / "genomes"))
    return report


def outcome(release: Path, folder: Path, out_name: str, penalty: int | None = None) -> dict[str, str]:
    """case -> what the run said: the product's Tm and window counts, or why there was none."""
    report = run(release, folder, out_name, penalty)
    found: dict[str, str] = {}
    for row in insilico._rows(report, is_report=True):
        said = (f"yes  F Tm {row['ForwardTm']} m={row['ForwardMismatches']} "
                f"3'w={row['Forward3PrimeMismatches']} 5'u={row['Forward5PrimeUnpaired']}")
        if row.get("Probe"):
            said += f"  P Tm {row['ProbeTm']} m={row['ProbeMismatches']}"
        found[row["Sample"]] = said
    for row in insilico._rows(insilico.beside(report, "no_product")):
        if row["Sample"] not in found:
            found[row["Sample"]] = f"no   {row['Oligo'].rsplit('-', 1)[-1]}: {row['Reason']}"
    return found


def table(title: str, cases: list[str], results: dict[str, dict[str, str]]) -> None:
    width = max(len(name) for name in cases + ["case"])
    print(f"\n== {title}")
    columns = list(results)
    print(f"{'case':<{width}} | " + " | ".join(f"{c:<7}" for c in columns) + " | first column says")
    for name in cases:
        marks = " | ".join(f"{results[c].get(name, '?')[:3].strip():<{max(7, len(c))}}" for c in columns)
        print(f"{name:<{width}} | {marks} | {results[columns[0]].get(name, '')[5:]}")
        for c in columns[1:]:
            if results[c].get(name, "")[:3] != results[columns[0]].get(name, "")[:3]:
                print(f"{'':<{width}} | {'':<7}   {c}: {results[c].get(name, '')[5:]}")


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    release = Path(argv[0]).resolve()
    work = Path(argv[1] if len(argv) > 1 else "tm_zones_work").resolve()
    print(f"sites need {REACTION.needed_tm:g} C: {' '.join(REACTION.options())}")

    # -- 1. one 3'-side mismatch per genome, strict and lenient --------------------------------------
    folder = work / "three_prime"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    cases.update({f"run{n}": list(range(1, n + 1)) for n in RUNS})
    (folder / "primers.fasta").parent.mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", FORWARD, cases))
    results = {"strict": outcome(release, folder, "strict"),
               "lenient": outcome(release, folder, "lenient", insilico.LENIENT_PENALTY)}
    table(f"1. one 3'-side mismatch (offset from the 3' end), --3p-penalty 0 and {insilico.LENIENT_PENALTY}",
          list(cases), results)

    # -- 2. the G/C rule ------------------------------------------------------------------------------
    # The forward primer with its 3'-terminal base made an A, and a G. One folder each: the two primers
    # differ in that base alone, so the one's mismatched template is the other's perfect one.
    results: dict[str, dict[str, str]] = {f"p={p}": {} for p in PENALTIES}
    for name, end, case in (("endsat", "A", "at"), ("endsgc", "G", "gc")):
        folder = work / "terminal_base" / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "primers.fasta").write_text(
            assay(folder, name, FORWARD[:-1] + end, {f"{case}_term": [1], f"{case}_perfect": []}))
        for p in PENALTIES:
            results[f"p={p}"].update(outcome(release, folder, f"p{p}", p))
    table("2. a terminal mismatch, primer ending A vs G, by --3p-penalty",
          ["at_perfect", "at_term", "gc_perfect", "gc_term"], results)

    # -- 3. one 5'-side mismatch per genome, at the defaults -----------------------------------------
    folder = work / "five_prime"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    cases.update({"off1+off2": [1, 2], "off1-3": [1, 2, 3], "off1-5": [1, 2, 3, 4, 5]})
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", FORWARD, cases, at_five=True))
    results = {"default": outcome(release, folder, "default")}
    table("3. mismatches at the 5' end (offset from the 5' end), at the defaults", list(cases), results)

    # -- 4. the 3'-end dependence check of a design --------------------------------------------------
    # The strict pass and the lenient one, as design.run_insilico_pcr runs them, and the reading
    # insilico.read_three_prime_held makes of the lenient report: a genome is "held" when it amplifies
    # only in the lenient pass and one 3'-window imperfection per primer explains every row of it.
    folder = work / "held"
    cases = {"perfect": []}
    cases.update({f"off{o}": [o] for o in OFFSETS})
    cases.update({f"off{a}+off{b}": [a, b] for a, b in
                  ((1, 2), (1, 5), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),   # both in the window
                   (1, 10), (3, 10), (5, 10),                               # one in it, one internal
                   (8, 12))})                                               # both internal
    cases["off3+off4+off5"] = [3, 4, 5]
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", FORWARD, cases))
    strict_report = run(release, folder, "strict")
    lenient_report = run(release, folder, "lenient", insilico.LENIENT_PENALTY)
    counted = insilico.read_report(strict_report).get("test", set())
    held = insilico.read_three_prime_held(lenient_report).get("test", set()) - counted
    results = {"strict": outcome(release, folder, "strict"),
               "lenient": outcome(release, folder, "lenient", insilico.LENIENT_PENALTY)}
    results["held"] = {name: results["lenient"][name] for name in held}
    lenient_only = set(insilico.read_report(lenient_report).get("test", set())) - counted
    windowless = [name for name in lenient_only if "3'w=0" in results["lenient"].get(name, "")]
    table("4. a design's exclusion passes: strict, lenient (--3p-penalty 4), and what is counted as held "
          "by the 3' end", list(cases), results)
    print(f"   lenient-only genomes: {len(lenient_only)}; of them without a 3'-window imperfection: "
          f"{len(windowless)} {sorted(windowless)}")

    # -- 5. one mismatch under the probe -------------------------------------------------------------
    folder = work / "probe"
    cases = {"perfect": []}
    cases.update({f"5p_off{o}": [o] for o in OFFSETS})
    cases.update({f"3p_off{o}": [-o] for o in OFFSETS})
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", FORWARD, cases, probe=True))
    results = {"default": outcome(release, folder, "default")}
    table("5. one mismatch under the probe (offset from its 5' or its 3' end), at the defaults",
          list(cases), results)

    # -- 6. two and three mismatches under the probe -------------------------------------------------
    folder = work / "probe_two"
    mid = PROBE // 2
    cases = {"middle": [mid], "middle+5p_off1": [mid, 1], "middle+5p_off3": [mid, 3],
             "middle+3p_off1": [mid, -1], "middle+3p_off3": [mid, -3], "5p_off1+3p_off1": [1, -1],
             "middle+middle4": [mid, mid + 4], "three_inside": [mid - 4, mid, mid + 4]}
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "primers.fasta").write_text(assay(folder, "test", FORWARD, cases, probe=True))
    results = {"default": outcome(release, folder, "default")}
    table("6. two or three mismatches under the probe, at the defaults", list(cases), results)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
