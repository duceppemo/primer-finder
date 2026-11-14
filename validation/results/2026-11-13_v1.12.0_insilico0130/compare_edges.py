#!/usr/bin/env python3
"""Two insilicoPCR releases on the contig-edge genomes, run and read the way `primer-finder design` does.

The genomes and assays are those of ../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py: 300 assays
on 600 kb of *Listeria monocytogenes*, a copy of it broken into contigs inside primer sites, a copy with
one base in 45 substituted, and a near neighbour. That record found two releases reporting the same
amplifications. 0.12.0 changes how a mismatch near the end of an oligo is counted, so here they are not
expected to: this says what moved, and checks that it moved only where the release says it should.

For each release and each primer file it makes the two kinds of pass a design makes, through
primer-finder's own `insilico` module and at the design defaults: the strict pass (`-m 1 -l 10000`), read
with `read_report`, and the lenient 3'-end pass (`-m 2 --probe-mismatches 1 --3p-penalty 4`), read with
`read_three_prime_held`. The primer files are rewritten first with the product size each assay was built
with on the forward primer's header, as primer-finder 1.12.0 writes them.

    python compare_edges.py <folder of make_contig_edges.py> <insilicoPCR-0.11.0> <insilicoPCR-0.12.0 or newer>

Runs the primer-finder this script sits in, and needs the commit this record was added in, which is
primer-finder 1.12.0: a later primer-finder requires insilicoPCR 0.12.0 and refuses both 0.11.0 and its
reports.

What 0.12.0 should leave alone, and this exits 1 if it does not:

- the strict pass of the primer pairs alone: `--3p-penalty 0` refuses any imperfection in the last five
  bases whichever way the trimmed ones are counted;
- every product both releases report: the same place, size and primer mismatch counts.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from primer_finder import insilico  # noqa: E402
from primer_finder.design import DEFAULT_MISMATCHES  # noqa: E402

THREADS = 8
PROBE = 24  # The length of every probe make_contig_edges.py writes
# What names one product in a report, and what must not move on a product both releases report
KEY = ("Sample", "Gene", "GenomeLocation", "Strand")
STABLE = ("AmpliconSize", "ForwardMismatches", "ReverseMismatches", "SizeCheck")
WINDOWS = ("Forward3PrimeMismatches", "Reverse3PrimeMismatches")


def sized(work: Path, plan: dict[str, dict[str, str]], kind: str) -> Path:
    """The primer file with `size=` on each forward primer's header."""
    out = work / "sized" / f"assays_{kind}.fasta"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for line in (work / f"assays_{kind}.fasta").read_text().splitlines():
        if line.startswith(">") and line.endswith("-F"):
            line += f" size={plan[line[1:-2]]['size']}"
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


def products(folder: Path) -> dict[tuple[str, ...], list[dict[str, str]]]:
    """The rows of a report by product. A product has one row, or one per place its probe binds inside
    it: a 6 kb product across two copies of a repeat has two."""
    found: dict[tuple[str, ...], list[dict[str, str]]] = {}
    with (folder / insilico.REPORT).open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            found.setdefault(tuple(row[column].strip() for column in KEY), []).append(row)
    return found


def probes(rows: list[dict[str, str]], column: str) -> tuple[str, ...]:
    """What a product's rows say of its probe, one value per place it binds."""
    return tuple(sorted(row[column].strip() for row in rows))


def pairs(amplified: dict[str, set[str]]) -> set[tuple[str, str]]:
    return {(assay, sample) for assay, samples in amplified.items() for sample in samples}


def main(work: Path, releases: list[Path]) -> int:
    with (work / "plan.tsv").open(newline="") as fh:
        plan = {row["assay"]: row for row in csv.DictReader(fh, delimiter="\t")}
    genomes = work / "genomes"
    tags: list[str] = []
    read, held, rows = {}, {}, {}
    for release in releases:
        launcher = insilico.resolve(release)
        tag = insilico.check_version(launcher).split()[-1]
        tags.append(tag)
        for kind in ("pcr", "qpcr"):
            primers = sized(work, plan, kind)
            strict = insilico.run(launcher, genomes, primers, work / tag / f"{kind}_strict", THREADS,
                                  DEFAULT_MISMATCHES)
            lenient = insilico.run(launcher, genomes, primers, work / tag / f"{kind}_3prime", THREADS,
                                   DEFAULT_MISMATCHES + 1, DEFAULT_MISMATCHES,
                                   three_prime_penalty=insilico.LENIENT_PENALTY)
            read[tag, kind] = pairs(insilico.read_report(strict))
            held[tag, kind] = pairs(insilico.read_three_prime_held(lenient, DEFAULT_MISMATCHES)) - read[tag, kind]
            rows[tag, kind, "strict"], rows[tag, kind, "3prime"] = products(strict), products(lenient)
        print(f"ran insilicoPCR {tag}")

    old, new = tags
    wrong = 0
    for kind in ("pcr", "qpcr"):
        print(f"\n== {kind}")
        for what, sets in (("amplified (assay, sample) pairs", read), ("held only by the 3' end", held)):
            before, after = sets[old, kind], sets[new, kind]
            print(f"   {what}: {len(before)} under {old}, {len(after)} under {new}; "
                  f"{len(before & after)} in both, {len(before - after)} only under {old}, "
                  f"{len(after - before)} only under {new}")
        if kind == "pcr" and read[old, kind] != read[new, kind]:
            wrong += 1
            print("   THE STRICT PASS OF THE PRIMER PAIRS MOVED")
        for name in ("strict", "3prime"):
            before, after = rows[old, kind, name], rows[new, kind, name]
            shared = before.keys() & after.keys()
            gone, come = before.keys() - after.keys(), after.keys() - before.keys()
            moved = Counter(column for key in shared for column in STABLE
                            if before[key][0][column] != after[key][0][column])
            wrong += bool(moved)
            print(f"   {kind}_{name} products: {len(before)} under {old}, {len(after)} under {new}; "
                  f"{len(shared)} in both" + (f", DIFFERING IN {dict(moved)}" if moved else ", the same size "
                  "and primer mismatch counts") + f"; {len(gone)} only under {old}, {len(come)} only under {new}")
            if kind == "qpcr":
                # A product in both whose probe is read differently: the bases blast trimmed, now compared
                relocated = sum(1 for key in shared if probes(before[key], "ProbeLocation")
                                != probes(after[key], "ProbeLocation"))
                recounted = Counter(
                    f"{'+'.join(probes(before[key], 'ProbeMismatches'))} -> "
                    f"{'+'.join(probes(after[key], 'ProbeMismatches'))}" for key in shared
                    if probes(before[key], "ProbeMismatches") != probes(after[key], "ProbeMismatches"))
                print(f"      in both: the probe's place differs on {relocated}; its mismatch count "
                      f"moved on {sum(recounted.values())} {dict(sorted(recounted.items())) or ''}")
            if gone:
                # What the older release said of the products the newer one no longer reports
                said = Counter(
                    f"probe aligned over {'+'.join(probes(before[key], 'ProbeSize'))} of {PROBE} bases, "
                    f"{'+'.join(probes(before[key], 'ProbeMismatches'))} mismatch(es) counted"
                    if kind == "qpcr"
                    else "primer mismatches " + "/".join(before[key][0][c].strip() for c in STABLE[1:3])
                    for key in gone)
                for line, count in sorted(said.items()):
                    print(f"      only under {old}: {count} x {line}")
            if come:
                said = Counter(
                    (f"probe {'+'.join(probes(after[key], 'ProbeMismatches'))} mismatch(es), "
                     if kind == "qpcr" else "")
                    + "3'-window imperfections " + "/".join(after[key][0][c].strip() for c in WINDOWS)
                    for key in come)
                for line, count in sorted(said.items()):
                    print(f"      only under {new}: {count} x {line}")
    print("\nNothing moved that the release does not announce." if not wrong
          else f"\n{wrong} UNEXPLAINED DIFFERENCE(S)")
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), [Path(sys.argv[2]), Path(sys.argv[3])]))
