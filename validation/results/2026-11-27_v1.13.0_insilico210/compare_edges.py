#!/usr/bin/env python3
"""The contig-edge genomes under insilicoPCR 0.14.0 and 2.1.0, each run the way its primer-finder runs it.

The genomes and assays are those of ../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py: 300 assays
on 600 kb of *Listeria monocytogenes*, a copy of it broken into contigs inside primer sites, a copy with
one base in 45 substituted, and a near neighbour. Under 1.x a primer running off a contig end had its
bases past the end free up to `--5p-truncation`; under 2.x they have no template and lower the site's Tm.
So these sets are where the two rules should disagree most, and this says on which pairs, and why.

- old: insilicoPCR 0.14.0 as primer-finder 1.13.0 runs it -- strict `-m 1 -l 10000`, lenient
  `-m 2 --probe-mismatches 1 --3p-penalty 4`;
- new: insilicoPCR 2.1.0 as this primer-finder runs it (`insilico.command`, the design defaults), strict
  and with `--3p-penalty 4`.

The primer files are rewritten first with the product size each assay was built with on the forward
primer's header, as a design writes them.

    python compare_edges.py <folder of make_contig_edges.py> <insilicoPCR-0.14.0> <insilicoPCR-2.1.0>
"""
from __future__ import annotations

import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path

import reports

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from primer_finder import insilico  # noqa: E402

THREADS = 8
REACTION = insilico.Reaction(anneal_temp=60.0, tm_margin=5.0, monovalent_mm=50.0, divalent_mm=3.0,
                             dntp_mm=0.8, oligo_nm=250.0)
OLD_STRICT = ["-m", "1", "-l", str(insilico.LONGEST_PRODUCT)]
OLD_LENIENT = ["-m", "2", "-l", str(insilico.LONGEST_PRODUCT), "--probe-mismatches", "1",
               "--3p-penalty", str(insilico.LENIENT_PENALTY)]


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


def launch(release: Path) -> list[str]:
    return [str(release / "runtime" / "linux" / "jdk" / "bin" / "java"), "-jar", str(release / "insilicoPCR.jar")]


def run_old(release: Path, genomes: Path, primers: Path, out: Path, options: list[str]) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run([*launch(release), "-i", str(genomes), "-p", str(primers), "-o", str(out), "-t", str(THREADS),
                    *options], check=True, capture_output=True)
    return out


def run_new(release: Path, genomes: Path, primers: Path, out: Path, penalty: int | None) -> Path:
    report = insilico.run(launch(release), genomes, primers, out, THREADS, REACTION, penalty)
    insilico.check_samples(report, insilico.count_genomes(genomes))
    return out


def held(lenient: Path, strict: set[reports.Pair], window_rule: bool) -> set[reports.Pair]:
    """The pairs only the lenient pass amplifies, and whose every row one 3'-window imperfection per
    primer explains (both arms' reading; the old one also asked the rest of each primer to fit -m 1)."""
    found: dict[reports.Pair, bool] = {}
    for pair, rows in reports.products(lenient).items():
        for row in rows:
            w = [int(row["Forward3PrimeMismatches"] or 0), int(row["Reverse3PrimeMismatches"] or 0)]
            m = [int(row["ForwardMismatches"] or 0), int(row["ReverseMismatches"] or 0)]
            inside = window_rule or all(a - b <= 1 for a, b in zip(m, w, strict=True))
            ok = any(w) and all(x <= 1 for x in w) and inside
            found[pair] = found.get(pair, True) and ok
    return {pair for pair, ok in found.items() if ok} - strict


def main(work: Path, old_release: Path, new_release: Path) -> int:
    with (work / "plan.tsv").open(newline="") as fh:
        plan = {row["assay"]: row for row in csv.DictReader(fh, delimiter="\t")}
    genomes = work / "genomes"
    for kind in ("pcr", "qpcr"):
        primers = sized(work, plan, kind)
        old_strict = run_old(old_release, genomes, primers, work / "old" / f"{kind}_strict", OLD_STRICT)
        old_lenient = run_old(old_release, genomes, primers, work / "old" / f"{kind}_3prime", OLD_LENIENT)
        new_strict = run_new(new_release, genomes, primers, work / "new" / f"{kind}_strict", None)
        new_lenient = run_new(new_release, genomes, primers, work / "new" / f"{kind}_3prime",
                              insilico.LENIENT_PENALTY)
        old_found, new_found = reports.products(old_strict), reports.products(new_strict)
        old_pairs, new_pairs = set(old_found), set(new_found)
        print(f"\n== {kind}: amplified (assay, genome) pairs: {len(old_pairs)} under 0.14.0, {len(new_pairs)} "
              f"under 2.1.0; {len(old_pairs & new_pairs)} in both")
        none = reports.no_product(new_strict)
        lost, lost_ex = Counter(), {}
        for pair in sorted(old_pairs - new_pairs):
            line = f"{reports.what_made(old_found[pair])} -> {reports.why_none(none.get(pair, []))}"
            lost[line] += 1
            lost_ex.setdefault(line, []).append(f"{pair[0]} in {pair[1]}: was {reports.example(old_found[pair])}; "
                                                f"now {reports.best_site(none.get(pair, []))}")
        gained, gained_ex = Counter(), {}
        for pair in sorted(new_pairs - old_pairs):
            line = reports.what_made(new_found[pair])
            five = max(int(row.get("Forward5PrimeUnpaired") or 0) + int(row.get("Reverse5PrimeUnpaired") or 0)
                       for row in new_found[pair])
            line += f", {five} unpaired 5' base(s)" if five else ""
            gained[line] += 1
            gained_ex.setdefault(line, []).append(f"{pair[0]} in {pair[1]}: {reports.example(new_found[pair])}")
        if lost:
            print(f"   only under 0.14.0: {sum(lost.values())}")
            reports.tally(lost, lost_ex)
        if gained:
            print(f"   only under 2.1.0: {sum(gained.values())}")
            reports.tally(gained, gained_ex)
        by_genome = Counter(sample for _, sample in old_pairs ^ new_pairs)
        print(f"   the pairs that moved, by genome: {dict(by_genome.most_common())}")
        old_held, new_held = held(old_lenient, old_pairs, False), held(new_lenient, new_pairs, True)
        print(f"   held only by the 3' end: {len(old_held)} under 0.14.0, {len(new_held)} under 2.1.0; "
              f"{len(old_held & new_held)} in both")
        # The set's oligos are stretches of the genome taken as they come, not oligos designed to a Tm,
        # and many melt below the needed Tm on their own perfect site: 2.x refuses those wherever they
        # are. What the set was built for -- a primer hanging off a contig end -- is read on the others.
        weak = {pair[0] for pair, why in none.items()
                if any("perfect match is" in row.get("Reason", "") for row in why)}
        strong = [name for name in plan if name not in weak]
        print(f"   assays whose oligos all reach {REACTION.needed_tm:g} C on a perfect site: {len(strong)} of "
              f"{len(plan)}; on lmono_cut, by how many primer bases lie past a contig end:")
        cells: dict[tuple[str, str], list[int]] = {}
        weakest: dict[tuple[str, str], list[float]] = {}
        for name in strong:
            where = (plan[name]["site"], plan[name]["bases_off"])
            cell = cells.setdefault(where, [0, 0, 0])
            cell[0] += 1
            cell[1] += (name, "lmono_cut") in old_pairs
            cell[2] += (name, "lmono_cut") in new_pairs
            weakest.setdefault(where, []).extend(
                tm for tm in (reports.weakest(row) for row in new_found.get((name, "lmono_cut"), []))
                if tm is not None)
        for where, (count, before, after) in sorted(cells.items(), key=lambda c: (c[0][0], int(c[0][1]))):
            tms = sorted(weakest[where])
            low = f", weakest site {tms[0]:.1f}-{tms[-1]:.1f} C" if tms else ""
            print(f"      {where[0]:<18} {where[1]} off: {count:3d} assays, amplified on lmono_cut by "
                  f"{before:3d} under 0.14.0, {after:3d} under 2.1.0{low}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])))
