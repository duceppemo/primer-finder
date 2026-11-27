#!/usr/bin/env python3
"""The example and the Xylella panel, designed by primer-finder 1.13.0 with insilicoPCR 0.14.0 (the
mismatch budget) and by this primer-finder with insilicoPCR 2.1.0 (the Tm rule), on the same find results.

Reads the output folder of run_arms.sh: `old/<target>/` and `new/<target>/`. Primer3 does not always propose
the same number of assays on a region, which renumbers that region's assays, so assays are joined by oligo
set (forward, reverse, probe), never by name.

Per target and mode (qpcr, pcr):

1. how the verdicts moved, and how the in silico counts moved over the assays both arms hold;
2. every (assay, genome) pair that moved, by why:
   - an inclusion genome lost, or an exclusion genome no longer amplified: the rule insilicoPCR 2.1.0
     says the pair fails, from `report_no_product.tsv` (which oligo, and why);
   - an inclusion genome gained, or an exclusion genome now amplified: what the new product rests on --
     how its primers came together (`Pairing`), the mismatches its primers and probe carry, and its
     weakest site's Tm;
3. the 3'-end check: how many selective assays it flags under each arm;
4. how close the new arm's selective assays are to failing: their weakest inclusion site against the Tm
   a site needs.

    python compare_arms.py <out of run_arms.sh>

Needs nothing of primer-finder: the reports are read by reports.py beside it, which takes both arms'.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import reports

TARGETS = ("example", "multiplex", "fastidiosa", "pauca", "morus")
OLIGOS = ("forward", "reverse", "probe")
MODES = ("qpcr", "pcr")
Key = tuple[str, ...]


def load(folder: Path) -> dict[Key, dict[str, str]]:
    """Every assay that has a verdict, by oligo set, in table order."""
    with (folder / "assays.tsv").open(newline="") as fh:
        return {tuple(row[c] for c in OLIGOS): row
                for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE) if row.get("pcr_selective")}


def pairs(found: dict[reports.Pair, list[dict[str, str]]]) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for assay, sample in found:
        out.setdefault(assay, set()).add(sample)
    return out


def main(out: Path) -> int:
    totals: dict[str, Counter[str]] = {}
    for target in TARGETS:
        old_dir, new_dir = out / "old" / target, out / "new" / target
        if not (old_dir / "assays.tsv").is_file() or not (new_dir / "assays.tsv").is_file():
            print(f"\n######## {target}: not designed under both arms, skipped")
            continue
        old, new = load(old_dir), load(new_dir)
        shared = [key for key in new if key in old]
        print(f"\n######## {target}: {len(old)} assays under the old arm, {len(new)} under the new, "
              f"{len(shared)} with the same oligos; sites need "
              f"{reports.needed_tm(new_dir / 'primer_finder_design.log')} C")
        for mode in MODES:
            judged = [key for key in shared if old[key][f"{mode}_selective"] and new[key][f"{mode}_selective"]]
            if not judged:
                continue
            moved = Counter(f"{old[key][f'{mode}_selective']} -> {new[key][f'{mode}_selective']}"
                            for key in judged if old[key][f"{mode}_selective"] != new[key][f"{mode}_selective"])
            totals.setdefault(mode, Counter()).update(moved)
            print(f"\n== {target} {mode}: {len(judged)} assays judged by both; verdicts moved: "
                  f"{dict(moved.most_common()) or 'none'}")
            for column in ("inclusion_amplified", "exclusion_amplified", "exclusion_held_by_3prime"):
                change = [int(new[key][f"{mode}_{column}"] or 0) - int(old[key][f"{mode}_{column}"] or 0)
                          for key in judged]
                up, down = sum(1 for c in change if c > 0), sum(1 for c in change if c < 0)
                print(f"   {mode}_{column}: up on {up}, down on {down}")

            # Every (assay, genome) pair that moved, by why
            runs = {(arm, group): (out / arm / target / "insilico_pcr" / f"{mode}_{group}")
                    for arm in ("old", "new") for group in ("inclusion", "exclusion", "exclusion_3prime")}
            old_incl = pairs(reports.products(runs["old", "inclusion"]))
            new_incl_found = reports.products(runs["new", "inclusion"], designed_only=True)
            new_incl = pairs(new_incl_found)
            old_excl_found = reports.products(runs["old", "exclusion"])
            old_excl = pairs(old_excl_found)
            new_excl_found = reports.products(runs["new", "exclusion"])
            new_excl = pairs(new_excl_found)
            incl_none, excl_none = (reports.no_product(runs["new", "inclusion"]),
                                    reports.no_product(runs["new", "exclusion"]))
            changes: dict[str, Counter[str]] = {k: Counter() for k in
                                                ("inclusion lost", "inclusion gained", "exclusion gained",
                                                 "exclusion lost")}
            examples: dict[str, dict[str, list[str]]] = {k: {} for k in changes}

            def note(kind: str, line: str, shown: str, changes=changes, examples=examples) -> None:
                changes[kind][line] += 1
                examples[kind].setdefault(line, []).append(shown)

            for key in judged:
                was, now = old[key]["assay"], new[key]["assay"]
                for sample in sorted(old_incl.get(was, set()) - new_incl.get(now, set())):
                    why = incl_none.get((now, sample), [])
                    note("inclusion lost", reports.why_none(why),
                         f"{now} in {sample}: {reports.best_site(why)}")
                for sample in sorted(new_incl.get(now, set()) - old_incl.get(was, set())):
                    found = new_incl_found[(now, sample)]
                    note("inclusion gained", reports.what_made(found),
                         f"{now} in {sample}: {reports.example(found)}")
                for sample in sorted(new_excl.get(now, set()) - old_excl.get(was, set())):
                    found = new_excl_found[(now, sample)]
                    note("exclusion gained", reports.what_made(found),
                         f"{now} in {sample}: {reports.example(found)}")
                for sample in sorted(old_excl.get(was, set()) - new_excl.get(now, set())):
                    why = excl_none.get((now, sample), [])
                    before = reports.example(old_excl_found[(was, sample)])
                    note("exclusion lost", f"{reports.what_made(old_excl_found[(was, sample)])}, now "
                         f"{reports.why_none(why)}",
                         f"{now} in {sample}: was {before}; now {reports.best_site(why)}")
            for kind, lines in changes.items():
                if lines:
                    print(f"   (assay, genome) pairs, {kind}: {sum(lines.values())}")
                    reports.tally(lines, examples[kind])

            # The 3'-end check
            for arm, table in (("old", old), ("new", new)):
                yes = [key for key in judged if table[key][f"{mode}_selective"] == "yes"]
                flagged = sum(1 for key in yes if int(table[key][f"{mode}_exclusion_held_by_3prime"] or 0))
                print(f"   3'-end check, {arm} arm: {flagged} of {len(yes)} `yes` assays held off an exclusion "
                      "genome only by a primer's last five bases")

            # How close the selective assays are to failing
            tms = [float(new[key][f"{mode}_inclusion_weakest_tm"]) for key in judged
                   if new[key][f"{mode}_selective"] == "yes" and new[key][f"{mode}_inclusion_weakest_tm"]]
            if tms:
                bands = Counter("under 56" if tm < 56 else "56-58" if tm < 58 else "58-60" if tm < 60
                                else "60+" for tm in tms)
                print(f"   new arm, `yes` assays by weakest inclusion site: "
                      f"{', '.join(f'{band} C: {bands[band]}' for band in ('under 56', '56-58', '58-60', '60+'))}")
            missing = Counter(new[key][f"{mode}_inclusion_no_product"].split(":")[0] for key in judged
                              if new[key][f"{mode}_inclusion_no_product"])
            if missing:
                print(f"   new arm, assays missing an inclusion genome, by the first oligo named: {dict(missing)}")

        firsts = [(next(iter(table)) if table else None) for table in (old, new)]
        print(f"   first assay of the table: {'the same' if firsts[0] == firsts[1] else 'differs'} "
              f"({old[firsts[0]]['assay'] if firsts[0] else '-'} / "
              f"{new[firsts[1]]['assay'] if firsts[1] else '-'})")
    print("\n######## all targets: verdicts moved")
    for mode, moved in totals.items():
        print(f"   {mode}: {dict(moved.most_common()) or 'none'}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
