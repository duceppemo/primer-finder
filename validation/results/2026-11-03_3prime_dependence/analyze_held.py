#!/usr/bin/env python3
"""What *_exclusion_held_by_3prime finds on the Xylella panel.

Joins the assays_held arm (current code, lenient 3'-end pass) against the 1.10.0 record's assays_new
(verdicts must be identical -- the annotation is not a verdict) and assays_old (where did the formerly
`undecided (3' end)` assays land in the held counts).

    python analyze_held.py <work_folder>
"""
import csv
import sys
from collections import Counter
from pathlib import Path

TARGETS = ("multiplex", "fastidiosa", "pauca", "morus")


def load(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        return {r["assay"]: r for r in csv.DictReader(fh, delimiter="\t") if r.get("pcr_selective")}


def main(work: Path) -> int:
    held_all, new_all, old_all = {}, {}, {}
    for t in TARGETS:
        for store, arm in ((held_all, "assays_held"), (new_all, "assays_new"), (old_all, "assays_old")):
            for n, r in load(work / t / arm / "assays.tsv").items():
                store[f"{t}/{n}"] = r

    # 1. the annotation moved no verdict
    moved = [n for n in held_all if n in new_all and any(
        held_all[n][c] != new_all[n][c] for c in ("pcr_selective", "qpcr_selective"))]
    print(f"verdicts moved vs assays_new: {len(moved)}" + (f" {moved[:5]}" if moved else ""))

    # 2. held counts among the yes assays, per mode
    for mode in ("pcr", "qpcr"):
        sel, held_col = f"{mode}_selective", f"{mode}_exclusion_held_by_3prime"
        yes = [n for n in held_all if held_all[n][sel] == "yes"]
        flagged = [n for n in yes if int(held_all[n][held_col] or 0) > 0]
        dist = Counter(int(held_all[n][held_col]) for n in flagged)
        by_target = Counter(n.split("/")[0] for n in flagged)
        print(f"\n== {mode}: {len(flagged)} of {len(yes)} yes assays held off >=1 exclusion genome "
              f"by the far 3' end alone")
        print(f"   genomes held per assay: {dict(sorted(dist.items()))}")
        print(f"   by target: {dict(by_target)}")

    # 3. the formerly-undecided against the flag
    for mode in ("pcr", "qpcr"):
        sel, held_col = f"{mode}_selective", f"{mode}_exclusion_held_by_3prime"
        undecided = [n for n in old_all if old_all[n][sel].startswith("undecided")]
        flagged = [n for n in undecided if int(held_all[n][held_col] or 0) > 0]
        print(f"\n== {mode}: of the {len(undecided)} formerly `undecided (3' end)`, "
              f"{len(flagged)} are flagged held_by_3prime")
        missed = [n for n in undecided if int(held_all[n][held_col] or 0) == 0]
        for n in missed:
            print(f"   not flagged: {n} (held={held_all[n][held_col]}, {sel}={held_all[n][sel]})")
        # and the converse: flagged yes assays the old model never parked
        yes_flagged = [n for n in held_all if held_all[n][sel] == "yes"
                       and int(held_all[n][held_col] or 0) > 0]
        newly = [n for n in yes_flagged if not old_all[n][sel].startswith("undecided")]
        print(f"   flagged yes assays the old model never parked as undecided: {len(newly)} "
              f"(old verdicts: {dict(Counter(old_all[n][sel].split(' (')[0] for n in newly))})")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
