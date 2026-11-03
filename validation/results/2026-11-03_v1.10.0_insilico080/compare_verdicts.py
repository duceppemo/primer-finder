#!/usr/bin/env python3
"""How the verdicts moved from primer-finder 1.9.3 + insilicoPCR 0.6.1 to 1.10.0 + 0.8.0.

Both designs were run on the SAME find results, so the assays are identical and the verdicts join by
name. For each mode (pcr, qpcr) this prints the old-verdict x new-verdict table -- the question of the
release being where the old `undecided (3' end)` assays land -- and, beside it, what --probe-mismatches 2
changes against the 1.10.0 defaults.

    python compare_verdicts.py <work_folder> [target ...]

<work_folder>/<target>/assays_old/assays.tsv, assays_new/assays.tsv and assays_probe2/assays.tsv are read.
Section 2 of SUMMARY.md is this output.
"""

from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

TARGETS = ("multiplex", "fastidiosa", "pauca", "morus")


def verdicts(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # Only the assays in silico PCR judged: the ones set aside by --min-primer-differences (or refused
    # earlier) carry empty verdict cells in both versions.
    return {row["assay"]: row for row in rows if row.get("pcr_selective")}


def simplify(label: str) -> str:
    """Fold `partial (88%)` into `partial` and `undecided (3' end, 88%)` into `undecided (3' end)`."""
    if label.startswith("partial"):
        return "partial"
    if label.startswith("undecided"):
        return "undecided (3' end)"
    return label


def crosstab(title: str, pairs: Counter) -> None:
    olds = sorted({old for old, _ in pairs}, key=lambda v: ("yes", "partial", "undecided (3' end)", "no").index(v))
    news = sorted({new for _, new in pairs}, key=lambda v: ("yes", "partial", "no").index(v))
    width = max([len(o) for o in olds] + [len("1.9.3 + 0.6.1")])
    print(f"\n== {title}")
    print(f"{'1.9.3 + 0.6.1':<{width}} | " + " | ".join(f"{n:>8}" for n in news) + " |    total")
    for old in olds:
        row = [pairs.get((old, new), 0) for new in news]
        print(f"{old:<{width}} | " + " | ".join(f"{n:>8}" for n in row) + f" | {sum(row):>8}")
    print(f"{'total':<{width}} | "
          + " | ".join(f"{sum(pairs.get((o, n), 0) for o in olds):>8}" for n in news)
          + f" | {sum(pairs.values()):>8}")


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    work = Path(argv[0])
    targets = argv[1:] or TARGETS

    old_all: dict[str, dict[str, str]] = {}
    new_all: dict[str, dict[str, str]] = {}
    probe2_all: dict[str, dict[str, str]] = {}
    for target in targets:
        old, new = verdicts(work / target / "assays_old" / "assays.tsv"), \
                   verdicts(work / target / "assays_new" / "assays.tsv")
        probe2 = verdicts(work / target / "assays_probe2" / "assays.tsv")
        if set(old) != set(new) or set(new) != set(probe2):
            print(f"WARNING: {target} judged different assays across runs "
                  f"(old {len(old)}, new {len(new)}, probe2 {len(probe2)})")
        for name, row in old.items():
            old_all[f"{target}/{name}"] = row
        for name, row in new.items():
            new_all[f"{target}/{name}"] = row
        for name, row in probe2.items():
            probe2_all[f"{target}/{name}"] = row

    shared = sorted(set(old_all) & set(new_all))
    print(f"{len(shared)} assays judged by both versions across {', '.join(targets)}")

    for mode in ("pcr", "qpcr"):
        column = f"{mode}_selective"
        judged = [name for name in shared if old_all[name].get(column) and new_all[name].get(column)]
        pairs = Counter((simplify(old_all[name][column]), simplify(new_all[name][column]))
                        for name in judged)
        crosstab(f"{column}: {len(judged)} assays", pairs)

    # Where exactly the formerly-undecided assays landed, with the inclusion side checked too
    print("\n== the formerly `undecided (3' end)` assays, one line each (pcr mode)")
    for name in shared:
        if not simplify(old_all[name].get("pcr_selective", "")).startswith("undecided"):
            continue
        new = new_all[name]
        print(f"{name}: now {new['pcr_selective']} "
              f"(inclusion {new['pcr_inclusion_amplified']}/{new['inclusion_total']}, "
              f"exclusion {new['pcr_exclusion_amplified']})")

    # --probe-mismatches 2 against the 1.10.0 default, qPCR mode only (the probe plays no part in pcr mode)
    moved = [(name, simplify(new_all[name].get("qpcr_selective", "")),
              simplify(probe2_all[name].get("qpcr_selective", "")))
             for name in sorted(set(new_all) & set(probe2_all))
             if new_all[name].get("qpcr_selective")
             and simplify(new_all[name]["qpcr_selective"]) != simplify(probe2_all[name].get("qpcr_selective", ""))]
    print(f"\n== --probe-mismatches 2 vs the default (qpcr mode): {len(moved)} verdict(s) moved")
    for name, default, loose in moved:
        print(f"{name}: {default} -> {loose}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
