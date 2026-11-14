#!/usr/bin/env python3
"""The Xylella panel designed under two insilicoPCR releases on the same find results: 0.11.0 and 0.13.0.

Reads the output folder of run_panel.sh: `<setting>_old/` and `<setting>_new/` for the settings `default`,
`probe2` (--probe-mismatches 2) and `m0` (--mismatches 0), each holding the four targets; each release is
read from the design logs. Primer3 does not always propose the same number of assays on a region, which
renumbers that region's assays, so assays are joined by oligo set (forward, reverse, probe), never by name.

1. Per setting, the two releases: which in silico columns moved, on how many assays, and which way; how
   many selective assays the 3'-end check flags; whether the order of the table moved.
2. Per release, what `--probe-mismatches 2` and `--mismatches 0` do to the verdicts of the default: the
   two comparisons the wiki argues from.

It opens with what the 3'-end check counts at the default settings, read out of the lenient reports
through primer-finder's own `insilico` module: how many of the held (assay, genome) pairs rest on one
imperfection in a primer's last five bases, one in each primer's, or two or more in one primer's -- and
how many `yes` assays would stay flagged if, as primer-finder now requires (MOST_IN_WINDOW), no primer
may carry two.

    python compare_panel.py <output_folder of run_panel.sh>

Runs the primer-finder this script sits in, and needs the commit this record was added in, which is
primer-finder 1.12.0: a later primer-finder refuses the reports insilicoPCR 0.11.0 wrote.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from primer_finder import insilico  # noqa: E402
from primer_finder.design import DEFAULT_MISMATCHES as MISMATCHES  # noqa: E402

TARGETS = ("multiplex", "fastidiosa", "pauca", "morus")
SETTINGS = ("default", "probe2", "m0")
OLIGOS = ("forward", "reverse", "probe")
MODES = ("pcr", "qpcr")
COUNTS = ("inclusion_amplified", "exclusion_amplified", "exclusion_held_by_3prime")
Key = tuple[str, ...]


def load(folder: Path) -> dict[Key, dict[str, str]]:
    """Every assay of the four targets that has a verdict, by target and oligo set, in table order."""
    found: dict[Key, dict[str, str]] = {}
    for target in TARGETS:
        with (folder / target / "assays.tsv").open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                if row.get("pcr_selective"):
                    found[(target, *(row[c] for c in OLIGOS))] = row
    return found


def number(row: dict[str, str], column: str) -> int:
    return int(row[column] or 0)


def moves(before: dict[Key, dict[str, str]], after: dict[Key, dict[str, str]], column: str) -> str:
    """How a count column moved over the assays both tables hold."""
    shared = [key for key in before if key in after and before[key][column] != "" and after[key][column] != ""]
    up = sum(1 for key in shared if number(after[key], column) > number(before[key], column))
    down = sum(1 for key in shared if number(after[key], column) < number(before[key], column))
    return f"up on {up}, down on {down}" if up or down else "the same on every assay"


def verdict_moves(before: dict[Key, dict[str, str]], after: dict[Key, dict[str, str]], mode: str) -> str:
    column = f"{mode}_selective"
    moved = Counter((before[key][column], after[key][column]) for key in before
                    if key in after and before[key][column] and before[key][column] != after[key][column])
    return ", ".join(f"{count} {was} -> {now}" for (was, now), count in sorted(moved.items())) or "none"


def flagged(table: dict[Key, dict[str, str]], mode: str) -> tuple[int, int, Counter[str]]:
    """Of the `yes` assays of a mode, how many the 3'-end check flags, and in which target."""
    yes = [key for key, row in table.items() if row[f"{mode}_selective"] == "yes"]
    held = [key for key in yes if number(table[key], f"{mode}_exclusion_held_by_3prime") > 0]
    return len(held), len(yes), Counter(key[0] for key in held)


def spread(table: dict[Key, dict[str, str]], mode: str) -> str:
    """How many exclusion genomes the flagged `yes` assays are held off: 1, 2, 3-9, 10 or more."""
    bins: Counter[str] = Counter()
    for row in table.values():
        held = number(row, f"{mode}_exclusion_held_by_3prime")
        if row[f"{mode}_selective"] == "yes" and held:
            bins["1" if held == 1 else "2" if held == 2 else "3-9" if held < 10 else "10+"] += 1
    return " / ".join(str(bins[name]) for name in ("1", "2", "3-9", "10+"))


def order(before: dict[Key, dict[str, str]], after: dict[Key, dict[str, str]]) -> str:
    """Whether the tables rank their shared assays in the same order, target by target."""
    said = []
    for target in TARGETS:
        was = [key for key in before if key[0] == target and key in after]
        now = [key for key in after if key[0] == target and key in before]
        displaced = sum(1 for a, b in zip(was, now, strict=True) if a != b)
        top = len(set(was[:20]) & set(now[:20]))
        said.append(f"{target} {displaced} of {len(was)} at another rank, {top} of the first 20 shared"
                    + ("" if was[:1] == now[:1] else ", ANOTHER ASSAY FIRST"))
    return "; ".join(said)


def release(folder: Path) -> str:
    """The insilicoPCR release a design folder was made with, from the version its log records."""
    for line in (folder / f"{TARGETS[0]}.log").read_text().splitlines():
        if "insilicoPCR: insilicoPCR " in line:
            return line.rsplit(" ", 1)[-1]
    return folder.name


def held_by(folder: Path, mode: str) -> tuple[Counter[str], int, int]:
    """The (assay, exclusion genome) pairs a design counts as held by the 3' end, over the four targets,
    by the 3'-window imperfections of the lightest of their products; and how many `yes` assays are
    flagged, as counted and with no primer allowed two imperfections in its window."""
    shapes: Counter[str] = Counter()
    flagged_all: set[tuple[str, str]] = set()
    flagged_one: set[tuple[str, str]] = set()
    for target in TARGETS:
        runs = folder / target / "insilico_pcr"
        strict = insilico.read_report(runs / f"{mode}_exclusion")
        lenient = runs / f"{mode}_exclusion_3prime"
        held = insilico.read_three_prime_held(lenient, MISMATCHES)
        windows: dict[tuple[str, str], list[tuple[int, int]]] = {}
        with (lenient / insilico.REPORT).open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                pair = (row["Gene"].strip(), row["Sample"].strip())
                windows.setdefault(pair, []).append(tuple(int(row[c] or 0) for c in insilico.WINDOW_COLUMNS))
        with (folder / target / "assays.tsv").open(newline="") as fh:
            yes = {row["assay"] for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
                   if row.get(f"{mode}_selective") == "yes"}
        for assay, samples in held.items():
            for sample in samples - strict.get(assay, set()):
                lightest = min(windows[assay, sample], key=sum)
                shapes["two or more in one primer" if max(lightest) >= 2
                       else "one in each primer" if sum(lightest) == 2 else "one"] += 1
                if assay in yes:
                    flagged_all.add((target, assay))
                    if all(max(window) <= 1 for window in windows[assay, sample]):
                        flagged_one.add((target, assay))
    return shapes, len(flagged_all), len(flagged_one)


def main(out: Path) -> int:
    tables = {(setting, arm): load(out / f"{setting}_{arm}") for setting in SETTINGS for arm in ("old", "new")}
    old, new = release(out / "default_old"), release(out / "default_new")

    print("== default: the (assay, exclusion genome) pairs held by the 3' end, by the 3'-window "
          "imperfections of their lightest product")
    for arm, tag in (("old", old), ("new", new)):
        for mode in MODES:
            shapes, flagged_all, flagged_one = held_by(out / f"default_{arm}", mode)
            print(f"   {mode}, {tag}: {sum(shapes.values())} pairs, {dict(shapes)}; `yes` assays flagged "
                  f"{flagged_all}, {flagged_one} when no primer may carry two")

    for setting in SETTINGS:
        before, after = tables[setting, "old"], tables[setting, "new"]
        shared = before.keys() & after.keys()
        print(f"\n== {setting}: {len(before)} assays under {old}, {len(after)} under {new}, "
              f"{len(shared)} with the same oligos")
        for mode in MODES:
            print(f"   {mode}_selective moved: {verdict_moves(before, after, mode)}")
            for column in COUNTS:
                print(f"   {mode}_{column}: {moves(before, after, f'{mode}_{column}')}")
            for tag, table in ((old, before), (new, after)):
                held, yes, by_target = flagged(table, mode)
                print(f"   {mode}, {tag}: {held} of {yes} `yes` assays ({held / yes:.0%}) held off an "
                      f"exclusion genome only by a 3' end; by genomes held (1 / 2 / 3-9 / 10+): "
                      f"{spread(table, mode)}; by target: {dict(by_target)}")
            both = [key for key in shared if before[key][f"{mode}_selective"] == "yes"
                    and after[key][f"{mode}_selective"] == "yes"]
            column = f"{mode}_exclusion_held_by_3prime"
            newly = sum(1 for key in both if not number(before[key], column) and number(after[key], column))
            no_longer = sum(1 for key in both if number(before[key], column) and not number(after[key], column))
            print(f"   {mode}: of the {len(both)} assays `yes` under both, {newly} are flagged under "
                  f"{new} only and {no_longer} under {old} only")
        print(f"   table order: {order(before, after)}")

    for arm, tag in (("old", old), ("new", new)):
        base = tables["default", arm]
        print(f"\n== {tag}: against the default")
        for setting, what in (("probe2", "--probe-mismatches 2"), ("m0", "--mismatches 0")):
            other = tables[setting, arm]
            for mode in MODES:
                print(f"   {what}, {mode}_selective moved: {verdict_moves(base, other, mode)}")
            for mode in MODES:
                column = f"{mode}_exclusion_amplified"
                print(f"   {what}, {column}: {moves(base, other, column)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1])))
