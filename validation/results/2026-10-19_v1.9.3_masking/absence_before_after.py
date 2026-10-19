#!/usr/bin/env python3
"""What the unmasked second look changes on a real panel.

    python absence_before_after.py <xylella_work_folder> [threads] [memory_GB]

The work folder is the one `validation/xylella/run.sh` fills: `<work>/<target>/{inclusion,exclusion}` for
each of the four published subspecies assays. Each target is run twice -- once as the code stands, and
once with the second, unmasked search disabled, which is what every version up to 1.9.2 did -- and the two
answers are compared: how many regions each reported, how many of them rest on a stretch missing from the
whole exclusion group, and how many bases of that there are in all.

Nothing in the package is changed to get the "before" run: `blast.unmasked_hits` is replaced with a
function that returns nothing, which is exactly what the code did before it existed.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

from primer_finder import blast, pipeline

TARGETS = ("multiplex", "fastidiosa", "pauca", "morus")


def run_one(inclusion: Path, exclusion: Path, output: Path, threads: int, memory: int) -> dict:
    """One `find` run. Returns what it reported, or {} when it reported nothing."""
    shutil.rmtree(output, ignore_errors=True)
    try:
        pipeline.run(pipeline.Settings(inclusion=inclusion, exclusion=exclusion, output=output,
                                       threads=threads, memory_gb=memory))
    except Exception as error:  # noqa: BLE001 - a target with nothing to report is an ordinary outcome
        print(f"    (nothing reported: {type(error).__name__}: {str(error).splitlines()[0][:90]})")
    info_file = output / "run_info.json"
    if not info_file.is_file():
        return {}
    info = json.loads(info_file.read_text())
    absent = info.get("absent_stretches") or {}
    regions = {}
    final = output / "final_kmers.fasta"
    if final.is_file():
        for block in final.read_text().split(">")[1:]:
            head, *rest = block.splitlines()
            regions[head.split()[0]] = len("".join(rest))
    return {
        "counts": info.get("counts", {}),
        "regions": regions,
        "absent_regions": sorted(absent),
        "absent_bases": sum(last - first + 1 for runs in absent.values() for first, last in runs),
    }


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__)
        return 2
    work = Path(argv[0]).resolve()
    threads = int(argv[1]) if len(argv) > 1 else 8
    memory = int(argv[2]) if len(argv) > 2 else 16
    real = blast.unmasked_hits
    rows = []
    for target in TARGETS:
        inclusion, exclusion = work / target / "inclusion", work / target / "exclusion"
        if not inclusion.is_dir() or not exclusion.is_dir():
            print(f"{target}: no such groups under {work}; skipped")
            continue
        print(f"=============== {target}")
        for label, patched in (("after", real), ("before", lambda *_args, **_kwargs: [])):
            blast.unmasked_hits = patched
            print(f"  {label}:")
            answer = run_one(inclusion, exclusion, work / target / f"masking_{label}", threads, memory)
            rows.append((target, label, answer))
        blast.unmasked_hits = real

    print(f"\n{'target':12} {'run':8} {'regions':>8} {'bases':>9} {'with an absent stretch':>23} "
          f"{'absent bases':>13}")
    for target, label, answer in rows:
        if not answer:
            print(f"{target:12} {label:8} {'nothing reported':>8}")
            continue
        print(f"{target:12} {label:8} {len(answer['regions']):>8} "
              f"{sum(answer['regions'].values()):>9} {len(answer['absent_regions']):>23} "
              f"{answer['absent_bases']:>13}")
    by_key = {(target, label): answer for target, label, answer in rows}
    print("\nRegions whose absent stretch the second look took away:")
    for target in TARGETS:
        before, after = by_key.get((target, "before")), by_key.get((target, "after"))
        if not before or not after:
            continue
        gone = sorted(set(before["absent_regions"]) - set(after["absent_regions"]))
        lost = sorted(set(before["regions"]) - set(after["regions"]))
        print(f"  {target}: {len(gone)} of {len(before['absent_regions'])} "
              f"({', '.join(gone[:6])}{' ...' if len(gone) > 6 else ''}); "
              f"regions no longer reported at all: {len(lost)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
