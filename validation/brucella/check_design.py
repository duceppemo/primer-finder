#!/usr/bin/env python3
"""Ask whether `primer-finder design` proposes the published oligos on the region that holds them.

    python check_design.py /path/to/work/runs/suis_corrected_p0.95/results --target suis -t 16

Finding the region is step one; the paper then designed primers on it with Primer3, and so does
`primer-finder design`. This takes the one region of a finished run that holds the published assay, makes a
results folder containing only that region, runs the design command on it, and compares what came back with
what was published:

- an **exact** match is a proposed oligo whose sequence is the published one;
- an **overlapping** match is one that covers part of the published oligo's footprint, which is what
  agreement on where the assay should sit looks like when Primer3 is given different constraints (the
  paper aimed at amplicons of about 100 bp; `design` is asked for 70-150 here, its own default).

Nothing here is a pass or a fail. Two tools pointed at the same 200-base region will not propose the same
20-mer unless their constraints coincide, and the published assay was designed on a different genome set
with a different tool. What the output answers is the weaker, fairer question: having found the region,
does the design step put an assay in the same place?
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.seqio import Record, iter_records, write_fasta  # noqa: E402

HERE = Path(__file__).resolve().parent
PARTS = ("forward", "reverse", "probe")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def footprint(oligo: str, sequence: str) -> tuple[int, int] | None:
    """Where an oligo sits in a region, on either strand, as a half-open span."""
    for query in (oligo.upper(), reverse_complement(oligo.upper())):
        at = sequence.find(query)
        if at != -1:
            return at, at + len(query)
    return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("results", type=Path, help="The output folder of a finished primer-finder run")
    parser.add_argument("--target", required=True)
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("--output", type=Path, help="Where to design (default: beside the run)")
    parser.add_argument("-t", "--threads", type=int, default=8)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)

    assay = {row["target"]: row for row in read_tsv(args.assays)}[args.target]
    regions = [record for record in iter_records(args.results / "final_kmers.fasta")]
    spans = {}
    chosen = None
    for record in regions:
        found = {part: footprint(assay[part], record.seq.upper()) for part in PARTS}
        if all(found.values()):
            chosen, spans = record, found
            break
    if chosen is None:
        print(f"{assay['assay']}: no region of this run holds all three published oligos", file=sys.stderr)
        return 1
    print(f"{assay['assay']}: designing on {chosen.name} ({len(chosen.seq)} bp), which holds the "
          f"published assay at {spans['forward'][0]}-{spans['reverse'][1]}")

    # A results folder holding that one region, so the design command sees nothing else
    work = args.output or args.results.parent / "design_one_region"
    shutil.rmtree(work, ignore_errors=True)
    (work / "run").mkdir(parents=True)
    write_fasta(work / "run" / "final_kmers.fasta", [Record(chosen.name, chosen.desc, chosen.seq)])
    info = json.loads((args.results / "run_info.json").read_text())
    (work / "run" / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")

    done = subprocess.run(["primer-finder", "design", str(work / "run"), "-o", str(work / "assays"),
                           "-t", str(args.threads)], capture_output=True, text=True)
    if done.returncode != 0:
        sys.stderr.write(done.stderr[-2000:])
        return 1
    print(re.sub(r"^.*\[INFO\] ", "  ", "\n".join(done.stdout.splitlines()[-6:]), flags=re.M))

    rows = read_tsv(work / "assays" / "assays.tsv")
    outcome = {"assay": assay["assay"], "target": args.target, "region": chosen.name,
               "region_bp": len(chosen.seq), "proposed": len(rows), "exact": {}, "overlapping": {}}
    print(f"\n  {len(rows)} assay(s) proposed on that region. Against the published oligos:")
    for part in PARTS:
        published = assay[part].upper()
        start, end = spans[part]
        exact, overlapping = [], []
        for row in rows:
            oligo = (row.get(part) or "").upper()
            if not oligo:
                continue
            at = footprint(oligo, chosen.seq.upper())
            if oligo == published or oligo == reverse_complement(published):
                exact.append(row["assay"])
            elif at and at[0] < end and start < at[1]:
                overlapping.append((row["assay"], at[0] - start, len(oligo)))
        outcome["exact"][part] = exact
        outcome["overlapping"][part] = [one[0] for one in overlapping]
        shown = f"{len(exact)} exact" if exact else "no exact match"
        if overlapping:
            offsets = ", ".join(f"{name} ({shift:+d} bp, {length} nt)"
                               for name, shift, length in overlapping[:3])
            shown += f"; {len(overlapping)} overlapping: {offsets}"
        print(f"    {part:>8s} {published:>28s}  {shown}")
    if args.json:
        args.json.write_text(json.dumps(outcome, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
