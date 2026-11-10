#!/usr/bin/env python3
"""The Xylella panel designed with insilicoPCR 0.10.1 against the same find results designed with 0.8.0.

The baseline is the `assays_held` arm of ../2026-11-03_3prime_dependence/ (primer-finder fdaf1c9, whose
design code 1.11.0 changes only by a docstring rewrap). The new arm is run_panel.sh's output.

For each target it says whether assays.tsv, both primer fastas and the six consolidated reports are
byte-identical. Primer3 does not always propose the same number of assays on a region (see the 3'-end
record), which renumbers that region's assays; where a region's assay count differs between the arms,
every file is compared again with that region's lines left out. Assays are then joined by oligo set
(forward, reverse, probe) and every in silico column compared.

    python compare_versions.py <baseline_work_folder> <new_work_folder>
"""
import csv
import filecmp
import sys
from collections import Counter
from pathlib import Path

TARGETS = ("multiplex", "fastidiosa", "pauca", "morus")
OLIGOS = ("forward", "reverse", "probe")
INSILICO = ("qpcr_inclusion_amplified", "qpcr_inclusion_percent", "qpcr_exclusion_amplified",
            "qpcr_exclusion_held_by_3prime", "qpcr_selective", "pcr_inclusion_amplified",
            "pcr_inclusion_percent", "pcr_exclusion_amplified", "pcr_exclusion_held_by_3prime",
            "pcr_selective")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def lines_without(path: Path, regions: set[str]) -> list[str]:
    """The file's lines (fasta: whole records), sorted, minus any naming one of the regions -- assay
    names embed the region."""
    lines = path.read_text().splitlines()
    if path.suffix == ".fasta":
        records, current = [], []
        for line in lines:
            if line.startswith(">") and current:
                records.append("\n".join(current))
                current = []
            current.append(line)
        if current:
            records.append("\n".join(current))
        lines = records
    return sorted(line for line in lines if not any(region in line for region in regions))


def main(base_work: Path, new_work: Path) -> int:
    totals = Counter()
    for target in TARGETS:
        base, new = base_work / target / "assays_held", new_work / target
        base_rows, new_rows = rows(base / "assays.tsv"), rows(new / "assays.tsv")
        base_counts = Counter(r["region"] for r in base_rows)
        new_counts = Counter(r["region"] for r in new_rows)
        wobble = {region for region in base_counts | new_counts
                  if base_counts[region] != new_counts[region]}
        print(f"== {target}: {len(base_rows)} assays with 0.8.0, {len(new_rows)} with 0.10.1")
        for region in sorted(wobble):
            print(f"   Primer3 proposed {base_counts[region]} then {new_counts[region]} assays on {region}")

        files = [Path("assays.tsv"), Path("assays_pcr.fasta"), Path("assays_qpcr.fasta")]
        files += sorted(p.relative_to(base) for p in base.glob("insilico_pcr/*/consolidated_report/report.tsv"))
        for rel in files:
            exact = (new / rel).is_file() and filecmp.cmp(base / rel, new / rel, shallow=False)
            outside = exact or ((new / rel).is_file()
                                and lines_without(base / rel, wobble) == lines_without(new / rel, wobble))
            label = rel.parts[1] if rel.parts[0] == "insilico_pcr" else str(rel)
            note = "byte-identical" if exact else (
                "identical outside " + ", ".join(sorted(wobble)) if outside else "DIFFERS")
            print(f"   {label:24s} {note}")
            totals["files"] += 1
            totals["files_byte_identical"] += exact
            totals["files_differing_outside_wobble"] += not outside

        base_by, new_by = ({tuple(r[c] for c in OLIGOS): r for r in table} for table in (base_rows, new_rows))
        shared = base_by.keys() & new_by.keys()
        moved = [k for k in shared if any(base_by[k][c] != new_by[k][c] for c in INSILICO)]
        print(f"   joined by oligo set: {len(shared)} shared, in silico columns differing on {len(moved)}")
        for k in moved:
            print(f"     {new_by[k]['assay']}: "
                  f"{ {c: (base_by[k][c], new_by[k][c]) for c in INSILICO if base_by[k][c] != new_by[k][c]} }")
        totals["assays_shared"] += len(shared)
        totals["assays_moved"] += len(moved)
        totals["assays_only_0.8.0"] += len(base_by.keys() - new_by.keys())
        totals["assays_only_0.10.1"] += len(new_by.keys() - base_by.keys())
    print("\n" + "  ".join(f"{k}={v}" for k, v in totals.items()))
    return 1 if totals["assays_moved"] or totals["files_differing_outside_wobble"] else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), Path(sys.argv[2])))
