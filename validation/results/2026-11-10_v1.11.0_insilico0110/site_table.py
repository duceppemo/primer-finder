#!/usr/bin/env python3
"""What each release made of the broken primer sites, from the reports compare_releases.py left.

For the PCR strict pass on `lmono_cut`, it counts per kind of site and per number of bases hanging off the
contig: the assays reported at the size they were built with, those reported only at a shorter size, and
those not reported at all.

    python site_table.py <folder of make_contig_edges.py>
"""
import csv
import sys
from collections import Counter
from pathlib import Path

REPORT = Path("pcr_strict") / "consolidated_report" / "report.tsv"
SAMPLE = "lmono_cut"


def main(work: Path) -> None:
    with (work / "plan.tsv").open(newline="") as fh:
        plan = list(csv.DictReader(fh, delimiter="\t"))
    releases = sorted(folder.name for folder in work.iterdir() if (folder / REPORT).is_file())
    counts: dict[str, Counter[tuple[str, int, str]]] = {}
    for release in releases:
        sizes: dict[str, set[int]] = {}
        with (work / release / REPORT).open(newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                if row["Sample"] == SAMPLE:
                    sizes.setdefault(row["Gene"], set()).add(int(row["AmpliconSize"]))
        counts[release] = Counter()
        for assay in plan:
            found = sizes.get(assay["assay"], set())
            outcome = "none" if not found else ("as built" if int(assay["size"]) in found else "shorter")
            counts[release][assay["site"], int(assay["bases_off"]), outcome] += 1
    print("site\tbases_off\tassays\t" + "\t".join(
        f"{release}: as built / shorter / none" for release in releases))
    for site, bases in sorted({(assay["site"], int(assay["bases_off"])) for assay in plan}):
        total = sum(1 for assay in plan if (assay["site"], int(assay["bases_off"])) == (site, bases))
        cells = [" / ".join(str(counts[release][site, bases, outcome])
                            for outcome in ("as built", "shorter", "none")) for release in releases]
        print(f"{site}\t{bases}\t{total}\t" + "\t".join(cells))


if __name__ == "__main__":
    main(Path(sys.argv[1]))
