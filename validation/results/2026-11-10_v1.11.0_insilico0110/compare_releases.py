#!/usr/bin/env python3
"""Two insilicoPCR releases on the contig-edge genomes, run and read the way `primer-finder design` does.

For each release and each primer file of make_contig_edges.py it makes the two kinds of pass a design
makes, through primer-finder's own `insilico` module and at the design defaults: the strict pass
(`-m 1 -l 10000`), read with `read_report`, and the lenient 3'-end pass (`-m 2 --probe-mismatches 1
--3p-penalty 4`), read with `read_three_prime_held`. It then compares what primer-finder reads -- which
samples each assay amplified, and which it amplified only at the 3' end -- and the reports themselves,
column by column.

A report row that differs is expected to differ in `AmpliconSize` alone, on a product that starts at the
first base of its contig, where a primer's 5' end hangs off: there the newer release must give a size
1-5 bases longer than the older one. Where that product is the one the assay was built for, the newer
size must be the size it was built with (plan.tsv).

    python compare_releases.py <folder of make_contig_edges.py> <insilicoPCR-0.10.1> <insilicoPCR-0.11.0>

Exits 1 if what primer-finder reads differs, or if a report differs in any other way.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from primer_finder import insilico  # noqa: E402
from primer_finder.design import DEFAULT_MISMATCHES  # noqa: E402

THREADS = 8


def pairs(amplified: dict[str, set[str]]) -> int:
    return sum(len(samples) for samples in amplified.values())


def main(work: Path, releases: list[Path]) -> int:
    with (work / "plan.tsv").open(newline="") as fh:
        plan = {row["assay"]: row for row in csv.DictReader(fh, delimiter="\t")}
    genomes = work / "genomes"
    tags: list[str] = []
    read, held, reports = {}, {}, {}
    for release in releases:
        launcher = insilico.resolve(release)
        tag = insilico.check_version(launcher).split()[-1]
        tags.append(tag)
        for kind in ("pcr", "qpcr"):
            primers = work / f"assays_{kind}.fasta"
            strict = insilico.run(launcher, genomes, primers, work / tag / f"{kind}_strict", THREADS,
                                  DEFAULT_MISMATCHES)
            lenient = insilico.run(launcher, genomes, primers, work / tag / f"{kind}_3prime", THREADS,
                                   DEFAULT_MISMATCHES + 1, DEFAULT_MISMATCHES,
                                   three_prime_penalty=insilico.LENIENT_PENALTY)
            read[tag, kind] = insilico.read_report(strict)
            held[tag, kind] = insilico.read_three_prime_held(lenient, DEFAULT_MISMATCHES)
            for name, folder in (("strict", strict), ("3prime", lenient)):
                reports[tag, kind, name] = (folder / insilico.REPORT).read_text().splitlines()
        print(f"ran insilicoPCR {tag}")

    old, new = tags
    wrong = 0
    for kind in ("pcr", "qpcr"):
        same_read = read[old, kind] == read[new, kind]
        same_held = held[old, kind] == held[new, kind]
        wrong += (not same_read) + (not same_held)
        print(f"== {kind}: amplified (assay, sample) pairs {pairs(read[old, kind])} and "
              f"{pairs(read[new, kind])}, {'the same' if same_read else 'DIFFERENT'}; held only by the "
              f"3' end {pairs(held[old, kind])} and {pairs(held[new, kind])}, "
              f"{'the same' if same_held else 'DIFFERENT'}")
        for name in ("strict", "3prime"):
            before, after = reports[old, kind, name], reports[new, kind, name]
            if before == after:
                print(f"   {kind}_{name}: byte-identical, {len(before) - 1} rows")
                continue
            if len(before) != len(after) or before[0] != after[0]:
                wrong += 1
                print(f"   {kind}_{name}: ROWS OR HEADER DIFFER ({len(before) - 1} and {len(after) - 1} rows)")
                continue
            header = before[0].split("\t")
            columns: Counter[str] = Counter()
            differing = at_start = built = not_as_built = 0
            for row_before, row_after in zip(before[1:], after[1:], strict=True):
                if row_before == row_after:
                    continue
                differing += 1
                was, now = (dict(zip(header, row.split("\t"), strict=True)) for row in (row_before, row_after))
                changed = [column for column in header if was[column] != now[column]]
                columns.update(changed)
                longer = int(now["AmpliconSize"]) - int(was["AmpliconSize"])
                at_contig_start = now["GenomeLocation"].startswith("1-")
                if changed != ["AmpliconSize"] or not at_contig_start or not 1 <= longer <= 5:
                    continue
                at_start += 1
                site = plan[now["Gene"]]
                # The product the assay was built for, as against an off-target one that happens to
                # start at the same broken site: the older release has it within 5 bases of the size
                if site["site"] == "forward_off_start" and 0 < int(site["size"]) - int(was["AmpliconSize"]) <= 5:
                    built += now["AmpliconSize"] == site["size"]
                    not_as_built += now["AmpliconSize"] != site["size"]
            print(f"   {kind}_{name}: {len(before) - 1} rows, {differing} differ, in {dict(columns)}; "
                  f"{at_start} of them start at the first base of a contig and are 1-5 bases longer under "
                  f"{new}; {built} of those are the product the assay was built for, which {new} sizes "
                  f"as built" + (f" -- AND {not_as_built} IT DOES NOT" if not_as_built else ""))
            wrong += (differing != at_start) + not_as_built
    print("\nWhat primer-finder reads is the same under both releases." if not wrong
          else f"\n{wrong} UNEXPLAINED DIFFERENCE(S)")
    return 1 if wrong else 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]), [Path(sys.argv[2]), Path(sys.argv[3])]))
