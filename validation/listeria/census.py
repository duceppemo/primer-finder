#!/usr/bin/env python3
"""Where each published oligo actually occurs in the genome panel, independently of primer-finder.

Every oligo of assays.tsv is blasted against every typed genome, and only a full-length, mismatch-free hit
counts. For each assay this gives the number of genomes of its own clonal complex that hold all three
oligos, and the genomes of *other* clonal complexes that hold them too.

That is the ceiling on what primer-finder could recover: with its defaults it reports a region only when the
region is in every inclusion genome and differs from every exclusion genome, so an assay missing from one
genome of its own group, or present in one genome of another, cannot be reported -- and should not be.

    python census.py /path/to/work          # needs genomes.tsv from type_genomes.py, and blast on PATH

Writes `census.tsv` beside it. One caveat: blast counts an IUPAC code in a primer as a mismatch, so the
count for SD_CC1, whose reverse primer carries an M, is a lower bound.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARTS = ("forward", "reverse", "probe")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def panel_fasta(genomes: list[dict[str, str]], out: Path) -> int:
    """Every genome in one file, each record named `accession|cc|n` so a hit says which group it is in."""
    if out.is_file():
        return len(genomes)
    with out.open("w") as fh:
        for row in genomes:
            number = 0
            for line in Path(row["path"]).read_text().splitlines():
                if line.startswith(">"):
                    number += 1
                    fh.write(f">{row['accession']}|{row['cc']}|{number}\n")
                else:
                    fh.write(line + "\n")
    return len(genomes)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding genomes.tsv")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("-t", "--threads", type=int, default=8)
    args = parser.parse_args(argv)

    genomes = [row for row in read_tsv(args.work / "genomes.tsv") if row["cc"]]
    assays = read_tsv(args.assays)
    by_cc: dict[str, int] = defaultdict(int)
    for row in genomes:
        by_cc[row["cc"]] += 1

    fasta = args.work / "panel.fasta"
    panel_fasta(genomes, fasta)
    database = args.work / "panel"
    if not database.with_suffix(".nin").is_file():
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(database)],
                       check=True, capture_output=True)
    query = args.work / "oligos.fasta"
    with query.open("w") as fh:
        for assay in assays:
            for part in PARTS:
                fh.write(f">{assay['target']}:{part}\n{assay[part]}\n")
    hits = args.work / "oligo_hits.tsv"
    if not hits.is_file():
        subprocess.run(["blastn", "-task", "blastn-short", "-query", str(query), "-db", str(database),
                        "-outfmt", "6 qseqid sseqid pident length qlen mismatch gapopen",
                        "-evalue", "1000", "-max_target_seqs", "1000000",
                        "-num_threads", str(args.threads), "-out", str(hits)], check=True)

    exact: dict[tuple[str, str], set[str]] = defaultdict(set)
    with hits.open() as fh:
        for query_id, subject, _pident, length, qlen, mismatch, gaps in csv.reader(fh, delimiter="\t"):
            if int(length) == int(qlen) and int(mismatch) == 0 and int(gaps) == 0:
                target, part = query_id.split(":")
                exact[(target, part)].add(subject.rsplit("|", 1)[0])

    out = args.work / "census.tsv"
    with out.open("w") as fh:
        fh.write("target\tcc\tcc_genomes\twith_all_three\tmissing_from\tother_cc_genomes\twhich\n")
        for assay in assays:
            target = assay["target"]
            cc = target.split("-")[0].replace("SD_", "")
            sets = [exact[(target, part)] for part in PARTS]
            have = set.intersection(*sets) if all(sets) else set()
            mine = sorted(g for g in have if g.endswith("|" + cc))
            others = sorted(g for g in have if not g.endswith("|" + cc))
            missing = sorted(row["accession"] for row in genomes
                             if row["cc"] == cc and f"{row['accession']}|{cc}" not in have)
            fh.write("\t".join([target, cc, str(by_cc.get(cc, 0)), str(len(mine)),
                                ";".join(missing) or "-", str(len(others)),
                                ";".join(others) or "-"]) + "\n")
    print(f"{out}: {len(assays)} assay(s) against {len(genomes)} genome(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
