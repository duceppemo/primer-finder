#!/usr/bin/env python3
"""Score any set of assays against the whole panel, on the same yardstick for every one of them.

    python score_assays.py /path/to/work --assays candidates.tsv --out scorecard.tsv

`census.py` answers this question for the published assays only, and keys its blast queries on the target,
so it cannot hold two assays for the same group. This takes an assays.tsv with any number of rows -- the
published ones, the ones `design` proposed, or both in one file -- keys everything on the assay name, and
reports for each:

- `inclusion_exact`: genomes of the assay's own group holding all three oligos full-length with no
  mismatch, and `inclusion_pct` of the group;
- `one_contig`: how many of those hold them on a single contig, which is what an amplicon needs;
- `exclusion_exact`: genomes of the exclusion group holding all three. For a genus assay that is the
  non-classical relatives; for a species assay, the other classical species as well;
- `held_out_exact`: genomes in neither group holding all three, reported but not counted either way.

Every column is about sequence. None of it says what a reaction would do: an assay whose oligos are
present may still fail, and one with a mismatch may still amplify. Those are wet-lab questions.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARTS = ("forward", "reverse", "probe")
CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")
CLASSICAL_EXTRA = ("pinnipedialis", "amazoniensis", "vulpis", "inopinata")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def panel_database(work: Path) -> Path:
    """The blast database of the whole panel, as census.py builds it."""
    fasta, database = work / "panel.fasta", work / "panel"
    if database.with_suffix(".nin").is_file():
        return database
    if not fasta.is_file():
        raise SystemExit(f"{fasta} is missing: run census.py first")
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(database)],
                   check=True, capture_output=True)
    return database


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv and panel.fasta")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("--out", type=Path, help="Where to write the scorecard (default: beside --assays)")
    parser.add_argument("-t", "--threads", type=int, default=8)
    args = parser.parse_args(argv)

    panel = read_tsv(args.work / "panel.tsv")
    assays = read_tsv(args.assays)
    where = {row["accession"]: row for row in panel}
    by_species: dict[str, int] = defaultdict(int)
    for row in panel:
        by_species[row["species"]] += 1

    database = panel_database(args.work)
    query = args.work / "score_oligos.fasta"
    with query.open("w") as fh:
        for assay in assays:
            for part in PARTS:
                if assay.get(part):
                    fh.write(f">{assay['assay']}:{part}\n{assay[part]}\n")
    hits = args.work / "score_hits.tsv"
    subprocess.run(["blastn", "-task", "blastn-short", "-query", str(query), "-db", str(database),
                    "-outfmt", "6 qseqid sseqid length qlen mismatch gapopen",
                    "-evalue", "1000", "-max_target_seqs", "10000000",
                    "-num_threads", str(args.threads), "-out", str(hits)], check=True)

    exact: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    with hits.open() as fh:
        for query_id, subject, length, qlen, mismatch, gaps in csv.reader(fh, delimiter="\t"):
            if int(length) == int(qlen) and int(mismatch) == 0 and int(gaps) == 0:
                name, part = query_id.rsplit(":", 1)
                accession, _species, contig = subject.split("|")
                exact[(name, part)][accession].add(contig)

    out = args.out or args.assays.with_name(args.assays.stem + "_scorecard.tsv")
    with out.open("w") as fh:
        fh.write("assay\ttarget\tprobe\tamplicon_bp\tgroup_genomes\tinclusion_exact\tinclusion_pct\t"
                 "one_contig\texclusion_exact\texclusion_species\theld_out_exact\n")
        for assay in assays:
            name, target = assay["assay"], assay["target"]
            mine = (CLASSICAL if target == "genus" else
                    CLASSICAL + CLASSICAL_EXTRA if target == "genus_wide" else (target,))
            theirs = () if target.startswith("genus") else tuple(n for n in CLASSICAL if n != target)
            parts = [part for part in PARTS if assay.get(part)]
            found: dict[str, set[str]] = {}
            for accession in where:
                contigs = [exact[(name, part)].get(accession, set()) for part in parts]
                if all(contigs):
                    found[accession] = set.intersection(*contigs)
            inclusion = [a for a in found if where[a]["species"] in mine]
            exclusion = [a for a in found if where[a]["group"] == "non_classical"
                         or where[a]["species"] in theirs]
            held = [a for a in found if a not in set(inclusion) | set(exclusion)]
            total = sum(by_species[one] for one in mine)
            together = [a for a in inclusion if found[a]]
            species = sorted({where[a]["species"] for a in exclusion})
            fh.write("\t".join([
                name, target, "yes" if assay.get("probe") else "no", assay.get("amplicon_bp", "-"),
                str(total), str(len(inclusion)),
                f"{100 * len(inclusion) / total:.2f}" if total else "-",
                str(len(together)), str(len(exclusion)), ";".join(species) or "-", str(len(held)),
            ]) + "\n")
    print(f"{out}: {len(assays)} assay(s) against {len(panel)} genome(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
