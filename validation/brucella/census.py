#!/usr/bin/env python3
"""Where each published Brucella oligo actually occurs in the panel, independently of primer-finder.

Every oligo of assays.tsv is blasted against every panel genome, and only a full-length, mismatch-free hit
counts. For each assay this gives the genomes of its own group that hold all three oligos, and the genomes
of any other group that hold them too.

That is the ceiling on what primer-finder could recover: with its defaults it reports a region only when
the region is in every inclusion genome and differs from every exclusion genome, so an assay missing from
one genome of its own group, or present in one genome of the exclusion group, cannot be reported -- and
should not be. The `-p` a run would need in order to report it anyway is computed here too.

    python census.py /path/to/work -t 16     # needs panel.tsv from panel.py, and blast on PATH

Writes `census.tsv` and `census_detail.tsv` beside it. The three oligos being present says nothing about
whether they would amplify: `census_detail.tsv` records whether they are on one contig, which is what the
paper's two over-fragmented B. suis assemblies are about.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from groups import exclusion_species, inclusion_species  # noqa: E402

PARTS = ("forward", "reverse", "probe")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def panel_stamp(genomes: list[dict[str, str]]) -> str:
    """What the cached fasta, database and hit table were built from.

    Every count below is read fresh from `panel.tsv` while the hits come from these three caches, so a
    panel that has changed since they were written would put a new denominator against an old numerator
    and quietly move every percentage. Written beside the fasta and compared on each run.
    """
    return f"{len(genomes)} genomes\n" + "".join(f"{row['accession']}\t{row['species']}\n"
                                                 for row in sorted(genomes,
                                                                   key=lambda row: row["accession"]))


def panel_fasta(genomes: list[dict[str, str]], out: Path) -> None:
    """Every genome in one file, each record named `accession|species|contig_number`."""
    stamp = out.with_suffix(".panel")
    if out.is_file() and stamp.is_file() and stamp.read_text() == panel_stamp(genomes):
        return
    stamp.write_text(panel_stamp(genomes))
    with out.open("w") as fh:
        for row in genomes:
            number = 0
            for line in Path(row["path"]).read_text().splitlines():
                if line.startswith(">"):
                    number += 1
                    fh.write(f">{row['accession']}|{row['species']}|{number}\n")
                elif line.strip():
                    fh.write(line + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv")
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("-t", "--threads", type=int, default=8)
    args = parser.parse_args(argv)

    genomes = read_tsv(args.work / "panel.tsv")
    assays = read_tsv(args.assays)
    by_species: dict[str, int] = defaultdict(int)
    for row in genomes:
        by_species[row["species"]] += 1
    where = {row["accession"]: row for row in genomes}

    fasta = args.work / "panel.fasta"
    before = fasta.stat().st_mtime if fasta.is_file() else 0.0
    panel_fasta(genomes, fasta)
    rebuilt = fasta.stat().st_mtime != before  # The panel changed, so the two caches below are stale
    database = args.work / "panel"
    if rebuilt or not database.with_suffix(".nin").is_file():
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(database)],
                       check=True, capture_output=True)
    query = args.work / "oligos.fasta"
    with query.open("w") as fh:
        for assay in assays:
            for part in PARTS:
                fh.write(f">{assay['target']}:{part}\n{assay[part]}\n")
    hits = args.work / "oligo_hits.tsv"
    if rebuilt or not hits.is_file():
        subprocess.run(["blastn", "-task", "blastn-short", "-query", str(query), "-db", str(database),
                        "-outfmt", "6 qseqid sseqid pident length qlen mismatch gapopen",
                        "-evalue", "1000", "-max_target_seqs", "10000000",
                        "-num_threads", str(args.threads), "-out", str(hits)], check=True)

    # (target, part) -> {accession -> {contig numbers holding it exactly}}
    exact: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    with hits.open() as fh:
        for query_id, subject, _pident, length, qlen, mismatch, gaps in csv.reader(fh, delimiter="\t"):
            if int(length) == int(qlen) and int(mismatch) == 0 and int(gaps) == 0:
                target, part = query_id.split(":")
                accession, _species, contig = subject.split("|")
                exact[(target, part)][accession].add(contig)

    summary = args.work / "census.tsv"
    detail = args.work / "census_detail.tsv"
    with summary.open("w") as out, detail.open("w") as full:
        out.write("target\tassay\tinclusion_genomes\twith_all_three\tone_contig\tmissing\t"
                  "exclusion_with_all_three\twhich_species\theld_out_with_all_three\tp_needed\n")
        full.write("target\taccession\tspecies\tgroup\tlevel\tparts_found\tone_contig\tnote\n")
        for assay in assays:
            target = assay["target"]
            mine = set(inclusion_species(target))
            found: dict[str, set[str]] = {}
            for accession in where:
                contigs = [exact[(target, part)].get(accession, set()) for part in PARTS]
                if all(contigs):
                    found[accession] = set.intersection(*contigs)
            theirs = exclusion_species(target)
            inclusion = [a for a in found if where[a]["species"] in mine]
            exclusion = [a for a in found if where[a]["group"] == "non_classical"
                         or where[a]["species"] in theirs]
            # Whatever is in neither group: unidentified Brucella sp., and for a species target the four
            # classical_extra species as well. Counted as the remainder so that nothing vanishes.
            held = [a for a in found if a not in set(inclusion) | set(exclusion)]
            total = sum(by_species[name] for name in mine)
            missing = sorted(a for a, row in where.items()
                             if row["species"] in mine and a not in found)
            together = [a for a in inclusion if found[a]]
            others = sorted({where[a]["species"] for a in exclusion})
            out.write("\t".join([
                target, assay["assay"], str(total), str(len(inclusion)), str(len(together)),
                ";".join(missing) if len(missing) <= 10 else f"{len(missing)} genomes",
                str(len(exclusion)), ";".join(others) or "-", str(len(held)),
                f"{len(inclusion) / total:.4f}" if total else "-",
            ]) + "\n")
            for accession in sorted(found):
                row = where[accession]
                full.write("\t".join([target, accession, row["species"], row["group"], row["level"],
                                      "3", "yes" if found[accession] else "no", row["note"]]) + "\n")
    print(f"{summary}: {len(assays)} assay(s) against {len(genomes)} genome(s)")
    print(f"{detail}: every genome holding all three oligos of an assay, and whether on one contig")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
