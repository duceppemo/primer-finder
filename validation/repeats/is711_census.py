#!/usr/bin/env python3
"""Where IS711 really is in the Brucella panel, before asking what primer-finder does with it.

    # the panel and its blast database come from the Brucella harness
    bash validation/brucella/run.sh /path/to/work 48 128        # or at least panel.py and census.py
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M94960&rettype=fasta" \
        -o /path/to/work/is711.fasta
    python validation/repeats/is711_census.py /path/to/work -t 32

Two things are counted for every genome of the panel: whole copies of the 842-base element (GenBank
M94960), and copies of the 63-base amplicon of the published assay. They answer different questions. The
element count is the copy number the literature reports, and says whether the panel looks like the
organism. The amplicon count is what the assay would see, and is the number `-d/--duplication` has to
reach -- as a *mean* over the group, because the ceiling is on the pooled count.

For a draft assembly the two part company, and that is the point of this script: an exact repeat breaks
the assembly at every copy, so the element is cut at a contig end and no whole copy is left, while the
piece holding the amplicon may survive. A marker lost that way is lost before primer-finder sees it.

Writes `is711_census.tsv` beside the panel. Prints the summaries the record quotes.
"""

from __future__ import annotations

import argparse
import collections
import csv
import statistics
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder.mapping import reverse_complement  # noqa: E402
from primer_finder.seqio import iter_records  # noqa: E402

HERE = Path(__file__).resolve().parent

CLASSICAL = ("abortus", "melitensis", "suis", "canis", "ceti", "ovis", "neotomae", "microti")
# A hit this close to the end of a contig was cut off by the assembly, not by biology
EDGE = 50
# What counts as a copy: near-identical, and nearly all of the query. The element's copies diverge by a
# base or two, and `-perc_identity 90` with most of the length aligned is how a copy is usually counted.
MIN_IDENTITY = 90
ELEMENT_COVERAGE = 80
AMPLICON_COVERAGE = 90


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


def amplicon_of(assay: dict[str, str], genome: Path) -> tuple[str, int]:
    """The published amplicon, read off a genome rather than taken on trust.

    Both primers are looked for on both strands, and an amplicon is a forward primer with the reverse
    primer's complement to its right and the probe between them. A published oligo that is nowhere in a
    genome of the organism is a typo in the table, which is why this runs before anything else.
    """
    forward, reverse, probe = assay["forward"], assay["reverse"], assay["probe"]
    found: list[str] = []
    for record in iter_records(genome):
        for sequence in (record.seq.upper(), reverse_complement(record.seq.upper())):
            start = sequence.find(forward)
            while start != -1:
                end = sequence.find(reverse_complement(reverse), start)
                if end != -1 and end - start < 400:
                    piece = sequence[start:end + len(reverse)]
                    if probe in piece or reverse_complement(probe) in piece:
                        found.append(piece)
                start = sequence.find(forward, start + 1)
    if not found:
        raise SystemExit(f"None of the published oligos pair up in {genome.name}: check the table "
                         f"against the genome, the way verify_amplicon.py does.")
    sizes = {len(piece) for piece in found}
    if len(sizes) > 1:
        print(f"note: the amplicon is not one length in this genome: {sorted(sizes)}", file=sys.stderr)
    return found[0], len(found[0])


def hits_of(query: Path, db: Path, out: Path, coverage: int, threads: int) -> list[list[str]]:
    """Every near-identical copy of `query` in the panel database."""
    subprocess.run([
        "blastn", "-task", "blastn", "-query", str(query), "-db", str(db),
        "-outfmt", "6 sseqid pident length qstart qend sstart send slen",
        "-perc_identity", str(MIN_IDENTITY), "-qcov_hsp_perc", str(coverage), "-dust", "no",
        "-max_target_seqs", "1000000", "-num_threads", str(threads), "-out", str(out),
    ], check=True)
    return [line.split("\t") for line in out.read_text().splitlines() if line.strip()]


def accession_of(subject: str) -> str:
    """The panel fasta names each record `accession|species|contig_number`."""
    return subject.split("|")[0]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder holding panel.tsv and the panel blast database")
    parser.add_argument("--assays", type=Path, default=HERE / "is711_assay.tsv")
    parser.add_argument("--element", type=Path, help="The element in fasta. Default: work/is711.fasta")
    parser.add_argument("-t", "--threads", type=int, default=16)
    args = parser.parse_args(argv)

    panel = read_tsv(args.work / "panel.tsv")
    database = args.work / "panel"
    if not (args.work / "panel.nin").is_file() and not (args.work / "panel.nal").is_file():
        raise SystemExit(f"No blast database at {database}: run validation/brucella/census.py first, "
                         f"which builds it from panel.tsv.")
    element = args.element or args.work / "is711.fasta"
    if not element.is_file():
        raise SystemExit(f"{element} is missing. The header of this script has the command that fetches "
                         f"M94960 from NCBI.")
    assay = read_tsv(args.assays)[0]

    reference = next(row for row in sorted(panel, key=lambda row: row["accession"])
                     if row["species"] in CLASSICAL and row["level"] == "Complete Genome")
    sequence, size = amplicon_of(assay, Path(reference["path"]))
    published = int(assay["amplicon_bp"])
    print(f"amplicon measured in {reference['accession']} ({reference['species']}): {size} bp, "
          f"the table says {published} bp -- {'as recorded' if size == published else 'DIFFERENT'}")
    if size != published:
        print("The table and the genome disagree; fix the table before reading anything below.",
              file=sys.stderr)
        return 1
    amplicon_file = args.work / "is711_amplicon.fasta"
    amplicon_file.write_text(f">is711_amplicon\n{sequence}\n")

    whole = collections.Counter()
    for hit in hits_of(element, database, args.work / "is711_element_hits.tsv",
                       ELEMENT_COVERAGE, args.threads):
        whole[accession_of(hit[0])] += 1
    amplicons = collections.Counter()
    for hit in hits_of(amplicon_file, database, args.work / "is711_amplicon_hits.tsv",
                       AMPLICON_COVERAGE, args.threads):
        amplicons[accession_of(hit[0])] += 1
    # Every piece of the element, however short, to tell "absent" from "cut at a contig end" apart
    pieces: dict[str, list[list[str]]] = collections.defaultdict(list)
    for hit in hits_of(element, database, args.work / "is711_element_pieces.tsv", 0, args.threads):
        pieces[accession_of(hit[0])].append(hit)

    out = args.work / "is711_census.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["accession", "species", "group", "level", "contigs", "element_copies",
                         "amplicon_copies", "longest_piece", "piece_at_contig_end"])
        for row in sorted(panel, key=lambda row: row["accession"]):
            here = pieces.get(row["accession"], [])
            longest = max((int(hit[2]) for hit in here), default=0)
            at_edge = any(min(int(hit[5]), int(hit[6])) <= EDGE
                          or max(int(hit[5]), int(hit[6])) >= int(hit[7]) - EDGE for hit in here)
            writer.writerow([row["accession"], row["species"], row["group"], row["level"], row["contigs"],
                             whole.get(row["accession"], 0), amplicons.get(row["accession"], 0),
                             longest, "yes" if at_edge else "no"])
    print(f"\nwrote {out}")

    classical = [row for row in panel if row["species"] in CLASSICAL]
    print(f"\nthe eight classical species: {len(classical)} genomes")
    for what, counter in (("element", whole), ("amplicon", amplicons)):
        counts = [counter.get(row["accession"], 0) for row in classical]
        print(f"  {what:<9} pooled {sum(counts):>6}, mean {sum(counts)/len(counts):>6.2f} per genome, "
              f"none in {sum(1 for c in counts if c == 0):>4} of them")
    mean = sum(amplicons.get(row["accession"], 0) for row in classical) / len(classical)
    print(f"  so -d has to be at least {mean:.2f} for a kmer of the amplicon to survive the ceiling, "
          f"which the whole group shares")

    print(f"\n{'assembly level':<18}{'genomes':>8}{'element=0':>11}{'element med':>12}{'amplicon med':>13}"
          f"{'contigs med':>12}{'cut at an end':>14}")
    by_level = collections.defaultdict(list)
    for row in classical:
        by_level[row["level"]].append(row)
    for level, rows in sorted(by_level.items(), key=lambda pair: -len(pair[1])):
        el = [whole.get(row["accession"], 0) for row in rows]
        am = [amplicons.get(row["accession"], 0) for row in rows]
        none = [row for row in rows if whole.get(row["accession"], 0) == 0]
        cut = sum(1 for row in none
                  if any(min(int(h[5]), int(h[6])) <= EDGE
                         or max(int(h[5]), int(h[6])) >= int(h[7]) - EDGE
                         for h in pieces.get(row["accession"], [])))
        print(f"{level:<18}{len(rows):>8}{len(none):>11}{statistics.median(el):>12}"
              f"{statistics.median(am):>13}{statistics.median(int(row['contigs']) for row in rows):>12}"
              f"{cut:>14}")

    print("\nelement copies in the complete classical genomes, by species "
          "(what the literature reports for the organism):")
    per_species = collections.defaultdict(list)
    for row in classical:
        if row["level"] == "Complete Genome":
            per_species[row["species"]].append(whole.get(row["accession"], 0))
    for species, counts in sorted(per_species.items()):
        print(f"  {species:<12} n={len(counts):<4} median {statistics.median(counts):<6} "
              f"range {min(counts)}-{max(counts)}")

    others = [row for row in panel if row["group"] == "non_classical"]
    print(f"\nthe exclusion group of the genus target: {len(others)} genomes, "
          f"{sum(whole.get(row['accession'], 0) for row in others)} element copies, "
          f"{sum(amplicons.get(row['accession'], 0) for row in others)} amplicon copies")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
