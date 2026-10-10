#!/usr/bin/env python3
"""Download the complete RefSeq Listeria monocytogenes genomes and give each one a clonal complex.

The clonal complex is what Felix et al. 2023 group by, and it is not in the NCBI metadata: it comes from
the 7-locus MLST profile. `mlst` reads the profile table of the Institut Pasteur scheme (`listeria_2`),
which carries the ST, its clonal complex and its lineage, so the assignment is the scheme's own rather than
anything computed here.

Writes `genomes.tsv`: accession, ST, clonal complex, lineage and the path of the fasta. Genomes whose ST is
new to the scheme, or whose ST the scheme leaves outside every clonal complex, are written with an empty
`cc` and are then left out of both groups by `group_genomes.py` -- they cannot be placed.

    conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
    python type_genomes.py /path/to/work

Takes about ten minutes: a 700 MB download and one `mlst` run per genome.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

SCHEME = "listeria_2"
TAXON = "Listeria monocytogenes"
ACCESSION = re.compile(r"(GCF_\d+\.\d+)")


def download(work: Path) -> Path:
    """Every complete RefSeq assembly of the species. Complete rather than draft because a region broken
    across contigs would look absent, which is a property of the assembly and not of the genome."""
    archive, folder = work / "lm_complete.zip", work / "lm_complete"
    if not folder.is_dir():
        if not archive.is_file():
            subprocess.run(["datasets", "download", "genome", "taxon", TAXON,
                            "--assembly-level", "complete", "--assembly-source", "RefSeq",
                            "--include", "genome", "--filename", str(archive), "--no-progressbar"],
                           check=True)
        shutil.unpack_archive(archive, folder)
    return folder


def profile_table() -> dict[str, tuple[str, str]]:
    """ST -> (clonal complex, lineage), from the scheme mlst has installed."""
    mlst = shutil.which("mlst")
    if mlst is None:
        raise SystemExit("mlst is not on PATH; see the module docstring for the environment")
    table = Path(mlst).resolve().parents[1] / "db" / "pubmlst" / SCHEME / f"{SCHEME}.txt"
    if not table.is_file():
        raise SystemExit(f"{table} is missing: the {SCHEME} scheme is not installed")
    with table.open() as fh:
        return {row["ST"]: (row.get("CC", ""), row.get("Lineage", ""))
                for row in csv.DictReader(fh, delimiter="\t")}


def type_genomes(genomes: list[Path], threads: int) -> dict[Path, str]:
    """The sequence type of each genome, as `mlst` calls it."""
    found: dict[Path, str] = {}
    batch = 20
    for start in range(0, len(genomes), batch * threads):
        chunk = genomes[start:start + batch * threads]
        result = subprocess.run(["mlst", "--scheme", SCHEME, "--quiet", "--threads", str(threads),
                                 *[str(path) for path in chunk]],
                                check=True, capture_output=True, text=True)
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) >= 3:
                found[Path(fields[0])] = fields[2]
        print(f"  typed {len(found)} of {len(genomes)}", file=sys.stderr)
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder for the download and the table.")
    parser.add_argument("-t", "--threads", type=int, default=8)
    args = parser.parse_args(argv)

    args.work.mkdir(parents=True, exist_ok=True)
    folder = download(args.work)
    genomes = sorted(folder.rglob("*.fna"))
    if not genomes:
        raise SystemExit(f"no genome fasta under {folder}")
    print(f"{len(genomes)} complete RefSeq genome(s)", file=sys.stderr)

    scheme = profile_table()
    types = type_genomes(genomes, args.threads)
    out = args.work / "genomes.tsv"
    placed = 0
    with out.open("w") as fh:
        fh.write("accession\tst\tcc\tlineage\tpath\n")
        for path in genomes:
            st = types.get(path, "-")
            cc, lineage = scheme.get(st, ("", ""))
            placed += bool(cc)
            match = ACCESSION.search(path.name)
            fh.write(f"{match.group(1) if match else path.stem}\t{st}\t{cc}\t{lineage}\t{path}\n")
    print(f"{out}: {len(genomes)} genome(s), {placed} with a clonal complex", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
