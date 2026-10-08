#!/usr/bin/env python3
"""Download the genomes of accessions.tsv into inclusion/ and exclusion/ folders.

Needs the NCBI datasets command line tool (conda install -c conda-forge ncbi-datasets-cli).
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_accessions(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        rows = [line for line in fh if not line.startswith("#")]
    return list(csv.DictReader(rows, delimiter="\t"))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="Folder to write inclusion/ and exclusion/ to.")
    parser.add_argument("--accessions", type=Path, default=HERE / "accessions.tsv")
    args = parser.parse_args(argv)

    rows = read_accessions(args.accessions)
    archive = args.output / "ncbi.zip"
    args.output.mkdir(parents=True, exist_ok=True)
    if not archive.exists():
        subprocess.run([
            "datasets", "download", "genome", "accession", *[row["accession"] for row in rows],
            "--include", "genome", "--filename", str(archive), "--no-progressbar",
        ], check=True)
    extracted = args.output / "ncbi"
    shutil.rmtree(extracted, ignore_errors=True)
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(extracted)

    missing = []
    for row in rows:
        folder = extracted / "ncbi_dataset" / "data" / row["accession"]
        fasta = next(folder.glob("*.fna"), None) if folder.is_dir() else None
        if fasta is None:
            missing.append(row["accession"])
            continue
        # One folder per subspecies: which of them is the inclusion group depends on the target, and
        # group_genomes.py builds the two folders a run needs out of these.
        folder_out = args.output / "genomes" / row["subspecies"]
        folder_out.mkdir(parents=True, exist_ok=True)
        # The file name carries the accession, the subspecies and the strain, so that the presence table
        # and the log of a run can be read without looking anything up.
        name = f"{row['accession']}_{row['subspecies']}_{row['strain']}.fasta"
        shutil.copyfile(fasta, folder_out / name)
    if missing:
        print(f"Not downloaded: {', '.join(missing)}", file=sys.stderr)
        return 1
    for folder in sorted((args.output / "genomes").iterdir()):
        print(f"{folder.name}: {len(list(folder.glob('*.fasta')))} genome(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
