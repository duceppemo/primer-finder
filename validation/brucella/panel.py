#!/usr/bin/env python3
"""Build the Brucella genome panel: download every RefSeq assembly of the relevant taxa and label it.

    python panel.py /path/to/work -t 16

The paper's in-silico panel is 1,160 classical Brucella assemblies and 194 of close non-Brucella relatives,
but **its accessions are not published** -- supplementary Table S2 gives only the count per species. The
panel here is therefore rebuilt by the same rule rather than reproduced: every RefSeq assembly, at any
assembly level, of the eight classical species the paper's panel holds, and of the four non-classical
groups it used as negatives. RefSeq has grown since October 2020, so the groups are larger than the
paper's; `panel.tsv` records what was actually used.

Two corrections are applied, both from the paper:

- the six assemblies of its Table S5, whose RefSeq species label disagrees with PubMLST, GTDB-Tk and the
  paper's phylogeny. Grouping on the RefSeq label alone would put them in the wrong group, which would
  make a correct assay look wrong;
- the two over-fragmented B. suis assemblies of its Table S4 are kept and flagged, not dropped. They hold
  the B. suis amplicon split across two contigs, so no PCR can amplify them, and the census is what should
  say so.

The classical-clade species the paper's panel does not hold -- B. pinnipedialis, B. amazoniensis,
B. vulpis and B. inopinata -- are labelled `classical_extra`. The `genus` target ignores them, as the
paper's panel did; the `genus_wide` target puts them in the inclusion group. "Brucella sp." is labelled
`held_out` and goes in neither group whatever the target: it is unidentified, and roughly half of those
assemblies are former Ochrobactrum. The census reports it separately.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

# The species lists live in groups.py, which also says what each target's two groups are made of. The four
# CLASSICAL_EXTRA species are labelled `classical_extra` here rather than held out, so that `genus_wide`
# can ask what spanning the whole clade costs.
from groups import CLASSICAL, CLASSICAL_EXTRA  # noqa: E402

# The non-classical groups it used as negatives. "Brucella (Ochrobactrum)" is the former genus
# Ochrobactrum, which NCBI now files under Brucella: every Brucella species name below is one of those.
OCHROBACTRUM = ("anthropi", "intermedia", "pseudogrignonensis", "pituitosa", "pseudintermedia", "cytisi",
                "tritici", "lupini", "rhizosphaerae", "ciceri", "oryzae", "haematophila", "gallinifaecis",
                "thiophenivorans", "pecoris", "grignonensis", "daejeonensis")
OTHER_GENERA = ("Ochrobactrum", "Pseudochrobactrum", "Paenochrobactrum", "Falsochrobactrum")

# Table S5: RefSeq label -> what the paper's phylogeny makes it. Grouping on the label would misplace these.
REVISED = {
    "GCF_000292125.1": "melitensis",  # labelled B. suis S2-30
    "GCF_009823695.1": "suis",        # labelled B. melitensis HN20190002
    "GCF_000250775.1": "suis",        # labelled B. melitensis S66
    "GCF_000250835.1": "abortus",     # labelled B. melitensis 16M13W
    "GCF_000292145.1": "melitensis",  # labelled B. abortus BCB027
    "GCF_000662015.2": "canis",       # labelled B. inopinata BO1
}

# Table S4: the amplicon is there, but split across two contigs, so no PCR can amplify them.
OVER_FRAGMENTED = ("GCF_000292005.1", "GCF_000292105.1")

TAXA = ("Brucella", "Ochrobactrum", "Pseudochrobactrum", "Paenochrobactrum", "Falsochrobactrum")


def run(command: list[str], **kwargs) -> subprocess.CompletedProcess:
    return subprocess.run(command, check=True, text=True, **kwargs)


def metadata(work: Path) -> list[dict[str, str]]:
    """Accession, organism and assembly level for every RefSeq assembly of the taxa above."""
    out = work / "metadata.tsv"
    if out.is_file():
        with out.open() as fh:
            return list(csv.DictReader(fh, delimiter="\t"))
    rows: list[dict[str, str]] = []
    for taxon in TAXA:
        done = run(["datasets", "summary", "genome", "taxon", taxon, "--assembly-source", "RefSeq",
                    "--as-json-lines"], capture_output=True)
        for line in done.stdout.splitlines():
            if not line.strip():
                continue
            record = json.loads(line)
            info = record.get("assembly_info", {})
            rows.append({
                "accession": record.get("accession", ""),
                "organism": record.get("organism", {}).get("organism_name", ""),
                "level": info.get("assembly_level", ""),
                "contigs": str(record.get("assembly_stats", {}).get("number_of_contigs", "")),
            })
    seen: dict[str, dict[str, str]] = {}
    for row in rows:
        seen.setdefault(row["accession"], row)
    rows = [seen[key] for key in sorted(seen)]
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["accession", "organism", "level", "contigs"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def label(row: dict[str, str]) -> tuple[str, str]:
    """The species and the group of one assembly: a classical species name, `non_classical`, or
    `held_out`."""
    genus, species = (row["organism"].split() + ["", ""])[:2]
    species = REVISED.get(row["accession"], species)
    if genus in OTHER_GENERA:
        return f"{genus.lower()}_{species}", "non_classical"
    if species in CLASSICAL:
        return species, species
    if species in CLASSICAL_EXTRA:
        return species, "classical_extra"
    if species in OCHROBACTRUM:
        return f"ochrobactrum_{species}", "non_classical"
    return species or "sp.", "held_out"


def download(accessions: list[str], work: Path) -> Path:
    """Download the assemblies NCBI will give us, and leave each one as its own fasta file."""
    genomes = work / "genomes"
    genomes.mkdir(parents=True, exist_ok=True)
    missing = [one for one in accessions if not (genomes / f"{one}.fasta").is_file()]
    if missing:
        listed = work / "wanted.txt"
        listed.write_text("".join(f"{one}\n" for one in missing))
        archive = work / "download.zip"
        print(f"Downloading {len(missing)} assembly/assemblies...", flush=True)
        run(["datasets", "download", "genome", "accession", "--inputfile", str(listed),
             "--include", "genome", "--filename", str(archive), "--no-progressbar"])
        with zipfile.ZipFile(archive) as bundle:
            for name in bundle.namelist():
                if not name.endswith((".fna", ".fna.gz")):
                    continue
                parts = Path(name).parts
                accession = next((part for part in parts if part.startswith("GC")), "")
                if not accession:
                    continue
                with bundle.open(name) as source, (genomes / f"{accession}.fasta").open("wb") as plain:
                    shutil.copyfileobj(source, plain)
        archive.unlink()
    return genomes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("work", type=Path, help="Folder to build the panel in")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Unused; kept for a uniform interface")
    args = parser.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=True)

    rows = metadata(args.work)
    print(f"{len(rows)} RefSeq assembly/assemblies in the taxa of the paper's panel")
    for row in rows:
        row["species"], row["group"] = label(row)
    genomes = download([row["accession"] for row in rows], args.work)

    kept, lost = [], []
    for row in rows:
        path = genomes / f"{row['accession']}.fasta"
        if path.is_file() and path.stat().st_size:
            row["path"] = str(path)
            row["note"] = "over_fragmented" if row["accession"] in OVER_FRAGMENTED else ""
            row["revised"] = "yes" if row["accession"] in REVISED else ""
            kept.append(row)
        else:
            lost.append(row["accession"])

    out = args.work / "panel.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "accession", "organism", "species", "group", "level", "contigs", "revised", "note", "path"])
        writer.writeheader()
        writer.writerows(kept)

    counts: dict[str, int] = {}
    for row in kept:
        counts[row["species"]] = counts.get(row["species"], 0) + 1
    classical = sum(counts.get(name, 0) for name in CLASSICAL)
    extra = sum(counts.get(name, 0) for name in CLASSICAL_EXTRA)
    negative = sum(1 for row in kept if row["group"] == "non_classical")
    held = sum(1 for row in kept if row["group"] == "held_out")
    print(f"{out}: {len(kept)} genome(s) -- {classical} classical Brucella of the paper's eight species, "
          f"{extra} of the other four classical-clade species, {negative} non-classical, {held} held out")
    for name in CLASSICAL:
        print(f"  B. {name}: {counts.get(name, 0)}")
    for name in CLASSICAL_EXTRA:
        print(f"  B. {name}: {counts.get(name, 0)} (classical_extra)")
    if lost:
        print(f"{len(lost)} accession(s) NCBI would not give us (suppressed since): {', '.join(lost[:5])}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
