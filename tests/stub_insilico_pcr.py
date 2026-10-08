#!/usr/bin/env python3
"""A stand-in for insilicoPCR: writes the consolidated report its real counterpart writes.

The `design` command only reads `consolidated_report/report.tsv`, so that is all this produces. What it
reports is read from PRIMER_FINDER_STUBS, as for the other stub programs:

    insilico_misses    {assay: [sample, ...]}  samples of the inclusion group that do not amplify
    insilico_extra     {assay: [sample, ...]}  samples of the exclusion group that do amplify
    insilico_fail      true to exit non-zero, as a failing run would

By default every assay amplifies every sample of the folder whose name holds "inclusion", and none of the
others, which is the answer a selective assay would get.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

COLUMNS = ("Sample", "Gene", "GenomeLocation", "AmpliconSize", "Contig", "Contig Description",
           "ForwardPrimers", "ReversePrimers", "ForwardMismatches", "ReverseMismatches",
           "ForwardEndMismatch", "ReverseEndMismatch", "Probe", "ProbeLocation", "ProbeSize",
           "ProbeMismatches")


def main(argv: list[str]) -> int:
    options = dict(zip(argv[::2], argv[1::2], strict=False))
    scenario = json.loads(os.environ.get("PRIMER_FINDER_STUBS", "{}"))
    calls = os.environ.get("PRIMER_FINDER_CALLS")
    if calls:
        with open(calls, "a") as fh:
            fh.write("insilico_pcr " + " ".join(argv) + "\n")
    if scenario.get("insilico_fail"):
        print("insilicoPCR: deliberate failure", file=sys.stderr)
        return 1

    genomes = Path(options["-i"])
    output = Path(options["-o"])
    assays = sorted({line[1:].rsplit("-", 1)[0] for line in Path(options["-p"]).read_text().splitlines()
                     if line.startswith(">")})
    samples = sorted(path.name.split(".")[0] for path in genomes.glob("*.fasta"))
    inclusion = "inclusion" in genomes.name

    report = output / "consolidated_report" / "report.tsv"
    report.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for assay in assays:
        misses = set(scenario.get("insilico_misses", {}).get(assay, []))
        extra = set(scenario.get("insilico_extra", {}).get(assay, []))
        for sample in samples:
            amplifies = (sample not in misses) if inclusion else (sample in extra)
            if amplifies:
                rows.append([sample, assay, "1-100", "100", "chr", "", f"{assay}-F", f"{assay}-R",
                             "0", "0", "0", "0", f"{assay}-P", "20-40", "20", "0"])
    with report.open("w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
