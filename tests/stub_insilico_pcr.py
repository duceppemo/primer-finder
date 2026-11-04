#!/usr/bin/env python3
"""A stand-in for insilicoPCR 0.8.0: writes the consolidated report its real counterpart writes.

`design` reads `consolidated_report/report.tsv`, whose rows are called amplicons since 0.8.0 -- the stub
writes one row per assay and amplifying sample, with the 0.8.0 columns. Like the real program it never
overwrites a report: an existing `report.tsv` makes it write `report(1).tsv` instead, which is what
`insilico.run` clears stale reports to avoid. What the rows say is read from PRIMER_FINDER_STUBS, as for
the other stub programs:

    insilico_misses    {assay: [sample, ...]}  samples of the inclusion group that do not amplify
    insilico_extra     {assay: [sample, ...]}  samples of the exclusion group that do amplify
    insilico_3prime_held {assay: [sample, ...]}  exclusion samples that amplify only in the lenient
                                               3'-end run (--3p-penalty on the command line), with a
                                               3'-window imperfection on their rows -- how the real
                                               program reports a genome held off by a terminal base
    insilico_fail      true to exit non-zero, as a failing run would
    insilico_no_report true to exit 0 having written no consolidated report at all

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
           "Forward3PrimeMismatches", "Reverse3PrimeMismatches", "Strand")
# Appended when any assay in the primer file has a probe; a probe-less assay in such a mixed file keeps
# them empty, as the real program pads them since 0.8.0.
PROBE_COLUMNS = ("Probe", "ProbeLocation", "ProbeSize", "ProbeMismatches")


def report_path(folder: Path) -> Path:
    """Where the real program puts the report: report.tsv, or report(1).tsv beside an existing one."""
    plain = folder / "report.tsv"
    if not plain.exists():
        return plain
    number = 1
    while (folder / f"report({number}).tsv").exists():
        number += 1
    return folder / f"report({number}).tsv"


def main(argv: list[str]) -> int:
    if "-v" in argv:
        # The version probe `design` sends before any of the work. Not logged to PRIMER_FINDER_CALLS:
        # the tests count the amplification runs, and this is not one.
        print("insilicoPCR 0.8.0")
        return 0
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
    primers = Path(options["-p"]).read_text()
    headers = [line[1:].rstrip() for line in primers.splitlines() if line.startswith(">")]
    assays = sorted({name.rsplit("-", 1)[0] for name in headers})
    with_probe = {name.rsplit("-", 1)[0] for name in headers if name.endswith("-P")}
    samples = sorted(path.name.split(".")[0] for path in genomes.glob("*.fasta"))
    inclusion = "inclusion" in genomes.name
    lenient = "--3p-penalty" in options and options["--3p-penalty"] != "0"
    columns = COLUMNS + PROBE_COLUMNS if with_probe else COLUMNS

    if scenario.get("insilico_no_report"):
        return 0
    folder = output / "consolidated_report"
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for assay in assays:
        misses = set(scenario.get("insilico_misses", {}).get(assay, []))
        extra = set(scenario.get("insilico_extra", {}).get(assay, []))
        held = set(scenario.get("insilico_3prime_held", {}).get(assay, [])) if lenient else set()
        for sample in samples:
            amplifies = (sample not in misses) if inclusion else (sample in extra or sample in held)
            if not amplifies:
                continue
            # A held sample's binding carries the tolerated 3'-terminal imperfection: one mismatch in
            # the forward primer's window, which is also the one unit of -m the trimmed base consumes
            forward, window = ("1", "1") if sample in held and sample not in extra else ("0", "0")
            row = [sample, assay, "1-100", "100", "chr", "",
                   f"{assay}-F", f"{assay}-R", forward, "0", window, "0", "+"]
            if with_probe:
                row += ([f"{assay}-P", "20-40", "20", "0"] if assay in with_probe else ["", "", "", ""])
            rows.append(row)
    with report_path(folder).open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
