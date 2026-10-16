#!/usr/bin/env python3
"""A stand-in for insilicoPCR: writes the consolidated report its real counterpart writes.

`design` reads `consolidated_report/report.tsv` and, for an amplicon with a trimmed primer end,
`detailed_report/<sample>/<sample>.tsv` -- the only place the real program says WHICH end it trimmed. Both
are produced here. What they report is read from PRIMER_FINDER_STUBS, as for the other stub programs:

    insilico_misses    {assay: [sample, ...]}  samples of the inclusion group that do not amplify
    insilico_extra     {assay: [sample, ...]}  samples of the exclusion group that do amplify
    insilico_terminal  {assay: [sample, ...]}  of those, the ones where blast had to trim a primer's 3'
                                               end, which is how the real one reports a difference in the
                                               last two bases: a mismatch it does not count. The reverse
                                               primer is trimmed on every second sample, by two bases, so
                                               that both end columns are exercised
    insilico_five_prime  {assay: [sample, ...]}  the same trim at the primer's 5' end instead, which the
                                               consolidated report cannot be told apart from the above --
                                               only the detailed report can
    insilico_counted   {assay: [sample, ...]}  of those, the ones that bound through a mismatch this check
                                               does count, i.e. one the tolerance let through
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

# The real program writes the probe columns only in qPCR mode: a primer file holding a probe puts the
# whole report in that mode, and the PCR report is the first twelve columns alone.
COLUMNS = ("Sample", "Gene", "GenomeLocation", "AmpliconSize", "Contig", "Contig Description",
           "ForwardPrimers", "ReversePrimers", "ForwardMismatches", "ReverseMismatches",
           "ForwardEndMismatch", "ReverseEndMismatch")
PROBE_COLUMNS = ("Probe", "ProbeLocation", "ProbeSize", "ProbeMismatches")
# The detailed report aligns each primer to the genome, so the genome is the query and the primer the
# subject: `slen` is the primer's length and `sstart`/`send` are positions along it.
DETAILED_COLUMNS = ("qseqid", "sseqid", "slen", "length", "qstart", "qend", "sstart", "send")
PRIMER = 20


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
    primers = Path(options["-p"]).read_text()
    assays = sorted({line[1:].rsplit("-", 1)[0] for line in primers.splitlines()
                     if line.startswith(">")})
    samples = sorted(path.name.split(".")[0] for path in genomes.glob("*.fasta"))
    inclusion = "inclusion" in genomes.name
    qpcr = any(line.startswith(">") and line.rstrip().endswith("-P") for line in primers.splitlines())
    columns = COLUMNS + PROBE_COLUMNS if qpcr else COLUMNS

    if scenario.get("insilico_no_report"):
        return 0
    report = output / "consolidated_report" / "report.tsv"
    report.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    detailed: dict[str, list[list[str]]] = {}
    for assay in assays:
        misses = set(scenario.get("insilico_misses", {}).get(assay, []))
        extra = set(scenario.get("insilico_extra", {}).get(assay, []))
        terminal = set(scenario.get("insilico_terminal", {}).get(assay, []))
        five_prime = set(scenario.get("insilico_five_prime", {}).get(assay, []))
        counted = set(scenario.get("insilico_counted", {}).get(assay, []))
        for number, sample in enumerate(samples):
            amplifies = (sample not in misses) if inclusion else (sample in extra)
            if not amplifies:
                continue
            forward_end, reverse_end = "0", "0"
            if sample in terminal or sample in five_prime:
                # Both columns carry a trim in the real reports, and -1 or -2
                forward_end, reverse_end = ("0", "-2") if number % 2 else ("-1", "0")
            mismatch = "1" if sample in counted else "0"
            at_five = sample in five_prime
            trims = (-int(forward_end), -int(reverse_end))
            # A 5' trim moves the outer end of the amplicon in by the bases it took; a 3' trim does not
            low_end, high_end = (1 + trims[0], 100 - trims[1]) if at_five else (1, 100)
            row = [sample, assay, f"{low_end}-{high_end}", str(high_end - low_end + 1), "chr", "",
                   f"{assay}-F", f"{assay}-R", mismatch, "0", forward_end, reverse_end]
            if qpcr:
                row += [f"{assay}-P", "20-40", "20", "0"]
            rows.append(row)
            # One alignment per primer, reaching its end of the amplicon, saying which end of the primer
            # was left out of it
            for primer, trimmed in ((f"{assay}-F", trims[0]), (f"{assay}-R", trims[1])):
                low, high = (1 + trimmed, PRIMER) if at_five else (1, PRIMER - trimmed)
                span = high - low
                reverse_primer = primer.endswith("-R")
                where = high_end - span if reverse_primer else low_end
                # A reverse primer aligns reversed, so its positions along the primer descend
                on_primer = (high, low) if reverse_primer else (low, high)
                detailed.setdefault(sample, []).append(
                    ["chr", primer, str(PRIMER), str(span + 1),
                     str(where), str(where + span), str(on_primer[0]), str(on_primer[1])])
    with report.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")
    for sample, hits in detailed.items():
        path = output / "detailed_report" / sample / f"{sample}.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as fh:
            fh.write("\t".join(DETAILED_COLUMNS) + "\n")
            for hit in hits:
                fh.write("\t".join(hit) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
