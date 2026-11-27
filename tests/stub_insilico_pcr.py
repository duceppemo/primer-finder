#!/usr/bin/env python3
"""A stand-in for insilicoPCR 2.1.0: writes the files of a run that its real counterpart writes.

`design` reads `consolidated_report/report.tsv`, whose rows are called amplicons; the stub writes one row
per assay and amplifying sample, with the 25 columns of 2.1.0 (`Assay`, `Pairing`, `Source`, `SizeCheck`,
the Tm columns). It ends a primer's name at the first blank of its header and reads a `size=` after it as
the product size the assay declares. Beside the report it writes `report_samples.tsv` (every sample,
analysed or not) and `report_no_product.tsv` (each assay without a product in a sample, and why), named
after the report; it prints `Report written to <report>` on standard output, as the real program does.
It refuses `-m`, `--probe-mismatches` and `--5p-truncation` with exit 2 and the real message, and so
any option it does not know and a `--tm-margin` not below `--anneal-temp`. A sample is named after its
file without the extension (`GCF_1.1_a.fasta` is `GCF_1.1_a`). Like the real program it never
overwrites a report: an existing `report.tsv` makes it write `report(1).tsv` and `report(1)_samples.tsv`
instead, and a numbered report is numbered after the highest there is. What the rows say is read from PRIMER_FINDER_STUBS, as for the other stub
programs:

    insilico_misses      {assay: [sample, ...]}  inclusion samples that do not amplify
    insilico_extra       {assay: [sample, ...]}  exclusion samples that do amplify
    insilico_3prime_held {assay: [sample, ...]}  exclusion samples that amplify only in the lenient
                                                 3'-end run (--3p-penalty above 0), with a 3'-window
                                                 imperfection on their rows -- how the real program
                                                 reports a genome held off by a terminal base
    insilico_one_primer  {assay: [sample, ...]}  samples of either group whose only product is one of
                                                 the forward primer alone (Pairing one-primer)
    insilico_tm          {assay: tm}             the primers' Tm on their sites (default 60.0)
    insilico_no_product_tm {assay: tm}           the Tm of the best reverse-primer site where the assay
                                                 gives no product (default 50.0); "3p" for a best site
                                                 that fails the 3' rule, which has no Tm
    insilico_not_analysed [sample, ...]          samples left out for a faulty file, as 1.0+ does
    insilico_fail        true to exit 1, as a failing run would
    insilico_no_report   true to exit 0 having written no consolidated report at all
    insilico_silent      true to print no "Report written to" line

By default every assay amplifies every sample of the folder whose name holds "inclusion", and none of the
others, which is the answer a selective assay would get.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

COLUMNS = ("Sample", "Assay", "GenomeLocation", "AmpliconSize", "Contig", "ContigDescription",
           "ForwardPrimer", "ReversePrimer", "ForwardMismatches", "ReverseMismatches",
           "Forward3PrimeMismatches", "Reverse3PrimeMismatches", "Strand", "Probe", "ProbeLocation",
           "ProbeSize", "ProbeMismatches", "Forward5PrimeUnpaired", "Reverse5PrimeUnpaired", "Pairing",
           "Source", "SizeCheck", "ForwardTm", "ReverseTm", "ProbeTm")
SAMPLE_COLUMNS = ("Sample", "Status", "Input", "Files", "Products", "Reason")
NO_PRODUCT_COLUMNS = ("Sample", "Assay", "Oligo", "Reason", "BestSite", "BestSiteTm", "NeededTm")
REMOVED = ("-m", "--mismatches", "--probe-mismatches", "--5p-truncation")
KNOWN = {"-i", "-p", "-o", "-t", "-l", "--anneal-temp", "--tm-margin", "--na", "--mg", "--dntp",
         "--oligo-conc", "--3p-window", "--3p-penalty"}


def report_path(folder: Path) -> Path:
    """Where the real program puts the report: report.tsv when neither it nor a numbered one is there,
    else the number after the highest numbered one (report(1).tsv beside report.tsv alone)."""
    numbers = [int(path.name[len("report("):-len(").tsv")]) for path in folder.glob("report(*).tsv")
               if path.name[len("report("):-len(").tsv")].isdigit()]
    plain = folder / "report.tsv"
    if not numbers and not plain.exists():
        return plain
    return folder / f"report({max(numbers, default=0) + 1}).tsv"


def write(path: Path, columns: tuple[str, ...], rows: list[list[str]]) -> None:
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def main(argv: list[str]) -> int:
    if "-v" in argv:
        # The version probe `design` sends before any of the work. Not logged to PRIMER_FINDER_CALLS:
        # the tests count the amplification runs, and this is not one.
        print("insilicoPCR 2.1.0")
        return 0
    calls = os.environ.get("PRIMER_FINDER_CALLS")
    if calls:
        with open(calls, "a") as fh:
            fh.write("insilico_pcr " + " ".join(argv) + "\n")
    for option in REMOVED:
        if option in argv:
            print(f"Error: {option} was removed in 2.0: a primer or probe site is accepted when its Tm "
                  "reaches --anneal-temp (default 60 °C) less --tm-margin (default 5 °C), whatever its "
                  f"mismatches. Leave {option} out; lower --anneal-temp or raise --tm-margin to accept "
                  "weaker sites\nRun insilicoPCR --help for the options.", file=sys.stderr)
            return 2
    options = dict(zip(argv[::2], argv[1::2], strict=False))
    unknown = sorted(set(options) - KNOWN)
    if unknown:
        print(f"Error: unknown option {unknown[0]}\nRun insilicoPCR --help for the options.", file=sys.stderr)
        return 2
    if float(options.get("--tm-margin", "5")) >= float(options.get("--anneal-temp", "60")):
        print("Error: --tm-margin must be below --anneal-temp", file=sys.stderr)
        return 2
    scenario = json.loads(os.environ.get("PRIMER_FINDER_STUBS", "{}"))
    if scenario.get("insilico_fail"):
        print("insilicoPCR: deliberate failure", file=sys.stderr)
        return 1

    genomes = Path(options["-i"])
    output = Path(options["-o"])
    primers = Path(options["-p"]).read_text()
    needed = float(options.get("--anneal-temp", "60")) - float(options.get("--tm-margin", "5"))
    headers = [line[1:].split(None, 1) for line in primers.splitlines() if line.startswith(">")]
    names = [header[0] for header in headers]
    assays = sorted({name.rsplit("-", 1)[0] for name in names})
    with_probe = {name.rsplit("-", 1)[0] for name in names if name.endswith("-P")}
    # The product size an assay declares, on the header of any of its primers
    declared = {header[0].rsplit("-", 1)[0]: word.split("=", 1)[1]
                for header in headers for word in header[1:] and header[1].split()
                if word.startswith("size=")}
    samples = sorted(path.stem for path in genomes.glob("*.fasta"))
    not_analysed = set(scenario.get("insilico_not_analysed", []))
    inclusion = "inclusion" in genomes.name
    lenient = "--3p-penalty" in options and options["--3p-penalty"] != "0"

    if scenario.get("insilico_no_report"):
        return 0
    folder = output / "consolidated_report"
    folder.mkdir(parents=True, exist_ok=True)
    rows, missing = [], []
    products = {sample: 0 for sample in samples}
    for assay in assays:
        misses = set(scenario.get("insilico_misses", {}).get(assay, []))
        extra = set(scenario.get("insilico_extra", {}).get(assay, []))
        held = set(scenario.get("insilico_3prime_held", {}).get(assay, [])) if lenient else set()
        one_primer = set(scenario.get("insilico_one_primer", {}).get(assay, []))
        tm = str(scenario.get("insilico_tm", {}).get(assay, "60.0"))
        for sample in samples:
            if sample in not_analysed:
                continue
            amplifies = (sample not in misses) if inclusion else (sample in extra or sample in held)
            alone = sample in one_primer
            if not amplifies:
                # The assay gave no product of its own, whatever its forward primer made alone: the
                # real program lists it, and says so at the end of the reason
                best = str(scenario.get("insilico_no_product_tm", {}).get(assay, "50.0"))
                why = ("no site is accepted; the best fails the 3' rule, with a 3' penalty of 4, more than "
                       "0 (--3p-penalty) (1 mismatch or indel)" if best == "3p" else
                       f"no site reaches {needed:g} °C; the best melts at {best} °C (2 mismatches and indels)")
                missing.append([sample, assay, f"{assay}-R",
                                why + ("; the report has a product of one of its primers alone" if alone else ""),
                                "chr:1-20(-)", "" if best == "3p" else best, f"{needed:.1f}"])
                if not alone:
                    continue
            # A held sample's binding carries the tolerated 3'-end imperfection: one mismatch in the
            # forward primer's window
            forward, window = ("1", "1") if sample in held and sample not in extra else ("0", "0")
            size = declared.get(assay, "100")
            probe = [f"{assay}-P", "20-40", "20", "0"] if assay in with_probe else ["", "", "", ""]
            pairing, reverse = ("one-primer", f"{assay}-F") if alone and not amplifies else ("designed",
                                                                                            f"{assay}-R")
            row = [sample, assay, f"1-{size}", size, "chr", "", f"{assay}-F", reverse, forward, "0",
                   window, "0", "+", *probe, "0", "0", pairing, "assembly",
                   "target" if assay in declared else "", tm, tm, "68.0" if assay in with_probe else ""]
            rows.append(row)
            products[sample] += 1
    report = report_path(folder)
    write(report, COLUMNS, rows)
    write(report.with_name(f"{report.stem}_samples.tsv"), SAMPLE_COLUMNS, [
        [sample, "not-analysed", "assembly", f"{sample}.fasta", "", "the file holds no sequence"]
        if sample in not_analysed else
        [sample, "analysed", "assembly", f"{sample}.fasta", str(products[sample]), ""]
        for sample in samples])
    write(report.with_name(f"{report.stem}_no_product.tsv"), NO_PRODUCT_COLUMNS, missing)
    if not scenario.get("insilico_silent"):
        print(f"Report written to {report}")
        print(f"HTML report written to {report.with_suffix('.html')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
