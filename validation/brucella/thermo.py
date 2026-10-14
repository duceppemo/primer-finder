#!/usr/bin/env python3
"""Melting temperatures and secondary structure of given oligos, under primer-finder's own conditions.

    python thermo.py --assays candidates.tsv --out thermo.tsv

`design` reports these numbers for the assays it proposes, which makes them useless for comparison on
their own: the published assays were designed elsewhere, under other conditions, and have no such row.
This runs Primer3's `check_primers` task on any set of oligos with the **same** `Conditions` the design
command uses by default, so that both sets are judged by one yardstick.

Reported per assay: each oligo's Tm and GC, the hairpin and self-dimer temperatures Primer3 predicts, the
pair's dimer temperatures, and the spread between the two primer Tms. A `0.0` structure temperature means
no structure was predicted.

Each oligo is checked in a **request of its own**, with every structure limit opened right up, for two
reasons: Primer3 returns nothing at all for a request it rejects, so one oligo over a limit would leave
the whole assay blank; and the point here is to measure the numbers and then say which exceed
primer-finder's defaults, not to let the defaults hide them. `over_limit` names the oligos whose hairpin
or dimer temperature is at or above what `design` would accept (`--max-hairpin-tm`, `--max-dimer-tm`,
47 °C each by default).

This is thermodynamic prediction under one assumed reaction, not a measurement. Primer3's model knows
nothing about the buffer, the polymerase or the cycling of a real assay.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from primer_finder import tools  # noqa: E402
from primer_finder.design import DEFAULT_GC_CLAMP, Conditions  # noqa: E402

HERE = Path(__file__).resolve().parent
PROGRAM = "primer3_core"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader([line for line in fh if not line.startswith("#")], delimiter="\t"))


OPEN = {  # Limits opened right up: this task is asked for numbers, not for a verdict
    "PRIMER_MAX_HAIRPIN_TH": 100.0,
    "PRIMER_MAX_SELF_ANY_TH": 100.0,
    "PRIMER_MAX_SELF_END_TH": 100.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 100.0,
    "PRIMER_PAIR_MAX_COMPL_END_TH": 100.0,
    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": 100.0,
    "PRIMER_INTERNAL_MAX_SELF_ANY_TH": 100.0,
    "PRIMER_INTERNAL_MAX_SELF_END_TH": 100.0,
    "PRIMER_MIN_SIZE": 15, "PRIMER_MAX_SIZE": 36,
    "PRIMER_MIN_TM": 20.0, "PRIMER_MAX_TM": 95.0,
    "PRIMER_MIN_GC": 5.0, "PRIMER_MAX_GC": 95.0,
    "PRIMER_INTERNAL_MIN_SIZE": 15, "PRIMER_INTERNAL_MAX_SIZE": 36,
    "PRIMER_INTERNAL_MIN_TM": 20.0, "PRIMER_INTERNAL_MAX_TM": 95.0,
    "PRIMER_INTERNAL_MIN_GC": 5.0, "PRIMER_INTERNAL_MAX_GC": 95.0,
    "PRIMER_MAX_NS_ACCEPTED": 5,
    "PRIMER_GC_CLAMP": 0,  # Not a constraint here: a published primer is what it is
    "PRIMER_MAX_POLY_X": 36,
    "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
    "PRIMER_EXPLAIN_FLAG": 1,
}


def conditions_tags(conditions: Conditions) -> dict[str, object]:
    """The reaction, taken from the design command's defaults so both sets of assays match."""
    return {
        "PRIMER_SALT_MONOVALENT": conditions.monovalent_mm,
        "PRIMER_SALT_DIVALENT": conditions.divalent_mm,
        "PRIMER_DNTP_CONC": conditions.dntp_mm,
        "PRIMER_DNA_CONC": conditions.primer_nm,
        "PRIMER_INTERNAL_SALT_MONOVALENT": conditions.monovalent_mm,
        "PRIMER_INTERNAL_SALT_DIVALENT": conditions.divalent_mm,
        "PRIMER_INTERNAL_DNTP_CONC": conditions.dntp_mm,
        "PRIMER_INTERNAL_DNA_CONC": conditions.probe_nm,
    }


def record(name: str, conditions: Conditions, extra: dict[str, object]) -> str:
    settings: dict[str, object] = {"SEQUENCE_ID": name, "PRIMER_TASK": "check_primers"}
    settings.update(OPEN)
    settings.update(conditions_tags(conditions))
    settings.update(extra)
    return "".join(f"{key}={value}\n" for key, value in settings.items()) + "=\n"


def requests(assay: dict[str, str], conditions: Conditions) -> list[tuple[str, str]]:
    """One request per oligo, plus one for the pair, each tagged with what it is for."""
    name = assay["assay"]
    out = [
        ("left", record(f"{name}:left", conditions, {
            "SEQUENCE_PRIMER": assay["forward"].upper(),
            "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0,
            "PRIMER_PICK_INTERNAL_OLIGO": 0})),
        ("right", record(f"{name}:right", conditions, {
            "SEQUENCE_PRIMER_REVCOMP": assay["reverse"].upper(),
            "PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0})),
        ("pair", record(f"{name}:pair", conditions, {
            "SEQUENCE_PRIMER": assay["forward"].upper(),
            "SEQUENCE_PRIMER_REVCOMP": assay["reverse"].upper(),
            "PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 1,
            "PRIMER_PICK_INTERNAL_OLIGO": 0})),
    ]
    if assay.get("probe"):
        out.append(("internal", record(f"{name}:internal", conditions, {
            "SEQUENCE_INTERNAL_OLIGO": assay["probe"].upper(),
            "PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 0,
            "PRIMER_PICK_INTERNAL_OLIGO": 1})))
    return out


def run_primer3(records: str) -> list[dict[str, str]]:
    """Every boulder-IO record in one call, parsed back into one dictionary each."""
    done = subprocess.run([PROGRAM], input=records, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        raise SystemExit(f"{PROGRAM} failed: {done.stderr[-500:]}")
    out: list[dict[str, str]] = []
    current: dict[str, str] = {}
    for line in done.stdout.splitlines():
        if line == "=":
            out.append(current)
            current = {}
        elif "=" in line:
            key, value = line.split("=", 1)
            current[key] = value
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--assays", type=Path, default=HERE / "assays.tsv")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    tools.require([PROGRAM])

    assays = read_tsv(args.assays)
    conditions = Conditions()
    limit_hairpin, limit_dimer = conditions.max_hairpin_tm, conditions.max_dimer_tm

    batch, kinds = [], []
    for assay in assays:
        for kind, text in requests(assay, conditions):
            kinds.append((assay["assay"], kind))
            batch.append(text)
    answers = run_primer3("".join(batch))
    if len(answers) != len(batch):
        raise SystemExit(f"{PROGRAM} answered {len(answers)} of {len(batch)} request(s)")
    replies: dict[tuple[str, str], dict[str, str]] = dict(zip(kinds, answers, strict=True))

    def number(name: str, kind: str, key: str) -> float | None:
        raw = replies.get((name, kind), {}).get(key, "")
        try:
            return float(raw)
        except ValueError:
            return None

    def shown(value: float | None) -> str:
        return "-" if value is None else f"{value:.1f}"

    out = args.out or args.assays.with_name(args.assays.stem + "_thermo.tsv")
    columns = ["assay", "target", "forward_tm", "reverse_tm", "probe_tm", "tm_spread",
               "forward_gc", "reverse_gc", "probe_gc",
               "forward_hairpin", "reverse_hairpin", "probe_hairpin",
               "forward_self_dimer", "reverse_self_dimer", "probe_self_dimer",
               "pair_dimer", "pair_dimer_end", "over_limit"]
    with out.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for assay in assays:
            name = assay["assay"]
            values = {
                "forward_tm": number(name, "left", "PRIMER_LEFT_0_TM"),
                "reverse_tm": number(name, "right", "PRIMER_RIGHT_0_TM"),
                "probe_tm": number(name, "internal", "PRIMER_INTERNAL_0_TM"),
                "forward_gc": number(name, "left", "PRIMER_LEFT_0_GC_PERCENT"),
                "reverse_gc": number(name, "right", "PRIMER_RIGHT_0_GC_PERCENT"),
                "probe_gc": number(name, "internal", "PRIMER_INTERNAL_0_GC_PERCENT"),
                "forward_hairpin": number(name, "left", "PRIMER_LEFT_0_HAIRPIN_TH"),
                "reverse_hairpin": number(name, "right", "PRIMER_RIGHT_0_HAIRPIN_TH"),
                "probe_hairpin": number(name, "internal", "PRIMER_INTERNAL_0_HAIRPIN_TH"),
                "forward_self_dimer": number(name, "left", "PRIMER_LEFT_0_SELF_ANY_TH"),
                "reverse_self_dimer": number(name, "right", "PRIMER_RIGHT_0_SELF_ANY_TH"),
                "probe_self_dimer": number(name, "internal", "PRIMER_INTERNAL_0_SELF_ANY_TH"),
                "pair_dimer": number(name, "pair", "PRIMER_PAIR_0_COMPL_ANY_TH"),
                "pair_dimer_end": number(name, "pair", "PRIMER_PAIR_0_COMPL_END_TH"),
            }
            left, right = values["forward_tm"], values["reverse_tm"]
            spread = f"{abs(left - right):.1f}" if left is not None and right is not None else "-"
            over = [key for key in ("forward_hairpin", "reverse_hairpin", "probe_hairpin")
                    if values[key] is not None and values[key] >= limit_hairpin]
            over += [key for key in ("forward_self_dimer", "reverse_self_dimer", "probe_self_dimer",
                                     "pair_dimer", "pair_dimer_end")
                     if values[key] is not None and values[key] >= limit_dimer]
            fh.write("\t".join([name, assay.get("target", "-")]
                                + [shown(values[key]) for key in
                                   ("forward_tm", "reverse_tm", "probe_tm")]
                                + [spread]
                                + [shown(values[key]) for key in columns[6:17]]
                                + [";".join(over) or "none"]) + "\n")
    print(f"{out}: {len(assays)} assay(s) under {conditions}; "
          f"hairpin and dimer limits of the design defaults are {limit_hairpin} and {limit_dimer} degC, "
          f"and the GC clamp ({DEFAULT_GC_CLAMP}) is not applied here")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
