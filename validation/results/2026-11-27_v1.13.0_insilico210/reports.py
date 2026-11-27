"""Reading the reports of both arms, and saying why an (assay, genome) pair moved between them.

primer-finder's own reader refuses the reports of insilicoPCR before 2.0, as it should, so the old arm's
(0.14.0: `Gene`, no Tm, no `Pairing`) are read here, by a reader that takes both. The new arm's reports of
the assays without a product (`report_no_product.tsv`, 2.0+) say why a pair has no product there; the old
arm's reports, and the new arm's Tm columns, say what the product rested on where there is one.
"""
from __future__ import annotations

import contextlib
import csv
import re
from collections import Counter
from pathlib import Path

REPORT = Path("consolidated_report") / "report.tsv"
Pair = tuple[str, str]


def _cell(value: str | None) -> str:
    """As primer_finder.insilico._cell, enough of it for the sample and assay names these runs hold."""
    text = (value or "").strip()
    if len(text) >= 2 and text[0] == text[-1] == '"':
        text = text[1:-1].replace('""', '"')
    if text.startswith("'") and len(text) > 1 and text[1] in "=@+-'":
        text = text[1:]
    return text


def rows(path: Path) -> list[dict[str, str]]:
    """A report's rows, the assay under `Assay` whichever release wrote it."""
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as fh:
        found = [{key: _cell(value) for key, value in row.items() if key is not None}
                 for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)]
    for row in found:
        row.setdefault("Assay", row.get("Gene", ""))
    return found


def products(folder: Path, designed_only: bool = False) -> dict[Pair, list[dict[str, str]]]:
    """A run's products by (assay, sample). Before 1.0 every product is the designed pair."""
    found: dict[Pair, list[dict[str, str]]] = {}
    for row in rows(folder / REPORT):
        if designed_only and row.get("Pairing", "designed") != "designed":
            continue
        found.setdefault((row["Assay"], row["Sample"]), []).append(row)
    return found


def no_product(folder: Path) -> dict[Pair, list[dict[str, str]]]:
    found: dict[Pair, list[dict[str, str]]] = {}
    for row in rows(folder / REPORT.with_name("report_no_product.tsv")):
        found.setdefault((row["Assay"], row["Sample"]), []).append(row)
    return found


def role(oligo: str) -> str:
    return {"F": "forward", "R": "reverse", "P": "probe"}.get(oligo.rpartition("-")[2][:1], oligo or "?")


def why_none(reasons: list[dict[str, str]]) -> str:
    """What insilicoPCR 2.x says of a pair without a product, as a category: per failing oligo, the rule."""
    if not reasons:
        return "no reason given"
    said = []
    for row in reasons:
        text = row.get("Reason", "")
        if not row.get("Oligo"):
            what = "no forward and reverse site facing each other"
        elif "perfect match is" in text:
            what = f"{role(row['Oligo'])}: below the needed Tm even on a perfect match"
        elif text.startswith("no site reaches"):
            what = f"{role(row['Oligo'])}: no site melts at the needed Tm"
        elif "fails the 3' rule" in text:
            what = f"{role(row['Oligo'])}: its best site fails the 3' rule"
        elif "BLAST found no match" in text:
            what = f"{role(row['Oligo'])}: no BLAST hit"
        else:
            what = f"{role(row.get('Oligo', ''))}: {text}"
        if "one of its primers alone" in text:
            what += " (a product of one primer alone instead)"
        said.append(what)
    return "; ".join(sorted(set(said)))


def best_site(reasons: list[dict[str, str]]) -> str:
    """The Tm of the failing oligos' best sites, for an example line."""
    return ", ".join(f"{role(row['Oligo'])} {row['BestSiteTm']} C" for row in reasons
                     if row.get("Oligo") and row.get("BestSiteTm"))


def _mismatches(row: dict[str, str]) -> str:
    probe = f" P{row['ProbeMismatches']}" if row.get("Probe") else ""
    return f"F{row.get('ForwardMismatches', '?')} R{row.get('ReverseMismatches', '?')}{probe}"


def weakest(row: dict[str, str]) -> float | None:
    tms = []
    for column in ("ForwardTm", "ReverseTm", "ProbeTm"):
        with contextlib.suppress(ValueError):
            tms.append(float(row.get(column, "")))
    return min(tms) if tms else None


def what_made(found: list[dict[str, str]]) -> str:
    """What a pair's best product rests on, as a category: how it paired, and the most mismatches its
    primers and its probe carry (the old arm's budget counted them; the new arm's Tm weighs them)."""
    best = max(found, key=lambda row: weakest(row) or 0.0)
    pairing = best.get("Pairing", "designed")
    primers = max(int(best.get("ForwardMismatches") or 0), int(best.get("ReverseMismatches") or 0))
    text = f"{pairing}, a primer with {primers} mismatch(es)"
    if best.get("Probe"):
        text += f", probe with {best.get('ProbeMismatches') or 0}"
    return text


def example(found: list[dict[str, str]]) -> str:
    best = max(found, key=lambda row: weakest(row) or 0.0)
    tm = weakest(best)
    return _mismatches(best) + (f", weakest site {tm:.1f} C" if tm is not None else "")


def tally(lines: Counter[str], examples: dict[str, list[str]], indent: str = "      ") -> None:
    for line, count in lines.most_common():
        print(f"{indent}{count:4d} x {line}")
        for shown in examples.get(line, [])[:2]:
            print(f"{indent}       e.g. {shown}")


def needed_tm(log: Path) -> str:
    match = re.search(r"reaches ([\d.]+) C \(--anneal-temp", log.read_text()) if log.is_file() else None
    return match.group(1) if match else "?"
