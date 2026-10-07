"""Turning an IDT order sheet into a fasta file of assays (forward primer, probe, reverse primer)."""

from __future__ import annotations

import csv
import logging
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree

from primer_finder import PrimerFinderError
from primer_finder.seqio import Record, write_fasta

log = logging.getLogger(__name__)

# The oligo types IDT uses, and the suffix each one gets in the fasta header.
SUFFIXES = {"forward primer": "-F", "probe": "-P", "reverse primer": "-R"}
# The columns the sheet must have (matched without regard to case or surrounding spaces).
TYPE_COLUMN, SEQUENCE_COLUMN, AMPLICON_COLUMN = "type", "sequence", "amplicon"

CELL = re.compile(r"([A-Z]+)(\d+)")


def convert(table: Path, output: Path, prefix: str = "") -> int:
    """Write the assays of an IDT sheet (.xlsx, .csv or .tsv) as a fasta file. Returns how many assays."""
    rows = read_table(table)
    records = to_records(rows, prefix)
    if not records:
        raise PrimerFinderError(
            f"No assay found in {table}. The sheet needs a \"Type\" column holding "
            f"{', '.join(sorted(SUFFIXES))} and a \"Sequence\" column."
        )
    write_fasta(output, records)
    assays = len({record.name.rsplit("-", 1)[0] for record in records})
    log.info("Wrote %d oligo(s) for %d assay(s) to %s", len(records), assays, output)
    return assays


def to_records(rows: list[dict[str, str]], prefix: str = "") -> list[Record]:
    """One record per oligo, named after the assay it belongs to and the amplicon size.

    IDT lists the oligos of an assay together (forward primer, probe, reverse primer), and the amplicon size
    may sit on any of those rows or on a summary row below them. An assay therefore ends when an oligo type
    comes round again, or at the end of the sheet.
    """
    records: list[Record] = []
    assay: dict[str, str] = {}
    amplicon = ""
    number = 0

    def flush() -> None:
        nonlocal assay, amplicon, number
        if not assay:
            amplicon = ""
            return
        name = f"{prefix}_{number}" if prefix else str(number)
        size = f"_{amplicon}bp" if amplicon else ""
        for kind in ("forward primer", "probe", "reverse primer"):
            if kind in assay:
                records.append(Record(f"{name}{size}{SUFFIXES[kind]}", "", assay[kind]))
        number += 1
        assay = {}
        amplicon = ""

    for row in rows:
        kind = row.get(TYPE_COLUMN, "").strip().lower()
        sequence = row.get(SEQUENCE_COLUMN, "").strip()
        size = row.get(AMPLICON_COLUMN, "").strip()
        if kind in SUFFIXES and sequence:
            if kind in assay:  # The oligos of the next assay start here
                flush()
            assay[kind] = sequence
        if size:
            amplicon = _whole_number(size)
    flush()
    return records


def _whole_number(value: str) -> str:
    """"123", "123.0" and "123 bp" all give "123"; anything else is kept as it is."""
    match = re.search(r"\d+(?:\.\d+)?", value)
    if not match:
        return value
    return str(int(float(match.group())))


def read_table(path: Path) -> list[dict[str, str]]:
    """Read a spreadsheet or a delimited text file as a list of rows keyed by lower-case column name."""
    if not path.is_file():
        raise PrimerFinderError(f"No such file: {path}")
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        grid = read_xlsx(path)
    elif path.suffix.lower() in (".csv", ".tsv", ".txt"):
        grid = read_delimited(path)
    else:
        raise PrimerFinderError(
            f"Cannot read {path}: give an .xlsx, .csv or .tsv file (in Excel: File, Save as, CSV)."
        )
    if not grid:
        raise PrimerFinderError(f"{path} is empty")
    header = [cell.strip().lower() for cell in grid[0]]
    return [dict(zip(header, row, strict=False)) for row in grid[1:] if any(cell.strip() for cell in row)]


def read_delimited(path: Path) -> list[list[str]]:
    """Read a csv or tab-separated file, guessing which of the two it is from the first line."""
    with path.open(newline="") as fh:
        first = fh.readline()
        fh.seek(0)
        delimiter = "\t" if first.count("\t") > first.count(",") else ","
        return [list(row) for row in csv.reader(fh, delimiter=delimiter)]


def read_xlsx(path: Path, sheet: str | None = None) -> list[list[str]]:
    """Read the cells of one worksheet of an .xlsx file, without any third-party library.

    An .xlsx file is a zip of XML parts: the shared strings, the workbook (which names the sheets) and one
    part per worksheet. Only what this converter needs is read: cell values, as text.
    """
    try:
        with zipfile.ZipFile(path) as archive:
            strings = _shared_strings(archive)
            part = _sheet_part(archive, sheet)
            with archive.open(part) as fh:
                return _cells(ElementTree.parse(fh).getroot(), strings)
    except (zipfile.BadZipFile, KeyError, OSError, ElementTree.ParseError) as exc:
        raise PrimerFinderError(f"Could not read {path} as an Excel file: {exc}") from exc


def _shared_strings(archive: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    with archive.open("xl/sharedStrings.xml") as fh:
        root = ElementTree.parse(fh).getroot()
    return ["".join(text.text or "" for text in item.iterfind(".//{*}t")) for item in root.iterfind(".//{*}si")]


def _sheet_part(archive: zipfile.ZipFile, sheet: str | None) -> str:
    """The path inside the archive of the wanted worksheet, or of the first one."""
    with archive.open("xl/workbook.xml") as fh:
        workbook = ElementTree.parse(fh).getroot()
    relations = _relations(archive)
    sheets = list(workbook.iterfind(".//{*}sheet"))
    for element in sheets:
        name = element.get("name", "")
        if sheet is None or name == sheet:
            target = relations.get(_relation_id(element))
            if target:
                return target
            break
    if sheet is not None:
        raise PrimerFinderError(f"No sheet named {sheet!r} in the file")
    parts = [name for name in archive.namelist() if name.startswith("xl/worksheets/") and name.endswith(".xml")]
    if not parts:
        raise PrimerFinderError("The Excel file holds no worksheet")
    return sorted(parts)[0]


def _relation_id(element: ElementTree.Element) -> str:
    for key, value in element.attrib.items():
        if key.endswith("}id") or key == "id":
            return value
    return ""


def _relations(archive: zipfile.ZipFile) -> dict[str, str]:
    if "xl/_rels/workbook.xml.rels" not in archive.namelist():
        return {}
    with archive.open("xl/_rels/workbook.xml.rels") as fh:
        root = ElementTree.parse(fh).getroot()
    return {
        element.get("Id", ""): _inside_xl(element.get("Target", ""))
        for element in root.iterfind(".//{*}Relationship")
    }


def _inside_xl(target: str) -> str:
    """A relationship target as a path in the archive: it may be absolute ("/xl/worksheets/sheet1.xml",
    which is what Excel and openpyxl write) or relative to the xl folder ("worksheets/sheet1.xml")."""
    target = target.lstrip("/")
    return target if target.startswith("xl/") else "xl/" + target


def _cells(sheet: ElementTree.Element, strings: list[str]) -> list[list[str]]:
    """The worksheet as a list of rows of text, with the blank cells kept in place."""
    grid: list[list[str]] = []
    for row in sheet.iterfind(".//{*}row"):
        values: list[str] = []
        for cell in row.iterfind("{*}c"):
            index = _column_index(cell.get("r", ""), len(values))
            values.extend([""] * (index - len(values)))
            values.append(_value(cell, strings))
        grid.append(values)
    return grid


def _column_index(reference: str, fallback: int) -> int:
    """The 0-based column of a cell reference such as "B3"."""
    match = CELL.match(reference)
    if not match:
        return fallback
    index = 0
    for character in match.group(1):
        index = index * 26 + (ord(character) - ord("A") + 1)
    return index - 1


def _value(cell: ElementTree.Element, strings: list[str]) -> str:
    kind = cell.get("t", "n")
    if kind == "inlineStr":
        return "".join(text.text or "" for text in cell.iterfind(".//{*}t"))
    value = cell.find("{*}v")
    text = value.text if value is not None and value.text is not None else ""
    if kind == "s":
        position = int(text) if text.isdigit() else -1
        return strings[position] if 0 <= position < len(strings) else ""
    if kind == "b":
        return "TRUE" if text == "1" else "FALSE"
    return text
