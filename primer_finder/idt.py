"""Turning an IDT order sheet into a fasta file of assays (forward primer, probe, reverse primer)."""

from __future__ import annotations

import codecs
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

# A number written with a grouping separator, which both the comma and the point can be: three digits
# after the separator is what tells "1,234" and "1.234" from the decimal "123.0", since an amplicon size
# has no thousandth of a base. A plain number, with or without a decimal part, is the second pattern.
GROUPED_NUMBER = re.compile(r"\d{1,3}(?:[.,]\d{3})+(?!\d)")
PLAIN_NUMBER = re.compile(r"\d+(?:[.,]\d+)?")

# What a base may be: the four, plus the IUPAC ambiguity codes and U for an RNA oligo.
BASES = frozenset("ACGTURYSWKMBDHVN")
# IDT writes a dye, a quencher or a modified base between slashes -- /56-FAM/, /3IABkFQ/, /i6-FAM/ -- and
# those are not bases. They are left out of the sequence, which is what an amplification check needs.
MODIFICATION = re.compile(r"/[^/]*/")


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

    for line, row in enumerate(rows, start=2):  # Row 1 of the sheet is the header
        kind = row.get(TYPE_COLUMN, "").strip().lower()
        sequence = row.get(SEQUENCE_COLUMN, "").strip()
        size = row.get(AMPLICON_COLUMN, "").strip()
        if kind in SUFFIXES and sequence:
            if kind in assay:  # The oligos of the next assay start here
                flush()
            assay[kind] = oligo_sequence(sequence, kind, line)
        if size:
            amplicon = _whole_number(size)
    flush()
    return records


def _whole_number(value: str) -> str:
    """"123", "123.0", "123 bp", "1,234" and "1.234" all give the number; anything else is kept as it is.

    The grouping separator matters: without it "1,234" matched only the "1" and the oligo was named after
    a 1 bp amplicon. "1.234" is the same 1,234 bases written by an Excel whose decimal separator is the
    comma -- the locale this converter already reads semicolon-separated CSV for -- and it used to give 1
    the same way, through `float("1.234")`. Three digits after the separator is what distinguishes them.
    """
    grouped = GROUPED_NUMBER.search(value)
    if grouped:
        return str(int(grouped.group().replace(",", "").replace(".", "")))
    plain = PLAIN_NUMBER.search(value)
    if not plain:
        return value
    return str(int(float(plain.group().replace(",", "."))))


def oligo_sequence(raw: str, kind: str, row: int) -> str:
    """The bases of one oligo as the fasta should hold them, or an error saying what the cell held.

    Nothing checked these: whatever the cell held was written into the fasta, and a sheet holding a note,
    a formula's error value or an order code in the sequence column produced a fasta that blast and
    insilicoPCR read as a sequence -- quietly, since neither refuses a record it cannot make sense of.
    """
    sequence = "".join(MODIFICATION.sub("", raw).split())  # Excel wraps a long sequence over lines
    bad = sorted({character for character in sequence.upper() if character not in BASES})
    if not sequence or bad:
        raise PrimerFinderError(
            f"Row {row} of the sheet gives the {kind} as {raw!r}, which is not a sequence"
            + (f": {', '.join(repr(character) for character in bad)} is not a base" if bad
               else " (it holds no base at all)")
            + ". An IDT modification between slashes, such as /56-FAM/, is not counted as one."
        )
    return sequence


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
    # The header is the first row that holds anything: Excel writes empty rows that carry only formatting.
    rows = [row for row in grid if any(cell.strip() for cell in row)]
    if not rows:
        raise PrimerFinderError(f"{path} is empty")
    header = [cell.strip().lower() for cell in rows[0]]
    return [dict(zip(header, row, strict=False)) for row in rows[1:]]


def decode(data: bytes) -> str:
    """Text out of whatever Excel saved: UTF-8 with or without a byte-order mark ("CSV UTF-8"), UTF-16
    ("Unicode Text"), or the legacy Windows code page."""
    for mark, encoding in ((codecs.BOM_UTF16_LE, "utf-16"), (codecs.BOM_UTF16_BE, "utf-16")):
        if data.startswith(mark):
            return data.decode(encoding)
    for encoding in ("utf-8-sig", "cp1252"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    return data.decode("utf-8", errors="replace")  # pragma: no cover - cp1252 decodes any byte


def read_delimited(path: Path) -> list[list[str]]:
    """Read a delimited text file, guessing the delimiter from the first line: a tab, a comma, or the
    semicolon that Excel writes where the comma is the decimal separator."""
    text = decode(path.read_bytes())
    first = text.splitlines()[0] if text else ""
    delimiter = max(("\t", ",", ";"), key=first.count)
    return [list(row) for row in csv.reader(text.splitlines(), delimiter=delimiter)]


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
    return [_shared_string(item) for item in root.iterfind(".//{*}si")]


def _shared_string(item: ElementTree.Element) -> str:
    """One shared string: its text, or the text of its rich-text runs. The phonetic runs (`rPh`) that East
    Asian versions of Excel add alongside are not part of the value."""
    parts: list[str] = []
    for child in item:
        tag = child.tag.rsplit("}", 1)[-1]
        if tag == "t":
            parts.append(child.text or "")
        elif tag == "r":  # A run of rich text: <r><rPr/><t>...</t></r>
            parts.extend(text.text or "" for text in child.iterfind("{*}t"))
    return "".join(parts)


def _sheet_part(archive: zipfile.ZipFile, sheet: str | None) -> str:
    """The path inside the archive of the wanted worksheet, or of the first one."""
    with archive.open("xl/workbook.xml") as fh:
        workbook = ElementTree.parse(fh).getroot()
    relations = _relations(archive)
    sheets = list(workbook.iterfind(".//{*}sheet"))
    if sheet is None and len(sheets) > 1:
        names = [element.get("name", "") for element in sheets]
        log.warning("The file holds several sheets (%s); reading the first one, %s",
                    ", ".join(names), names[0])
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
