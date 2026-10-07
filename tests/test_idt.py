"""Converting an IDT order sheet into a fasta file."""

from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.idt import convert, read_table, to_records

ROWS = [
    {"type": "Forward Primer", "sequence": "AAAACCCC", "amplicon": ""},
    {"type": "Probe", "sequence": "TTTTGGGG", "amplicon": ""},
    {"type": "Reverse Primer", "sequence": "CCCCAAAA", "amplicon": ""},
    {"type": "", "sequence": "", "amplicon": "120.0"},
]


def test_one_assay_becomes_three_records():
    assert [(record.name, record.seq) for record in to_records(ROWS)] == [
        ("0_120bp-F", "AAAACCCC"),
        ("0_120bp-P", "TTTTGGGG"),
        ("0_120bp-R", "CCCCAAAA"),
    ]


def test_the_prefix_goes_in_front_of_the_assay_number():
    assert to_records(ROWS, "mytarget")[0].name == "mytarget_0_120bp-F"


def test_a_second_assay_starts_when_a_type_comes_round_again():
    rows = ROWS + [
        {"type": "Forward Primer", "sequence": "GGGG", "amplicon": "95"},
        {"type": "Reverse Primer", "sequence": "TTTT", "amplicon": ""},
    ]
    names = [record.name for record in to_records(rows)]
    assert names == ["0_120bp-F", "0_120bp-P", "0_120bp-R", "1_95bp-F", "1_95bp-R"]


def test_an_assay_without_an_amplicon_size():
    rows = [row | {"amplicon": ""} for row in ROWS]
    assert [record.name for record in to_records(rows)] == ["0-F", "0-P", "0-R"]


def test_rows_that_are_not_oligos_are_ignored():
    rows = [{"type": "Plate", "sequence": "ACGT", "amplicon": ""}, *ROWS]
    assert len(to_records(rows)) == 3


MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
RELS_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACKAGE_RELS_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def write_xlsx(path: Path, grid: list[list[str]], namespaces: bool = True) -> Path:
    """An .xlsx holding one worksheet, as Excel writes it: namespaced XML, shared strings and an absolute
    relationship target. With `namespaces` off, the same file without any namespace, which some exporters
    produce. openpyxl instead writes inline strings and no shared-strings part, which
    test_inline_strings covers.
    """
    strings: list[str] = []
    rows_xml = []
    for row_number, row in enumerate(grid, start=1):
        cells = []
        for column, value in enumerate(row):
            letter = chr(ord("A") + column)
            if value == "":
                continue
            if value.replace(".", "", 1).isdigit():
                cells.append(f'<c r="{letter}{row_number}"><v>{value}</v></c>')
            else:
                if value not in strings:
                    strings.append(value)
                cells.append(f'<c r="{letter}{row_number}" t="s"><v>{strings.index(value)}</v></c>')
        rows_xml.append(f'<row r="{row_number}">{"".join(cells)}</row>')
    shared = "".join(f"<si><t>{text}</t></si>" for text in strings)
    main = f' xmlns="{MAIN_NS}"' if namespaces else ""
    rels_attr = f' xmlns:r="{RELS_NS}"' if namespaces else ""
    package_attr = f' xmlns="{PACKAGE_RELS_NS}"' if namespaces else ""
    id_attr = "r:id" if namespaces else "id"
    target = "/xl/worksheets/sheet1.xml" if namespaces else "worksheets/sheet1.xml"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("[Content_Types].xml", "<Types/>")
        archive.writestr(
            "xl/workbook.xml",
            f"<workbook{main}{rels_attr}><sheets>"
            f'<sheet name="Sheet1" sheetId="1" {id_attr}="rId1"/></sheets></workbook>',
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            f"<Relationships{package_attr}>"
            f'<Relationship Id="rId1" Target="{target}"/></Relationships>',
        )
        archive.writestr("xl/sharedStrings.xml", f"<sst{main}>{shared}</sst>")
        archive.writestr(
            "xl/worksheets/sheet1.xml",
            f'<worksheet{main}><sheetData>{"".join(rows_xml)}</sheetData></worksheet>',
        )
    return path


GRID = [
    ["Type", "Sequence", "Amplicon"],
    ["Forward Primer", "AAAACCCC", ""],
    ["Probe", "TTTTGGGG", ""],
    ["Reverse Primer", "CCCCAAAA", ""],
    ["", "", "120"],
]


def test_convert_an_xlsx_file(tmp_path):
    table = write_xlsx(tmp_path / "order.xlsx", GRID)
    output = tmp_path / "assays.fasta"
    assert convert(table, output, "target") == 1
    assert output.read_text() == (
        ">target_0_120bp-F\nAAAACCCC\n>target_0_120bp-P\nTTTTGGGG\n>target_0_120bp-R\nCCCCAAAA\n"
    )


def test_the_columns_may_be_in_any_order_and_any_case(tmp_path):
    grid = [["sequence", "AMPLICON", "Type"], ["AAAA", "75", "forward primer"], ["TTTT", "", "Reverse primer"]]
    table = write_xlsx(tmp_path / "order.xlsx", grid)
    assert convert(table, tmp_path / "out.fasta") == 1
    assert ">0_75bp-F" in (tmp_path / "out.fasta").read_text()


@pytest.mark.parametrize("suffix,separator", [(".csv", ","), (".tsv", "\t")])
def test_convert_a_delimited_file(tmp_path, suffix, separator):
    table = tmp_path / f"order{suffix}"
    table.write_text("".join(separator.join(row) + "\n" for row in GRID))
    assert convert(table, tmp_path / "out.fasta") == 1
    assert ">0_120bp-P" in (tmp_path / "out.fasta").read_text()


def test_an_unknown_file_type(tmp_path):
    table = tmp_path / "order.xls"
    table.write_text("x")
    with pytest.raises(PrimerFinderError, match="give an .xlsx, .csv or .tsv file"):
        convert(table, tmp_path / "out.fasta")


def test_a_sheet_without_assays(tmp_path):
    table = write_xlsx(tmp_path / "order.xlsx", [["Type", "Sequence"], ["Plate", "ACGT"]])
    with pytest.raises(PrimerFinderError, match="No assay found"):
        convert(table, tmp_path / "out.fasta")


def test_an_empty_sheet(tmp_path):
    table = write_xlsx(tmp_path / "order.xlsx", [])
    with pytest.raises(PrimerFinderError, match="is empty"):
        read_table(table)


def test_a_file_that_is_not_a_spreadsheet(tmp_path):
    table = tmp_path / "order.xlsx"
    table.write_text("not a zip file")
    with pytest.raises(PrimerFinderError, match="Could not read"):
        read_table(table)


def test_blank_cells_keep_their_column(tmp_path):
    """A row whose first cells are empty must still line up with the header."""
    grid = [["Type", "Sequence", "Amplicon"], ["", "", "95"], ["Forward Primer", "AAAA", ""],
            ["Reverse Primer", "TTTT", ""]]
    table = write_xlsx(tmp_path / "order.xlsx", grid)
    rows = read_table(table)
    assert rows[0] == {"type": "", "sequence": "", "amplicon": "95"}
    assert rows[1]["type"] == "Forward Primer"


def test_a_sheet_without_namespaces(tmp_path):
    """Not every exporter writes the namespaces Excel does."""
    table = write_xlsx(tmp_path / "order.xlsx", GRID, namespaces=False)
    assert convert(table, tmp_path / "out.fasta") == 1
    assert ">0_120bp-R" in (tmp_path / "out.fasta").read_text()


def test_a_sheet_can_be_chosen_by_name(tmp_path):
    from primer_finder.idt import read_xlsx

    table = write_xlsx(tmp_path / "order.xlsx", GRID)
    assert read_xlsx(table, "Sheet1")[0] == ["Type", "Sequence", "Amplicon"]
    with pytest.raises(PrimerFinderError, match="No sheet named"):
        read_xlsx(table, "Nope")


def write_parts(path: Path, parts: dict[str, str], sheet_name: str = "Sheet1") -> Path:
    """An .xlsx whose worksheet and shared strings are given as raw XML, for the shapes a spreadsheet
    program produces that write_xlsx does not."""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("[Content_Types].xml", "<Types/>")
        archive.writestr(
            "xl/workbook.xml",
            f'<workbook xmlns="{MAIN_NS}" xmlns:r="{RELS_NS}"><sheets>'
            f'<sheet name="{sheet_name}" sheetId="1" r:id="rId1"/></sheets></workbook>',
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            f'<Relationships xmlns="{PACKAGE_RELS_NS}">'
            '<Relationship Id="rId1" Target="/xl/worksheets/sheet1.xml"/></Relationships>',
        )
        for name, xml in parts.items():
            archive.writestr(name, xml)
    return path


def sheet_xml(rows: str) -> str:
    return f'<worksheet xmlns="{MAIN_NS}"><sheetData>{rows}</sheetData></worksheet>'


ASSAY_CSV = "Type,Sequence,Amplicon\nForward Primer,AAAA,95\nReverse Primer,TTTT,\n"


@pytest.mark.parametrize(
    ("name", "data"),
    [
        # Excel's "CSV UTF-8", the default CSV format, starts with a byte-order mark
        ("bom.csv", b"\xef\xbb\xbf" + ASSAY_CSV.encode()),
        # Excel writes semicolons where the comma is the decimal separator (fr, de, ...)
        ("semicolon.csv", ASSAY_CSV.replace(",", ";").encode()),
        # The legacy "CSV (Comma delimited)" on Windows is cp1252, here with a degree sign
        ("cp1252.csv", ASSAY_CSV.replace("95", "95,60\xb0C").encode("cp1252")),
        # Excel's "Unicode Text (*.txt)" is UTF-16 with tabs
        ("unicode.txt", ASSAY_CSV.replace(",", "\t").encode("utf-16")),
    ],
)
def test_the_text_formats_excel_saves(tmp_path, name, data):
    table = tmp_path / name
    table.write_bytes(data)
    assert convert(table, tmp_path / "out.fasta") == 1
    assert (tmp_path / "out.fasta").read_text() == ">0_95bp-F\nAAAA\n>0_95bp-R\nTTTT\n"


def test_inline_strings(tmp_path):
    """openpyxl writes the text inside the cell and no shared-strings part."""
    def row(number, values):
        cells = "".join(
            f'<c r="{chr(ord("A") + column)}{number}" t="inlineStr"><is><t>{value}</t></is></c>'
            for column, value in enumerate(values)
        )
        return f'<row r="{number}">{cells}</row>'

    table = write_parts(tmp_path / "order.xlsx", {"xl/worksheets/sheet1.xml": sheet_xml(
        row(1, ["Type", "Sequence", "Amplicon"]) + row(2, ["Forward Primer", "AAAA", "95"])
        + row(3, ["Reverse Primer", "TTTT"])
    )})
    assert convert(table, tmp_path / "out.fasta") == 1


def test_phonetic_runs_are_not_part_of_the_value(tmp_path):
    """East Asian versions of Excel store a phonetic run (rPh) beside the text; it is not the value."""
    shared = (
        f'<sst xmlns="{MAIN_NS}">'
        '<si><r><t>Forward </t></r><r><t>Primer</t></r><rPh sb="0" eb="7"><t>FW</t></rPh></si>'
        "<si><t>Sequence</t></si><si><t>Type</t></si>"
        "<si><r><t>Reverse Primer</t></r><rPh sb=\"0\" eb=\"7\"><t>RV</t></rPh></si>"
        "</sst>"
    )
    rows = (
        '<row r="1"><c r="A1" t="s"><v>2</v></c><c r="B1" t="s"><v>1</v></c></row>'
        '<row r="2"><c r="A2" t="s"><v>0</v></c><c r="B2" t="s"><v>1</v></c></row>'
        '<row r="3"><c r="A3" t="s"><v>3</v></c><c r="B3" t="s"><v>1</v></c></row>'
    )
    table = write_parts(tmp_path / "order.xlsx",
                        {"xl/sharedStrings.xml": shared, "xl/worksheets/sheet1.xml": sheet_xml(rows)})
    rows_read = read_table(table)
    assert [row["type"] for row in rows_read] == ["Forward Primer", "Reverse Primer"]


def test_a_first_row_that_only_carries_formatting_is_not_the_header(tmp_path):
    """Excel writes a styled but empty row as cells with no value at all."""
    grid = [["", "", ""], ["Type", "Sequence", "Amplicon"], ["Forward Primer", "AAAA", "95"],
            ["Reverse Primer", "TTTT", ""]]
    table = write_xlsx(tmp_path / "order.xlsx", grid)
    assert read_table(table)[0]["type"] == "Forward Primer"
    assert convert(table, tmp_path / "out.fasta") == 1


def test_columns_past_z(tmp_path):
    rows = (
        '<row r="1"><c r="AA1" t="inlineStr"><is><t>Type</t></is></c>'
        '<c r="AB1" t="inlineStr"><is><t>Sequence</t></is></c></row>'
        '<row r="2"><c r="AA2" t="inlineStr"><is><t>Forward Primer</t></is></c>'
        '<c r="AB2" t="inlineStr"><is><t>AAAA</t></is></c></row>'
        '<row r="3"><c r="AA3" t="inlineStr"><is><t>Reverse Primer</t></is></c>'
        '<c r="AB3" t="inlineStr"><is><t>TTTT</t></is></c></row>'
    )
    table = write_parts(tmp_path / "order.xlsx", {"xl/worksheets/sheet1.xml": sheet_xml(rows)})
    assert convert(table, tmp_path / "out.fasta") == 1


def test_an_xlsm_file_is_read(tmp_path):
    table = write_xlsx(tmp_path / "order.xlsm", GRID)
    assert convert(table, tmp_path / "out.fasta") == 1


def test_only_the_first_sheet_is_read_and_it_is_said(tmp_path, caplog):
    import zipfile as zf

    table = write_xlsx(tmp_path / "order.xlsx", GRID)
    # Add a second sheet, so that the workbook lists two
    with zf.ZipFile(table) as archive:
        parts = {name: archive.read(name) for name in archive.namelist()}
    parts["xl/workbook.xml"] = parts["xl/workbook.xml"].replace(
        b"</sheets>", b'<sheet name="Notes" sheetId="2" r:id="rId2"/></sheets>')
    with zf.ZipFile(table, "w") as archive:
        for name, data in parts.items():
            archive.writestr(name, data)
    assert convert(table, tmp_path / "out.fasta") == 1
    assert "several sheets (Sheet1, Notes); reading the first one" in caplog.text
