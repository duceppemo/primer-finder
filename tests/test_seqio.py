"""Reading and writing fasta files."""

from __future__ import annotations

import gzip

import pytest

from primer_finder import PrimerFinderError
from primer_finder.seqio import (
    MOST_LISTED,
    Record,
    base_name,
    count_records,
    find_genomes,
    is_fasta,
    iter_records,
    read_fasta,
    require_genomes,
    write_fasta,
)


def test_iter_records_splits_the_header_and_joins_the_lines(fasta):
    path = fasta("x.fasta", {})
    path.write_text(">one some description\nACGT\nTTTT\n\n>two\nGG\n")
    records = list(iter_records(path))
    assert [(record.name, record.desc, record.seq) for record in records] == [
        ("one", "some description", "ACGTTTTT"),
        ("two", "", "GG"),
    ]
    assert records[0].header == "one some description"
    assert count_records(path) == 2


def test_an_empty_file_holds_no_record(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    assert list(iter_records(path)) == []
    assert count_records(path) == 0
    assert not is_fasta(path)


def test_gzipped_files_are_read(tmp_path):
    path = tmp_path / "g.fasta.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">one\nACGT\n")
    assert read_fasta(path)["one"].seq == "ACGT"
    assert is_fasta(path)
    assert base_name(path) == "g"


def test_a_file_that_is_not_gzipped_but_named_gz(tmp_path):
    path = tmp_path / "broken.fasta.gz"
    path.write_text(">one\nACGT\n")
    assert not is_fasta(path)


def test_write_fasta(tmp_path):
    path = tmp_path / "out.fasta"
    assert write_fasta(path, [Record("a", "note", "ACGT"), Record("b", "", "TT")]) == 2
    assert path.read_text() == ">a note\nACGT\n>b\nTT\n"


def test_find_genomes_is_recursive_sorted_and_follows_links(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "b.fasta").write_text(">b\nA\n")
    (tmp_path / "sub" / "a.fna").write_text(">a\nA\n")
    (tmp_path / "notes.txt").write_text("ignored")
    (tmp_path / "c.fasta").symlink_to(tmp_path / "b.fasta")
    assert [path.name for path in find_genomes(tmp_path)] == ["b.fasta", "c.fasta", "a.fna"]


def test_require_genomes_reports_what_is_wrong(tmp_path):
    with pytest.raises(PrimerFinderError, match="inclusion folder does not exist"):
        require_genomes(tmp_path / "gone", "inclusion")
    (tmp_path / "empty").mkdir()
    with pytest.raises(PrimerFinderError, match="no file with an accepted extension"):
        require_genomes(tmp_path / "empty", "inclusion")
    (tmp_path / "empty" / "a.fasta").write_text("not a fasta\n")
    with pytest.raises(PrimerFinderError, match="1 file.* are not fasta"):
        require_genomes(tmp_path / "empty", "inclusion")


def test_base_name():
    from pathlib import Path

    assert base_name(Path("/a/b/genome.fasta")) == "genome"
    assert base_name(Path("genome.fna.gz")) == "genome"
    assert base_name(Path("genome.v2.fasta")) == "genome.v2"


def test_a_header_that_is_not_utf8_is_read_anyway(tmp_path, caplog):
    """Older pipelines wrote accented strain names; a header is only a label here."""
    path = tmp_path / "latin1.fasta"
    path.write_bytes(b">Mycobacterium bovis caf\xe9 strain\nACGT\n")
    assert is_fasta(path)
    record = next(iter(iter_records(path)))
    assert record.name == "Mycobacterium"
    assert record.seq == "ACGT"


def test_gzip_data_in_a_file_that_is_not_named_gz_is_not_a_fasta(tmp_path):
    import gzip as gziplib

    path = tmp_path / "mislabelled.fasta"
    with gziplib.open(path, "wt") as fh:
        fh.write(">chr\nACGT\n")
    assert not is_fasta(path)
    with pytest.raises(PrimerFinderError, match="are not fasta"):
        require_genomes(tmp_path, "inclusion")


def test_duplicate_record_names_are_reported(tmp_path, caplog):
    path = tmp_path / "dup.fasta"
    path.write_text(">ctg1\nAAAA\n>ctg1\nTTTT\n>ctg2\nGG\n")
    records = read_fasta(path)
    assert set(records) == {"ctg1", "ctg2"}
    assert records["ctg1"].seq == "TTTT"
    assert "name is not unique (ctg1)" in caplog.text
    assert count_records(path) == 3


def test_a_folder_full_of_files_that_are_not_fasta_gives_a_short_error(tmp_path):
    """Pointing -i at a large folder by mistake must not print a thousand paths."""
    for number in range(40):
        (tmp_path / f"g{number:02d}.fasta").write_text("not a fasta\n")
    with pytest.raises(PrimerFinderError) as error:
        require_genomes(tmp_path, "inclusion")
    message = str(error.value)
    assert "40 file(s)" in message
    assert "and 35 more" in message
    # The length of the message depends on how long the paths are, so count the names it lists
    assert message.count(".fasta") == MOST_LISTED
