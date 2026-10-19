"""Reading and writing fasta files."""

from __future__ import annotations

import bz2
import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.seqio import (
    MOST_LISTED,
    Record,
    base_name,
    compressed_as,
    count_records,
    find_genomes,
    is_fasta,
    iter_records,
    open_text,
    read_fasta,
    require_genomes,
    walk_genomes,
    write_fasta,
)


def test_a_truncated_gzip_says_which_file_it_is(tmp_path):
    """The decompressor raises EOFError, which names nothing; a genome folder holds hundreds of files, so
    the name is the whole of the useful information."""
    good = tmp_path / "whole.fasta.gz"
    with gzip.open(good, "wt") as fh:
        fh.write(">chr\n" + "ACGT" * 2000 + "\n")
    cut = tmp_path / "cut.fasta.gz"
    cut.write_bytes(good.read_bytes()[: len(good.read_bytes()) // 2])

    with pytest.raises(PrimerFinderError, match="ends before its data does"):
        list(iter_records(cut))
    with pytest.raises(PrimerFinderError, match="cut.fasta.gz"):
        count_records(cut)
    assert len(list(iter_records(good))) == 1   # and a whole one still reads


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
    # A link to a folder, which the walk descends -- so it finds what is in it a second time
    (tmp_path / "d.fasta").symlink_to(tmp_path / "sub")
    assert [path.name for path in walk_genomes(tmp_path)] == ["b.fasta", "a.fna", "a.fna"]
    assert [path.name for path in find_genomes(tmp_path)] == ["b.fasta", "a.fna"]


def test_find_genomes_counts_a_genome_once_however_many_names_it_has(tmp_path):
    """The staging and the genome counts the in-silico verdicts are measured against call this directly,
    so a second name for a genome here makes the denominator of those verdicts larger than the panel the
    kmers were counted from."""
    (tmp_path / "b.fasta").write_text(">b\nA\n")
    (tmp_path / "c.fasta").symlink_to(tmp_path / "b.fasta")
    assert [path.name for path in walk_genomes(tmp_path)] == ["b.fasta", "c.fasta"]
    assert [path.name for path in find_genomes(tmp_path)] == ["b.fasta"]


def test_a_compressed_genome_with_an_unaccepted_name_is_reported_not_dropped(tmp_path, caplog):
    """The dangerous direction. A genome left out of the exclusion group has its kmers read as specific to
    the inclusion group, and a compressed file used to leave in total silence: `is_fasta` opened it as
    text, found no `>`, and the "looks like fasta" warning never named it. `.GZ` comes from Windows,
    `.bz2` and `.zip` from older archives."""
    folder = tmp_path / "exclusion"
    folder.mkdir()
    for name in ("g1.fasta", "g2.fasta"):
        (folder / name).write_text(">chr\nACGT\n")
    with gzip.open(folder / "g3.fasta.GZ", "wt") as fh:      # accepted extension, wrong case
        fh.write(">chr\nACGT\n")
    (folder / "g4.fasta.bz2").write_bytes(bz2.compress(b">chr\nACGT\n"))
    (folder / "notes.txt").write_text("nothing to see\n")    # not a genome: still no warning
    # Neither of these is compressed data, and both used to be named as genomes: a zip is how a
    # spreadsheet of strain metadata looks, and bzip2's magic is "BZh" followed by a block-size digit.
    (folder / "strains.xlsx").write_bytes(b"PK\x03\x04" + b"\x00" * 20)
    (folder / "codes.txt").write_text("BZh is a strain code\n")
    with caplog.at_level("WARNING"):
        found = require_genomes(folder, "exclusion")
    assert [path.name for path in found] == ["g1.fasta", "g2.fasta"]
    for name in ("g3.fasta.GZ", "g4.fasta.bz2"):
        assert name in caplog.text, name
    for name in ("notes.txt", "strains.xlsx", "codes.txt"):
        assert name not in caplog.text, name


def test_a_gzipped_file_named_with_capitals_is_read_as_gzip(tmp_path):
    """`.GZ` comes from Windows. Nothing proved the reader handles it: the file is named in the "not part
    of the group" warning either way, which is a different claim."""
    path = tmp_path / "g.FASTA.GZ"
    with gzip.open(path, "wt") as fh:
        fh.write(">chr\nACGTACGT\n")
    with open_text(path) as fh:
        assert fh.readline().startswith(">chr")
    assert is_fasta(path)


def test_compressed_as_names_the_format_from_the_first_bytes(tmp_path):
    path = tmp_path / "x"
    for data, expected in ((bz2.compress(b"hello"), "bzip2"),
                           (b"BZh is a strain code", None),   # the magic needs its block-size digit
                           (b"PK\x03\x04junk", None),          # a zip is not a genome this can read
                           (b"\xfd7zXZ\x00junk", "xz"),
                           (b"\x28\xb5/\xfdjunk", "zstandard"),
                           (b">chr\nACGT\n", None)):
        path.write_bytes(data)
        assert compressed_as(path) == expected, data[:4]
    with gzip.open(path, "wb") as fh:
        fh.write(b">chr\nACGT\n")
    assert compressed_as(path) == "gzip"
    assert compressed_as(tmp_path / "gone") is None


def test_a_named_pipe_with_an_accepted_extension_is_refused_not_waited_on(tmp_path):
    """`require_genomes` reads every path with an accepted extension to check it is fasta, and opening a
    named pipe with nothing on the other end never returns: the run used to hang at discovery with no
    output at all, not even its version banner. A subprocess with a timeout, since the failure is a hang."""
    folder = tmp_path / "inclusion"
    folder.mkdir()
    (folder / "g1.fasta").write_text(">chr\nACGT\n")
    os.mkfifo(folder / "pipe.fasta")
    done = subprocess.run(
        [sys.executable, "-c",
         "import sys; from primer_finder import PrimerFinderError;"
         "from primer_finder.seqio import require_genomes;"
         "from pathlib import Path;"
         "\ntry:\n print(len(require_genomes(Path(sys.argv[1]), 'inclusion')))"
         "\nexcept PrimerFinderError as exc:\n print('refused:', exc)",
         str(folder)],
        capture_output=True, text=True, timeout=30, cwd=Path(__file__).resolve().parent.parent)
    assert done.returncode == 0, done.stderr
    assert "not ordinary files" in done.stdout, done.stdout
    assert "pipe.fasta" in done.stdout


def test_a_named_pipe_in_a_genome_folder_does_not_stop_the_run(tmp_path):
    """Looking for fasta files with an unaccepted extension means opening whatever is in the folder, and
    opening a named pipe with nothing on the other end never returns."""
    folder = tmp_path / "inclusion"
    folder.mkdir()
    (folder / "g1.fasta").write_text(">chr\nACGT\n")
    os.mkfifo(folder / "pipe")
    done = subprocess.run(
        [sys.executable, "-c",
         "import sys; from primer_finder.seqio import require_genomes;"
         "print(len(require_genomes(__import__('pathlib').Path(sys.argv[1]), 'inclusion')))",
         str(folder)],
        capture_output=True, text=True, timeout=30, cwd=Path(__file__).resolve().parent.parent)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == "1"


def test_a_symlink_loop_in_a_genome_folder_is_not_a_crash(tmp_path):
    """`resolve()` on a loop raises -- RuntimeError on some versions, OSError on others -- and the command
    line turns neither into a message."""
    folder = tmp_path / "inclusion"
    folder.mkdir()
    (folder / "g1.fasta").write_text(">chr\nACGT\n")
    (folder / "loop_a").symlink_to(folder / "loop_b")
    (folder / "loop_b").symlink_to(folder / "loop_a")
    assert len(require_genomes(folder, "inclusion")) == 1


def test_a_fasta_with_an_unaccepted_extension_is_reported_not_ignored(tmp_path, caplog):
    """A genome left out of the exclusion group is the dangerous direction: its kmers then read as
    specific to the inclusion group. Being quietly dropped for its extension is how that happens."""
    folder = tmp_path / "exclusion"
    folder.mkdir()
    for name in ("g1.fasta", "g2.fasta"):
        (folder / name).write_text(">chr\nACGT\n")
    (folder / "g3.fas").write_text(">chr\nACGT\n")      # fasta, unaccepted extension
    (folder / "notes.txt").write_text("nothing to see\n")  # not fasta: no warning for this one
    with caplog.at_level("WARNING"):
        found = require_genomes(folder, "exclusion")
    assert [path.name for path in found] == ["g1.fasta", "g2.fasta"]
    assert "g3.fas" in caplog.text
    assert "notes.txt" not in caplog.text


def test_two_paths_to_one_genome_count_once(tmp_path, caplog):
    """`find_genomes` follows links, so a link to the folder itself used to multiply the whole group by
    the depth it descended."""
    folder = tmp_path / "inclusion"
    folder.mkdir()
    (folder / "g1.fasta").write_text(">chr\nACGT\n")
    (folder / "also_g1.fasta").symlink_to(folder / "g1.fasta")
    with caplog.at_level("WARNING"):
        found = require_genomes(folder, "inclusion")
    assert len(found) == 1
    assert "already counted" in caplog.text


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
