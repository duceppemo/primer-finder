"""Counting kmers with KMC."""

from __future__ import annotations

import pytest

from primer_finder import PrimerFinderError
from primer_finder.kmers import (
    NO_MAX_COUNT,
    count,
    dump,
    dump_to_fasta,
    kmc_memory,
    require_db,
    subtract,
    write_file_list,
)


def test_kmc_memory_stays_in_the_range_kmc_accepts():
    assert kmc_memory(0) == 1
    assert kmc_memory(64) == 64
    assert kmc_memory(5000) == 1024


def test_write_file_list(tmp_path):
    from pathlib import Path

    path = write_file_list([Path("/a/b.fasta"), Path("/c/d.fasta")], tmp_path / "list.txt")
    assert path.read_text() == "/a/b.fasta\n/c/d.fasta\n"


def test_dump_to_fasta_numbers_the_kmers(tmp_path):
    dump_file = tmp_path / "dump.txt"
    dump_file.write_text("ACGT\t4\n\nTTGA\t4\n")
    fasta = tmp_path / "kmers.fasta"
    assert dump_to_fasta(dump_file, fasta) == 2
    assert fasta.read_text() == ">kmer_0\nACGT\n>kmer_1\nTTGA\n"


def test_require_db_names_the_missing_file(tmp_path):
    (tmp_path / "db.kmc_pre").touch()
    with pytest.raises(PrimerFinderError, match="broke .db.kmc_suf missing."):
        require_db(tmp_path / "db", "broke")
    (tmp_path / "db.kmc_suf").touch()
    assert require_db(tmp_path / "db", "broke") == tmp_path / "db"


def test_count_and_subtract_and_dump(stubs, tmp_path):
    calls = stubs(kmers=["ACGT", "TTGA"])
    list_file = write_file_list([tmp_path / "g.fasta"], tmp_path / "list.txt")
    first = count(list_file, tmp_path / "a", tmp_path / "work", 31, 2, 4, min_count=2, max_count=NO_MAX_COUNT)
    second = count(list_file, tmp_path / "b", tmp_path / "work", 31, 2, 4, min_count=1, max_count=1)
    specific = subtract(first, second, tmp_path / "c", threads=2)
    dump_file = dump(specific, tmp_path / "dump.txt", threads=2)
    assert dump_file.read_text().splitlines() == ["ACGT\t2", "TTGA\t2"]
    assert "kmc -k31 -t2 -m4 -fm -ci2 -cx1000000000" in calls.read_text()
    assert (tmp_path / "work").is_dir()


def test_an_empty_dump_is_an_error(stubs, tmp_path):
    stubs(kmers=[])
    (tmp_path / "db.kmc_pre").touch()
    (tmp_path / "db.kmc_suf").touch()
    with pytest.raises(PrimerFinderError, match="No inclusion-specific kmer"):
        dump(tmp_path / "db", tmp_path / "dump.txt", threads=1)


def test_the_counter_ceiling_is_raised_for_a_large_inclusion_group(stubs, tmp_path):
    """KMC counts to 255 unless -cs says otherwise, so -ci300 would match nothing in 300 genomes."""
    calls = stubs()
    list_file = write_file_list([tmp_path / "g.fasta"], tmp_path / "list.txt")
    count(list_file, tmp_path / "small", tmp_path / "work", 31, 1, 2, min_count=4, max_count=4)
    count(list_file, tmp_path / "big", tmp_path / "work", 31, 1, 2, min_count=300, max_count=600)
    text = calls.read_text()
    assert "-fm -ci4 -cx4" in text  # Up to 255, KMC's default counter is enough
    assert "-fm -cs600 -ci300 -cx600" in text
