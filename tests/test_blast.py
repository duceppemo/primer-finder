"""Reading blast output, and what the exclusion genomes say about a contig."""

from __future__ import annotations

from pathlib import Path

import pytest

from primer_finder.blast import (
    ALIGNMENT_FIELDS,
    PRESENCE_FIELDS,
    ExclusionResult,
    Hit,
    exclusion_variants,
    has_close_variants,
    make_db,
    parse_hits,
    presence_in_genomes,
    shared_variants,
    variant_positions,
)


def test_parse_hits_ignores_malformed_lines(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("ctg1\t1\t20\t1e-30\tACGT\tACGA\nbroken\n")
    hits = parse_hits(path, ALIGNMENT_FIELDS)
    assert len(hits) == 1
    assert (hits[0].query, hits[0].qstart, hits[0].qend, hits[0].evalue) == ("ctg1", 1, 20, 1e-30)


def test_parse_hits_without_the_alignment_fields(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("ctg1\t1e-30\n")
    hits = parse_hits(path, PRESENCE_FIELDS)
    assert (hits[0].query, hits[0].qseq, hits[0].qstart) == ("ctg1", "", 0)


def test_variant_positions_are_contig_coordinates():
    # The hit starts at base 5 of the contig (1-based), so the first aligned base is position 4
    hit = Hit("ctg1", 5, 12, 1e-30, "ACGTACGT", "ACCTACGA")
    assert variant_positions(hit) == [6, 11]


def test_variant_positions_of_an_insertion_in_the_contig():
    hit = Hit("ctg1", 1, 8, 1e-30, "ACGTTTACGT", "ACGT--ACGT")
    assert variant_positions(hit) == [4, 5]


def test_variant_positions_of_a_deletion_in_the_contig():
    """The exclusion genome has bases the contig lacks: the position where they would fit is a variant."""
    hit = Hit("ctg1", 1, 8, 1e-30, "ACGT--ACGT", "ACGTTTACGT")
    assert variant_positions(hit) == [4]


def test_lower_case_sequences_are_not_variants():
    assert variant_positions(Hit("ctg1", 1, 4, 1e-30, "acgt", "ACGT")) == []


def test_has_close_variants():
    assert has_close_variants([5, 26]) is False  # 21 bases apart: no primer covers both
    assert has_close_variants([5, 25]) is True   # 20 bases apart
    assert has_close_variants([]) is False
    assert has_close_variants([1, 100, 110]) is True


def test_shared_variants_needs_most_of_the_genomes():
    result = ExclusionResult(genomes_hit=3, positions=[[1, 2], [1, 2], [1]])
    assert shared_variants(result) == [1]  # Position 2 differs in only two of the three genomes
    assert shared_variants(result, fraction=0.5) == [1, 2]
    assert shared_variants(ExclusionResult()) == []


def test_make_db_decompresses_a_gzipped_genome(stubs, tmp_path):
    import gzip

    genome = tmp_path / "g.fasta.gz"
    with gzip.open(genome, "wt") as fh:
        fh.write(">chr\nACGT\n")
    db = make_db(genome, tmp_path / "db")
    assert db.name == "g"
    assert (tmp_path / "db" / "g.fasta").read_text() == ">chr\nACGT\n"


def test_presence_in_genomes(stubs, tmp_path):
    stubs(presence={"g1": ["ctg1"], "g2": []})
    query = tmp_path / "q.fasta"
    query.write_text(">ctg1\nACGT\n")
    genomes = []
    for name in ("g1", "g2"):
        genome = tmp_path / f"{name}.fasta"
        genome.write_text(">chr\nACGT\n")
        genomes.append(genome)
    presence = presence_in_genomes(query, genomes, tmp_path / "work", threads=2)
    assert presence == {"ctg1": {"g1.fasta": True, "g2.fasta": False}}


def test_exclusion_variants_keeps_the_best_hit_of_each_genome(stubs, tmp_path):
    stubs(exclusion_hits={"g1": [
        ["ctg1", 1, 4, 1e-30, "ACGT", "ACGA"],   # The best hit: one variant at position 3
        ["ctg1", 1, 4, 1e-11, "ACGT", "AAAA"],   # A weaker hit, ignored
    ]})
    query = tmp_path / "q.fasta"
    query.write_text(">ctg1\nACGT\n")
    genome = tmp_path / "g1.fasta"
    genome.write_text(">chr\nACGT\n")
    results = exclusion_variants(query, [genome], tmp_path / "work", threads=1)
    assert results["ctg1"].genomes_hit == 1
    assert results["ctg1"].positions == [[3]]


def test_a_contig_no_exclusion_genome_hits(stubs, tmp_path):
    stubs(exclusion_hits={"g1": []})
    query = tmp_path / "q.fasta"
    query.write_text(">ctg1\nACGT\n")
    genome = tmp_path / "g1.fasta"
    genome.write_text(">chr\nACGT\n")
    results = exclusion_variants(query, [genome], tmp_path / "work", threads=1)
    assert results["ctg1"] == ExclusionResult(genomes_hit=0, positions=[])


def test_no_genome_at_all(stubs, tmp_path):
    query = tmp_path / "q.fasta"
    query.write_text(">ctg1\nACGT\n")
    assert presence_in_genomes(query, [], tmp_path / "work", threads=1) == {"ctg1": {}}


@pytest.mark.parametrize("fields", [ALIGNMENT_FIELDS, PRESENCE_FIELDS])
def test_blastn_asks_for_the_fields_it_parses(stubs, tmp_path, fields):
    from primer_finder.blast import blastn

    calls = stubs()
    query = tmp_path / "q.fasta"
    query.write_text(">ctg1\nACGT\n")
    out = blastn(Path("db"), query, tmp_path / "out.tsv", fields)
    assert f"-outfmt '6 {' '.join(fields)}'" in calls.read_text()
    assert out.exists()
