"""Reading blast output, and what the exclusion genomes say about a contig."""

from __future__ import annotations

from pathlib import Path

import pytest

from primer_finder.blast import (
    ALIGNMENT_FIELDS,
    DB_NAME,
    LOCAL_FASTA,
    PRESENCE_FIELDS,
    Alignments,
    ExclusionResult,
    Hit,
    blastn,
    exclusion_variants,
    genome_folder,
    has_close_variants,
    make_db,
    parse_hits,
    presence_in_genomes,
    shared_variants,
    variant_bases,
    variant_positions,
)


def genome(folder: Path, name: str, sequence: str = "ACGT" * 10) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    path.write_text(f">chr\n{sequence}\n")
    return path


def query_file(folder: Path, names: tuple[str, ...] = ("ctg1",)) -> Path:
    path = folder / "q.fasta"
    path.write_text("".join(f">{name}\nACGT\n" for name in names))
    return path


def hit(start: int, end: int, qseq: str, sseq: str, evalue: float = 1e-30) -> Hit:
    return Hit("ctg1", start, end, evalue, qseq, sseq)


def alignments(*hits: Hit) -> Alignments:
    result = Alignments()
    for one in hits:
        result.add(one)
    return result


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


def test_variant_bases_give_what_the_exclusion_genome_has_there():
    """How much a mismatch costs depends on which bases are involved, so the base is carried, not only the
    position. The hit starts at base 5 of the contig (1-based), so the first aligned base is position 4."""
    assert variant_bases(hit(5, 12, "ACGTACGT", "ACCTACGA")) == {6: "C", 11: "A"}


def test_variant_bases_of_an_insertion_in_the_contig():
    """The contig has bases the exclusion genome does not: there is nothing to pair with, which is a
    difference of its own and is marked as such."""
    assert variant_bases(hit(1, 8, "ACGTTTACGT", "ACGT--ACGT")) == {4: "-", 5: "-"}


def test_variant_bases_of_a_deletion_in_the_contig_take_the_first_of_them():
    """The exclusion genome has two bases the contig lacks, and they both fall at the same contig
    position. The first is the one an oligo would run into, so it is the one kept."""
    assert variant_bases(hit(1, 8, "ACGT--ACGT", "ACGTTTACGT")) == {4: "T"}
    assert variant_bases(hit(1, 8, "ACGT--ACGT", "ACGTGCACGT")) == {4: "G"}  # not the C


def test_variant_bases_are_upper_case_whatever_the_alignment_says():
    """blast lower-cases what it masks, and a difference is a difference whatever the case."""
    assert variant_bases(hit(1, 4, "ACGT", "acga")) == {3: "A"}


def test_variant_positions_are_contig_coordinates():
    # The hit starts at base 5 of the contig (1-based), so the first aligned base is position 4
    assert variant_positions(hit(5, 12, "ACGTACGT", "ACCTACGA")) == [6, 11]


def test_variant_positions_of_an_insertion_in_the_contig():
    assert variant_positions(hit(1, 8, "ACGTTTACGT", "ACGT--ACGT")) == [4, 5]


def test_variant_positions_of_a_deletion_in_the_contig():
    """The exclusion genome has bases the contig lacks: the position where they would fit is a variant."""
    assert variant_positions(hit(1, 8, "ACGT--ACGT", "ACGTTTACGT")) == [4]


def test_lower_case_sequences_are_not_variants():
    assert variant_positions(hit(1, 4, "acgt", "ACGT")) == []


def test_aligned_sequences_of_unequal_length_do_not_crash():
    assert variant_positions(hit(1, 4, "ACGT", "AC")) == []


def test_has_close_variants():
    assert has_close_variants([5, 26]) is False  # 21 bases apart: no primer covers both
    assert has_close_variants([5, 25]) is True   # 20 bases apart
    assert has_close_variants([]) is False
    assert has_close_variants([1, 100, 110]) is True


def test_shared_variants_needs_most_of_the_genomes_that_hold_the_region():
    """A position must differ in at least 90% of the exclusion genomes that align it."""
    differs = alignments(hit(1, 10, "ACGTACGTAC", "ACGTAAGTAC"))  # variant at position 5
    identical = alignments(hit(1, 10, "ACGTACGTAC", "ACGTACGTAC"))
    assert shared_variants(ExclusionResult([differs] * 9 + [identical])) == [5]   # 9 of 10
    assert shared_variants(ExclusionResult([differs] * 8 + [identical] * 2)) == []  # 8 of 10
    assert shared_variants(ExclusionResult([differs, identical])) == []  # 1 of 2
    assert shared_variants(ExclusionResult()) == []


def test_a_genome_that_does_not_hold_the_region_does_not_count():
    """Its hit covers the first half only; it says nothing about a position in the second half, which is
    not the same as saying that the position is identical."""
    partial = alignments(hit(1, 200, "A" * 200, "A" * 200))
    whole = alignments(hit(1, 400, "A" * 400, "A" * 299 + "C" + "A" * 10 + "C" + "A" * 89))
    assert shared_variants(ExclusionResult([partial, whole])) == [299, 310]


def test_a_second_perfect_copy_in_the_same_genome_rules_the_position_out():
    """An assay would amplify that copy, so a position matching any copy is not specific."""
    diverged = hit(1, 600, "A" * 600, "A" * 299 + "C" + "A" * 300)
    perfect = hit(100, 450, "A" * 351, "A" * 351)
    assert shared_variants(ExclusionResult([alignments(diverged)])) == [299]
    assert shared_variants(ExclusionResult([alignments(diverged, perfect)])) == []


def test_alignments_cover_and_differ():
    one = alignments(hit(11, 20, "ACGTACGTAC", "ACGTAAGTAC"))
    assert one.covers(10) and one.covers(19) and not one.covers(20) and not one.covers(9)
    assert one.differs_at(15)  # The sixth base of the alignment, which starts at position 10
    assert not one.differs_at(10)
    assert not one.differs_at(100)  # Outside every alignment
    assert one.variants == {15}


def test_make_db_links_the_genome_and_builds_the_database_in_its_folder(stubs, tmp_path):
    path = genome(tmp_path / "in", "g.fasta")
    db = make_db(path, tmp_path / "db")
    assert db == tmp_path / "db" / DB_NAME
    assert (tmp_path / "db" / LOCAL_FASTA).resolve() == path.resolve()


def test_make_db_decompresses_a_gzipped_genome(stubs, tmp_path):
    import gzip

    path = tmp_path / "g.fasta.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(">chr\nACGT\n")
    make_db(path, tmp_path / "db")
    assert (tmp_path / "db" / LOCAL_FASTA).read_text() == ">chr\nACGT\n"


def test_make_db_can_be_run_again_over_an_existing_folder(stubs, tmp_path):
    path = genome(tmp_path / "in", "g.fasta")
    make_db(path, tmp_path / "db")
    make_db(path, tmp_path / "db")  # The link is already there
    assert (tmp_path / "db" / LOCAL_FASTA).is_symlink()


def test_blast_is_given_no_path_that_could_hold_a_space(stubs, tmp_path):
    """BLAST splits its path arguments on whitespace, so it runs in the database's folder."""
    calls = stubs()
    folder = tmp_path / "with space" / "db"
    path = genome(tmp_path / "with space" / "in", "g.fasta")
    db = make_db(path, folder)
    blastn(db, query_file(tmp_path / "with space"), folder / "hits.tsv", PRESENCE_FIELDS)
    for line in calls.read_text().splitlines():
        arguments = line.split()[1:]
        assert "space" not in line, line
        assert all(not part.startswith("/") for part in arguments), line


def test_genome_folder_is_unique_per_genome(tmp_path):
    first = genome(tmp_path / "a", "contigs.fasta")
    second = genome(tmp_path / "b", "contigs.fasta")
    assert genome_folder(tmp_path, 0, first) != genome_folder(tmp_path, 1, second)
    assert genome_folder(tmp_path, 3, first).name == "0003_contigs"


def test_presence_in_genomes(stubs, tmp_path):
    stubs(presence={"g1": ["ctg1"], "g2": []})
    genomes = [genome(tmp_path / "in", f"{name}.fasta") for name in ("g1", "g2")]
    presence = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=2)
    assert presence == {"ctg1": {str(genomes[0]): True, str(genomes[1]): False}}


def test_two_genomes_with_the_same_name_stay_apart(stubs, tmp_path):
    stubs(presence={"a/g": ["ctg1"], "b/g": []})
    genomes = [genome(tmp_path / "in" / sample, "g.fasta") for sample in ("a", "b")]
    presence = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=2)
    assert presence == {"ctg1": {str(genomes[0]): True, str(genomes[1]): False}}


def test_a_hit_that_is_not_good_enough_is_not_a_presence(stubs, tmp_path):
    """blastn is told -evalue 1e-10, and the value is checked again when the output is read."""
    stubs(presence_evalue="1e-3")
    genomes = [genome(tmp_path / "in", "g1.fasta")]
    presence = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=1)
    assert presence == {"ctg1": {str(genomes[0]): False}}


def test_exclusion_variants_reads_every_hit_of_each_genome(stubs, tmp_path):
    stubs(exclusion_hits={"g1": [
        ["ctg1", 1, 4, 1e-30, "ACGT", "ACGA"],
        ["ctg1", 1, 4, 1e-11, "ACGT", "AAAA"],
        ["ctg1", 1, 4, 1e-3, "ACGT", "TTTT"],  # Not good enough: ignored
    ]})
    genomes = [genome(tmp_path / "in", "g1.fasta")]
    results = exclusion_variants(query_file(tmp_path), genomes, tmp_path / "work", threads=1)
    assert results["ctg1"].genomes_hit == 1
    assert [sorted(variants) for _, _, variants in results["ctg1"].genomes[0].spans] == [[3], [1, 2, 3]]
    assert shared_variants(results["ctg1"]) == [3]  # Position 3 differs in both copies, 1 and 2 do not


def test_a_contig_no_exclusion_genome_hits(stubs, tmp_path):
    stubs(exclusion_hits={"g1": []})
    genomes = [genome(tmp_path / "in", "g1.fasta")]
    results = exclusion_variants(query_file(tmp_path), genomes, tmp_path / "work", threads=1)
    assert results["ctg1"].genomes_hit == 0
    assert shared_variants(results["ctg1"]) == []


def test_no_genome_at_all(stubs, tmp_path):
    assert presence_in_genomes(query_file(tmp_path), [], tmp_path / "work", threads=1) == {"ctg1": {}}


@pytest.mark.parametrize("fields", [ALIGNMENT_FIELDS, PRESENCE_FIELDS])
def test_blastn_asks_for_the_fields_it_parses(stubs, tmp_path, fields):
    calls = stubs()
    path = genome(tmp_path / "in", "g.fasta")
    db = make_db(path, tmp_path / "db")
    out = blastn(db, query_file(tmp_path), tmp_path / "db" / "hits.tsv", fields)
    assert f"-outfmt '6 {' '.join(fields)}'" in calls.read_text()
    assert out.exists()
