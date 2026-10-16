"""Reading blast output, and what the exclusion genomes say about a contig."""

from __future__ import annotations

from pathlib import Path

import pytest

from primer_finder.blast import (
    ALIGNMENT_FIELDS,
    DB_NAME,
    LOCAL_FASTA,
    PRESENCE_FIELDS,
    SEARCH_TASK,
    SHORTEST_ABSENT_STRETCH,
    Alignments,
    ExclusionResult,
    Hit,
    blastn,
    enough_variants,
    exclusion_variants,
    genome_folder,
    has_close_variants,
    make_db,
    parse_hits,
    presence_in_genomes,
    shared_variants,
    unaligned_positions,
    variant_bases,
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


def test_a_masked_contig_base_is_not_a_difference():
    """blast lower-cases what it masks on either side; the bases still match."""
    assert variant_bases(hit(1, 4, "acgt", "ACGT")) == {}


def test_aligned_sequences_of_unequal_length_do_not_crash():
    assert variant_bases(hit(1, 4, "ACGT", "AC")) == {}


def test_has_close_variants():
    """The window is the widest primer, so two positions share one when the second is fewer than that many
    bases past the first -- a span of 25 for a 25-base primer."""
    assert has_close_variants([5, 30]) is False  # span 26: wider than any primer
    assert has_close_variants([5, 29]) is True   # span 25: the widest primer covers both
    assert has_close_variants([5, 25]) is True
    assert has_close_variants([]) is False
    assert has_close_variants([1, 100, 110]) is True


def test_enough_variants_at_each_stringency():
    """The strict rule wants two close enough to share a primer; the last-resort rule wants one, anywhere.
    Neither accepts a region with no shared variant at all, which would be no assay."""
    assert enough_variants([5, 25]) and enough_variants([5, 25], 1)
    assert not enough_variants([5, 30]) and enough_variants([5, 30], 1)  # span 26: too wide for a primer
    assert not enough_variants([5]) and enough_variants([5], 1)
    assert not enough_variants([]) and not enough_variants([], 1)


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


def test_the_search_is_sensitive_enough_to_find_a_copy_that_differs(stubs, tmp_path):
    """blastn's default task is megablast, which seeds on an exact 28-base match and therefore finds
    *nothing* in a copy whose differences are closer together than that. Every search here exists to find
    homology an assay would amplify, so missing it is the dangerous direction: an amplicon an exclusion
    genome holds at 97% identity was reported absent, and the assay on it labelled specific by absence.
    """
    calls = stubs()
    path = genome(tmp_path / "in", "g.fasta")
    db = make_db(path, tmp_path / "db")
    blastn(db, query_file(tmp_path), tmp_path / "db" / "hits.tsv", PRESENCE_FIELDS)
    asked = [line for line in calls.read_text().splitlines() if line.startswith("blastn ")]
    assert asked and all(f"-task {SEARCH_TASK}" in line for line in asked), asked
    assert SEARCH_TASK == "blastn", "megablast cannot see a copy with differences under 28 bases apart"


def test_a_stretch_no_exclusion_genome_aligns_counts_as_differing():
    """The shape a repeat-targeted marker has. A multi-copy element is a novel block in shared backbone,
    and the junctions either side of it are inclusion-specific too, so the assembly is flank + block +
    flank. Every exclusion genome then matches the flanks -- enough for `genomes_hit` to be non-zero --
    and nothing in between. Measured on a 444-base contig carrying a 300-base marker: hits at `q1-72` and
    `q373-444` in every exclusion genome, nothing between, and the region used to be dropped for having
    no differences."""
    flanks = ExclusionResult(genomes=[])
    for _ in range(4):
        alignments = Alignments()
        alignments.spans.append((0, 72, {}))      # q1-72, matching exactly
        alignments.spans.append((372, 444, {}))   # q373-444, the same
        flanks.genomes.append(alignments)
    assert shared_variants(flanks) == []          # Nothing differs where they do align
    missing = unaligned_positions(flanks, 444)
    assert missing == list(range(72, 372))        # The 300 bases of the marker
    assert enough_variants(missing, 2)


def one_genome_aligning(length: int, *spans: tuple[int, int]) -> ExclusionResult:
    """One exclusion genome that aligns exactly the given half-open spans of a contig."""
    alignments = Alignments()
    for start, end in spans:
        alignments.spans.append((start, end, {}))
    return ExclusionResult(genomes=[alignments])


@pytest.mark.parametrize(("gap", "counted"),
                         [(1, False), (6, False), (32, False),
                          (SHORTEST_ABSENT_STRETCH - 1, False), (SHORTEST_ABSENT_STRETCH, True)])
def test_how_long_an_unaligned_run_has_to_be_before_it_is_an_absence(gap, counted):
    """The gap between one contig's local alignments is not necessarily small -- measured per genome
    across five *Listeria* complexes it runs to 137 bases and beyond. What has to decide this is whether
    the search would have reported the stretch had the genome held it: against a 3 Mb genome at this
    module's `-evalue 1e-10`, an exact island of 25 or 30 bases is not reported at all and one of 33 is,
    so anything shorter than the floor could be present and merely invisible."""
    result = one_genome_aligning(400, (0, 50), (50 + gap, 400))
    assert (unaligned_positions(result, 400) == list(range(50, 50 + gap))) is counted


def test_the_floor_is_above_what_the_search_can_miss():
    """Not a round number chosen for comfort: 33 bases is where an exact island becomes reportable at
    `-evalue 1e-10` against a 3 Mb genome, and the floor leaves room for a larger one."""
    assert SHORTEST_ABSENT_STRETCH >= 33


def test_a_contig_no_exclusion_genome_holds_at_all_has_no_unaligned_positions():
    """`genomes_hit == 0` already keeps such a contig, so this reports nothing rather than every base:
    the two paths must not both claim it."""
    assert unaligned_positions(ExclusionResult(genomes=[]), 200) == []


def test_a_position_one_exclusion_genome_aligns_is_not_unaligned():
    """Absent from the group means absent from all of it: one genome holding that stretch is enough to
    make it shared."""
    result = ExclusionResult(genomes=[])
    first, second = Alignments(), Alignments()
    first.spans.append((0, 50, {}))     # This one holds the first half
    second.spans.append((50, 100, {}))  # and this one the second
    result.genomes.extend([first, second])
    assert unaligned_positions(result, 100) == []
