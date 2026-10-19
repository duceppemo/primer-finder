"""Reading blast output, and what the exclusion genomes say about a contig."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from primer_finder.blast import (
    ALIGNMENT_FIELDS,
    DB_NAME,
    LOCAL_FASTA,
    MOST_SUBJECTS,
    PRESENCE_FIELDS,
    SEARCH_TASK,
    SHORTEST_ABSENT_STRETCH,
    Alignments,
    ExclusionResult,
    Hit,
    blastn,
    covered_fraction,
    enough_variants,
    exclusion_variants,
    genome_folder,
    has_close_variants,
    in_one_piece,
    make_db,
    masked_intervals,
    overlaps_masked,
    parse_hits,
    presence_in_genomes,
    safe_stem,
    shared_variants,
    unaligned_positions,
    unaligned_stretches,
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
    """The presence search asks for the two spans and nothing more: no aligned sequences to read."""
    path = tmp_path / "hits.tsv"
    path.write_text("ctg1\t1\t20\t1e-30\tchr\t900\t919\n")
    hits = parse_hits(path, PRESENCE_FIELDS)
    assert (hits[0].query, hits[0].qseq, hits[0].qstart, hits[0].qend) == ("ctg1", "", 1, 20)
    assert (hits[0].subject, hits[0].sstart, hits[0].send) == ("chr", 900, 919)


def test_a_genome_name_blast_cannot_use_does_not_reach_the_database_folder():
    """makeblastdb reopens its own output by absolute path and splits it on whitespace, so a space in the
    folder it runs in makes it fail -- and the folder used to carry the genome's file stem. Naming an
    assembly after a strain is ordinary, and it cost the whole run at step 4."""
    assert safe_stem(Path("/x/Brucella abortus 2308.fasta")) == "Brucella_abortus_2308"
    assert safe_stem(Path('/x/B"abortus.fasta')) == "B_abortus"
    assert safe_stem(Path("/x/B'abortus.fasta")) == "B'abortus"       # a quote blast does accept
    assert safe_stem(Path("/x/plain.fasta")) == "plain"
    for name in ("a b.fasta", 'a"b.fasta', "a\tb.fasta"):
        folder = genome_folder(Path("/work"), 3, Path("/x") / name)
        assert folder.name == "0003_a_b", folder
        assert not any(character in folder.name for character in '" \t\n'), folder


@pytest.mark.tools
@pytest.mark.skipif(not shutil.which("makeblastdb"), reason="needs makeblastdb on PATH")
def test_a_database_is_built_for_a_genome_whose_name_holds_a_space(tmp_path):
    """The real program, which is the only thing that can show this: the stub cannot fail on its cwd."""
    genome = tmp_path / "Brucella abortus 2308.fasta"
    genome.write_text(">chr\n" + "ACGTACCGGTTA" * 40 + "\n")
    database = make_db(genome, genome_folder(tmp_path / "work", 0, genome))
    assert database.parent.name == "0000_Brucella_abortus_2308"
    assert (database.parent / f"{DB_NAME}.nin").is_file() or (database.parent / f"{DB_NAME}.nal").is_file()


def test_covered_fraction_merges_the_pieces_instead_of_adding_them():
    """A repeated region is hit once per copy, over the same columns of the query. Added up, a genome
    would look as though it held the contig several times over -- and the figure would be worst exactly
    where repeats are, which is what it is there to describe."""
    three_copies = [hit(1, 50, "", ""), hit(1, 50, "", ""), hit(1, 50, "", "")]
    assert covered_fraction(three_copies, 100) == 0.5
    assert covered_fraction([hit(1, 50, "", ""), hit(41, 60, "", "")], 100) == 0.6   # overlapping
    assert covered_fraction([hit(1, 20, "", ""), hit(81, 100, "", "")], 100) == 0.4  # two pieces
    assert covered_fraction([hit(1, 100, "", "")], 100) == 1.0
    assert covered_fraction([], 100) == 0.0


def test_covered_fraction_counts_pieces_a_draft_genome_split_between_its_contigs():
    """The case the figure exists for: a genome that broke the region in two holds all of it, and the two
    halves come back as two hits because they are on two of its sequences."""
    halves = [hit(1, 60, "", ""), hit(61, 120, "", "")]
    assert covered_fraction(halves, 120) == 1.0


def test_a_hit_reaching_past_the_contig_is_not_more_than_all_of_it():
    """Nothing blast reports should, but a fraction above 1 would read as a sure thing."""
    assert covered_fraction([hit(1, 150, "", "")], 100) == 1.0


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
    """The fraction of the contig each genome holds: all of it, or none."""
    stubs(presence={"g1": ["ctg1"], "g2": []})
    genomes = [genome(tmp_path / "in", f"{name}.fasta") for name in ("g1", "g2")]
    presence, _ = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=2)
    assert presence == {"ctg1": {str(genomes[0]): 1.0, str(genomes[1]): 0.0}}


def test_presence_in_genomes_reports_the_part_a_genome_holds(stubs, tmp_path):
    """A genome holding some of a contig is told from one holding all of it -- not a yes or no. Asserted
    on the number, not on its truthiness, since 1.0 and 0.0 satisfy an assertion written for True and
    False and so would a return to booleans."""
    contig = tmp_path / "long.fasta"
    contig.write_text(">ctg1\n" + "ACGT" * 50 + "\n")          # 200 bases
    stubs(presence={"g1": ["ctg1"], "g2": ["ctg1"]},
          presence_spans={"g2": {"ctg1": [[1, 24], [101, 130]]}})  # 24 + 30 of the 200
    genomes = [genome(tmp_path / "in", f"{name}.fasta") for name in ("g1", "g2")]
    presence, _ = presence_in_genomes(contig, genomes, tmp_path / "work", threads=2)
    assert presence["ctg1"][str(genomes[0])] == 1.0
    assert presence["ctg1"][str(genomes[1])] == pytest.approx(54 / 200)


def test_two_genomes_with_the_same_name_stay_apart(stubs, tmp_path):
    stubs(presence={"a/g": ["ctg1"], "b/g": []})
    genomes = [genome(tmp_path / "in" / sample, "g.fasta") for sample in ("a", "b")]
    presence, _ = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=2)
    assert presence == {"ctg1": {str(genomes[0]): 1.0, str(genomes[1]): 0.0}}


def test_a_hit_that_is_not_good_enough_is_not_a_presence(stubs, tmp_path):
    """blastn is told -evalue 1e-10, and the value is checked again when the output is read."""
    stubs(presence_evalue="1e-3")
    genomes = [genome(tmp_path / "in", "g1.fasta")]
    presence, _ = presence_in_genomes(query_file(tmp_path), genomes, tmp_path / "work", threads=1)
    assert presence == {"ctg1": {str(genomes[0]): 0.0}}


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
    assert presence_in_genomes(query_file(tmp_path), [], tmp_path / "work", threads=1) == (
        {"ctg1": {}}, {"ctg1": {}})


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


def test_the_exclusion_search_asks_for_more_than_ten_subjects(stubs, tmp_path, fasta):
    """A draft assembly can put a short stretch of a contig on a sequence of its own, and blast's default
    of ten subjects can crowd that sequence out -- which reads as the genome lacking a stretch it has.
    The stand-in ignores `-max_target_seqs`, so what the search was asked for is what can be checked."""
    calls = stubs()
    query = tmp_path / "contigs.fasta"
    query.write_text(">ctg1\n" + "ACGT" * 25 + "\n")
    exclusion_variants(query, [fasta("g1.fasta", {"chr": "ACGT" * 50})], tmp_path / "work", threads=1)
    asked = [line for line in calls.read_text().splitlines() if line.startswith("blastn")]
    assert asked, "no search was run"
    assert all(f"-max_target_seqs {MOST_SUBJECTS}" in line for line in asked), asked


def test_the_unaligned_stretches_of_a_contig():
    """What the second, unmasked pass is given to look at: what the first pass left untouched, in runs
    long enough for anything downstream to read as an absence."""
    sizes = {"ctg1": 200, "ctg2": 100}
    # One hit over the middle of ctg1, nothing at all on ctg2
    stretches = unaligned_stretches([hit(51, 150, "", "")], sizes)
    assert stretches == [("ctg1", 0, 50), ("ctg1", 150, 200), ("ctg2", 0, 100)]
    # A gap shorter than the floor is not one: nothing downstream would read it as an absence
    short = unaligned_stretches([hit(1, 90, "", ""), hit(111, 200, "", "")], {"ctg1": 200})
    assert short == []
    assert SHORTEST_ABSENT_STRETCH == 40
    assert unaligned_stretches([hit(1, 200, "", "")], {"ctg1": 200}) == []


def test_overlaps_masked():
    assert overlaps_masked([(10, 20)], 0, 11)
    assert overlaps_masked([(10, 20)], 19, 100)
    assert not overlaps_masked([(10, 20)], 20, 100)
    assert not overlaps_masked([(10, 20)], 0, 10)
    assert not overlaps_masked([], 0, 100)


# A 200-base element with a three-base period: low-complexity, and what blast masks out of a query.
LOW_COMPLEXITY = ("AAT" * 70)[:200]


@pytest.mark.tools
@pytest.mark.skipif(not all(shutil.which(program) for program in ("blastn", "makeblastdb", "dustmasker")),
                    reason="needs blastn, makeblastdb and dustmasker on PATH")
def test_a_low_complexity_contig_a_genome_holds_exactly_is_not_reported_absent(tmp_path):
    """The real programs, which are the only thing that can show this. blast masks a low-complexity
    query by default, and a masked stretch is not found however exactly the genome holds it -- so the
    contig was reported absent from every exclusion genome, which is the strongest claim this tool
    makes, and missing from every inclusion genome, which drops the region altogether.
    """
    import random

    rng = random.Random(29)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    genome = tmp_path / "genome.fasta"
    genome.write_text(f">chr\n{dna(3000)}{LOW_COMPLEXITY}{dna(3000)}\n")
    query = tmp_path / "contigs.fasta"
    query.write_text(f">ctg1\n{LOW_COMPLEXITY}\n")

    # What the masker says, and what the default search can therefore see: nothing
    work = tmp_path / "work"
    assert masked_intervals(query, work)["ctg1"] == [(0, 200)]
    folder = genome_folder(work, 0, genome)
    db = make_db(genome, folder)
    assert parse_hits(blastn(db, query, folder / "masked.tsv", PRESENCE_FIELDS), PRESENCE_FIELDS) == []
    assert parse_hits(blastn(db, query, folder / "unmasked.tsv", PRESENCE_FIELDS, dust=False),
                      PRESENCE_FIELDS), "the genome holds it exactly, masking or not"

    # And what the two search steps now say about it
    held, whole = presence_in_genomes(query, [genome], tmp_path / "presence", threads=1)
    assert held["ctg1"][str(genome)] == 1.0
    assert whole["ctg1"][str(genome)] == 1.0
    found = exclusion_variants(query, [genome], tmp_path / "exclusion", threads=1)
    assert found["ctg1"].genomes_hit == 1
    assert unaligned_positions(found["ctg1"], len(LOW_COMPLEXITY)) == []


def piece(qstart: int, qend: int, sstart: int, send: int, subject: str = "chr") -> Hit:
    """One hit of a contig in a genome, with where it sits in the genome."""
    return Hit("ctg1", qstart, qend, 1e-30, subject=subject, sstart=sstart, send=send)


def test_how_much_of_a_contig_a_genome_holds_in_one_piece():
    """`covered_fraction` beside this asks nothing about where in the genome the pieces came from, so two
    halves of a contig landing 100 bases apart -- which is what a region spliced across a tandem repeat
    does -- read as the whole contig being there."""
    halves = [piece(1, 100, 1000, 1099), piece(51, 150, 1150, 1249)]
    assert covered_fraction(halves, 150) == 1.0        # 150 of the 150 bases, in two pieces
    assert in_one_piece(halves, 150) == 100 / 150      # and no more than 100 of them in one
    # A genome that holds the contig whole, in one alignment or in two that follow each other
    assert in_one_piece([piece(1, 150, 1000, 1149)], 150) == 1.0
    assert in_one_piece([piece(1, 60, 1000, 1059), piece(61, 150, 1060, 1149)], 150) == 1.0
    # Two copies of the contig side by side in the genome are two pieces, not one longer one
    assert in_one_piece([piece(1, 150, 1000, 1149), piece(1, 150, 1200, 1349)], 150) == 1.0
    # On two sequences of a draft assembly, which is the case the merged fraction exists for
    assert covered_fraction([piece(1, 60, 1, 60), piece(61, 150, 1, 90, "ctg_b")], 150) == 1.0
    assert in_one_piece([piece(1, 60, 1, 60), piece(61, 150, 1, 90, "ctg_b")], 150) == 90 / 150
    assert in_one_piece([], 150) == 0.0
