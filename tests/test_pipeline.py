"""The pipeline from end to end, with the external programs replaced by stubs."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.pipeline import (
    HITS_NAME,
    Settings,
    lower_positions,
    required_genomes,
    run,
    write_presence_table,
)
from primer_finder.seqio import iter_records


def final_records(output: Path) -> dict[str, str]:
    from primer_finder.seqio import read_fasta

    return {name: record.seq for name, record in read_fasta(output / "final_kmers.fasta").items()}


def test_full_run(stubs, settings):
    scenario = settings()
    assert run(scenario) == 0
    assert final_records(scenario.output) == {"ctg1": "ACGTAcGTACGTaCGTACGT" + "ACGT" * 20}
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["counts"] == {"kmers": 3, "contigs": 1, "candidates": 1, "in_all_inclusion": 1, "final": 1}
    assert Path(info["reference"]).name == "exclusion_1.fasta"
    assert info["parameters"]["assembler"] == "skesa"
    # Every program used must have reported a version, including the two KMC ones (no --version flag)
    assert set(info["programs"]) == {"blastn", "kmc", "kmc_tools", "makeblastdb", "minimap2", "skesa"}
    assert all(version == f"{name} stub 1.0" for name, version in info["programs"].items()), info["programs"]
    # The variant positions are in the header, in contig coordinates
    assert "[5, 12]" in (scenario.output / "final_kmers.fasta").read_text()


def test_intermediate_files_are_removed(stubs, settings):
    scenario = settings()
    run(scenario)
    assert not (scenario.output / "1_kmers" / "kmers.txt").exists()
    assert not (scenario.output / "1_kmers" / "inclusion.kmc_pre").exists()
    assert not (scenario.output / "3_candidates" / "mapping.sam").exists()
    assert not (scenario.output / "4_blast" / "inclusion_db").exists()
    # What the user looks at is kept
    assert (scenario.output / "1_kmers" / "inclusion_specific_99mers.fasta").exists()
    assert (scenario.output / "2_assembly" / "assembly.fasta").exists()
    assert (scenario.output / "3_candidates" / "best_kmers.fasta").exists()
    assert (scenario.output / "4_blast" / "inclusion_blast_hits.tsv").exists()
    assert (scenario.output / "primer_finder.log").read_text().count("[INFO]") > 5


def test_keep_intermediate(stubs, settings):
    scenario = settings(keep_intermediate=True)
    run(scenario)
    assert (scenario.output / "1_kmers" / "kmers.txt").exists()
    assert (scenario.output / "3_candidates" / "mapping.sam").exists()
    assert (scenario.output / "4_blast" / "exclusion_db").is_dir()
    assert (scenario.output / "4_blast" / "exclusion_db" / "0000_exclusion_1" / "hits.tsv").exists()


def test_programs_are_called_with_the_right_options(stubs, settings):
    calls = stubs()
    scenario = settings(kmer_size=31, duplication=2)
    run(scenario)
    text = calls.read_text()
    assert "kmc -k31 -t1 -m2 -fm -ci2 -cx4" in text  # Two inclusion genomes, -d 2
    assert "-ci1 -cx1000000000" in text  # The exclusion group: every kmer, however often it occurs
    assert "kmc_tools -t1 simple" in text and "kmers_subtract" in text
    assert "minimap2 -t 1 -a --eqx" in text
    assert "skesa --cores 1 --mem 2" in text
    assert "spades.py" not in text


def test_spades_is_used_when_asked(stubs, settings):
    calls = stubs()
    run(settings(assembler="spades"))
    assert "spades.py --s 1" in calls.read_text()
    assert "--only-assembler" in calls.read_text()


def test_gzipped_genomes_are_accepted(stubs, settings, tmp_path):
    with gzip.open(tmp_path / "inclusion" / "inclusion_3.fasta.gz", "wt") as fh:
        fh.write(">chr\nACGT\n")
    scenario = settings()
    assert run(scenario) == 0
    header = (scenario.output / "4_blast" / "inclusion_blast_hits.tsv").read_text().splitlines()[0]
    assert "inclusion_3.fasta.gz" in header


def test_contig_missing_from_an_inclusion_genome_is_dropped(stubs, settings):
    stubs(contigs={"ctg1": "ACGT" * 25, "ctg2": "TTGG" * 25},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["ctg2", 0, "40=1X7=1X50="]],
          presence={"inclusion_1": ["ctg1", "ctg2"], "inclusion_2": ["ctg1"]})
    scenario = settings()
    assert run(scenario) == 0
    rows = (scenario.output / "4_blast" / "inclusion_blast_hits.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["contig", "inclusion_1.fasta", "inclusion_2.fasta"]
    assert rows[1:] == ["ctg1\t1\t1", "ctg2\t1\t0"]
    assert set(final_records(scenario.output)) == {"ctg1"}


def test_contig_absent_from_every_exclusion_genome_is_kept_as_it_is(stubs, settings):
    stubs(sam=[["ctg1", 4, "*"]], exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = settings()
    assert run(scenario) == 0
    records = final_records(scenario.output)
    assert records["ctg1"] == ("ACGT" * 25).lower()  # All of it differs from the exclusion group
    assert "100I" in (scenario.output / "final_kmers.fasta").read_text()


def test_variants_too_far_apart_are_dropped(stubs, settings):
    # One mismatch at offset 0 and one at offset 30: no primer or probe can cover both
    stubs(exclusion_hits={genome: [["ctg1", 1, 40, 1e-30, "A" * 40, "C" + "A" * 28 + "C" + "A" * 10]]
                          for genome in ("exclusion_1", "exclusion_2")})
    scenario = settings()
    assert run(scenario) == 0
    assert final_records(scenario.output) == {}


def test_variant_in_only_one_exclusion_genome_is_dropped(stubs, settings):
    """A position has to differ in at least 90% of the exclusion genomes the contig hits."""
    stubs(exclusion_hits={
        "exclusion_1": [["ctg1", 1, 20, 1e-30, "ACGTACGTACGTACGTACGT", "ACGTAAGTACGTTCGTACGT"]],
        "exclusion_2": [["ctg1", 1, 20, 1e-30, "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"]],
    })
    scenario = settings()
    assert run(scenario) == 0
    assert final_records(scenario.output) == {}


def test_no_candidate_after_mapping(stubs, settings):
    stubs(sam=[["ctg1", 0, "100="]])  # Identical to the exclusion genome
    with pytest.raises(PrimerFinderError, match="No contig differs enough"):
        run(settings())


def test_no_contig_in_all_inclusion_genomes(stubs, settings):
    stubs(presence={"inclusion_1": ["ctg1"], "inclusion_2": []})
    with pytest.raises(PrimerFinderError, match="No candidate contig is present in all"):
        run(settings())


def test_no_specific_kmer(stubs, settings):
    stubs(kmers=[])
    with pytest.raises(PrimerFinderError, match="No inclusion-specific kmer"):
        run(settings())


def test_a_failing_program_is_reported(stubs, settings):
    stubs(fail="kmc")
    with pytest.raises(PrimerFinderError, match="(?s)kmc failed .exit code 1.*deliberate"):
        run(settings())


def test_empty_assembly(stubs, settings):
    stubs(contigs={})
    with pytest.raises(PrimerFinderError, match="could not assemble"):
        run(settings())


def test_missing_program(settings, monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", str(tmp_path / "nothing"))
    with pytest.raises(PrimerFinderError, match="Required program.*conda install"):
        run(settings())


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"kmer_size": 0}, "between 1 and 256"),
        ({"kmer_size": 300}, "between 1 and 256"),
        ({"duplication": 0}, "must be 1 or more"),
        ({"assembler": "velvet"}, "Unknown assembler"),
    ],
)
def test_rejected_settings(stubs, settings, changes, message):
    with pytest.raises(PrimerFinderError, match=message):
        run(settings(**changes))


def test_reference_must_be_an_exclusion_genome(stubs, settings, genomes):
    inclusion, _ = genomes
    with pytest.raises(PrimerFinderError, match="must be one of the genomes in the exclusion folder"):
        run(settings(reference=inclusion / "inclusion_1.fasta"))


def test_reference_can_be_chosen(stubs, settings, genomes):
    _, exclusion = genomes
    scenario = settings(reference=Path(exclusion / "exclusion_2.fasta"))
    run(scenario)
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert Path(info["reference"]).name == "exclusion_2.fasta"


def test_the_same_genome_in_both_groups_is_refused(stubs, settings, genomes, tmp_path):
    inclusion, exclusion = genomes
    (exclusion / "inclusion_1.fasta").symlink_to(inclusion / "inclusion_1.fasta")
    with pytest.raises(PrimerFinderError, match="in both groups"):
        run(settings())


def test_empty_input_folder(stubs, tmp_path, genomes):
    inclusion, exclusion = genomes
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(PrimerFinderError, match="no file with an accepted extension"):
        run(Settings(inclusion=empty, exclusion=exclusion, output=tmp_path / "out", threads=1, memory_gb=2))


def test_missing_input_folder(stubs, tmp_path, genomes):
    _, exclusion = genomes
    with pytest.raises(PrimerFinderError, match="inclusion folder does not exist"):
        run(Settings(inclusion=tmp_path / "gone", exclusion=exclusion, output=tmp_path / "out",
                     threads=1, memory_gb=2))


def test_a_file_that_is_not_a_fasta_is_refused(stubs, settings, genomes):
    inclusion, _ = genomes
    (inclusion / "notes.fasta").write_text("this is not a fasta file\n")
    with pytest.raises(PrimerFinderError, match="1 file.* in the inclusion folder are not fasta"):
        run(settings())


@pytest.mark.parametrize("where", ["same", "below", "above"])
def test_the_output_folder_may_not_overlap_an_input_folder(stubs, settings, genomes, tmp_path, where):
    """The input folders are searched recursively, so results written inside one would come back as input."""
    inclusion, _ = genomes
    output = {"same": inclusion, "below": inclusion / "results", "above": tmp_path}[where]
    with pytest.raises(PrimerFinderError, match="output folder and the inclusion folder overlap"):
        run(settings(output=output))


def test_lower_positions():
    assert lower_positions("ACGTACGT", [0, 3, 99]) == "aCGtACGT"


def test_write_presence_table(tmp_path):
    path = tmp_path / "hits.tsv"
    write_presence_table(path, {"b": {"g2.fasta": True}, "a": {"g1.fasta": False, "g2.fasta": True}})
    assert path.read_text() == "contig\tg1.fasta\tg2.fasta\na\t0\t1\nb\t0\t1\n"


def test_several_genomes_are_blasted_in_parallel(stubs, settings, tmp_path):
    """The blast steps run one genome per thread; the results must not be mixed up."""
    inclusion = tmp_path / "inclusion"
    for number in range(3, 9):
        (inclusion / f"inclusion_{number}.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
    stubs(presence={f"inclusion_{number}": ["ctg1"] for number in range(1, 9)})
    scenario = settings(threads=4)
    assert run(scenario) == 0
    rows = (scenario.output / "4_blast" / "inclusion_blast_hits.tsv").read_text().splitlines()
    assert len(rows[0].split("\t")) == 1 + 8  # The contig column and one per genome
    assert rows[1] == "ctg1\t" + "\t".join(["1"] * 8)


def test_two_genomes_with_the_same_file_name(stubs, genomes, tmp_path):
    """Per-sample assembler folders give several genomes the same file name; they must stay apart."""
    _, exclusion = genomes
    inclusion = tmp_path / "samples"
    for sample in ("sampleA", "sampleB"):
        (inclusion / sample).mkdir(parents=True)
        (inclusion / sample / "contigs.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
    stubs(contigs={"ctg1": "ACGT" * 25, "ctg2": "TTGG" * 25},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["ctg2", 0, "40=1X7=1X50="]],
          presence={"sampleA/contigs": ["ctg1", "ctg2"], "sampleB/contigs": ["ctg1"]})
    output = tmp_path / "out"
    assert run(Settings(inclusion=inclusion, exclusion=exclusion, output=output,
                        threads=2, memory_gb=2, keep_intermediate=True)) == 0
    rows = (output / "4_blast" / "inclusion_blast_hits.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["contig", "sampleA/contigs.fasta", "sampleB/contigs.fasta"]
    assert rows[1:] == ["ctg1\t1\t1", "ctg2\t1\t0"]
    # Each genome got a folder of its own, named after its position in the list
    folders = sorted(path.name for path in (output / "4_blast" / "inclusion_db").iterdir())
    assert folders == ["0000_contigs", "0001_contigs"]


def test_genome_labels(tmp_path):
    from primer_finder.pipeline import genome_labels

    root = tmp_path / "inclusion"
    labels = genome_labels([root / "a.fasta", root / "sub" / "a.fasta"], root)
    assert sorted(labels.values()) == ["a.fasta", "sub/a.fasta"]


def test_an_output_path_with_a_space_is_refused_before_anything_runs(stubs, settings, tmp_path):
    """BLAST cannot open a database whose path holds a space; failing late would waste the whole run."""
    with pytest.raises(PrimerFinderError, match="output folder path contains a space"):
        run(settings(output=tmp_path / "with space" / "out"))


def test_input_folders_may_hold_a_space(stubs, genomes, tmp_path):
    """Each genome is linked into the output folder, so blast never sees the input path."""
    import shutil

    inclusion, exclusion = genomes
    spaced = tmp_path / "my genomes"
    shutil.copytree(inclusion, spaced)
    output = tmp_path / "out"
    assert run(Settings(inclusion=spaced, exclusion=exclusion, output=output,
                        threads=1, memory_gb=2)) == 0
    assert (output / "final_kmers.fasta").exists()


@pytest.fixture
def four_inclusion(stubs, genomes, tmp_path):
    """Four inclusion genomes, and stubs where ctg1 is missing from one of them and ctg2 from two."""
    _, exclusion = genomes
    inclusion = tmp_path / "four"
    inclusion.mkdir()
    for number in range(1, 5):
        (inclusion / f"g{number}.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
    stubs(contigs={"ctg1": "ACGT" * 25, "ctg2": "TTGG" * 25},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["ctg2", 0, "40=1X7=1X50="]],
          presence={"g1": ["ctg1", "ctg2"], "g2": ["ctg1", "ctg2"], "g3": ["ctg1"], "g4": []})
    return lambda **changes: Settings(inclusion=inclusion, exclusion=exclusion,
                                      output=tmp_path / "out", threads=1, memory_gb=2, **changes)


@pytest.mark.parametrize(
    ("fraction", "genomes_given", "expected"),
    [(1.0, 8, 8), (0.9, 8, 8), (0.75, 8, 6), (0.5, 3, 2), (0.01, 8, 1), (1.0, 1, 1)],
)
def test_required_genomes_rounds_up(tmp_path, fraction, genomes_given, expected):
    settings = Settings(inclusion=tmp_path, exclusion=tmp_path, output=tmp_path, threads=1, memory_gb=2,
                        min_inclusion=fraction)
    assert required_genomes(settings, [tmp_path] * genomes_given) == expected


def test_a_contig_some_inclusion_genomes_lack_is_kept_with_a_threshold(stubs, four_inclusion):
    """ctg1 is in 3 of the 4 genomes, ctg2 in 2: -p 0.75 keeps ctg1 only."""
    scenario = four_inclusion(min_inclusion=0.75)
    assert run(scenario) == 0
    records = {r.name: r for r in iter_records(scenario.output / "final_kmers.fasta")}
    assert set(records) == {"ctg1"}
    assert "inclusion=3/4" in records["ctg1"].desc
    assert "[5, 12]" in records["ctg1"].desc  # The variant positions are still there


def test_the_default_needs_every_inclusion_genome(stubs, four_inclusion):
    with pytest.raises(PrimerFinderError, match="present in all the inclusion genomes"):
        run(four_inclusion())


def test_a_lower_threshold_keeps_more_and_puts_the_widest_first(stubs, four_inclusion):
    """ctg_a is in 2 of the 4 genomes and ctg_b in 3. The mapping step sorts equal cigars by name, so
    ctg_a comes first until the presence counts reorder them.
    """
    stubs(contigs={"ctg_a": "ACGT" * 25, "ctg_b": "TTGG" * 25},
          sam=[["ctg_a", 0, "40=1X7=1X50="], ["ctg_b", 0, "40=1X7=1X50="]],
          presence={"g1": ["ctg_a", "ctg_b"], "g2": ["ctg_a", "ctg_b"], "g3": ["ctg_b"], "g4": []})
    scenario = four_inclusion(min_inclusion=0.5)
    assert run(scenario) == 0
    candidates = [r.name for r in iter_records(scenario.output / "3_candidates" / "best_kmers.fasta")]
    assert candidates == ["ctg_a", "ctg_b"]  # Before this step, the narrower one is first
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    assert [r.name for r in records] == ["ctg_b", "ctg_a"]  # Most inclusion genomes first
    assert "inclusion=3/4" in records[0].desc
    assert "inclusion=2/4" in records[1].desc


def test_the_threshold_sets_kmc_minimum_count(stubs, four_inclusion):
    calls = stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
                  presence={f"g{n}": ["ctg1"] for n in range(1, 5)})
    scenario = four_inclusion(min_inclusion=0.75)
    run(scenario)
    assert "-ci3 -cx4" in calls.read_text()  # Three of four genomes, and -d 1 still caps at four


def test_run_info_records_the_threshold(stubs, four_inclusion):
    scenario = four_inclusion(min_inclusion=0.75)
    run(scenario)
    parameters = json.loads((scenario.output / "run_info.json").read_text())["parameters"]
    assert parameters["min_inclusion"] == 0.75
    assert parameters["inclusion_genomes_required"] == 3


def test_rounding_up_can_still_mean_every_genome(stubs, four_inclusion):
    """0.9 of 4 genomes rounds up to 4, so the run asks for all of them and says so."""
    with pytest.raises(PrimerFinderError, match="present in all the inclusion genomes"):
        run(four_inclusion(min_inclusion=0.9))


def test_the_error_names_the_threshold_that_was_asked_for(stubs, four_inclusion):
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          presence={"g1": ["ctg1"], "g2": [], "g3": [], "g4": []})
    with pytest.raises(PrimerFinderError, match="at least 3 of the inclusion genomes"):
        run(four_inclusion(min_inclusion=0.75))


@pytest.mark.parametrize("fraction", [0, -0.5, 1.5])
def test_a_threshold_outside_the_range_is_refused(stubs, settings, fraction):
    with pytest.raises(PrimerFinderError, match="must be a fraction greater than 0 and at most 1"):
        run(settings(min_inclusion=fraction))


@pytest.mark.parametrize(
    ("fraction", "genomes_given", "expected"),
    [
        # A fraction is a little more than itself in binary: 0.14 * 50 is 7.000000000000001, and asking
        # for 8 of 50 genomes instead of 7 would be wrong.
        (0.14, 50, 7),
        (0.07, 100, 7),
        (0.29, 100, 29),
        (0.1, 10, 1),
        (1 / 3, 3, 1),
        (0.999, 1000, 999),
    ],
)
def test_required_genomes_is_not_thrown_off_by_floating_point(tmp_path, fraction, genomes_given, expected):
    settings = Settings(inclusion=tmp_path, exclusion=tmp_path, output=tmp_path, threads=1, memory_gb=2,
                        min_inclusion=fraction)
    assert required_genomes(settings, [tmp_path] * genomes_given) == expected


def test_required_genomes_matches_exact_arithmetic(tmp_path):
    """Every two-decimal fraction, against the same sum done without binary floating point."""
    import math
    from fractions import Fraction

    for cents in range(1, 101):
        settings = Settings(inclusion=tmp_path, exclusion=tmp_path, output=tmp_path, threads=1,
                            memory_gb=2, min_inclusion=cents / 100)
        for number in (1, 3, 7, 8, 12, 50, 97, 100, 250):
            expected = max(1, math.ceil(Fraction(cents, 100) * number))
            assert required_genomes(settings, [tmp_path] * number) == expected, (cents, number)


def test_a_fraction_too_small_to_round_still_asks_for_one_genome(tmp_path):
    """-p 1e-12 of 4 genomes rounds to nothing; one genome is the floor, or a region in none of them
    would pass."""
    settings = Settings(inclusion=tmp_path, exclusion=tmp_path, output=tmp_path, threads=1, memory_gb=2,
                        min_inclusion=1e-12)
    assert required_genomes(settings, [tmp_path] * 4) == 1
    assert required_genomes(settings, [tmp_path]) == 1


def test_a_contig_no_inclusion_genome_holds_is_dropped_whatever_the_threshold(stubs, four_inclusion):
    """At least one genome has to hold it: a region in none of them is not a region of this group."""
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          presence={f"g{number}": [] for number in range(1, 5)})
    with pytest.raises(PrimerFinderError, match="at least 1 of the inclusion genomes"):
        run(four_inclusion(min_inclusion=1e-12))


def test_no_inclusion_tag_at_the_default(stubs, four_inclusion):
    """With every genome required, the records are equal and nothing is added to their description."""
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          presence={f"g{number}": ["ctg1"] for number in range(1, 5)})
    scenario = four_inclusion()
    assert run(scenario) == 0
    for name in ("4_blast/all_inclusion_contigs.fasta", "final_kmers.fasta"):
        assert "inclusion=" not in (scenario.output / name).read_text()
    parameters = json.loads((scenario.output / "run_info.json").read_text())["parameters"]
    assert parameters["min_inclusion"] == 1.0
    assert parameters["inclusion_genomes_required"] == 4


def test_one_inclusion_genome_and_a_low_threshold(stubs, genomes, tmp_path):
    """Half of one genome is still one genome, so nothing is partial and no tag is written."""
    _, exclusion = genomes
    inclusion = tmp_path / "one"
    inclusion.mkdir()
    (inclusion / "only.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]], presence={"only": ["ctg1"]})
    output = tmp_path / "out"
    assert run(Settings(inclusion=inclusion, exclusion=exclusion, output=output, threads=1, memory_gb=2,
                        min_inclusion=0.5)) == 0
    assert "inclusion=" not in (output / "final_kmers.fasta").read_text()


def test_contigs_in_the_same_number_of_genomes_keep_the_order_of_the_mapping_step(stubs, four_inclusion):
    """The sort is stable: contigs that cover as much as each other are not shuffled."""
    stubs(contigs={"ctg_a": "ACGT" * 25, "ctg_b": "TTGG" * 25, "ctg_c": "GGAA" * 25},
          sam=[["ctg_a", 0, "40=1X7=1X50="], ["ctg_b", 0, "40=1X7=1X50="], ["ctg_c", 0, "40=1X7=1X50="]],
          presence={"g1": ["ctg_a", "ctg_b", "ctg_c"], "g2": ["ctg_a", "ctg_b", "ctg_c"],
                    "g3": ["ctg_a", "ctg_c"], "g4": []})
    scenario = four_inclusion(min_inclusion=0.5)
    assert run(scenario) == 0
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    # ctg_a and ctg_c are in three genomes, ctg_b in two; a and c keep their relative order
    assert [record.name for record in records] == ["ctg_a", "ctg_c", "ctg_b"]


def test_the_tag_survives_on_a_contig_no_exclusion_genome_hits(stubs, four_inclusion):
    """That contig keeps its own description instead of a list of positions; the tag must still be there."""
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 4, "*"]],
          presence={"g1": ["ctg1"], "g2": ["ctg1"], "g3": ["ctg1"], "g4": []},
          exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = four_inclusion(min_inclusion=0.75)
    assert run(scenario) == 0
    record = next(iter_records(scenario.output / "final_kmers.fasta"))
    assert "inclusion=3/4" in record.desc
    assert "100I" in record.desc  # And what the mapping step said about it


def test_the_presence_table_lists_the_candidates_that_were_dropped(stubs, four_inclusion):
    """The table is the place to look when a region disappears, so it is written before the filtering."""
    scenario = four_inclusion(min_inclusion=0.75)
    run(scenario)
    rows = (scenario.output / "4_blast" / HITS_NAME).read_text().splitlines()
    assert rows[0].split("\t") == ["contig", "g1.fasta", "g2.fasta", "g3.fasta", "g4.fasta"]
    assert rows[1:] == ["ctg1\t1\t1\t1\t0", "ctg2\t1\t1\t0\t0"]  # ctg2 was dropped, and is still here
    assert "ctg2" not in (scenario.output / "final_kmers.fasta").read_text()


def test_the_threshold_and_the_duplication_limit_together(stubs, four_inclusion):
    calls = stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
                  presence={f"g{number}": ["ctg1"] for number in range(1, 5)})
    run(four_inclusion(min_inclusion=0.5, duplication=2))
    assert "-ci2 -cx8" in calls.read_text()  # Two of four genomes, up to twice in each of the four
