"""The pipeline from end to end, with the external programs replaced by stubs."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.blast import MOST_SUBJECTS, PRESENCE_FIELDS
from primer_finder.design import MIN_AMPLICON_COVERAGE
from primer_finder.pipeline import (
    ABSENCE_WANTS_GENOMES,
    HITS_NAME,
    MOSTLY_PRESENT,
    Settings,
    lower_positions,
    report_partial_presence,
    required_genomes,
    run,
    write_presence_table,
)
from primer_finder.seqio import Record, iter_records

# The whole of the default contig, differing from it at offsets 5 and 12: two differences
# close enough to share one primer, which is what the strict rule asks for.
DIFFERS_AT_5_AND_12 = "ACGTA" + "A" + "GTACGT" + "T" + "CGTACGT" + "ACGT" * 20


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


def test_a_run_that_gives_up_still_removes_its_databases(stubs, settings):
    """A KMC database for a few hundred exclusion genomes is a gigabyte, and the run that stopped short is
    exactly the one whose folder nobody goes back to tidy."""
    stubs(kmers=[], shared_kmers=["ACGT" * 5])
    scenario = settings()
    with pytest.raises(PrimerFinderError):
        run(scenario)
    folder = scenario.output / "1_kmers"
    assert not (folder / "exclusion.kmc_pre").exists()
    assert not (folder / "inclusion_specific.kmc_suf").exists()
    assert not (folder / "kmc_work").exists()
    assert not (folder / "inclusion_kmers.txt").exists()


def test_a_collapsed_assembly_is_recorded_as_such(stubs, settings):
    """A run whose assembler got nothing out of the kmers still reports regions, and run_info.json says
    which assembler the contigs actually came from -- otherwise the output looks like a SKESA run."""
    import json
    import random

    # Kmers that tile one 100-base stretch, so that assembling them directly gives a contig the rest of
    # the run can work with. The assembler writes nothing, so the fallback is what produces it.
    rng = random.Random(11)
    stretch = "".join(rng.choice("ACGT") for _ in range(100))
    stubs(contigs={}, kmers=[stretch[at:at + 21] for at in range(len(stretch) - 20)],
          sam=[["unitig_1", 0, "40=1X7=1X51="]])
    scenario = settings(kmer_size=21)
    assert run(scenario) == 0
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["parameters"]["assembler"] == "skesa"
    assert info["parameters"]["assembler_used"] == "unitig"


def test_a_run_that_cannot_move_the_old_results_aside_says_so_and_keeps_everything(stubs, settings):
    """Moving the earlier results aside can fail -- a `.previous` path that is a directory, a folder that
    cannot be written to -- and it used to fail from inside the block that runs `clean_up` on the way out,
    so the error took that earlier run's `--keep-intermediate` files with it."""
    first = settings(keep_intermediate=True)
    assert run(first) == 0
    kept = [first.output / "1_kmers" / "inclusion.kmc_pre",
            first.output / "3_candidates" / "mapping.sam",
            first.output / "4_blast" / "exclusion_db"]
    assert all(path.exists() for path in kept)
    (first.output / "final_kmers.fasta.previous").mkdir()
    # The second run does not ask to keep them, so `clean_up` would delete the first run's if it ran
    with pytest.raises(PrimerFinderError, match="Cannot move the earlier run's final_kmers.fasta aside"):
        run(settings(output=first.output))
    assert all(path.exists() for path in kept), "the earlier run's files were tidied away"
    assert (first.output / "final_kmers.fasta").is_file()  # And its answer is still there


def test_a_refused_run_leaves_an_earlier_runs_files_alone(stubs, settings, tmp_path):
    """Clean-up is armed only once the settings have been accepted. It deletes by fixed path under the
    output folder, so a typo in `-i` used to take a previous run's `--keep-intermediate` files with it."""
    first = settings(keep_intermediate=True)
    assert run(first) == 0
    kept = [first.output / "3_candidates" / "mapping.sam", first.output / "4_blast" / "exclusion_db"]
    assert all(path.exists() for path in kept)
    with pytest.raises(PrimerFinderError, match="inclusion folder does not exist"):
        run(settings(inclusion=tmp_path / "typo", output=first.output))
    assert all(path.exists() for path in kept), "a refused run deleted the previous run's files"


def test_an_output_path_that_is_a_file_says_so(stubs, settings, tmp_path):
    """mkdir's own FileExistsError, and then the clean-up tripping over it, used to surface as
    "Not a directory: .../out.txt/1_kmers/inclusion.kmc_pre"."""
    path = tmp_path / "out.txt"
    path.write_text("not a folder\n")
    with pytest.raises(PrimerFinderError, match="output path is a file, or sits under one"):
        run(settings(output=path))
    # And a path that runs *through* a file, which mkdir reports as NotADirectoryError instead
    with pytest.raises(PrimerFinderError, match="output path is a file, or sits under one"):
        run(settings(output=path / "under"))


def test_a_dangling_symlink_where_the_results_go_is_moved_aside(stubs, settings):
    """`exists()` is False for a dangling symlink, so one named like the results file used to survive
    this step and then break the final write, after the whole run had been done."""
    scenario = settings()
    scenario.output.mkdir(parents=True)
    link = scenario.output / "final_kmers.fasta"
    link.symlink_to(scenario.output / "nothing_here.fasta")
    assert link.is_symlink() and not link.exists()
    assert run(scenario) == 0
    assert (scenario.output / "final_kmers.fasta").is_file()
    assert (scenario.output / "final_kmers.fasta.previous").is_symlink()


def test_a_failed_rerun_does_not_leave_the_previous_answer_behind(stubs, settings, caplog):
    """run_info.json is only written by a run that reaches the end, so a folder was left holding the
    previous run's regions with no record of where they came from -- and `design` pointed at it would
    design on them."""
    first = settings()
    assert run(first) == 0
    assert (first.output / "final_kmers.fasta").is_file()
    stubs(kmers=[], shared_kmers=["ACGT" * 5])  # the second run stops at the kmer step
    with pytest.raises(PrimerFinderError, match="No inclusion-specific kmer"):
        run(settings(output=first.output))
    assert not (first.output / "final_kmers.fasta").exists()
    assert not (first.output / "run_info.json").exists()
    assert "An earlier run's results are in" in caplog.text
    # Moved aside, not destroyed: preventing `design` from using them as this run's does not require
    # throwing them away
    assert (first.output / "final_kmers.fasta.previous").is_file()
    assert (first.output / "run_info.json.previous").is_file()


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


def test_the_assembler_is_told_to_report_a_region_as_short_as_one_kmer(stubs, settings):
    """Both assemblers stop at 200 bases of their own accord, which is longer than a kmer: a
    specific stretch shorter than that would never be reported, and those are commonest where the two
    groups are closely related. The floor is the kmer size instead."""
    calls = stubs()
    run(settings(kmer_size=71))
    assert "--min_contig 71" in calls.read_text()


def test_the_shortest_region_can_be_asked_for(stubs, settings):
    import json

    calls = stubs(contigs={"ctg1": "ACGT" * 50}, sam=[["ctg1", 0, "40=1X7=1X151="]])
    scenario = settings(min_region=150)
    run(scenario)
    assert "--min_contig 150" in calls.read_text()
    assert json.loads((scenario.output / "run_info.json").read_text())["parameters"]["min_region"] == 150


def test_a_region_shorter_than_two_hundred_bases_is_reported(stubs, settings):
    """The whole point of the floor: before it, nothing shorter than 200 bases came out of the assembler
    at all, so no amount of looking downstream would have found it."""
    short = "ACGTT" * 24 + "ACGTACGTACGTACGTACGTACGT"   # 144 bases
    calls = stubs(contigs={"ctg1": short}, sam=[["ctg1", 0, "40=1X7=1X95="]])
    scenario = settings()
    run(scenario)
    assert "--min_contig 99" in calls.read_text()
    assert "ctg1" in (scenario.output / "final_kmers.fasta").read_text()


def test_spades_is_used_when_asked(stubs, settings):
    calls = stubs()
    run(settings(assembler="spades"))
    assert "spades.py --s 1" in calls.read_text()
    assert "--only-assembler" in calls.read_text()


def test_spades_says_so_when_it_cannot_go_as_low_as_asked(stubs, settings, caplog):
    """SPAdes loses a short stretch in graph simplification and has no option for it, so asking for a
    floor below that gets fewer regions than asked for -- with a warning rather than in silence."""
    stubs()
    with caplog.at_level("WARNING"):
        run(settings(assembler="spades", kmer_size=99))
    assert "SPAdes does not report a contig shorter" in caplog.text


def test_a_negative_shortest_region_is_refused(stubs, settings):
    with pytest.raises(PrimerFinderError, match="min-region cannot be negative"):
        run(settings(min_region=-1))


def test_spades_contigs_are_held_to_the_same_floor(stubs, settings):
    """SPAdes has no option for a shortest contig, so the command applies the floor to its output
    instead: either assembler reports the same lengths."""
    stubs(contigs={"ctg1": "ACGT" * 25, "tiny": "ACGTACGTAC"},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["tiny", 0, "10="]])
    scenario = settings(assembler="spades")
    run(scenario)
    assembly = (scenario.output / "2_assembly" / "assembly.fasta").read_text()
    assert "ctg1" in assembly and "tiny" not in assembly


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
    assert rows[1:] == ["ctg1\t1.000\t1.000", "ctg2\t1.000\t0.000"]
    assert set(final_records(scenario.output)) == {"ctg1"}


def test_contig_absent_from_every_exclusion_genome_is_kept_as_it_is(stubs, settings):
    stubs(sam=[["ctg1", 4, "*"]], exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = settings()
    assert run(scenario) == 0
    records = final_records(scenario.output)
    assert records["ctg1"] == ("ACGT" * 25).lower()  # All of it differs from the exclusion group
    assert "100I" in (scenario.output / "final_kmers.fasta").read_text()


def test_variants_too_far_apart_are_dropped(stubs, settings):
    """One mismatch at offset 0 and one at offset 30: no primer or probe can cover both, so the strict
    rule turns the region down. `--no-last-resort` is what asks for the strict answer and nothing else."""
    stubs(exclusion_hits={genome: [["ctg1", 1, 100, 1e-30, "A" * 100, "C" + "A" * 29 + "C" + "A" * 69]]
                          for genome in ("exclusion_1", "exclusion_2")})
    scenario = settings(last_resort=False)
    assert run(scenario) == 0
    assert final_records(scenario.output) == {}


def test_variants_too_far_apart_come_back_as_a_last_resort(stubs, settings, caplog):
    """The same region, with the last-resort pass left on: it is reported, and every way of reading the
    output says it only got through on the weaker rule."""
    import json

    stubs(exclusion_hits={genome: [["ctg1", 1, 100, 1e-30, "A" * 100, "C" + "A" * 29 + "C" + "A" * 69]]
                          for genome in ("exclusion_1", "exclusion_2")})
    scenario = settings()
    assert run(scenario) == 0
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    assert len(records) == 1
    assert "last_resort" in records[0].desc
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["last_resort"] is True
    assert info["parameters"]["min_differences"] == 1
    assert "last resort" in caplog.text.lower()


def test_no_last_resort_keeps_the_error_the_run_would_have_given(stubs, settings):
    """With the pass turned off, a run that finds nothing under the strict rule says so and stops, which
    is what it did before the pass existed."""
    stubs(sam=[["ctg1", 0, "50=1X69="]])  # one mismatch: not two within 25 bases
    with pytest.raises(PrimerFinderError, match="No contig differs enough"):
        run(settings(last_resort=False))


def test_the_last_resort_pass_rescues_a_run_that_step_three_emptied(stubs, settings, caplog):
    """The same single mismatch, with the pass left on: step 3 keeps the contig the second time round."""
    stubs(sam=[["ctg1", 0, "50=1X69="]])
    scenario = settings()
    assert run(scenario) == 0
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    assert len(records) == 1 and "last_resort" in records[0].desc
    assert "single difference" in caplog.text


def test_the_last_resort_pass_does_not_run_when_the_strict_rule_found_something(stubs, settings, caplog):
    """It is a fallback, not a second opinion: a run that works must not be relaxed, and nothing it
    reports may carry the mark."""
    scenario = settings()  # the default stubs give a contig with two differences 7 bases apart
    assert run(scenario) == 0
    import json
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["parameters"]["min_differences"] == 2
    assert info["last_resort"] is False
    assert "last_resort" not in (scenario.output / "final_kmers.fasta").read_text()
    assert "last resort" not in caplog.text.lower()


def test_a_failed_program_is_not_answered_by_asking_for_less(stubs, settings):
    """The reason `NothingSelectiveEnough` has a class of its own. A failed program is a PrimerFinderError
    too, so a retry that caught those as well would run the whole of step 3 a second time before failing
    again. That it failed is not enough to check: minimap2 must have been run once."""
    calls = stubs(fail="minimap2")
    with pytest.raises(PrimerFinderError, match="minimap2 failed"):
        run(settings())
    attempts = [line for line in calls.read_text().splitlines()
                if line.startswith("minimap2 ") and "--version" not in line]
    assert len(attempts) == 1, attempts


def test_the_last_resort_pass_still_needs_a_difference(stubs, settings):
    """It lowers how many differences are wanted, not whether any are. A contig identical to the
    exclusion genome is not an assay at either stringency."""
    stubs(sam=[["ctg1", 0, "120="]])
    with pytest.raises(PrimerFinderError, match="No contig differs enough"):
        run(settings())


def test_the_pass_also_rescues_a_run_that_step_four_emptied(stubs, settings, caplog):
    """The retry is not only for step 3. Here contigs do meet the two-differences rule but are not in
    every inclusion genome, while a single-difference contig is -- so the second pass has something to
    report, and the warning must not claim the rule was what failed."""
    stubs(contigs={"ctg1": "ACGT" * 25, "ctg2": "TTGG" * 25},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["ctg2", 0, "50=1X69="]],
          presence={"inclusion_1": ["ctg1", "ctg2"], "inclusion_2": ["ctg2"]},
          exclusion_hits={genome: [["ctg2", 1, 100, 1e-30, "ACGT" * 25,
                                    "ACGTAAGTACGTACGTACGT"]]
                          for genome in ("exclusion_1", "exclusion_2")})
    scenario = settings()
    assert run(scenario) == 0
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    assert [record.name for record in records] == ["ctg2"]
    assert "last_resort" in records[0].desc
    # The reason quoted is step 4's, and the sentence around it does not blame the span rule. Step 3's
    # own INFO line does mention the rule, so this looks at the warnings only.
    warnings = " ".join(record.getMessage() for record in caplog.records
                        if record.levelname == "WARNING")
    assert "No candidate contig is present in all the inclusion genomes" in warnings
    assert "nothing to report at the usual stringency" in warnings
    assert "two differences within" not in warnings


def test_a_last_resort_region_absent_from_every_exclusion_genome_is_marked_too(stubs, settings):
    """The other branch of step 5: a contig no exclusion genome holds is specific on its own, and in a
    last-resort pass it still only got through step 3 on one difference. It carries the mark."""
    stubs(sam=[["ctg1", 0, "50=1X69="]], exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = settings()
    assert run(scenario) == 0
    records = list(iter_records(scenario.output / "final_kmers.fasta"))
    assert len(records) == 1
    assert "last_resort" in records[0].desc.split()


def test_a_last_resort_pass_that_finds_nothing_says_so_and_claims_nothing(stubs, settings, caplog):
    """`min_differences` records the rule the run ended on; `last_resort` records whether any region came
    back under it. A pass that reports nothing leaves the second false and the file empty."""
    import json

    stubs(sam=[["ctg1", 0, "50=1X69="]],
          exclusion_hits={genome: [["ctg1", 1, 100, 1e-30, "ACGT" * 25,
                                    "ACGT" * 25]]
                          for genome in ("exclusion_1", "exclusion_2")})
    scenario = settings()
    assert run(scenario) == 0
    assert list(iter_records(scenario.output / "final_kmers.fasta")) == []
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["parameters"]["min_differences"] == 1
    assert info["last_resort"] is False
    assert "found nothing either" in caplog.text


def test_a_second_pass_that_fails_says_a_weaker_rule_was_already_tried(stubs, settings):
    """Otherwise the error is word for word the one the strict pass gives, and the user cannot tell that
    anything else was attempted."""
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "50=1X69="]],
          presence={"inclusion_1": ["ctg1"], "inclusion_2": []})
    with pytest.raises(PrimerFinderError, match="last-resort pass found nothing either"):
        run(settings())


def test_variant_in_only_one_exclusion_genome_is_dropped(stubs, settings):
    """A position has to differ in at least 90% of the exclusion genomes the contig hits."""
    stubs(exclusion_hits={
        "exclusion_1": [["ctg1", 1, 100, 1e-30, "ACGT" * 25, DIFFERS_AT_5_AND_12]],
        "exclusion_2": [["ctg1", 1, 100, 1e-30, "ACGT" * 25, "ACGT" * 25]],
    })
    scenario = settings(last_resort=False)  # The rule under test is the strict one
    assert run(scenario) == 0
    assert final_records(scenario.output) == {}


def test_a_contig_one_exclusion_genome_holds_identically_is_dropped(stubs, settings):
    """`genomes_hit == 0` is what makes a contig specific on its own, and it has to mean *every* exclusion
    genome. Here one holds it with no difference at all and another does not hold it; relaxing the rule to
    "all but one" would report it."""
    stubs(exclusion_hits={
        "exclusion_1": [["ctg1", 1, 100, 1e-30, "ACGT" * 25, "ACGT" * 25]],
        "exclusion_2": [],
    })
    scenario = settings(last_resort=False)
    assert run(scenario) == 0
    assert final_records(scenario.output) == {}


def test_a_region_whose_specificity_is_a_block_no_exclusion_genome_holds_is_reported(stubs, settings):
    """End to end: the exclusion genomes align the contig's two ends exactly and nothing in between, which
    is what a repeat-targeted marker looks like once its junctions are assembled with it. The aligned ends
    hold no differences, so the run used to drop the region -- while the mapping step two steps earlier
    had called the same contig a candidate on the very block that is missing."""
    ends = [["ctg1", 1, 20, 1e-30, "ACGT" * 5, "ACGT" * 5],
            ["ctg1", 81, 100, 1e-30, "ACGT" * 5, "ACGT" * 5]]  # 60 bases unaligned, over the floor
    stubs(exclusion_hits={"exclusion_1": ends, "exclusion_2": ends})
    scenario = settings(last_resort=False)
    assert run(scenario) == 0
    records = final_records(scenario.output)
    assert list(records) == ["ctg1"]
    marked = sum(1 for base in records["ctg1"] if base.islower())
    assert marked == 60, marked  # The 60 bases between the two aligned ends


def test_a_region_no_exclusion_genome_holds_records_the_whole_of_itself_as_absent(stubs, settings):
    """`design` reads a region's marked bases and cannot tell a difference from a base that is simply not
    there, so `find` records which are which. A region no exclusion genome hits at all is absent in its
    entirety, and until this was recorded `design` read every one of its bases as a mismatch: on a
    289-base region of the IS711 panel, a 20-base 3'-terminal mismatch run over a primer whose template
    is in none of the genomes."""
    stubs(sam=[["ctg1", 4, "*"]], exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = settings()
    assert run(scenario) == 0
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["absent_stretches"] == {"ctg1": [[0, 99]]}   # the whole 100-base contig
    assert "ctg1" in final_records(scenario.output)


def test_only_a_region_the_absence_carries_is_called_one_that_rests_on_it(stubs, settings, caplog):
    """The warning says the regions "are specific because a stretch of them is missing". A region whose
    differences alone meet the rule is not one of them, however much of it is also unaligned -- counting
    it named the wrong regions as the fragile ones."""
    whole = DIFFERS_AT_5_AND_12
    aligned = whole[:50]
    changed = list(aligned)
    changed[5], changed[12] = "G", "G"        # two differences, 7 apart: enough on their own
    hits = [["ctg1", 1, 50, 1e-30, aligned, "".join(changed)]]
    stubs(contigs={"ctg1": whole}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          exclusion_hits={"exclusion_1": hits, "exclusion_2": hits})
    scenario = settings()
    with caplog.at_level("WARNING"):
        assert run(scenario) == 0
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["absent_stretches"] == {"ctg1": [[50, 99]]}   # still recorded, for design
    assert "specific because a stretch of them is missing" not in caplog.text


def test_a_region_the_absence_does_carry_is_warned_about(stubs, settings, caplog):
    """The same shape with one difference instead of two: now the absent tail is what keeps the region,
    so the exclusion group's size is worth saying."""
    whole = DIFFERS_AT_5_AND_12
    aligned = whole[:50]
    changed = list(aligned)
    changed[5] = "G"                           # one difference: not enough on its own
    hits = [["ctg1", 1, 50, 1e-30, aligned, "".join(changed)]]
    stubs(contigs={"ctg1": whole}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          exclusion_hits={"exclusion_1": hits, "exclusion_2": hits})
    with caplog.at_level("WARNING"):
        assert run(settings()) == 0
    assert "specific because a stretch of them is missing from all 2 exclusion genome(s)" in caplog.text


def test_the_absent_stretches_are_recorded_for_the_design_command(stubs, settings):
    """`design` reads the marked positions off the sequence and cannot tell a base that differs from one
    that is missing. This is the only channel that tells it, so the run has to write it."""
    ends = [["ctg1", 1, 20, 1e-30, "ACGT" * 5, "ACGT" * 5],
            ["ctg1", 81, 100, 1e-30, "ACGT" * 5, "ACGT" * 5]]  # 60 bases unaligned between them
    stubs(exclusion_hits={"exclusion_1": ends, "exclusion_2": ends})
    scenario = settings(last_resort=False)
    assert run(scenario) == 0
    info = json.loads((scenario.output / "run_info.json").read_text())
    assert info["absent_stretches"] == {"ctg1": [[20, 79]]}


def test_a_run_with_nothing_absent_records_nothing(stubs, settings):
    """An ordinary region leaves the record empty rather than absent, so `design` can tell the two apart
    from an older run that never wrote it."""
    assert run(settings(last_resort=False)) == 0
    info = json.loads((settings().output / "run_info.json").read_text())
    assert info["absent_stretches"] == {}


def test_a_region_resting_on_absence_says_what_a_small_group_is_worth(stubs, settings, caplog):
    """The differences rule needs 90% of the covering genomes to agree, so no one genome can carry it. A
    stretch missing from the group is credited on whatever the group says -- with a handful of genomes, one
    alignment gap becomes an absence. Measured on a Listeria complex: against 3 genomes, 95% of the
    stretches reported absent were held by others of its 716."""
    ends = [["ctg1", 1, 20, 1e-30, "ACGT" * 5, "ACGT" * 5],
            ["ctg1", 81, 100, 1e-30, "ACGT" * 5, "ACGT" * 5]]  # 60 bases unaligned between them
    stubs(exclusion_hits={"exclusion_1": ends, "exclusion_2": ends})
    with caplog.at_level("WARNING"):
        assert run(settings(last_resort=False)) == 0
    assert "missing from all 2 exclusion genome(s)" in caplog.text
    assert "Widen the exclusion group" in caplog.text


def test_a_region_resting_on_differences_alone_says_nothing_about_the_group(stubs, settings, caplog):
    """The warning is about absence, so an ordinary region must not draw it."""
    with caplog.at_level("WARNING"):
        assert run(settings(last_resort=False)) == 0
    assert "missing from all" not in caplog.text


def test_an_exclusion_group_big_enough_to_mean_it_draws_no_warning(stubs, settings, caplog, tmp_path):
    """The same absence, said by enough genomes to be worth something, passes without comment: the point
    is what a *small* group's say-so is worth, not that absence is doubtful in itself."""
    exclusion = tmp_path / "many"
    ends = [["ctg1", 1, 20, 1e-30, "ACGT" * 5, "ACGT" * 5],
            ["ctg1", 81, 100, 1e-30, "ACGT" * 5, "ACGT" * 5]]
    hits = {}
    for number in range(1, ABSENCE_WANTS_GENOMES + 1):
        (exclusion / f"exclusion_{number}.fasta").parent.mkdir(parents=True, exist_ok=True)
        (exclusion / f"exclusion_{number}.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
        hits[f"exclusion_{number}"] = ends
    stubs(exclusion_hits=hits)
    with caplog.at_level("WARNING"):
        assert run(settings(exclusion=exclusion, last_resort=False)) == 0
    assert list(final_records(settings(exclusion=exclusion).output)) == ["ctg1"]
    assert "missing from all" not in caplog.text


def test_a_contig_no_exclusion_genome_holds_needs_no_difference(stubs, settings):
    """The mirror: hit by none of them, so it is kept whatever it carries."""
    stubs(exclusion_hits={"exclusion_1": [], "exclusion_2": []})
    scenario = settings(last_resort=False)
    assert run(scenario) == 0
    assert list(final_records(scenario.output)) == ["ctg1"]
    assert "last_resort" not in (scenario.output / "final_kmers.fasta").read_text()


def test_no_candidate_after_mapping(stubs, settings):
    stubs(sam=[["ctg1", 0, "100="]])  # Identical to the exclusion genome
    with pytest.raises(PrimerFinderError, match="No contig differs enough"):
        run(settings())


def test_no_contig_in_all_inclusion_genomes(stubs, settings):
    stubs(presence={"inclusion_1": ["ctg1"], "inclusion_2": []})
    with pytest.raises(PrimerFinderError, match="No candidate contig is present in all"):
        run(settings())


def test_no_specific_kmer(stubs, settings):
    """The inclusion genomes share kmers, but the exclusion group has every one of them."""
    stubs(kmers=[], shared_kmers=["ACGT" * 5])
    with pytest.raises(PrimerFinderError, match="No inclusion-specific kmer"):
        run(settings())


def test_inclusion_genomes_that_share_no_kmer_at_all_say_so(stubs, settings):
    """A different problem from the one above, and the advice that fits it is different too: there was
    nothing for the subtraction to remove."""
    stubs(kmers=[], shared_kmers=[])
    with pytest.raises(PrimerFinderError, match="No kmer at all is shared by the inclusion genomes"):
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
    write_presence_table(path, {"b": {"g2.fasta": 1.0}, "a": {"g1.fasta": 0.0, "g2.fasta": 0.4}})
    assert path.read_text() == "contig\tg1.fasta\tg2.fasta\na\t0.000\t0.400\nb\t0.000\t1.000\n"


def test_several_genomes_are_blasted_in_parallel(stubs, settings, tmp_path):
    """The blast steps run one genome per thread; the results must not be mixed up."""
    inclusion = tmp_path / "inclusion"
    for number in range(3, 9):
        (inclusion / f"inclusion_{number}.fasta").write_text(">chr\n" + "ACGT" * 50 + "\n")
    # A second contig that one genome -- and only that genome -- does not hold. Were every genome given
    # the same answer, the columns would read the same in any order and a mix-up would be invisible.
    stubs(contigs={"ctg1": "ACGT" * 25, "ctg2": "TTGA" * 25},
          sam=[["ctg1", 0, "40=1X7=1X50="], ["ctg2", 0, "40=1X7=1X50="]],
          presence={f"inclusion_{number}": ["ctg1"] if number == 5 else ["ctg1", "ctg2"]
                    for number in range(1, 9)})
    scenario = settings(threads=4)
    assert run(scenario) == 0
    rows = (scenario.output / "4_blast" / "inclusion_blast_hits.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    assert len(header) == 1 + 8  # The contig column and one per genome
    counts = dict(zip(header[1:], rows[2].split("\t")[1:], strict=True))
    assert rows[2].startswith("ctg2\t")
    assert counts["inclusion_5.fasta"] == "0.000"  # The one that lacks it, in its own column
    assert sorted(counts.values()) == ["0.000"] + ["1.000"] * 7
    assert rows[1] == "ctg1\t" + "\t".join(["1.000"] * 8)


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
    assert rows[1:] == ["ctg1\t1.000\t1.000", "ctg2\t1.000\t0.000"]
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
    """Each genome is linked into the output folder, so blast never sees the input path.

    The stubs do not mind a space, so a run that merely finishes proves nothing: what is checked is what
    the programs were handed. BLAST splits its path arguments on whitespace and would read `my` and
    `genomes` as two files."""
    import shutil

    calls = stubs()
    inclusion, exclusion = genomes
    spaced = tmp_path / "my genomes"
    shutil.copytree(inclusion, spaced)
    output = tmp_path / "out"
    assert run(Settings(inclusion=spaced, exclusion=exclusion, output=output,
                        threads=1, memory_gb=2)) == 0
    assert (output / "final_kmers.fasta").exists()
    blasted = [line for line in calls.read_text().splitlines()
               if line.startswith(("blastn", "makeblastdb"))]
    assert blasted
    for line in blasted:
        assert "my genomes" not in line, line


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


def test_a_contig_an_inclusion_genome_holds_part_of_is_kept_and_said_to_be_partial(
        stubs, settings, caplog):
    """Presence is any hit at all, and stays that way: the design step checks each amplicon it proposes
    against every genome one at a time, so a region whose shared part carries a usable assay is worth
    reporting even where one genome holds only some of it. What the run has to do is say so, because
    holding a region and holding an eighth of it were the same tick in the table.

    Measured case behind this: on a draft Brucella group a 338-base region was reported as present in all
    thirty inclusion genomes while one of them held 41 bases of it at the end of a contig.
    """
    stubs(contigs={"ctg1": DIFFERS_AT_5_AND_12},
          sam=[["ctg1", 0, "40=1X7=1X50="]],
          presence_spans={"inclusion_2": {"ctg1": [[1, 12]]}})   # 12 of its 100 bases
    scenario = settings()
    with caplog.at_level("INFO"):
        assert run(scenario) == 0
    assert "ctg1" in final_records(scenario.output)          # kept, which is the point
    rows = (scenario.output / "4_blast" / HITS_NAME).read_text().splitlines()
    held = dict(zip(rows[0].split("\t")[1:], rows[1].split("\t")[1:], strict=True))
    assert held["inclusion_2.fasta"] == "0.120"
    assert held["inclusion_1.fasta"] == "1.000"
    assert "held only in part by an inclusion genome" in caplog.text
    assert "as little as 12 of the 100 bases of ctg1" in caplog.text


def test_nothing_is_said_about_partial_presence_when_every_genome_holds_the_whole_region(
        stubs, settings, caplog):
    """The line is worth printing because it is unusual; on an ordinary run it would be noise."""
    stubs(contigs={"ctg1": DIFFERS_AT_5_AND_12}, sam=[["ctg1", 0, "40=1X7=1X50="]])
    with caplog.at_level("INFO"):
        assert run(settings()) == 0
    assert "held only in part" not in caplog.text


def test_a_genome_that_holds_none_of_a_region_is_not_called_a_partial_one(caplog):
    """Below -p 1.0 a region is kept although some genomes have none of it, and the fraction in its header
    says so. This line is about the genomes that do hold it: counting the others would report every such
    region as held "only in part", which is the regime the filter exists for and the one the end-to-end
    test above cannot reach, since it runs at -p 1.0."""
    records = [Record("r1", "", "A" * 100), Record("r2", "", "A" * 100)]
    presence = {
        "r1": {"g1": 1.0, "g2": 0.0, "g3": 0.3},   # one genome lacks it, another holds 30 of its 100
        "r2": {"g1": 1.0, "g2": 0.0, "g3": 1.0},   # one genome lacks it; the rest hold all of it
    }
    with caplog.at_level("INFO"):
        assert report_partial_presence(presence, records, Path("table.tsv")) == 1
    assert "as little as 30 of the 100 bases of r1" in caplog.text
    assert "r2" not in caplog.text


def test_nothing_is_reported_when_no_kept_region_is_partly_held(caplog):
    with caplog.at_level("INFO"):
        assert report_partial_presence({"r1": {"g1": 1.0}}, [Record("r1", "", "AC")], Path("t")) == 0
    assert "held only in part" not in caplog.text


def test_the_partial_presence_threshold_is_the_one_the_design_step_requires():
    """`MOSTLY_PRESENT` is reported against, never applied, and it is worth reporting against the figure
    the design step does apply: what `design` asks of an amplicon is what makes part of a region worth
    telling from the whole of it. Two independent literals until one of them moves."""
    assert MOSTLY_PRESENT == MIN_AMPLICON_COVERAGE


def test_the_presence_search_asks_for_every_subject_not_the_best_one(stubs, settings):
    """A draft assembly breaks a region across two of its own sequences exactly where a repeat is. With
    one subject the other piece is missing from the fraction, which is the case the fraction is for."""
    calls = stubs(contigs={"ctg1": DIFFERS_AT_5_AND_12}, sam=[["ctg1", 0, "40=1X7=1X50="]])
    assert run(settings()) == 0
    # blast runs in the database's own folder with relative paths, so the search is told apart by the
    # fields it asks for rather than by where it was run. The whole list, quoted as the log writes it:
    # the exclusion search asks for these four and two more, and its first four read the same.
    wanted = "-outfmt '6 " + " ".join(PRESENCE_FIELDS) + "'"
    presence_calls = [line for line in calls.read_text().splitlines()
                      if line.startswith("blastn ") and wanted in line]
    assert len(presence_calls) == 2, presence_calls  # One per inclusion genome
    assert all(f"-max_target_seqs {MOST_SUBJECTS}" in line for line in presence_calls), presence_calls


def test_the_presence_table_lists_the_candidates_that_were_dropped(stubs, four_inclusion):
    """The table is the place to look when a region disappears, so it is written before the filtering."""
    scenario = four_inclusion(min_inclusion=0.75)
    run(scenario)
    rows = (scenario.output / "4_blast" / HITS_NAME).read_text().splitlines()
    assert rows[0].split("\t") == ["contig", "g1.fasta", "g2.fasta", "g3.fasta", "g4.fasta"]
    assert rows[1:] == ["ctg1\t1.000\t1.000\t1.000\t0.000",
                        "ctg2\t1.000\t1.000\t0.000\t0.000"]  # ctg2 was dropped, and is still here
    assert "ctg2" not in (scenario.output / "final_kmers.fasta").read_text()


def test_the_threshold_and_the_duplication_limit_together(stubs, four_inclusion):
    calls = stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
                  presence={f"g{number}": ["ctg1"] for number in range(1, 5)})
    run(four_inclusion(min_inclusion=0.5, duplication=2))
    assert "-ci2 -cx8" in calls.read_text()  # Two of four genomes, up to twice in each of the four


def test_the_exclusion_step_still_matches_contigs_when_a_tag_is_in_the_header(stubs, four_inclusion):
    """With a threshold, the query headers of the exclusion step carry `inclusion=x/y` after the name.

    blast reports the first word of a header as the query name, so the hits still line up with the contigs
    (checked against the real blastn, which returns `ctg1` for `>ctg1 40=1X7=1X50= inclusion=3/4`). If they
    did not, a contig would look as though no exclusion genome holds it and would be kept as fully
    specific, which is the opposite of the truth. Here the two variants are 30 bases apart, so the contig
    must be dropped.
    """
    far_apart = ["ctg1", 1, 100, 1e-30, "A" * 100, "C" + "A" * 29 + "C" + "A" * 69]
    stubs(contigs={"ctg1": "ACGT" * 25}, sam=[["ctg1", 0, "40=1X7=1X50="]],
          presence={"g1": ["ctg1"], "g2": ["ctg1"], "g3": ["ctg1"], "g4": []},
          exclusion_hits={"exclusion_1": [far_apart], "exclusion_2": [far_apart]})
    scenario = four_inclusion(min_inclusion=0.75, last_resort=False)
    assert run(scenario) == 0
    # It reached the exclusion step with the tag on its header, and was dropped there
    assert "inclusion=3/4" in (scenario.output / "4_blast" / "all_inclusion_contigs.fasta").read_text()
    assert list(iter_records(scenario.output / "final_kmers.fasta")) == []
