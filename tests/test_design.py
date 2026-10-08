"""Designing assays with Primer3, and deciding which of them could tell the two groups apart."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError
from primer_finder.design import (
    DEFAULT_PRODUCT_SIZE,
    STRONG_BASE_WEIGHT,
    WEAK_BASE_WEIGHT,
    Assay,
    Conditions,
    DesignSettings,
    Oligo,
    amplicon_copies,
    assays_of,
    build_request,
    check_against_exclusion,
    count_inclusion_copies,
    design,
    forced_requests,
    genome_folders,
    parse_records,
    rank,
    requests_for,
    run,
    specific_by,
    variant_runs,
    write_assays,
)
from primer_finder.seqio import Record


def oligo(start=100, length=20, variants=(), reverse=False) -> Oligo:
    return Oligo("x", "A" * length, start=start, length=length, tm=60.0, gc=50.0,
                 variants=list(variants), reverse=reverse)


def assay(**changes) -> Assay:
    defaults = dict(region="r", number=0, forward=oligo(), reverse=oligo(200, reverse=True), probe=None,
                    product_size=120, penalty=0.5)
    defaults.update(changes)
    return Assay(**defaults)


# ---------------------------------------------------------------- where the differences sit


def test_the_three_prime_end_of_each_primer():
    """The left primer reads left to right, the right primer the other way."""
    assert oligo(100, 20).three_prime == 119
    assert oligo(100, 20, reverse=True).three_prime == 100


@pytest.mark.parametrize(
    ("variants", "reverse", "run_length", "near"),
    [
        ((118, 119), False, 2, 2),      # the last two bases of a left primer
        ((119,), False, 1, 1),
        ((100, 101), True, 2, 2),       # the 3' end of a right primer is its lowest coordinate
        ((119, 118, 117), False, 3, 3),
        ((110, 111), False, 0, 0),      # in the middle: no run at the end, and too far from it
        ((115, 119), False, 1, 2),      # one at the end, one four bases away
        ((), False, 0, 0),
    ],
)
def test_runs_of_differences_at_the_three_prime_end(variants, reverse, run_length, near):
    one = oligo(100, 20, variants=variants, reverse=reverse)
    assert one.terminal_run == run_length
    assert one.near_three_prime == near


def test_variant_runs():
    assert variant_runs("AAacgtAAAAAggAAAAt") == [(2, 4), (11, 2), (17, 1)]
    assert variant_runs("ACGT") == []
    assert variant_runs("acgt") == [(0, 4)]


# ---------------------------------------------------------------- talking to Primer3


def test_build_request_asks_for_a_probe_and_ends_the_record():
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 3)
    assert "SEQUENCE_ID=ctg1" in request
    assert "PRIMER_PICK_INTERNAL_OLIGO=1" in request
    assert f"PRIMER_PRODUCT_SIZE_RANGE={DEFAULT_PRODUCT_SIZE}" in request
    assert "PRIMER_INTERNAL_MIN_TM=62.0" in request  # A probe has to be a usable probe
    assert request.endswith("=\n")


def test_the_gc_clamp_is_asked_for_where_primer3_picks_the_3_prime_end():
    """A G or C at the 3' end holds the primer where extension starts, so it is asked for -- except where
    the 3' end is pinned on a differing base, which is whatever the genomes made it."""
    region = Record("ctg1", "", "A" * 700 + "acg" + "A" * 700)
    requests = dict(requests_for(region, DEFAULT_PRODUCT_SIZE, 1, force_ends=True))
    for request in requests:
        pins_a_primer = "SEQUENCE_FORCE_LEFT_END" in request or "SEQUENCE_FORCE_RIGHT_END" in request
        clamp = next(line for line in request.splitlines() if line.startswith("PRIMER_GC_CLAMP="))
        assert clamp == ("PRIMER_GC_CLAMP=0" if pins_a_primer else "PRIMER_GC_CLAMP=1"), request[:80]


def test_the_gc_clamp_can_be_turned_off():
    plain = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1, gc_clamp=0)
    assert "PRIMER_GC_CLAMP=0" in plain
    assert "PRIMER_GC_CLAMP=2" in build_request(Record("c", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1,
                                                gc_clamp=2)


def test_the_structures_primer3_must_not_leave_in_are_asked_for():
    """An oligo that folds on itself, or pairs with itself or its partner, is spent before it ever reaches
    the template. Primer3 rejects those itself, so the limits have to be in every request."""
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1)
    for tag in ("PRIMER_MAX_HAIRPIN_TH", "PRIMER_INTERNAL_MAX_HAIRPIN_TH"):
        assert f"{tag}=47.0" in request
    for tag in ("PRIMER_MAX_SELF_ANY_TH", "PRIMER_MAX_SELF_END_TH",
                "PRIMER_PAIR_MAX_COMPL_ANY_TH", "PRIMER_PAIR_MAX_COMPL_END_TH",
                "PRIMER_INTERNAL_MAX_SELF_ANY_TH", "PRIMER_INTERNAL_MAX_SELF_END_TH"):
        assert f"{tag}=47.0" in request


def test_the_reaction_the_temperatures_are_predicted_for_is_asked_for():
    """A melting temperature is only meaningful for a given reaction, so the salt, the dNTPs and the
    oligo concentrations go in the request -- the probe at its own concentration, not the primers'."""
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1)
    assert "PRIMER_SALT_MONOVALENT=50.0" in request
    assert "PRIMER_SALT_DIVALENT=3.0" in request
    assert "PRIMER_DNTP_CONC=0.8" in request
    assert "PRIMER_DNA_CONC=250.0" in request
    assert "PRIMER_INTERNAL_SALT_MONOVALENT=50.0" in request
    assert "PRIMER_INTERNAL_DNA_CONC=200.0" in request


def test_the_reaction_can_be_described_differently():
    conditions = Conditions(max_hairpin_tm=40.0, max_dimer_tm=35.0, monovalent_mm=60.0, divalent_mm=5.0,
                            dntp_mm=1.2, primer_nm=500.0, probe_nm=100.0)
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1,
                            conditions=conditions)
    assert "PRIMER_MAX_HAIRPIN_TH=40.0" in request
    assert "PRIMER_MAX_SELF_ANY_TH=35.0" in request
    assert "PRIMER_SALT_MONOVALENT=60.0" in request
    assert "PRIMER_SALT_DIVALENT=5.0" in request
    assert "PRIMER_DNTP_CONC=1.2" in request
    assert "PRIMER_DNA_CONC=500.0" in request
    assert "PRIMER_INTERNAL_DNA_CONC=100.0" in request


def test_the_reaction_reaches_the_forced_requests_too():
    """The requests aimed at a run of differences are the ones most likely to return a marginal oligo, so
    they are the ones that most need the same limits as the rest."""
    region = Record("ctg1", "", "A" * 200 + "cgt" + "A" * 200)
    conditions = Conditions(max_hairpin_tm=42.0, monovalent_mm=70.0)
    requests = forced_requests(region, DEFAULT_PRODUCT_SIZE, 1, conditions=conditions)
    assert requests
    for request in requests:
        assert "PRIMER_MAX_HAIRPIN_TH=42.0" in request
        assert "PRIMER_SALT_MONOVALENT=70.0" in request


def test_forced_requests_aim_at_the_runs_of_differences():
    region = Record("ctg1", "", "A" * 100 + "cgt" + "A" * 100)
    requests = forced_requests(region, DEFAULT_PRODUCT_SIZE, 1)
    joined = "\n".join(requests)
    assert "SEQUENCE_FORCE_LEFT_END=102" in joined   # the last base of the run
    assert "SEQUENCE_FORCE_RIGHT_END=100" in joined  # the first, for the primer that reads backwards
    assert "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST=101" in joined
    assert all("SEQUENCE_ID=ctg1~force_" in request for request in requests)


def test_nothing_is_forced_where_the_rest_of_the_assay_would_not_fit():
    """Primer3 does not fail fast on an impossible constraint: it searches the whole template first.

    A probe forced at base 15 of a region, with no room for a primer before it, kept Primer3 busy for 45
    minutes in the run that led to these guards.
    """
    at_the_start = Record("ctg1", "", "ac" + "A" * 1500)
    assert forced_requests(at_the_start, DEFAULT_PRODUCT_SIZE, 1) == []

    at_the_end = Record("ctg2", "", "A" * 1500 + "ac")
    assert forced_requests(at_the_end, DEFAULT_PRODUCT_SIZE, 1) == []

    # In the middle there is room for all three
    in_the_middle = Record("ctg3", "", "A" * 700 + "acg" + "A" * 700)
    joined = "\n".join(forced_requests(in_the_middle, DEFAULT_PRODUCT_SIZE, 1))
    assert "SEQUENCE_FORCE_LEFT_END=702" in joined
    assert "SEQUENCE_FORCE_RIGHT_END=700" in joined
    assert "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST=701" in joined


def test_a_request_primer3_does_not_answer_is_dropped(stubs, tmp_path, caplog):
    """A slow request must not hold up the rest: it is dropped and counted."""
    import logging

    import primer_finder.design as design_module

    stubs(primer3_sleep=5)
    originals = design_module.PRIMER3_TIMEOUT, design_module.PROBE_TIMEOUT
    design_module.PRIMER3_TIMEOUT = design_module.PROBE_TIMEOUT = 0.5
    try:
        with caplog.at_level(logging.INFO):
            assert design([Record("ctg1", "", "A" * 200)], tmp_path / "work") == []
    finally:
        design_module.PRIMER3_TIMEOUT, design_module.PROBE_TIMEOUT = originals
    assert "ran out of time" in caplog.text


def test_pinning_the_probe_gets_less_time_than_the_rest(stubs, tmp_path):
    """It is the expensive constraint, and the one worth abandoning first."""
    from primer_finder.design import PRIMER3_TIMEOUT, PROBE_TIMEOUT, requests_for

    assert PROBE_TIMEOUT < PRIMER3_TIMEOUT
    region = Record("ctg1", "", "A" * 700 + "acg" + "A" * 700)
    requests = requests_for(region, DEFAULT_PRODUCT_SIZE, 1, force_ends=True)
    probe = [timeout for request, timeout in requests
             if "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST" in request]
    others = [timeout for request, timeout in requests
              if "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST" not in request]
    assert probe and all(timeout == PROBE_TIMEOUT for timeout in probe)
    assert others and all(timeout == PRIMER3_TIMEOUT for timeout in others)


def test_parse_records():
    records = parse_records("SEQUENCE_ID=a\nPRIMER_PAIR_NUM_RETURNED=0\n=\nSEQUENCE_ID=b\n=\n")
    assert [record["SEQUENCE_ID"] for record in records] == ["a", "b"]


def test_assays_of_converts_the_right_primer_coordinates():
    """Primer3 gives the rightmost base of a right primer; everything here counts from the leftmost."""
    region = Record("ctg1", "", "ACGT" * 50)
    fields = {
        "SEQUENCE_ID": "ctg1", "PRIMER_PAIR_NUM_RETURNED": "1",
        "PRIMER_LEFT_0": "10,20", "PRIMER_LEFT_0_SEQUENCE": "A" * 20,
        "PRIMER_RIGHT_0": "119,20", "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
        "PRIMER_PAIR_0_PRODUCT_SIZE": "110", "PRIMER_PAIR_0_PENALTY": "0.25",
    }
    designed = assays_of(region, fields)[0]
    assert designed.forward.start == 10 and not designed.forward.reverse
    assert designed.reverse.start == 100 and designed.reverse.end == 120 and designed.reverse.reverse
    assert designed.probe is None
    assert designed.product_size == 110


def test_assays_of_counts_the_differences_each_oligo_covers():
    region = Record("ctg1", "", "A" * 10 + "cg" + "A" * 108 + "t" + "A" * 50)
    fields = {
        "SEQUENCE_ID": "ctg1", "PRIMER_PAIR_NUM_RETURNED": "1",
        "PRIMER_LEFT_0": "0,20", "PRIMER_LEFT_0_SEQUENCE": "A" * 20,
        "PRIMER_RIGHT_0": "139,20", "PRIMER_RIGHT_0_SEQUENCE": "T" * 20,
        "PRIMER_PAIR_0_PRODUCT_SIZE": "140", "PRIMER_PAIR_0_PENALTY": "0.25",
    }
    designed = assays_of(region, fields)[0]
    assert designed.forward.variants == [10, 11]
    assert designed.reverse.variants == [120]
    assert designed.primer_variants_covered == 3


def test_design_runs_both_kinds_of_request_and_numbers_what_comes_back(stubs, tmp_path):
    region = Record("ctg1", "", "A" * 100 + "cgt" + "A" * 100)
    designed = design([region], tmp_path / "work", how_many=1)
    names = [one.name for one in designed]
    assert names == [f"ctg1_assay{number}" for number in range(len(designed))]
    assert len(designed) > 1  # The plain request plus the forced ones, minus the duplicates
    assert all(one.forward.name.endswith("-F") for one in designed)
    # One of them has its left primer's 3' end on the last base of the run
    assert any(one.forward.end - 1 == 102 for one in designed)


def test_design_drops_the_duplicates_the_two_kinds_of_request_produce(stubs, tmp_path):
    stubs(primer3_pairs=1)
    region = Record("ctg1", "", "A" * 100 + "cgt" + "A" * 100)
    designed = design([region], tmp_path / "work", how_many=1)
    # Two assays may share their primers and differ in their probe, which is not a duplicate
    seen = [(one.forward.sequence, one.reverse.sequence, one.probe.sequence if one.probe else "")
            for one in designed]
    assert len(seen) == len(set(seen))


def test_a_primer3_error_is_reported(stubs, tmp_path):
    stubs(primer3_error="SEQUENCE_TEMPLATE is too short")
    with pytest.raises(PrimerFinderError, match="Primer3 failed on ctg1.*too short"):
        design([Record("ctg1", "", "ACGT" * 50)], tmp_path / "work")


# ---------------------------------------------------------------- the amplicon against the genomes


def test_amplicon_copies_counts_each_full_length_hit(stubs, tmp_path, fasta):
    stubs(amplicon_hits={"g1": 2, "g2": 0})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50}), fasta("g2.fasta", {"chr": "TTTT" * 50})]
    one = assay(amplicon="ACGT" * 10)
    one.region, one.number = "ctg1", 0
    counts, _ = amplicon_copies([one], genomes, tmp_path / "work", threads=2)
    assert counts[one.name] == [2, 0]


def test_a_difference_over_a_g_or_c_weighs_more_than_one_over_an_a_or_t():
    """Breaking a G:C pair costs three hydrogen bonds against two, so at equal numbers it counts for more.
    This only ever breaks a tie; it is reasoning from the chemistry, not a measurement."""
    strong = assay(region="gc", exclusion_genomes_with_amplicon=1,
                   forward=oligo(variants=[118, 119], ))
    strong.forward.variant_bases = {118: "G", 119: "C"}
    weak = assay(region="at", exclusion_genomes_with_amplicon=1, forward=oligo(variants=[118, 119]))
    weak.forward.variant_bases = {118: "A", 119: "T"}
    assert strong.primer_variant_weight == 2.0
    assert weak.primer_variant_weight == 1.0
    assert [one.region for one in rank([weak, strong])] == ["gc", "at"]


def test_the_number_of_differences_still_comes_first():
    """The weight breaks ties; it does not outrank having more differences."""
    three_weak = assay(region="three", exclusion_genomes_with_amplicon=1,
                       forward=oligo(variants=[115, 117, 119]))
    three_weak.forward.variant_bases = {115: "A", 117: "T", 119: "A"}
    two_strong = assay(region="two", exclusion_genomes_with_amplicon=1, forward=oligo(variants=[118, 119]))
    two_strong.forward.variant_bases = {118: "G", 119: "C"}
    assert two_strong.primer_variant_weight > three_weak.primer_variant_weight
    assert [one.region for one in rank([two_strong, three_weak])] == ["three", "two"]


def test_a_difference_whose_base_is_unknown_weighs_in_between():
    unknown = oligo(variants=[119])
    assert unknown.variant_weight == 0.75


def test_the_bases_come_from_the_exclusion_alignments(stubs, tmp_path, fasta):
    """What an oligo has to mismatch is what the exclusion genomes have there, so it is read from their
    own alignments of the amplicon."""
    stubs(amplicon_hits={"g1": 1}, amplicon_mismatch={"3": "G", "7": "A"})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGTACGTACGT", forward=oligo(start=0, length=12, variants=[3, 7]))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.forward.variant_bases == {3: "G", 7: "A"}
    assert one.forward.strong_variants == 1
    assert one.primer_variant_weight == STRONG_BASE_WEIGHT + WEAK_BASE_WEIGHT


def test_the_inclusion_group_is_not_asked_for_the_aligned_sequences(stubs, tmp_path, fasta):
    """Only the exclusion genomes' bases are weighed, so counting copies in the inclusion group does not
    ask blast for the alignments or read them."""
    calls = stubs(amplicon_hits={"g1": 2})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGT" * 10)
    count_inclusion_copies([one], genomes, tmp_path / "incl", threads=1)
    assert one.inclusion_copies == [2]
    blastn = [line for line in calls.read_text().splitlines() if line.startswith("blastn")]
    assert blastn and all("qseq sseq" not in line for line in blastn)  # "qseqid" is not "qseq"

    check_against_exclusion([one], genomes, tmp_path / "excl", threads=1)
    blastn = [line for line in calls.read_text().splitlines() if line.startswith("blastn")]
    assert any("qseq sseq" in line for line in blastn)  # the exclusion group is


def test_check_against_exclusion_and_inclusion_copies(stubs, tmp_path, fasta):
    stubs(amplicon_hits={"g1": 1, "g2": 3})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50}), fasta("g2.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGT" * 10)
    check_against_exclusion([one], genomes, tmp_path / "excl", threads=1)
    assert one.exclusion_genomes_with_amplicon == 2
    count_inclusion_copies([one], genomes, tmp_path / "incl", threads=1)
    assert one.inclusion_copies == [1, 3]
    assert one.min_inclusion_copies == 1


# ---------------------------------------------------------------- which assays are worth having


def test_an_assay_is_usable_when_nothing_can_amplify_or_an_oligo_sits_on_a_difference():
    nothing_to_amplify = assay(exclusion_genomes_with_amplicon=0)
    assert nothing_to_amplify.usable and nothing_to_amplify.specific_by_absence
    assert specific_by(nothing_to_amplify) == "absence"

    on_a_difference = assay(forward=oligo(variants=[118, 119]), exclusion_genomes_with_amplicon=3)
    assert on_a_difference.usable and not on_a_difference.specific_by_absence
    assert specific_by(on_a_difference) == "difference"

    neither = assay(exclusion_genomes_with_amplicon=3)
    assert not neither.usable
    assert specific_by(neither) == "nothing"


def test_the_order_puts_the_assays_that_cannot_amplify_the_wrong_group_first():
    absence = assay(region="a", exclusion_genomes_with_amplicon=0, penalty=2.5)
    many_differences = assay(region="b", exclusion_genomes_with_amplicon=2, penalty=0.1,
                             forward=oligo(variants=[110, 111, 118, 119]))
    one_difference = assay(region="c", exclusion_genomes_with_amplicon=2, penalty=0.1,
                           forward=oligo(variants=[119]))
    ordered = [one.region for one in rank([one_difference, many_differences, absence])]
    assert ordered == ["a", "b", "c"]


def test_more_differences_beat_a_longer_run_at_the_end():
    """The order this tool uses: more differences first, where they sit second. See
    docs/wiki/Designing-assays.md for what that is and is not based on."""
    four = assay(region="four", exclusion_genomes_with_amplicon=1,
                 forward=oligo(variants=[105, 110, 115, 119]))
    two_at_the_end = assay(region="two", exclusion_genomes_with_amplicon=1,
                           forward=oligo(variants=[118, 119]))
    assert [one.region for one in rank([two_at_the_end, four])] == ["four", "two"]


def test_at_equal_numbers_the_run_at_the_end_wins():
    at_the_end = assay(region="end", exclusion_genomes_with_amplicon=1,
                       forward=oligo(variants=[118, 119]))
    in_the_middle = assay(region="middle", exclusion_genomes_with_amplicon=1,
                          forward=oligo(variants=[108, 109]))
    assert [one.region for one in rank([in_the_middle, at_the_end])] == ["end", "middle"]


def test_a_probe_covering_more_differences_does_not_beat_better_chemistry():
    """There is no point forcing a probe onto differences if the assay it belongs to is a poor one."""
    good_chemistry = assay(region="good", exclusion_genomes_with_amplicon=0, penalty=0.2,
                           probe=oligo(130, variants=[]))
    poor_chemistry = assay(region="poor", exclusion_genomes_with_amplicon=0, penalty=3.5,
                           probe=oligo(130, variants=[130, 131, 132]))
    assert [one.region for one in rank([poor_chemistry, good_chemistry])] == ["good", "poor"]
    # Within the same band, the probe that covers differences wins
    same_band = assay(region="same", exclusion_genomes_with_amplicon=0, penalty=0.3,
                      probe=oligo(130, variants=[130, 131]))
    assert [one.region for one in rank([good_chemistry, same_band])] == ["same", "good"]


def test_a_repeated_target_ranks_above_a_single_copy_one():
    """More copies per genome usually means a better limit of detection."""
    repeated = assay(region="many", exclusion_genomes_with_amplicon=0, inclusion_copies=[3, 4])
    single = assay(region="one", exclusion_genomes_with_amplicon=0, inclusion_copies=[1, 1])
    assert [one.region for one in rank([single, repeated])] == ["many", "one"]


# ---------------------------------------------------------------- the command from end to end


@pytest.fixture
def finished_run(genomes, tmp_path):
    """A folder that looks like a finished `find` run: the regions it reported and what it recorded."""
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    region = "A" * 100 + "cgt" + "A" * 100  # one run of differences, in the middle
    (results / "final_kmers.fasta").write_text(f">ctg1 [100, 101, 102]\n{region}\n>ctg2 203I\n{region}\n")
    (results / "run_info.json").write_text(json.dumps({
        "version": "1.1.0",
        "parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)},
    }))
    return results


def design_settings(results: Path, output: Path, **changes) -> DesignSettings:
    defaults = dict(results=results, output=output, threads=1)
    defaults.update(changes)
    return DesignSettings(**defaults)


def test_the_genome_folders_come_from_what_the_run_recorded(finished_run, genomes, tmp_path):
    inclusion, exclusion = genomes
    assert genome_folders(design_settings(finished_run, tmp_path / "out")) == (inclusion, exclusion)


def test_the_genome_folders_can_be_given_instead(finished_run, genomes, tmp_path):
    inclusion, exclusion = genomes
    settings = design_settings(finished_run, tmp_path / "out", inclusion=exclusion, exclusion=inclusion)
    assert genome_folders(settings) == (exclusion, inclusion)


def test_a_run_whose_genomes_have_moved(finished_run, tmp_path):
    (finished_run / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": "/gone/incl", "exclusion": "/gone/excl"}}))
    with pytest.raises(PrimerFinderError, match="inclusion folder of that run is gone"):
        genome_folders(design_settings(finished_run, tmp_path / "out"))


def test_a_run_without_run_info(finished_run, tmp_path):
    (finished_run / "run_info.json").unlink()
    with pytest.raises(PrimerFinderError, match="Give them with -i and -e"):
        genome_folders(design_settings(finished_run, tmp_path / "out"))


def test_design_writes_the_assays_and_the_primer_files(stubs, finished_run, tmp_path):
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output)) == 0
    rows = (output / "assays.tsv").read_text().splitlines()
    assert rows[0].split("\t")[:3] == ["assay", "region", "specific_by"]
    assert len(rows) > 1
    # insilicoPCR cannot read both kinds in one file, so both are written
    assert ">ctg1_assay0-F" in (output / "assays_qpcr.fasta").read_text()
    assert ">ctg1_assay0-P" in (output / "assays_qpcr.fasta").read_text()
    assert ">ctg1_assay0-P" not in (output / "assays_pcr.fasta").read_text()
    # and a script to run it later
    script = (output / "run_insilico_pcr.sh").read_text()
    assert "INSILICO_PCR" in script and "assays_qpcr.fasta" in script
    info = json.loads((output / "design_info.json").read_text())
    assert info["counts"]["regions"] == 2
    assert info["counts"]["selective"] is None  # In silico PCR was not run


def test_design_refuses_a_folder_that_is_not_a_finished_run(stubs, tmp_path):
    with pytest.raises(PrimerFinderError, match="give the output folder of a finished"):
        run(design_settings(tmp_path / "nowhere", tmp_path / "out"))


def test_design_with_no_regions(stubs, finished_run, tmp_path):
    (finished_run / "final_kmers.fasta").write_text("")
    with pytest.raises(PrimerFinderError, match="holds no candidate region"):
        run(design_settings(finished_run, tmp_path / "out"))


def test_design_when_primer3_finds_nothing(stubs, finished_run, tmp_path):
    stubs(primer3_pairs=0)
    with pytest.raises(PrimerFinderError, match="no assay in any region"):
        run(design_settings(finished_run, tmp_path / "out"))


def test_design_when_every_assay_would_amplify_both_groups(stubs, finished_run, tmp_path, genomes):
    """The exclusion genomes hold the amplicon and no oligo sits on a difference: nothing to design on."""
    _, exclusion = genomes
    (finished_run / "final_kmers.fasta").write_text(">ctg1 203I\n" + "A" * 203 + "\n")
    stubs(amplicon_hits={path.stem: 1 for path in exclusion.glob("*.fasta")})
    with pytest.raises(PrimerFinderError, match="would amplify the exclusion genomes too"):
        run(design_settings(finished_run, tmp_path / "out"))


def test_design_runs_in_silico_pcr_and_reports_a_verdict_per_mode(stubs, finished_run, tmp_path,
                                                                  insilico_pcr):
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    header, *rows = (output / "assays.tsv").read_text().splitlines()
    columns = header.split("\t")
    for column in ("qpcr_selective", "pcr_selective", "qpcr_inclusion_amplified",
                   "pcr_exclusion_amplified", "inclusion_total"):
        assert column in columns
    first = dict(zip(columns, rows[0].split("\t"), strict=True))
    assert first["qpcr_selective"] == "yes" and first["pcr_selective"] == "yes"
    assert first["inclusion_total"] == "2" and first["pcr_exclusion_amplified"] == "0"
    info = json.loads((output / "design_info.json").read_text())
    assert info["counts"]["selective"]["qpcr"]["complete"] > 0
    assert info["counts"]["selective"]["pcr"]["complete"] > 0


def test_an_assay_that_amplifies_an_exclusion_genome_is_not_selective(stubs, finished_run, tmp_path,
                                                                      genomes, insilico_pcr):
    inclusion, exclusion = genomes
    sample = sorted(path.name.split(".")[0] for path in exclusion.glob("*.fasta"))[0]
    stubs(insilico_extra={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    rows = list(iter_rows(output / "assays.tsv"))
    bad = next(row for row in rows if row["assay"] == "ctg1_assay0")
    assert bad["pcr_selective"] == "no" and bad["pcr_exclusion_amplified"] == "1"
    # and it is no longer at the top of the list
    assert rows[0]["assay"] != "ctg1_assay0"


def test_an_assay_that_misses_an_inclusion_genome_is_not_selective(stubs, finished_run, tmp_path,
                                                                   genomes, insilico_pcr):
    inclusion, _ = genomes
    sample = sorted(path.name.split(".")[0] for path in inclusion.glob("*.fasta"))[0]
    stubs(insilico_misses={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    run(design_settings(finished_run, output, insilico_pcr=insilico_pcr))
    bad = next(row for row in iter_rows(output / "assays.tsv") if row["assay"] == "ctg1_assay0")
    assert bad["qpcr_inclusion_amplified"] == "1" and bad["qpcr_selective"] == "no"


def iter_rows(path: Path):
    import csv

    with path.open() as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def test_the_structures_of_each_oligo_are_reported(tmp_path):
    """What Primer3 predicted is written out, so that a user who runs a reaction under other conditions
    can see how much room an assay had."""
    one = assay(probe=oligo(140, 22), pair_dimer_tm=8.7, pair_dimer_end_tm=4.1)
    one.forward.hairpin_tm = 35.7
    one.forward.self_dimer_tm = 13.6
    one.probe.hairpin_tm = 33.3
    path = tmp_path / "assays.tsv"
    write_assays(path, [one])
    row = next(iter_rows(path))
    assert row["pair_dimer_tm"] == "8.7" and row["pair_dimer_end_tm"] == "4.1"
    assert row["forward_hairpin_tm"] == "35.7" and row["forward_self_dimer_tm"] == "13.6"
    assert row["probe_hairpin_tm"] == "33.3"
    assert row["reverse_hairpin_tm"] == "0.0"


def test_an_assay_without_a_probe_reports_no_probe_structures(tmp_path):
    path = tmp_path / "assays.tsv"
    write_assays(path, [assay()])
    row = next(iter_rows(path))
    assert row["probe_hairpin_tm"] == "" and row["probe_self_dimer_tm"] == ""


def test_the_structures_primer3_returned_are_read_back():
    fields = {
        "PRIMER_PAIR_NUM_RETURNED": "1",
        "PRIMER_LEFT_0_SEQUENCE": "A" * 20, "PRIMER_LEFT_0": "0,20",
        "PRIMER_LEFT_0_TM": "60.0", "PRIMER_LEFT_0_GC_PERCENT": "50.0",
        "PRIMER_LEFT_0_HAIRPIN_TH": "35.7", "PRIMER_LEFT_0_SELF_ANY_TH": "13.6",
        "PRIMER_RIGHT_0_SEQUENCE": "T" * 20, "PRIMER_RIGHT_0": "119,20",
        "PRIMER_RIGHT_0_TM": "60.0", "PRIMER_RIGHT_0_GC_PERCENT": "50.0",
        "PRIMER_RIGHT_0_HAIRPIN_TH": "0.00", "PRIMER_RIGHT_0_SELF_ANY_TH": "0.00",
        "PRIMER_PAIR_0_PRODUCT_SIZE": "120",
        "PRIMER_PAIR_0_COMPL_ANY_TH": "8.70", "PRIMER_PAIR_0_COMPL_END_TH": "4.10",
    }
    one, = assays_of(Record("ctg1", "", "ACGT" * 50), fields)
    assert one.forward.hairpin_tm == 35.7 and one.forward.self_dimer_tm == 13.6
    assert one.reverse.hairpin_tm == 0.0
    assert one.pair_dimer_tm == 8.7 and one.pair_dimer_end_tm == 4.1


def test_an_old_primer3_that_reports_no_structures_is_not_an_error():
    """Those fields are only there when Primer3 does the thermodynamic alignment; without them an assay
    still has primers, and a zero means nothing was predicted rather than nothing was checked."""
    fields = {
        "PRIMER_PAIR_NUM_RETURNED": "1",
        "PRIMER_LEFT_0_SEQUENCE": "A" * 20, "PRIMER_LEFT_0": "0,20",
        "PRIMER_LEFT_0_TM": "60.0", "PRIMER_LEFT_0_GC_PERCENT": "50.0",
        "PRIMER_RIGHT_0_SEQUENCE": "T" * 20, "PRIMER_RIGHT_0": "119,20",
        "PRIMER_RIGHT_0_TM": "60.0", "PRIMER_RIGHT_0_GC_PERCENT": "50.0",
        "PRIMER_PAIR_0_PRODUCT_SIZE": "120",
    }
    one, = assays_of(Record("ctg1", "", "ACGT" * 50), fields)
    assert one.forward.hairpin_tm == 0.0 and one.pair_dimer_tm == 0.0


def test_write_assays_without_a_verdict(tmp_path):
    one = assay(exclusion_genomes_with_amplicon=0, inclusion_copies=[2, 2])
    path = tmp_path / "assays.tsv"
    write_assays(path, [one])
    row = next(iter_rows(path))
    assert row["specific_by"] == "absence"
    assert row["inclusion_copies_min"] == "2"
    assert "qpcr_selective" not in row


def test_the_threshold_of_a_run_is_read_from_what_it_recorded(finished_run, tmp_path):
    from primer_finder.design import threshold_of

    assert threshold_of(design_settings(finished_run, tmp_path / "out")) == 1.0
    (finished_run / "run_info.json").write_text(json.dumps(
        {"parameters": {"min_inclusion": 0.75, "inclusion": "x", "exclusion": "y"}}))
    assert threshold_of(design_settings(finished_run, tmp_path / "out")) == 0.75


def test_a_run_without_run_info_has_the_default_threshold(finished_run, tmp_path):
    from primer_finder.design import threshold_of

    (finished_run / "run_info.json").unlink()
    assert threshold_of(design_settings(finished_run, tmp_path / "out")) == 1.0


def test_an_assay_missing_an_inclusion_genome_counts_when_the_run_allowed_it(stubs, finished_run,
                                                                             tmp_path, genomes,
                                                                             insilico_pcr):
    """-p 0.5 on the find run means an assay need only reach half the inclusion genomes here."""
    inclusion, exclusion = genomes
    (finished_run / "run_info.json").write_text(json.dumps({"parameters": {
        "inclusion": str(inclusion), "exclusion": str(exclusion), "min_inclusion": 0.5}}))
    sample = sorted(path.name.split(".")[0] for path in inclusion.glob("*.fasta"))[0]
    stubs(insilico_misses={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    row = next(one for one in iter_rows(output / "assays.tsv") if one["assay"] == "ctg1_assay0")
    assert row["qpcr_selective"] == "partial (50%)"
    assert row["qpcr_inclusion_percent"] == "50"
    # and an assay that amplifies them all still ranks above it
    assert next(iter_rows(output / "assays.tsv"))["qpcr_selective"] == "yes"
    info = json.loads((output / "design_info.json").read_text())
    assert info["parameters"]["min_inclusion"] == 0.5
    assert info["counts"]["selective"]["qpcr"]["at_threshold"] >= info["counts"]["selective"]["qpcr"]["complete"]
