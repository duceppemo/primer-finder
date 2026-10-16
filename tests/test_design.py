"""Designing assays with Primer3, and deciding which of them could tell the two groups apart."""

from __future__ import annotations

import json
from itertools import permutations
from pathlib import Path

import pytest

from primer_finder import PrimerFinderError, seqio
from primer_finder.blast import Hit
from primer_finder.design import (
    DEFAULT_MIN_PRIMER_DIFFERENCES,
    DEFAULT_MISMATCHES,
    DEFAULT_PRODUCT_SIZE,
    MIN_AMPLICON_IDENTITY,
    STRONG_BASE_WEIGHT,
    WEAK_BASE_WEIGHT,
    Assay,
    Conditions,
    DesignSettings,
    Oligo,
    amplicon_copies,
    assays_of,
    build_request,
    chains_of,
    check_against_exclusion,
    check_product_size,
    count_inclusion_copies,
    design,
    difference_targets,
    forced_requests,
    genome_folders,
    has_the_amplicon,
    parse_records,
    rank,
    requests_for,
    run,
    specific_by,
    unambiguous_contigs,
    with_flanks,
    write_assays,
)
from primer_finder.seqio import Record, reverse_complement


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


def test_the_targets_are_stretches_one_oligo_could_cover():
    """`find` keeps a region because its differences could sit in one oligo -- a run of them, or two of
    them fewer than a primer's length apart. A pair seven bases apart is not a run, so aiming only at runs
    would ignore the evidence the region was kept for."""
    # start, width, how many differences
    assert difference_targets("A" * 60 + "cgt" + "A" * 60) == [(60, 3, 3)]          # a run of three
    assert difference_targets("A" * 60 + "c" + "A" * 6 + "t" + "A" * 60) == [(60, 8, 2)]   # 7 apart
    assert difference_targets("A" * 60 + "c" + "A" * 8 + "g" + "A" * 9 + "t" + "A" * 60) == [(60, 20, 3)]
    assert difference_targets("ACGT") == []


def test_a_pair_beside_a_richer_target_still_gets_one():
    """The windows are measured again over what is left after each pick. Taking them from one list would
    orphan a pair whose only window overlapped a richer one, and those are the pairs a region was kept
    for."""
    seq = list("A" * 40)
    for position in (0, 5, 24, 30, 35):
        seq[position] = "a"
    targets = difference_targets("".join(seq))
    assert targets == [(24, 12, 3), (0, 6, 2)]      # and not the (24, 12, 3) alone
    assert difference_targets("a" + "A" * 4 + "a") == [(0, 6, 2)]   # the pair on its own, for comparison


def test_a_lone_difference_beside_a_richer_target_is_still_aimed_at():
    seq = list("A" * 40)
    for position in (0, 24, 30):
        seq[position] = "a"
    assert difference_targets("".join(seq)) == [(24, 7, 2), (0, 1, 1)]


def test_only_the_targets_that_will_be_used_are_measured():
    """Measuring the rest costs real time on a region dense with differences, and nothing reads it. The
    limit must not change which targets come first, and the whole thing must stay quick: measuring every
    target of a region like this one took over a minute before the window scan became a bisect."""
    import random
    import time

    random.seed(1)
    bases = list("A" * 20000)
    for position in random.sample(range(20000), 5000):
        bases[position] = "a"
    seq = "".join(bases)
    started = time.monotonic()
    four = difference_targets(seq, limit=4)
    elapsed = time.monotonic() - started
    assert four == difference_targets(seq, limit=0)[:4]
    assert elapsed < 2.0, f"{elapsed:.1f}s for four targets of a dense region"


def test_a_difference_too_far_from_the_next_is_its_own_target():
    """Past the width of the widest oligo there is no placement that covers both, so they are separate
    targets and each is worth only itself."""
    seq = "A" * 60 + "c" + "A" * 30 + "t" + "A" * 60
    assert difference_targets(seq) == [(60, 1, 1), (91, 1, 1)]


def test_the_richest_target_comes_first_and_they_do_not_overlap():
    """Only the first few are aimed at, so the order decides which differences get a request; and two
    targets over the same bases would spend those requests twice on the same thing."""
    seq = "A" * 40 + "c" + "A" * 5 + "g" + "A" * 5 + "t" + "A" * 40 + "cg" + "A" * 40
    targets = difference_targets(seq)
    assert [count for _, _, count in targets] == sorted((c for _, _, c in targets), reverse=True)
    assert targets[0][2] == 3
    covered = [range(start, start + width) for start, width, _ in targets]
    for first, second in zip(covered, covered[1:], strict=False):
        assert not set(first) & set(second)


def test_a_window_wider_than_one_oligo_is_not_a_target():
    """The span is the widest primer Primer3 may return: a stretch wider than that is not one an oligo
    could cover, however many differences are in it."""
    seq = "A" * 60 + "c" + "A" * 23 + "t" + "A" * 60      # 24 apart: a 25-mer covers both
    assert difference_targets(seq) == [(60, 25, 2)]
    seq = "A" * 60 + "c" + "A" * 24 + "t" + "A" * 60      # 25 apart: nothing covers both
    assert difference_targets(seq) == [(60, 1, 1), (85, 1, 1)]


def test_a_short_region_still_gets_its_pinned_requests():
    """The amplicon starts where the primer starts, up to a primer's length before the 3' end that is
    pinned -- not at the 3' end. Demanding the whole amplicon beyond it skipped the pins on exactly the
    short regions the kmer-size floor made reportable."""
    # 106 bases with its differences at 48-50. A left primer pinned at 50 starts at 26 at the earliest,
    # which leaves 80 bases for a 70-base product -- but the old guard asked for 70 bases *past* 50
    region = Record("ctg1", "", "A" * 48 + "cgt" + "A" * 55)
    joined = "\n".join(forced_requests(region, DEFAULT_PRODUCT_SIZE, 1))
    assert "SEQUENCE_FORCE_LEFT_END=50" in joined
    assert "SEQUENCE_FORCE_RIGHT_END=48" in joined


def test_a_pin_with_no_room_for_an_amplicon_is_still_skipped():
    """The guard is looser, not gone: a region only as long as the smallest product has no room for one
    once a primer has to fit behind the pin as well."""
    region = Record("ctg1", "", "acg" + "A" * 40)
    assert forced_requests(region, "120-200", 1) == []


def test_a_right_pin_cannot_run_past_the_end_of_the_region():
    """The amplicon may start up to a primer's length past the pinned 3' end, but not past the region:
    without the clamp this asks Primer3 for a product longer than the template it was given."""
    region = Record("ctg1", "", "A" * 95 + "cg" + "A" * 18)   # 115 bases, differences at 95 and 96
    assert len(region.seq) == 115
    assert forced_requests(region, "120-200", 1) == []


def test_a_region_that_is_all_differences_gets_no_pinned_request():
    """Nothing to aim at: minimap2 could not place it in the exclusion genomes at all, so it is specific
    by absence and pinning an end on an arbitrary base of it only spends requests."""
    assert forced_requests(Record("ctg1", "", "acgt" * 60), DEFAULT_PRODUCT_SIZE, 1) == []


def test_a_forced_request_pins_the_end_past_a_whole_pair_not_just_one_base():
    """The point is to spend both differences on one primer: a left primer whose 3' end is the *second* of
    them has both behind it, under the primer."""
    region = Record("ctg1", "", "A" * 200 + "c" + "A" * 6 + "t" + "A" * 200)
    joined = "\n".join(forced_requests(region, DEFAULT_PRODUCT_SIZE, 1))
    assert "SEQUENCE_FORCE_LEFT_END=207" in joined    # the second difference
    assert "SEQUENCE_FORCE_RIGHT_END=200" in joined   # the first, for the primer that reads backwards
    assert "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST=204" in joined   # the middle of the pair


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
    requests = [request for request, _ in requests_for(region, DEFAULT_PRODUCT_SIZE, 1, force_ends=True)]
    assert len(requests) > 1
    for request in requests:
        pins_a_primer = "SEQUENCE_FORCE_LEFT_END" in request or "SEQUENCE_FORCE_RIGHT_END" in request
        clamp = next(line for line in request.splitlines() if line.startswith("PRIMER_GC_CLAMP="))
        assert clamp == ("PRIMER_GC_CLAMP=0" if pins_a_primer else "PRIMER_GC_CLAMP=1"), request[:80]


def test_the_gc_clamp_can_be_turned_off():
    plain = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1, gc_clamp=0)
    assert "PRIMER_GC_CLAMP=0" in plain
    assert "PRIMER_GC_CLAMP=2" in build_request(Record("c", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1,
                                                gc_clamp=2)


def test_in_silico_pcr_allows_one_mismatch_by_default():
    """insilicoPCR does not count the last two bases of a primer at all, so this governs the mismatches
    further in -- and one of those is often not enough to stop a reaction. Assuming it is would flatter
    every assay that rests on a single difference."""
    assert DEFAULT_MISMATCHES == 1
    assert DesignSettings(results=Path("r"), output=Path("o"), threads=1).mismatches == 1


def test_the_default_tolerance_reaches_insilico_pcr(stubs, finished_run, tmp_path, genomes, insilico_pcr):
    """`insilico.run` defaults to 0 of its own, so a settings value that never arrives would quietly
    tighten the check instead of failing."""
    calls = stubs()
    assert run(design_settings(finished_run, tmp_path / "assays", insilico_pcr=insilico_pcr)) == 0
    asked = [line for line in calls.read_text().splitlines() if line.startswith("insilico_pcr")]
    assert len(asked) == 4  # two modes, two groups
    assert all(" -m 1" in line for line in asked), asked


def test_another_tolerance_reaches_insilico_pcr(stubs, finished_run, tmp_path, genomes, insilico_pcr):
    calls = stubs()
    settings = design_settings(finished_run, tmp_path / "assays", insilico_pcr=insilico_pcr, mismatches=2)
    assert run(settings) == 0
    asked = [line for line in calls.read_text().splitlines() if line.startswith("insilico_pcr")]
    assert all(" -m 2" in line for line in asked), asked


def test_a_failed_check_keeps_the_assays_but_reports_the_failure(stubs, finished_run, tmp_path,
                                                                 insilico_pcr, caplog):
    """Everything before the check is the expensive part of the command. Throwing it away because an
    optional check failed would be the wrong trade, so the table is written and the exit code says so."""
    stubs(insilico_fail=True)
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 1
    assert "without its verdicts" in caplog.text
    rows = list(iter_rows(output / "assays.tsv"))
    assert rows and "pcr_inclusion_amplified" not in rows[0]  # No verdict to report, so no column


# Every field distinct, so that a tag fed from the wrong one cannot pass. The dimer limit covers six tags
# and the hairpin limit two, which is why they are kept apart here: at the defaults both are 47 °C, and a
# test written against the defaults cannot tell them apart.
DISTINCT = Conditions(max_hairpin_tm=41.0, max_dimer_tm=35.0, monovalent_mm=62.0, divalent_mm=5.0,
                      dntp_mm=1.2, primer_nm=500.0, probe_nm=100.0)
DISTINCT_TAGS = {
    "PRIMER_GC_CLAMP": "2",
    "PRIMER_MAX_HAIRPIN_TH": "41.0",
    "PRIMER_INTERNAL_MAX_HAIRPIN_TH": "41.0",
    "PRIMER_MAX_SELF_ANY_TH": "35.0",
    "PRIMER_MAX_SELF_END_TH": "35.0",
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": "35.0",
    "PRIMER_PAIR_MAX_COMPL_END_TH": "35.0",
    "PRIMER_INTERNAL_MAX_SELF_ANY_TH": "35.0",
    "PRIMER_INTERNAL_MAX_SELF_END_TH": "35.0",
    "PRIMER_SALT_MONOVALENT": "62.0",
    "PRIMER_INTERNAL_SALT_MONOVALENT": "62.0",
    "PRIMER_SALT_DIVALENT": "5.0",
    "PRIMER_INTERNAL_SALT_DIVALENT": "5.0",
    "PRIMER_DNTP_CONC": "1.2",
    "PRIMER_INTERNAL_DNTP_CONC": "1.2",
    "PRIMER_DNA_CONC": "500.0",
    "PRIMER_INTERNAL_DNA_CONC": "100.0",   # the probe's own concentration, not the primers'
}


def test_every_tag_the_reaction_decides_carries_the_right_field():
    """An oligo that folds on itself, or pairs with itself or its partner, is spent before it reaches the
    template, and a melting temperature only means anything for a given reaction. Primer3 is told all of
    it, and each tag has to carry the field it belongs to -- a probe fed the primers' concentration, or a
    hairpin limit fed the dimer limit, would be silently wrong."""
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1,
                            gc_clamp=2, conditions=DISTINCT)
    lines = request.splitlines()
    for tag, value in DISTINCT_TAGS.items():
        assert f"{tag}={value}" in lines, tag


def test_the_defaults_are_primer3s_own_structure_limits_and_a_taqman_reaction():
    request = build_request(Record("ctg1", "", "ACGT" * 50), DEFAULT_PRODUCT_SIZE, 1).splitlines()
    assert "PRIMER_MAX_HAIRPIN_TH=47.0" in request
    assert "PRIMER_MAX_SELF_ANY_TH=47.0" in request
    assert "PRIMER_SALT_MONOVALENT=50.0" in request
    assert "PRIMER_SALT_DIVALENT=3.0" in request
    assert "PRIMER_DNTP_CONC=0.8" in request
    assert "PRIMER_DNA_CONC=250.0" in request
    assert "PRIMER_INTERNAL_DNA_CONC=200.0" in request


def test_every_request_of_a_region_carries_the_whole_reaction():
    """Including the ones aimed at a run of differences, which are the most likely to return a marginal
    oligo and so the ones that most need the same limits as the rest. Only the GC clamp differs."""
    region = Record("ctg1", "", "A" * 300 + "cgt" + "A" * 300)
    requests = [request for request, _ in
                requests_for(region, DEFAULT_PRODUCT_SIZE, 1, force_ends=True,
                             gc_clamp=2, conditions=DISTINCT)]
    assert len(requests) > 1
    for request in requests:
        lines = request.splitlines()
        for tag, value in DISTINCT_TAGS.items():
            if tag == "PRIMER_GC_CLAMP":
                continue
            assert f"{tag}={value}" in lines, (tag, request[:60])


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
    # The amplicon is what gets blasted to decide specificity, so it has to be the whole product: from
    # the forward primer's first base to the reverse primer's last, both included.
    assert designed.amplicon == region.seq[10:120]
    assert len(designed.amplicon) == designed.product_size


def test_a_right_primer_keeps_the_strand_primer3_reports_it_on():
    """Primer3 reports a right primer as the oligo to order, which is the reverse complement of the
    template under it. That sequence is what goes into the order sheet and the insilicoPCR files, so it
    is passed through untouched -- while the coordinates still describe the template."""
    # Not a palindrome over the right primer: "ACGT" * n is its own reverse complement, which would make
    # the test pass whichever strand the sequence came back on
    region = Record("ctg1", "", "ACGT" * 25 + "AAAACCCCGGGGTTTT" * 2 + "ACGT" * 25)
    template = region.seq[100:120]
    fields = {
        "SEQUENCE_ID": "ctg1", "PRIMER_PAIR_NUM_RETURNED": "1",
        "PRIMER_LEFT_0": "10,20", "PRIMER_LEFT_0_SEQUENCE": region.seq[10:30],
        "PRIMER_RIGHT_0": "119,20", "PRIMER_RIGHT_0_SEQUENCE": reverse_complement(template),
        "PRIMER_PAIR_0_PRODUCT_SIZE": "110", "PRIMER_PAIR_0_PENALTY": "0.25",
    }
    designed = assays_of(region, fields)[0]
    assert designed.reverse.sequence == reverse_complement(template) != template
    assert region.seq[designed.reverse.start:designed.reverse.end] == template


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
    """The plain request and a forced one can propose the same assay, and only one of them is kept.

    The stub places its left primer so that its 3' end lands on base 29, which is where a run of
    differences at 27-29 forces it, so the plain request and that forced one really do come back with the
    same oligos -- Primer3 asked twice about overlapping neighbourhoods of the same region. Without a raw
    duplicate among the records there would be nothing for the de-duplication to remove, and the
    assertion below would hold with it taken out."""
    stubs(primer3_pairs=1, primer3_probe=False)
    region = Record("ctg1", "", "A" * 27 + "cgt" + "A" * 170)
    designed = design([region], tmp_path / "work", how_many=1)
    seen = [(one.forward.sequence, one.reverse.sequence, one.probe.sequence if one.probe else "")
            for one in designed]
    assert len(seen) == len(set(seen))
    # The duplicate was in what Primer3 returned: more records came back than assays came out
    records = sum(text.count("PRIMER_PAIR_NUM_RETURNED=1")
                  for text in [path.read_text()
                               for path in sorted((tmp_path / "work").glob("*_output.txt"))])
    assert records > len(designed)


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
    assert WEAK_BASE_WEIGHT < unknown.variant_weight < STRONG_BASE_WEIGHT


def hit(size, *, qstart=1, qend=None, gaps=0, mismatches=0, subject="chr", sstart=None, send=None):
    """One blast row the way blast reports it: `length` counts gap columns, `nident` does not."""
    qend = size if qend is None else qend
    aligned = qend - qstart + 1 - gaps
    sstart = qstart if sstart is None else sstart
    send = sstart + (qend - qstart) if send is None else send
    return Hit(query="a", qstart=qstart, qend=qend, evalue=1e-30,
               length=qend - qstart + 1, identical=aligned - mismatches, gaps=gaps,
               identity=100.0 * (aligned - mismatches) / (qend - qstart + 1),
               subject=subject, sstart=sstart, send=send)


def test_a_genome_holding_a_copy_blast_reports_in_two_pieces_is_counted(stubs, tmp_path, fasta):
    """Through blast's own output and `amplicon_copies`, not just the chaining on its own: the genome
    holds the whole amplicon with 300 of its own bases in the middle, which blast reports as two
    pieces."""
    stubs(amplicon_hits={"g1": 1}, amplicon_split=[60, 300])
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 200})]
    one = assay(amplicon="ACGT" * 30, forward=oligo(start=0, length=20))
    counts, _ = amplicon_copies([one], genomes, tmp_path / "work", threads=1)
    assert counts[one.name] == [1]  # One copy, not two and not none


def test_two_hits_starting_at_the_same_base_of_the_amplicon_are_not_one_copy():
    """Two copies, one of which blast aligns further than the other. They start at the same base of the
    amplicon, so the second does not continue the first however far apart they sit in the genome."""
    chains = chains_of([hit(200, qstart=1, qend=30, sstart=1, send=30),
                        hit(200, qstart=1, qend=200, sstart=500, send=699)])
    assert len(chains) == 2
    assert [has_the_amplicon(chain, 200) for chain in chains] == [False, True]


def test_pieces_on_opposite_strands_are_not_one_copy():
    """An inverted repeat: the second piece continues the amplicon and lies close by in the genome, but
    on the other strand, so no single primer pair reads through both."""
    # Again coordinates that would chain were both on one strand: 99 bases along, continuing the amplicon
    chains = chains_of([hit(200, qstart=1, qend=100, sstart=1000, send=901),   # Minus
                        hit(200, qstart=101, qend=200, sstart=1100, send=1199)])  # Plus
    assert len(chains) == 2

@pytest.mark.parametrize("overlap", [1, 4])
def test_the_pieces_of_a_split_copy_usually_overlap_by_a_base_or_two(overlap):
    """Where blast breaks a copy the two pieces rarely abut exactly. Over fifteen real cases -- five
    sequences at three insertion sizes -- the break came back as q1-120/q121-207 (abutting), but also
    q1-120/q120-207, q1-121/q120-207, q1-120/q118-207 and q1-123/q120-207: overlaps of one to four
    bases. Refusing any overlap would leave those copies counted as absent, which is the thing this is
    all here to avoid."""
    pieces = [hit(207, qstart=1, qend=120, sstart=11, send=130),
              hit(207, qstart=121 - overlap, qend=207, sstart=231, send=317 - overlap)]
    chains = chains_of(pieces)
    assert len(chains) == 1
    assert has_the_amplicon(chains[0], 207)


def test_a_chain_counts_a_base_of_the_amplicon_once_however_many_pieces_reach_it():
    """Two pieces overlapping by 30 bases, the first matching throughout and the second carrying 10
    mismatches among its own 30 columns: 80 of the 90 covered columns match, which is under the 90% rule.
    Adding the pieces up without regard to the overlap gives 110 of 120, which passes it."""
    pieces = [hit(100, qstart=1, qend=60, sstart=1, send=60),
              hit(100, qstart=31, qend=90, sstart=31, send=90, mismatches=10)]
    chains = chains_of(pieces)
    assert len(chains) == 1
    assert not has_the_amplicon(chains[0], 100)
    # And the same pieces without the mismatches are held, so it is the identity that decided it
    clean = chains_of([hit(100, qstart=1, qend=60, sstart=1, send=60),
                       hit(100, qstart=31, qend=90, sstart=31, send=90)])
    assert has_the_amplicon(clean[0], 100)


def test_two_split_copies_near_each_other_are_counted_whatever_order_blast_prints_them_in():
    """Blast prints its hits by score, not by position, and each piece joins the first chain that will
    take it. These four pieces are two copies 1,200 bases apart, each split by an insertion; in the order
    blast actually emitted them -- the downstream copy's second piece scored highest -- the second piece
    of one copy joined the chain of the other, leaving one countable copy where there are two."""
    first = [hit(207, qstart=1, qend=120, sstart=501, send=620),
             hit(207, qstart=120, qend=207, sstart=720, send=807, mismatches=1)]
    second = [hit(207, qstart=1, qend=120, sstart=1808, send=1927),
              hit(207, qstart=120, qend=207, sstart=2027, send=2114)]
    blast_order = [first[0], second[0], second[1], first[1]]  # What blastn printed
    counted = [sum(1 for chain in chains_of(list(order)) if has_the_amplicon(chain, 207))
               for order in permutations(blast_order)]
    assert set(counted) == {2}, "the count depends on the order blast printed its hits in"


def test_the_bases_come_from_every_piece_of_a_split_copy(stubs, tmp_path, fasta):
    """What an oligo has to mismatch is read off these alignments, so a copy blast split in two has to be
    read in both halves -- including where the second half is the one carrying the difference."""
    # Bases that really differ: "ACGT" * 30 has G at 10 and A at 80
    stubs(amplicon_hits={"g1": 1}, amplicon_split=[60, 300], amplicon_mismatch={"10": "T", "80": "T"})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 200})]
    one = assay(amplicon="ACGT" * 30, forward=oligo(start=0, length=20, variants=[10]),
                reverse=oligo(start=70, length=20, variants=[80]))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.forward.variant_bases == {10: "T"}   # In the first piece
    assert one.reverse.variant_bases == {80: "T"}   # In the second

def test_a_chain_of_three_pieces_counts_each_base_once():
    """Three pieces of a 100-base amplicon overlapping by 10 each, the last carrying 11 mismatches: 89 of
    the 100 covered columns match, which is under the rule. Counting the overlapping columns twice gives
    109 of 120, which is over it."""
    pieces = [hit(100, qstart=1, qend=40, sstart=1, send=40),
              hit(100, qstart=31, qend=70, sstart=31, send=70),
              hit(100, qstart=61, qend=100, sstart=61, send=100, mismatches=11)]
    chains = chains_of(pieces)
    assert len(chains) == 1
    assert not has_the_amplicon(chains[0], 100)
    # The same three pieces matching throughout are held, so it is the identity that decided it
    pieces[-1] = hit(100, qstart=61, qend=100, sstart=61, send=100)
    assert has_the_amplicon(chains_of(pieces)[0], 100)


def test_the_reach_of_a_chain_is_the_furthest_any_piece_got():
    """Not the last piece's end: a piece that reaches past the one after it must not make the overlap
    with that one read as negative, which would count its columns twice."""
    chain = [hit(100, qstart=1, qend=60, sstart=1, send=60),
             hit(100, qstart=21, qend=50, sstart=21, send=50),    # Inside the reach already
             hit(100, qstart=51, qend=100, sstart=51, send=100, mismatches=11)]
    # 100 columns covered, 11 of them mismatching: under the rule however the pieces are added up
    assert not has_the_amplicon(chain, 100)


@pytest.mark.parametrize(("overlap", "chained"), [(30, True), (31, False)])
def test_how_much_two_pieces_may_overlap_before_they_are_two_copies(overlap, chained):
    """Blast overlaps its pieces by a base or two where it breaks a copy; an overlap much larger than that
    is two things aligning to the same stretch rather than one copy reported in parts."""
    pieces = [hit(207, qstart=1, qend=120, sstart=1, send=120),
              hit(207, qstart=121 - overlap, qend=207, sstart=121 - overlap, send=207)]
    assert (len(chains_of(pieces)) == 1) is chained


@pytest.mark.parametrize(("gap", "chained"), [(2000, True), (2001, False)])
def test_how_far_apart_two_pieces_may_be_before_they_are_two_copies(gap, chained):
    """`MAX_PRODUCT_GAP`, measured from the last base of one piece to the first of the next."""
    pieces = [hit(200, qstart=1, qend=100, sstart=1, send=100),
              hit(200, qstart=101, qend=200, sstart=101 + gap, send=200 + gap)]
    assert (len(chains_of(pieces)) == 1) is chained

def contig_holding(flank: str, body: str, before: int = 100) -> dict:
    """A genome whose amplicon sits at `before` + len(flank), with `flank` as the genome spells it."""
    return {"chr": Record("chr", "", "A" * before + flank + body)}


def test_a_flank_blast_would_not_pay_to_align_is_still_counted():
    """Measured on real output: a 72-base amplicon against a *B. suis* genome came back as `q9-72`, 64 of
    72 bases, because its first eight sit either side of a one-base insertion and a one-base deletion and
    blastn will not pay its gap costs to cross that for four matching bases. megablast aligned all 72.
    64 of 72 is 88.9%, just under the coverage rule, so the genome read as not holding an amplicon it
    holds -- which hands the assay an unearned claim that there is nothing to amplify."""
    body = "ACGTTGCA" * 8                      # The 64 bases blast did align
    amplicon = "TGCCCATA" + body               # ... and the eight it did not
    contigs = contig_holding("TGCCGCAT", body)  # As the genome spells them: one insertion, one deletion
    one = hit(72, qstart=9, qend=72, sstart=109, send=172, mismatches=1)
    assert not has_the_amplicon([one], 72)      # 64 of 72 on its own
    whole = with_flanks([one], 72, amplicon, contigs)
    assert [(piece.qstart, piece.qend) for piece in whole] == [(1, 8), (9, 72)]
    assert whole[0].length == 8  # The flank is eight bases wherever they match
    assert has_the_amplicon(whole, 72)


def test_a_genome_holding_a_flank_blast_left_out_is_counted(stubs, tmp_path, fasta):
    """Through `amplicon_copies` and the genome beside the blast database, not just the two functions on
    their own: a hit covering 85% of the amplicon is under the coverage rule, and the genome holds the
    rest exactly where it would have to be."""
    amplicon = "ACGTTGCA" * 9  # 72 bases
    stubs(amplicon_hits={"g1": 1}, amplicon_covers=0.85)
    genomes = [fasta("g1.fasta", {"chr": amplicon + "TTTTTTTTTT"})]
    one = assay(amplicon=amplicon, forward=oligo(start=0, length=20))
    counts, _ = amplicon_copies([one], genomes, tmp_path / "work", threads=1)
    assert counts[one.name] == [1]


def test_a_flank_can_rescue_a_copy_but_never_take_one_away(stubs, tmp_path, fasta):
    """The guarantee the flank check rests on. It is consulted only for a chain that fell short, so a copy
    blast's own evidence reaches is never reconsidered -- an `absence` is the claim that costs, ranking
    above every difference-based assay and exempt from the differences filter, so the check may only ever
    shed one. Measured over four Brucella targets and 958 assays: two absences shed, none gained."""
    amplicon = "ACGTTGCA" * 12  # 96 bases
    # A hit reaching 95% of the amplicon at 91.2% identity -- over the rule, but only just, so that an
    # unheld flank folded in would take it under (86.5%) if the guard were not there.
    mismatches = {str(position): "A" for position in range(1, 16, 2)}
    stubs(amplicon_hits={"g1": 1}, amplicon_covers=0.95, amplicon_mismatch=mismatches)
    genomes = [fasta("g1.fasta", {"chr": amplicon[:91] + "GGGGGGGGGGGGGGG"})]
    one = assay(amplicon=amplicon, forward=oligo(start=0, length=20))
    counts, _ = amplicon_copies([one], genomes, tmp_path / "work", threads=1)
    assert counts[one.name] == [1], "a chain that reaches far enough on its own must stay a copy"


def test_a_flank_cut_off_by_the_start_of_a_contig_is_compared_in_place():
    """A window the contig's start cuts short holds the bases nearest the alignment -- the *last* of the
    flank. Comparing it against the flank's first bases instead lines up the wrong pairs, and could credit
    a flank the genome does not hold at all."""
    body = "ACGTTGCA" * 8                 # 64 bases the hit covers
    amplicon = "AAAAACCC" + body           # An 8-base flank whose halves are plainly different
    contig = {"chr": seqio.Record("chr", "", "CCC" + body)}  # Only its last three bases are there
    one = hit(72, qstart=9, qend=72, sstart=4, send=67)
    flank = with_flanks([one], 72, amplicon, contig)[0]
    assert flank.length == 8
    assert flank.identical == 3, flank.identical  # "CCC" against "CCC", not "AAA" against "CCC"


def test_a_flank_is_not_taken_from_a_contig_whose_name_is_ambiguous(tmp_path):
    """Blast reports the name, so two records sharing one cannot be told apart: a hit in the first would
    be given a window out of the second. The rescue declines rather than guessing."""
    path = tmp_path / "genome.fasta"
    path.write_text(">c\nACGTACGTACGT\n>c\nTTTTTTTTTTTT\n")
    assert unambiguous_contigs(path) == {}
    path.write_text(">c\nACGTACGTACGT\n>d\nTTTTTTTTTTTT\n")
    assert sorted(unambiguous_contigs(path)) == ["c", "d"]


def test_the_tail_flank_on_the_minus_strand_is_taken_from_before_the_alignment():
    """The mirror of the head case: on the minus strand the amplicon's tail lies *before* the alignment in
    the genome, reverse complemented."""
    body = "ACGTTGCA" * 8
    amplicon = body + "TTTTCCCC"
    # The whole amplicon on the minus strand: the genome read forwards is its reverse complement, so the
    # amplicon's tail sits at the genome's *start* and the genome coordinate falls as the amplicon rises.
    genome = seqio.reverse_complement(amplicon)
    contigs = {"chr": seqio.Record("chr", "", genome)}
    one = hit(72, qstart=1, qend=64, sstart=72, send=9)   # The body alone, 72 down to 9
    flank = with_flanks([one], 72, amplicon, contigs)[-1]
    assert (flank.qstart, flank.qend) == (65, 72)
    assert flank.identical == 8, flank.identical          # Genome 8..1, reverse complemented


def test_the_flanks_of_a_chain_are_taken_from_its_own_ends():
    """A chain of two pieces: the head flank belongs to the first piece and the tail to the last. Using
    one piece's bounds for both would take a window from the wrong place."""
    body = "ACGTTGCA" * 4                 # 32 bases
    amplicon = "AAAA" + body + body + "CCCC"   # 72: flank, two pieces, flank
    genome = "AAAA" + body + "G" * 50 + body + "CCCC"
    contigs = {"chr": seqio.Record("chr", "", genome)}
    first = hit(72, qstart=5, qend=36, sstart=5, send=36)
    second = hit(72, qstart=37, qend=68, sstart=87, send=118)
    whole = with_flanks([first, second], 72, amplicon, contigs)
    assert [(piece.qstart, piece.qend) for piece in whole] == [(1, 4), (5, 36), (37, 68), (69, 72)]
    assert whole[0].identical == 4 and whole[-1].identical == 4


def test_a_flank_is_compared_whatever_case_the_genome_is_written_in():
    """Soft-masked genomes are ordinary input, and a lower-case flank is still the flank."""
    body = "ACGTTGCA" * 8
    amplicon = "TTTTCCCC" + body
    contigs = {"chr": seqio.Record("chr", "", "ttttcccc" + body)}
    one = hit(72, qstart=9, qend=72, sstart=9, send=72)
    assert with_flanks([one], 72, amplicon, contigs)[0].identical == 8

def test_a_flank_the_genome_does_not_hold_earns_little_and_stays_absent():
    """What keeps the above from turning every short hit into a copy: the credit is the bases that match
    where they would have to be, so a flank that is not there does not buy coverage."""
    body = "ACGTTGCA" * 8
    amplicon = "TTTTTTTTTTTTTTTTTTTT" + body   # 20 bases the genome has nothing like
    contigs = contig_holding("CGACGACGACGACGACGACG", body)
    one = hit(84, qstart=21, qend=84, sstart=121, send=184)
    whole = with_flanks([one], 84, amplicon, contigs)
    assert not has_the_amplicon(whole, 84)


def test_a_difference_in_a_flank_is_not_held_against_the_copy():
    """A flank is the amplicon's outermost bases, so a difference among them is a difference at a primer's
    5' end -- the end whose mismatches do not stop the polymerase extending. Counting them against the
    copy's identity denied amplifications that are real, on seven assays of the Brucella genus panel."""
    body = "ACGTTGCA" * 10
    amplicon = "TTTTTTTT" + body
    contigs = contig_holding("TTTTTTTT", body)
    one = hit(88, qstart=9, qend=88, sstart=109, send=188, mismatches=7)   # 73 of 80: 91.3%
    whole = with_flanks([one], 88, amplicon, contigs)
    assert (whole[0].identical, whole[0].length) == (8, 8)  # The flank is there exactly
    assert has_the_amplicon(whole, 88)


def test_a_flank_on_the_minus_strand_is_looked_for_the_other_way_round():
    """On the minus strand the amplicon runs against the genome, so the head flank lies after the
    alignment rather than before it, reverse complemented."""
    body = "ACGTTGCA" * 8
    amplicon = "TGCCCATA" + body
    genome = reverse_complement("A" * 100 + "TGCCGCAT" + body)
    contigs = {"chr": Record("chr", "", genome)}
    # The body occupies the genome's first 64 bases, descending as the amplicon ascends
    one = hit(72, qstart=9, qend=72, sstart=64, send=1)
    whole = with_flanks([one], 72, amplicon, contigs)
    assert whole[0].length == 8 and whole[0].identical == 4  # TGCCCATA against TGCCGCAT, in place

def test_the_pieces_blast_splits_a_copy_into_are_read_as_one_copy():
    """Measured with real blastn 2.17: a 207-base amplicon against a genome carrying 100 extra bases in
    the middle of it comes back as two pieces, q1-120 at s11-130 and q121-207 at s231-317, each 100%
    identical. Judged one piece at a time neither reaches 90% of the amplicon -- 58% and 42% -- so the
    genome counted as holding nothing, while holding every base of the amplicon and both primer sites."""
    pieces = [hit(207, qstart=1, qend=120, sstart=11, send=130),
              hit(207, qstart=121, qend=207, sstart=231, send=317)]
    chains = chains_of(pieces)
    assert len(chains) == 1
    assert has_the_amplicon(chains[0], 207)
    assert not any(has_the_amplicon([piece], 207) for piece in pieces)  # Which is the whole point


def test_a_copy_split_by_a_deletion_is_read_as_one_copy():
    """The mirror: the genome is missing 30 bases of the amplicon, so the gap is in the query instead."""
    chains = chains_of([hit(150, qstart=1, qend=60, sstart=1, send=60),
                        hit(150, qstart=91, qend=150, sstart=61, send=120)])
    assert len(chains) == 1
    assert has_the_amplicon(chains[0], 150)


def test_a_chain_that_stops_short_of_the_amplicon_is_not_a_copy():
    """What a primer site inside an insertion looks like: the pieces are collinear and chain, but the
    chain does not reach the end of the amplicon, so there is nothing to amplify from."""
    chains = chains_of([hit(200, qstart=1, qend=80, sstart=1, send=80),
                        hit(200, qstart=81, qend=150, sstart=300, send=369)])
    assert len(chains) == 1
    assert not has_the_amplicon(chains[0], 200)


def test_two_copies_side_by_side_in_the_genome_stay_two_copies():
    """A tandem repeat of the amplicon. The second piece repeats the query rather than continuing it, so
    chaining it would turn two copies into one."""
    chains = chains_of([hit(100, sstart=1000, send=1099), hit(100, sstart=1100, send=1199)])
    assert len(chains) == 2
    assert all(has_the_amplicon(chain, 100) for chain in chains)


def test_pieces_too_far_apart_in_the_genome_are_not_one_copy():
    chains = chains_of([hit(200, qstart=1, qend=100, sstart=1, send=100),
                        hit(200, qstart=101, qend=200, sstart=50_000, send=50_099)])
    assert len(chains) == 2


def test_pieces_on_different_sequences_of_the_genome_are_not_one_copy():
    # Coordinates that would chain happily were they on one sequence: the second piece continues the
    # amplicon and sits 99 bases along. Only the name of the sequence says they are two copies.
    chains = chains_of([hit(200, qstart=1, qend=100, subject="ctg_a", sstart=1, send=100),
                        hit(200, qstart=101, qend=200, subject="ctg_b", sstart=200, send=299)])
    assert len(chains) == 2


def test_the_pieces_of_a_copy_on_the_minus_strand_chain_the_other_way_round():
    """On the minus strand the genome coordinates run down as the amplicon runs up, so a chain has to be
    read in the opposite direction -- and two pieces that ascend together are not one copy."""
    down = chains_of([hit(207, qstart=1, qend=120, sstart=1000, send=881),
                      hit(207, qstart=121, qend=207, sstart=780, send=694)])
    assert len(down) == 1
    assert has_the_amplicon(down[0], 207)
    up = chains_of([hit(207, qstart=1, qend=120, sstart=1000, send=881),
                    hit(207, qstart=121, qend=207, sstart=1100, send=1014)])
    assert len(up) == 2

def test_an_amplicon_the_genome_has_apart_from_an_insertion_is_a_copy():
    """The numbers are a real hit from the bundled example: assay8's 148-base amplicon against
    exclusion_1, which holds all of it but the planted 15-base insertion. blast calls that 87.2% identity
    because the 15 gap columns count as mismatches, so judging it on `pident` said the exclusion genome
    did not have the amplicon at all -- while it holds both primer sites exactly, and in silico PCR
    amplifies it. 129 of the 133 aligned bases match."""
    one = hit(148, gaps=15, mismatches=4)
    assert one.identity < MIN_AMPLICON_IDENTITY  # What blast reports, and what this used to be judged on
    assert has_the_amplicon([one], 148)


def test_an_amplicon_whose_hit_stops_short_of_it_is_not_a_copy():
    """The other half of that example: where a primer site is inside the insertion, blast stops the
    alignment at it instead of gapping through. Those 11 assays are specific by absence, and stay so --
    98% identity over the 93 bases it does reach, which is why identity alone cannot tell the two apart."""
    one = hit(110, qend=93, mismatches=2)
    assert one.identity > MIN_AMPLICON_IDENTITY
    assert not has_the_amplicon([one], 110)


@pytest.mark.parametrize(("reaches", "expected"), [(90, True), (89, False)])
def test_the_coverage_boundary_is_how_much_of_the_amplicon_the_hit_reaches(reaches, expected):
    assert has_the_amplicon([hit(100, qend=reaches)], 100) is expected


@pytest.mark.parametrize(("mismatches", "expected"), [(10, True), (11, False)])
def test_the_identity_boundary_counts_only_the_columns_that_align(mismatches, expected):
    assert has_the_amplicon([hit(100, mismatches=mismatches)], 100) is expected


@pytest.mark.parametrize(("mismatches", "expected"), [(9, True), (10, False)])
def test_the_gap_columns_are_left_out_of_the_identity(mismatches, expected):
    """With 10 gap columns in the alignment, 90 columns have a base in both: 9 mismatches among them is
    90% and 10 is 89%. Counting the gaps as columns -- as blast's own `pident` does -- would make both
    of them 100 columns and put the second over the line."""
    assert has_the_amplicon([hit(100, gaps=10, mismatches=mismatches)], 100) is expected


def test_a_hit_of_nothing_but_gaps_is_not_a_copy():
    """Blast does not report one, but dividing by the aligned columns would raise if it did."""
    all_gaps = Hit(query="a", qstart=1, qend=20, evalue=1e-30, length=20, gaps=20)
    assert not has_the_amplicon([all_gaps], 20)
    assert not has_the_amplicon([hit(20)], 0)
    assert not has_the_amplicon([], 20)


def test_an_exclusion_genome_that_has_the_amplicon_apart_from_an_insertion_is_counted(stubs, tmp_path,
                                                                                     fasta):
    """End to end, through blast's own output: the label has to say the genome has it."""
    stubs(amplicon_hits={"g1": 1}, amplicon_insertion={"40": 15})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGT" * 25, forward=oligo(start=0, length=20))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.exclusion_genomes_with_amplicon == 1
    assert not one.specific_by_absence


def test_the_bases_come_from_the_exclusion_alignments(stubs, tmp_path, fasta):
    """What an oligo has to mismatch is what the exclusion genomes have there, so it is read from their
    own alignments of the amplicon."""
    stubs(amplicon_hits={"g1": 1}, amplicon_mismatch={"3": "G", "7": "A"})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    # Long enough that two mismatches leave the alignment well above `MIN_AMPLICON_IDENTITY`, which is
    # what this test needs: a genome that holds the amplicon and differs from the primer inside it.
    one = assay(amplicon="ACGT" * 10, forward=oligo(start=0, length=12, variants=[3, 7]))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.forward.variant_bases == {3: "G", 7: "A"}
    assert one.forward.strong_variants == 1
    assert one.primer_variant_weight == STRONG_BASE_WEIGHT + WEAK_BASE_WEIGHT


def test_a_base_is_put_where_it_sits_in_the_region_not_in_the_amplicon(stubs, tmp_path, fasta):
    """blast reports a position in the amplicon, and the alignment need not start at the amplicon's first
    base. Both shifts have to be applied or every base lands under the wrong oligo -- or under none."""
    # The hit starts at base 3 of the amplicon, and the amplicon starts at base 40 of the region, so a
    # difference at offset 1 of the alignment is region position 40 + 2 + 1 = 43
    stubs(amplicon_hits={"g1": 1}, amplicon_qstart=3, amplicon_mismatch={"1": "G", "25": "G"})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGTACGT" * 10, forward=oligo(start=40, length=20, variants=[43]),
                reverse=oligo(start=100, length=20, reverse=True))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.forward.variant_bases == {43: "G"}
    # 40 + 2 + 25 = 67 is inside the amplicon but under neither oligo, so neither of them carries it
    assert one.reverse.variant_bases == {}


def test_a_base_outside_an_oligo_is_not_given_to_it(stubs, tmp_path, fasta):
    """A difference in the middle of the amplicon is real, but it is not something an oligo has to
    mismatch: only the bases under the oligo itself count."""
    stubs(amplicon_hits={"g1": 1}, amplicon_mismatch={"2": "T", "60": "C"})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGTACGT" * 10, forward=oligo(start=0, length=20, variants=[2]))
    check_against_exclusion([one], genomes, tmp_path / "work", threads=1)
    assert one.forward.variant_bases == {2: "T"}   # and nothing from position 60, in neither oligo


def test_the_inclusion_group_is_not_asked_for_the_aligned_sequences(stubs, tmp_path, fasta):
    """Only the exclusion genomes' bases are weighed, so counting copies in the inclusion group does not
    ask blast for the alignments or read them."""
    calls = stubs(amplicon_hits={"g1": 2})
    genomes = [fasta("g1.fasta", {"chr": "ACGT" * 50})]
    one = assay(amplicon="ACGT" * 10)
    count_inclusion_copies([one], genomes, tmp_path / "incl", threads=1)
    assert one.inclusion_copies == [2]
    assert all(not (asked_fields(line) & {"qseq", "sseq"}) for line in blast_calls(calls))

    check_against_exclusion([one], genomes, tmp_path / "excl", threads=1)
    assert any({"qseq", "sseq"} <= asked_fields(line) for line in blast_calls(calls))  # the exclusion group is


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


def test_designing_on_a_last_resort_region_says_so(stubs, finished_run, tmp_path, caplog):
    """Someone running `design` days later on a results folder should not have to remember how the
    regions were found. The find command marks them; this repeats it."""
    final = finished_run / "final_kmers.fasta"
    final.write_text(final.read_text().replace(">ctg1 [100, 101, 102]",
                                               ">ctg1 [100, 101, 102] last_resort"))
    assert run(design_settings(finished_run, tmp_path / "assays")) == 0
    assert "last_resort" in caplog.text
    assert "rests on one mismatch" in caplog.text


def test_designing_on_ordinary_regions_says_nothing_about_a_last_resort(stubs, finished_run, tmp_path,
                                                                        caplog):
    assert run(design_settings(finished_run, tmp_path / "assays")) == 0
    assert "last_resort" not in caplog.text


def test_an_output_folder_inside_an_input_folder_is_refused(stubs, finished_run, genomes, tmp_path):
    """Everything the command writes under its output folder -- the amplicon databases, a copy of each
    genome, the staged fasta files -- is read back as an input genome if the output sits inside an input
    folder. One run counted 16 inclusion genomes where there were 4, which changes every in-silico
    denominator. `find` refuses this; so does this.
    """
    inclusion, _ = genomes
    with pytest.raises(PrimerFinderError, match="output folder is inside the inclusion folder"):
        run(design_settings(finished_run, inclusion / "assays"))


def test_an_input_folder_inside_the_output_folder_is_refused(stubs, finished_run, tmp_path, genomes):
    """The staged folders are emptied before they are filled, so pointing -i at something under the
    output folder deletes those files. Reached by pointing a second run at the first run's staged folder,
    which the documentation describes as a flat folder of plain fasta -- and that run wiped it, including
    the one genome that was a real file rather than a symlink."""
    inclusion, exclusion = genomes
    staged = tmp_path / "assays" / "insilico_genomes" / "inclusion"
    staged.mkdir(parents=True)
    (staged / "g.fasta").write_text(">chr\nACGT\n")
    with pytest.raises(PrimerFinderError, match="inclusion folder is inside the output folder"):
        run(design_settings(finished_run, tmp_path / "assays", inclusion=staged, exclusion=exclusion))
    assert (staged / "g.fasta").is_file(), "the refusal has to come before anything is deleted"


def test_an_output_path_with_a_space_is_refused(stubs, finished_run, tmp_path):
    """blast reopens a database by absolute path and cannot quote whitespace. The design command builds
    databases under its output folder, so it refuses the path up front rather than failing part way in."""
    with pytest.raises(PrimerFinderError, match="contains a space"):
        run(design_settings(finished_run, tmp_path / "with space"))


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


def one_difference_run(tmp_path, genomes, stubs) -> Path:
    """A finished run whose region offers two differences too far apart for one oligo, so every assay in
    it can only ever spend one of them. The stub's left primer covers bases 10-29."""
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    region = "A" * 15 + "c" + "A" * 80 + "t" + "A" * 100
    (results / "final_kmers.fasta").write_text(f">ctg1 [15, 96]\n{region}\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    stubs(amplicon_hits={path.name.split(".")[0]: 1 for path in exclusion.glob("*.fasta")})
    return results


def test_an_assay_resting_on_one_difference_is_set_aside(stubs, tmp_path, genomes):
    """A region is kept because two differences could sit in one oligo. An assay inside it that spends
    only one of them is the weak case that rule was there to avoid, so it is not carried further -- and
    when that is every assay, the error says which option accepts them."""
    results = one_difference_run(tmp_path, genomes, stubs)
    with pytest.raises(PrimerFinderError, match="min-primer-differences 1"):
        run(design_settings(results, tmp_path / "strict"))


def test_the_weaker_assay_can_be_asked_for(stubs, tmp_path, genomes):
    results = one_difference_run(tmp_path, genomes, stubs)
    output = tmp_path / "loose"
    assert run(design_settings(results, output, min_primer_differences=1)) == 0
    rows = list(iter_rows(output / "assays.tsv"))
    assert rows and all(row["best_primer_variants"] == "1" for row in rows)
    assert rows[0]["assay"] in (output / "assays_pcr.fasta").read_text()


def probe_only_run(tmp_path, genomes, stubs) -> Path:
    """A finished run whose one difference falls under the stub's probe (35-54) and outside its plain
    left primer (10-29), so the plain request yields assays no primer can discriminate with."""
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    region = "A" * 40 + "c" + "A" * 159
    (results / "final_kmers.fasta").write_text(f">ctg1 [40]\n{region}\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    stubs(amplicon_hits={path.name.split(".")[0]: 1 for path in exclusion.glob("*.fasta")})
    return results


def test_an_assay_with_no_primer_difference_is_set_aside_at_the_lowest_setting(stubs, tmp_path, genomes):
    """`--min-primer-differences 1` means one difference under a *primer*, not one under any oligo: an
    assay whose only differences are under its probe is set aside however low the setting goes."""
    results = probe_only_run(tmp_path, genomes, stubs)
    output = tmp_path / "loose"
    assert run(design_settings(results, output, min_primer_differences=1)) == 0
    rows = list(iter_rows(output / "assays.tsv"))
    probe_only = [row for row in rows if row["best_primer_variants"] == "0"]
    assert probe_only, "the scenario is pointless without an assay no primer can discriminate with"
    assert any(row["best_primer_variants"] != "0" for row in rows), "nor without one that a primer can"
    carried = (output / "assays_pcr.fasta").read_text()
    for row in probe_only:
        assert row["assay"] not in carried, row["assay"]


def test_the_advice_when_everything_is_set_aside_fits_the_setting(stubs, tmp_path, genomes):
    """At 2 there is a looser setting worth suggesting; at 1 the only way back is 0, and suggesting 1
    again would be no help at all."""
    results = probe_only_run(tmp_path, genomes, stubs)
    with pytest.raises(PrimerFinderError, match="min-primer-differences 1 accepts"):
        run(design_settings(results, tmp_path / "two", min_primer_differences=2))

    # 80 bases: there is room for the probe pin but not for either primer pin, so no primer can be put
    # over the difference at all and every assay is probe-only
    inclusion, exclusion = genomes
    narrow = tmp_path / "narrow"
    narrow.mkdir()
    (narrow / "final_kmers.fasta").write_text(">ctg1 [40]\n" + "A" * 40 + "c" + "A" * 39 + "\n")
    (narrow / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    stubs(amplicon_hits={path.name.split(".")[0]: 1 for path in exclusion.glob("*.fasta")})
    with pytest.raises(PrimerFinderError, match="min-primer-differences 0 turns this off"):
        run(design_settings(narrow, tmp_path / "one", min_primer_differences=1))


def test_turning_the_filter_off_accepts_them(stubs, tmp_path, genomes):
    """0 is the only way back to the older, looser rule, and the error says so when everything is set
    aside."""
    results = probe_only_run(tmp_path, genomes, stubs)
    output = tmp_path / "off"
    assert run(design_settings(results, output, min_primer_differences=0)) == 0
    carried = (output / "assays_pcr.fasta").read_text()
    probe_only = [row for row in iter_rows(output / "assays.tsv")
                  if row["best_primer_variants"] == "0"]
    assert probe_only and all(row["assay"] in carried for row in probe_only)


def test_a_region_too_short_for_the_product_range_is_skipped_not_fatal(stubs, tmp_path, genomes, caplog):
    """Primer3 answers a request whose template is shorter than the smallest product with an error, not
    with an empty result, and the command treats a Primer3 error as fatal. One short region would
    therefore end a run over hundreds of good ones -- and short regions became common when `find` started
    reporting them down to one kmer."""
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text(
        ">tiny [40]\n" + "A" * 40 + "c" + "A" * 59 + "\n"          # 100 bases, under the 120 floor
        ">roomy [40]\n" + "A" * 40 + "c" + "A" * 259 + "\n")       # 300 bases, over it
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    output = tmp_path / "assays"
    with caplog.at_level("INFO"):
        assert run(design_settings(results, output, product_size="120-200")) == 0
    assert "1 of 2 region(s) are shorter than the 120-base smallest product" in caplog.text
    assert {row["region"] for row in iter_rows(output / "assays.tsv")} == {"roomy"}


def test_the_regions_with_no_room_are_left_out_before_max_regions(stubs, tmp_path, genomes, caplog):
    """--max-regions says how many regions to design on. Spending one of those places on a region that
    cannot hold a product would leave out one that can, so the length filter comes first."""
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    # the short ones first: at --max-regions 2 the old order spent both places on them and designed on none
    regions = [("tiny1", 100), ("tiny2", 100), ("roomy1", 300), ("roomy2", 300), ("roomy3", 300)]
    results.joinpath("final_kmers.fasta").write_text("".join(
        f">{name} [40]\n{'A' * 40}c{'A' * (length - 41)}\n" for name, length in regions))
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    output = tmp_path / "assays"
    with caplog.at_level("INFO"):
        assert run(design_settings(results, output, product_size="120-200", max_regions=2)) == 0
    assert "2 of 5 region(s) are shorter than the 120-base smallest product" in caplog.text
    assert "the first 2 of 3 regions with room in them" in caplog.text
    designed = {row["region"] for row in iter_rows(output / "assays.tsv")}
    assert designed == {"roomy1", "roomy2"}
    counts = json.loads((output / "design_info.json").read_text())["counts"]
    assert counts["regions"] == 2    # what it designed on, not what it was handed


def test_every_region_too_short_says_what_to_change(stubs, tmp_path, genomes):
    inclusion, exclusion = genomes
    results = tmp_path / "results"
    results.mkdir()
    (results / "final_kmers.fasta").write_text(">tiny [40]\n" + "A" * 40 + "c" + "A" * 59 + "\n")
    (results / "run_info.json").write_text(json.dumps(
        {"parameters": {"inclusion": str(inclusion), "exclusion": str(exclusion)}}))
    with pytest.raises(PrimerFinderError, match="Every region is shorter than"):
        run(design_settings(results, tmp_path / "assays", product_size="120-200"))


def test_a_product_size_range_primer3_cannot_read_is_refused():
    """The design step reads the smallest product before Primer3 ever runs, to leave out the regions with
    no room for one. A range it cannot read would give a wrong reason for leaving them out."""
    assert check_product_size("70-150") == "70-150"
    assert check_product_size("70-150 200-300") == "70-150 200-300"
    assert "runs backwards" in check_product_size("150-70")
    assert "two whole numbers" in check_product_size("abc")
    assert "two whole numbers" in check_product_size("70_150")
    assert "below 1" in check_product_size("0-100")
    assert "empty" in check_product_size("   ")
    assert "two whole numbers" in check_product_size("70-150 oops")


def test_the_default_asks_for_two_differences_under_one_primer():
    assert DEFAULT_MIN_PRIMER_DIFFERENCES == 2
    assert DesignSettings(results=Path("r"), output=Path("o"), threads=1).min_primer_differences == 2


def test_a_probes_differences_do_not_qualify_an_assay():
    """The probe is longer and hotter than a primer, so two mismatches under it are no evidence it will
    not bind: a conventional TaqMan probe still gave signal through five (Yao et al. 2006). Its
    differences are reported and break ties in the ranking, but they cannot carry an assay."""
    probe_only = assay(probe=oligo(140, 22, variants=[150, 151]))
    assert probe_only.probe_variants_covered == 2
    assert probe_only.best_primer_variants == 0

    on_a_primer = assay(forward=oligo(start=100, length=20, variants=[118, 119]))
    assert on_a_primer.best_primer_variants == 2


def test_one_difference_on_each_primer_does_not_qualify_either():
    """Two differences, but neither primer carries two: the region was kept because two could sit in one
    oligo, and this assay does not do that."""
    split = assay(forward=oligo(start=100, length=20, variants=[119]),
                  reverse=oligo(200, variants=[200], reverse=True))
    assert split.primer_variants_covered == 2
    assert split.best_primer_variants == 1


def test_an_assay_specific_by_absence_is_not_asked_for_differences(stubs, finished_run, tmp_path, genomes):
    """It does not rest on a difference at all: there is nothing in the exclusion genomes to amplify, so
    the number of differences under its primers is beside the point."""
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, min_primer_differences=4)) == 0
    rows = list(iter_rows(output / "assays.tsv"))
    kept = [row for row in rows if row["assay"] in (output / "assays_pcr.fasta").read_text()]
    assert kept and all(row["specific_by"] == "absence" for row in kept)


def test_how_many_differences_the_best_primer_covers_is_reported(tmp_path):
    one = assay(forward=oligo(start=100, length=20, variants=[115, 119]),
                probe=oligo(start=140, length=22, variants=[150, 151, 152]))
    path = tmp_path / "assays.tsv"
    write_assays(path, [one])
    row = next(iter_rows(path))
    assert row["best_primer_variants"] == "2"   # the forward primer's two, not the probe's three
    assert row["probe_variants"] == "3"         # which is still reported, just not what qualifies


def test_the_structure_temperatures_survive_a_whole_run(stubs, finished_run, tmp_path, genomes):
    """Primer3 reports these whenever it does its thermodynamic alignment, so a real run fills them. They
    are all zero if anything between the output and the table drops them."""
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output)) == 0
    row = next(iter_rows(output / "assays.tsv"))
    assert row["forward_hairpin_tm"] == "35.7" and row["forward_self_dimer_tm"] == "13.6"
    assert row["reverse_hairpin_tm"] == "24.3" and row["reverse_self_dimer_tm"] == "11.2"
    assert row["probe_hairpin_tm"] == "33.3" and row["probe_self_dimer_tm"] == "19.6"
    assert row["pair_dimer_tm"] == "8.7" and row["pair_dimer_end_tm"] == "4.1"


def test_a_primer3_that_reports_no_structures_still_runs(stubs, finished_run, tmp_path, genomes):
    """Those fields need Primer3's thermodynamic alignment. Without them an assay still has primers, and
    the columns say zero rather than stopping the run."""
    stubs(primer3_structures=False)
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output)) == 0
    row = next(iter_rows(output / "assays.tsv"))
    assert row["forward_hairpin_tm"] == "0.0" and row["pair_dimer_tm"] == "0.0"


def test_an_exclusion_amplification_that_rests_on_a_trimmed_primer_end_is_reported(stubs, finished_run,
                                                                                   tmp_path, genomes,
                                                                                   insilico_pcr):
    """insilicoPCR said it amplified an exclusion genome, so it is not selective -- but the amplification
    rests on bases insilicoPCR does not count, which it could not have refused. That gets its own label
    rather than being lumped in with a cross-reaction."""
    inclusion, exclusion = genomes
    sample = sorted(path.name.split(".")[0] for path in exclusion.glob("*.fasta"))[0]
    stubs(insilico_extra={"ctg1_assay0": [sample]}, insilico_terminal={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    rows = list(iter_rows(output / "assays.tsv"))
    row = next(r for r in rows if r["assay"] == "ctg1_assay0")
    for kind in ("pcr", "qpcr"):
        assert row[f"{kind}_selective"] == "undecided (3' end)", kind
        assert row[f"{kind}_exclusion_amplified"] == "1", kind
        assert row[f"{kind}_exclusion_terminal_only"] == "1", kind
    # below the assays nothing was found against, and above any the check refused on its own evidence
    assert rows[0]["assay"] != "ctg1_assay0"
    assert rows[0]["pcr_selective"] == "yes"


def test_an_exclusion_amplification_through_a_primers_five_prime_end_is_a_cross_reaction(
        stubs, finished_run, tmp_path, genomes, insilico_pcr):
    """The same consolidated report, with the trim at the other end of the primer. A 5'-end difference
    does not stop the polymerase extending, so that amplification is real and gets the plain `no` -- not
    the benefit of the doubt that the 3'-end blind spot earns."""
    inclusion, exclusion = genomes
    sample = sorted(path.name.split(".")[0] for path in exclusion.glob("*.fasta"))[0]
    stubs(insilico_extra={"ctg1_assay0": [sample]}, insilico_five_prime={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    row = next(r for r in iter_rows(output / "assays.tsv") if r["assay"] == "ctg1_assay0")
    for kind in ("pcr", "qpcr"):
        assert row[f"{kind}_selective"] == "no", kind
        assert row[f"{kind}_exclusion_amplified"] == "1", kind
        assert row[f"{kind}_exclusion_terminal_only"] == "0", kind


def test_a_real_exclusion_amplification_is_not_counted_as_terminal(stubs, finished_run, tmp_path, genomes,
                                                                   insilico_pcr):
    inclusion, exclusion = genomes
    sample = sorted(path.name.split(".")[0] for path in exclusion.glob("*.fasta"))[0]
    stubs(insilico_extra={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    assert run(design_settings(finished_run, output, insilico_pcr=insilico_pcr)) == 0
    row = next(r for r in iter_rows(output / "assays.tsv") if r["assay"] == "ctg1_assay0")
    assert row["pcr_exclusion_amplified"] == "1" and row["pcr_exclusion_terminal_only"] == "0"


def test_an_assay_that_misses_an_inclusion_genome_is_not_selective(stubs, finished_run, tmp_path,
                                                                   genomes, insilico_pcr):
    inclusion, _ = genomes
    sample = sorted(path.name.split(".")[0] for path in inclusion.glob("*.fasta"))[0]
    stubs(insilico_misses={"ctg1_assay0": [sample]})
    output = tmp_path / "assays"
    run(design_settings(finished_run, output, insilico_pcr=insilico_pcr))
    bad = next(row for row in iter_rows(output / "assays.tsv") if row["assay"] == "ctg1_assay0")
    assert bad["qpcr_inclusion_amplified"] == "1" and bad["qpcr_selective"] == "no"


def blast_calls(calls: Path) -> list[str]:
    lines = [line for line in calls.read_text().splitlines() if line.startswith("blastn")]
    assert lines, "no blastn call was made"
    return lines


def asked_fields(call: str) -> set[str]:
    """The fields of a `-outfmt "6 ..."` argument, as a set: looking for a substring would make `qseqid`
    answer for `qseq`, and would depend on the order they happen to be written in."""
    import shlex

    words = shlex.split(call)
    outfmt = words[words.index("-outfmt") + 1].split()
    assert outfmt and outfmt[0] == "6", call
    return set(outfmt[1:])


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


def test_the_weighting_columns_are_written_from_their_own_numbers(tmp_path):
    """Three columns of the same row come from three different properties, and all three are small
    integers: a column fed from the wrong one would look plausible."""
    one = assay(forward=oligo(start=100, length=20, variants=[105, 117, 118, 119]))
    one.forward.variant_bases = {105: "G", 117: "A", 118: "A", 119: "A"}
    path = tmp_path / "assays.tsv"
    write_assays(path, [one])
    row = next(iter_rows(path))
    assert row["forward_variants"] == "4"
    assert row["forward_strong_variants"] == "1"        # only the G, and it is the one furthest from the end
    assert row["forward_near_3prime"] == "3"            # 117, 118 and 119, within five bases of the end
    assert row["forward_terminal_run"] == "3"           # and all three in a row, ending at it
    assert row["primer_variant_weight"] == "2.5"        # 1.0 for the G and 0.5 for each A


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
