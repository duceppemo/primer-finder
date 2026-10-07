"""Cigar strings, and which contigs could carry a selective assay."""

from __future__ import annotations

import pytest

from primer_finder.mapping import (
    Alignment,
    cigar_string,
    is_candidate,
    mark_differences,
    parse_cigar,
    parse_sam,
    reverse_complement,
    select_candidates,
    sort_candidates,
)
from primer_finder.seqio import Record


def test_parse_and_write_cigar():
    assert parse_cigar("10=1X5I3D2S") == [(10, "="), (1, "X"), (5, "I"), (3, "D"), (2, "S")]
    assert cigar_string(parse_cigar("10=1X")) == "10=1X"
    assert parse_cigar("*") == []


def test_reverse_complement_keeps_the_case():
    assert reverse_complement("AAccGGtN") == "NaCCggTT"


@pytest.mark.parametrize(
    ("cigar", "expected"),
    [
        ("100=", False),                      # Identical to the exclusion genome
        ("50=1X49=", False),                  # One mismatch: a primer cannot be selective on its own
        ("50=1X21=1X27=", False),             # 21 matching bases between them: too far apart
        ("50=1X20=1X28=", True),              # 20 matching bases between them: one primer covers both
        ("50=2X48=", True),                   # Two mismatches side by side
        ("50=3I50=", True),                   # An insertion
        ("50=3D50=", True),                   # A deletion
        ("50=1I50=", False),                  # A single inserted base
        ("50=1X5=1I44=", True),               # A mismatch and an insertion, close together
        ("1X5=1X", True),                     # Two mismatches at the very end of the cigar
        ("150=1I1X149=", True),               # An inserted base and a mismatch that touch
        ("50=1I1X48=", True),
        ("100=1X1D1X197=", True),             # Two mismatches on either side of a deleted base
        ("250S150=1X149=", True),             # 250 bases the exclusion genome does not hold at all
        ("5S150=1X149=", False),              # A clipped end too short for a primer
        ("150=1X149=250S", True),
    ],
)
def test_is_candidate(cigar, expected):
    assert is_candidate(parse_cigar(cigar)) is expected


def test_mark_differences_marks_mismatches_and_insertions_only():
    # 4 matches, 2 mismatches, 3 deleted bases (not in the contig), 2 matches, 2 inserted bases
    assert mark_differences("AAAACCGGTT", parse_cigar("4=2X3D2=2I")) == "AAAAccGGtt"


def test_mark_differences_keeps_what_the_cigar_does_not_cover():
    assert mark_differences("AAAACCCC", parse_cigar("4=")) == "AAAACCCC"


def test_mark_differences_marks_clipped_ends():
    """The aligner could not place those bases in the exclusion genome: they are not matches."""
    assert mark_differences("AAAACCCC", parse_cigar("2S4=2S")) == "aaAACCcc"


def test_parse_sam_skips_headers_and_secondary_alignments(tmp_path):
    sam = tmp_path / "x.sam"
    sam.write_text(
        "@HD\tVN:1.6\n"
        "ctg1\t0\tref\t1\t60\t10=\t*\t0\t0\tACGTACGTAC\t*\n"
        "ctg1\t256\tref\t1\t60\t10=\t*\t0\t0\tACGTACGTAC\t*\n"   # secondary
        "ctg2\t2064\tref\t1\t60\t10=\t*\t0\t0\tACGTACGTAC\t*\n"  # supplementary
        "ctg3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
        "truncated\t0\tref\n"
    )
    alignments = list(parse_sam(sam))
    assert [(a.name, a.mapped, a.reverse) for a in alignments] == [
        ("ctg1", True, False), ("ctg3", False, False)
    ]
    assert alignments[0].cigar == [(10, "=")]


def test_an_unmapped_contig_is_entirely_specific():
    assembly = {"ctg1": Record("ctg1", "", "ACGTACGT")}
    alignment = Alignment("ctg1", mapped=False, reverse=False, cigar=[])
    candidates = select_candidates(assembly, [alignment])
    assert candidates["ctg1"].desc == "8I"
    assert candidates["ctg1"].seq == "acgtacgt"


def test_a_reverse_alignment_marks_the_bases_of_the_forward_contig():
    # The contig maps to the minus strand: the cigar describes its reverse complement, whose first base
    # is a mismatch. That base is the last one of the contig as it was assembled.
    assembly = {"ctg1": Record("ctg1", "", "AAAACCCCGG")}
    alignment = Alignment("ctg1", mapped=True, reverse=True, cigar=[(2, "X"), (8, "=")])
    candidate = select_candidates(assembly, [alignment])["ctg1"]
    assert candidate.seq == "AAAACCCCgg"
    assert candidate.desc == "2X8="


def test_a_contig_the_assembly_does_not_hold_is_ignored(caplog):
    alignment = Alignment("ghost", mapped=True, reverse=False, cigar=[(2, "X")])
    assert select_candidates({}, [alignment]) == {}
    assert "not in the assembly" in caplog.text


def test_contigs_with_the_most_differing_bases_come_first():
    from primer_finder.mapping import Candidate

    candidates = {
        "two_marks": Candidate("two_marks", "10=1X5=1X10=", "AAAAAAAAAAaAAAAAaAAAAA"),
        "one_mark_b": Candidate("one_mark_b", "10=1X10=", "AAAAAAAAAAaAAAAAAAAAA"),
        "one_mark_a": Candidate("one_mark_a", "10=1X10=", "AAAAAAAAAAaAAAAAAAAAA"),
        "all_specific": Candidate("all_specific", "412I", "a" * 412),
    }
    assert [candidate.name for candidate in sort_candidates(candidates)] == [
        "all_specific", "two_marks", "one_mark_a", "one_mark_b"
    ]


def test_a_mapped_record_without_a_cigar_is_treated_as_unmapped():
    """A SAM record can carry a position and no cigar string at all."""
    assembly = {"ctg1": Record("ctg1", "", "ACGTACGT")}
    alignment = Alignment("ctg1", mapped=True, reverse=False, cigar=[])
    assert select_candidates(assembly, [alignment])["ctg1"].desc == "8I"


def test_reverse_complement_of_the_ambiguity_codes():
    assert reverse_complement("RYKMBDHVN") == "NBDHVKMRY"  # B<->V and D<->H
    assert reverse_complement(reverse_complement("ACGTRYKMBDHVNacgt")) == "ACGTRYKMBDHVNacgt"
