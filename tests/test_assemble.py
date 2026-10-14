"""Assembling the inclusion-specific kmers: the unitig walk, and what each assembler is given."""

from __future__ import annotations

import random

import pytest

from primer_finder import PrimerFinderError
from primer_finder.assemble import SPADES_FLOOR, assemble, programs_for, unitigs
from primer_finder.seqio import Record, iter_records, reverse_complement, write_fasta

SEQUENCE = "ACGTTGCAAGGCTTACCGGATTCAGGCTAACTTGCCAATGCCGTTAAGCCTTAGGCAATTGGCCAATTGG"


def kmers_of(sequence: str, k: int = 21) -> list[str]:
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def test_every_kmer_of_one_stretch_rebuilds_that_stretch():
    built = unitigs(kmers_of(SEQUENCE))
    assert len(built) == 1
    # A kmer set does not say which strand it came from, so either orientation is the right answer
    assert built[0] in (SEQUENCE, reverse_complement(SEQUENCE))


def test_the_answer_does_not_depend_on_the_order_or_the_strand_of_the_input():
    kmers = kmers_of(SEQUENCE)
    shuffled = kmers[:]
    random.Random(0).shuffle(shuffled)
    flipped = [reverse_complement(kmer) for kmer in kmers]
    assert unitigs(shuffled) == unitigs(kmers) == unitigs(flipped)


def test_a_branch_ends_a_unitig():
    """Where the sequence could carry on two ways the stretch stops: a unitig is what the kmers
    determine, never a choice between alternatives."""
    kmers = kmers_of(SEQUENCE)
    fork = SEQUENCE[20:40] + "A"  # a 21-mer that leaves the stretch at base 40
    built = unitigs([*kmers, fork])
    assert len(built) > 1
    assert all(len(one) < len(SEQUENCE) for one in built)


OTHER = "TTGACCTAGCATCGGATACGTTTCAGGCATAGCCTTGAACCGTATCAGGTT"


def test_two_stretches_running_into_one_keep_the_shared_part_separate():
    """Two different sequences that end in the same 30 bases. The kmers do not say which of the two the
    shared part belongs to, so it comes back as a stretch of its own rather than being tacked onto one of
    them."""
    left_one = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACT"
    left_two = "AAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCAC"
    shared = "TGTTGGCCCAGTGTGAATCGCTTAAGGGTT"
    built = unitigs(kmers_of(left_one + shared) + kmers_of(left_two + shared))
    assert [len(one) for one in built] == [60, 60, 30]
    assert sum(1 for one in built if shared in one or reverse_complement(shared) in one) == 1


def test_backing_up_stops_where_the_sequence_could_have_gone_two_ways():
    """Walking back to the start of a stretch has to stop not only where two sequences run together, but
    also one step earlier than that -- where the step back could itself have gone on two ways. Without
    that, a stretch comes back in pieces: this pair of sequences gives 32, 29 and 26 bases, and gives
    32, 29, 21 and 15 if the walk backs up one node too far. Found by searching random graphs rather than
    by reasoning, which is why it is pinned here with its own sequences."""
    one = "CATTGGTGCATAGAGCCTGTACTAGCTTAATGGCATAATCGGTCCACCGGA"
    two = "GCGTTAACGCCCTTTATACTAGCTTAATGGCATAATCGGTCCACCGGA"
    built = unitigs(kmers_of(one, k=11) + kmers_of(two, k=11))
    assert [len(piece) for piece in built] == [32, 29, 26]


def test_the_stretches_come_back_longest_first():
    """`find` reports regions in the order it gets them, so a stable order here is what makes a run
    repeatable. Equal lengths are broken by sequence, which nothing else would settle."""
    built = unitigs(kmers_of(SEQUENCE) + kmers_of(OTHER))
    assert [len(one) for one in built] == sorted((len(one) for one in built), reverse=True)
    same_length = unitigs(kmers_of("ACGTTGCAAGGCTTACCGGATTCAGGCTAAC")
                          + kmers_of("TTTTCCCCGGGGAAATTTGGGCCCAAATTTG"))
    assert [len(one) for one in same_length] == [31, 31]
    assert same_length == sorted(same_length)  # the tie-break, which is the sequence itself


def test_two_unrelated_stretches_stay_two():
    built = unitigs(kmers_of(SEQUENCE) + kmers_of(OTHER))
    assert sorted(len(one) for one in built) == sorted([len(SEQUENCE), len(OTHER)])


def test_a_repeat_ends_the_stretches_around_it():
    """Not a defect. Where the same 30 bases occur twice with different sequence after each copy, the
    kmers no longer say which way the sequence carried on, and a unitig that picked one would be an
    invention. It is the reason to prefer a longer -k, and the reason `-d` exists for repeated regions.

    A, B and C are distinct 40-base stretches and R a 30-base repeat, laid out as A-R-B-R-C."""
    a = "GCTAAAGACAATTACATAACATACACGTCAGCACGAAACT"
    repeat = "TGTTGGCCCAGTGTGAATCGCTTAAGGGTT"
    b = "AAGTAAGTGTGATGCATACGCCTTTACTTGCTGTGTCCAC"
    c = "CCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAAC"
    built = unitigs(kmers_of(a + repeat + b + repeat + c))
    # The repeat comes back once, as a stretch of its own, and the three unique parts keep the 20 bases
    # of overlap they share with it on either side
    assert sorted(len(one) for one in built) == [30, 60, 60, 80]
    assert any(one in (repeat, reverse_complement(repeat)) for one in built)


def test_a_circular_stretch_comes_back_once_instead_of_going_round_for_ever():
    """Every kmer of a circle has one way in and one way out, so the walk would never meet a reason to
    stop. A plasmid, or a tandem array, is such a circle. The stretch is reported once, the length of the
    circle plus the overlap a linear reading of it needs."""
    circle = SEQUENCE[:45]
    k = 21
    wrapped = circle + circle[:k - 1]  # the kmer set of the circle, origin included
    built = unitigs([wrapped[i:i + k] for i in range(len(wrapped) - k + 1)])
    assert len(built) == 1
    assert len(built[0]) == len(circle) + k - 1


def test_a_single_kmer_is_a_unitig_of_its_own_length():
    assert unitigs([SEQUENCE[:21]]) in ([SEQUENCE[:21]], [reverse_complement(SEQUENCE[:21])])


def test_no_kmers_at_all():
    assert unitigs([]) == []


def test_fifteen_kmers_are_enough():
    """The case this was written for: SKESA refuses a handful of kmers, which is exactly the situation
    where the signal is a single short region and the run has nothing else to report."""
    fifteen = kmers_of(SEQUENCE)[:15]
    built = unitigs(fifteen)
    assert len(built) == 1
    assert len(built[0]) == 21 + 14  # fifteen kmers each one base further on


def test_kmers_of_different_lengths_are_refused():
    with pytest.raises(PrimerFinderError, match="not all the same length"):
        unitigs([SEQUENCE[:21], SEQUENCE[:31]])


def test_lower_case_kmers_are_read_the_same(tmp_path):
    assert unitigs([one.lower() for one in kmers_of(SEQUENCE)]) == unitigs(kmers_of(SEQUENCE))


def test_the_unitig_assembler_needs_no_external_program():
    assert programs_for("unitig") == set()
    assert programs_for("skesa") == {"skesa"}
    assert programs_for("spades") == {"spades.py"}


def test_assembling_with_the_unitig_walk_writes_the_contigs(tmp_path):
    kmer_fasta = tmp_path / "kmers.fasta"
    write_fasta(kmer_fasta, [Record(f"kmer_{n}", "", one)
                             for n, one in enumerate(kmers_of(SEQUENCE))])
    out = assemble(kmer_fasta, tmp_path / "assembly.fasta", "unitig", threads=1, memory_gb=1)
    written = list(iter_records(out))
    assert len(written) == 1
    assert len(written[0].seq) == len(SEQUENCE)
    assert written[0].name == "unitig_1"


def test_the_unitig_assembler_honours_a_floor_below_what_spades_can_do(tmp_path):
    """SPAdes cannot report a contig shorter than about 200 bases whatever it is asked; this can, down to
    a single kmer, which is the whole point of `--min-region`."""
    short = SEQUENCE[:40]
    kmer_fasta = tmp_path / "kmers.fasta"
    write_fasta(kmer_fasta, [Record(f"kmer_{n}", "", one) for n, one in enumerate(kmers_of(short))])
    assert len(short) < SPADES_FLOOR
    out = assemble(kmer_fasta, tmp_path / "assembly.fasta", "unitig", threads=1, memory_gb=1,
                   min_contig=len(short))
    assert len(next(iter(iter_records(out))).seq) == len(short)


def test_a_floor_above_every_unitig_is_an_error_not_an_empty_file(tmp_path):
    kmer_fasta = tmp_path / "kmers.fasta"
    write_fasta(kmer_fasta, [Record(f"kmer_{n}", "", one)
                             for n, one in enumerate(kmers_of(SEQUENCE))])
    with pytest.raises(PrimerFinderError, match="could not assemble"):
        assemble(kmer_fasta, tmp_path / "assembly.fasta", "unitig", threads=1, memory_gb=1,
                 min_contig=len(SEQUENCE) + 1)


def test_the_error_from_an_assembler_for_reads_points_at_the_unitig_walk(stubs, tmp_path):
    """The old advice was to lower -k, which runs into the same assembler from the other end."""
    kmer_fasta = tmp_path / "kmers.fasta"
    write_fasta(kmer_fasta, [Record("kmer_0", "", SEQUENCE[:21])])
    stubs(contigs={})
    with pytest.raises(PrimerFinderError, match="-a unitig"):
        assemble(kmer_fasta, tmp_path / "assembly.fasta", "skesa", threads=1, memory_gb=1)

def every_kmer_once(kmers: list[str], built: list[str]) -> None:
    """The two invariants that define a correct decomposition, which no single example pins down:

    - every kmer given appears in exactly one stretch, in one orientation or the other;
    - every k-base window of every stretch is a kmer that was given.

    The first says nothing was lost or duplicated, the second that nothing was invented.
    """
    canonical = {min(one, reverse_complement(one)) for one in (kmer.upper() for kmer in kmers)}
    k = len(next(iter(canonical)))
    for key in canonical:
        holding = [one for one in built if key in one or reverse_complement(key) in one]
        assert len(holding) == 1, f"{key} is in {len(holding)} stretch(es), not 1"
    for one in built:
        for at in range(len(one) - k + 1):
            window = min(one[at:at + k], reverse_complement(one[at:at + k]))
            assert window in canonical, f"{one} spells out {window}, which was never given"


def test_a_hairpin_does_not_lose_kmers():
    """A kmer whose last k-1 bases are their own reverse complement leads back into itself the other way
    round. Backing up through it used to leave the seed in no stretch at all: of these four kmers -- two,
    once reverse complements are collapsed -- the walk returned one, and the other was simply gone.
    `ACGTAT` is the palindrome."""
    kmers = ["ACGTATC", "TACGTAT", "ATACGTA", "GATACGT"]
    built = unitigs(kmers)
    every_kmer_once(kmers, built)
    assert built == ["TACGTATC"]


def test_a_palindromic_kmer_is_its_own_stretch():
    """With an even k a kmer can be its own reverse complement. It has to come back, once."""
    kmers = ["ACGT", "CGTA", "GTAC"]
    built = unitigs(kmers)
    every_kmer_once(kmers, built)


# Twelve arbitrary seeds for general cover, and three found by searching 400 of them for a graph the
# pre-fix walk got wrong: 203, 245 and 309 each lose a kmer without the fix, and 19, 162 and 249 change
# the answer without losing one. Pinned because the shape needs a hairpin *and* a tip in the right place,
# which an arbitrary seed produces about once in seventy.
@pytest.mark.parametrize("seed", [*range(12), 19, 162, 203, 245, 249, 309])
def test_the_invariants_hold_on_random_graphs_including_palindromes(seed):
    """Branches, tips, repeats and hairpins in combination, which no hand-written example covers. The
    palindrome is spliced in on purpose: that is the case the walk got wrong."""
    rng = random.Random(seed)
    k = rng.choice([5, 7, 9, 11, 21])

    def chunk(length: int) -> str:
        return "".join(rng.choice("ACGT") for _ in range(length))

    # Built to guarantee the shape that broke the walk rather than hope for it: a palindrome of exactly
    # k - 1 bases means the kmer starting one base before it ends on a palindrome, and backing up through
    # such a kmer leads into itself the other way round. An odd length cannot be a palindrome -- the
    # middle base would have to pair with itself -- so with an even k the palindrome is k long instead,
    # which makes a kmer that is its own reverse complement.
    half = (k - 1) // 2 if (k - 1) % 2 == 0 else k // 2
    palindrome = chunk(half)
    palindrome += reverse_complement(palindrome)
    assert palindrome == reverse_complement(palindrome)
    pieces = [chunk(rng.randint(k, 3 * k)) for _ in range(3)]
    sequences = [
        pieces[0] + palindrome + pieces[1],
        pieces[2] + palindrome + pieces[0],          # the palindrome reached from two directions
        pieces[1] + pieces[2],                        # and a join that makes a branch
    ]
    kmers = [one[at:at + k] for one in sequences for at in range(len(one) - k + 1)]
    kmers = [one for one in kmers if rng.random() > 0.1]  # tips, from dropping a tenth of them
    if not kmers:
        return
    built = unitigs(kmers)
    every_kmer_once(kmers, built)
    # And the answer may not depend on the order the kmers arrive in
    shuffled = kmers[:]
    rng.shuffle(shuffled)
    assert unitigs(shuffled) == built


def test_the_invariants_hold_on_the_hand_written_fixtures():
    """The same two invariants on the examples the other tests use, which assert lengths and sequences
    but never that the decomposition is complete."""
    for kmers in (kmers_of(SEQUENCE),
                  kmers_of(SEQUENCE) + kmers_of(OTHER),
                  kmers_of(SEQUENCE) + [SEQUENCE[20:40] + "A"],
                  kmers_of("CATTGGTGCATAGAGCCTGTACTAGCTTAATGGCATAATCGGTCCACCGGA", k=11)
                  + kmers_of("GCGTTAACGCCCTTTATACTAGCTTAATGGCATAATCGGTCCACCGGA", k=11)):
        every_kmer_once(kmers, unitigs(kmers))
