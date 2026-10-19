"""Assembling the inclusion-specific kmers into contigs: by unitig, or with SKESA or SPAdes."""

from __future__ import annotations

import logging
import shutil
from collections.abc import Iterable
from pathlib import Path

from primer_finder import PrimerFinderError, seqio, tools
from primer_finder.seqio import reverse_complement

log = logging.getLogger(__name__)

ASSEMBLERS = ("skesa", "spades", "unitig")
# The external program each one needs. `unitig` needs none: it is this module.
PROGRAMS: dict[str, str | None] = {"skesa": "skesa", "spades": "spades.py", "unitig": None}

BASES = "ACGT"


def programs_for(assembler: str) -> set[str]:
    """The external programs an assembler needs, which for `unitig` is none at all."""
    program = PROGRAMS.get(assembler)
    return {program} if program else set()

# Where SPAdes stops, measured rather than read off a manual: assembled from every 99-mer of a synthetic
# stretch, stretches of 150, 180 and 199 bases are already absent from its simplified assembly graph, while
# 210 and longer come through. The floor is in graph simplification, not in what it chooses to write, so
# unlike SKESA's --min_contig it cannot be lowered by asking.
SPADES_FLOOR = 200

# The least memory SKESA will start on, measured: `--mem 1` and `--mem 2` both exit 1 with "Memory
# provided is insufficient to do runs in 10 cycles for the read coverage", whatever the input -- it gives
# its kmer counting a third of what it is given and wants a gigabyte of that -- while at `--mem 3` the
# same 1,949-kmer input reports "Memory available (GB): 0.999959" and assembles. So `-a skesa` with `-m 2`
# always ends in the fallback below, and `default_memory_gb()` has its own floor of 2: on a small machine,
# or under a tight control group, the default assembler cannot run at all. The run is told so at the
# start, where it can still be answered, rather than finding out at the assembly step.
SKESA_MIN_MEMORY_GB = 3

# Below how many assembled bases per kmer the chosen assembler is treated as having collapsed rather than
# merely done less well, and the unitig walk is used instead. Measured over 21 runs on three validation
# panels (validation/results/2026-10-14_v1.7.0_unitig/bases_per_kmer.txt): SKESA assembled 0.000, 0.000
# and 0.019 bases per kmer on the three runs where it returned nothing usable, and 0.702 to 1.388 on the
# other eighteen. A 37-fold gap with nothing in it, so the exact value matters little; what matters is
# that no healthy run comes near it.
#
# It deliberately does not try to catch an assembler that is merely leaving a lot behind. On the same
# panels SKESA recovered between 10.7% and 69.8% of the bases the walk found while staying inside the
# healthy band above, so a run cannot tell from its own numbers that it is losing most of the signal.
COLLAPSE_BASES_PER_KMER = 0.1

# Exit codes that mean the machine could not run the assembler, rather than that the assembler could make
# nothing of the kmers. Answering either of those by starting a second assembly in this process -- which
# honours no memory budget at all -- would be the wrong move.
#
# A negative code is a signal, which covers a kernel OOM kill of a native program like SKESA. SPAdes is a
# Python wrapper and never exits by a signal, so the two codes its own machinery produces are listed
# instead. Read out of the installed `spades_pipeline/support.py`:
#
#   * 68, `ErrorCode.MemoryLimitExceeded`: its own --memory check, "I could not stay inside the budget".
#   * 255, `ErrorCode.GeneralError` (-1, and a negative status comes back from the shell as 255). Its
#     `ErrorCode._missing_` maps every value it does not recognise to this one, and a child killed by a
#     signal arrives as such a value -- `ErrorCode(-9)` is `GeneralError` -- so a SPAdes whose assembler
#     was OOM-killed exits 255, not 68. The codes in 64..127, which SPAdes calls `user_end_error`, are
#     the ones that mean it did not like the input, and those still reach the fallback.
FATAL_EXIT_CODES = frozenset({68, 255})


def unitigs(kmers: Iterable[str]) -> list[str]:
    """The maximal unambiguous stretches a set of exact kmers spells out.

    A genome assembler is built for noisy, redundant reads: it estimates coverage, treats what it sees
    once as error, and reconciles disagreement. None of that applies here. These kmers are exact
    substrings of real genomes, each one present or absent with nothing in between, so there is nothing to
    correct and nothing to vote on -- and an assembler's coverage model actively throws them away.
    Measured on one run's 51-mers, SKESA returned 5 contigs where the walk below returns 615, and on 15
    kmers it returned nothing at all.

    What is left once the noise handling goes is the de Bruijn graph on (k-1)-mers and a rule for walking
    it: extend while the step is forced, which means the current kmer has exactly one successor and that
    successor exactly one predecessor. Anywhere the sequence could continue two ways, the unitig ends --
    so a unitig is a stretch the kmers determine, not a guess between alternatives.

    KMC reports a kmer in canonical form, standing for itself and its reverse complement, so the set is
    held that way and both orientations are tried at every step. Successors are found by appending each of
    the four bases and asking whether that kmer is present, which keeps the memory to the kmer set itself.

    Returns the stretches longest first. A kmer set says nothing about which strand a stretch was read
    from, so the orientation each one comes back in is decided here: the walk starts from the stretch's
    lexicographically smallest canonical kmer and reads outward from it, which depends on that stretch's
    own kmers and nothing else. It is therefore stable across runs and unaffected by what else is in the
    set, but it need not match the strand another assembler would choose. Nothing downstream cares: the
    aligner reports either strand and marks the bases accordingly.
    """
    present = {min(kmer, reverse_complement(kmer)) for kmer in (one.upper() for one in kmers)}
    if not present:
        return []
    sizes = {len(kmer) for kmer in present}
    if len(sizes) != 1:
        raise PrimerFinderError(
            f"The kmers are not all the same length ({', '.join(str(one) for one in sorted(sizes))}): "
            "they do not come from one KMC run."
        )

    def known(kmer: str) -> bool:
        return min(kmer, reverse_complement(kmer)) in present

    def after(kmer: str) -> list[str]:
        return [kmer[1:] + base for base in BASES if known(kmer[1:] + base)]

    def before(kmer: str) -> list[str]:
        return [base + kmer[:-1] for base in BASES if known(base + kmer[:-1])]

    def forced_step(step: str, back: bool) -> bool:
        """Whether a step is the only way in and the only way on: the condition for staying in a unitig."""
        return len(after(step) if back else before(step)) == 1

    used: set[str] = set()
    found: list[str] = []
    for seed in sorted(present):  # sorted so that the same kmers always give the same stretches
        if seed in used:
            continue
        # Back up to where the unitig starts, then walk forward building the sequence
        start = seed
        passed = {seed}  # every kmer this walk has backed through, as a canonical key
        while True:
            coming = before(start)
            if len(coming) != 1 or not forced_step(coming[0], back=True):
                break
            step = coming[0]
            key = min(step, reverse_complement(step))
            # Stop rather than go round for ever, which happens in two ways: a circular stretch comes
            # back to the seed, and a hairpin -- a kmer whose last k-1 bases are their own reverse
            # complement -- leads back into the same kmer the other way round. Remembering only the seed
            # caught the first and not the second, which backed the walk through a kmer it had already
            # passed and left the seed itself in no unitig at all.
            if key in used or key in passed:
                break
            passed.add(key)
            start = step
        # Only the base each step adds is kept, not the whole kmer: a stretch of a million kmers would
        # otherwise be held twice over
        grown: list[str] = []
        used.add(min(start, reverse_complement(start)))
        current = start
        while True:
            going = after(current)
            if len(going) != 1 or not forced_step(going[0], back=False):
                break
            step = going[0]
            key = min(step, reverse_complement(step))
            if key in used:
                break
            used.add(key)
            grown.append(step[-1])
            current = step
        found.append(start + "".join(grown))

    # Longest first, and by sequence within a length, so that the order never depends on the order of the
    # input. No deduplication is needed: every kmer is claimed by its canonical key, so two stretches
    # cannot be each other's reverse complement -- they would have to share every key.
    return sorted(found, key=lambda one: (-len(one), one))


def bases_per_kmer(kmer_fasta: Path, output: Path) -> float:
    """How many assembled bases the assembler got out of each kmer it was given.

    A kmer set of N separate stretches spells out `kmers + N(k-1)` bases, so this sits a little above 1
    when an assembler is working on this kind of input and collapses to nearly nothing when it is not.
    It needs no second assembly to compute, which is what makes it usable as a trigger.
    """
    kmers = seqio.count_records(kmer_fasta)
    if not kmers:
        return 0.0
    if not output.is_file():
        return 0.0
    return sum(len(record.seq) for record in seqio.iter_records(output)) / kmers


def assemble(kmer_fasta: Path, output: Path, assembler: str, threads: int, memory_gb: int,
             min_contig: int = 0) -> tuple[Path, str]:
    """Assemble a fasta file of kmers into `output`. Returns `output` and the assembler that filled it.

    `min_contig` is the shortest contig to report. Left alone both assemblers stop at 200 bases, which is
    longer than a kmer: an inclusion-specific stretch shorter than that is then never reported, and nothing
    downstream can see it. Those stretches are short exactly where the two groups are closely related, so
    the caller sets this to the kmer size. 0 leaves each assembler's own default alone.

    **`unitig` and SKESA honour it; SPAdes cannot.** `unitig` builds every stretch the kmers determine and
    then drops the short ones, so the floor is exactly what is asked for and the shortest it can be is a
    kmer. `--min_contig` tells SKESA what to report, so asking for less gets less. SPAdes has no such
    option and its floor is in graph simplification (`SPADES_FLOOR`), so a stretch shorter than that never
    reaches its output and asking cannot change it; here the value is only applied to what SPAdes did
    write, which keeps the assemblers from reporting *longer* contigs than asked for. `-a spades` with a
    floor below `SPADES_FLOOR` is warned about rather than silently ignored.

    **If the chosen assembler collapses, the unitig walk is used instead and the return says so.** An
    assembler built for sequencing reads can get nothing at all out of a kmer set -- it fails outright, or
    writes a handful of bases out of tens of thousands of kmers -- and the error it used to raise advised
    a smaller `-k`, which runs into the same assembler from the other end. `COLLAPSE_BASES_PER_KMER` says
    where "collapsed" starts, and the measurements behind it are beside that constant. This is not a
    judgement about which assembler is better on a run that worked: the fallback is only for the case
    where there is nothing to lose.
    """
    used = assembler
    failure: tools.ToolError | None = None
    try:
        if assembler == "unitig":
            _unitig(kmer_fasta, output, min_contig)
        elif assembler == "skesa":
            _skesa(kmer_fasta, output, threads, memory_gb, min_contig)
        elif assembler == "spades":
            if 0 < min_contig < SPADES_FLOOR:
                log.warning("SPAdes does not report a contig shorter than about %d bases and has no "
                            "option for it, so regions down to %d will be missing. Use the default "
                            "assembler to get them", SPADES_FLOOR, min_contig)
            _spades(kmer_fasta, output, threads, memory_gb, min_contig)
        else:  # pragma: no cover - the command line only accepts the three of ASSEMBLERS
            raise PrimerFinderError(
                f'Unknown assembler "{assembler}". Choose one of: {", ".join(ASSEMBLERS)}')
    except tools.ToolError as error:
        # An assembler that exits non-zero on a kmer set has collapsed in the same way as one that writes
        # nothing, and is answered the same way rather than ending the run -- unless it was killed, or
        # could not fit in the memory it was given. See FATAL_EXIT_CODES.
        if error.returncode is not None and (error.returncode < 0
                                             or error.returncode in FATAL_EXIT_CODES):
            raise
        failure = error

    if assembler != "unitig":
        recovered = 0.0 if failure is not None else bases_per_kmer(kmer_fasta, output)
        if failure is not None or recovered < COLLAPSE_BASES_PER_KMER:
            reason = (f"it failed: {failure}" if failure is not None
                      else f"it assembled {recovered:.3f} bases per kmer, where a working run gives "
                           f"about 1")
            # Deliberately not "it got nothing out of the kmers": a floor higher than every contig
            # empties a perfectly healthy assembler too, and from here the two look the same
            log.warning("Nothing usable came out of %s -- %s", assembler, reason)
            log.warning("Assembling the kmers directly instead (as -a unitig does), which needs no "
                        "minimum number of them and no coverage model. If --min-region is what emptied "
                        "the assembly rather than the assembler, this will say so and stop. "
                        "run_info.json records which assembler the contigs came from.")
            _unitig(kmer_fasta, output, min_contig)
            used = "unitig"

    if not output.exists() or output.stat().st_size == 0:
        raise PrimerFinderError(
            f"{used} could not assemble the inclusion-specific kmers into contigs. A single kmer is "
            "already a unitig, so with a floor of the kmer size this happens only when --min-region is "
            "longer than every stretch the kmers spell out."
        )
    return output, used


def _unitig(kmer_fasta: Path, output: Path, min_contig: int) -> None:
    """Walk the kmers into unitigs in this process, with no external assembler."""
    built = unitigs(record.seq for record in seqio.iter_records(kmer_fasta))
    keep = [one for one in built if len(one) >= min_contig] if min_contig else built
    if len(keep) < len(built):
        log.info("%d of %d unitig(s) are shorter than %d bases and are left out",
                 len(built) - len(keep), len(built), min_contig)
    seqio.write_fasta(output, [seqio.Record(f"unitig_{number}", "", one)
                               for number, one in enumerate(keep, start=1)])


def _skesa(kmer_fasta: Path, output: Path, threads: int, memory_gb: int, min_contig: int) -> None:
    tools.run([
        "skesa",
        "--cores", str(threads),
        "--mem", str(memory_gb),
        "--fasta", kmer_fasta,
        "--contigs_out", output,
        *(["--min_contig", str(min_contig)] if min_contig else []),
    ])


def _spades(kmer_fasta: Path, output: Path, threads: int, memory_gb: int, min_contig: int) -> None:
    work_dir = output.parent / "spades"
    # Cleared first, not only on the way out: the folder is only removed when SPAdes succeeded, so a run
    # that failed and was answered by the fallback used to leave `spades/contigs.fasta` sitting beside the
    # walk's assembly -- contigs this run did not use, named as though it had. SPAdes would also have
    # taken the leftovers as a half-finished run of its own to continue from.
    shutil.rmtree(work_dir, ignore_errors=True)
    tools.run([
        "spades.py",
        "--s", "1", kmer_fasta,
        "--isolate",
        "--only-assembler",  # The kmers carry no quality values: there is nothing to correct
        "--threads", str(threads),
        "--memory", str(memory_gb),
        "-o", work_dir,
    ])
    contigs = work_dir / "contigs.fasta"
    if not contigs.exists():
        # A ToolError, not a plain error, so that it reaches the fallback: SPAdes exiting 0 without
        # writing contigs is the same situation as SPAdes exiting 1
        raise tools.ToolError(f"SPAdes wrote no contigs ({contigs} missing)")
    shutil.move(str(contigs), output)
    shutil.rmtree(work_dir, ignore_errors=True)
    if min_contig:
        # SPAdes has no option for this. Applying it to the output cannot bring a short contig back -- it
        # was lost in graph simplification -- but it does keep `-a spades` from reporting contigs the
        # command was asked to leave out
        keep = [record for record in seqio.read_fasta(output).values() if len(record.seq) >= min_contig]
        seqio.write_fasta(output, keep)
