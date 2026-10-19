"""The pipeline: count kmers, assemble them, then check the contigs against every genome."""

from __future__ import annotations

import json
import logging
import math
import shutil
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from primer_finder import (
    NothingSelectiveEnough,
    PrimerFinderError,
    __version__,
    assemble,
    blast,
    kmers,
    mapping,
    seqio,
    system,
    tools,
)
from primer_finder.mapping import MAX_PRIMER_SPAN

log = logging.getLogger(__name__)

LOG_NAME = "primer_finder.log"
RUN_INFO_NAME = "run_info.json"
FINAL_NAME = "final_kmers.fasta"
BEST_NAME = "best_kmers.fasta"
ALL_INCLUSION_NAME = "all_inclusion_contigs.fasta"
HITS_NAME = "inclusion_blast_hits.tsv"
# What a record's description says, with -p/--min-inclusion, about how many inclusion genomes hold it.
INCLUSION_TAG = "inclusion="

# What marks a region that only the last-resort rule accepted, in its header in final_kmers.fasta. The
# whole point is that it cannot be missed by anything reading the file, this tool included.
LAST_RESORT_TAG = "last_resort"

# What an earlier run's results are renamed to when a new run replaces them.
PREVIOUS_SUFFIX = ".previous"

KMER_FOLDER = "1_kmers"
ASSEMBLY_FOLDER = "2_assembly"
CANDIDATE_FOLDER = "3_candidates"
BLAST_FOLDER = "4_blast"

BASE_PROGRAMS = ("kmc", "kmc_tools", "minimap2", "makeblastdb", "blastn")


def required_genomes(settings: Settings, inclusion: list[Path]) -> int:
    """How many of the inclusion genomes a kmer, and later a contig, has to be in: every one of them by
    default, or the fraction -p/--min-inclusion asks for, rounded up.

    The product is rounded to nine decimals before rounding up, because a fraction such as 0.14 is a little
    more than 0.14 in binary: 0.14 * 50 is 7.000000000000001, which would ask for 8 genomes out of 50.
    """
    return max(1, math.ceil(round(settings.min_inclusion * len(inclusion), 9)))


@dataclass
class Settings:
    """Everything the pipeline needs, as the command line gave it."""

    inclusion: Path
    exclusion: Path
    output: Path
    threads: int
    memory_gb: int
    kmer_size: int = 99
    duplication: int = 1
    min_inclusion: float = 1.0  # The fraction of the inclusion genomes a kmer has to be in
    min_region: int = 0  # The shortest region to report; 0 for the kmer size
    reference: Path | None = None
    assembler: str = "skesa"
    # Whether to try again accepting a single difference when the strict rule reports nothing at all
    last_resort: bool = True
    keep_intermediate: bool = False
    command_line: list[str] = field(default_factory=list)


def run(settings: Settings) -> int:
    """Run the whole pipeline. Returns 0, or raises PrimerFinderError.

    The settings are checked and the output folder made **before** the clean-up is armed. A run refused
    before it has written anything must not tidy the folder: the same `-o` may hold the results of an
    earlier run, and `clean_up` deletes by fixed path, so a typo in `-i` used to take the previous run's
    `--keep-intermediate` files with it.

    The folder is then claimed for this run, before anything in it is written or moved, so that two runs
    cannot share it; and the earlier log is moved aside before this run's log replaces it.
    """
    inclusion, exclusion, reference = check(settings)
    try:
        settings.output.mkdir(parents=True, exist_ok=True)
    except (FileExistsError, NotADirectoryError) as error:
        # mkdir's own message is not obvious: a file where the folder should be gives FileExistsError,
        # and a file part way along the path gives NotADirectoryError
        raise PrimerFinderError(
            f"The output path is a file, or sits under one, rather than a folder: {settings.output}"
        ) from error
    # Both of these can fail on a folder an earlier run wrote -- a read-only folder, or a `.previous`
    # path that is a directory -- and failing inside the block below would run `clean_up` on that earlier
    # run's `--keep-intermediate` files, which is the very thing the ordering here exists to prevent.
    if settings.assembler == "skesa" and settings.memory_gb < assemble.SKESA_MIN_MEMORY_GB:
        # Said here, before any of the work, rather than left to the assembly step: there the fallback
        # reported "nothing usable came out of skesa", which reads as a judgement on the kmers
        log.warning("SKESA refuses to start on the %d GB of -m/--memory, whatever the kmers are, so the "
                    "assembly will fall back to the unitig walk (it needs at least %d GB). Raise -m if "
                    "the machine has the memory, or -a unitig to ask for the walk outright.",
                    settings.memory_gb, assemble.SKESA_MIN_MEMORY_GB)
    with system.lock_output(settings.output):
        keep_previous_log(settings.output / LOG_NAME)
        add_log_file(settings.output / LOG_NAME)
        replace_previous_results(settings)
        try:
            return _run(settings, inclusion, exclusion, reference)
        finally:
            # On the way out either way: a run that ends at "no contig differs enough" used to leave its
            # KMC databases behind, and an exclusion database for a few hundred genomes is a gigabyte
            clean_up(settings)


def _run(settings: Settings, inclusion: list[Path], exclusion: list[Path], reference: Path) -> int:
    started = time.time()
    log.info("primer-finder %s", __version__)
    log.info("%d inclusion genome(s), %d exclusion genome(s)", len(inclusion), len(exclusion))

    info: dict[str, object] = {
        "version": __version__,
        "command_line": settings.command_line,
        "parameters": {
            # Absolute, because `design` reads these back and resolves them against *its own* working
            # directory. Recorded as typed, a run started with `-i inclusion -e exclusion` left `design`
            # either stopping ("the inclusion folder of that run is gone") or, where a folder of that
            # name happened to exist beside it, designing against the wrong genomes with no warning at
            # all -- checking amplicons and running in silico PCR on a panel the run never saw.
            "inclusion": str(settings.inclusion.resolve()),
            "exclusion": str(settings.exclusion.resolve()),
            "output": str(settings.output.resolve()),
            "kmer_size": settings.kmer_size,
            "duplication": settings.duplication,
            "min_inclusion": settings.min_inclusion,
            "min_region": settings.min_region or settings.kmer_size,
            "inclusion_genomes_required": required_genomes(settings, inclusion),
            "assembler": settings.assembler,
            "threads": settings.threads,
            "memory_gb": settings.memory_gb,
        },
        # Absolute, but not resolved: a panel of links into a shared assembly store would otherwise be
        # recorded under the store's names, where the log, the presence table's columns and KMC's own file
        # list all show the panel's. Only the folders above need resolving, for `design` to find them.
        "inclusion_genomes": [str(path.absolute()) for path in inclusion],
        "exclusion_genomes": [str(path.absolute()) for path in exclusion],
        "reference": str(reference.absolute()),
        "programs": tools.versions(sorted({*BASE_PROGRAMS, *assemble.programs_for(settings.assembler)})),
        "counts": {},
    }
    counts: dict[str, int] = info["counts"]  # type: ignore[assignment]

    kmer_fasta, counts["kmers"] = find_specific_kmers(settings, inclusion, exclusion)
    assembly, assembler_used = assemble_kmers(settings, kmer_fasta)
    info["parameters"]["assembler_used"] = assembler_used  # type: ignore[index]
    counts["contigs"] = seqio.count_records(assembly)

    absent: dict[str, list[list[int]]] = {}
    final, min_differences = check_the_contigs(settings, assembly, reference, inclusion, exclusion,
                                               counts, absent)
    # Which marked stretches are missing from the exclusion group rather than differing from it. `design`
    # reads the marks off the sequence and cannot tell the two apart on its own.
    info["absent_stretches"] = {name: runs for name, runs in sorted(absent.items())}
    info["parameters"]["min_differences"] = min_differences  # type: ignore[index]
    info["last_resort"] = min_differences <= 1 and counts["final"] > 0

    info["seconds"] = round(time.time() - started, 1)
    (settings.output / RUN_INFO_NAME).write_text(json.dumps(info, indent=2) + "\n")
    log.info("Final number of contigs: %d", counts["final"])
    log.info("Results: %s", final)
    if counts["final"] == 0:
        log.warning("No contig passed every filter: there is no inclusion-specific assay to design here.")
    return 0


def check(settings: Settings) -> tuple[list[Path], list[Path], Path]:
    """Check the arguments and the input folders, and choose the exclusion genome to map the contigs to."""
    if not kmers.MIN_KMER_SIZE <= settings.kmer_size <= kmers.MAX_KMER_SIZE:
        raise PrimerFinderError(
            f"-k/--kmer_size must be between {kmers.MIN_KMER_SIZE} and {kmers.MAX_KMER_SIZE} "
            f"(KMC's range), not {settings.kmer_size}"
        )
    if settings.duplication < 1:
        raise PrimerFinderError(f"-d/--duplication must be 1 or more, not {settings.duplication}")
    if not 0 < settings.min_inclusion <= 1:
        raise PrimerFinderError(
            f"-p/--min-inclusion must be a fraction greater than 0 and at most 1, not {settings.min_inclusion}"
        )
    if settings.min_region < 0:
        raise PrimerFinderError("--min-region cannot be negative; 0 means the kmer size")
    if settings.assembler not in assemble.ASSEMBLERS:
        raise PrimerFinderError(
            f'Unknown assembler "{settings.assembler}". Choose one of: {", ".join(assemble.ASSEMBLERS)}'
        )
    inclusion = seqio.require_genomes(settings.inclusion, "inclusion")
    exclusion = seqio.require_genomes(settings.exclusion, "exclusion")
    shared = {path.resolve() for path in inclusion} & {path.resolve() for path in exclusion}
    if shared:
        raise PrimerFinderError(
            "The same genome is in both groups: " + ", ".join(sorted(str(path) for path in shared))
        )
    output = settings.output.resolve()
    # makeblastdb reopens the database it has just written by its absolute path, and BLAST splits paths on
    # whitespace, so a blast database cannot live under a folder whose name holds a space. The input
    # folders may: every genome is linked into the output folder before blast sees it.
    if any(character in blast.UNSAFE_IN_DB_PATH for character in str(output)):
        raise PrimerFinderError(
            f"The output folder path contains a space or a quote, which BLAST cannot use: "
            f"{settings.output}. Choose an output folder whose path has neither (the input paths and the "
            "genome file names may have both: each genome is linked into a folder named after its "
            "position in the list)."
        )
    for group, folder in (("inclusion", settings.inclusion), ("exclusion", settings.exclusion)):
        root = folder.resolve()
        if _inside(output, root) or _inside(root, output):
            raise PrimerFinderError(
                f"The output folder and the {group} folder overlap ({settings.output} and {folder}). "
                "primer-finder searches the input folders recursively, so it would read its own results as "
                "input genomes. Choose an output folder outside both input folders."
            )

    if settings.reference is not None:
        wanted = settings.reference.resolve()
        match = [path for path in exclusion if path.resolve() == wanted]
        if not match:
            raise PrimerFinderError(
                f"-r/--reference must be one of the genomes in the exclusion folder, which {settings.reference} "
                "is not"
            )
        reference = match[0]
    else:
        reference = exclusion[0]  # The first one in alphabetical order, so that a run can be repeated

    tools.require(sorted({*BASE_PROGRAMS, *assemble.programs_for(settings.assembler)}))
    return inclusion, exclusion, reference


def _inside(child: Path, parent: Path) -> bool:
    """True if `child` is `parent` or sits below it. Both paths must be resolved."""
    return child == parent or parent in child.parents


def add_log_file(path: Path) -> None:
    """Also write the log to a file in the output folder."""
    root = logging.getLogger()
    for existing in [h for h in root.handlers if getattr(h, "_primer_finder", False)]:
        root.removeHandler(existing)  # A second run in the same process logs to its own file
        existing.close()
    handler = logging.FileHandler(path, mode="w")
    handler._primer_finder = True  # type: ignore[attr-defined]
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S"))
    root.addHandler(handler)
    if root.level == logging.NOTSET or root.level > logging.INFO:
        root.setLevel(logging.INFO)  # The log file is part of the results, whatever the caller set up


def find_specific_kmers(settings: Settings, inclusion: list[Path], exclusion: list[Path]) -> tuple[Path, int]:
    """Count the kmers of both groups with KMC, subtract the exclusion ones, and write what is left as a
    fasta file of kmers."""
    folder = settings.output / KMER_FOLDER
    work_dir = folder / "kmc_work"
    folder.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    required = required_genomes(settings, inclusion)
    if required < len(inclusion):
        log.info("Counting the %d-mers shared by at least %d of the %d inclusion genomes (-p %g)...",
                 settings.kmer_size, required, len(inclusion), settings.min_inclusion)
    else:
        log.info("Counting the %d-mers shared by the %d inclusion genomes...",
                 settings.kmer_size, len(inclusion))
    inclusion_db = kmers.count(
        kmers.write_file_list(inclusion, folder / "inclusion_list.txt"),
        folder / "inclusion", work_dir, settings.kmer_size, settings.threads, settings.memory_gb,
        min_count=required, max_count=len(inclusion) * settings.duplication,
    )
    log.info("Counting the %d-mers of the %d exclusion genomes...", settings.kmer_size, len(exclusion))
    exclusion_db = kmers.count(
        kmers.write_file_list(exclusion, folder / "exclusion_list.txt"),
        folder / "exclusion", work_dir, settings.kmer_size, settings.threads, settings.memory_gb,
        min_count=1, max_count=kmers.NO_MAX_COUNT,
    )
    log.info("Subtracting the exclusion kmers from the inclusion ones...")
    specific_db = kmers.subtract(inclusion_db, exclusion_db, folder / "inclusion_specific", settings.threads)

    # An empty result can be the fault of either group, so if nothing is left, count what the inclusion
    # group shared before the subtraction to say which of the two it was
    dump_file = kmers.dump(
        specific_db, folder / "kmers.txt", settings.threads,
        count_shared=lambda: kmers.count_in_db(inclusion_db, folder / "inclusion_kmers.txt",
                                               settings.threads),
    )
    kmer_fasta = folder / f"inclusion_specific_{settings.kmer_size}mers.fasta"
    count = kmers.dump_to_fasta(dump_file, kmer_fasta)
    # Here rather than inside `dump`, so that the kmers are counted once: a dump of the inclusion group
    # of a large panel is gigabytes of text, and reading it twice to check its length is reading it twice
    kmers.require_whole_dump(dump_file, count)
    log.info("Found %d inclusion-specific %d-mers", count, settings.kmer_size)
    return kmer_fasta, count


def assemble_kmers(settings: Settings, kmer_fasta: Path) -> tuple[Path, str]:
    """Assemble the inclusion-specific kmers, to have fewer and longer sequences to check.

    Returns the assembly and the assembler that filled it, which is not always the one that was asked
    for: an assembler that gets nothing out of a kmer set is answered by assembling them directly.
    """
    folder = settings.output / ASSEMBLY_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    assembly = folder / "assembly.fasta"
    log.info("Assembling the kmers with %s...", settings.assembler)
    _, used = assemble.assemble(kmer_fasta, assembly, settings.assembler, settings.threads,
                                settings.memory_gb,
                                min_contig=settings.min_region or settings.kmer_size)
    log.info("Assembled %d contig(s)%s", seqio.count_records(assembly),
             "" if used == settings.assembler else f" (with {used}, not {settings.assembler})")
    return assembly, used


def check_the_contigs(settings: Settings, assembly: Path, reference: Path, inclusion: list[Path],
                      exclusion: list[Path], counts: dict[str, int],
                      absent: dict[str, list[list[int]]] | None = None) -> tuple[Path, int]:
    """Steps 3 to 5 at the strict rule, and again accepting a single difference if that found nothing.

    The three steps between the assembly and the answer all turn on how much a region has to differ, and
    any of them can empty the output. When they do, the run has a choice: stop, or say what it would have
    found had it asked for less. Asking for less is worth doing *and* worth admitting to, so the second
    pass marks everything it reports -- in the headers, in `run_info.json` and in the log -- and the
    design command repeats the warning. `--no-last-resort` turns it off for a run that would rather have
    the error.

    Only `NothingSelectiveEnough` is answered this way. A failed program is a `PrimerFinderError` too and
    must stop the run, which is why that error has a class of its own.

    Returns the regions file and the rule they were found under.
    """
    def pass_at(min_differences: int) -> Path:
        best, counts["candidates"] = find_candidates(settings, assembly, reference, min_differences)
        all_inclusion, counts["in_all_inclusion"] = keep_shared_by_inclusion(settings, best, inclusion)
        final, counts["final"], found = keep_absent_from_exclusion(settings, all_inclusion, exclusion,
                                                                    min_differences)
        if absent is not None:
            absent.clear()          # A second pass replaces what the first one found
            absent.update(found)
        return final

    strict = mapping.STRICT_DIFFERENCES
    try:
        final = pass_at(strict)
        if counts["final"] > 0 or not settings.last_resort:
            return final, strict
        reason = "no region's differences held up against every exclusion genome"
    except NothingSelectiveEnough as error:
        if not settings.last_resort:
            raise
        reason = str(error)

    # Deliberately not "nothing met the two-differences rule": any of the three steps can empty the run,
    # and when it is step 4 the contigs did meet that rule and were simply not in every inclusion genome.
    # The reason says which, so the sentence around it must not claim to know.
    log.warning("The run has nothing to report at the usual stringency. %s", reason)
    log.warning("Trying again as a last resort, accepting a region with a single difference instead of "
                "two close enough to share a primer. An assay resting on one mismatch is the harder one "
                "to make work, so everything this reports is marked %r and should be treated as a lead "
                "rather than a result. --no-last-resort skips this.", LAST_RESORT_TAG)
    try:
        final = pass_at(mapping.LAST_RESORT_DIFFERENCES)
    except NothingSelectiveEnough as error:
        # Otherwise this is word for word the error the strict pass gives, and nothing in it says that a
        # second, weaker attempt has already been made
        raise NothingSelectiveEnough(f"The last-resort pass found nothing either. {error}") from error
    if counts["final"] > 0:
        log.warning("%d region(s) found only by the last-resort rule, each marked %r in %s",
                    counts["final"], LAST_RESORT_TAG, final.name)
    else:
        log.warning("The last-resort pass found nothing either: no region is left to design on.")
    return final, mapping.LAST_RESORT_DIFFERENCES


def find_candidates(settings: Settings, assembly: Path, reference: Path,
                    min_differences: int = mapping.STRICT_DIFFERENCES) -> tuple[Path, int]:
    """Map the contigs to one exclusion genome and keep those whose differences could fit in a primer."""
    folder = settings.output / CANDIDATE_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    log.info("Mapping the contigs to the exclusion genome %s...", reference.name)
    sam_file = mapping.map_contigs(assembly, reference, folder / "mapping.sam", settings.threads)

    if min_differences > 1:
        log.info("Keeping the contigs with at least two differences within %d bases...", MAX_PRIMER_SPAN)
    else:
        log.info("Last resort: keeping the contigs with even a single difference...")
    candidates = mapping.select_candidates(seqio.read_fasta(assembly), mapping.parse_sam(sam_file),
                                           min_differences)
    best = folder / BEST_NAME
    count = seqio.write_fasta(best, (candidate.as_record() for candidate in mapping.sort_candidates(candidates)))
    log.info("%d contig(s) could carry a selective assay", count)
    if count == 0:
        asked = (f"at least two differences within {MAX_PRIMER_SPAN} bases" if min_differences > 1
                 else "even a single difference")
        raise NothingSelectiveEnough(
            "No contig differs enough from the exclusion genome to carry a selective assay "
            f"({asked})."
        )
    return best, count


def keep_shared_by_inclusion(settings: Settings, best: Path, inclusion: list[Path]) -> tuple[Path, int]:
    """Drop the candidates that too few inclusion genomes hold: the mapping step only compared them with one
    exclusion genome, and the assembly may hold kmers that come from some of the genomes only.

    By default a candidate has to be in every inclusion genome. With -p/--min-inclusion below 1, it has to
    be in that fraction of them, the ones that hold it come first, and each record says how many hold it.
    """
    folder = settings.output / BLAST_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    required = required_genomes(settings, inclusion)
    log.info("Checking the candidates against the %d inclusion genomes...", len(inclusion))
    presence = blast.presence_in_genomes(best, inclusion, folder / "inclusion_db", settings.threads)
    write_presence_table(folder / HITS_NAME, presence, genome_labels(inclusion, settings.inclusion))

    partial = required < len(inclusion)
    kept: list[tuple[int, seqio.Record]] = []
    for record in seqio.iter_records(best):
        # Any hit at all is the contig being there, which is the rule this has always applied. See
        # `report_partial_presence` for why it is not a fraction of the contig.
        found = sum(1 for covered in presence[record.name].values() if covered > 0)
        if found < required:
            continue
        if partial:  # The candidates are no longer equal: say what each one covers
            record = seqio.Record(record.name, f"{record.desc} {INCLUSION_TAG}{found}/{len(inclusion)}".strip(),
                                  record.seq)
        kept.append((found, record))
    if partial:  # Most inclusion genomes first, keeping the order of the mapping step within a count
        kept.sort(key=lambda pair: -pair[0])
    all_inclusion = folder / ALL_INCLUSION_NAME
    count = seqio.write_fasta(all_inclusion, (record for _, record in kept))

    if partial:
        log.info("%d contig(s) are present in at least %d of the %d inclusion genomes",
                 count, required, len(inclusion))
    else:
        log.info("%d contig(s) are present in all inclusion genomes", count)
    report_partial_presence(presence, (record for _, record in kept), folder / HITS_NAME)
    if count == 0:
        how_many = f"at least {required} of the inclusion genomes" if partial else "all the inclusion genomes"
        hint = (f" A -p/--min-inclusion below {settings.min_inclusion:g} would ask for fewer than {required}."
                if partial else
                " A lower -p/--min-inclusion would accept a contig that some of them lack.")
        raise NothingSelectiveEnough(
            f"No candidate contig is present in {how_many}. The presence of each one in each genome is in "
            f"{folder / HITS_NAME}." + hint
        )
    return all_inclusion, count


def genome_labels(genomes: list[Path], root: Path) -> dict[str, str]:
    """A short name for each genome, for the columns of the presence table: its path inside the input
    folder, so that two genomes with the same file name in different subfolders stay apart."""
    labels: dict[str, str] = {}
    for genome in genomes:
        try:
            label = str(genome.relative_to(root))
        except ValueError:  # pragma: no cover - find_genomes returns paths inside the folder
            label = genome.name
        labels[str(genome)] = label
    return labels


# How many exclusion genomes make "absent from the group" worth much. Not a limit -- nothing is refused
# below it -- only the point under which the log says what the evidence is worth. Measured on one
# *Listeria* complex by running it against exclusion groups sampled from its own 716 genomes and testing
# every absence against all of them: with 1 to 3 genomes, 95 to 97% of the stretches reported absent are
# held somewhere in the full set; with 10, 44%; with 25, 17%; with all 716, none. The real absences were
# found at every group size, so what grows as the group shrinks is the false ones.
ABSENCE_WANTS_GENOMES = 25
# An inclusion genome that aligns less than this much of a region holds it only partly, and the run says
# so. The same figure as `design.MIN_AMPLICON_COVERAGE`, which is the rule this one reports against: what
# the design step requires of an amplicon is what makes a part of a region worth telling apart from the
# whole of it. It is a reporting threshold only -- nothing is kept or dropped by it.
MOSTLY_PRESENT = 0.9


def warn_about_absence(regions: int, beside_differences: int, genomes: int) -> None:
    """Say what a region resting on an absent stretch is worth, given how many genomes said so.

    The differences rule needs 90% of the genomes covering a position to agree, so no single genome can
    carry it. A stretch missing from the group is credited on whatever the group says -- with one genome,
    one alignment gap of 6,888 bases becomes an absence. It is a true statement about the genomes given
    and it generalises only as far as they do, which is worth saying when they are few.
    """
    if genomes >= ABSENCE_WANTS_GENOMES or not (regions or beside_differences):
        return
    # Both counts are said. A region whose own differences meet the rule does not *rest* on its absent
    # stretch -- but `design` aims its oligos at the marked bases, and an absent block is marked in full,
    # so it is usually the richest target and the assays designed there rest on it even though the region
    # did not. Naming only the first kind made the warning disappear for a region with two differences
    # beside a 400-base absent block, on which every assay then came out specific by absence.
    log.warning(
        "%d region(s) are specific because a stretch of them is missing from all %d exclusion genome(s), "
        "and %d more hold such a stretch beside differences of their own, which `design` will aim at "
        "first. That is a statement about those %d genomes: with a group this small such a stretch is "
        "often present elsewhere in the species -- measured on one Listeria complex, 95%% of the stretches "
        "reported absent against 3 genomes were held by others of its 716. Widen the exclusion group, or "
        "check these regions against more genomes before trusting the specificity.",
        regions, genomes, beside_differences, genomes)


def runs_of(positions: Sequence[int]) -> list[list[int]]:
    """Consecutive positions as [first, last] pairs."""
    runs: list[list[int]] = []
    for position in positions:
        if runs and position == runs[-1][1] + 1:
            runs[-1][1] = position
        else:
            runs.append([position, position])
    return runs


def write_presence_table(path: Path, presence: dict[str, dict[str, float]],
                         labels: dict[str, str] | None = None) -> None:
    """A table of contigs (rows) against genomes (columns), holding the fraction of the contig each genome
    aligns: 0 for absent, 1 for all of it.

    It was 1 and 0 until 1.9.2. Presence is still decided on anything above 0, so the rows mean what they
    meant; what the fraction adds is that a genome holding a part of a region can be told from one holding
    the region, which the tick could not say and which is the difference that matters where an assembly
    has broken a repeat in two.

    Three decimals, so a genome holding under a two-thousandth of a region reads as `0.000` like one
    holding none of it. That needs a region of some 66,000 bases hit over the 33 the search can see, which
    no run here has produced; the figure to go by in that case is the log line, which counts bases.
    """
    labels = labels or {}
    genomes = sorted({genome for hits in presence.values() for genome in hits},
                     key=lambda genome: labels.get(genome, genome))
    with path.open("w") as fh:
        fh.write("contig\t" + "\t".join(labels.get(genome, genome) for genome in genomes) + "\n")
        for contig in sorted(presence):
            row = [f"{presence[contig].get(genome, 0.0):.3f}" for genome in genomes]
            fh.write(f"{contig}\t" + "\t".join(row) + "\n")


def report_partial_presence(presence: dict[str, dict[str, float]], kept: Iterable[seqio.Record],
                            table: Path) -> int:
    """Say how many of the kept contigs an inclusion genome holds only part of. Returns that count.

    A contig counts as present in a genome on any hit at all, with no requirement on how much of it that
    hit covers, and that is left alone on purpose: the design step checks each amplicon it proposes
    against every genome separately, so a region whose shared part carries a usable assay is worth
    keeping, and asking for most of the region here would drop it before anything looked at the part the
    genomes agree on. What was missing is that the output could not tell the two apart.

    Measured on the draft-heavy *Brucella* group of `validation/results/2026-10-16_v1.9.1_is711/`: a
    338-base region was reported as present in all thirty inclusion genomes while one of them, a
    contig-level *B. canis* assembly, held 41 identical bases of it at the start of a contig and not the
    published amplicon at all. The region was worth reporting -- the other genomes carry it -- and the
    run had no way of saying which case it was.
    """
    worst: tuple[float, str, int] | None = None
    partial = 0
    for record in kept:
        # Only the genomes that hold some of it. Below -p 1.0 a region is kept although some genomes have
        # none of it, and that is what the fraction in its header is for; this line is about the genomes
        # that do hold it, which is a different thing and the one the table cannot be read without.
        held = [covered for covered in presence.get(record.name, {}).values() if covered > 0]
        if not held:  # pragma: no cover - a contig held by no genome is never kept
            continue
        least = min(held)
        if least < MOSTLY_PRESENT:
            partial += 1
            if worst is None or least < worst[0]:
                worst = (least, record.name, len(record.seq))
    if worst is not None:
        # Bases rather than a percentage: a rounded percentage of a long region reads as a contradiction
        # -- "as little as 0%" for a genome that does hold 33 bases of a 7 kb region, or "90%" for one
        # just under the threshold this line is about.
        log.info("%d of them are held only in part by an inclusion genome that holds them at all: as "
                 "little as %d of the %d bases of %s. A region can be worth reporting on the part its "
                 "genomes share -- the design step checks each amplicon against every genome -- but see "
                 "%s for which genome holds how much.",
                 partial, round(worst[0] * worst[2]), worst[2], worst[1], table)
    return partial


def keep_absent_from_exclusion(
    settings: Settings, all_inclusion: Path, exclusion: list[Path],
    min_differences: int = mapping.STRICT_DIFFERENCES,
) -> tuple[Path, int, dict[str, list[list[int]]]]:
    """Keep the contigs whose differences hold up against every exclusion genome, not just the one they were
    mapped to. A contig that no exclusion genome hits at all is kept as it is.

    In a last-resort pass (`min_differences` of 1) every region reported is tagged `LAST_RESORT_TAG` in its
    header. Nothing in such a pass would have passed the strict rule -- the strict pass ran first and
    reported nothing -- so the tag belongs on all of them, including a region specific by absence.
    """
    mark = LAST_RESORT_TAG if min_differences <= 1 else ""
    folder = settings.output / BLAST_FOLDER
    log.info("Checking the differences against the %d exclusion genomes...", len(exclusion))
    results = blast.exclusion_variants(all_inclusion, exclusion, folder / "exclusion_db", settings.threads)
    candidates = seqio.read_fasta(all_inclusion)

    records: list[seqio.Record] = []
    resting_on_absence: list[str] = []
    absence_beside_differences: list[str] = []
    absent_by_region: dict[str, list[list[int]]] = {}
    for name, record in candidates.items():
        result = results.get(name, blast.ExclusionResult())
        if result.genomes_hit == 0:  # Absent from every exclusion genome: specific on its own
            # The whole of it is the absent stretch, and saying so is what stops `design` reading the
            # region's marked bases as mismatches. Until 1.9.2 this branch recorded nothing, so a region
            # no exclusion genome holds came out of `design` claiming a difference at every base: on a
            # 289-base region of the IS711 validation panel, `forward_variants 20`, `forward_terminal_run
            # 20` and a `primer_variant_weight` of 30 over a primer whose template is not there at all.
            # That is the misreport 1.9.1 set out to end; it only covered the partly-absent regions.
            absent_by_region[name] = [[0, len(record.seq) - 1]]
            resting_on_absence.append(name)
            records.append(seqio.Record(record.name, f"{record.desc} {mark}".strip(), record.seq)
                           if mark else record)
            continue
        # A stretch no exclusion genome aligns counts as differing, which is what it does: the whole of
        # it is missing from the group. Without this a region whose specificity is a novel block in shared
        # backbone -- a repeat-targeted marker, most often -- was dropped for having no differences, while
        # the mapping step two steps earlier had called the same contig a candidate on that very block.
        absent = blast.unaligned_positions(result, len(record.seq))
        differences = sorted(blast.shared_variants(result))
        if absent:
            absent_by_region[name] = runs_of(absent)
            # "Resting on" means the absence is what carries it. A region whose differences alone meet the
            # rule does not rest on anything, and counting it made the warning name the wrong regions as
            # the fragile ones -- but it is still worth counting separately, since the assays designed on
            # it will aim at the absent block.
            if blast.enough_variants(differences, min_differences):
                absence_beside_differences.append(name)
            else:
                resting_on_absence.append(name)
        positions = sorted(set(differences) | set(absent))
        if blast.enough_variants(positions, min_differences):
            tag = next((part for part in record.desc.split() if part.startswith(INCLUSION_TAG)), "")
            desc = " ".join(part for part in (f"{positions}", tag, mark) if part)
            records.append(seqio.Record(name, desc, lower_positions(record.seq, positions)))
    kept = {record.name for record in records}
    # Returned for `design`, which reads the marked positions off the sequence and cannot otherwise tell a
    # base that *differs* from one that is *missing* -- and counted them alike, so an assay inside an
    # absent block came out with twenty "mismatches" under a primer carrying none.
    absent_kept = {name: runs for name, runs in absent_by_region.items() if name in kept}
    warn_about_absence(sum(1 for name in resting_on_absence if name in kept),
                       sum(1 for name in absence_beside_differences if name in kept), len(exclusion))
    final = settings.output / FINAL_NAME
    count = seqio.write_fasta(final, records)
    return final, count, absent_kept


def lower_positions(seq: str, positions: list[int]) -> str:
    """The sequence in upper case, with the given positions in lower case."""
    bases = list(seq.upper())
    for position in positions:
        if 0 <= position < len(bases):
            bases[position] = bases[position].lower()
    return "".join(bases)


def keep_previous_log(path: Path) -> None:
    """Move an earlier run's log aside before this run's log replaces it, one generation.

    `add_log_file` opens the log for writing, which empties it, so until this the record of whatever ran
    in this folder before was destroyed by the act of starting -- including by a run that then went no
    further, since `replace_previous_results` can refuse one. The results are kept one generation
    (`PREVIOUS_SUFFIX`) and the log was not, which left the previous `final_kmers.fasta.previous` with no
    log to say what produced it.
    """
    if not (path.exists() or path.is_symlink()):
        return
    aside = path.with_name(path.name + PREVIOUS_SUFFIX)
    try:
        path.replace(aside)
    except OSError as error:  # A directory in the way, or a folder that cannot be written to
        raise PrimerFinderError(
            f"Cannot move the earlier run's {path.name} aside to {aside}: {error}"
        ) from error


def replace_previous_results(settings: Settings) -> None:
    """Remove the answer an earlier run left in this folder, before producing a new one.

    A run is not resumable: a repeated command starts again from the kmers. But it can also fail part way
    through, and `run_info.json` is only written by a run that reaches the end -- so the folder was left
    holding the *previous* run's `final_kmers.fasta` with no record of which run it came from, and
    `primer-finder design` pointed at it would happily design on those stale regions. Both are removed
    here, at the point where this run commits to replacing them.
    """
    stale = [settings.output / FINAL_NAME, settings.output / RUN_INFO_NAME]
    # `exists()` is False for a dangling symlink, which would then survive this and break the final write
    found = [path for path in stale if path.exists() or path.is_symlink()]
    if found:
        log.warning("An earlier run's results are in %s (%s). This run replaces them; the old ones are "
                    "kept with a %s suffix, one generation only, its log included. A run is not "
                    "resumable: this one starts again from the kmers", settings.output,
                    ", ".join(path.name for path in found), PREVIOUS_SUFFIX)
    for path in found:
        # Moved aside rather than deleted. What this has to prevent is a *failed* re-run leaving an older
        # answer behind for `design` to use as though it were this run's; it does not have to destroy
        # the answer to do that. Only one generation is kept: a third run overwrites what the second
        # moved aside. The log goes the same way, in `keep_previous_log`, before it is reopened.
        aside = path.with_name(path.name + PREVIOUS_SUFFIX)
        try:
            path.replace(aside)
        except OSError as error:  # A directory in the way, or a folder that cannot be written to
            raise PrimerFinderError(
                f"Cannot move the earlier run's {path.name} aside to {aside}: {error}"
            ) from error


def clean_up(settings: Settings) -> None:
    """Remove what was only needed while running, unless the run asked to keep it."""
    if settings.keep_intermediate:
        return
    folder = settings.output
    shutil.rmtree(folder / KMER_FOLDER / "kmc_work", ignore_errors=True)
    for name in ("inclusion", "exclusion", "inclusion_specific"):
        for suffix in (".kmc_pre", ".kmc_suf"):
            (folder / KMER_FOLDER / (name + suffix)).unlink(missing_ok=True)
    for name in ("kmers.txt", "kmers_histogram.txt", "inclusion_kmers.txt"):
        (folder / KMER_FOLDER / name).unlink(missing_ok=True)
    (folder / CANDIDATE_FOLDER / "mapping.sam").unlink(missing_ok=True)
    shutil.rmtree(folder / BLAST_FOLDER / "inclusion_db", ignore_errors=True)
    shutil.rmtree(folder / BLAST_FOLDER / "exclusion_db", ignore_errors=True)
