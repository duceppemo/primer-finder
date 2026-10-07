"""The pipeline: count kmers, assemble them, then check the contigs against every genome."""

from __future__ import annotations

import json
import logging
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

from primer_finder import PrimerFinderError, __version__, assemble, blast, kmers, mapping, seqio, tools
from primer_finder.mapping import PRIMER_LENGTH

log = logging.getLogger(__name__)

LOG_NAME = "primer_finder.log"
RUN_INFO_NAME = "run_info.json"
FINAL_NAME = "final_kmers.fasta"
BEST_NAME = "best_kmers.fasta"
ALL_INCLUSION_NAME = "all_inclusion_contigs.fasta"
HITS_NAME = "inclusion_blast_hits.tsv"

KMER_FOLDER = "1_kmers"
ASSEMBLY_FOLDER = "2_assembly"
CANDIDATE_FOLDER = "3_candidates"
BLAST_FOLDER = "4_blast"

BASE_PROGRAMS = ("kmc", "kmc_tools", "minimap2", "makeblastdb", "blastn")


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
    reference: Path | None = None
    assembler: str = "skesa"
    keep_intermediate: bool = False
    command_line: list[str] = field(default_factory=list)


def run(settings: Settings) -> int:
    """Run the whole pipeline. Returns 0, or raises PrimerFinderError."""
    started = time.time()
    inclusion, exclusion, reference = check(settings)
    settings.output.mkdir(parents=True, exist_ok=True)
    add_log_file(settings.output / LOG_NAME)
    log.info("primer-finder %s", __version__)
    log.info("%d inclusion genome(s), %d exclusion genome(s)", len(inclusion), len(exclusion))

    info: dict[str, object] = {
        "version": __version__,
        "command_line": settings.command_line,
        "parameters": {
            "inclusion": str(settings.inclusion),
            "exclusion": str(settings.exclusion),
            "output": str(settings.output),
            "kmer_size": settings.kmer_size,
            "duplication": settings.duplication,
            "assembler": settings.assembler,
            "threads": settings.threads,
            "memory_gb": settings.memory_gb,
        },
        "inclusion_genomes": [str(path) for path in inclusion],
        "exclusion_genomes": [str(path) for path in exclusion],
        "reference": str(reference),
        "programs": tools.versions(sorted({*BASE_PROGRAMS, assemble.PROGRAMS[settings.assembler]})),
        "counts": {},
    }
    counts: dict[str, int] = info["counts"]  # type: ignore[assignment]

    kmer_fasta, counts["kmers"] = find_specific_kmers(settings, inclusion, exclusion)
    assembly = assemble_kmers(settings, kmer_fasta)
    counts["contigs"] = seqio.count_records(assembly)

    best, counts["candidates"] = find_candidates(settings, assembly, reference)
    all_inclusion, counts["in_all_inclusion"] = keep_shared_by_inclusion(settings, best, inclusion)
    final, counts["final"] = keep_absent_from_exclusion(settings, all_inclusion, exclusion)

    info["seconds"] = round(time.time() - started, 1)
    (settings.output / RUN_INFO_NAME).write_text(json.dumps(info, indent=2) + "\n")
    clean_up(settings)
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
    if settings.output.resolve() in {path.resolve().parent for path in inclusion + exclusion}:
        raise PrimerFinderError(f"The output folder holds input genomes: {settings.output}")

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

    tools.require(sorted({*BASE_PROGRAMS, assemble.PROGRAMS[settings.assembler]}))
    return inclusion, exclusion, reference


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

    log.info("Counting the %d-mers shared by the %d inclusion genomes...", settings.kmer_size, len(inclusion))
    inclusion_db = kmers.count(
        kmers.write_file_list(inclusion, folder / "inclusion_list.txt"),
        folder / "inclusion", work_dir, settings.kmer_size, settings.threads, settings.memory_gb,
        min_count=len(inclusion), max_count=len(inclusion) * settings.duplication,
    )
    log.info("Counting the %d-mers of the %d exclusion genomes...", settings.kmer_size, len(exclusion))
    exclusion_db = kmers.count(
        kmers.write_file_list(exclusion, folder / "exclusion_list.txt"),
        folder / "exclusion", work_dir, settings.kmer_size, settings.threads, settings.memory_gb,
        min_count=1, max_count=kmers.NO_MAX_COUNT,
    )
    log.info("Subtracting the exclusion kmers from the inclusion ones...")
    specific_db = kmers.subtract(inclusion_db, exclusion_db, folder / "inclusion_specific", settings.threads)

    dump_file = kmers.dump(specific_db, folder / "kmers.txt", settings.threads)
    kmer_fasta = folder / f"inclusion_specific_{settings.kmer_size}mers.fasta"
    count = kmers.dump_to_fasta(dump_file, kmer_fasta)
    log.info("Found %d inclusion-specific %d-mers", count, settings.kmer_size)
    return kmer_fasta, count


def assemble_kmers(settings: Settings, kmer_fasta: Path) -> Path:
    """Assemble the inclusion-specific kmers, to have fewer and longer sequences to check."""
    folder = settings.output / ASSEMBLY_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    assembly = folder / "assembly.fasta"
    log.info("Assembling the kmers with %s...", settings.assembler)
    assemble.assemble(kmer_fasta, assembly, settings.assembler, settings.threads, settings.memory_gb)
    log.info("Assembled %d contig(s)", seqio.count_records(assembly))
    return assembly


def find_candidates(settings: Settings, assembly: Path, reference: Path) -> tuple[Path, int]:
    """Map the contigs to one exclusion genome and keep those whose differences could fit in a primer."""
    folder = settings.output / CANDIDATE_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    log.info("Mapping the contigs to the exclusion genome %s...", reference.name)
    sam_file = mapping.map_contigs(assembly, reference, folder / "mapping.sam", settings.threads)

    log.info("Keeping the contigs with at least two differences within %d bases...", PRIMER_LENGTH)
    candidates = mapping.select_candidates(seqio.read_fasta(assembly), mapping.parse_sam(sam_file))
    best = folder / BEST_NAME
    count = seqio.write_fasta(best, (candidate.as_record() for candidate in mapping.sort_candidates(candidates)))
    log.info("%d contig(s) could carry a selective assay", count)
    if count == 0:
        raise PrimerFinderError(
            "No contig differs enough from the exclusion genome to carry a selective assay "
            f"(at least two differences within {PRIMER_LENGTH} bases)."
        )
    return best, count


def keep_shared_by_inclusion(settings: Settings, best: Path, inclusion: list[Path]) -> tuple[Path, int]:
    """Drop the candidates that are missing from any inclusion genome: the mapping step only compared them
    with one exclusion genome, and the assembly may hold kmers from several genomes."""
    folder = settings.output / BLAST_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    log.info("Checking the candidates against the %d inclusion genomes...", len(inclusion))
    presence = blast.presence_in_genomes(best, inclusion, folder / "inclusion_db", settings.threads)
    write_presence_table(folder / HITS_NAME, presence)

    all_inclusion = folder / ALL_INCLUSION_NAME
    count = seqio.write_fasta(all_inclusion, (
        record for record in seqio.iter_records(best) if all(presence[record.name].values())
    ))
    log.info("%d contig(s) are present in all inclusion genomes", count)
    if count == 0:
        raise PrimerFinderError(
            "No candidate contig is present in all the inclusion genomes. The presence of each one in each "
            f"genome is in {folder / HITS_NAME}."
        )
    return all_inclusion, count


def write_presence_table(path: Path, presence: dict[str, dict[str, bool]]) -> None:
    """A table of contigs (rows) against genomes (columns), 1 when the contig is present."""
    genomes = sorted({genome for hits in presence.values() for genome in hits})
    with path.open("w") as fh:
        fh.write("contig\t" + "\t".join(genomes) + "\n")
        for contig in sorted(presence):
            row = ["1" if presence[contig].get(genome) else "0" for genome in genomes]
            fh.write(f"{contig}\t" + "\t".join(row) + "\n")


def keep_absent_from_exclusion(
    settings: Settings, all_inclusion: Path, exclusion: list[Path]
) -> tuple[Path, int]:
    """Keep the contigs whose differences hold up against every exclusion genome, not just the one they were
    mapped to. A contig that no exclusion genome hits at all is kept as it is."""
    folder = settings.output / BLAST_FOLDER
    log.info("Checking the differences against the %d exclusion genomes...", len(exclusion))
    results = blast.exclusion_variants(all_inclusion, exclusion, folder / "exclusion_db", settings.threads)
    candidates = seqio.read_fasta(all_inclusion)

    records: list[seqio.Record] = []
    for name, record in candidates.items():
        result = results.get(name, blast.ExclusionResult())
        if result.genomes_hit == 0:  # Absent from every exclusion genome: specific on its own
            records.append(record)
            continue
        positions = blast.shared_variants(result)
        if len(positions) > 1 and blast.has_close_variants(positions):
            records.append(seqio.Record(name, str(positions), lower_positions(record.seq, positions)))
    final = settings.output / FINAL_NAME
    count = seqio.write_fasta(final, records)
    return final, count


def lower_positions(seq: str, positions: list[int]) -> str:
    """The sequence in upper case, with the given positions in lower case."""
    bases = list(seq.upper())
    for position in positions:
        if 0 <= position < len(bases):
            bases[position] = bases[position].lower()
    return "".join(bases)


def clean_up(settings: Settings) -> None:
    """Remove what was only needed while running, unless the run asked to keep it."""
    if settings.keep_intermediate:
        return
    folder = settings.output
    shutil.rmtree(folder / KMER_FOLDER / "kmc_work", ignore_errors=True)
    for name in ("inclusion", "exclusion", "inclusion_specific"):
        for suffix in (".kmc_pre", ".kmc_suf"):
            (folder / KMER_FOLDER / (name + suffix)).unlink(missing_ok=True)
    (folder / KMER_FOLDER / "kmers.txt").unlink(missing_ok=True)
    (folder / CANDIDATE_FOLDER / "mapping.sam").unlink(missing_ok=True)
    shutil.rmtree(folder / BLAST_FOLDER / "inclusion_db", ignore_errors=True)
    shutil.rmtree(folder / BLAST_FOLDER / "exclusion_db", ignore_errors=True)
