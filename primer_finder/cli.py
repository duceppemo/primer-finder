"""The command line: `primer-finder` runs the pipeline, `primer-finder idt` converts an IDT order sheet."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from primer_finder import PrimerFinderError, __version__, assemble, design, idt, kmers
from primer_finder.pipeline import Settings, run
from primer_finder.system import default_memory_gb, usable_cpus

log = logging.getLogger(__name__)

COMMANDS = ("find", "design", "idt")
DESCRIPTION = (
    "Find group-specific kmers to design selective qPCR assays: kmers shared by every inclusion genome "
    "and absent from every exclusion genome."
)


def _fraction(value: str) -> float:
    """A fraction greater than 0 and at most 1, as -p/--min-inclusion takes."""
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not a number") from None
    if not 0 < number <= 1:
        raise argparse.ArgumentTypeError(f"{value} is not a fraction greater than 0 and at most 1")
    return number


def build_parser(max_cpu: int, max_mem: int) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="primer-finder", description=DESCRIPTION)
    parser.add_argument("-v", "--version", action="version", version=f"primer-finder {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="{find,design,idt}")

    find = commands.add_parser(
        "find", help="Find inclusion-specific contigs (the default command).", description=DESCRIPTION
    )
    find.add_argument("-i", "--inclusion", metavar="/inclusion_folder/", required=True, type=Path,
                      help="Folder that holds the genomes the assay should amplify.")
    find.add_argument("-e", "--exclusion", metavar="/exclusion_folder/", required=True, type=Path,
                      help="Folder that holds the genomes the assay should not amplify.")
    find.add_argument("-o", "--output", metavar="/output_folder/", required=True, type=Path,
                      help="Folder to hold the result files.")
    find.add_argument("-t", "--threads", metavar=str(max_cpu), type=int, default=max_cpu,
                      help=f"Number of CPU. Default is every CPU available ({max_cpu}).")
    find.add_argument("-m", "--memory", metavar=str(max_mem), type=int, default=max_mem,
                      help=f"Memory in GB. Default is 85%% of the total memory ({max_mem}).")
    find.add_argument("-k", "--kmer_size", "--kmer-size", metavar="99", type=int, default=99,
                      help=f"Kmer size for KMC ({kmers.MIN_KMER_SIZE}-{kmers.MAX_KMER_SIZE}). Default 99.")
    find.add_argument("-d", "--duplication", metavar="1", type=int, default=1,
                      help="Maximum number of times a kmer can be found in each inclusion genome. Default 1, "
                           "meaning that repeated regions are discarded.")
    find.add_argument("-p", "--min-inclusion", metavar="1.0", type=_fraction, default=1.0,
                      help="Fraction of the inclusion genomes a kmer must be in, between 0 and 1. Default "
                           "1.0, meaning every one of them. 0.9 keeps what 90%% of them share, rounded up, "
                           "and the assay will not amplify the genomes that lack it.")
    find.add_argument("-r", "--reference", metavar="/exclusion_folder/genome.fasta", type=Path,
                      help="Exclusion genome to map the assembled kmers to. Default is the first genome of "
                           "the exclusion folder in alphabetical order.")
    find.add_argument("-a", "--assembler", metavar="skesa", default="skesa", choices=assemble.ASSEMBLERS,
                      help='Assembler for the kmers: "skesa" or "spades". Default skesa.')
    find.add_argument("--keep-intermediate", action="store_true",
                      help="Keep the KMC databases, the SAM file and the blast databases.")
    find.add_argument("--debug", action="store_true", help="Verbose logging.")
    find.add_argument("-v", "--version", action="version", version=f"primer-finder {__version__}")

    design_parser = commands.add_parser(
        "design", help="Design assays on the regions of a finished run, with Primer3.",
        description="Design primers and probes on the candidate regions of a finished primer-finder run, "
                    "keep the assays that could tell the two groups apart, and write the primer files "
                    "insilicoPCR reads.",
    )
    design_parser.add_argument("results", type=Path, help="The output folder of a finished run.")
    design_parser.add_argument("-o", "--output", metavar="/assay_folder/", required=True, type=Path,
                               help="Folder to hold the designed assays.")
    design_parser.add_argument("-i", "--inclusion", metavar="/inclusion_folder/", type=Path,
                               help="Inclusion genomes. Default: what the run recorded.")
    design_parser.add_argument("-e", "--exclusion", metavar="/exclusion_folder/", type=Path,
                               help="Exclusion genomes. Default: what the run recorded.")
    design_parser.add_argument("-t", "--threads", metavar=str(max_cpu), type=int, default=max_cpu,
                               help=f"Number of CPU. Default is every CPU available ({max_cpu}).")
    design_parser.add_argument("--product-size", metavar=design.DEFAULT_PRODUCT_SIZE,
                               default=design.DEFAULT_PRODUCT_SIZE,
                               help=f"Amplicon size range for Primer3. Default {design.DEFAULT_PRODUCT_SIZE}.")
    design_parser.add_argument("--assays-per-region", metavar=str(design.DEFAULT_ASSAYS_PER_REGION),
                               type=int, default=design.DEFAULT_ASSAYS_PER_REGION,
                               help="How many assays Primer3 should propose per region. Default "
                                    f"{design.DEFAULT_ASSAYS_PER_REGION}.")
    design_parser.add_argument("--max-regions", metavar=str(design.DEFAULT_MAX_REGIONS), type=int,
                               default=design.DEFAULT_MAX_REGIONS,
                               help="Design on this many regions, the most promising first. 0 for every "
                                    f"region. Default {design.DEFAULT_MAX_REGIONS}.")
    design_parser.add_argument("--insilico-pcr", metavar="/insilicoPCR/", type=Path,
                               help="Run in silico PCR of the designed assays against both groups with "
                                    "insilicoPCR, and report which assays amplify every inclusion genome "
                                    "and no exclusion genome. Give the folder of an extracted portable "
                                    "release, its jar, or a launcher script.")
    design_parser.add_argument("--mismatches", metavar="0", type=int, default=0,
                               help="Primer mismatches insilicoPCR should allow. Default 0.")
    design_parser.add_argument("--debug", action="store_true", help="Verbose logging.")

    convert = commands.add_parser(
        "idt", help="Convert an IDT order sheet into a fasta file of assays.",
        description="Convert an IDT order sheet (.xlsx, .csv or .tsv) into a fasta file holding one record "
                    "per oligo, named after its assay and amplicon size.",
    )
    convert.add_argument("table", type=Path, help="The IDT sheet (.xlsx, .csv or .tsv).")
    convert.add_argument("output", type=Path, help="The fasta file to write.")
    convert.add_argument("prefix", nargs="?", default="", help="Optional prefix for the assay names.")
    convert.add_argument("--debug", action="store_true", help="Verbose logging.")
    return parser


def with_default_command(argv: list[str]) -> list[str]:
    """`primer-finder -i ... -e ... -o ...` runs the `find` command, as the previous versions did."""
    if argv and argv[0] not in COMMANDS and not argv[0].startswith("-"):
        return argv  # An unknown command: let argparse report it
    if argv and argv[0] not in COMMANDS and argv[0] not in ("-h", "--help", "-v", "--version"):
        return ["find", *argv]
    return argv


def setup_logging(debug: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )


def main(argv: list[str] | None = None) -> int:
    given = list(sys.argv[1:] if argv is None else argv)
    max_cpu, max_mem = usable_cpus(), default_memory_gb()
    parser = build_parser(max_cpu, max_mem)
    args = parser.parse_args(with_default_command(given))
    if args.command is None:
        parser.print_help()
        return 1
    setup_logging(args.debug)

    try:
        if args.command == "idt":
            idt.convert(args.table, args.output, args.prefix)
            return 0
        if args.command == "design":
            return design.run(design.DesignSettings(
                results=args.results,
                output=args.output,
                threads=args.threads,
                inclusion=args.inclusion,
                exclusion=args.exclusion,
                product_size=args.product_size,
                assays_per_region=args.assays_per_region,
                max_regions=args.max_regions,
                insilico_pcr=args.insilico_pcr,
                mismatches=args.mismatches,
                command_line=["primer-finder", *given],
            ))
        threads = args.threads
        if threads < 1:
            parser.error("-t/--threads must be 1 or more")
        if threads > max_cpu:
            log.warning("Asked for %d threads but only %d CPU(s) are available; using %d",
                        threads, max_cpu, max_cpu)
            threads = max_cpu
        memory = args.memory
        if memory < 1:
            parser.error("-m/--memory must be 1 GB or more")
        if memory > max_mem:
            log.warning("Asked for %d GB but only about %d GB are available; using %d GB",
                        memory, max_mem, max_mem)
            memory = max_mem
        settings = Settings(
            inclusion=args.inclusion,
            exclusion=args.exclusion,
            output=args.output,
            threads=threads,
            memory_gb=memory,
            kmer_size=args.kmer_size,
            duplication=args.duplication,
            min_inclusion=args.min_inclusion,
            reference=args.reference,
            assembler=args.assembler,
            keep_intermediate=args.keep_intermediate,
            command_line=["primer-finder", *given],
        )
        return run(settings)
    except PrimerFinderError as exc:
        log.error("%s", exc)
        return 1
    except OSError as exc:  # An unreadable input, a read-only or missing output folder, a full disk
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        log.error("Interrupted")
        return 130
