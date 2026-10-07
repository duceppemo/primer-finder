"""The command line: `primer-finder` runs the pipeline, `primer-finder idt` converts an IDT order sheet."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from primer_finder import PrimerFinderError, __version__, assemble, idt, kmers
from primer_finder.pipeline import Settings, run
from primer_finder.system import default_memory_gb, usable_cpus

log = logging.getLogger(__name__)

COMMANDS = ("find", "idt")
DESCRIPTION = (
    "Find group-specific kmers to design selective qPCR assays: kmers shared by every inclusion genome "
    "and absent from every exclusion genome."
)


def build_parser(max_cpu: int, max_mem: int) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="primer-finder", description=DESCRIPTION)
    parser.add_argument("-v", "--version", action="version", version=f"primer-finder {__version__}")
    commands = parser.add_subparsers(dest="command", metavar="{find,idt}")

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
    find.add_argument("-r", "--reference", metavar="/exclusion_folder/genome.fasta", type=Path,
                      help="Exclusion genome to map the assembled kmers to. Default is the first genome of "
                           "the exclusion folder in alphabetical order.")
    find.add_argument("-a", "--assembler", metavar="skesa", default="skesa", choices=assemble.ASSEMBLERS,
                      help='Assembler for the kmers: "skesa" or "spades". Default skesa.')
    find.add_argument("--keep-intermediate", action="store_true",
                      help="Keep the KMC databases, the SAM file and the blast databases.")
    find.add_argument("--debug", action="store_true", help="Verbose logging.")
    find.add_argument("-v", "--version", action="version", version=f"primer-finder {__version__}")

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
            reference=args.reference,
            assembler=args.assembler,
            keep_intermediate=args.keep_intermediate,
            command_line=["primer-finder", *given],
        )
        return run(settings)
    except PrimerFinderError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        log.error("Interrupted")
        return 130
