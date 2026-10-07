"""Assembling the inclusion-specific kmers into contigs, with SKESA or SPAdes."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from primer_finder import PrimerFinderError, tools

log = logging.getLogger(__name__)

ASSEMBLERS = ("skesa", "spades")
PROGRAMS = {"skesa": "skesa", "spades": "spades.py"}


def assemble(kmer_fasta: Path, output: Path, assembler: str, threads: int, memory_gb: int) -> Path:
    """Assemble a fasta file of kmers into `output`. Returns `output`."""
    if assembler == "skesa":
        _skesa(kmer_fasta, output, threads, memory_gb)
    elif assembler == "spades":
        _spades(kmer_fasta, output, threads, memory_gb)
    else:  # pragma: no cover - the command line only accepts the two
        raise PrimerFinderError(f'Unknown assembler "{assembler}". Choose one of: {", ".join(ASSEMBLERS)}')
    if not output.exists() or output.stat().st_size == 0:
        raise PrimerFinderError(
            f"{assembler} could not assemble the inclusion-specific kmers into contigs. "
            "There may be too few of them; a smaller kmer size (-k) may help."
        )
    return output


def _skesa(kmer_fasta: Path, output: Path, threads: int, memory_gb: int) -> None:
    tools.run([
        "skesa",
        "--cores", str(threads),
        "--mem", str(memory_gb),
        "--fasta", kmer_fasta,
        "--contigs_out", output,
    ])


def _spades(kmer_fasta: Path, output: Path, threads: int, memory_gb: int) -> None:
    work_dir = output.parent / "spades"
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
        raise PrimerFinderError(f"SPAdes wrote no contigs ({contigs} missing)")
    shutil.move(str(contigs), output)
    shutil.rmtree(work_dir, ignore_errors=True)
