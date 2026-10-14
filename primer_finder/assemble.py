"""Assembling the inclusion-specific kmers into contigs, with SKESA or SPAdes."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from primer_finder import PrimerFinderError, seqio, tools

log = logging.getLogger(__name__)

ASSEMBLERS = ("skesa", "spades")
PROGRAMS = {"skesa": "skesa", "spades": "spades.py"}


def assemble(kmer_fasta: Path, output: Path, assembler: str, threads: int, memory_gb: int,
             min_contig: int = 0) -> Path:
    """Assemble a fasta file of kmers into `output`. Returns `output`.

    `min_contig` is the shortest contig to report. Both assemblers stop at 200 bases of their own accord,
    which is longer than a kmer: an inclusion-specific stretch shorter than that is then never reported, and
    nothing downstream can see it. The stretches are short exactly where the two groups are closely
    related, so the caller sets this to the kmer size and lets a single kmer's worth of sequence through.
    0 leaves each assembler's own default alone.
    """
    if assembler == "skesa":
        _skesa(kmer_fasta, output, threads, memory_gb, min_contig)
    elif assembler == "spades":
        _spades(kmer_fasta, output, threads, memory_gb, min_contig)
    else:  # pragma: no cover - the command line only accepts the two
        raise PrimerFinderError(f'Unknown assembler "{assembler}". Choose one of: {", ".join(ASSEMBLERS)}')
    if not output.exists() or output.stat().st_size == 0:
        raise PrimerFinderError(
            f"{assembler} could not assemble the inclusion-specific kmers into contigs. "
            "There may be too few of them; a smaller kmer size (-k) may help."
        )
    return output


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
    if min_contig:
        # SPAdes has no option for this, so the same floor is applied here instead: whatever the two
        # assemblers do of their own accord, the command reports the same lengths with either
        keep = [record for record in seqio.read_fasta(output).values() if len(record.seq) >= min_contig]
        seqio.write_fasta(output, keep)
