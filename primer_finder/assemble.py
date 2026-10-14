"""Assembling the inclusion-specific kmers into contigs, with SKESA or SPAdes."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from primer_finder import PrimerFinderError, seqio, tools

log = logging.getLogger(__name__)

ASSEMBLERS = ("skesa", "spades")
PROGRAMS = {"skesa": "skesa", "spades": "spades.py"}

# Where SPAdes stops, measured rather than read off a manual: assembled from every 99-mer of a synthetic
# stretch, stretches of 150, 180 and 199 bases are already absent from its simplified assembly graph, while
# 210 and longer come through. The floor is in graph simplification, not in what it chooses to write, so
# unlike SKESA's --min_contig it cannot be lowered by asking.
SPADES_FLOOR = 200


def assemble(kmer_fasta: Path, output: Path, assembler: str, threads: int, memory_gb: int,
             min_contig: int = 0) -> Path:
    """Assemble a fasta file of kmers into `output`. Returns `output`.

    `min_contig` is the shortest contig to report. Left alone both assemblers stop at 200 bases, which is
    longer than a kmer: an inclusion-specific stretch shorter than that is then never reported, and nothing
    downstream can see it. Those stretches are short exactly where the two groups are closely related, so
    the caller sets this to the kmer size. 0 leaves each assembler's own default alone.

    **It only lowers the floor for SKESA.** `--min_contig` tells SKESA what to report, so asking for less
    gets less. SPAdes has no such option and its floor is in graph simplification (`SPADES_FLOOR`), so a
    stretch shorter than that never reaches its output and asking cannot change it; here the value is only
    applied to what SPAdes did write, which keeps the two assemblers from reporting *longer* contigs than
    asked for. `-a spades` with a floor below `SPADES_FLOOR` is warned about rather than silently ignored.
    """
    if assembler == "skesa":
        _skesa(kmer_fasta, output, threads, memory_gb, min_contig)
    elif assembler == "spades":
        if 0 < min_contig < SPADES_FLOOR:
            log.warning("SPAdes does not report a contig shorter than about %d bases and has no option "
                        "for it, so regions down to %d will be missing. Use the default assembler to get "
                        "them", SPADES_FLOOR, min_contig)
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
        # SPAdes has no option for this. Applying it to the output cannot bring a short contig back -- it
        # was lost in graph simplification -- but it does keep `-a spades` from reporting contigs the
        # command was asked to leave out
        keep = [record for record in seqio.read_fasta(output).values() if len(record.seq) >= min_contig]
        seqio.write_fasta(output, keep)
