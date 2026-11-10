#!/usr/bin/env python3
"""Genomes whose contigs break inside primer sites, and assays in the files primer-finder writes.

insilicoPCR 0.11.0 changed one thing in assembly mode: the size of a product whose primer hangs off the
start of a contig. This builds the case on purpose, on real sequence. From the first 600 kb of *Listeria
monocytogenes* EGD-e it takes 300 assays (a 20-base primer pair 80-150 bases apart, and a 24-base probe
between them) and writes three genomes:

- `lmono_cut`: the region, broken into contigs. For one assay in six the break sits 1-7 bases inside the
  forward primer, so the contig starts with the primer's 5' end missing; for another one in six it sits
  1-7 bases before the 5' end of the reverse primer, so the contig ends inside that primer. insilicoPCR
  lets up to 5 unaligned 5' bases through (`--5p-truncation`), so 1-5 should amplify and 6-7 should not.
- `lmono_mutated`: the same region with one base in 45 substituted and half of the breaks, for primers
  that bind through a mismatch, at the 3' end or elsewhere -- what the lenient 3'-end pass is about.
- `linnocua`: the first 600 kb of *L. innocua* Clip11262, a near neighbour.

The region is 600 kb and not the genome because insilicoPCR blasts each genome against the primers, and
with 1,500 assays on the whole genomes one pass had not finished after ten minutes.

    python make_contig_edges.py <lmono_EGDe.fasta> <linnocua_Clip11262.fasta> <output_folder>

The genomes are GCF_000196035.1 and GCF_000195795.1. The seed is fixed, so the output is the same every
time. Writes genomes/, assays_pcr.fasta, assays_qpcr.fasta and plan.tsv (the size each assay was built
with, and where its site was broken).
"""
import random
import sys
from pathlib import Path

REGION = 600_000
ASSAYS = 300
SEED = 11


def read(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip().upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGT", "TGCA"))


def wrap(sequence: str) -> str:
    return "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80))


def contigs(sequence: str, breaks: list[int], prefix: str) -> str:
    edges = [0, *sorted(breaks), len(sequence)]
    return "".join(f">{prefix}_{number + 1}\n{wrap(sequence[start:end])}\n"
                   for number, (start, end) in enumerate(zip(edges, edges[1:], strict=False)))


def main(lmono: Path, linnocua: Path, out: Path) -> None:
    rng = random.Random(SEED)
    mono = read(lmono)[0][1][:REGION]
    innocua_name, innocua = read(linnocua)[0]
    (out / "genomes").mkdir(parents=True, exist_ok=True)

    breaks: list[int] = []
    pcr: list[str] = []
    qpcr: list[str] = []
    plan: list[str] = []
    step = (len(mono) - 4000) // ASSAYS
    for i in range(ASSAYS):
        start = 1000 + i * step + rng.randrange(200)
        size = rng.randrange(80, 151)
        forward = mono[start:start + 20]
        reverse = reverse_complement(mono[start + size - 20:start + size])
        probe = mono[start + 28:start + 52]
        if set(forward + reverse + probe) - set("ACGT"):
            continue
        name = f"lm_{start}_{i}"
        pcr += [f">{name}-F\n{forward}", f">{name}-R\n{reverse}"]
        qpcr += [f">{name}-F\n{forward}", f">{name}-R\n{reverse}", f">{name}-P\n{probe}"]
        kind, bases = "whole", 0
        if i % 6 == 1:  # The contig starts inside the forward primer: its 5' end hangs off
            bases = 1 + (i // 6) % 7
            breaks.append(start + bases)
            kind = "forward_off_start"
        elif i % 6 == 3:  # The contig ends inside the reverse primer: its 5' end hangs off
            bases = 1 + (i // 6) % 7
            breaks.append(start + size - bases)
            kind = "reverse_off_end"
        plan.append(f"{name}\t{size}\t{kind}\t{bases}")

    (out / "assays_pcr.fasta").write_text("\n".join(pcr) + "\n")
    (out / "assays_qpcr.fasta").write_text("\n".join(qpcr) + "\n")
    (out / "plan.tsv").write_text("assay\tsize\tsite\tbases_off\n" + "\n".join(plan) + "\n")

    (out / "genomes" / "lmono_cut.fasta").write_text(contigs(mono, breaks, "lmono_contig"))
    mutated = list(mono)
    for position, base in enumerate(mutated):
        if rng.random() < 1 / 45:
            mutated[position] = rng.choice([other for other in "ACGT" if other != base])
    (out / "genomes" / "lmono_mutated.fasta").write_text(
        contigs("".join(mutated), breaks[::2], "lmono_mutated_contig"))
    (out / "genomes" / "linnocua.fasta").write_text(f">{innocua_name}\n{wrap(innocua[:REGION])}\n")
    print(f"{len(plan)} assays, {len(breaks)} sites broken")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
