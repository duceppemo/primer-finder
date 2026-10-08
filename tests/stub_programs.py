"""Stand-ins for the external programs, so that the pipeline can be tested without installing them.

Each stub writes the files the real program would write. The `stubs` fixture puts a two-line launcher for
every program on PATH; what they produce is read from the PRIMER_FINDER_STUBS environment variable (JSON):

    kmers            the kmers "kmc_tools transform ... dump" writes (an empty list writes an empty file)
    contigs          {name: sequence} for the assembler to write
    sam              [[query name, flag, cigar], ...] for minimap2 to write
    presence         {genome: [contig, ...]} the contigs blastn finds in each inclusion genome
    exclusion_hits   {genome: [[query, qstart, qend, evalue, qseq, sseq], ...]} for the exclusion genomes
    fail             the name of a program that should exit 1

Every call is appended to the file named by PRIMER_FINDER_CALLS, one command line per line.
"""

from __future__ import annotations

import json
import os
import shlex
import sys
from pathlib import Path

DEFAULT_KMERS = ["ACGT" * 5, "TTGA" * 5, "GGCA" * 5]
DEFAULT_CONTIGS = {"ctg1": "ACGT" * 25}
# A contig that differs from the exclusion genome by two mismatches 7 bases apart: a candidate.
DEFAULT_SAM = [["ctg1", 0, "40=1X7=1X50="]]
DEFAULT_QSEQ = "ACGTACGTACGTACGTACGT"
DEFAULT_SSEQ = "ACGTAAGTACGTTCGTACGT"  # Mismatches at offsets 5 and 12


def scenario() -> dict:
    return json.loads(os.environ.get("PRIMER_FINDER_STUBS", "{}"))


def value(key: str, default):
    return scenario().get(key, default)


def argument(argv: list[str], flag: str, default: str = "") -> str:
    """The value that follows `flag`."""
    return argv[argv.index(flag) + 1] if flag in argv else default


def kmc(argv: list[str]) -> int:
    """kmc -k.. -t.. -m.. -fm -ci.. -cx.. @list <db> <work_dir>"""
    positional = [part for part in argv if not part.startswith(("-", "@"))]
    database = Path(positional[0])
    for suffix in (".kmc_pre", ".kmc_suf"):
        database.with_name(database.name + suffix).touch()
    return 0


def kmc_tools(argv: list[str]) -> int:
    """kmc_tools -t.. simple <a> <b> kmers_subtract <db>   |   kmc_tools -t.. transform <db> dump <file>"""
    if "kmers_subtract" in argv:
        database = Path(argv[argv.index("kmers_subtract") + 1])
        for suffix in (".kmc_pre", ".kmc_suf"):
            database.with_name(database.name + suffix).touch()
        return 0
    if "dump" in argv:
        out_file = Path(argv[argv.index("dump") + 1])
        kmers = value("kmers", DEFAULT_KMERS)
        out_file.write_text("".join(f"{kmer}\t2\n" for kmer in kmers))
        return 0
    return 1  # pragma: no cover


def _write_contigs(path: Path) -> None:
    contigs: dict[str, str] = value("contigs", DEFAULT_CONTIGS)
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in contigs.items()))


def skesa(argv: list[str]) -> int:
    _write_contigs(Path(argument(argv, "--contigs_out")))
    return 0


def spades(argv: list[str]) -> int:
    work_dir = Path(argument(argv, "-o"))
    work_dir.mkdir(parents=True, exist_ok=True)
    _write_contigs(work_dir / "contigs.fasta")
    return 0


def minimap2(argv: list[str]) -> int:
    """Write a SAM file on standard output, one record per entry of `sam`."""
    contigs: dict[str, str] = value("contigs", DEFAULT_CONTIGS)
    print("@HD\tVN:1.6\tSO:unsorted")
    for name, flag, cigar in value("sam", DEFAULT_SAM):
        sequence = contigs.get(name, "A" * 100)
        fields = [name, str(flag), "reference", "1", "60", cigar, "*", "0", "0", sequence, "*"]
        print("\t".join(fields))
    return 0


def makeblastdb(argv: list[str]) -> int:
    database = Path(argument(argv, "-out"))
    database.parent.mkdir(parents=True, exist_ok=True)
    database.with_name(database.name + ".nin").touch()
    return 0


def stub_programs_local_fasta() -> str:
    """The name make_db() links each genome to, next to the database."""
    from primer_finder.blast import LOCAL_FASTA

    return LOCAL_FASTA


def blastn(argv: list[str]) -> int:
    """Write the tabular output the pipeline asked for, for the genome the database was made from."""
    fields = argument(argv, "-outfmt").split()[1:]
    out_file = Path(argument(argv, "-out"))
    # blastn runs in the database's folder, with the genome linked there under a fixed name. A scenario can
    # name a genome by its file stem, by "<parent folder>/<stem>" when two genomes share a name, or by the
    # folder blast is running in.
    local = Path(stub_programs_local_fasta())
    real = local.resolve() if local.exists() else local
    keys = [real.stem, f"{real.parent.name}/{real.stem}", Path.cwd().name]
    genome = next((key for key in keys if key in value("presence", {}) or key in value("exclusion_hits", {})),
                  keys[0])
    queries = [line[1:].split()[0] for line in Path(argument(argv, "-query")).read_text().splitlines()
               if line.startswith(">")]
    rows: list[list[str]] = []
    if "length" in fields and "pident" in fields:
        # Counting the copies of an amplicon: one row per copy, each covering the whole query
        lengths, name = {}, None
        for line in Path(argument(argv, "-query")).read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name:
                lengths[name] += len(line.strip())
        copies = value("amplicon_hits", {}).get(genome, 0)
        for query in queries:
            rows.extend([[query, "1e-30", str(lengths.get(query, 100)), "100.0"]] * copies)
    elif "qseq" in fields:  # The exclusion step: aligned sequences
        hits = value("exclusion_hits", {})
        default = [[query, 1, len(DEFAULT_QSEQ), 1e-30, DEFAULT_QSEQ, DEFAULT_SSEQ] for query in queries]
        for query, qstart, qend, evalue, qseq, sseq in hits.get(genome, default):
            rows.append([str(query), str(qstart), str(qend), str(evalue), qseq, sseq])
    else:  # The inclusion step: is the contig there at all?
        present = value("presence", {}).get(genome, queries)
        evalue = value("presence_evalue", "1e-30")
        rows = [[query, evalue] for query in queries if query in present]
    out_file.write_text("".join("\t".join(row) + "\n" for row in rows))
    return 0


def primer3_core(argv: list[str]) -> int:
    """Answer a boulder-IO file the way Primer3 does, placing the oligos where it was told to.

    One pair per record, unless the scenario says otherwise:
      primer3_pairs   how many pairs to return (0 for "nothing could be designed")
      primer3_error   a message to put in PRIMER_ERROR
    A forced 3' end or probe junction is honoured, so that the caller's coordinate handling is tested.
    """
    records, current = [], {}
    for line in Path(argv[0]).read_text().splitlines():
        if line == "=":
            records.append(current)
            current = {}
        elif "=" in line:
            key, value = line.split("=", 1)
            current[key] = value

    pairs = value_of("primer3_pairs", 1)
    error = value_of("primer3_error", "")
    out = []
    for record in records:
        fields = [f"SEQUENCE_ID={record.get('SEQUENCE_ID', '')}"]
        template = record.get("SEQUENCE_TEMPLATE", "")
        if error:
            fields.append(f"PRIMER_ERROR={error}")
            out.append("\n".join(fields) + "\n=\n")
            continue
        fields.append(f"PRIMER_PAIR_NUM_RETURNED={pairs}")
        for number in range(pairs):
            size = 20
            if "SEQUENCE_FORCE_LEFT_END" in record:
                left = int(record["SEQUENCE_FORCE_LEFT_END"]) - size + 1
            else:
                left = 10 + number
            if "SEQUENCE_FORCE_RIGHT_END" in record:
                right_start = int(record["SEQUENCE_FORCE_RIGHT_END"])
            else:
                right_start = left + 80
            if "SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST" in record:
                probe_start = max(left + size, int(record["SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST"]) - 5)
            else:
                probe_start = left + size + 5
            right_end = right_start + size - 1  # Primer3 reports the rightmost base of a right primer
            fields += [
                f"PRIMER_LEFT_{number}={left},{size}",
                f"PRIMER_LEFT_{number}_SEQUENCE={template[left:left + size] or 'A' * size}",
                f"PRIMER_LEFT_{number}_TM=60.0",
                f"PRIMER_LEFT_{number}_GC_PERCENT=50.0",
                f"PRIMER_RIGHT_{number}={right_end},{size}",
                f"PRIMER_RIGHT_{number}_SEQUENCE={template[right_start:right_end + 1] or 'T' * size}",
                f"PRIMER_RIGHT_{number}_TM=60.0",
                f"PRIMER_RIGHT_{number}_GC_PERCENT=50.0",
                f"PRIMER_PAIR_{number}_PRODUCT_SIZE={right_end - left + 1}",
                f"PRIMER_PAIR_{number}_PENALTY={0.5 + number}",
            ]
            if value_of("primer3_probe", True):
                fields += [
                    f"PRIMER_INTERNAL_{number}={probe_start},{size}",
                    f"PRIMER_INTERNAL_{number}_SEQUENCE={template[probe_start:probe_start + size] or 'C' * size}",
                    f"PRIMER_INTERNAL_{number}_TM=68.0",
                    f"PRIMER_INTERNAL_{number}_GC_PERCENT=55.0",
                ]
        out.append("\n".join(fields) + "\n=\n")
    print("".join(out), end="")
    return 0


def value_of(key: str, default):
    return scenario().get(key, default)


PROGRAMS = {
    "kmc": kmc,
    "primer3_core": primer3_core,
    "kmc_tools": kmc_tools,
    "skesa": skesa,
    "spades.py": spades,
    "minimap2": minimap2,
    "makeblastdb": makeblastdb,
    "blastn": blastn,
}


def main(name: str, argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    calls = os.environ.get("PRIMER_FINDER_CALLS")
    if calls:
        with open(calls, "a") as fh:
            fh.write(f"{name} {shlex.join(argv)}\n")
    if argv in (["--version"], ["-version"]) or (not argv and name in ("kmc", "kmc_tools")):
        # KMC has no version flag: it prints its banner when run with no argument (see tools.VERSION_FLAGS)
        print(f"{name} stub 1.0")
        return 0
    if value("fail", "") == name:
        print(f"{name}: deliberate failure", file=sys.stderr)
        return 1
    return PROGRAMS[name](argv)
