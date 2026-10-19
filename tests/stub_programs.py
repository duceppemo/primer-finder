"""Stand-ins for the external programs, so that the pipeline can be tested without installing them.

Each stub writes the files the real program would write. The `stubs` fixture puts a two-line launcher for
every program on PATH; what they produce is read from the PRIMER_FINDER_STUBS environment variable (JSON):

    kmers            the kmers "kmc_tools transform ... dump" writes (an empty list writes an empty file)
    dump_short       how many of them to leave out of the dump, as a full disk does, reporting success
    shared_kmers     the kmers the inclusion database holds, if not the same as `kmers`
    contigs          {name: sequence} for the assembler to write
    sam              [[query name, flag, cigar], ...] for minimap2 to write
    genome_kmers     {genome stem: [kmer, ...]} what each genome holds, so that KMC's operands are read
    presence         {genome: [contig, ...]} the contigs blastn finds in each inclusion genome
    presence_spans   {genome: {contig: [[qstart, qend], ...]}} which part of a contig it finds there
    exclusion_hits   {genome: [[query, qstart, qend, evalue, qseq, sseq], ...]} for the exclusion genomes
    fail             the name of a program that should exit 1
    kill             the name of a program that should die by a signal, as an OOM kill looks
    exit_code        {program: code} for a program that should exit with that code

For the amplicon searches of the design step, where blastn is asked for one query per assay:

    amplicon_hits    {genome: how many copies of each amplicon it holds}
    amplicon_mismatch   {position: base} the genome has instead, within the amplicon
    amplicon_qstart  the base of the amplicon the alignment starts at (1 by default)
    amplicon_covers  the fraction of the amplicon the alignment reaches (all of it by default)
    amplicon_insertion  {position: size} bases of the amplicon the genome has nothing for
    amplicon_split   [where to break, how much genome between] report each copy as two pieces
    amplicon_evalue  the evalue every amplicon row carries (default 1e-30)

And for Primer3:

    primer3_pairs    how many pairs to return (0 for "nothing could be designed")
    primer3_probe    False to return no internal oligo
    primer3_structures  False to report no hairpin or dimer temperatures
    primer3_error    the PRIMER_ERROR to report instead of an assay
    primer3_sleep    seconds to sleep, for the time-limit tests
    presence_evalue  the evalue the inclusion search reports

Every call is appended to the file named by PRIMER_FINDER_CALLS, one command line per line.
"""

from __future__ import annotations

import json
import os
import shlex
import signal
import sys
from pathlib import Path


def reverse_complement(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))


DEFAULT_KMERS = ["ACGT" * 5, "TTGA" * 5, "GGCA" * 5]
DEFAULT_CONTIGS = {"ctg1": "ACGT" * 25}
# A contig that differs from the exclusion genome by two mismatches 7 bases apart: a candidate.
DEFAULT_SAM = [["ctg1", 0, "40=1X7=1X50="]]
# The default exclusion alignment covers the whole of DEFAULT_CONTIGS["ctg1"], as a real blast would: a
# contig differing at two positions aligns end to end. An alignment covering only part of it would leave
# the rest *unaligned*, which the pipeline now reads as the exclusion group not holding that stretch --
# true of a real partial alignment, and not what a test about two close mismatches means to say.
DEFAULT_QSEQ = "ACGT" * 25
DEFAULT_SSEQ = ("ACGTA" + "A" + "GTACGT" + "T" + "CGTACGT" + "ACGT" * 20)  # Mismatches at 5 and 12


# Every key the stubs understand. A scenario asking for anything else is a mistyped fixture, and silently
# ignoring it has produced tests that passed while exercising nothing -- so it is an error instead.
KNOWN_KEYS = frozenset({
    "kmers", "shared_kmers", "dump_short", "genome_kmers", "contigs", "sam", "presence",
    "presence_evalue",
    "presence_spans", "exclusion_hits",
    "fail", "kill", "exit_code",
    "amplicon_hits", "amplicon_mismatch", "amplicon_qstart", "amplicon_covers", "amplicon_insertion",
    "amplicon_split", "amplicon_evalue",
    "primer3_pairs", "primer3_probe", "primer3_structures", "primer3_error", "primer3_sleep",
    # `stub_insilico_pcr.py` reads the same variable, so its keys are known here too
    "insilico_fail", "insilico_misses", "insilico_extra", "insilico_terminal", "insilico_counted",
    "insilico_five_prime", "insilico_no_report",
})


def scenario() -> dict:
    asked = json.loads(os.environ.get("PRIMER_FINDER_STUBS", "{}"))
    unknown = sorted(set(asked) - KNOWN_KEYS)
    if unknown:
        raise SystemExit(f"stub_programs: no such scenario key: {', '.join(unknown)}. "
                         f"Known keys: {', '.join(sorted(KNOWN_KEYS))}")
    return asked


def value(key: str, default):
    return scenario().get(key, default)


def argument(argv: list[str], flag: str, default: str = "") -> str:
    """The value that follows `flag`."""
    return argv[argv.index(flag) + 1] if flag in argv else default


MARKER = "stub kmers:"


def _joined(argv: list[str], prefix: str, default: int) -> int:
    """KMC writes its numbers onto the flag: `-ci2`, not `-ci 2`."""
    part = next((part for part in argv if part.startswith(prefix)), None)
    return int(part[len(prefix):]) if part else default


def _write_db(database: Path, kmers: list[str] | None) -> None:
    """The two files KMC writes. With `genome_kmers` in the scenario they carry what the database holds,
    so that the operands of the subtraction are read rather than guessed from the file name."""
    body = "" if kmers is None else MARKER + " " + " ".join(kmers) + "\n"
    for suffix in (".kmc_pre", ".kmc_suf"):
        database.with_name(database.name + suffix).write_text(body)


def _read_db(database: Path) -> list[str] | None:
    """What a stub database holds, or None for one written without `genome_kmers`."""
    path = database.with_name(database.name + ".kmc_pre")
    text = path.read_text() if path.is_file() else ""
    return text[len(MARKER):].split() if text.startswith(MARKER) else None


def kmc(argv: list[str]) -> int:
    """kmc -k.. -t.. -m.. -fm -ci.. -cx.. @list <db> <work_dir>

    With `genome_kmers` {genome stem: [kmer, ...]} in the scenario, the @list is read and `-ci`/`-cx` are
    applied to the pooled count across the listed files, as KMC applies them. Without it the database is
    empty and the dump comes from the `kmers`/`shared_kmers` keys, as before.
    """
    positional = [part for part in argv if not part.startswith(("-", "@"))]
    database = Path(positional[0])
    per_genome = value("genome_kmers", None)
    if per_genome is None:
        _write_db(database, None)
        return 0
    listed = next((part[1:] for part in argv if part.startswith("@")), "")
    paths = [line for line in Path(listed).read_text().split("\n") if line.strip()] if listed else []
    counts: dict[str, int] = {}
    for path in paths:
        for kmer in per_genome.get(Path(path).name.split(".")[0], []):
            counts[kmer] = counts.get(kmer, 0) + 1
    least, most = _joined(argv, "-ci", 1), _joined(argv, "-cx", 1_000_000_000)
    _write_db(database, sorted(kmer for kmer, count in counts.items() if least <= count <= most))
    return 0


def kmc_tools(argv: list[str]) -> int:
    """kmc_tools -t.. simple <a> <b> kmers_subtract <db>
    |   kmc_tools -t.. transform <db> [histogram <file>] dump <file>"""
    if "kmers_subtract" in argv:
        database = Path(argv[argv.index("kmers_subtract") + 1])
        minuend, subtrahend = Path(argv[argv.index("simple") + 1]), Path(argv[argv.index("simple") + 2])
        # The operands, in the order given: kmers_subtract keeps what the first holds and the second does
        # not. Reading them is what makes a swap visible -- it used to leave every test passing.
        held, removed = _read_db(minuend), _read_db(subtrahend)
        _write_db(database, None if held is None else
                  sorted(set(held) - set(removed or [])))
        return 0
    if "dump" in argv or "histogram" in argv:
        # The inclusion database is dumped too, to tell the two causes of an empty result apart. Unless a
        # scenario says otherwise, both dumps hold the same kmers.
        database = Path(argv[argv.index("transform") + 1])
        recorded = _read_db(database)
        if recorded is None:
            key = "kmers" if database.name.endswith("_specific") else "shared_kmers"
            recorded = value(key, value("kmers", DEFAULT_KMERS))
        if "histogram" in argv:
            # One line per count, as the real one writes it: every stubbed kmer is seen twice
            histogram = Path(argv[argv.index("histogram") + 1])
            histogram.write_text(f"2\t{len(recorded)}\n" if recorded else "")
        if "dump" in argv:
            # `dump_short` is what a full disk does: the real kmc_tools writes what it can, says nothing
            # and exits 0
            keep = len(recorded) - value("dump_short", 0)
            out_file = Path(argv[argv.index("dump") + 1])
            out_file.write_text("".join(f"{kmer}\t2\n" for kmer in recorded[:max(0, keep)]))
        return 0
    return 1  # pragma: no cover


def _write_contigs(path: Path, min_contig: int = 0) -> None:
    contigs: dict[str, str] = value("contigs", DEFAULT_CONTIGS)
    kept = {name: seq for name, seq in contigs.items() if len(seq) >= min_contig}
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in kept.items()))


def skesa(argv: list[str]) -> int:
    # The real one drops anything shorter than --min_contig, 200 by default
    asked = argument(argv, "--min_contig")
    _write_contigs(Path(argument(argv, "--contigs_out")), int(asked) if asked else 200)
    return 0


def spades(argv: list[str]) -> int:
    # SPAdes has no such option, so it writes everything and the caller does the filtering
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


def query_sequences(argv: list[str]) -> dict[str, str]:
    """The query fasta of a blastn call, as {name: sequence}."""
    sequences: dict[str, str] = {}
    name = None
    for line in Path(argument(argv, "-query")).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        elif name:
            sequences[name] += line.strip()
    return sequences


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
        # Counting the copies of an amplicon: one row per copy, each covering the whole query.
        # `amplicon_mismatch` puts a base of the genome's own at a position of the amplicon, which is what
        # the design step weighs a difference by.
        sequences, name = {}, None
        for line in Path(argument(argv, "-query")).read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = ""
            elif name:
                sequences[name] += line.strip()
        copies = value("amplicon_hits", {}).get(genome, 0)
        mismatches = {int(position): base
                      for position, base in value("amplicon_mismatch", {}).items()}
        # `amplicon_qstart` makes the alignment start partway into the amplicon, as a real hit may: the
        # design step has to shift every reported position by it to land on the region's coordinates.
        qstart = int(value("amplicon_qstart", 1))
        # `amplicon_covers` ends the alignment partway through the amplicon, as a real hit does where a
        # primer site is missing, and `amplicon_insertion` gives the genome a gap of that many bases at
        # that position of the amplicon, as a real hit does where the amplicon has bases the genome lacks.
        # Both are needed to reach either side of `design.has_the_amplicon`. The counts below are worked
        # out from the two aligned strings the way blast works them out, so they stay consistent.
        covers = float(value("amplicon_covers", 1.0))
        insertion = {int(position): int(size)
                     for position, size in value("amplicon_insertion", {}).items()}
        for query in queries:
            whole = sequences.get(query, "A" * 100)
            sequence = whole[qstart - 1:qstart - 1 + max(1, round(covers * len(whole)))]
            subject = list(sequence)
            for position, base in mismatches.items():
                if position < len(subject):
                    subject[position] = base
            for position, size in insertion.items():
                if position < len(subject):  # The genome has nothing where the amplicon has these bases
                    subject[position:position + size] = ["-"] * min(size, len(subject) - position)
            aligned_subject = "".join(subject)
            # Each copy sits somewhere else in the genome, far enough apart that nothing reads two of
            # them as one interrupted copy. `amplicon_split` reports a copy as two pieces instead of one,
            # the way blast does once gapping through an indel costs more than the rest of the alignment
            # is worth: [where in the amplicon to break, how much genome lies between the pieces].
            split = value("amplicon_split", None)
            for copy in range(copies):
                start = 1 + copy * 100_000
                pieces = [(0, len(sequence), 0)]
                if split:
                    pieces = [(0, int(split[0]), 0), (int(split[0]), len(sequence), int(split[1]))]
                between = 0
                for begin, end, gap in pieces:
                    piece, piece_subject = sequence[begin:end], aligned_subject[begin:end]
                    if not piece:
                        continue
                    between += gap
                    at = start + between + begin
                    identical_here = sum(1 for one, other in zip(piece, piece_subject, strict=True)
                                         if one != "-" and one.upper() == other.upper())
                    row = {
                        "qseqid": query, "qstart": str(qstart + begin),
                        "qend": str(qstart + end - 1),
                        "evalue": str(value("amplicon_evalue", "1e-30")),
                        "length": str(len(piece)),
                        "pident": f"{100.0 * identical_here / len(piece):.3f}",
                        "nident": str(identical_here),
                        "gaps": str(piece.count("-") + piece_subject.count("-")),
                        "sseqid": "chr", "sstart": str(at),
                        "send": str(at + len(piece_subject.replace("-", "")) - 1),
                        "qseq": piece, "sseq": piece_subject,
                    }
                    rows.append([row[name] for name in fields])
    elif "qseq" in fields:  # The exclusion step: aligned sequences
        hits = value("exclusion_hits", {})
        default = [[query, 1, len(DEFAULT_QSEQ), 1e-30, DEFAULT_QSEQ, DEFAULT_SSEQ] for query in queries]
        for query, qstart, qend, evalue, qseq, sseq in hits.get(genome, default):
            rows.append([str(query), str(qstart), str(qend), str(evalue), qseq, sseq])
    else:  # The inclusion step: is the contig there, and how much of it?
        present = value("presence", {}).get(genome, queries)
        evalue = value("presence_evalue", "1e-30")
        # `presence_spans` {genome: {contig: [[qstart, qend], ...]}} says which part of a contig the
        # genome aligns, as blast reports one row per high-scoring pair. The default is all of it, which
        # is what a genome holding the contig looks like.
        spans = value("presence_spans", {}).get(genome, {})
        sizes = {name: len(sequence) for name, sequence in query_sequences(argv).items()}
        for query in queries:
            if query not in present:
                continue
            for start, end in spans.get(query, [[1, sizes.get(query, 100)]]):
                row = {"qseqid": query, "qstart": str(start), "qend": str(end), "evalue": evalue}
                rows.append([row[name] for name in fields])
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

    sleep = value_of("primer3_sleep", 0)
    if sleep:  # To test what happens when Primer3 does not answer in time
        import time

        time.sleep(sleep)
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
                # Primer3 reports a right primer as the oligo to order, which is the reverse complement
                # of the template it sits on -- checked against primer3_core, whose PRIMER_RIGHT_0_SEQUENCE
                # is the reverse complement of template[start:end + 1] and not that slice itself.
                f"PRIMER_RIGHT_{number}_SEQUENCE="
                f"{reverse_complement(template[right_start:right_end + 1]) or 'T' * size}",
                f"PRIMER_RIGHT_{number}_TM=60.0",
                f"PRIMER_RIGHT_{number}_GC_PERCENT=50.0",
                f"PRIMER_PAIR_{number}_PRODUCT_SIZE={right_end - left + 1}",
                f"PRIMER_PAIR_{number}_PENALTY={0.5 + number}",
            ]
            if value_of("primer3_structures", True):
                # Real Primer3 reports these whenever it does the thermodynamic alignment, which is its
                # default. Distinct values, so that a column fed from the wrong one is visible.
                fields += [
                    f"PRIMER_LEFT_{number}_HAIRPIN_TH=35.70",
                    f"PRIMER_LEFT_{number}_SELF_ANY_TH=13.60",
                    f"PRIMER_LEFT_{number}_SELF_END_TH=2.10",
                    f"PRIMER_RIGHT_{number}_HAIRPIN_TH=24.30",
                    f"PRIMER_RIGHT_{number}_SELF_ANY_TH=11.20",
                    f"PRIMER_RIGHT_{number}_SELF_END_TH=1.40",
                    f"PRIMER_PAIR_{number}_COMPL_ANY_TH=8.70",
                    f"PRIMER_PAIR_{number}_COMPL_END_TH=4.10",
                ]
            if value_of("primer3_probe", True):
                fields += [
                    f"PRIMER_INTERNAL_{number}={probe_start},{size}",
                    f"PRIMER_INTERNAL_{number}_SEQUENCE={template[probe_start:probe_start + size] or 'C' * size}",
                    f"PRIMER_INTERNAL_{number}_TM=68.0",
                    f"PRIMER_INTERNAL_{number}_GC_PERCENT=55.0",
                ]
                if value_of("primer3_structures", True):
                    fields += [
                        f"PRIMER_INTERNAL_{number}_HAIRPIN_TH=33.30",
                        f"PRIMER_INTERNAL_{number}_SELF_ANY_TH=19.60",
                        f"PRIMER_INTERNAL_{number}_SELF_END_TH=3.20",
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
    if value("kill", None) == name:
        os.kill(os.getpid(), signal.SIGKILL)
    chosen = value("exit_code", {}).get(name)
    if chosen is not None:
        print(f"{name}: exiting {chosen} as the scenario asked", file=sys.stderr)
        return int(chosen)
    if argv in (["--version"], ["-version"]) or (not argv and name in ("kmc", "kmc_tools")):
        # KMC has no version flag: it prints its banner when run with no argument (see tools.VERSION_FLAGS)
        print(f"{name} stub 1.0")
        return 0
    if value("fail", "") == name:
        print(f"{name}: deliberate failure", file=sys.stderr)
        return 1
    return PROGRAMS[name](argv)
