# Changelog

## 1.3.0

### Added

- The `design` command rejects the oligos that cannot work, rather than proposing them: a hairpin or a dimer
  — with itself or with the other primer, over the whole oligo and at the 3' end where a polymerase could
  extend it — that melts at or above 47 °C now disqualifies a candidate, for the probe as well as the
  primers. The limits are Primer3's own thermodynamic model, so nothing extra has to be installed, and they
  are adjustable with `--max-hairpin-tm` and `--max-dimer-tm`. What was predicted is reported in
  `assays.tsv` as `pair_dimer_tm`, `pair_dimer_end_tm` and the per-oligo `*_hairpin_tm` and
  `*_self_dimer_tm`.
- The reaction those temperatures are predicted for is now stated instead of left at Primer3's bare
  defaults: 50 mM monovalent cations, 3 mM Mg²⁺, 0.8 mM dNTPs, 250 nM primers and 200 nM probe, which is an
  ordinary TaqMan qPCR. Each is an option (`--monovalent`, `--divalent`, `--dntp`, `--primer-conc`,
  `--probe-conc`) and all of them are recorded in `design_info.json`. They change both the temperatures
  reported and which oligos come back, so a run can be matched to the master mix it is for.
- `--gc-clamp`: how many G or C bases a primer must carry at its 3' end, where extension starts. One by
  default. Primer3 takes one clamp per request, so it is dropped from the requests that pin a 3' end on a
  differing base — that base is whatever the genomes made it, and with the clamp left on there, seven of
  eight forced-end requests on the *Xylella* set returned nothing. Those requests therefore ignore this
  option, and no clamp applies to the probe, which Primer3 has no tag for.
- A verdict of its own, `undecided (3' end)`, for an assay that covers the inclusion group but amplifies
  exclusion genomes only where a difference sits in the last two bases of a primer. `no` would claim the
  check had refused it on evidence worth believing, and it has not: those bases are not counted at any
  tolerance, so there was nothing it could refuse. Those assays rank above the ones refused on evidence the
  check can defend and below the ones it cleared, and `design_info.json` counts them.
- `qpcr_exclusion_terminal_only` and `pcr_exclusion_terminal_only` in `assays.tsv`: of the exclusion genomes
  an assay amplified, how many did so only where a difference sits in the last two bases of a primer — the
  bases insilicoPCR does not count. When it equals `*_exclusion_amplified`, every amplification behind the
  verdict is one the tool could not have refused, so it is the model's blind spot rather than a
  cross-reaction. Only the amplicons that needed no counted mismatch are weighed, which makes the count the
  same at any `--mismatches`, and a genome counts only when none of those is an exact match. On the
  *Xylella* set this separates 96 unjudgeable assays from the 149 the check refuses on its own evidence, and
  the run says so as it goes.

### Changed

- In silico PCR now lets each primer bind through **one** mismatch by default, instead of none
  (`--mismatches`). What that tolerance governs was measured rather than assumed: insilicoPCR does not count
  the last two bases of a primer at all — blast trims an unmatched base off its alignment — never amplifies
  through a mismatch 3 or 4 bases from the 3' end whatever the tolerance, and only counts the mismatches 5
  or more bases in, per primer. A single internal mismatch frequently does not stop a real reaction, so
  assuming it does flattered every assay resting on one difference, and it hid 12 of the 30 absence-based
  assays in the *Xylella* set that amplify something off-target in an exclusion genome. `--mismatches 2` is
  a stricter second pass; 3 measured nothing 2 did not.
- The claim that assays resting on a run of differences at a 3' end were never selective in silico is
  **withdrawn as evidence**. It is still what happens, at every tolerance, but the cause is that the
  alignment cannot see those bases — not anything about PCR. The sweep, the controlled test behind it and
  the numbers that replace it are in `validation/results/2026-10-08_v1.3.0/`.

- A difference is weighed by the base it replaces as well as by where it is: a position where the exclusion
  genomes have a G or a C counts for more than one where they have an A or a T, since G:C holds with three
  hydrogen bonds and A:T with two. The bases come from the exclusion genomes' own alignments of each
  amplicon, and the weight separates only assays covering the same number of differences. `assays.tsv`
  reports it (`primer_variant_weight`, `forward_strong_variants`, `reverse_strong_variants`).

## 1.2.0

### Added

- `primer-finder design`: Primer3 on the regions of a finished run, and the assays it proposes ordered by a
  heuristic for which to try first — not a prediction of bench behaviour, and nothing here is wet-lab
  tested. An assay works either because no exclusion genome holds its amplicon or
  because an oligo sits on bases that differ; one that is neither would amplify both groups and is set
  aside. Primer3 is also asked for primers whose 3' end is forced onto a run of differing bases, and for
  probes straddling one, so that allele-specific assays are proposed as well as convenient ones.
- The assays are written as the primer files
  [insilicoPCR](https://github.com/duceppemo/insilicoPCR) reads, twice over — with the probe and without,
  since it cannot report both kinds at once — together with a script that runs it over both groups.
  `--insilico-pcr` runs it directly and records, per assay and per mode, whether it amplifies every
  inclusion genome and no exclusion genome. When the regions were found with `-p` below 1, an assay that
  reaches that share counts too and is reported as `partial (88%)`: the regions themselves were allowed to
  be missing from some inclusion genomes, so the assays on them are judged on the same terms.
- The copies of each amplicon in the inclusion genomes are counted and preferred, since a repeated target
  usually improves the limit of detection. Note that `find -d 1`, the default, discards repeated regions
  before the design step can see them.
- A logo, in `docs/images/`.

### Fixed

- A Primer3 request whose answer is hard to find no longer holds up the command. Each request is run on its
  own with a time limit, and pinning the probe over a position — the expensive constraint, and the one that
  ran for 45 minutes on a real region — gets a shorter one. A request that runs out of time is dropped and
  counted; its region is still designed on from the others.
- Nothing is forced where the rest of the assay would not fit, and a forced request only searches the
  neighbourhood of what it aims at rather than the whole region.
- `primer3_core --version` prints nothing, so `run_info.json` recorded a copyright line as its version;
  `--about` is used instead.

## 1.1.0

### Added

- `-p/--min-inclusion`: the fraction of the inclusion genomes a region has to be in, instead of all of them.
  The fraction is rounded up to a whole number of genomes, and both the kmer counting and the blast check
  use it. Below 1, each record in `all_inclusion_contigs.fasta` and `final_kmers.fasta` carries
  `inclusion=<holding>/<total>`, and the regions that cover the most inclusion genomes come first.
  The exclusion side is unchanged: a region must still be absent from every exclusion genome. This was the
  first of the two things the 2022 README listed as wanted.
- The package is published to PyPI as well as bioconda (`pip install primer-finder` gives the command, not
  the programs it runs).
- A validation suite (`validation/xylella/`) that recovers four published *Xylella fastidiosa* subspecies
  qPCR assays from public genomes, with a dated record of the run.

### Changed

- A file that is not a fasta is reported with the count and the first few names, instead of every path: a
  `-i` pointing at the wrong folder could print thousands of them.

## 1.0.0

First packaged release. The two scripts (`primer_finder.py` and `IDT_results_converter.py`) became the
`primer_finder` package with a `primer-finder` command, a test suite and documentation. The pipeline itself
is the same one: shared inclusion kmers (KMC), minus the exclusion kmers, assembled (SKESA or SPAdes),
mapped to one exclusion genome (minimap2) and checked against every genome (blast).

### Added

- A `primer-finder` command (`primer-finder -i inclusion/ -e exclusion/ -o results/`), installable from
  bioconda, plus `primer-finder idt` for IDT order sheets. `primer_finder.py` and
  `IDT_results_converter.py` still work.
- `run_info.json` with the parameters, the genomes, the version of every program used and the counts of
  each step, and `primer_finder.log` with the whole log, in the output folder.
- `--keep-intermediate` keeps the KMC databases, the SAM file and the blast databases; the run removes them
  by default.
- Gzipped genomes (`.fa.gz`, `.fasta.gz`, `.fna.gz`) work throughout: they are decompressed where blast
  needs plain fasta.
- A bundled example with simulated genomes and a known answer (`example/run_example.sh`), run by the CI.

### Changed

- No Python dependency at all: pysam, Biopython, pandas, psutil, pyahocorasick and xlrd are gone.
  Biopython's `NcbiblastnCommandline` wrapper, which the previous version used, no longer exists.
- bowtie2 and samtools are no longer needed. The contigs are mapped once with minimap2, whose SAM output is
  read directly, instead of being split into one file per contig and mapped one by one.
- The exclusion genome the contigs are mapped to is the first one in alphabetical order, not a random one,
  so that a run can be repeated. `-r` still forces a particular genome.
- The candidates are checked against each exclusion genome separately. A variant position has to differ in
  at least 90% of the exclusion genomes that the contig hits, as before, but the genomes are now counted
  instead of the blast alignments, and the inclusion and exclusion genomes are searched in parallel.
- The results are in subfolders (`1_kmers/`, `2_assembly/`, `3_candidates/`, `4_blast/`), with
  `final_kmers.fasta` at the top. `inclusion_blast_hits.tsv` holds 1 and 0 instead of `True` and `False`.
- Every hit of a contig in an exclusion genome is examined, not only the best one, and a position counts as
  different in that genome only if every copy of the region differs there: one matching copy is enough for
  an assay to amplify it. A genome whose hits do not reach a position no longer counts towards the 90%
  either, since it does not hold that region at all.
- SPAdes is called with `--only-assembler` (kmers carry no quality values to correct) and `--memory`, which
  the previous version did not pass.
- The blast databases are built in a folder of their own per genome, named after the genome's position in
  the list, and blast runs in that folder with relative paths. Two genomes with the same file name in
  different subfolders (`sampleA/contigs.fasta`, `sampleB/contigs.fasta`) therefore no longer collide, and
  the columns of `inclusion_blast_hits.tsv` are named after the path inside the input folder.
- An output folder that overlaps an input folder, or whose path contains a space (which BLAST cannot use in
  a database path), is refused before any work is done. The input folders may contain spaces.
- Python 3.10 or later (was 3.6).

### Fixed

- The variant positions in the headers of `final_kmers.fasta` are positions in the contig. They used to be
  positions in the blast alignment, which differ as soon as the hit does not start at the first base of the
  contig or holds a gap, so the lower-case bases could be the wrong ones.
- A contig that maps to the minus strand of the exclusion genome has its differences marked on the right
  bases: the cigar string describes the reverse complement, which the previous version applied to the
  forward sequence.
- Deletions no longer shift the lower-case marks: a deletion consumes no base of the contig.
- The rule "two differences within 21 bases" is applied to insertions and deletions as well as mismatches,
  and no longer misses the differences in the last three operations of the cigar string.
- `-k`/`--kmer_size` is checked (1-256): the previous check could not fail.
- `makeblastdb` no longer writes its index files into the inclusion and exclusion folders, which fails when
  they are read only.
- A genome listed in both groups, an output folder holding the input genomes, a file that is not a fasta and
  an unreadable kmer size are reported as errors instead of running anyway.
- A missing program is reported once, with the conda package that provides it, instead of failing midway.
- All the exclusion kmers are now subtracted. The old `-cx1e9` was read by KMC as `1`, so only kmers seen
  exactly once in the whole exclusion group were removed, and candidates in repeated exclusion regions
  survived. This changes the result of every run.
- With more than 255 inclusion genomes, `-ci <N>` could not match anything, because KMC counts to 255
  unless `-cs` says otherwise: such a run reported that no inclusion-specific kmer existed. The ceiling is
  now raised when the group needs it.
- Two differences that touch each other (a mismatch next to an inserted base, or two mismatches on either
  side of a deletion) count as two differences; the rule used to need a run of matching bases between them.
  A clipped end of at least 21 bases — a whole primer's worth of sequence the exclusion genome does not hold
  — now keeps the contig too.
- Two shared variant positions are enough to keep a contig, as documented; an off-by-three in the previous
  loop meant four were needed.
- `best_kmers.fasta` lists the contigs with the most differing bases first, so a contig the exclusion genome
  does not hold at all comes first instead of last.
- Fasta files whose headers are not UTF-8 (an accented strain name) are read instead of ending the run with
  a `UnicodeDecodeError` traceback; duplicate record names in an assembly are reported.
- An unreadable input, a read-only output folder or an `-o` that names an existing file is reported as an
  error instead of a traceback.
- The memory limit of a control group (Slurm, systemd, a container) is found wherever the job's cgroup sits,
  not only at the root of the cgroup filesystem, so the default `-m` is 85% of what the job may use rather
  than of the whole machine.
- `primer-finder idt` reads every format Excel saves: "CSV UTF-8" (with its byte-order mark), the legacy
  Windows CSV, "Unicode Text" (UTF-16) and the semicolon-separated CSV of comma-decimal locales. It no
  longer takes a styled but empty first row as the header, nor glues a phonetic run onto a cell value, and
  it says which sheet it read when a workbook holds several.
