# Changelog

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
