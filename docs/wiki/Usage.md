# Usage

```bash
primer-finder -i inclusion/ -e exclusion/ -o results/
```

`primer-finder find ...` is the same command written in full; `primer-finder idt ...` is the
[order-sheet converter](#converting-an-idt-order-sheet).

## Input

Two folders of **assembled** genomes, one file per genome:

- accepted extensions: `.fasta`, `.fa`, `.fna`, and the gzipped `.fasta.gz`, `.fa.gz`, `.fna.gz`;
- subfolders are searched too, and symbolic links are followed, so the folders can hold links to a central
  collection of genomes;
- anything else in the folders is ignored, but a file with an accepted extension that is not a fasta stops
  the run;
- the same genome cannot be in both groups (primer-finder compares the resolved paths, so a link to it
  counts).

The two groups must be curated. primer-finder keeps only perfect matches, so a contaminated assembly, a
misassembly or a genome put in the wrong group is enough to lose every candidate region. A tool such as
[genome_comparator](https://github.com/duceppemo/genome_comparator) helps to check that the inclusion
genomes really belong together.

### Splitting a collection into the two groups

With a folder of genomes and a list of the ones that belong to the inclusion group:

```bash
# Every genome, and the inclusion ones, both sorted
find -L /all_genomes -type f -name "*.fasta" | sort > all.list   # -L also follows symbolic links
sort inclusion.list -o inclusion.list

# The rest is the exclusion group
comm -3 all.list inclusion.list | sort > exclusion.list

# Folders of symbolic links, which is what primer-finder reads
mkdir -p inclusion exclusion
while read -r genome; do ln -sf "$genome" inclusion/; done < inclusion.list
while read -r genome; do ln -sf "$genome" exclusion/; done < exclusion.list
```

## Options

| Option | Default | What it does |
|---|---|---|
| `-i`, `--inclusion` | required | Folder of genomes the assay should amplify. |
| `-e`, `--exclusion` | required | Folder of genomes the assay must not amplify. |
| `-o`, `--output` | required | Folder for the results; it is created if needed. It must be outside both input folders, and its path must not contain a space (a BLAST limitation). |
| `-t`, `--threads` | every CPU | CPUs for KMC, the assembler, minimap2, and how many genomes are blasted at a time. |
| `-m`, `--memory` | 85% of the memory | Memory in GB, passed to KMC and the assembler. |
| `-k`, `--kmer_size`, `--kmer-size` | 99 | Kmer size for KMC, 1-256. Shorter kmers find shorter specific regions, and more of them; longer kmers are more specific. |
| `-d`, `--duplication` | 1 | How many times a kmer may occur in each inclusion genome. 1 discards repeated regions, which make poor assays. |
| `-p`, `--min-inclusion` | 1.0 | The fraction of the inclusion genomes a kmer has to be in, rounded up to a whole number of genomes. 1.0 means every one of them. |
| `-r`, `--reference` | the first exclusion genome, alphabetically | Which exclusion genome the contigs are mapped to, to find their differences. It must be one of the genomes in the exclusion folder. |
| `-a`, `--assembler` | `skesa` | `skesa` or `spades`, to assemble the inclusion-specific kmers. |
| `--keep-intermediate` | off | Keep the KMC databases, the kmer dump, the SAM file and the blast databases. |
| `--debug` | off | Log every command line and the output of every program. |
| `-v`, `--version` | | Print the version. |

`-t` and `-m` are capped at what the machine (or the control group a job scheduler put the run in) actually
offers, with a warning.

### Accepting a region that some inclusion genomes lack

By default a region has to be in **every** inclusion genome, which is what makes the assay amplify all of
them — and which one incomplete assembly is enough to prevent. `-p/--min-inclusion` lowers that bar:

```bash
primer-finder -i inclusion/ -e exclusion/ -o results/ -p 0.9
```

With 12 inclusion genomes, `-p 0.9` asks for the regions that at least 11 of them share (0.9 × 12 = 10.8,
rounded up); the run says which number it settled on. The exclusion side is untouched: a region must still
be absent from every exclusion genome.

What this buys and what it costs:

- it recovers regions lost to a gap, a misassembly or genuine diversity in one or two genomes;
- the assay designed on such a region **will not amplify the genomes that lack it**. Each record then says
  how many hold it (`inclusion=11/12`), the widest come first, and `4_blast/inclusion_blast_hits.tsv` names
  the genomes that are missing it.

Use it to find out *why* a run came up empty before using it to design on: if a region appears at `-p 0.9`
but not at 1.0, the table tells you which genome is responsible, and that genome may be the problem rather
than the threshold.

## What to expect

A run prints what it is doing and ends with the number of candidate regions:

```
10:21:33 [INFO] primer-finder 1.0.0
10:21:33 [INFO] 12 inclusion genome(s), 40 exclusion genome(s)
10:21:33 [INFO] Counting the 99-mers shared by the 12 inclusion genomes...
10:22:04 [INFO] Counting the 99-mers of the 40 exclusion genomes...
10:23:41 [INFO] Subtracting the exclusion kmers from the inclusion ones...
10:23:48 [INFO] Found 18422 inclusion-specific 99-mers
10:23:48 [INFO] Assembling the kmers with skesa...
10:23:55 [INFO] Assembled 34 contig(s)
10:23:55 [INFO] Mapping the contigs to the exclusion genome exclusion_01.fasta...
10:23:56 [INFO] Keeping the contigs with at least two differences within 21 bases...
10:23:56 [INFO] 21 contig(s) could carry a selective assay
10:23:56 [INFO] Checking the candidates against the 12 inclusion genomes...
10:24:09 [INFO] 17 contig(s) are present in all inclusion genomes
10:24:09 [INFO] Checking the differences against the 40 exclusion genomes...
10:24:38 [INFO] Final number of contigs: 9
10:24:38 [INFO] Results: results/final_kmers.fasta
```

The same log is written to `results/primer_finder.log`, and the parameters, the genomes, the version of every
program and the counts of each step to `results/run_info.json`. See [Outputs](Outputs).

A run that finds nothing is not a failure of the program: it means there is no region that every inclusion
genome shares and no exclusion genome has. [FAQ](FAQ) lists what to try.

## Performance

The slow steps are KMC on the exclusion group and the two blast rounds. As a rough guide, 50 bacterial
genomes with `-t 16` take a few minutes. Memory is KMC's main need: `-m` is what it is allowed to use before
it starts writing temporary files to the output folder.

Runs are not resumable: a repeated command starts again from the kmers. The intermediate files of a finished
run are deleted unless `--keep-intermediate` is given.

## Designing assays on the regions

`primer-finder design` runs Primer3 on the regions of a finished run, keeps the assays that could tell the
two groups apart, and writes the primer files [insilicoPCR](https://github.com/duceppemo/insilicoPCR) reads:

```bash
primer-finder design results/ -o assays/ -t 16
primer-finder design results/ -o assays/ -t 16 --insilico-pcr /path/to/insilicoPCR-linux-x64
```

| Option | Default | What it does |
|---|---|---|
| `results` | required | The output folder of a finished `find` run. |
| `-o`, `--output` | required | Folder for the designed assays. |
| `-i`, `--inclusion`, `-e`, `--exclusion` | what the run recorded | The genome folders, if they have moved. |
| `-t`, `--threads` | every CPU | For Primer3's blast checks and for insilicoPCR. |
| `--product-size` | `70-150` | Amplicon size range for Primer3. |
| `--assays-per-region` | 3 | How many assays Primer3 proposes per request. |
| `--gc-clamp` | 1 | G or C bases required at a primer's 3' end. 0 for none. Primer3 takes one clamp per request, so this does not reach the requests that pin a 3' end on a differing base. |
| `--max-hairpin-tm` | 47 °C | Reject an oligo whose hairpin melts at or above this. |
| `--max-dimer-tm` | 47 °C | Reject an oligo that pairs with itself, or with its partner, at or above this. |
| `--monovalent` | 50 mM | Monovalent cations in the reaction, for the melting temperatures. |
| `--divalent` | 3 mM | Magnesium in the reaction. |
| `--dntp` | 0.8 mM | Total dNTPs (0.2 mM of each). |
| `--primer-conc` | 250 nM | Primer concentration. |
| `--probe-conc` | 200 nM | Probe concentration. |
| `--max-regions` | 50 | Design on this many regions, most promising first. 0 for all of them. |
| `--insilico-pcr` | off | Run in silico PCR of the assays against both groups. Takes the folder of an extracted portable release, its jar, or a launcher script. |
| `--mismatches` | 1 | Mismatches insilicoPCR lets each primer bind through (`-m`). One is often not enough to stop a reaction; 0 assumes any mismatch is, 2 is a stricter check. It never applies to the last two bases of a primer. |

The defaults of `--max-hairpin-tm` through `--probe-conc` are an ordinary TaqMan qPCR; set them to your own
master mix if it differs, since they change both the melting temperatures Primer3 predicts and which oligos
it returns.

The answer is `assays/assays.tsv`, the most promising assays first. [Designing assays](Designing-assays)
explains what an oligo has to satisfy, what the order means and how an assay is judged;
[Outputs](Outputs) lists the columns.

## Converting an IDT order sheet

Once an assay has been designed on a candidate region and ordered, the IDT order sheet can be turned into a
fasta file of oligos:

```bash
primer-finder idt order.xlsx assays.fasta my_target
```

It reads `.xlsx` and `.xlsm` (no Excel or extra library needed) and `.csv`, `.tsv` and `.txt` files with a
`Type` column (`Forward Primer`, `Probe`, `Reverse Primer`), a `Sequence` column and, optionally, an
`Amplicon` column; the columns may be in any order and any case, and rows that are not oligos are ignored.
Every format Excel saves is accepted: "CSV UTF-8" (with its byte-order mark), the legacy Windows "CSV",
"Unicode Text" (UTF-16), and the semicolon-separated CSV of locales where the comma is the decimal
separator. Of a workbook, the **first sheet** is read; if there are several, the run says which one it took.
Each assay becomes up to three records:

```
>my_target_0_120bp-F
AAAACCCC...
>my_target_0_120bp-P
TTTTGGGG...
>my_target_0_120bp-R
CCCCAAAA...
```

The prefix is optional. `python IDT_results_converter.py sheet.xlsx out.fasta prefix` still works and does
the same thing.
