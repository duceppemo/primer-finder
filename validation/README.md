# Validation

The unit tests and the bundled `example/` use simulated data: they check that the code does what it says.
This folder checks something else — that primer-finder, run on public genomes, finds the regions that
published diagnostic assays were designed on.

Each record of a validation run is kept in `validation/results/<date>_v<version>/`, with the numbers it
produced. Records are never rewritten: a later run that contradicts one gets its own record, and the old one
keeps an erratum pointing at it.

 re-designs the Brucella regions under 1.8.0 to measure what the
corrected amplicon rule changes, and reports a defect it found in that rule.

Not every record is about an assay. `2026-10-14_v1.7.0_unitig/` measures what a change to the code does to
a panel whose answer is already known — the question any change to the find or assembly rules has to
answer before it can become a default.

## xylella: four published subspecies assays

`validation/xylella/` reproduces the target selection of

> Dupas E., Briand M., Jacques M.-A., Cesbron S. (2019) *Novel tetraplex quantitative PCR assays for
> simultaneous detection and identification of Xylella fastidiosa subspecies in plant tissues.*
> Frontiers in Plant Science 10:1732. [doi:10.3389/fpls.2019.01732](https://doi.org/10.3389/fpls.2019.01732)

Their primers and probes were designed with SkIf, a kmer-based signature tool, on 58 *Xylella fastidiosa*
genomes split into an ingroup and an outgroup — the same question primer-finder answers. The paper publishes
the oligo sequences and the specific regions they sit in, which makes it a ready-made answer key.

The test: for each of four subspecies, take that subspecies as the inclusion group and all the other
subspecies as the exclusion group, run primer-finder, and look for the published primers and probe in the
regions it reports. An assay counts as recovered when all three oligos are found exactly, in one region, with
the primers facing each other, the probe between them, and the amplicon the published length.

```bash
bash validation/xylella/run.sh /path/to/work 16 32
```

It needs primer-finder and its programs, plus the NCBI datasets command line tool
(`conda install -c conda-forge ncbi-datasets-cli`) to fetch the genomes. It takes a few minutes, most of it
downloading. The pieces:

| File | What it holds |
|---|---|
| `accessions.tsv` | the 25 complete RefSeq genomes used, with subspecies and strain, pinned by accession |
| `assays.tsv` | the published primers and probes (their Table 3) |
| `regions.tsv` | the specific regions the paper reports, with coordinates in its reference genomes (their Table 2) |
| `get_genomes.py` | downloads the genomes into one folder per subspecies |
| `group_genomes.py` | builds the inclusion and exclusion folders of one target, as symbolic links |
| `check_assay.py` | looks for the oligos in `final_kmers.fasta` and decides whether the assay was recovered |
| `check_region.py` | says where a recovered region sits in the paper's reference genome, as context |

`check_region.py` is deliberately not a pass or fail: the region primer-finder reports is not the same object
as the paper's long-mer. The exclusion genomes are not the same ones, and a region here is a contig assembled
from inclusion-specific kmers rather than a sequence bounded by the outgroup. Whether the assay sits inside
the region is the question that matters, and `check_assay.py` answers it.

## listeria: 30 published clonal-complex assays, and much larger groups

`validation/listeria/` asks the same question of groups one to two orders of magnitude larger, against

> Félix B. et al. (2023) *Identification by High-Throughput Real-Time PCR of 30 Major Circulating Listeria
> monocytogenes Clonal Complexes in Europe.* Microbiology Spectrum 11(3):e03954-22.
> [doi:10.1128/spectrum.03954-22](https://doi.org/10.1128/spectrum.03954-22)

Their 34 TaqMan sets were designed on kmers found in 954 genomes, at least 15 per clonal complex, and
checked against a further 2,388 — again the same question primer-finder answers. Clonal complexes within one
species are a harder target than subspecies, and the panel is far bigger: 18 to 86 inclusion genomes against
677 to 745 exclusion ones.

```bash
conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
conda activate primer-finder_listeria
bash validation/listeria/run.sh /path/to/work 32 64
```

| File | What it holds |
|---|---|
| `assays.tsv` | the 34 published primer-and-probe sets (their Table 1), without the dyes and quencher |
| `panel.tsv` | the 784 complete RefSeq genomes the record used, with the ST, clonal complex and lineage of each |
| `census.tsv` | where each published assay actually occurs in that panel, produced by `census.py` |
| `type_genomes.py` | downloads the genomes and gives each one a clonal complex from its MLST profile |
| `census.py` | blasts every published oligo against every genome: the ceiling on what could be recovered |
| `group_genomes.py` | builds the inclusion and exclusion folders of one clonal complex |
| `check_assay.py` | looks for the oligos in `final_kmers.fasta` and decides whether the assay was recovered |

The clonal complex is not in the NCBI metadata: it comes from the 7-locus MLST profile, and `type_genomes.py`
takes the ST → complex assignment from the Institut Pasteur scheme that `mlst` installs, rather than
computing one.

**Read `census.tsv` beside the results.** primer-finder reports a region only when it is in every inclusion
genome and differs from every exclusion genome, so an assay the census shows is missing from one genome of
its own complex, or present in one genome of another, cannot be reported — and should not be. On this panel
that is true of four of the nine complex-wide assays, which is why they are not recovered. The record sets
out which, and why.

## brucella: a genus assay and a three-species multiplex, on the tool's own ancestor

`validation/brucella/` is the first record where the published assays were found with **this software**
(one of the paper's eleven authors maintains it):

> Hoover E., Chmara J., Duceppe M., Charron P., Lloyd D., Huang H., Surujballi O.,
> Thanthrige-Don N., Arnold K., Hazelwood J., Kang M. (2026) *Two novel real-time PCR assays for
> Brucella detection: a species-specific multiplex and a genus-level singleplex developed via
> large-scale whole genome sequencing data analysis.* Microbiology Spectrum 14:e02974-25.
> [doi:10.1128/spectrum.02974-25](https://doi.org/10.1128/spectrum.02974-25)

Its four assays — one genus-level singleplex and a *B. abortus* / *B. melitensis* / *B. suis* multiplex —
were designed on unique core sequences found with `selective_primer_finder 0.2.1` and RUCS 1.0.3, on 1,597
assemblies the authors built from SRA reads. primer-finder is the rewrite of that 0.2.1, so this record
asks whether the rewrite still finds what its ancestor found, on a panel rebuilt from RefSeq.

Two things make it a harder test than the first two records:

- **the groups are closely related and very unequal.** *B. canis* sits inside *B. suis* phylogenetically,
  and the exclusion group of a species target is every other classical species — 1,571 genomes against 120
  for *B. suis*;
- **the panel is draft-heavy.** 873 of the 1,376 classical genomes are contig-level assemblies, which is
  the hostile case for a rule that asks for an exact match in *every* inclusion genome.

```bash
conda install -c conda-forge ncbi-datasets-cli
bash validation/brucella/run.sh /path/to/work 48 128
```

| File | What it holds |
|---|---|
| `assays.tsv` | the four published primer-and-probe sets (their Table 1), without the dyes and quencher |
| `groups.py` | what each target's inclusion and exclusion group is made of, and why |
| `panel.py` | rebuilds the panel from RefSeq and labels each genome; the paper's accessions are not published |
| `verify_amplicon.py` | measures each published amplicon in a real genome before it is trusted as a key |
| `census.py` | blasts every published oligo against every genome: the ceiling on what could be recovered |
| `mislabelled.txt` | the assemblies whose RefSeq species label the census shows cannot be right |
| `group_genomes.py` | builds the inclusion and exclusion folders of one target, `genus_wide` included |
| `check_assay.py` | looks for the oligos in `final_kmers.fasta` and decides whether the assay was recovered |
| `check_design.py` | asks whether `design` proposes the published oligos on the region that holds them |
| `score_assays.py` | any set of assays against the whole panel on one yardstick: published and proposed alike |
| `thermo.py` | Tm and secondary structure of any set of oligos, under the design command's own conditions |

**The panel is rebuilt, not reproduced.** Supplementary Table S2 of the paper gives only the number of
assemblies per species, not their accessions, so `panel.py` applies the same rule to today's RefSeq: every
assembly of the eight classical species the paper's panel holds, and of the four non-classical groups it
used as negatives. RefSeq has grown since October 2020, so the groups are larger than the paper's — 1,376
classical against 315 non-classical, where the paper had 1,160 against 194.

There are six targets, defined in `groups.py`: `genus` (the eight classical species the paper's panel
holds), `genus_wide` (those plus *B. pinnipedialis*, *B. amazoniensis*, *B. vulpis* and *B. inopinata*,
which the paper's panel had none of), one for each of the three species, and `suis_lato` (*B. suis*
together with *B. canis*, which lies inside it phylogenetically). The eight-species target is kept as it
was so that the record of the published assays stays reproducible.

**Read `census.tsv` beside the results**, as in the *Listeria* record, and `mislabelled.txt` beside that.
Eight assemblies of this panel carry all three oligos of exactly one species assay and it is not the assay
of their own RefSeq label, while no genome of the panel carries two. Left in the exclusion group they
remove a region that is specific to the species they really belong to, so every target is run twice: once
on the panel as labelled, once with those eight held out.
