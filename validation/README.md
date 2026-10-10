# Validation

The unit tests and the bundled `example/` use simulated data: they check that the code does what it says.
This folder checks something else — that primer-finder, run on public genomes, finds the regions that
published diagnostic assays were designed on.

Each record of a validation run is kept in `validation/results/<date>_v<version>/`, with the numbers it
produced. Records are never rewritten: a later run that contradicts one gets its own record, and the old one
keeps an erratum pointing at it.

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

`check_region.py` is deliberately not a pass or fail: the region primer-finder reports is not the same object
as the paper's long-mer. The exclusion genomes are not the same ones, and a region here is a contig assembled
from inclusion-specific kmers rather than a sequence bounded by the outgroup. Whether the assay sits inside
the region is the question that matters, and `check_assay.py` answers it.
