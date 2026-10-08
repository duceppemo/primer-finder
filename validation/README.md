# Validation

The unit tests and the bundled `example/` use simulated data: they check that the code does what it says.
This folder checks something else — that primer-finder, run on public genomes, finds the regions that
published diagnostic assays were designed on.

Each record of a validation run is kept in `validation/results/<date>_v<version>/`, with the numbers it
produced. Records are never rewritten: a later run that contradicts one gets its own record, and the old one
keeps an erratum pointing at it.

## xylella: three published subspecies assays

`validation/xylella/` reproduces the target selection of

> Dupas E., Briand M., Jacques M.-A., Cesbron S. (2019) *Novel tetraplex quantitative PCR assays for
> simultaneous detection and identification of Xylella fastidiosa subspecies in plant tissues.*
> Frontiers in Plant Science 10:1732. [doi:10.3389/fpls.2019.01732](https://doi.org/10.3389/fpls.2019.01732)

Their primers and probes were designed with SkIf, a kmer-based signature tool, on 58 *Xylella fastidiosa*
genomes split into an ingroup and an outgroup — the same question primer-finder answers. The paper publishes
the oligo sequences and the specific regions they sit in, which makes it a ready-made answer key.

The test: for each of three subspecies, take that subspecies as the inclusion group and all the other
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
| `accessions.tsv` | the 22 complete RefSeq genomes used, with subspecies and strain, pinned by accession |
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
