<p align="center">
  <img src="docs/images/primer-finder-logo.svg" alt="primer-finder" width="680">
</p>

<p align="center">
  <a href="https://github.com/duceppemo/primer-finder/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/primer-finder/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://codecov.io/gh/duceppemo/primer-finder"><img src="https://codecov.io/gh/duceppemo/primer-finder/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/primer-finder/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/primer-finder?label=release&cacheSeconds=3600" alt="Latest release"></a>
  <a href="https://anaconda.org/bioconda/primer-finder"><img src="https://img.shields.io/conda/vn/bioconda/primer-finder?label=bioconda" alt="Bioconda"></a>
  <a href="https://pypi.org/project/primer-finder/"><img src="https://img.shields.io/pypi/v/primer-finder?label=pypi" alt="PyPI"></a>
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <a href="https://github.com/duceppemo/primer-finder/blob/main/LICENSE"><img src="https://img.shields.io/github/license/duceppemo/primer-finder" alt="License: MIT"></a>
  <a href="https://github.com/duceppemo/primer-finder/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23226411"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23226411.svg" alt="DOI"></a>
</p>

primer-finder finds the sequences that tell one group of genomes from another, so that a selective (q)PCR
assay can be designed on them. Give it two folders of assembled genomes — the ones the assay should amplify
(inclusion) and the ones it must not (exclusion) — and it reports the regions that every inclusion genome
carries, that no exclusion genome carries, and whose differences are close enough together to sit in one
primer or probe.

```
inclusion/ ──┐                                                        ┌─► final_kmers.fasta
             ├─► kmers ─► subtract ─► assemble ─► map ─────► blast ───┤   the candidate regions
exclusion/ ──┘   (KMC)     (KMC)     (SKESA or   (minimap2)  (every   │
                                      SPAdes)                genome)  └─► primer-finder design
                                                                          Primer3 + in silico PCR
                                                                          ─► assays.tsv
```

## Quick start

```bash
conda create -n primer-finder -c conda-forge -c bioconda primer-finder
conda activate primer-finder

primer-finder -i inclusion/ -e exclusion/ -o results/
```

> **The bioconda package is awaiting review** ([bioconda-recipes#70034](https://github.com/bioconda/bioconda-recipes/pull/70034)).
> Until it is merged, use `pip install primer-finder` in an environment that already holds KMC, SKESA or
> SPAdes, minimap2 and BLAST, or install from the source code with `environment.yml`, which brings them.

`inclusion/` and `exclusion/` hold one assembled genome per file (`.fasta`, `.fna`, `.fa`, gzipped or not;
subfolders and symbolic links are followed). The answer is `results/final_kmers.fasta`: one record per
candidate region, the specific bases in lower case and their positions in the header, the most promising
first. `results/run_info.json` records the parameters, the genomes and the version of every program used.

To turn those regions into assays, with Primer3 and a check of what would actually make each one selective:

```bash
primer-finder design results/ -o assays/ --insilico-pcr /path/to/insilicoPCR-linux-x64
```

`--insilico-pcr` is optional and points at [insilicoPCR](https://github.com/duceppemo/insilicoPCR), a
separate program: given it, the designed assays are amplified in silico against both groups and each one is
marked with whether it amplifies every inclusion genome and no exclusion genome.

Every oligo proposed satisfies Primer3's length, melting temperature and GC limits and folds or dimerises —
with itself and with its partner — below 47 °C; and unless the request pinned a 3' end on a differing base,
each primer carries a G or a C at that end. All of it is predicted for an ordinary TaqMan reaction whose salt, dNTP and oligo
concentrations are options.
So a doomed assay is never proposed. See
[Designing assays](https://github.com/duceppemo/primer-finder/wiki/Designing-assays).

To install from the source code instead, see
[Installation](https://github.com/duceppemo/primer-finder/wiki/Installation).

Given the genomes of five *Xylella fastidiosa* subspecies, it reports the regions that four independently
published subspecies-specific qPCR assays were designed on, three of them in the top five candidates — see
[Validation](https://github.com/duceppemo/primer-finder/wiki/Validation).

primer-finder only keeps perfect matches: a kmer must be in **all** the inclusion genomes with no mismatch,
and in **none** of the exclusion genomes (`-p/--min-inclusion` relaxes the first half when some of the
inclusion genomes are incomplete). It is therefore very sensitive to the quality of the assemblies and
to how the genomes were assigned to the two groups. Curate the input genomes;
[genome_comparator](https://github.com/duceppemo/genome_comparator) helps with that.

To check an installation, run the bundled example (simulated genomes with a known answer, a few seconds).
It is in the repository, not in the conda package:

```bash
curl -sL https://github.com/duceppemo/primer-finder/archive/refs/tags/v1.4.0.tar.gz | tar -xz --strip-components=1 primer-finder-1.4.0/example
bash example/run_example.sh
```

## Ordering the assays

Once an assay has been designed from a candidate region and ordered, `primer-finder idt` turns the IDT order
sheet into a fasta file of oligos, one record per primer and probe:

```bash
primer-finder idt order.xlsx assays.fasta my_target
```

## Documentation

Everything else is in the [wiki](https://github.com/duceppemo/primer-finder/wiki), whose sources are
maintained in [`docs/wiki`](https://github.com/duceppemo/primer-finder/blob/main/docs/wiki):

| Page | Contents |
|---|---|
| [Installation](https://github.com/duceppemo/primer-finder/wiki/Installation) | conda, bioconda, from source, checking the installation |
| [Usage](https://github.com/duceppemo/primer-finder/wiki/Usage) | inputs, every option, choosing the two groups, performance |
| [Designing assays](https://github.com/duceppemo/primer-finder/wiki/Designing-assays) | Primer3 on the regions, how assays are scored, in silico PCR |
| [Methods](https://github.com/duceppemo/primer-finder/wiki/Methods) | what each step does, the filtering rules, limits |
| [Outputs](https://github.com/duceppemo/primer-finder/wiki/Outputs) | every file and field |
| [Example](https://github.com/duceppemo/primer-finder/wiki/Example) | the simulated dataset and its expected result |
| [Validation](https://github.com/duceppemo/primer-finder/wiki/Validation) | finding four published qPCR assays in public genomes |
| [FAQ](https://github.com/duceppemo/primer-finder/wiki/FAQ) | troubleshooting, "no contig passed" |
| [Development](https://github.com/duceppemo/primer-finder/wiki/Development) | tests, continuous integration, releases |

## Citation

If primer-finder helped your work, please cite it — [doi:10.5281/zenodo.23226411](https://doi.org/10.5281/zenodo.23226411),
which always resolves to the latest version (see [CITATION.cff](https://github.com/duceppemo/primer-finder/blob/main/CITATION.cff)) — together with the
programs it runs: [KMC](https://github.com/refresh-bio/KMC),
[SKESA](https://github.com/ncbi/SKESA) or [SPAdes](https://github.com/ablab/spades),
[minimap2](https://github.com/lh3/minimap2), [BLAST](https://blast.ncbi.nlm.nih.gov/) and, for the design
step, [Primer3](https://primer3.org/) and [insilicoPCR](https://github.com/duceppemo/insilicoPCR).

## License

[MIT](https://github.com/duceppemo/primer-finder/blob/main/LICENSE)
