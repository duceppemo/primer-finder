# Installation

primer-finder is a Python package with no Python dependency; what it needs are the programs it runs: KMC,
SKESA or SPAdes, minimap2 and BLAST. Python 3.10 or later.

## From bioconda (recommended)

> **Not available yet**: the recipe is awaiting review ([bioconda-recipes#70034](https://github.com/bioconda/bioconda-recipes/pull/70034)). Until it is merged,
> install [from PyPI](#from-pypi) into an environment that already holds the programs, or
> [from the source code](#from-the-source-code) with `environment.yml`, which brings them.

```bash
conda create -n primer-finder -c conda-forge -c bioconda primer-finder
conda activate primer-finder
primer-finder --version
```

Mamba works the same way (`mamba create ...`) and is much faster.

## From PyPI

`pip install primer-finder` installs the command, but **not** the programs it runs: KMC, SKESA or SPAdes,
minimap2 and BLAST are not Python packages. Use it inside an environment that already holds them (or install
them separately); otherwise prefer the bioconda package above, which brings everything.

```bash
pip install primer-finder   # or: pipx install primer-finder
```

## From the source code

The conda environment file brings the external programs in, then pip installs primer-finder itself:

```bash
git clone https://github.com/duceppemo/primer-finder
cd primer-finder
conda env create -f environment.yml
conda activate primer-finder
pip install .
```

Without installing anything, `python primer_finder.py ...` works from a clone and takes the same options as
`primer-finder`.

## The programs it runs

| Program | Used for | conda package |
|---|---|---|
| `kmc`, `kmc_tools` | counting and subtracting kmers | `kmc` >=3.2 |
| `skesa` | assembling the kmers (`-a skesa`, the default) | `skesa` >=2.4 |
| `spades.py` | assembling the kmers (`-a spades`) | `spades` >=3.15 |
| `minimap2` | mapping the contigs to an exclusion genome | `minimap2` >=2.24 |
| `blastn`, `makeblastdb` | checking the contigs against every genome | `blast` >=2.14 |
| `primer3_core` | designing assays on the regions (`primer-finder design`) | `primer3` >=2.6 |

A missing program is reported before anything runs, with the conda package that provides it. Only the
assembler you ask for has to be installed.

[insilicoPCR](https://github.com/duceppemo/insilicoPCR) is **not** installed with primer-finder and is not
needed to run it. It is a separate program, and `primer-finder design --insilico-pcr` points at an extracted
portable release of it (which brings its own Java, BBMap and BLAST+). See
[Designing assays](Designing-assays).

## Checking the installation

The repository holds a small example: simulated genomes whose specific region is known, which the whole
pipeline must find. It takes a few seconds.

```bash
# In a clone
bash example/run_example.sh

# Or, with primer-finder installed from conda
curl -sL https://github.com/duceppemo/primer-finder/archive/refs/tags/v1.4.0.tar.gz | tar -xz --strip-components=1 primer-finder-1.4.0/example
bash example/run_example.sh
```

It ends with `All checks passed`. See [Example](Example) for what it checks.

## Upgrading from the 2022 scripts

The two scripts became a package, and the old command lines still work:

| Before | Now |
|---|---|
| `python primer_finder.py -i incl/ -e excl/ -o out/` | `primer-finder -i incl/ -e excl/ -o out/` (or the same `primer_finder.py` command) |
| `python IDT_results_converter.py sheet.xlsx out.fasta prefix` | `primer-finder idt sheet.xlsx out.fasta prefix` |

The options are unchanged. What moved is where the results are written: `final_kmers.fasta` is still at the
top of the output folder, and the intermediate files are now in numbered subfolders
(see [Outputs](Outputs)). pysam, Biopython, pandas, psutil, pyahocorasick and xlrd are no longer needed, and
neither are bowtie2 and samtools; `requirements.txt` is gone, replaced by `environment.yml`.
