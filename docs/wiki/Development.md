# Development

```bash
git clone https://github.com/duceppemo/primer-finder && cd primer-finder
conda env create -f environment.yml && conda activate primer-finder
pip install -e ".[test]"
pre-commit install
```

## Layout

```
primer_finder/
├── cli.py        the command line: `find` (the default), `design` and `idt`
├── pipeline.py   the steps of `find`, in order, and what each writes
├── kmers.py      KMC: counting, subtracting, dumping
├── assemble.py   SKESA and SPAdes
├── mapping.py    minimap2, cigar strings, which contigs are candidates
├── blast.py      makeblastdb/blastn, presence, and what differs where
├── design.py     the steps of `design`: Primer3, the amplicon checks, the order
├── insilico.py   the primer files insilicoPCR reads, running it, and its report
├── idt.py        IDT order sheets (.xlsx without any library, .csv, .tsv)
├── seqio.py      fasta reading and writing, gzipped or not
├── tools.py      running external programs, their versions, and their time limits
└── system.py     how many CPUs and how much memory may be used
example/          the simulated dataset, the runner and the checks
validation/       the published assays this is checked against, and dated records
tests/            the test suite, with stand-ins for the external programs
docs/wiki/        the sources of this wiki
recipe/           a copy of the bioconda recipe
```

primer-finder uses the **Python standard library only**. Everything else is an external program called on the
command line, which keeps the conda environment small and the tests fast. Please keep it that way.

## Tests

```bash
pytest -q                                  # the whole suite, a few seconds
pytest --cov=primer_finder --cov-report=term-missing
pytest -m tools                            # only the example, with the real programs
ruff check .
```

`tests/stub_programs.py` holds a stand-in for every external program, and `tests/stub_insilico_pcr.py` one
for insilicoPCR; the `stubs` fixture puts them first on `PATH` and lets each test say what they should
produce (the kmers KMC dumps, the contigs the assembler writes, the SAM records, the blast hits, the assays
Primer3 proposes, which samples amplify). The whole pipeline is therefore tested without conda, in seconds,
including its error paths. The test marked `tools` runs `example/run_example.sh` with the real programs and
is skipped when they are not installed.

Add a test with every fix, and check that it fails without the fix.

## Continuous integration

`.github/workflows/`:

| Workflow | When | What |
|---|---|---|
| `ci.yml` | push to `main`, pull requests | `ruff check`, the test suite on Python 3.10-3.13 and on macOS, coverage to Codecov, and the example with the real programs in a micromamba environment built from `environment.yml` |
| `wiki.yml` | push touching `docs/wiki/**` | copies `docs/wiki/` to the GitHub wiki |
| `release.yml` | a `v*` tag | builds the wheel and the sdist, checks that the tag matches the version, and creates the GitHub release with the notes taken from `CHANGELOG.md` |
| `publish.yml` | a `v*` tag, or by hand with a tag | builds the wheel and the sdist again and uploads them to PyPI through a trusted publisher (OpenID Connect, no API token). It runs from the tag rather than from the release, because a release created by `release.yml` cannot trigger another workflow. The `pypi` environment of this repository must match the environment of the publisher registered on PyPI |

The wiki is maintained in `docs/wiki/` in this repository and published by `wiki.yml`. Do not edit the pages
on GitHub directly: the next push overwrites them.

## Releasing

1. Bump the version in `pyproject.toml`, `primer_finder/__init__.py`, `CITATION.cff` (`version` and
   `date-released`) and `recipe/meta.yaml`; add the section to `CHANGELOG.md`. `tests/test_cli.py` checks
   that all five agree.
2. `pytest` and `bash example/run_example.sh` green; update the tarball version in the README and in
   `docs/wiki/Installation.md`.
3. Tag and push: `git tag v1.0.1 && git push --tags`. `release.yml` publishes the release, and
   `publish.yml` then uploads the wheel and the sdist to PyPI.
4. Zenodo archives the release and mints a DOI; add the new version DOI to `CITATION.cff`.
5. Update the bioconda recipe: bump `version`, recompute the sha256 of the release tarball
   (`curl -sL <url> | sha256sum`), reset `number` to 0, and open a pull request against
   [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) (its autobump bot often does it first).
   Keep `recipe/meta.yaml` in this repository the same as the one there.
