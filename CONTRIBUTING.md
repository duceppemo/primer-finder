# Contributing

Thanks for helping improve primer-finder. Bug reports, assay results and documentation fixes are as welcome
as code.

## Reporting a bug

Open an [issue](https://github.com/duceppemo/primer-finder/issues/new/choose) with the exact command and the
`primer_finder.log` and `run_info.json` files from the output folder: they record the version of
primer-finder and of every program it ran. Rerunning with `--debug --keep-intermediate` keeps the files of
each step, which usually shows where it went wrong.

## Changing code

```bash
git clone https://github.com/duceppemo/primer-finder && cd primer-finder
conda env create -f environment.yml && conda activate primer-finder
pip install -e ".[test]"
pre-commit install
```

- Keep `ruff check .` and `pytest` green; add a test for every fix or feature. The tests replace the external
  programs with stubs (`tests/stub_programs.py`), so the whole suite runs in seconds without conda.
- primer-finder itself uses the Python standard library only; external programs are called on the command
  line. Please do not add a Python dependency.
- A new or replacement program must be maintained, installable from bioconda together with the others, and
  compared with the current one on real data before it replaces it.
- User-visible changes go in `CHANGELOG.md`; options and outputs are documented in `docs/wiki/`.
- `bash example/run_example.sh` must still end with "All checks passed".
- Open the pull request against `main`.
