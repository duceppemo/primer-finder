# Validation

Two different things are checked, and they answer different questions.

**The tests and the example** (`pytest`, `example/run_example.sh`) check that the code does what the
documentation says: 182 tests with stand-ins for the external programs, and a simulated dataset whose
specific region is known, down to which bases must come out in lower case. They run in seconds and the CI
runs them on every push. See [Development](Development) and [Example](Example).

**The validation suite** (`validation/`) checks something the tests cannot: that primer-finder, run on public
genomes, finds the regions that published diagnostic assays were actually designed on. Records of each run
live in `validation/results/<date>_v<version>/` and are never rewritten.

## Four published *Xylella fastidiosa* subspecies assays

The answer key is:

> Dupas E., Briand M., Jacques M.-A., Cesbron S. (2019) *Novel tetraplex quantitative PCR assays for
> simultaneous detection and identification of Xylella fastidiosa subspecies in plant tissues.*
> Frontiers in Plant Science 10:1732. [doi:10.3389/fpls.2019.01732](https://doi.org/10.3389/fpls.2019.01732)

Their primers and probes were designed with SkIf, a kmer-based signature tool, on 58 *X. fastidiosa* genomes
split into an ingroup and an outgroup — the same question primer-finder answers, asked with different code.
The paper publishes the oligo sequences and the regions they sit in.

The test takes 25 complete RefSeq genomes covering five subspecies and, for each of four subspecies, runs
primer-finder with that subspecies as the inclusion group and every other subspecies as the exclusion group,
with default options. An assay counts as recovered when both primers and the probe are found exactly, in one
reported region, with the primers facing each other, the probe between them, and the amplicon the published
length.

### Result, primer-finder 1.1.0

| Assay | Inclusion group | Final regions | Recovered | Rank of its region |
|---|---|---|---|---|
| XFM | subsp. *multiplex* (8 genomes) | 711 | yes | 1 |
| XFF | subsp. *fastidiosa* (6 genomes) | 377 | yes | 5 |
| XFP | subsp. *pauca* (4 genomes) | 1,413 | yes | 105 |
| XFMO | subsp. *morus* (5 genomes) | 286 | yes | 5 |

Every amplicon came out exactly the published length, and each recovered region aligned to the paper's
reference genome at 100% identity. The full record, including how far each region overlaps the one the paper
reports and what the comparison does and does not prove, is in
[`validation/results/2026-10-08_v1.1.0/SUMMARY.md`](https://github.com/duceppemo/primer-finder/blob/master/validation/results/2026-10-08_v1.1.0/SUMMARY.md). Earlier records are kept beside it.

What it does not show: that the other hundreds of regions in each run would make working assays, or anything
about specificity outside those 22 genomes. Those regions are candidates to design on, which is what the
[Methods](Methods) page means by a candidate.

### Running it

```bash
bash validation/xylella/run.sh /path/to/work 16 32
```

It needs primer-finder and its programs, plus the NCBI datasets command line tool
(`conda install -c conda-forge ncbi-datasets-cli`), downloads the genomes once, and exits non-zero if any of
the four assays is not recovered. A few minutes, most of it downloading.
