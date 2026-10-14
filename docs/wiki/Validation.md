# Validation

Two different things are checked, and they answer different questions.

**The tests and the example** (`pytest`, `example/run_example.sh`) check that the code does what the
documentation says: 380 tests with stand-ins for the external programs, and a simulated dataset whose
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
[`validation/results/2026-10-08_v1.1.0/SUMMARY.md`](https://github.com/duceppemo/primer-finder/blob/main/validation/results/2026-10-08_v1.1.0/SUMMARY.md). Earlier records are kept beside it.

What it does not show: that the other hundreds of regions in each run would make working assays, or anything
about specificity outside those 25 genomes. Those regions are candidates to design on, which is what the
[Methods](Methods) page means by a candidate.

## Designing assays on those regions

`validation/results/2026-10-08_v1.2.0/SUMMARY.md` records a run of
[`primer-finder design`](Designing-assays) on the subsp. *multiplex* regions: 20 regions, 432 assays, 417
that could tell the groups apart, and in silico PCR of all of them against the 25 genomes in both modes.

That record is a dry run of the machinery, not evidence about PCR: nothing in it has been tested in a
laboratory, and in silico PCR partly restates its own matching rule. It is kept because it shows the command
working end to end on real data, and because the difference between its two modes is informative — assays
whose only differences sit under the probe are selective in qPCR mode and not in standard PCR mode.

### How many mismatches the check should allow

`validation/results/2026-10-08_v1.3.0/SUMMARY.md` asks what insilicoPCR's `-m` tolerance does, by measuring
it: the same 465 assays scored at 0, 1, 2 and 3, plus a controlled test with a synthetic template mutated
one base at a time.

The measured result is that a primer has three zones. The last two bases are **free** — a mismatch there is
not counted at any tolerance, because blast trims it off the alignment. Three and four bases in **never**
amplify, at any tolerance. Five or more bases in is what `-m` governs, counted per primer. So the tolerance
cannot reward a difference at the very 3' end, which is what the ranking cares most about, and one sentence
of the 1.2.0 page — that assays resting on a terminal run were selective in none of the cases seen — turned
out to be the alignment's behaviour rather than anything about PCR. It has been withdrawn as evidence.

The default is now `--mismatches 1`: 0 assumes any internal mismatch stops a primer, which hid 12 of the 30
absence-based assays that amplify something off-target in an exclusion genome; 2 assumes two mismatches per
primer still bind, which is what a deliberate allele-specific design uses to discriminate, so it is offered
as a stricter second pass; 3 measured nothing 2 did not.

## Thirty published *Listeria monocytogenes* clonal-complex assays

The *Xylella* set asks the right question of small groups: 8 inclusion genomes against 17. A second answer
key asks it of much larger ones, and of a harder target — clonal complexes *within* one species:

> Félix B. et al. (2023) *Identification by High-Throughput Real-Time PCR of 30 Major Circulating Listeria
> monocytogenes Clonal Complexes in Europe.* Microbiology Spectrum 11(3):e03954-22.
> [doi:10.1128/spectrum.03954-22](https://doi.org/10.1128/spectrum.03954-22)

Their 34 TaqMan primer-and-probe sets were designed on kmers found in 954 genomes, at least 15 per clonal
complex, and checked in silico against 2,388 more. The panel here is every complete RefSeq
*L. monocytogenes* assembly — 784 genomes, 763 of which the Institut Pasteur MLST scheme places in one of 79
clonal complexes — which gives inclusion groups of 18 to 86 genomes against exclusion groups of 677 to 745.

### Result, primer-finder 1.6.0

Nine of the 34 assays target a whole clonal complex that has at least 15 complete genomes. (CC1, CC14, CC37
and CC121 have two assays each, resolving subdivisions *within* the complex, so a complex-wide group is not
what they detect.)

These are what 1.6.0 reports; the region counts and CC5's rank have moved since the first run, as the
records below set out.

| Target | Inclusion | Exclusion | Regions | Recovered | Rank |
|---|---|---|---|---|---|
| CC4 (at `-p 0.94`) | 19 | 744 | 10 | yes | 1 |
| CC5 | 22 | 741 | 57 | yes | 10 |
| CC6 | 36 | 727 | 53 | yes | 5 |
| CC7 | 47 | 716 | 60 | yes | 1 |
| CC8 | 52 | 711 | 62 | yes | 2 |
| CC2, CC3, CC9, CC224 | 18–86 | 677–745 | 1–40 | no | — |

**Five of the nine were recovered**, four at the defaults and CC4 once `-p 0.94` allowed one of its 19
genomes to lack the region. In all four that were not, an independent blast census of the panel shows the
published assay is not both present in every genome of its own complex and absent from every other — so
primer-finder's rule is right to decline them, and for CC3 and CC224 the cross-reaction it declines on is
one the paper itself predicts in silico.

The full record, including the two limitations this panel exposed, is in
[`validation/results/2026-10-10_v1.5.0/SUMMARY.md`](https://github.com/duceppemo/primer-finder/blob/main/validation/results/2026-10-10_v1.5.0/SUMMARY.md).
The first of them — an inclusion-specific stretch shorter than 200 bases was never reported, because that is
where both assemblers stop of their own accord — has since been **fixed**: the floor is now the kmer size.
[`2026-10-13_v1.6.0/SUMMARY.md`](https://github.com/duceppemo/primer-finder/blob/main/validation/results/2026-10-13_v1.6.0/SUMMARY.md)
re-runs the nine targets with it removed. It puts 6 to 10 times as many stretches in front of the filters
and 2 to 3 times as many regions in the output, changes no verdict here, and turns one run that failed
outright into one that completes. The second — a published amplicon may straddle the edge of such a
stretch — stands.

It needs one more environment, because a clonal complex is not in the NCBI metadata and has to come from the
MLST profile:

```bash
conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
conda activate primer-finder_listeria
bash validation/listeria/run.sh /path/to/work 32 64
```

### Running it

```bash
bash validation/xylella/run.sh /path/to/work 16 32
```

It needs primer-finder and its programs, plus the NCBI datasets command line tool
(`conda install -c conda-forge ncbi-datasets-cli`), downloads the genomes once, and exits non-zero if any of
the four assays is not recovered. A few minutes, most of it downloading.
