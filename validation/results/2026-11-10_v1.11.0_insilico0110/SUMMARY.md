# 2026-11-10 — primer-finder 1.11.0 + insilicoPCR 0.11.0: nothing primer-finder reads moves

insilicoPCR 0.11.0 (2026-11-10) is mostly about FASTQ mode and reporting, like the three releases
before it, but it also changes one thing in assembly mode, the only mode primer-finder uses: the size it
reports for a product whose primer hangs off the start of a contig. This record asks whether
primer-finder 1.11.0 needs anything changed to run on 0.11.0, and answers it three ways: by reading what
the release changed, by designing the example under 0.10.1 and 0.11.0, and by building the changed case
on purpose and running both releases on it the way a design does.

**Answer: no change is needed.** The example's assays, primer files and six reports are byte-identical
under the two releases. On genomes whose contigs break inside primer sites, the two releases report the
same products, row for row; 0.11.0 gives 63 to 68 of them per report a size 1 to 5 bases longer, which
is the correction it announces, and primer-finder never reads the size. Which samples each assay
amplifies, and which only at the 3' end, are the same.

primer-finder 1.11.0 (`bb6a236`), Primer3 2.6.1, blast 2.17.0+, insilicoPCR 0.10.1 and 0.11.0 (the
`v0.11.0` linux-x64 release, sha256-checked), every arm in the same `primer-finder_dev` environment. The
two releases bundle the same runtime, file for file: blast 2.17.0+, BBTools and JDK.

## 1. What the release changed

From insilicoPCR's changelog and `git diff v0.10.1 v0.11.0`:

- **FASTQ mode**: read-level calls are made per product and count only the reads at the product's
  locus, the search after blast went from hours to seconds, and the diagnostics say more. primer-finder
  hands insilicoPCR plain fasta genomes only, so none of it runs on its behalf.
- **Product-size labels**: a primer header may carry `size=` and `other=` after the primer ID, and the
  report then ends with a `SizeCheck` column. Rows are labelled, never removed, and the column is absent
  when no assay declares a size. A header's ID now ends at the first blank.
- **HTML report**: every run writes `consolidated_report/report.html` beside `report.tsv`, and prints
  one more line saying so.
- **Product size at a contig start** (assembly mode too): when up to `--5p-truncation` 5' bases of a
  primer lie before the first base of a contig, `AmpliconSize` now counts them, as it already did at
  the other end of a contig. 0.10.1 clamped the product at the contig's first base and reported it
  short by those bases. `GenomeLocation` is unchanged, and `--min-length`/`-l` apply to the full size.
- The lookup of hits while pairing primers was rewritten for speed, with the same output.

What that means for what primer-finder relies on:

| what primer-finder relies on | 0.10.1 → 0.11.0 |
|---|---|
| `-i -p -o -t -m -l --probe-mismatches --3p-penalty` | unchanged; 0.11.0 adds no option |
| `-v` prints `insilicoPCR <x.y.z>` | unchanged |
| the primer files it writes (`>name-F`, no description) | read as before; with nothing declared there is no `SizeCheck` column |
| `consolidated_report/report.tsv` columns | unchanged |
| which hits are accepted (`-m`, the 3' window, `--5p-truncation`) | unchanged: the hit filters are not in the diff |
| which pairs of hits are a product | unchanged, except that the size held against `-l 10000` is the full one: a product within 5 bases of 10,000 at a contig start could now be refused. Nothing primer-finder designs is near that. |
| `AmpliconSize` | larger by 1-5 at a contig start; **not read** by primer-finder, which counts rows per sample and assay and reads the four mismatch columns |
| `report(n).tsv` beside an existing report | unchanged; `report.html` is new in the folder and is overwritten by a rerun |
| bundled blast 2.17.0+ | unchanged |

## 2. The example, designed under both releases

`compare_example.sh` runs the example search, designs its region under 0.10.1 and under 0.11.0 (17
assays, each run making six insilicoPCR passes over four inclusion and four exclusion genomes), and
compares what the in silico check wrote. Under 0.11.0 it then designs again into the same folder, and
runs the `run_insilico_pcr.sh` the design wrote.

| compared | files | result |
|---|---|---|
| 0.10.1 against 0.11.0: `assays.tsv`, both primer fastas | 3 | byte-identical |
| 0.10.1 against 0.11.0: the six `report.tsv` (164 rows in all) | 6 | byte-identical |
| 0.11.0 designed twice into one folder: `assays.tsv` and the six reports | 7 | byte-identical; the first report is set aside as `report.tsv.previous`, no `report(n)` file is left |
| 0.11.0 through `run_insilico_pcr.sh`: the six reports | 6 | byte-identical to the design's own |

The example has no contig that starts inside a primer site, so it does not meet the one assembly-mode
change. Section 3 does.

## 3. Contigs that break inside primer sites

`make_contig_edges.py` takes 300 assays from the first 600 kb of *Listeria monocytogenes* EGD-e (a
20-base primer pair 80-150 bases apart, a 24-base probe between them) and writes three genomes: the
region broken into 101 contigs, 50 of the breaks 1-7 bases inside a forward primer and 50 of them 1-7
bases before the 5' end of a reverse primer (`lmono_cut`); the same region with one base in 45
substituted and half of the breaks (`lmono_mutated`); and the first 600 kb of *L. innocua* Clip11262.
`compare_releases.py` runs both releases on them through primer-finder's own `insilico` module, at the
design defaults: the strict pass (`-m 1 -l 10000`) read with `read_report`, and the lenient 3'-end pass
(`-m 2 --probe-mismatches 1 --3p-penalty 4`) read with `read_three_prime_held`.

What primer-finder reads:

| primer file | amplified (assay, sample) pairs, 0.10.1 / 0.11.0 | held only by the 3' end, 0.10.1 / 0.11.0 |
|---|---|---|
| `assays_pcr.fasta` | 560 / 560, the same pairs | 37 / 37, the same pairs |
| `assays_qpcr.fasta` | 516 / 516, the same pairs | 30 / 30, the same pairs |

The reports themselves:

| report | rows (both releases) | rows that differ | differing in |
|---|---|---|---|
| PCR, strict | 586 | 65 | `AmpliconSize` only |
| PCR, lenient | 656 | 68 | `AmpliconSize` only |
| qPCR, strict | 544 | 63 | `AmpliconSize` only |
| qPCR, lenient | 598 | 66 | `AmpliconSize` only |

Every differing row is a product that starts at the first base of its contig and is 1 to 5 bases longer
under 0.11.0. All but one per lenient report are the product the assay was built for, and there 0.11.0
gives exactly the size it was built with. The exception is an off-target product of 8.8 kb in
`lmono_mutated` that starts at one of the broken sites: 8,806 bases under 0.10.1, 8,807 under 0.11.0.

`site_table.py` reads the PCR strict reports for `lmono_cut`, site by site (assays reported at the size
they were built with / only shorter / not at all):

| site | bases of the primer off the contig | assays | 0.10.1 | 0.11.0 |
|---|---|---|---|---|
| unbroken | 0 | 200 | 200 / 0 / 0 | 200 / 0 / 0 |
| contig starts inside the forward primer | 1-5 | 36 | 0 / 36 / 0 | 36 / 0 / 0 |
| contig starts inside the forward primer | 6-7 | 14 | 1 / 0 / 13 | 1 / 0 / 13 |
| contig ends inside the reverse primer | 1-5 | 36 | 36 / 0 / 0 | 36 / 0 / 0 |
| contig ends inside the reverse primer | 6-7 | 14 | 1 / 0 / 13 | 1 / 0 / 13 |

So the limit of five unaligned 5' bases holds in both releases and at both ends of a contig; the only
thing that moved is the size at a contig start, from short to right. The two assays reported with six
or seven bases off are not exceptions: their sites sit in a stretch of the genome that occurs twice,
6 kb apart, and they amplify from the copy that was not broken.

The example took 24 seconds, the contig-edge comparison 7 minutes (eight insilicoPCR runs on 8 threads).

## What this does and does not show

- It shows primer-finder 1.11.0 runs on insilicoPCR 0.11.0 unchanged, and that 0.11.0 reports the same
  amplifications as 0.10.1 on assemblies, including where a contig breaks inside a primer site. The
  verdicts of the earlier records hold as measured, and none needs an erratum.
- It does **not** rerun the *Xylella* panel: the work folder of the earlier records no longer exists,
  and the example is 17 assays on eight synthetic genomes. What stands in for the panel here is the
  reading of the diff, which leaves the hit filters untouched, and the contig-edge set, which meets the
  one changed case in 63 to 68 rows of each report, on real sequence. The contig-edge set is 600 kb of each genome and not the
  whole of it: with 1,500 assays on the whole genomes one insilicoPCR pass had not finished after ten
  minutes.
- It does not cover insilicoPCR's FASTQ mode, which primer-finder does not use, and says nothing of
  releases after 0.11.0.

## Reproducing

```bash
# kmc, skesa, minimap2, blast and primer3_core on PATH; both insilicoPCR releases extracted
bash compare_example.sh <out>/example /path/to/insilicoPCR-0.10.1 /path/to/insilicoPCR-0.11.0 8

# the two Listeria genomes: GCF_000196035.1 (EGD-e) and GCF_000195795.1 (Clip11262), unzipped
python make_contig_edges.py lmono_EGDe.fasta linnocua_Clip11262.fasta <out>/edges
python compare_releases.py <out>/edges /path/to/insilicoPCR-0.10.1 /path/to/insilicoPCR-0.11.0
python site_table.py <out>/edges
```

Both comparisons exit 1 if anything differs that this record does not describe.
