# 2026-11-27 — insilicoPCR 2.1.0: a site is accepted by its Tm, not by a mismatch count

insilicoPCR 2.0 (2026-11-26) replaced the mismatch budgets of 1.x (`-m`, `--probe-mismatches`,
`--5p-truncation`) with one rule. A primer or probe site is accepted when its melting temperature on the
template reaches `--anneal-temp` (60 °C) less `--tm-margin` (5 °C), and a primer's site also has to pass
the weighted 3' rule (`--3p-window`, `--3p-penalty`, unchanged). The Tm is Primer3's thermodynamic
alignment, worked out on the bases under the site, in the reaction `--na`, `--mg`, `--dntp` and
`--oligo-conc` describe. 2.1.0 (2026-11-27) changes FASTQ mode only, and fixes the count of 5' mismatches
in `report_no_product.tsv`.

primer-finder 1.13.0 could not run 2.x at all. It passes `-m`, which 2.0 refuses with exit 2, and it reads
`Gene` and `Forward5PrimeUnaligned`, which 1.0 renamed `Assay` and `Forward5PrimeUnpaired`. This record
asks what the Tm rule changes for primer-finder's verdicts, and says what primer-finder changed in
response.

**Answer.**

- **Requirement.** primer-finder now requires insilicoPCR 2.1.0. It passes the annealing temperature, the
  margin and the reaction Primer3 designed in, and no mismatch budget. `--mismatches` and
  `--probe-mismatches` are gone.
- **Verdicts.** On the *Xylella* panel (2,430 assays both arms hold, at the defaults) 117 PCR verdicts
  move `no` to `yes` and 8 `yes` to `no`; in qPCR mode 78 and 15. Each moved for a reason
  `report_no_product.tsv` or the Tm columns give, listed below.
- **Why most moved.** Most are exclusion genomes whose primer site carries one mismatch, or a few
  unpaired 5' bases. The budget of 1 let that site through; its Tm is below 55 °C.
- **The 3'-end check** keeps its meaning and flags far fewer assays: PCR 501 of 2,047 `yes` assays (24%)
  before, 146 of 2,156 (7%) now; qPCR 330 of 2,124 (16%), now 132 of 2,187 (6%).
- **The example.** Its 17 verdicts do not move.

primer-finder 1.13.0 (`v1.13.0`) with insilicoPCR 0.14.0, against the primer-finder of this record with
insilicoPCR 2.1.0. 0.14.0 gives byte-identical reports to 0.12.0 and 0.13.0
(`../2026-11-13_v1.12.0_insilico0130/`). Both are the linux-x64 releases, sha256-checked, run with Primer3
2.6.1 and blast 2.17.0+ in the `primer-finder_dev` environment.

## 1. What 2.0 and 2.1 changed for primer-finder

From insilicoPCR's CHANGELOG (1.0.0, 2.0.0 and 2.1.0), its release notes and its wiki (Reports,
Command-line Interface, Amplicon Calling, Why Tm, Primer FASTA Format):

| | 1.x and earlier | 2.1.0 | primer-finder now |
|---|---|---|---|
| site rule | mismatches within `-m` (probe: `--probe-mismatches`), up to 5 unpaired 5' bases free | Tm on the template ≥ `--anneal-temp` − `--tm-margin` | passes `--anneal-temp 60 --tm-margin 5` (new options), always explicitly |
| reaction | — | `--na --mg --dntp --oligo-conc` (defaults 50, 3, 0.8, 250) | passes `--monovalent --divalent --dntp --primer-conc`; the probe is judged at the primers' concentration, since insilicoPCR takes one |
| `-m`, `--probe-mismatches` | the budget | refused, exit 2 | never passed; refused with a message if given |
| 3' rule | `--3p-penalty`, default 0 | unchanged | strict runs at 0, the lenient 3'-end run at 4 and nothing else changed |
| assay column | `Gene` | `Assay` (1.0) | reads `Assay`; a report with `Gene` is refused as one from before 1.0 |
| 5' columns | `*5PrimeUnaligned` | `*5PrimeUnpaired` (1.0) | required, not read |
| `Pairing` | — | `designed`, `one-primer`, `two-assays` (1.0) | an inclusion genome counts on a `designed` product only; an exclusion genome on any |
| Tm columns | — | `ForwardTm`, `ReverseTm`, `ProbeTm` (2.0) | required; give `*_inclusion_weakest_tm` |
| samples | the log only | `report_samples.tsv`: every sample, `analysed` or `not-analysed` (1.0) | required; a `not-analysed` genome, or a count unlike the staged one, fails the check |
| no product | — | `report_no_product.tsv` (2.0): the failing oligo, its best site and its Tm | gives `*_inclusion_no_product` |
| the report | `report.tsv` or `report(n).tsv` | `Report written to <path>` on standard output (1.0) | read from there |
| exit codes | — | 0, 1, 2 for a command line it cannot use (1.0) | exit 2 is named |
| licence | MIT | GPL-3.0-or-later | run as a separate program, no code copied: primer-finder's licence is unchanged |

2.0.0 is refused as well as 1.x. Its products are the same as 2.1.0's (the 2.1.0 notes), but 2.1.0 is
what this record checks, and 2.0.0 counts a mismatch among a primer's first five bases twice in the
reasons primer-finder quotes.

## 2. The rule, measured one base at a time (`tm_zones.py`)

The same method as `mismatch_zones.py` in the earlier records, run on 2.1.0. An assay is placed in a
synthetic template, and the template is mutated at known distances from an oligo's end, one genome per
case. The oligos are the first assay primer-finder designs on its example, as Primer3 gave them:
`ACTACTGTCCCATGTGAACT` / `CACAGGCATTACGTTAACCG` (59.8 / 61.0 °C) and probe `TCGGAGTACTAGACCTCAGAATGT`
(63.5 °C). On their perfect sites insilicoPCR gives 59.7, 61.1 and 63.7 °C; sites need 55 °C.

Random oligos of 50-55% G+C were tried first and melt 5-7 °C hotter: on them every single mismatch looks
tolerable, which is not what a designed assay sees.

**The forward primer, one mismatch by its distance from the 3' end** (Tm of the site where it is
accepted):

| offset from the 3' end | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| strict (`--3p-penalty 0`) | no | no | no | no | no | yes, 56.9 | yes, 56.8 | no, 50.5 | no, 50.7 |
| lenient (`--3p-penalty 4`) | yes, 59.2 | yes, 57.6 | yes, 55.2 | yes, 56.0 | no | yes | yes | no | no |

**More cases:**

- **Mismatched 3' window.** Runs of 2-4 mismatches ending at the 3' end are refused by both passes
  (3' penalties 6 to 9). A terminal mismatch costs 3 under an A-ending primer and 4 under a G-ending one,
  as insilicoPCR's table says, so `--3p-penalty 4` is still the price of the worst single mismatch.
- **The 5' end is no longer free.**
  - Mismatches at the first 1-3 bases, and 1-3 unpaired bases together, still pass (59.5 to 56.2 °C).
  - A mismatch 4 or 5 bases in fails (54.3 and 52.5 °C), as do 5 unpaired bases (51.8 °C).
  - Under 1.x all of those were free up to 5 bases.
- **The probe.** One mismatch anywhere is accepted (56.9 to 63.3 °C). So are most pairs of mismatches
  (56.2 to 62.6 °C), but not one in the middle with another 3 bases from the 5' end (53.7 °C), nor three
  inside (50.5 °C).

**The 3'-end check of a design.** The same primer is run as a design's strict and lenient exclusion
passes, and the lenient report is read through `insilico.read_three_prime_held`:

- **What only the lenient pass amplifies.** 6 genomes, and every one carries a 3'-window imperfection.
  The pass differs from the strict one only by the 3' rule, as primer-finder now assumes.
- **One mismatch in the window.** Counted as held where the site still melts at the minimum: offsets 1-4,
  but not 5, whose site melts below 55 °C.
- **Two mismatches in one window.** Still readmitted at penalty 4 when they weigh no more than a terminal
  one, here at 2+4 and 3+4 (both 56.0 °C). They are not counted (`MOST_IN_WINDOW` = 1), so the rule is
  still needed.

## 3. The example

Its 17 assays get the same verdicts and the same counts in both modes, and the first assay of the table
stays. The weakest inclusion site of every assay is 58 °C or more. What moves is the held count:

- **Before.** 3 qPCR and 4 PCR assays were held off an exclusion genome by their 3' end.
- **Now.** None are. Under 0.14.0 the lenient run passed `-m 2`, which let through a primer site with two
  mismatches, one of them in the 3' window. Under 2.1.0 that site melts at 48.7 °C, and
  `report_no_product.tsv` says so.

## 4. The contig-edge set (`compare_edges.py`)

The 300 assays of `../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py` are run as each primer-finder
runs them. The set's oligos are stretches of *L. monocytogenes* taken as they come, not designed to a Tm.
Most of them melt below 55 °C even on their own perfect site: only 105 PCR assays (98 qPCR) have oligos
that all reach it. So 2.1.0 amplifies 163 (assay, genome) pairs in PCR mode where 0.14.0 amplified 560
(qPCR 146 against 529). Of the 398 pairs lost in PCR mode, 279 have no primer mismatch at all and fail
only because an oligo is "below the needed Tm even on a perfect match". The rest fail on a mismatch, or
on primer bases past a contig end, that the Tm does not survive. This is the rule doing what it says on oligos no one would design.

On the assays whose oligos are strong enough, it shows what the set was built for: a primer hanging off a
contig end.

| primer bases past a contig end (`lmono_cut`) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| assays (PCR) | 71 | 5 | 7 | 3 | 5 | 3 | 5 | 6 |
| amplified, 0.14.0 | 71 | 5 | 7 | 3 | 5 | 3 | 0 | 2 |
| amplified, 2.1.0 | 71 | 5 | 4 | 2 | 0 | 1 | 0 | 2 |

1.x let up to 5 such bases through for nothing (`--5p-truncation`). 2.x counts them against the site's
Tm: 1 base off is always survived here, 2-3 sometimes, 4-5 almost never. Two assays 7 bases off amplify
in both releases, which the plan's break does not explain; they were not examined further.

The held-only-by-the-3'-end pairs fall from 93 to 13 (PCR) and from 70 to 9 (qPCR), all 13 and 9 among
the old ones. The difference is the lenient run's `-m 2` of 1.13.0 again.

## 5. The *Xylella* panel (`run_arms.sh`, `compare_arms.py`)

The four subspecies' find results are designed by both arms at the defaults: 2,430 assays, the same
oligos in both. Pauca's Primer3 wobble did not occur this time. Every (assay, genome) pair that moved is
put to a reason. A pair the new arm no longer amplifies takes its reason from `report_no_product.tsv` of
the new run; a pair it now amplifies, from the new row's `Pairing`, mismatch counts and Tm columns.

| | verdicts moved | exclusion pairs lost | exclusion pairs gained | inclusion pairs lost | inclusion pairs gained |
|---|---|---|---|---|---|
| PCR | 117 `no` → `yes`, 8 `yes` → `no` | 374 | 44 | 14 | 0 |
| qPCR | 78 `no` → `yes`, 15 `yes` → `no` | 159 | 44 | 10 | 1 |

**Why exclusion genomes stopped amplifying** (the `no` → `yes` moves):

- **One primer mismatch.** In 341 of the 374 lost PCR pairs, and 128 of the 159 qPCR ones, a primer site
  in the exclusion genome carries one mismatch, which the old budget of 1 let through. Its Tm is now
  below 55 °C (`no site melts at the needed Tm`). For example, `Contig_2496_67.4597_assay18` in pauca
  De Donno had a forward site with one mismatch, whose best Tm is 53.8 °C.
- **No mismatch, but unpaired 5' bases.** In 29 PCR pairs (27 qPCR) the old product had no primer
  mismatch but 1, 2 or 5 unpaired 5' bases, free under 1.x. They now cost the site its Tm.
- **The sites no longer face each other.** In 4 pairs a primer's mismatched site is refused, and the
  genome's remaining sites do not face each other within 10 kb.

**Why exclusion genomes started amplifying** (most of the `yes` → `no` moves):

- **A primer alone, binding both strands** (`Pairing` `one-primer`). 32 PCR pairs, 16 qPCR, among them
  `Contig_2344_86.6818_assay3` in morus MUL0034 and Riv16, at 58.5 °C. insilicoPCR has looked for such
  products since 0.15.3. It is a real product of the tube, so it counts on the exclusion group.
- **A probe site with 2 or 3 mismatches** that still melts at 55 °C or more (16 qPCR pairs, the probe
  at 55.4 to 65.0 °C). The old probe budget of 1 refused it. This is the Why-Tm argument: a probe designed hotter
  than its primers binds through a couple of mismatches.
- **A primer site with one mismatch and 5 unpaired 5' bases** (12 pairs in each mode, every one a
  reverse primer, every product of an unexpected size), which 1.x refused on the two together. One example is `Contig_2395_61.3198_assay6` on pauca: its reverse primer's 14
  paired 3' bases melt at 58.6 °C and make a 62-base off-target product (`SizeCheck` `unexpected`).

**Why inclusion genomes were lost** (the other `yes` → `no` moves): 14 PCR pairs, 10 qPCR. A primer whose
inclusion site carries a mismatch now melts below 55 °C. For example, `Contig_1757_60.384_assay5`'s
forward primer melts at 52.5 °C in morus Riv19 and Riv25. `*_inclusion_no_product` names the oligo and
the Tm. These are findings about the assays, and no setting was changed to keep their old verdict.

**The 3'-end check.** It still asks one question, whether an exclusion genome is held off only by one
mismatch in a primer's last five bases, and it flags far fewer assays:

| `yes` assays flagged | 1.13.0 + 0.14.0 | this + 2.1.0 |
|---|---|---|
| PCR | 501 of 2,047 (24%) | 146 of 2,156 (7%) |
| qPCR | 330 of 2,124 (16%) | 132 of 2,187 (6%) |

The reason is the Tm rule's:

- **The old lenient run was lenient twice.** It passed `-m 2` to pay for the 3'-window mismatch out of
  the count, so a primer with a second mismatch further in was readmitted too. The Tm usually refuses
  that site.
- **The evidence.** All 16 of the example's lenient PCR rows under 0.14.0 are such sites, a primer with
  two mismatches of which one is in the 3' window (section 3), and
  on the contig-edge set the held pairs fall from 93 to 13 (section 4).

**How close the new `yes` assays are to failing.** By their weakest inclusion site:

| | under 56 °C | 56-58 °C | 58-60 °C | 60 °C or more |
|---|---|---|---|---|
| PCR | 0 | 7 | 1,602 | 547 |
| qPCR | 0 | 7 | 1,626 | 554 |

None rest on a site within 1 °C of the minimum. Primer3 designs the primers at 58-63 °C on the sequence
every inclusion genome shares, so their inclusion sites are perfect.

**The table's order.** The first assay of every target is the same in both arms.

## What this does and does not show

- It shows primer-finder had to change to run 2.x at all, and that the Tm rule moves verdicts both ways
  for reasons the reports name: 125 PCR and 93 qPCR verdicts out of 2,430 assays, each pair put to a rule.
- It shows the 3'-end check still means what it says: the lenient pass differs from the strict one by
  the 3' rule alone (section 2).
- It shows the new `inclusion_no_product` and `inclusion_weakest_tm` columns say why an assay misses an
  inclusion genome, and how close it comes to doing so.
- It does not show which rule a polymerase agrees with. The Tm rule is insilicoPCR's, argued on its wiki
  (Why Tm) from data of its own. Whether an exclusion genome whose primer site melts at 53 °C stays silent
  at 60 °C is a question for the bench.
- It is one organism, at one annealing temperature and margin, with the reaction primer-finder designs
  for.
- It does not cover FASTQ mode, which primer-finder does not use.

## Reproducing

```bash
# kmc, skesa, minimap2, blast and primer3_core on PATH (the primer-finder_dev environment), the NCBI
# datasets tool for the downloads, and the releases extracted. Scripts that run insilicoPCR 2.1.0 use the
# primer-finder of this record; the old arm is the v1.13.0 tag (git worktree add <dir> v1.13.0).
python tm_zones.py /path/to/insilicoPCR-2.1.0 work_zones > tm_zones.txt

python ../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py lmono_EGDe.fasta linnocua_Clip11262.fasta <out>/edges
python compare_edges.py <out>/edges /path/to/insilicoPCR-0.14.0 /path/to/insilicoPCR-2.1.0

bash run_arms.sh <primer-finder-v1.13.0> /path/to/insilicoPCR-0.14.0 /path/to/insilicoPCR-2.1.0 <panel_work> <out>/arms 6
python compare_arms.py <out>/arms
```

`tm_zones.py` takes about 2 minutes and `compare_edges.py` about 15, both arms included. The ten designs
of `run_arms.sh` run side by side in about 30 minutes on 64 cores and take about 19 GB: 18 GB for the
old arm, almost all of it insilicoPCR 0.14.0's `insilico_pcr/*/detailed_report/` (3.3-5.8 GB a target),
and 0.4 GB for the new arm, whose 2.1.0 writes far less there. `amplicons/` is about 31 MB a target. The panel's work folder is rebuilt from NCBI in about 5 minutes when it is missing
(`../2026-11-13_v1.12.0_insilico0130/run_panel.sh` describes it).
