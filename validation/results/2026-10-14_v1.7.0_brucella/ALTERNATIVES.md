# Could primer-finder 1.7.0 propose *better* Brucella assays from the same genomes?

2026-10-14, primer-finder 1.7.0, same panel as [`SUMMARY.md`](SUMMARY.md) (1,376 classical *Brucella*,
315 non-classical relatives, 103 held out). That record asked whether the rewrite finds the **published**
assays of Hoover et al. 2026. This one asks the other question: given the same genomes, what does it
propose of its own, and is any of it better?

**Everything below is in silico.** The published assays carry wet-lab validation on 235 classical
*Brucella* and 192 non-target strains. The candidates here have none, and no number in this file
substitutes for that. What is compared is sequence coverage, sequence specificity, and thermodynamic
prediction under one assumed reaction.

## How both sets were measured, on one yardstick

`design` reports Tm and structure for the assays it proposes, which makes those columns useless for
comparison on their own — the published assays were designed elsewhere and have no such row. So:

- **`score_assays.py`** blasts every oligo of every assay, published and proposed alike, against all 1,794
  panel genomes and counts only full-length, mismatch-free hits. Same rule as the census in `SUMMARY.md`;
- **`thermo.py`** runs Primer3's `check_primers` on every assay with **primer-finder's own default
  `Conditions`** (50 mM monovalent, 3 mM divalent, 0.8 mM dNTP, 250 nM primer, 200 nM probe), each oligo
  in a request of its own with the structure limits opened right up, so that a value over a limit is
  reported rather than causing a rejection.

`design` was run on the regions each `find` run reported: 398 regions for the genus target at the default
`-p 1.0`, and the corrected-panel runs for the species targets. 520, 15 and 98 assays came back for
genus, *B. abortus* and *B. melitensis*. The shortlist scored here is the best-ranked assay of each
distinct region, 12 + 3 + 11.

## *B. suis*: nothing to compare

No region, at any setting tried — `-p` 1.0, 0.95, 0.90 and 0.85, and `-k` 99, 71, 51 and 31. The reason is
phylogenetic rather than technical: **B. canis sits inside B. suis**, and the exclusion group of this
target holds 74 *B. canis* genomes, so a stretch present in every *B. suis* and absent from every
exclusion genome does not exist on this panel. The published BS assay is in the same position — the census
shows it is carried by one *B. canis*-labelled assembly — and it exists because its design did not demand
absence from *B. canis*.

To get a *B. suis* assay out of primer-finder, the groups have to say what is wanted: either put *B. canis*
in the inclusion group (an assay for *B. suis* sensu lato, which is what the phylogeny supports), or leave
it out of both groups and accept that the assay says nothing about *B. canis*.

## The genus assay: five candidates reach 100% where the published one reaches 99.85%

| assay | bp | panel species | exclusion | held out: sp./pinni/amazo/vulpis/inopi | hairpins F/R/P | Tm F/R/P |
|---|---|---|---|---|---|---|
| **Bru_genus** (published) | 93 | 1,374/1,376 = **99.85%** | 0 | 49 / 8 / 4 / **2** / 0 = 63 | 0.0 / 0.0 / 0.0 | 62.0 / 62.2 / 68.3 |
| genus:Contig_193_65.5561_assay0 | 101 | **1,376/1,376 = 100%** | 0 | 39 / 8 / 4 / 0 / 0 = 51 | **0.0 / 0.0 / 0.0** | 60.0 / 59.9 / 68.0 |
| genus:Contig_228_33.6114_assay0 | 125 | **100%** | 0 | 44 / 8 / 4 / 0 / **1** = 57 | 0.0 / 36.4 / 0.0 | 59.9 / 60.3 / 67.8 |
| genus:Contig_156_34.646_assay0 | 91 | **100%** | 0 | 36 / 8 / 4 / 0 / 0 = 48 | 0.0 / 0.0 / 39.9 | 59.9 / 60.0 / 67.4 |
| genus:Contig_429_56.0727_assay1 | 85 | **100%** | 0 | 39 / 8 / 4 / 0 / 0 = 51 | 36.6 / 42.6 / 37.9 | 60.0 / 60.1 / 65.8 |
| genus:Contig_80_59.5616_assay0 | 96 | **100%** | 0 | 34 / 8 / 4 / 0 / 0 = 46 | 38.2 / 46.4 / 36.3 | 60.0 / 60.0 / 67.9 |
| genus:Contig_590_36.7573_assay1 | 141 | 1,373/1,376 = 99.78% | 0 | **54** / 8 / 0 / **2** / **3** = **67** | 0.0 / 0.0 / 34.5 | 60.4 / 60.0 / 62.1 |
| genus:Contig_678_54.2745_assay0 | 108 | 1,375/1,376 = 99.93% | 0 | 52 / 8 / 0 / **2** / **3** = 65 | 45.3 / 42.7 / 0.0 | 60.0 / 60.0 / 68.1 |

Every genus candidate is specific **by absence**: `design` found the amplicon in 0 of the 315
non-classical genomes, as the published assay's is.

Read the held-out columns before calling any of this an improvement. The 103 held-out genomes are the ones
the paper's panel did not contain — unidentified `Brucella sp.` and the classical-clade species
*B. pinnipedialis*, *B. amazoniensis*, *B. vulpis* and *B. inopinata* — and for a **genus** assay, covering
more of them is better, not worse. No candidate covers all four named species, and neither does the
published assay:

- the five 100% candidates all **lose the two *B. vulpis* genomes** the published assay covers;
- the only two candidates that cover all three *B. inopinata* assemblies, which the published assay does
  not cover at all, give up *B. amazoniensis* and do not reach 100% on the panel species.
  `Contig_228` covers one of the three and keeps *B. amazoniensis*, which makes it the broadest of the
  five that do reach 100%.

So on the eight species of the panel a new assay is strictly better — the two genomes Bru_genus misses,
the ones the paper itself documents as its only in-silico failures, are covered at 100%. Across the genus
as a whole it is a lateral move. **`genus:Contig_193_65.5561_assay0` is the one worth looking at**: 100%
of the panel species, nothing in the exclusion group, and the only candidate with no predicted structure
at all on any of its three oligos.

```
Bru_genus (published)   F CGCGCATGATACCACTTCT    R GTGTTCGCGCATGTTGAAG    P AGGTGCGGTCGATAAACTTCTGGAAC
Contig_193_assay0 (new) F GATGTCGAGATAATGCACCC   R AAGGATTTGCACTGGTTCTC   P CCGTTACTGAAAATGGCCCGGC
```

## The species assays: slightly better coverage, the same specificity, a stronger discrimination basis

| assay | bp | own species | exclusion hits | which genomes | differences under one primer |
|---|---|---|---|---|---|
| **BA** (published) | 75 | 391/398 = 98.24% | 1 | GCF_055389475.1 | **one**, and 114 bp from nothing else |
| abortus:Contig_30_48.3778_assay3 | 87 | **392/398 = 98.49%** | 1 | **the same genome** | **two**, weight 1.0 |
| **BM_10** (published) | 147 | 728/734 = 99.18% | 4 | GCF_016091955.1, …65.1, GCF_018604785.1, GCF_041861985.1 | two, **114 bp apart** |
| melitensis:Contig_387_57.8_assay3 | 92 | **730/734 = 99.46%** | 4 | **the same four genomes** | **two, under one primer**, weight 2.0 |

Three things are worth drawing out of that table.

**The exclusion hits are identical, and in both cases they are the mislabelled assemblies.** The genome the
*B. abortus* assays "cross-react" with is GCF_055389475.1, labelled *melitensis* but carrying only the
*abortus* assay; the four for *melitensis* are the four labelled *abortus* but carrying only the
*melitensis* assay (`mislabelled.txt`). Neither the published nor the proposed assays have a demonstrated
cross-reaction on this panel — they have the same RefSeq label errors.

**The proposed assays rest on two differences under a single primer, which is what primer-finder is for.**
The published BA rests on one mismatch and BM_10 on two 114 bases apart, one under each primer — the very
configuration that stops `find` reporting their regions in the first place (`SUMMARY.md`). Whether two
mismatches under one primer really discriminate better than one under each of two is the judgement the
tool is built around; it is not settled by this record.

**The coverage gain is small and in the same direction for every candidate**: +1 genome for *B. abortus*,
+2 for *B. melitensis*, out of 398 and 734. That is a fraction of a percent and well inside what a
different panel would move.

## Where the published assays sit against primer-finder's own limits

`thermo.py` on the four published assays, under the design command's default conditions:

| assay | forward hairpin | probe hairpin | over the 47 °C default? |
|---|---|---|---|
| Bru_genus | 0.0 | 0.0 | no — and no dimer above 14.6 either |
| BA | 45.4 | 34.6 | no, but the forward primer is 1.6 °C under the limit |
| BM_10 | **76.4** | **61.7** | **yes, both** |
| BS | 0.0 | **63.9** | **yes, the probe** |

**Three of the four published assays exceed `--max-hairpin-tm` 47, so `design` would never have proposed
them** — and all four are wet-lab validated on 235 target and 192 non-target strains. The most useful
thing in this file may therefore be about the tool rather than the assays: the 47 °C default is
conservative enough to reject working assays, and this is the only wet-lab evidence in the whole
validation suite that bears on it. It is a reason to make the limit easy to relax, which it is
(`--max-hairpin-tm`), and a reason not to read a rejection as a verdict.

## The rating asked for

| axis | genus | *B. abortus* | *B. melitensis* | *B. suis* |
|---|---|---|---|---|
| in-silico coverage of the panel species | **new better** (100% vs 99.85%) | new better by 1 genome | new better by 2 genomes | no candidate |
| in-silico specificity | tied (0 vs 0) | tied, same label errors | tied, same label errors | — |
| breadth across the whole genus | **published at least as good** | n/a | n/a | — |
| predicted structure under one reaction | tied at best (both clean) | new better (no 45.4 °C hairpin) | **new better** (no 76.4 °C hairpin) | — |
| discrimination basis | both by absence | **new stronger** (2 under one primer vs 1) | **new stronger** (2 under one primer vs 2 at 114 bp) | — |
| evidence it works | **published only** | **published only** | **published only** | — |

**The answer: no, not on this evidence.** The candidates match or slightly beat the published assays on
every sequence measure, and `genus:Contig_193_65.5561_assay0` is a genuinely attractive design — but the
margins are one or two genomes out of hundreds, the genus candidates trade panel coverage for genus
breadth, and the published assays have the one thing none of the candidates has. A sequence-level
improvement of 0.15 percentage points does not justify replacing a validated assay.

Where this does point somewhere useful:

1. **If the two genomes Bru_genus misses matter** — GCF_002290125.1 and GCF_023651695.1, both with more
   than three mismatches at primer sites — then a second genus assay is available that covers them, and
   `genus:Contig_193_65.5561_assay0` is the candidate with the cleanest predicted structure.
2. **If full genus coverage matters**, neither the published assay nor any single candidate reaches it.
   all three *B. inopinata* assemblies are covered only by `Contig_590` and `Contig_678` (and one of the
   three by `Contig_228`); *B. vulpis* only by those two candidates and the published assay. Putting
   *B. pinnipedialis*, *B. vulpis*, *B. amazoniensis* and *B. inopinata* into the inclusion group and
   re-running would ask the question properly, and is the obvious next run.
   **Done: [`WHOLE_CLADE.md`](WHOLE_CLADE.md).** It changes this file's conclusion for the genus assay —
   with all twelve species in the group, one proposal covers every one of the 1,393 genomes with no
   predicted structure at all, and gives up nothing measured here.
3. **A *B. suis* assay needs the groups redrawn**, because *B. canis* is inside *B. suis*.
4. **The hairpin default rejects working assays.** Worth remembering before trusting `design` to have
   found everything.

## Files

| File | What it holds |
|---|---|
| `candidates.tsv` | the 30 assays compared: 4 published and 26 proposed, one per region, with the design columns |
| `scorecard.tsv` | each assay against all 1,794 genomes: coverage of its own group, exclusion hits, held out |
| `thermo.tsv` | Tm, GC, hairpin and dimer temperatures for all 30, under one set of conditions |
