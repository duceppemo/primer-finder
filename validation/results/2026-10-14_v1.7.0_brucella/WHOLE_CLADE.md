# A genus assay for the whole classical clade, not just the paper's eight species

2026-10-14, primer-finder 1.7.0. Follows [`ALTERNATIVES.md`](ALTERNATIVES.md), which found that every
candidate reaching 100% in-silico coverage of the paper's eight species gave up species outside them, and
that neither the published assay nor any candidate spanned the whole classical clade. So the four
remaining classical-clade species were put into the inclusion group and the question asked again.

**The answer is yes, and without a trade-off this time.** One proposed assay carries all three oligos
exactly in **every one of the 1,393 genomes of all twelve species**, in none of the 315 non-classical
relatives, and has **no predicted secondary structure of any kind**. It is still in silico: it has no
wet-lab data, and the published assay has it on 235 target and 192 non-target strains.

## The group

`panel.py` now labels *B. pinnipedialis*, *B. amazoniensis*, *B. vulpis* and *B. inopinata*
`classical_extra` instead of holding them out, and `group_genomes.py` has a **`genus_wide`** target whose
inclusion group is all twelve classical-clade species. The eight-species `genus` target is untouched, so
the record of the published assays stays reproducible.

| | genomes |
|---|---|
| the paper's eight species | 1,376 |
| *B. pinnipedialis* 8, *B. amazoniensis* 4, *B. inopinata* 3, *B. vulpis* 2 | **17** |
| `genus_wide` inclusion group | **1,393** |
| non-classical exclusion group | 315 |
| still held out: unidentified `Brucella sp.` | 86 |

Unidentified `Brucella sp.` stays out of both groups whatever the target: roughly half of those assemblies
are former *Ochrobactrum*, so they would corrupt either group.

## What the wider group costs, and what the published assay scores on it

| `-p` | required | regions | published assay |
|---|---|---|---|
| **1.0** (default) | 1,393 of 1,393 | **56** | not among them |
| 0.996 | 1,388 | 2,103 | not among them |
| 0.99 | 1,380 | 5,137 | recovered, rank 3,954 |

Fifty-six regions are present in every genome of all twelve species and absent from every relative, where
the eight-species group gives 398 — **adding 17 genomes of four species removes 86% of the candidates.**
Most of that is *B. inopinata*, the most divergent of the twelve.

The published Bru_genus scores **1,388 of 1,393 = 99.64%** on this group against 99.85% on the eight
species: the three *B. inopinata* assemblies are added to the two genomes it already missed. In all three
the amplicon is there and the forward primer site is intact; what differs is one base under the probe and,
in two of them, the last base of the reverse primer (`SUMMARY.md`).

## The assay that spans the clade

`design` proposed 334 assays on those 56 regions. Scoring the structurally cleanest twenty of them, plus
the published assay, on one yardstick (`wide_clean_scorecard.tsv`, `wide_clean_thermo.tsv`):

| assay | bp | all twelve species | relatives | hairpins F/R/P | dimers F/R/P, pair | Tm F/R/P |
|---|---|---|---|---|---|---|
| **Bru_genus** (published) | 93 | 1,388/1,393 = 99.64% | 0/315 | 0.0 / 0.0 / 0.0 | 14.6 / 14.6 / 0.0, 14.6 | 62.0 / 62.2 / 68.3 |
| **wide:Contig_83_60.1233_assay19** | 109 | **1,393/1,393 = 100%** | **0/315** | **0.0 / 0.0 / 0.0** | **0.0 / 0.0 / 0.0, 0.0** | 59.8 / 60.5 / 65.7 |
| wide:Contig_83_60.1233_assay20 | 145 | **100%** | 0/315 | 0.0 / 0.0 / 0.0 | 0.0 / 0.0 / 0.0, 0.0 | 59.8 / 59.0 / 65.7 |
| wide:Contig_83_60.1233_assay17 | — | **100%** | 0/315 | 0.0 / 0.0 / 0.0 | 9.4 / 0.0 / 0.0, 0.0 | 60.4 / 59.7 / 68.0 |
| wide:Contig_75_18.6543_assay0 | 79 | **100%** | 0/315 | 0.0 / 44.5 / 45.8 | 0.0 / 36.9 / 23.9, 0.0 | 60.2 / 60.3 / 68.0 |
| wide:Contig_59_21.6165_assay5 | 89 | 1,388 = 99.64% | 0/315 | 0.0 / 0.0 / 0.0 | all 0.0 | — |
| 15 others scored | 71-145 | 1,387-1,390 | 0/315 | — | — | — |

Four assays reach 100%, three of them on one region, `Contig_83_60.1233`. The pick of them:

```
Bru_genus (published)      F CGCGCATGATACCACTTCT   R GTGTTCGCGCATGTTGAAG   P AGGTGCGGTCGATAAACTTCTGGAAC    93 bp
Contig_83_assay19 (new)    F ATTTCATAAATCCCCTGCCG  R AATGCCTTTTGTTCCCTTCC  P TGTGTGGTAGGTTGACGCATCA       109 bp
```

Species by species, `Contig_83_60.1233_assay19` holds all three oligos exactly, on a single contig:

| | | | |
|---|---|---|---|
| *B. abortus* 398/398 | *B. melitensis* 734/734 | *B. suis* 120/120 | *B. canis* 74/74 |
| *B. ceti* 21/21 | *B. ovis* 18/18 | *B. neotomae* 6/6 | *B. microti* 5/5 |
| *B. pinnipedialis* 8/8 | *B. amazoniensis* 4/4 | *B. vulpis* 2/2 | ***B. inopinata* 3/3** |
| non-classical relatives **0/315** | unidentified `Brucella sp.` 57/86 | | |

It also covers GCF_002290125.1 and GCF_023651695.1, the two assemblies the paper reports as its only
in-silico misses, and the 57 of 86 unidentified `Brucella sp.` it hits are more than the published assay's
49 — consistent with it reaching further across the clade rather than less far.

## How this differs from the earlier comparison

`ALTERNATIVES.md` concluded that nothing proposed was worth acting on: the gains were a fraction of a
percent and every 100% candidate traded breadth for them. Widening the group changes that, because the
trade-off was an artefact of asking the eight-species question:

| | eight species | twelve species |
|---|---|---|
| best candidate's coverage | 100% of 1,376, but 0 of 2 *B. vulpis*, 0 of 3 *B. inopinata* | **100% of 1,393, every species** |
| published assay | 99.85% of 1,376, 0 of 3 *B. inopinata* | 99.64% of 1,393 |
| cost of the best candidate | loses species the published assay covers | **loses nothing measured here** |

So the recommendation changes. On sequence, `Contig_83_60.1233_assay19` is better than the published assay
on every axis measured: it covers five genomes the published assay does not, including all three
*B. inopinata*, nothing in the exclusion group, and no predicted structure at all where the published
assay has 14.6 °C dimers. The two caveats that matter:

- **no wet-lab data.** The published assay has 235 target and 192 non-target strains behind it. Nothing
  here substitutes for that, and a design that looks better in silico can fail in a reaction;
- **its primers are AT-rich at the 5' end** (`ATTTCATAAATCCCCTGCCG`, 45% GC, with a four-base C run) where
  the published primers are 52.6% GC and better balanced. Primer3 accepted it under primer-finder's
  defaults and predicts no structure, but predicted structure is not the only thing that makes a primer
  behave, and this is the kind of difference a wet-lab test would settle.

If a genus assay covering *B. inopinata* is wanted, this is where to start, and the sensible step is to
test it beside the published one rather than instead of it.

## Reproducing

```bash
bash validation/brucella/run.sh /path/to/work 48 128     # now runs genus_wide as well
```

## Files

| File | What it holds |
|---|---|
| `wide_runs.tsv` | the three `genus_wide` runs, with region counts and what became of the published assay |
| `wide_scorecard.tsv` | the first shortlist — one assay per region — against all 1,794 genomes |
| `wide_clean.tsv`, `wide_clean_scorecard.tsv`, `wide_clean_thermo.tsv` | the twenty structurally cleanest proposals, scored and measured |
| `wide_candidates.tsv` | the shortlist's sequences and the design columns behind them |
