# 2026-10-08 — primer-finder 1.4.0, spending a region's differences on one oligo

`primer-finder find` keeps a region when its differences from an exclusion genome **could sit in one
oligo**: a run of them, or two of them fewer than 21 matching bases apart
([Methods](../../../docs/wiki/Methods.md)). That is the only reason a region is a candidate. This record
asks whether the design step was actually placing an oligo over them, and the answer was: often not.

Same genome set and same regions as [`../2026-10-08_v1.3.0/SUMMARY.md`](../2026-10-08_v1.3.0/SUMMARY.md) —
25 *Xylella fastidiosa* genomes, subsp. *multiplex* as the inclusion group (8 genomes) against the other
four subspecies (17), 711 regions, the first 20 designed on. primer-finder 1.4.0, Primer3 2.6.1, blast
2.17.0+, insilicoPCR 0.6.1.

## 1. What the regions actually offer

For each of the 20 regions designed on, the marked (lower-case) differences, asked two ways: is there a run
of two or more consecutive ones, and is there a pair that is **not** consecutive but still closer than a
primer's length?

| | regions |
|---|---|
| hold a close pair that is **not** contiguous | **19** of 20 |
| hold a run of 2 or more consecutive differences | 12 of 20 |
| hold close pairs **only** — no run at all to aim at | **8** of 20 |

Up to 1.3.0 the forced Primer3 requests aimed at runs of consecutive differences only. For 8 of these 20
regions there was no such run, so every forced request pinned a primer's 3' end on a **single** base and the
pair the region was kept for went unused.

## 2. Aiming at a window instead

A target is now the stretch of at most 25 bases — the widest primer Primer3 may return — holding the most
differences; the four richest non-overlapping ones get requests. A run of three is such a window, so nothing
is lost; a pair seven bases apart is one too, which is what is gained.

Same regions, same options, only the target changed:

| | 1.3.0 (runs) | 1.4.0 (windows) |
|---|---|---|
| assays Primer3 proposed | 472 | **511** |
| of them resting on a difference | 287 | 309 |
| best single oligo covers 1 difference | 163 | **68** |
| …2 | 62 | 109 |
| …3 | 28 | 68 |
| …4 or more | 34 | 64 |
| **2 or more under one oligo** | 124 (43%) | **241 (77%)** |

That is the whole point of the change: the number of assays that use the evidence their region was selected
for, rather than one base of it.

## 3. Setting aside the rest

`--min-oligo-differences`, 2 by default, now carries no further an assay that rests on differences and
covers fewer than two under any one oligo. Assays specific by absence are exempt — they do not rest on a
difference at all.

```bash
primer-finder design multiplex/results -o multiplex/assays -t 32 --max-regions 20 \
    --insilico-pcr /path/to/insilicoPCR-linux-x64
```

| | 1.3.0 | 1.4.0 |
|---|---|---|
| proposed | 472 | 511 |
| carried forward | 465 | **436** (195 by absence, 241 by difference) |
| set aside for one difference under one oligo | — | 68 |
| would amplify both groups | 7 | 7 |
| in silico PCR, qPCR mode: selective | 311 / 465 (67%) | **383 / 436 (88%)** |
| in silico PCR, qPCR mode: `undecided (3' end)` | 53 | 3 |
| in silico PCR, PCR mode: selective | 220 / 465 (47%) | **305 / 436 (70%)** |
| in silico PCR, PCR mode: `undecided (3' end)` | 96 | 17 |

2 min 24 s for 20 regions on 32 threads, including all four in silico PCR runs.

**Do not read the selectivity rise as better assays.** in silico PCR decides a primer match by counting
mismatches, so demanding more differences under one oligo raises its pass rate close to by construction —
the same circularity as everywhere else in these records. Two things in that table are *not* circular:

- the **`undecided (3' end)` count collapses** (96 → 17 in PCR mode). Those are the assays whose only
  exclusion amplifications rest on bases insilicoPCR cannot see, which is what a single difference at a 3'
  end gives you. An assay with two differences under one oligo almost always has one of them more than two
  bases in, where the check can judge it. So the share of the assays it was given that it can say anything
  about at all — cleared or refused, rather than `undecided` — rose from 369 of 465 (79%) to 419 of 436
  (96%);
- the 68 assays set aside are not borderline — they rest on **one** difference, which is the case the
  `find` rule was written to exclude and which the literature on single mismatches does not support
  ([Lefever et al. 2019](https://doi.org/10.1038/s41598-019-38581-z)).

## What this does and does not show

- It shows the design step now uses the property the regions were selected for, measured on the same
  regions with everything else held fixed: 43% → 77% of difference-based assays carry two or more
  differences on one oligo.
- It shows the filter costs candidates rather than regions: all 20 regions still yield assays, 436 of them.
- It does **not** show that any of these assays works, and the in silico pass rates are not evidence that
  they do. No assay here has been near a bench.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/assays \
    -t 32 --max-regions 20 --insilico-pcr /path/to/insilicoPCR-linux-x64
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/loose \
    -t 32 --max-regions 20 --min-oligo-differences 1    # what 1.3.0 carried forward
```

The region shapes of section 1 come from `difference_targets` in `primer_finder/design.py`, applied to
`multiplex/results/final_kmers.fasta`.
