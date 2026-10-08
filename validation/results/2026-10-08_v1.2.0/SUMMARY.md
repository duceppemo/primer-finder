# 2026-10-08 — primer-finder 1.2.0, designing assays on the *Xylella* regions

> **Note, 2026-10-08.** The in silico PCR numbers in this record were taken with insilicoPCR's mismatch
> tolerance at 0, which was primer-finder's default at the time. 1.3.0 measured what that tolerance governs
> and made 1 the default, so these counts are the optimistic end of the range. Nothing here is withdrawn —
> the run happened as described — but read it beside
> [`../2026-10-08_v1.3.0/SUMMARY.md`](../2026-10-08_v1.3.0/SUMMARY.md), which also explains why a difference
> in the last two bases of a primer cannot be seen by this check at all.

Two things were run on the genome set of
[`../2026-10-08_v1.1.0/SUMMARY.md`](../2026-10-08_v1.1.0/SUMMARY.md): the four published assays were
recovered again, unchanged, and the new `design` command was run on the subsp. *multiplex* regions.

**Read the second part as a dry run of the machinery, not as evidence about PCR.** Nothing here was tested
in a laboratory, and in silico PCR judges a primer by counting mismatches against its own tolerance, so its
verdicts partly restate its own rule. What this record shows is that the command works end to end on real
data and what it produced on one dataset.

## Designing on the subsp. *multiplex* regions

```bash
primer-finder design multiplex/results -o multiplex/assays -t 16 --max-regions 20 \
    --insilico-pcr /path/to/insilicoPCR-linux-x64
```

primer-finder 1.2.0, Primer3 2.6.1 (libprimer3 release 2.6.1), blast 2.17.0+, insilicoPCR 0.6.1.
3 min 49 s for 20 regions, including both in silico PCR runs against all 25 genomes.

| | |
|---|---|
| Regions designed on | 20 (of 711, the most promising first) |
| Assays proposed by Primer3 | 432 |
| Primer3 requests dropped on time | 9 of 260, all of them probe-pinning ones; their regions still produced assays |
| Could tell the groups apart | 417 — 137 because no exclusion genome holds the amplicon, 280 because an oligo sits on differing bases |
| Would amplify both groups | 15, listed but carried no further |
| In silico PCR, qPCR mode | 374 of 417 amplified every inclusion genome and no exclusion genome |
| In silico PCR, standard PCR mode | 278 of 417 |

### What the two modes show

The gap between the modes is the point of running both. Of the 280 assays that rest on differences, 47 have
their differences under the **probe only**:

| | qPCR mode (with probe) | PCR mode (primers only) |
|---|---|---|
| differences under the primers: 0 (probe only) | 47 / 47 | 1 / 47 |
| 1 | 86 / 117 | 52 / 117 |
| 2 | 50 / 60 | 44 / 60 |
| 3 | 26 / 26 | 24 / 26 |
| 4 or more | 30 / 30 | 30 / 30 |

Assays whose only differences are under the probe are reported selective in qPCR mode and not in PCR mode,
which is what the two files are for: `pcr_selective` says whether the primer pair discriminates on its own.

The rise with the number of differences is consistent with how in silico PCR decides a primer match
(`--mismatches 0` by default), and is **not** independent evidence that more mismatches make a better assay.
The ranking's reasons are set out in `docs/wiki/Designing-assays.md`, with the published work they rest on.

### The region holding the published assay

18 assays were designed on `Contig_1563`, the region that contains the published XFM assay. All of them were
reported selective in both modes, and all are specific by absence — the exclusion genomes do not hold that
stretch at all, which is also why the published assay's own oligos sit on no differing base.

## What this does and does not show

- The command runs end to end on real data: Primer3, the amplicon checks against 25 genomes, both in silico
  PCR modes, and a ranked table, in under four minutes for 20 regions.
- A pathological Primer3 request can no longer stop it. Before the fixes in this version, one probe-pinning
  request ran for 45 minutes; 8 of 20 regions were then lost to the first, cruder timeout. Now 9 requests
  of 260 are dropped and every region still yields assays.
- It does **not** show that any of these assays works, or that the ranking predicts which will. No assay
  here has been near a bench.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32     # the four published assays
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/assays \
    -t 16 --max-regions 20 --insilico-pcr /path/to/insilicoPCR-linux-x64
```

insilicoPCR is a separate program: https://github.com/duceppemo/insilicoPCR
