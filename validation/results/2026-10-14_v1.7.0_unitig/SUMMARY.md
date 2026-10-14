# What the unitig assembler and the last-resort pass change, on a panel with an answer key

2026-10-14, primer-finder 1.7.0 plus the two unreleased features. Both were prompted by the *Brucella*
records, where SKESA refused the kmers of a real locus and the find rule turned down three published
assays. This measures them where there is a published answer to check against: the *Listeria
monocytogenes* clonal-complex panel of
[`2026-10-10_v1.5.0`](../2026-10-10_v1.5.0/SUMMARY.md), rebuilt from scratch — 784 complete RefSeq
genomes, 763 of which the Institut Pasteur scheme places in a clonal complex, identical to the original.

**The question this answers** is whether either feature is safe to make the default. The answer is that
the last-resort pass is inert here and the unitig walk is an improvement that costs one rank.

## The baseline reproduces, which is what makes the rest worth reading

The nine targets with SKESA at the defaults give **332 regions** — CC2 1, CC3 40, CC4 6, CC5 57, CC6 53,
CC7 60, CC8 62, CC9 19, CC224 34 — and recover four published assays at ranks 10, 5, 1 and 2. Those are
the figures of [`2026-10-14_v1.6.0`](../2026-10-14_v1.6.0/SUMMARY.md) to the region, on a panel downloaded
again and typed again. Nothing in the two features has moved the existing behaviour.

## `-a unitig`: 40% more regions, the same four assays, one rank worse

| target | SKESA contigs → regions | unitig contigs → regions | published assay |
|---|---|---|---|
| CC2 | 32 → 1 | 50 → **4** | not recovered either way |
| CC3 | 485 → 40 | 681 → **59** | not recovered either way |
| CC4 | 99 → 6 | 142 → **8** | not recovered either way |
| CC5 | 468 → 57 | 652 → **84** | recovered, **rank 10 → 15** |
| CC6 | 354 → 53 | 490 → **67** | recovered, rank 5 → 5 |
| CC7 | 533 → 60 | 741 → **85** | recovered, rank 1 → 1 |
| CC8 | 414 → 62 | 584 → **83** | recovered, rank 2 → 2 |
| CC9 | 158 → 19 | 243 → **29** | not recovered either way |
| CC224 | 288 → 34 | 387 → **47** | not recovered either way |
| **total** | **2,831 → 332** | **3,970 → 466** | **4 of 9 either way** |

**+40% regions, every target, and not one verdict changed.** The four targets whose published assay is
recoverable still recover it; the five the 2026-10-10 record explains as not recoverable on this panel
still do not.

The one cost is CC5's rank, 10 of 57 to 15 of 84. That is what more candidates does to a rank rather than
a worse answer — the region itself is the same one, 752 bases with SKESA and 804 with the unitig walk. The
other three ranks are unchanged.

Two things the region lengths show about how the walk differs:

- it usually recovers **more** of a stretch (CC5 752 → 804, CC6 1,449 → 1,501, CC8 2,208 → 2,260), because
  SKESA trims ends it is unsure of;
- sometimes **less** (CC7 8,142 → 6,591), because a unitig ends at a branch where SKESA resolves one. That
  is the deliberate difference: a unitig is what the kmers determine, and SKESA will pick between
  alternatives. Neither is wrong; they answer different questions.

## The last-resort pass: inert on this panel

**It did not fire once in eighteen runs.** The `last_resort` column of `compare.tsv` is `no` throughout.
That is the result worth having: on a panel where the strict rule works, a run is never relaxed and no
region carries the mark. The pass can only change a run that would otherwise have reported nothing.

Where it does fire — the *Brucella* `suis_lato` target, which has exactly one candidate locus — the two
features together carry it through a pipeline that could not report it before:

```
>unitig_1 [178] inclusion=193/193 last_resort        215 bases
```

`-a unitig` builds that stretch from the 15 kmers SKESA refused, and the last-resort pass keeps it on its
single difference — at position 178, the same base measured by hand against *B. abortus* in
[`SUIS_LATO.md`](../2026-10-14_v1.7.0_brucella/SUIS_LATO.md). It is still a lead and not a result, which is
what the mark says, and it still does not hold the published BS assay (that assay is *B. suis* sensu
stricto and covers 1 of 74 *B. canis*).

## Should either become the default?

**The last-resort pass already is on by default**, and this record is the argument for that being safe: it
cannot affect a run that finds anything, and everything it does report is marked in four places.
`--no-last-resort` restores the old error.

**`-a unitig` is not the default and this record does not settle whether it should be.** In its favour:
40% more candidates, no verdict changed, no external program, and it honours `--min-region` down to a
single kmer where SPAdes cannot go below 200 bases. Against: SKESA has years of use behind its output and
this has one panel; it reports a shorter region where a branch is resolvable; and a rank moved. The
cautious reading is that it is the right choice when a run ends at "could not assemble", when the specific
stretches are short, or when the kmer count is small — which is what the documentation says — and that
making it the default wants a second panel.

## A defect this record found before it was published

The walk as first written lost kmers where a kmer's last `k-1` bases were their own reverse complement: a
hairpin led the back-up into the same kmer the other way round, and the seed ended in no stretch at all.
It was found by a review, not by this panel — at `-k 99` the trigger needs a palindrome of 98 bases or
more, which these genomes do not have. All nine targets were re-run after the fix and report the same
numbers to the region (`recheck.tsv`), which is why the table above is still good.

## Reproducing

```bash
conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
bash validation/listeria/run.sh /path/to/work 48 128     # the baseline, with SKESA
# then, per target, the same groups with -a unitig
```

## Files

| File | What it holds |
|---|---|
| `compare.tsv` | the eighteen runs: contigs, regions, recovery, rank, region length, whether `last_resort` fired |
| `recheck.tsv` | the nine targets re-run after the hairpin fix, against the numbers above |
| `panel.tsv.gz` | the 784 genomes and their clonal complexes, as this record used them |
