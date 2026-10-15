# What the unitig assembler and the last-resort pass change, on three panels with answer keys

2026-10-14, primer-finder 1.7.0 plus the two unreleased features. Both were prompted by the *Brucella*
records, where SKESA refused the kmers of a real locus and the find rule turned down three published
assays. This measures them where there is a published answer to check against: the *Listeria
monocytogenes* clonal-complex panel of
[`2026-10-10_v1.5.0`](../2026-10-10_v1.5.0/SUMMARY.md), rebuilt from scratch — 784 complete RefSeq
genomes, 763 of which the Institut Pasteur scheme places in a clonal complex, identical to the original.

**The question this answers** is whether either feature is safe to make the default. It was first written
from the *Listeria* panel alone; two more were added afterwards, and the three together are what the
recommendation at the end rests on: **19 targets, 51 runs, 21 of them paired on the same groups with both
assemblers.**

| panel | genomes | targets | regions, SKESA → unitig | published assays recovered |
|---|---|---|---|---|
| *Listeria* clonal complexes | 784 complete | 9 | 332 → **466** (+40%) | 4 of 9, **the same 4 either way** |
| *Xylella* subspecies | 25 complete | 4 | 3,985 → **5,990** (+50%) | **4 of 4 either way** |
| *Brucella* | 1,794, draft-heavy | 6 | 5,811 → **13,873** (+139%) | 1 of 1 testable, either way |

**Across all three, the unitig walk never lost a recovery, and in one case it gained a run outright.**
`last_resort` fired once in 51 runs, on the one target that has a single candidate locus.

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

**It did not fire once in eighteen runs.** The `last_resort` column of `listeria_compare.tsv` is `no` throughout.
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

## The two panels added afterwards

### *Xylella fastidiosa*: four published assays, recovered either way

25 complete RefSeq genomes, each subspecies against the others
([`xylella_compare.tsv`](xylella_compare.tsv)). This is the easiest panel in the suite — the groups are
distinct enough that thousands of regions come through — so it tests the opposite end from *Brucella*.

| target | SKESA contigs → regions | unitig contigs → regions | rank |
|---|---|---|---|
| multiplex | 3,494 → 923 | 7,132 → **1,223** | 1 → 1 |
| fastidiosa | 3,205 → 612 | 5,506 → **782** | 5 → **2** |
| pauca | 3,150 → 2,040 | 7,246 → **3,471** | 173 → 263 |
| morus | 2,956 → 410 | 4,974 → **514** | 5 → **3** |
| **total** | **3,985** | **5,990** | all four recovered either way |

Two ranks improve, one is unchanged, one is worse. The rank moves in both directions, which is what a
larger candidate pool does to a position in it — not a change in whether the assay is found.

### *Brucella*: the panel that prompted the feature

The six targets of the *Brucella* records, at the settings those records used
([`brucella_compare.tsv`](brucella_compare.tsv)):

| target | settings | SKESA → unitig regions | published assay |
|---|---|---|---|
| genus | `-p 0.99` | 5,344 → **11,156** | recovered either way, rank 4,249 → 8,497 |
| genus | default | 398 → **2,159** (5.4×) | not recoverable at this `-p` either way |
| genus_wide | default | 56 → **537** (9.6×) | — |
| abortus | `-p 0.994` | 2 → 2 | not recovered either way |
| melitensis | `-p 0.994` | 11 → **18** | not recovered either way |
| suis_lato | `-k 201 -p 0.95` | **0, "could not assemble"** → **1** | — |

The gains are far larger here than on the other two panels, and the reason is the panel itself: 873 of
the 1,376 classical genomes are contig-level drafts, so the specific stretches are short and the kmer
counts per stretch small, which is exactly where an assembler built for reads discards the most. The
whole-clade `genus_wide` pool going from 56 regions to 537 is the most consequential single number in this
file: that is the set a genus assay spanning all twelve classical species has to be designed from.

The `suis_lato` row is the only run on any panel where SKESA fails outright, and the only one where
`last_resort` fired.

## Where the unitig walk outperforms SKESA, and by how much

[`bases_per_kmer.txt`](bases_per_kmer.txt) holds every paired run: the kmer size, how many
inclusion-specific kmers the kmer step found, how many bases each assembler turned them into, and SKESA's
bases per kmer. Twenty-one pairs. Two things come out of it.

**SKESA either works or collapses, with nothing in between.** Measured as bases assembled per kmer given:

| | SKESA bases per kmer | share of what the unitig walk recovers |
|---|---|---|
| `-k 201`, 15 kmers | **0.000** | 0% |
| `-k 31`, 31,256 kmers | **0.000** | 0% |
| `-k 51`, 23,413 kmers | **0.019** | 0.8% |
| the other **eighteen** runs, all `-k 99` | **0.702 – 1.388** | 10.7% – 69.8% |

A 37-fold gap separates the three collapses from the eighteen working runs, and every collapse is at a
kmer size away from the default. That is a quantity a run can measure on itself, without a second
assembly, and no threshold between 0.05 and 0.5 could pick wrongly on this evidence.

**But "SKESA is underperforming" is not measurable the same way.** Its share of the recoverable bases
ranges from **10.7% to 69.8%** across the eighteen healthy runs — the worst of them, *Brucella*
`genus_wide`, assembles 13 kb where the walk assembles 121 kb, nine times less, while its bases per kmer
(0.702) sits inside the healthy band. So a run cannot tell from its own numbers that it is leaving most of
the signal behind. Any rule that tried would either fire on healthy runs or miss this one.

## Should either become the default?

**The last-resort pass already is on by default**, and this record is the argument for that being safe: it
cannot affect a run that finds anything, and everything it does report is marked in four places.
`--no-last-resort` restores the old error.

**`-a unitig` should stay a choice, and a narrow automatic fallback should be added.** That is what the
three panels support, and it is less than making it the default.

In its favour: 40% to 139% more candidates on every panel, never a recovery lost over 19 targets, one run
rescued from an error, no external program, and `--min-region` honoured down to a single kmer where SPAdes
cannot go below 200 bases. Against: SKESA has years of use behind its output and this has three panels of
one day's work; a unitig ends at a branch SKESA would resolve, so a region can come back shorter (*Listeria*
CC7, 8,142 bases to 6,591); and a rank moves in both directions as the pool grows, which costs nothing if
the assay is still found but makes `--max-regions 50` cover a smaller share of a bigger pool.

**The fallback the evidence does support** is to use the unitig walk when the chosen assembler assembles
essentially nothing — below roughly 0.1 bases per kmer, against a healthy 0.7 to 1.4. It would have fired
on 3 of these 21 runs, every one of them a total failure, and on none of the eighteen that worked. Today
those three runs end in an error that tells the user to lower `-k`, which runs into the same assembler from
the other end; with the fallback they would report regions and say which assembler produced them.

**The fallback the evidence does not support** is anything cleverer than that. A run that is quietly
leaving 90% of its signal behind looks, from the inside, exactly like one that is not.

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
| `listeria_compare.tsv` | the eighteen *Listeria* runs: contigs, regions, recovery, rank, region length, whether `last_resort` fired |
| `xylella_compare.tsv` | the eight *Xylella* runs, the same columns |
| `brucella_compare.tsv` | the twelve *Brucella* runs, the same columns |
| `bases_per_kmer.txt` | all 21 paired runs: kmers, bases assembled by each, SKESA's share and its bases per kmer |
| `recheck.tsv` | the nine *Listeria* targets re-run after the hairpin fix, against the numbers above |
| `panel.tsv.gz` | the 784 *Listeria* genomes and their clonal complexes, as this record used them |
