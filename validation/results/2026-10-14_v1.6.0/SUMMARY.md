# 2026-10-14 — two defaults put to the *Listeria* panel: the kmer size, and the two-differences rule

Both of these were chosen by argument rather than by measurement. This record measures them on the nine
clonal-complex targets of [`../2026-10-10_v1.5.0/SUMMARY.md`](../2026-10-10_v1.5.0/SUMMARY.md) — the same
784-genome panel, 18 to 86 inclusion genomes against 677 to 745 — and changes one documented claim that
turns out to be false. Raw results in [`sweeps.tsv`](sweeps.tsv).

primer-finder 1.6.0. Recovery is judged as in the earlier record: all three published oligos in one reported
region, primers facing, probe between, amplicon the published length.

## The kmer size, `-k 99`

| `-k` | regions over the nine runs | runs with no region at all | published assays recovered |
|---|---|---|---|
| 31 | 70 | 2 (CC2, CC4) | 4 of 9 |
| 51 | 160 | 1 (CC2) | 4 of 9 |
| **99** (shipped) | **306** | **0** | **4 of 9** |
| 127 | 319 | 0 | 4 of 9 |

**A shorter kmer finds fewer regions, not more.** A kmer survives the subtraction only if no exclusion
genome holds it, so for a region that does exist in the exclusion genomes the kmer has to span a difference.
A long kmer spans one more often, so more kmers survive and the stretches they assemble into are longer; a
short kmer survives only when it sits almost on top of a difference. At `-k 31` two of the nine runs produce
no candidate region at all.

The wiki said the opposite — "a shorter `-k` finds shorter specific regions **and more of them**" — and that
sentence has been corrected to the table above. It was a plausible guess that nobody had checked.

**Which assays are recovered does not move**: the same four (CC5, CC6, CC7, CC8) at every kmer size, at much
the same ranks. The kmer size buys candidates, not answers. 99 sits near the top of the yield with no empty
runs; 127 adds 4% more candidates and nothing else, at more exposure to assembly error, since one wrong base
spoils every kmer covering it.

## The two-differences rule, 21 bases

`find` keeps a contig when two of its differences are fewer than `mapping.PRIMER_LENGTH` = 21 matching bases
apart, on the grounds that both could then sit in one oligo. But the design step's widest primer is
`MAX_PRIMER_SIZE` = 25, and `difference_targets` aims at windows of up to 25 bases. The two numbers disagree,
and the find step is the stricter of the two:

| matching bases between the pair | the pair spans | `find` keeps it | a 25-base primer covers it |
|---|---|---|---|
| 20 | 22 | yes | yes |
| 21 | 23 | **no** | yes |
| 22 | 24 | **no** | yes |
| 23 | 25 | **no** | yes |
| 24 | 26 | no | no |

So pairs spanning 23, 24 or 25 bases are thrown away by `find` although a primer could carry both. The value
that expresses "a 25-base primer covers both" is **24**, not 25: the rule counts the *matching* bases
between the two differences, so 24 admits spans up to 25 bases and 25 would admit a span of 26, which no
primer can cover.

Run at 24, everything else unchanged:

| | rule 21 (shipped) | rule 24 |
|---|---|---|
| regions over the nine runs | 306 | **322** (+5.2%) |
| published assays recovered | 4 of 9 | 4 of 9 |
| ranks of those four | 10, 5, 1, 2 | 10, 5, 1, 2 |

Per target the change is +0 to +4 regions; no target loses any, and no rank moves.

## What shipped, and why its numbers differ from the experiment above

The experiment was crude: it set `mapping.PRIMER_LENGTH` to 24 and re-ran. That one constant did **three**
jobs, so the experiment moved all three at once — the pair rule (the thing being tested), the clipped-end
rule (from "at least 21 unplaceable bases" to 24, *stricter*), and the window `blast.has_close_variants`
uses in step 5 (to 24). The +5.2% was the net of a loosening and two tightenings.

What shipped separates them. `MAX_PRIMER_SPAN = 25` is the pair rule and the step-5 window, stated as a
span — two differences qualify when the pair, counting both differences and the matches between them, is at
most 25 bases wide, which is the widest primer Primer3 may return. `MIN_CLIPPED = 21` keeps the value the
clipped-end rule has always had, because how much unplaceable sequence is worth keeping is a different
quantity from how far apart two differences may be, and raising it would only throw regions away.

Re-run on the same nine targets with that:

| target | span 22 (shipped before) | span 25 (now) | change |
|---|---|---|---|
| CC2 | 1, not recovered | 1, not recovered | +0 |
| CC3 | 38, not recovered | 40, not recovered | +2 |
| CC4 | 5, not recovered | 6, not recovered | +1 |
| CC5 | 49, **rank 10** | 57, **rank 10** | +8 |
| CC6 | 46, **rank 5** | 53, **rank 5** | +7 |
| CC7 | 57, **rank 1** | 60, **rank 1** | +3 |
| CC8 | 59, **rank 2** | 62, **rank 2** | +3 |
| CC9 | 19, not recovered | 19, not recovered | +0 |
| CC224 | 32, not recovered | 34, not recovered | +2 |
| **total** | **306 regions** | **332 regions** | **+26 (+8.5%)** |

Four of nine recovered either way, at identical ranks, and CC4 at `-p 0.94` is still recovered at rank 1 (10
regions instead of 9) — so five of nine counting that one, as in the earlier records. CC121 at `-p 0.93`
gives 71 regions instead of 66 and is still not recovered, which is the island-edge case the 1.5.0 record
describes and not something a wider span can fix.

So the shipped change gains 8.5% rather than the 5.2% the experiment suggested, because the experiment was
paying for two tightenings it did not intend. No target loses a region and no rank moves.

## What this does and does not show

- It shows the shipped kmer size is a good choice on this panel, and that the documented reason for
  preferring a shorter one was wrong.
- It shows the find step's 21 was inconsistent with the design step's 25; expressing the rule as a 25-base
  span instead gains 8.5% more candidate regions on this panel at no cost to recovery, and the clipped-end
  rule is a separate quantity that stays where it was.
- It does **not** show that the extra regions are good ones: nothing here runs `design` on them, and a
  region kept because two differences sit 25 bases apart gives Primer3 a harder placement problem than one
  where they touch.
- Recovery is 4 of 9 throughout because the other five are explained in the earlier record: four of the
  published assays are not specific on this panel, and CC4 needs `-p 0.94`. None of that is affected by
  either default.
- Nothing here has been near a bench.

## Reproducing

```bash
bash validation/listeria/run.sh /path/to/work 32 64          # the panel and the shipped settings
for k in 31 51 127; do
  primer-finder -i <work>/CC7/inclusion -e <work>/CC7/exclusion -o <work>/k$k -t 32 -m 64 -k $k
done
# the rule was first measured by setting the old mapping.PRIMER_LENGTH to 24; what shipped is
# MAX_PRIMER_SPAN = 25 with MIN_CLIPPED left at 21, which is the "span 25" column above
```
