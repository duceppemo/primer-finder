# Brucella, re-designed under 1.8.0: what the corrected amplicon rule changes

2026-10-16, primer-finder 1.8.0, on the regions the
[`2026-10-14_v1.7.0_brucella/`](../2026-10-14_v1.7.0_brucella/) runs found. `find` was not re-run: the
three changes in question all act at the design step, so the same `final_kmers.fasta` and the same
inclusion and exclusion folders were designed on again with the same settings
(`--product-size 70-150 --assays-per-region 3 --max-regions 50`).

What changed between the two versions and could move these numbers:

- whether an exclusion genome holds an amplicon is judged on the span a hit reaches and the identity among
  its aligned columns, not on blast's `length` and `pident`;
- a copy blast reports in several pieces is chained before that rule is applied;
- the amplicon searches use `-task blastn` rather than megablast.

## The totals

[`totals.tsv`](totals.tsv), per-assay in [`label_changes.tsv`](label_changes.tsv). Primer3 proposed the
same assays on the same regions in every case, so these are the same populations judged twice.

| target | assays | `absence` 1.7.0 → 1.8.0 | lost it | gained it |
|---|---|---|---|---|
| genus, as labelled, `-p 1.0` | 520 | 417 → 377 | 45 | 5 |
| genus_wide, as labelled, `-p 1.0` | 334 | 227 → 207 | 23 | 3 |
| melitensis, corrected, `-p 0.994` | 98 | 0 → 2 | 0 | 2 |
| abortus, corrected, `-p 1.0` | 6 | 0 → 0 | 0 | 0 |

**68 assays lost an unearned `absence`** and are now judged on how their primers differ, which is what the
corrected rule was for: `specific_by=absence` ranks above every difference-based assay and is exempt from
`--min-primer-differences`, so an assay that claimed it wrongly was both over-ranked and under-examined.
`best_primer_variants` did not change for a single assay in any run — only the question of what the
exclusion genomes hold.

**The conclusions of the earlier record stand.** The top of the genus ranking is unchanged, candidate for
candidate, including all five named in `ALTERNATIVES.md`; `abortus` is untouched; and the census figures
behind "100.00% against the published 99.85%" come from oligo presence rather than from this rule, so they
are unaffected.

## Ten assays moved the wrong way, and why

Ten assays **gained** `absence` — the unsafe direction, since it exempts them from the differences filter.
All ten have short amplicons (71 to 150 bases, median 98) and the cause is the blast task, not the rule.
Worked case, `Contig_120_57.4688_assay0`, a 72-base amplicon against *B. suis* GCF_000007505.1:

| | hit | coverage | verdict |
|---|---|---|---|
| 1.7.0, megablast | `q1-72`, length 73, one indel each side, 95.9% | 73 of 72 | holds it |
| 1.8.0, `-task blastn` | `q9-72`, span 64, ungapped, 98.4% | 64 of 72 = 88.9% | **absent** |

The genome holds the whole amplicon: megablast aligns all 72 bases through a one-base insertion and a
one-base deletion in the first eight. blastn will not pay its gap costs to extend an HSP across that for
the sake of four matching bases, so it starts the alignment eight bases in — and eight bases of 72 is 11%,
just past the 90% coverage rule. There is one HSP, so the chaining that exists for exactly this shape of
problem has nothing to join.

Eight of the ten concern one or two exclusion genomes each, where the label flips but little else does.
The other two are `Contig_120_57.4688_assay0` and `_assay1`, where **630 of 951** exclusion genomes went
from holding the amplicon to not holding it.

This is a defect in 1.8.0, found by this re-run and not yet fixed. It is the same failure the release set
out to correct, arriving from the other side: a genome that holds an amplicon is recorded as not holding
it, and the assay is handed a specificity claim it has not earned. Neither version is right in general —
1.7.0's own count rested on megablast, whose blind spot for a copy that differs every ~25 bases is what
`-task blastn` was adopted to close.

Fixing it means deciding what to do when a hit stops short of an amplicon's end by a few bases: the
missing flank can be compared against the subject directly, which is exact and needs the genome to hand,
or the coverage rule can be given an absolute allowance beside its fraction, which is simpler and
arbitrary. The first is the honest one. Until then, an `absence` verdict on an amplicon under about 150
bases is worth checking against the genome by hand.

## Running it again

```bash
# with the 1.7.0 work folder in place, holding runs/<target>/{results,inclusion,exclusion}
for target in genus_as_labelled_p1.0 genus_wide_as_labelled_p1.0 \
              abortus_corrected_p1.0 melitensis_corrected_p0.994; do
    primer-finder design "work/runs/${target}/results" -o "rerun/${target}" \
        -i "work/runs/${target}/inclusion" -e "work/runs/${target}/exclusion" \
        --product-size 70-150 --assays-per-region 3 --max-regions 50 -t 48
done
```
