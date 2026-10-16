# Markers in repeated regions: where they were lost, and what was done

2026-10-16, primer-finder 1.8.1. Prompted by a report that assays targeting repeated regions — chosen for
the sensitivity several copies per genome buys — had given false negatives in earlier testing.

Repeats are worth targeting and worth being careful about: a multi-copy element gives an assay more
template to find, and it is also the shape most likely to be dropped by a pipeline built around
"where does this region differ from the exclusion group".

## The test

[`validation/repeats/make_repeat_example.py`](../../repeats/make_repeat_example.py) writes eight genomes
on one 20 kb backbone. Each inclusion genome carries a 300-base marker **two to four times**, at spread
positions; no exclusion genome carries it at all. Copy number varies between genomes, as a transposon's
would, because that is what decides whether the kmer counting keeps the marker.

## Two barriers, one by design and one a defect

**`-d/--duplication` governs whether the marker's kmers survive at all**, and that is documented
behaviour working as intended. `-d` sets the ceiling on how often a kmer may occur per inclusion genome,
and the default of 1 discards anything repeated. On this dataset:

| `-d` | inclusion-specific kmers | what is assembled |
|---|---|---|
| 1 | 215 | fragments only: two regions of about 100 bases |
| 2 | 216 | the same |
| 3 | 419 | the whole 444-base contig, marker included |
| 4, 8 | 419 | the same |

**`-d` is a ceiling on the pooled count, not a per-genome one.** `pipeline` passes
`max_count = (number of inclusion genomes) x d`, and KMC tests that against the kmer's total count across
the whole group — so what `-d` has to reach is the marker's **mean** copy number, not its highest.
Measured, after an earlier draft of this record claimed the opposite: regenerate with copies
`(2, 2, 2, 6)` and no private changes, so the mean is 3 and the highest 6, and `-d 3` keeps the marker
while `-d 2` leaves only fragments. A run targeting a repeat should set `-d` deliberately; left at the
default it reports fragments of the marker or nothing.

**The defect: a region whose specificity is a block the exclusion group simply lacks was dropped.** With
`-d 3` the marker assembled correctly and the run still reported nothing. The marker's junctions are
inclusion-specific too, so the assembly is flank + marker + flank — a 444-base contig holding the 300-base
marker. Every exclusion genome then hit the two flanks and nothing between:

```
q1-72      evalue 2.87e-33     the leading backbone flank
q373-444   evalue 2.87e-33     the trailing backbone flank
q73-372    no hit at all       the marker
```

Two flanking hits are enough for `genomes_hit` to be non-zero, so the contig was not treated as absent
from the exclusion group. It was sent down the "where does it differ" path instead, the flanks align
exactly, there were no differences, and the region was dropped — while the mapping step two steps earlier
had called the same contig a candidate on the very block that is missing. `blast.shared_variants` says in
its own docstring that a position no genome aligns "only makes the assay more selective", and then cannot
credit it: only positions that *differ somewhere* become candidates there.

**Fixed in 1.8.1** by `blast.unaligned_positions`: a stretch no exclusion genome aligns counts as
differing, because that is what it does. With `-d 3` the marker is now reported as one 444-base region
with its 300 marker bases marked, `design` proposes 28 assays on it, all specific by absence. **17** of them
report `inclusion_copies 2-4` — the copy number that makes a repeat worth targeting — while the other 11
report 1, their amplicon straddling a junction between the marker and the backbone, which is unique to
each copy site.

## What the threshold is for, and how it was chosen

A run has to be at least **40 bases** before it is called an absence, and the reason is not the one an
earlier draft of this record gave. Two things that do *not* justify a floor, both measured:

* **The gaps between one contig's local alignments are not small.** Per genome, across five *Listeria*
  complexes, they run to 137 bases and beyond — contig insertions, one stretch replaced by another, and
  HSP pairs megabases apart. The 1-to-6-base runs quoted earlier were the marked runs in the *output*,
  which is a union over every exclusion genome; what keeps the per-genome gaps from being read as
  absences is that a stretch has to be missing from all ~700 of them, not the floor. With a small
  exclusion group they would be credited in full.
* **Being long enough to hold a primer is not the test.** The binding constraint is whether the search
  would have reported the stretch had the genome held it, and that is longer than a primer.

What the floor is for, measured against a 3 Mb genome at the exclusion search's own `-evalue 1e-10`, with
an exact island flanked by unrelated sequence:

| island | evalue | reported? |
|---|---|---|
| 25 bases | 1.16e-06 | no |
| 30 bases | 8.15e-09 | no |
| 33 bases | 5.62e-11 | yes |
| 40 bases | 3.27e-14 | yes |

So a run shorter than 33 bases could be present and merely invisible, and an earlier floor of 25 sat
inside that blind spot. 40 leaves room for a larger genome, whose evalues are worse in proportion. The
real absences in these panels run from 384 to 6,348 bases, so the margin costs nothing.

One blind spot the floor does not close: `-task blastn` seeds on 11 contiguous matches, so a block the
exclusion genome holds at 90% identity is reported in full when its mismatches fall at random, and **not
at all** when they fall every tenth base. Such a block would be marked absent in its entirety. Real
divergence is not periodic, but the shape exists.

With the floor in place, what the rule does on real panels is narrower than "unchanged counts" suggests.
Five *Listeria* complexes, each run twice — once as shipped and once with the rule removed — and the
per-genome and union figures measured from the kept exclusion hits
([`listeria_absence.tsv`](listeria_absence.tsv)):

| complex | exclusion genomes | longest per-genome gap | gaps over 40 | union absences ≥ 40 | regions, with / without |
|---|---|---|---|---|---|
| CC1 | 681 | 21 | 0 | 0 | 8 / 8 |
| CC3 | 745 | 48 | 83 | 0 | 43 / 43 |
| CC7 | 716 | **1,584** | 15 | 5 (445 to 6,348 bases) | 60 / 60 |
| CC224 | 715 | 919 | 58 | 1 (384 bases) | 36 / 36 |
| CC37 | 742 | 137 | **970** | 0 | 48 / 48 |

Three things follow, and the first two matter more than the third.

**The union is what protects, not the floor.** CC37 has 970 per-genome gaps longer than the 40-base floor
and not one union absence; CC3 has 83 and none. A gap in one genome says nothing, because the stretch has
to be missing from all 700-odd of them. With a small exclusion group those gaps would be credited in full,
and that is the case to be careful about — not these panels.

**It adds and removes no regions on any of the five.** On CC1, CC3 and CC37 there is nothing for it to act
on, so their unchanged counts are not evidence that it behaves well. On CC7 and CC224 it does fire, and the
counts are still identical, because the contigs carrying those absences were already being kept by the
differences rule. What changes is the **marking** inside them, which is what `design` reads to place
primers:

| region | length | marked without the rule | with it |
|---|---|---|---|
| CC7 `Contig_207_78.6509` | 8,142 bp | 165 | **6,513** |
| CC7 `Contig_450_75.8999` | 909 bp | 13 | 501 |
| CC224 `Contig_277_74.5974` | 646 bp | 11 | 395 |

**Where it changes an outcome is the case it was built for.** The repeat dataset goes from no region to one
holding its marker, and the bundled example — 1 region, 11 assays specific by absence — is untouched.

The *Listeria* figures are against 1.8.1 with the change reverted, measured for this purpose. The numbers
recorded in `2026-10-14_v1.7.0_unitig/` for the same complexes are 1 and 16, which is a separate matter:
those predate `-task blastn`, whose effect on this panel that record already describes.

## An exact repeat does not assemble at all

Measured while checking the dataset, and worth knowing on its own. Give every copy the identical sequence
— no divergence between them — and the run reports **one region of 99 bases, which is `-k` exactly**,
instead of the 444-base contig above. An exactly repeated block branches the de Bruijn graph at each of
its junctions, so the assembler cannot extend out of the repeat and returns a single kmer's worth of it.

So a repeat-targeted marker needs its copies to differ *somewhere* for the assembler to resolve more than
one kmer of it. Real copies of an element do diverge, which is why the dataset leaves the stray private
change that lands inside one copy in place rather than scrubbing it. It also means `-a unitig` is worth
trying on a repeat target: the walk stops at a branch for the same reason, but it is reached when SKESA
returns almost nothing, and the two fail differently.

## What is still worth knowing

The bundled example's contig is aligned end to end by every exclusion genome — blast gaps through its
15-base insertion — so it has no unaligned runs and never exercised this path. A novel block large enough
to break the alignment in two is what does, and nothing in the test suite or the example had that shape
before this. The repeat dataset is kept for that reason.

One caveat about what gets reported. The marked positions become the region's "variants", and `design`
counts them like any other difference — so an assay inside a 300-base absent block comes out with
`best_primer_variants 20`, `forward_terminal_run 20` and a large `primer_variant_weight`, which are
counts over bases that are *absent* rather than differing. The assay is labelled `absence` and the
numbers order such assays sensibly among themselves, but `best_terminal_run 20` is not a claim about
twenty 3'-end mismatches and should not be read as one. `strong_variants` and the per-oligo variant bases
are empty for these, which is correct.

This says nothing about how such an assay behaves in a reaction. What it establishes is that a
repeat-targeted marker now survives the pipeline that was discarding it, and that `-d` must be set to the
copy number for its kmers to reach the assembler in the first place.
