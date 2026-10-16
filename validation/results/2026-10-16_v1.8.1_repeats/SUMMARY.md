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

So **`-d` must be at least the highest copy number in any inclusion genome**, which is 4 here and was
satisfied at 3 because the marker's interior kmers occur at most that often in the genomes that carry it.
A run targeting a repeat should set `-d` deliberately; left at the default it reports fragments of the
marker or nothing.

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
with its 300 marker bases marked, and `design` proposes 28 assays on it, all specific by absence, every
one reporting `inclusion_copies 2-4` — the copy number that makes a repeat worth targeting.

## What the threshold is for, and how it was chosen

A run has to be at least 25 bases — the widest a primer may be, so that a primer could sit wholly inside
it — before it is called an absence. Measured, not guessed: blast's local alignments of one contig do not
quite abut, and crediting every gap between them turned a *Listeria* CC3 run from 43 regions into 44 and
marked runs of **1 to 6 bases** across both complexes tested. None of those was a real absence. 25 is far
clear of the largest gap seen and is the shortest run that could carry a primer on its own.

With the threshold in place the change is confined to what it was for:

| | before | after |
|---|---|---|
| the repeat dataset, `-d 3` | 0 regions | **1**, holding the marker |
| *Listeria* CC1 | 8 regions | 8 |
| *Listeria* CC3 | 43 regions | 43 |
| the bundled example | 1 region, 11 absence assays | unchanged |

The *Listeria* figures are against 1.8.1 with the change reverted, measured for this purpose. The numbers
recorded in `2026-10-14_v1.7.0_unitig/` for the same complexes are 1 and 16, which is a separate matter:
those predate `-task blastn`, whose effect on this panel that record already describes.

## What is still worth knowing

The bundled example's contig is aligned end to end by every exclusion genome — blast gaps through its
15-base insertion — so it has no unaligned runs and never exercised this path. A novel block large enough
to break the alignment in two is what does, and nothing in the test suite or the example had that shape
before this. The repeat dataset is kept for that reason.

This says nothing about how such an assay behaves in a reaction. What it establishes is that a
repeat-targeted marker now survives the pipeline that was discarding it, and that `-d` must be set to the
copy number for its kmers to reach the assembler in the first place.
