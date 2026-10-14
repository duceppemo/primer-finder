# 2026-10-13 — the 200-base floor removed, and the *Listeria* set re-run

The [*Listeria* record of three days ago](../2026-10-10_v1.5.0/SUMMARY.md) found that neither assembler
reports a contig shorter than 200 bases unless told to, and that primer-finder was not telling them. A
candidate region shorter than that was therefore dropped before any later step could see it, silently.

`--min_contig` is now set to the kmer size — the shortest a region can be — and SPAdes, which has no such
option, has the same floor applied to its output. `--min-region` overrides it. This record is the same nine
targets re-run, to say what that changed.

Everything else is the panel of the previous record, unchanged: 784 complete RefSeq *Listeria monocytogenes*
genomes, 763 placed in 79 clonal complexes by the Institut Pasteur `listeria_2` scheme.

## Before and after

Contigs the assembler returned / candidates after the two-differences filter / regions reported:

| Target | before | after | assay before | assay after |
|---|---|---|---|---|
| CC2 | 2 / 0 / **run failed** | 32 / 6 / 1 | — | no |
| CC3 | 49 / 35 / 16 | 485 / 154 / 38 | no | no |
| CC4 | 5 / 4 / 2 | 99 / 30 / 5 | no | no |
| CC5 | 46 / 31 / 22 | 468 / 176 / 49 | yes, rank 9 | yes, rank 10 |
| CC6 | 51 / 42 / 24 | 354 / 130 / 46 | yes, rank 5 | yes, rank 5 |
| CC7 | 52 / 50 / 15 | 533 / 473 / 57 | yes, rank 1 | yes, rank 1 |
| CC8 | 49 / 48 / 22 | 414 / 376 / 59 | yes, rank 2 | yes, rank 2 |
| CC9 | 16 / 16 / 11 | 158 / 147 / 19 | no | no |
| CC224 | 30 / 25 / 12 | 288 / 103 / 32 | no | no |
| CC4 at `-p 0.94` | 13 / 12 / 6 | 107 / 40 / 9 | yes, rank 1 | yes, rank 1 |
| CC121 at `-p 0.93` | 67 / 64 / 25 | 398 / 365 / 66 | no | no |

**Which published assays are recovered did not change**: the same five, at the same ranks but one. CC5
slipped from 9 to 10 because there are now more regions above it, which is what more candidates means. So
this is a blind spot removed, not a result improved — and that is worth saying plainly, because the
temptation with a change like this is to present a bigger number as a better one.

Two things did change:

- **CC2 no longer fails.** Before, its 3,564 inclusion-specific kmers assembled into 2 contigs, neither of
  which carried two differences within 21 bases, so the run ended with the "no contig differs enough" error
  and produced nothing at all. It now assembles 32 contigs, 6 of them candidates, and reports one region.
  The published CC2 assay is still not recovered — the census says it is missing from one of its 44 genomes
  and present in a CC1 genome — but the run now has an answer instead of an error.
- **Short regions are a large share of what is reported.** For CC121 at `-p 0.93`, 41 of the 66 regions are
  shorter than 200 bases. Those were all invisible before.

## What happened to the CC121 island

The previous record used CC121 to demonstrate the floor: the published probe sits in 77 of its
inclusion-specific 99-mers and in none of the 67 contigs the assembler returned. Traced through the fixed
pipeline, that island now:

- **is assembled** — `Contig_101_56.952`, 145 bases, holding the published probe and reverse primer;
- **passes** the two-differences filter and the inclusion presence check, appearing in `3_candidates/` and
  `4_blast/all_inclusion_contigs.fasta`;
- **is turned down by the exclusion check**, and so is not among the 66 reported regions.

That is the right shape of answer: the floor was hiding the stretch from every later step, and now the later
steps get to decide. This one they decide against, on the sequence rather than on its length. The second
limitation the previous record describes is untouched and still stands — the island holds the probe and the
reverse primer but not the forward primer, so that published amplicon straddles the edge of the specific
stretch, which is a difference between what primer-finder reports and what a published assay needs.

## What this does and does not show

- It shows the floor was real, and that removing it puts 6 to 10 times as many stretches in front of the
  filters, 2 to 3 times as many of which survive to be reported.
- It shows removing it changed no verdict on this panel: five of nine recovered, before and after.
- It shows one run that previously failed outright now completes.
- It does **not** show that the extra regions are good ones. They are shorter, so they have less room for a
  primer pair, and nothing here checks whether an assay can be designed in them — `design` was not run.
- Nothing here has been near a bench.

## Reproducing

```bash
bash validation/listeria/run.sh /path/to/work 32 64          # the nine targets, with the floor at the kmer size
primer-finder -i <work>/CC121/inclusion -e <work>/CC121/exclusion -o <work>/CC121_trace \
    -t 32 -m 64 -p 0.93 --keep-intermediate                  # to follow the 145-base island through
primer-finder ... --min-region 200                           # to get the old behaviour back
```
