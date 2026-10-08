# 2026-10-08 — primer-finder 1.1.0, Xylella fastidiosa subspecies assays

Ran `validation/xylella/run.sh` on primer-finder 1.1.0 (tag v1.1.0). All four published subspecies assays
were recovered.

This supersedes [`../2026-10-08_v1.0.0/SUMMARY.md`](../2026-10-08_v1.0.0/SUMMARY.md), which covered three
assays: the subsp. *morus* assay was added here, and with it the three other complete *morus* genomes, so
that *morus* is a credible inclusion group rather than a pair. That also enlarged the exclusion group of the
other three targets by three genomes, which is why their numbers differ a little from the earlier record.
Nothing in the earlier record was wrong; it describes a smaller genome set.

## What was run

25 complete RefSeq genomes (`validation/xylella/accessions.tsv`): 8 subsp. *multiplex*, 6 *fastidiosa*,
5 *morus*, 4 *pauca*, 2 *sandyi*, including the strains the paper used as references. Four runs, each with
one subspecies as the inclusion group and every other subspecies as the exclusion group, default options
(`-k 99 -d 1 -p 1.0 -a skesa`), `-t 16 -m 32`. 2 min 4 s for all four on a 64-core machine.

Programs: kmc 3.2.4, SKESA 2.5.1, minimap2 2.31-r1302, blast 2.17.0+.

## Result

| Assay | Inclusion group | Exclusion | Specific kmers | Contigs | Candidates | In all inclusion | Final regions | Recovered | Rank | Region | Amplicon |
|---|---|---|---|---|---|---|---|---|---|---|---|
| XFM | *multiplex* (8) | 17 | 845,820 | 1,832 | 1,546 | 1,540 | 711 | yes | 1 of 711 | 4,010 bp | 88 bp (published: 88) |
| XFF | *fastidiosa* (6) | 19 | 557,035 | 1,179 | 1,049 | 1,048 | 377 | yes | 5 of 377 | 1,111 bp | 100 bp (published: 100) |
| XFP | *pauca* (4) | 21 | 901,845 | 1,891 | 1,755 | 1,745 | 1,413 | yes | 105 of 1,413 | 488 bp | 154 bp (published: 154) |
| XFMO | *morus* (5) | 20 | 486,296 | 1,031 | 701 | 701 | 286 | yes | 5 of 286 | 1,171 bp | 123 bp (published: 123) |

In every run both primers and the probe of the published assay were found exactly, in a single reported
region, with the primers facing each other, the probe between them, and the reconstructed amplicon exactly
the published length. Three of the four regions were in the top five of the list.

## How the regions compare with the ones the paper reports

Context, not a verdict — see `validation/README.md`. Each recovered region was aligned to the reference
genome the paper gives coordinates in; all four aligned at 100% identity.

| Assay | Reported region, in the paper's reference | The paper's region | Overlap |
|---|---|---|---|
| XFM | M12 (NC_010513.1) 1,823,814–1,827,823 | 1,825,046–1,826,705 (1,660 bp) | the whole of it |
| XFF | M23 (NC_010577.1) 2,476,520–2,477,630 | 2,477,123–2,477,638 (516 bp) | 508 bp of 516 |
| XFP | De Donno (NZ_CP020870.1) 338,130–338,617 | 337,676–338,551 (876 bp) | 422 bp of 876 |
| XFMO | MUL0034 (NZ_CP006740.1) 1,907,980–1,909,150 | 1,908,250–1,908,548 (288 bp) | the whole of it |

The XFMO row needed one reading of the paper: its Table 2 prints the end of that region as "908,548",
which cannot be right, since the assay itself sits at 1,908,399–1,908,502 in the same genome. Taken as
1,908,548 the region is covered whole. The paper states 288 bp for a span of 299, which is left as printed;
`validation/xylella/regions.tsv` records both the reading and the discrepancy.

## What this does and does not show

- The pipeline, with default options and on genomes nobody tuned it for, puts four independently published
  diagnostic assays inside the regions it reports, three of them near the top of the list.
- It says nothing about the hundreds of other regions reported in each run, and nothing about specificity
  beyond these 25 genomes.
- It is not a wet-lab result. The published papers did that part.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32
```

Exits non-zero if any of the four assays is not recovered.
