# 2026-10-08 — primer-finder 1.0.0, Xylella fastidiosa subspecies assays

> Superseded by [`../2026-10-08_v1.1.0/SUMMARY.md`](../2026-10-08_v1.1.0/SUMMARY.md), which adds the subsp.
> *morus* assay and three more *morus* genomes. Nothing below is wrong; it describes a 22-genome set, and
> the accession list in the repository now holds 25, so the numbers here are no longer what a fresh run
> produces.

Ran `validation/xylella/run.sh` on primer-finder 1.0.0 (commit 371a455, tag v1.0.0). All three published
assays were recovered.

## What was run

22 complete RefSeq genomes (`validation/xylella/accessions.tsv`): 8 subsp. *multiplex*, 6 *fastidiosa*,
4 *pauca*, 2 *morus*, 2 *sandyi*, including the strains the paper used as references. Three runs, each with
one subspecies as the inclusion group and every other subspecies as the exclusion group, default options
(`-k 99 -d 1 -a skesa`), `-t 16 -m 32`. 80 s for all three on a 64-core machine.

Programs: kmc 3.2.4, SKESA 2.5.1, minimap2 2.31-r1302, blast 2.17.0+.

## Result

| Assay | Inclusion group | Specific kmers | Contigs | Candidates | In all inclusion | Final regions | Assay recovered | Rank of its region | Region size | Amplicon |
|---|---|---|---|---|---|---|---|---|---|---|
| XFM | subsp. *multiplex* (8) | 845,952 | 1,831 | 1,546 | 1,540 | 711 | yes | 1 of 711 | 4,010 bp | 88 bp (published: 88) |
| XFF | subsp. *fastidiosa* (6) | 557,214 | 1,179 | 1,049 | 1,048 | 370 | yes | 5 of 370 | 1,111 bp | 100 bp (published: 100) |
| XFP | subsp. *pauca* (4) | 901,845 | 1,891 | 1,755 | 1,745 | 1,412 | yes | 105 of 1,412 | 488 bp | 154 bp (published: 154) |

In each run the two primers and the probe of the published assay were found exactly, in a single reported
region, with the primers facing each other, the probe between them, and the reconstructed amplicon exactly
the published length. Two of the three regions were in the top five of the list, which is what the sorting
(most differing bases first) is meant to do.

## How the regions compare with the ones the paper reports

Context, not a verdict — see `validation/README.md` for why. Each recovered region was aligned to the
reference genome the paper gives coordinates in; all three aligned at 100% identity.

| Assay | Reported region, in the paper's reference | The paper's region | Overlap |
|---|---|---|---|
| XFM | M12 (NC_010513.1) 1,823,814–1,827,823 | 1,825,046–1,826,705 (1,660 bp) | 1,660 bp, the whole of it |
| XFF | M23 (NC_010577.1) 2,476,520–2,477,630 | 2,477,123–2,477,638 (516 bp) | 508 bp of 516 (98%) |
| XFP | De Donno (NZ_CP020870.1) 338,130–338,617 | 337,676–338,551 (876 bp) | 422 bp of 876 (48%) |

The regions differ in extent from the published long-mers, in both directions, which is expected: the
exclusion genomes are not the same 58 sequences, and a region here is a contig assembled from
inclusion-specific kmers, not a sequence trimmed to where the outgroup starts matching. The XFP region is
the clearest case — it holds the assay and 422 bp of the published region, and stops where this exclusion
set stops being different.

## What this does and does not show

- It shows that the pipeline, with its default options and on genomes nobody tuned it for, puts three
  independently published diagnostic assays inside the regions it reports, and ranks two of them near the
  top of the list.
- It does not show that the regions are good assay targets in the wet lab: that is what the published
  papers did. Nor does it test the assays' specificity against anything outside these 22 genomes.
- Hundreds of other regions are reported per run. They are candidates, not assays; this validation says
  nothing about how many of them would work.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32
```

Exits non-zero if any of the three assays is not recovered.
