# 2026-10-10 — primer-finder 1.5.0, a second answer key with large groups

The *Xylella* validation asks the right question of small groups: 8 inclusion genomes against 17. This
record asks it of groups one to two orders of magnitude larger, and of a harder target — clonal complexes
*within* one species rather than subspecies.

The answer key is

> Félix B., Capitaine K., Te S., Felten A., Gillot G., Feurer C., van den Bosch T., Torresi M.,
> Sréterné Lancz Z., Delannoy S., Brauge T., Midelet G., Leblanc J.-C., Roussel S. (2023)
> *Identification by High-Throughput Real-Time PCR of 30 Major Circulating Listeria monocytogenes Clonal
> Complexes in Europe.* Microbiology Spectrum 11(3):e03954-22.
> [doi:10.1128/spectrum.03954-22](https://doi.org/10.1128/spectrum.03954-22)

Their 34 TaqMan primer-and-probe sets were designed on kmers and point mutations found in 954
*L. monocytogenes* genomes — at least 15 per clonal complex — and checked in silico against a further
2,388. That is the same question primer-finder answers, asked with different code, which is what makes it an
answer key rather than a comparison.

primer-finder 1.5.0, KMC 3.2.4, SKESA 2.5.1, minimap2 2.31, blast 2.17.0+; mlst 2.33.1 for the typing.

## The panel

Every complete RefSeq *L. monocytogenes* assembly on 2026-10-10: **784 genomes**. Each was given a sequence
type by `mlst` against the Institut Pasteur `listeria_2` scheme, and the clonal complex and lineage that
scheme assigns to that ST. **763 genomes** fell in one of **79 clonal complexes**; the other 21 carry an ST
the scheme places in no complex, or does not know, and are left out of both groups — they might be the
target, and putting them in the exclusion group would throw away real regions.

Complete rather than draft assemblies, for the same reason as the *Xylella* set: a region broken across
contigs looks absent, which is a property of the assembly and not of the genome.

13 clonal complexes have at least 15 complete genomes. Four of them — CC1, CC14, CC37, CC121 — are left out
of the recovery test on purpose: the paper gives each of them **two** assays, which resolve subdivisions
*within* the complex, so a complex-wide inclusion group is not what those assays detect. That leaves **nine
complex-wide assays** to test. The group sizes dwarf the *Xylella* ones:

| | *Xylella* (2026-10-08) | *Listeria* (here) |
|---|---|---|
| inclusion genomes | 8 | **18 to 86** |
| exclusion genomes | 17 | **677 to 745** |
| groups to tell apart | 5 subspecies | **79 clonal complexes of one species** |

## What primer-finder had to beat: the census

Before asking whether primer-finder reports the right region, `census.py` asks where each published assay
actually *is*, by blasting all 102 oligos against all 763 genomes and counting only full-length,
mismatch-free hits. This is the ceiling on what could be recovered: primer-finder's rule is "in every
inclusion genome, different from every exclusion genome", so an assay missing from one genome of its own
complex, or present in one genome of another, cannot be reported — and should not be.

| Assay | in its own complex | in another complex |
|---|---|---|
| CC2 | 43 / 44 | 1 (GCF_051797645.1, CC1/ST1) |
| CC3 | 18 / 18 | 1 (GCF_002101275.1, CC489) |
| CC4 | 18 / 19 | none |
| CC5 | 22 / 22 | none |
| CC6 | 36 / 36 | none |
| CC7 | 47 / 47 | none |
| CC8 | 52 / 52 | none |
| CC9 | 85 / 86 | 1 (GCF_049078365.1, CC121/ST121) |
| CC224 | 48 / 48 | 1 (GCF_041764955.1, CC226/ST226) |

**Three of those four cross-reactions are ones the paper itself reports or predicts.** Its Table 2 lists
CC489 among the in-silico-predicted cross-reactions of the CC3 assay, and ST226a among those of the CC224
assay; the CC9 assay is the one whose analytical sensitivity the paper puts at **0.93**, the lowest of its
complex-wide assays. Finding the same things from a different panel is the strongest check available here
that the census is measuring what it claims to.

Two other observations from the census, reported because they are the kind of thing a panel this size
turns up:

- all 13 non-CC1 genomes that carry the **CC1** assay are **CC183**, which is exactly the cross-reaction the
  paper reports for CC1 SD_2;
- **GCF_051797645.1** is typed CC1 but lacks the CC1 assay and carries the CC2 assay. One genome behaving
  like that is worth no conclusion, but it is where I would look first for a mistyped or recombinant
  assembly.

## The result

primer-finder was run with its defaults for each complex: inclusion = that complex, exclusion = all 78
others.

| Target | inclusion | exclusion | regions reported | published assay recovered | its rank |
|---|---|---|---|---|---|
| CC2 | 44 | 719 | — | no | — |
| CC3 | 18 | 745 | 16 | no | — |
| CC4 | 19 | 744 | 2 | no at the default | — |
| **CC4** at `-p 0.94` | 19 | 744 | 6 | **yes** | **1 of 6** |
| **CC5** | 22 | 741 | 22 | **yes** | 9 of 22 |
| **CC6** | 36 | 727 | 24 | **yes** | 5 of 24 |
| **CC7** | 47 | 716 | 15 | **yes** | **1 of 15** |
| **CC8** | 52 | 711 | 22 | **yes** | 2 of 22 |
| CC9 | 86 | 677 | 11 | no | — |
| CC224 | 48 | 715 | 12 | no | — |

**Five of the nine published assays were recovered**: four at the defaults, and CC4 once `-p 0.94` allowed
one of its 19 genomes to lack the region. Recovered means all three oligos in one reported region, the
primers facing each other, the probe between them, and the amplicon the published length.

In every one of the four that were not recovered, the census says the published assay is **not**
simultaneously present in all of its own complex and absent from all others on this panel. The rejections
are therefore the rule working, not the rule failing:

- **CC3, CC224** — the assay is in every genome of its complex but also in one genome of another. The
  exclusion check drops the region, and the paper predicts both of those cross-reactions.
- **CC9** — missing from one of 86, and present in one CC121 genome. No `-p` can rescue it: the exclusion
  hit is fatal on its own.
- **CC2** — missing from one of 44 and present in one CC1 genome. Separately, no region survived at all:
  only 3,564 inclusion-specific 99-mers assembled into 2 contigs, and neither held two differences within
  21 bases of each other, so the run ended with the documented "no contig differs enough" error rather than
  an empty table.

Each run took between 33 seconds and about two minutes on 32 threads, the largest being 86 inclusion against
677 exclusion genomes.

## Two limitations this panel exposed

The *Xylella* set never showed either, because its regions are long. Both are about primer-finder, not about
the paper.

### An inclusion-specific island shorter than 200 bases is invisible

SKESA will not report a contig under **200 bp**, and primer-finder does not override that. At a kmer size
of 99 that is a real floor: an island covered by fewer than about 100 overlapping kmers is dropped before
anything else sees it. Clonal-complex islands are short, so this bites here.

Demonstrated on CC121 at `-p 0.93`: the published probe and primers **are** among the inclusion-specific
99-mers (77 of the 49,487 contain the probe), and are **not** in the assembly. Re-running SKESA on the same
kmers with `--min_contig 120` gave **349 contigs instead of 67**, the shortest 121 bp, including a 145 bp
island carrying the CC121 probe and reverse primer. Nothing downstream could have recovered that assay
because the region never reached the candidate step.

### A published amplicon may straddle the edge of an island

That same 145 bp island holds the CC121 probe and reverse primer but **not** its forward primer: the
forward primer sits in sequence that is not CC121-specific. primer-finder reports only stretches where
*every* kmer is inclusion-specific, whereas a published assay needs only its amplicon to be present in the
group and distinguishable outside it — the two definitions do not have to agree, and at this level they
often will not.

## What this does and does not show

- It shows primer-finder works unchanged at a scale the *Xylella* set never tested: 86 inclusion against
  677 exclusion genomes, in a species where the groups differ by very little, in under two minutes.
- It shows its accept/reject decisions agree with an independent blast census of the panel in **all nine**
  cases, and that three of its four rejections match cross-reactions the paper reports or predicts.
- It shows `-p/--min-inclusion` doing the job it exists for: one divergent genome out of 19 is the
  difference between no result and the published region at rank 1.
- It does **not** show that primer-finder would have designed those assays — only that it reports the
  regions they sit in. The design step is not exercised here at all.
- It does **not** show that the four unrecovered assays are bad. They were validated at a bench on 597
  strains; this panel is a different, smaller, differently biased sample (complete RefSeq assemblies
  over-represent reference and outbreak strains), and a disagreement about specificity between two genome
  panels is evidence about the panels.
- Nothing here has been near a bench.

One caveat on the census itself: blast counts an IUPAC code as a mismatch, so the count for SD_CC1 — whose
reverse primer carries an M — is a lower bound, and reads as zero. It is a subdivision assay and outside the
test set in any case.

## Reproducing

```bash
conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
conda activate primer-finder_listeria
bash validation/listeria/run.sh /path/to/work 32 64
```

`run.sh` downloads and types the genomes, runs the census, and then for each of the nine complex-wide
targets builds the groups, runs primer-finder and looks for the published oligos. The panel this record used
is pinned in [`validation/listeria/panel.tsv`](../../listeria/panel.tsv) and the census it produced in
[`census.tsv`](../../listeria/census.tsv), so a later run can be diffed against it — the set of public
genomes grows, and these numbers will move.

The `-p` run and the SKESA demonstration were by hand:

```bash
primer-finder -i <work>/CC4/inclusion -e <work>/CC4/exclusion -o <work>/CC4_p/results -t 32 -m 64 -p 0.94
skesa --fasta <work>/CC121_p0.93/results/1_kmers/inclusion_specific_99mers.fasta \
      --contigs_out /tmp/cc121.fasta --cores 16 --memory 32 --min_contig 120
```
