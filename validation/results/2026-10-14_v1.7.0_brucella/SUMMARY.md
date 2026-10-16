# Brucella: the published genus and species assays, primer-finder 1.7.0

2026-10-14, primer-finder 1.7.0, against

> Hoover E., Chmara J., Duceppe M., Charron P., Lloyd D., Huang H., Surujballi O.,
> Thanthrige-Don N., Arnold K., Hazelwood J., Kang M. (2026) *Two novel real-time PCR assays for
> Brucella detection: a species-specific multiplex and a genus-level singleplex developed via
> large-scale whole genome sequencing data analysis.* Microbiology Spectrum 14:e02974-25.
> doi:10.1128/spectrum.02974-25

This is the first record where the published assays were found with **this software**: the paper's unique
core sequences came from `selective_primer_finder 0.2.1` and RUCS 1.0.3, one of its eleven authors being
this package's maintainer. primer-finder is the rewrite of that 0.2.1, so the question is whether the
rewrite still finds what its ancestor found.

**Short answer.** The genus assay is recovered; the three species assays are not, and the reason is a rule
the rewrite has and 0.2.1 did not. Every number below is in silico — about sequence, not about how any
assay behaves in a reaction.

## The panel

The paper's accessions are not published (supplementary Table S2 gives only a count per species), so the
panel is **rebuilt by the same rule** from today's RefSeq rather than reproduced: every assembly of the
eight classical species its panel holds, and of the four non-classical groups it used as negatives.

| | this record | the paper |
|---|---|---|
| classical *Brucella* | **1,376** (398 abortus, 734 melitensis, 120 suis, 74 canis, 21 ceti, 18 ovis, 6 neotomae, 5 microti) | 1,160 |
| non-classical | **315** (former *Ochrobactrum*, *Pseudochrobactrum*, *Paenochrobactrum*, *Falsochrobactrum*) | 194 |
| held out | 103 (86 `Brucella sp.`, 8 pinnipedialis, 4 amazoniensis, 3 inopinata, 2 vulpis) | — |
| assembly level, classical | 201 complete, 132 chromosome, 170 scaffold, **873 contig** | drafts included |

RefSeq has grown since October 2020, so the groups are larger. The panel is draft-heavy in the same way
the paper's was, which matters: 873 of 1,376 are contig-level.

The paper's Table S5 revisions are applied before grouping — grouping on a species label the paper itself
corrected would put a genome in the wrong group. `panel.tsv.gz` records what was used.

## What the census says, before primer-finder runs at all

All four published amplicons were first measured in a reference genome of their own group: **93, 75, 147
and 102 bp, every one exactly as published** (`verify_amplicon.py`). The answer key is sound.

Then every oligo was blasted against all 1,794 genomes, counting only full-length, mismatch-free hits
(`census.tsv`):

| assay | holds all three oligos | of its own group | exclusion genomes holding it |
|---|---|---|---|
| Bru_genus | **1,374 of 1,376** | 99.85% | **0 of 315** |
| BA | 391 of 398 | 98.2% | 1 |
| BM_10 | 728 of 734 | 99.2% | 4 |
| BS | 114 of 120 | 95.0% | 3 |

Three independent confirmations of the paper fall out of this:

1. **The two genomes that lack the genus assay are the two the paper names** — GCF_002290125.1 and
   GCF_023651695.1, its only in-silico misses (99.8% sensitivity, "more than three mismatches at primer
   binding sites"). Found here on a panel assembled from scratch, 1,376 genomes instead of 1,160.
2. **Both over-fragmented *B. suis* assemblies reproduce Table S4 contig for contig**: in GCF_000292005.1
   the probe (24/24) and reverse primer (20/20) sit on NZ_ALOK02000260.1 at 539-562 and 497-516 while the
   forward primer is cut by a contig edge at 21 of its 22 bases on NZ_ALOK02000046.1; GCF_000292105.1 is
   the same picture on NZ_ALON02000226.1 and NZ_ALON02000030.1. Identical to the published table.
3. **No genome of the 1,794 carries two species assays.** The multiplex is mutually exclusive across the
   whole panel.

### Eight more mislabelled assemblies, by the same argument the paper used

The eight "exclusion" hits above are not cross-reactions. Each of those genomes carries all three oligos
of **exactly one** species assay, and it is not the assay of its own RefSeq label, while carrying the
genus assay as well. A cross-reacting assay would show a genome carrying its own species' assay *and*
another's; none does.

| accession | RefSeq label | the assay it carries | level |
|---|---|---|---|
| GCF_055389475.1 | melitensis | abortus | complete genome |
| GCF_016091955.1 | abortus | melitensis | contig |
| GCF_016091965.1 | abortus | melitensis | contig |
| GCF_018604785.1 | abortus | melitensis | contig |
| GCF_041861985.1 | abortus | melitensis | contig |
| GCF_000292045.1 | abortus | suis | contig |
| GCF_055388155.1 | melitensis | suis | complete genome |
| GCF_000662015.2 | canis (per Table S5) | suis | contig |

That extends the paper's Table S5 from six assemblies to fourteen, two of them complete genomes deposited
since.<a id="correction-1"></a>[^corr1]

[^corr1]: **Correction, 2026-10-16.** Fourteen double-counts `GCF_000662015.2`, which is listed among the
    eight above *and* is one of the paper's own six Table S5 revisions, as this record itself says two
    paragraphs down. Seven of the eight are new to the table, so the total is **thirteen**. Everything
    else here stands; the eight-assembly list and every count derived from it are unaffected, since the
    assembly belongs in the list either way -- it is why the *B. suis* assay shows one exclusion hit on
    the panel as labelled. Four of the paper's own six revisions are independently confirmed by which assay the genome
carries (suis→melitensis, melitensis→suis ×2, melitensis→abortus). Of the other two:

- **GCF_000292145.1** (abortus→melitensis, 292 contigs) carries no species assay at all: the melitensis
  forward primer and probe are intact and adjacent on NZ_ALOL02000099.1, but the reverse primer site is
  not, so the reassignment is supported as far as this evidence goes and the assay simply cannot be read
  off this assembly;
- **GCF_000662015.2** (inopinata BO1→canis) carries the *B. suis* assay. A 102 bp amplicon is far weaker
  evidence than the Mash, SNP, PubMLST and phylogenetic analysis behind Table S5, so this is recorded as
  an observation, not a correction. It is also worth knowing that RefSeq holds a **second** assembly of
  strain BO1, GCF_000182725.1 (55 contigs, 3.37 Mb against 3.21 Mb), and the two disagree about this
  locus — one of them is not what its label says.

Because leaving these eight in the exclusion group removes regions that belong to the species they really
come from, every target below was run twice: on the panel **as labelled**, and with the eight **held out**.

### Two things the paper's panel could not show

- The genus assay's oligos are present without mismatch in **all 8 *B. pinnipedialis*, all 4
  *B. amazoniensis* and both *B. vulpis* assemblies**, and in 49 of the 86 unidentified `Brucella sp.` —
  none of which the paper's panel held. It separates that unidentified set roughly in half, which is what
  one would expect of a genus assay given that about half of them are former *Ochrobactrum*.
- In **all three *B. inopinata* assemblies** the amplicon is present, the forward primer site is intact
  (19/19, no mismatch), the reverse primer site is intact or one base short — and the probe carries
  **exactly one mismatch**. Whether an assay tolerates that is a wet-lab question this says nothing about;
  the sequence observation is what is recorded.

## What primer-finder reports

Every run used the default `find` settings except where stated, `-t 48 -m 128`. `runs.tsv` and
`ksweep.tsv` hold all 41 runs.

| target | panel | `-p` | regions | published assay |
|---|---|---|---|---|
| genus | as labelled | 1.0 (default) | 398 | not among them |
| genus | as labelled | 0.998 (≥1374) | 2,139 | not among them |
| genus | as labelled | **0.99** (≥1363) | 5,344 | **recovered**, rank 4,249, in a 127 bp region |
| genus | corrected | 1.0 | 412 | not among them |
| genus | corrected | 0.998 | 2,167 | not among them |
| abortus | as labelled | 1.0 | — | 41 specific 99-mers, too few for the assembler |
| abortus | corrected | 1.0 / 0.994 | 1 / 2 | not among them |
| melitensis | as labelled | 1.0 / 0.991 | — | no inclusion-specific kmer at all |
| melitensis | corrected | 1.0 / 0.994 | 1 / 11 | not among them |
| suis | as labelled | 1.0 / 0.95 | — | no inclusion-specific kmer at all |
| suis | corrected | 0.95 | 0 | — |

### The genus assay: recovered, and only because of a change made in 1.6.0

`Contig_1467_53.243`, **127 bp**, holding the published amplicon at exactly 93 bp. Until 1.6.0 both
assemblers were left at their own 200 bp minimum contig length, so a 127 bp region could not be reported
at all: this assay would have been unrecoverable in 1.5.0 whatever the settings. The *Listeria* panel
showed that removing that floor changed no verdict; here it is the difference between recovering a
published assay and not.

Why `-p 0.99` and not the 0.9985 the census implies: the census asks whether the three oligos are present,
while the kmer step asks for an intact 99-mer across the locus. About eleven further genomes differ
somewhere in that window without touching an oligo. Lowering `-k` to 71, 51 or 31 does not substitute for
the relaxation (`ksweep.tsv`) — the variation is inside the amplicon, not in its flanks.

Having found the region, `design` was run on it alone (`check_design.py`, `design_genus.json`). It
proposed 16 assays, and the closest:

```
published   F CGCGCATGATACCACTTCT   R GTGTTCGCGCATGTTGAAG   P AGGTGCGGTCGATAAACTTCTGGAAC    93 bp
assay5      F  GCGCATGATACCACTTCTC  R  ATGTGTTCGCGCATGTTG   P CTTCTGGAACCGGGTCTGGTCG        94 bp
assay3      F   GCATGATACCACTTCTCGAA R  ATGTGTTCGCGCATGTTG  P CTTCTGGAACCGGGTCTGGTCG        92 bp
```

Not identical — two Primer3 runs under different constraints will not agree on a 19-mer — but the forward
primer shares 18 of 19 bases with the published one and the amplicon comes out at 92-94 bp against the
published 93. The design step independently puts the assay in the same place.

### The three species assays: not recovered, and the reason is a rule 0.2.1 did not have

The kmer step **does** find these signatures. At `-k 51` on the corrected panel, 16 of the 25 51-mer
windows of the *B. abortus* amplicon are among the 23,413 inclusion-specific 51-mers; at `-k 31`, 19 of
the melitensis amplicon's 117 windows are among 31,256. The specificity is there and primer-finder sees
it. What stops it is the next filter.

primer-finder keeps a candidate region only when **two differences could sit under one primer** — at most
25 bases apart — or when at least 21 bases of a contig end do not align at all. Measured against the
nearest homologue in each assay's own exclusion group:

| assay | against | differences in the amplicon | spacing | the find rule |
|---|---|---|---|---|
| BA | *B. melitensis* 16M | **one**, at amplicon position 16 | — | cannot keep it: one difference |
| BM_10 | *B. abortus* 2308 | **two**, at positions 15 and 129 | **114 bases apart** | cannot keep it: needs ≤25 |
| BS | *B. abortus* 2308 | **none** in the 81 aligned bases | the first 21 bases do not align | cannot keep it as an internal indel |

So all three published species assays rest on discrimination that primer-finder's find rule is designed to
reject: a single mismatch under one primer, or one mismatch under each of two primers 114 bases apart.
That rule is deliberate — two differences under one primer discriminate far better than one each under two
— and it is the rule `--min-primer-differences` defaults to 2 for. **It did not exist in 0.2.1**, which
retained any unique core sequence longer than 28 bp (the paper reports a mean of 59.5 bp). The rewrite is
stricter than the version that produced these assays, and that is the whole of the difference.

Two secondary observations from the same runs:

- **On the panel as labelled, the species targets return nothing at all** — "no inclusion-specific kmer"
  for melitensis and suis at every `-k` and `-p` tried. Holding out the eight mislabelled assemblies turns
  that into 1-11 regions. A single genome filed under the wrong species removes every kmer of the region
  that distinguishes it, which is worth knowing before blaming a tool for finding nothing.
- **At `-k 31` SKESA fails outright** ("Was not able to assemble anything", exit 1) on 15,000-31,000
  specific 31-mers, with `--min_count 1`, `--kmer 21` and `--kmer 25` all making no difference. The advice
  the "too few kmers" error gives — lower `-k` — therefore runs into a different wall. A region shorter
  than roughly 50 bases cannot currently come out of this pipeline whatever the kmer size, because the
  assembly step cannot build it.

  **Since fixed, by `-a unitig`** (`unitig_demo.tsv`): the same 31,256 31-mers give 1,237 contigs and 81
  reported regions, and the 23,413 51-mers give 615 contigs and 4 regions where SKESA gave 5 and none.
  None of those regions holds a published assay — the find rule is still the blocker, as the rest of this
  record sets out — so the verdicts above stand unchanged.

## What would have to change to recover the species assays

Nothing about the panel or the settings: the limit is the find rule, and no `-k`, `-p` or panel correction
gets past it. Recovering this class of assay needs the find rule to be able to keep a region on **one**
difference, or on two differences further apart than one primer — which is to say, the `find` side needs
the equivalent of the `--min-primer-differences 1` that `design` already has, and the span the pair rule
uses would need to become an option rather than a constant. Both are in `mapping.py` as `MAX_PRIMER_SPAN`
and the two-difference test in `is_candidate`, and in `blast.has_close_variants` at the exclusion step.

That is a change to the tool, not to this record, and it is not made here. It is also not obviously the
right default: the rule exists because an assay resting on a single mismatch is the harder one to make
work, which is a judgement about reactions and not about sequence.

## Reproducing this

```bash
conda install -c conda-forge ncbi-datasets-cli
bash validation/brucella/run.sh /path/to/work 48 128
```

About 6 GB of genomes and, on 48 cores, a minute or two per target plus ten minutes for the census. The
panel will differ from this one as RefSeq grows; `panel.tsv.gz` pins what was used here.

## Files

| File | What it holds |
|---|---|
| `panel.tsv.gz` | the 1,794 assemblies used, with species, group, assembly level and contig count |
| `census.tsv` | per assay: the genomes of its own group that hold all three oligos, and of any other |
| `census_detail.tsv.gz` | every genome holding all three oligos of an assay, and whether on one contig |
| `mislabelled.txt` | the eight assemblies whose RefSeq species label the census contradicts |
| `runs.tsv` | the 17 runs at the default kmer size, with region counts and whether the assay came out |
| `ksweep.tsv` | the 24 runs sweeping `-k` over 71, 51 and 31 |
| `design_genus.json` | what `design` proposed on the recovered genus region, against the published oligos |
| `unitig_demo.tsv` | both assemblers on the three runs where they differ most, which is what `-a unitig` was added for |

[`ALTERNATIVES.md`](ALTERNATIVES.md) beside this file asks the other question — what primer-finder proposes
of its own on the same genomes, and whether any of it beats the published assays.
[`WHOLE_CLADE.md`](WHOLE_CLADE.md) then puts the four remaining classical-clade species into the genus
group and asks it again, and [`SUIS_LATO.md`](SUIS_LATO.md) puts *B. canis* into the *B. suis* group.
