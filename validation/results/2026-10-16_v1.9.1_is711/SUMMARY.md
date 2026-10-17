# A published repeat-targeted assay as an answer key: IS711

2026-10-16, primer-finder 1.9.1. The repeats record of the same day
([`2026-10-16_v1.8.1_repeats/`](../2026-10-16_v1.8.1_repeats/SUMMARY.md)) rests entirely on a synthetic
dataset: eight genomes written to have the shape of a repeat-targeted marker. This one asks the same
question of a real published assay in 1,691 real genomes.

## The answer key

> Bounaadja L., Albert D., Chenais B., Hénault S., Zygmunt M.S., Poliak S., Garin-Bastuji B. (2009)
> *Real-time PCR for identification of Brucella spp.: a comparative study of IS711, bcsp31 and per target
> genes.* Veterinary Microbiology 137:156-164.
> [doi:10.1016/j.vetmic.2008.12.023](https://doi.org/10.1016/j.vetmic.2008.12.023)

That paper is the right answer key for this question because its subject *is* the question: it compares a
multi-copy target, the insertion sequence IS711, with two single-copy ones, bcsp31 and per. Its own
results are behind the paywall and are not quoted here; what is on the record publicly is that its IS711
set is the one later work adopts as the genus assay, and the papers reproducing it describe it as the
sensitive one. IS711 is 842 bases (GenBank [M94960](https://www.ncbi.nlm.nih.gov/nuccore/M94960)),
carried several times over by a *Brucella* genome and by nothing in the exclusion group below — which is
what a genus assay wants, and exactly the shape a pipeline built around "where does this region differ"
is most likely to drop.

The panel is the one the *Brucella* record already uses: 1,376 genomes of the eight classical species as
the inclusion group of the `genus` target, against 315 *Ochrobactrum*, *Pseudochrobactrum*,
*Paenochrobactrum* and *Falsochrobactrum* genomes. **No genome of that exclusion group carries IS711 at
all**, so unlike four of the nine *Listeria* complex assays, this one is not ruled out by the panel: it is
specific here, and a run should be able to report it.

**The oligos were checked against a genome before anything was run.** The paper is paywalled, so the three
oligos are taken from the open-access papers that reproduce them and cite it. Eighteen of those in PMC
carry the reverse primer as `GGCCTACCGCTGCGAAT` and one as `CGCCTACCGCTGCGAAT`; the second spelling occurs
in no *Brucella* genome, so it is a typo in that one paper. Locating all three oligos in
GCF_000008145.1 (*B. abortus* bv. 1 str. 9-941) before trusting the table is what caught it — the same job
[`verify_amplicon.py`](../../brucella/verify_amplicon.py) does for the *Brucella* key. The amplicon is
**63 bases**, measured here between the outer 5' ends of the two primers; the paper's own figure could not
be read, so where `check_assay.py` prints "the paper reports 63 bp" it means "measured".

## What is actually in the panel: `is711_census.py`

Whole copies of the element and copies of the 63-base amplicon, counted by blast
(`-perc_identity 90`, 80% and 90% of the query aligned) over every genome of the panel:

| | genomes | element copies: pooled | mean per genome | none at all |
|---|---|---|---|---|
| the eight classical species | 1,376 | 3,977 | 2.89 | 671 |
| the exclusion group | 315 | 0 | 0 | 315 |

and for the 63-base amplicon the assay would see: pooled **5,197** over the classical genomes, mean
**3.78**, none at all in **226** of them. That mean is a lower bound on what `-d` must reach and not the
figure to set it by — see recommendation 1 below, corrected after review: counted exactly, the
highest-counting 99-mer of the element needs 3.11 on this panel, and on the draft group the mean is short.

The copy numbers are right where the literature puts them, in the genomes that can show them:

| species | complete genomes | median element copies | range |
|---|---|---|---|
| abortus | 82 | 6 | 4-8 |
| melitensis | 78 | 8 | 5-10 |
| suis | 26 | 7.5 | 0-13 |
| canis | 9 | 6 | 4-6 |
| ceti | 4 | 40 | 12-40 |
| microti | 1 | 13 | - |
| ovis | 1 | 38 | - |

## Barrier zero, which is not primer-finder's: the assemblies do not hold the repeat

| assembly level | genomes | no whole element copy | median element copies | median amplicon copies |
|---|---|---|---|---|
| Complete Genome | 201 | 2 | 7 | 7 |
| Contig | 873 | 557 | 0 | 1 |
| Scaffold | 170 | 58 | 2 | 3 |
| Chromosome | 132 | 54 | 2 | 2 |

**556 of the 557** contig-level assemblies with no whole copy of the element do hold *pieces* of it, and at
least one piece ends within 50 bases of a contig end. The median longest piece is 422 of the 842 bases.
That is the de Bruijn graph branching at a repeat: the contig stops where the copies become
indistinguishable. Only one of the 873 has no trace of IS711 at all.

So in a draft-heavy public panel — 873 of these 1,376 genomes are contig-level, and 557 of those 873 have
no whole copy — a multi-copy marker is mostly gone before any tool is run, and the assemblies that lost it
are not marked as having lost anything. Across the whole classical group, 671 of 1,376 assemblies hold no
whole copy of the element and 226 hold no copy of the 63-base amplicon. This is the mechanism behind the
report that started this work: an assay targeting a repeat, tested against draft assemblies, gives false
negatives, and the assemblies rather than the assay are what changed.

## Barrier one and two: `-d` and the assembler, each necessary and neither sufficient

Thirty complete classical genomes (up to six per species, IS711 copies 4 to 40, mean 11.63) against thirty
complete non-classical ones, by `is711_panels.py`. On complete genomes the element is intact in every
genome, so nothing but primer-finder's own handling of repeats is in the way. `-d 12` is the mean copy
number of this group rounded up, which is what the ceiling has to reach.

| | `-a skesa` | `-a unitig` |
|---|---|---|
| `-d 1` (default) | 5,074 regions, **not recovered** | 10,064 regions, **not recovered** |
| `-d 12` | 5,443 regions, **not recovered** | 10,414 regions, **recovered in three of them** |

**Both changes are needed and neither is enough on its own**, and the two runs at `-d 12` say which step
is responsible, because they differ in nothing else. Both counted the same 1,329,878 inclusion-specific
99-mers, IS711's among them at 100% identity, and then:

| `-d 12` | longest piece of the element in `2_assembly/assembly.fasta` |
|---|---|
| `-a skesa` | 19 bases — nothing of it, the assembly has 5,616 records |
| `-a unitig` | 289, 225 and 197 bases -- the longest of 40 contigs holding 100 or more -- in 10,923 |

So at `-d 1` the ceiling discards the marker's kmers whichever assembler follows, and at `-d 12` the
kmers are there and SKESA does not assemble them: near-identical copies branch the graph at every
junction, the pieces are short, and an assembler built for sequencing reads treats short and
low-coverage as error. `unitig` builds every stretch the kmers determine and keeps them.

In the cell that works the three published oligos are found in five regions and are together, facing, with
the probe between them, in three:

| region | length | rank | amplicon |
|---|---|---|---|
| unitig_2208 | 289 bp | 1,227 of 10,414 | 63 bp, primers facing, probe between |
| unitig_6072 | 168 bp | 4,199 of 10,414 | the same |
| unitig_9340 | 116 bp | 7,069 of 10,414 | the same |

Three regions rather than one, and all three are the same stretch of the element: they align to M94960
positions 433-739, two of them at 100% identity and unitig_2208 at 98.6%. They differ from each other —
unitig_2208 against unitig_6072 is three mismatches over the 168 bases they share, while unitig_9340
matches both exactly over the 98 bases it overlaps and ends elsewhere. That is the graph branching where
the copies of the element diverge, which they do: in B. abortus 2308 the reverse primer's site reads
`GGCCTACCGCTGCGAAT` in one chromosome-I copy of *B. abortus* 9-941 and `GGCCTACCGCTGAGAAT` in three
others. A repeat-targeted
region therefore arrives as several short pieces rather than one long one, which is why `--min-region 99`
matters here as much as `-d`.

`design` on those three regions, `--product-size 50-150` because the published amplicon is 63 bases and
the default floor is 70: **every assay reports `inclusion_copies_min 4` and `inclusion_copies_max 40`**,
which is exactly what blast counts in those thirty genomes (4 to 40 whole element copies), with
`exclusion_genomes_with_amplicon 0` and `specific_by absence`. The copy reporting added in 1.8 and 1.9 had
only ever been checked against the synthetic dataset; it is right on a real repeat.

## On a draft-heavy group: recovered, and `design` is the one that says what it cost

The second group is fifteen complete genomes and fifteen contig-level assemblies, which is what a public
panel looks like. (It is sampled from the same pool rather than from the thirty above, so only five of its
complete genomes are among them — corrected after review, where this said "fifteen of those".) Whole element copies: pooled 148, mean 4.93, and **ten of the thirty hold no whole
copy at all**, so `-d 5`.

At `-d 5 -a unitig -p 1.0` the run reports 10,131 regions and the assay **is recovered**, in unitig_1828
(338 bp, rank 922), amplicon 63 bp, primers facing, probe between. That is a weaker result than it looks,
and the design step is what shows why:

| group | `inclusion_copies_min` | `inclusion_copies_max` |
|---|---|---|
| complete genomes, `-d 12` | 4 | 40 |
| the draft-heavy group, `-d 5` | **0** | 38 |

A minimum of 0 means at least one inclusion genome does not hold the amplicon at all. It does not: in
GCF_964187125.1 (*B. canis*, contig-level) the published amplicon is absent, and the whole 338-base region
is represented by **41 identical bases at the very start of a contig** — the only hit under the 1e-10
cut-off, at evalue 2.32e-14, covering 12% of the region. The find stage counted that as the region being
present in that genome, because that is all it asks for.

So on a draft panel the find stage says a repeat-targeted region is in every inclusion genome, and the
design step's copy count is what says it is not. Both were run here; a run that stopped at `find` would
have been misled.

## The whole panel: the assay is assembled and then correctly refused

All 1,376 classical genomes against all 315 non-classical ones, `-k 99 --min-region 99 -a unitig`:

| run | wall time | regions reported | the published assay |
|---|---|---|---|
| `-d 1 -p 1.0` | 4 min, 48 threads | 2,028 | not recovered |
| `-d 6 -p 1.0` | 6.5 min, 48 threads | 2,888 | not recovered |
| `-d 6 -p 0.85` | 29 min, 40 threads | 7,014 | **recovered in three of them** |

The last one is slower because lowering `-p` lowers KMC's `-ci` with it, so more kmers survive and more
candidates go through the per-genome checks: 7,341 contigs against 1,376 genomes rather than 5,542.

`-d 6` is above both the 3.78 the amplicon census suggests and the 3.11 the element's highest-counting
99-mer actually needs here. Following the `-d 6 -p 1.0` run through its own
intermediate folders says exactly where the marker goes:

| stage | records | holding 100+ bases of IS711 |
|---|---|---|
| `2_assembly/assembly.fasta` | 5,561 | 14 |
| `3_candidates/best_kmers.fasta` | 5,542 | 14 |
| `4_blast/all_inclusion_contigs.fasta` | 3,184 | **0** |

So with the ceiling raised the kmers survive, assemble, and reach the candidate stage — four of those
contigs hold the published 63-base amplicon at 100% identity — and the **"present in every inclusion
genome" check drops them**. That is the correct answer to the question as asked. The run's own
`inclusion_blast_hits.tsv` says how close the four amplicon-bearing ones came:

| contig | inclusion genomes holding it |
|---|---|
| unitig_1596 | 1,213 of 1,376 (0.882) |
| unitig_4695 | 1,181 of 1,376 (0.858) |
| unitig_3892 | 1,176 of 1,376 (0.855) |
| unitig_5498 | 1,168 of 1,376 (0.849) |

**Corrected after review**: this table first listed unitig_1751 at 1,119 in place of unitig_5498.
unitig_1751 holds 163 bases of the element but not the amplicon, so the four that carry the assay span
1,168 to 1,213, not 1,119 to 1,213. The 14 contigs holding 100 or more bases of the element span 1,119 to
1,359, eight of them above 96% of the panel — a weaker version of "most assemblies have lost the repeat"
than the four amplicon contigs alone suggest. unitig_5498 at 1,168 is also why the `-p 0.85` run reports
three regions and not four: the quorum there is 1,170.

which is the census number from the other direction: 1,150 of the 1,376 assemblies hold the amplicon
(0.836). A panel most of whose draft assemblies have lost the repeat cannot support a marker inside it at
`-p 1.0`, and nothing in the tool could rescue that — the sequence is not in the files.

At `-p 0.85` — which asks for 1,170 of the 1,376 — the published assay **is** recovered on the whole
panel, in unitig_3509 (338 bp), unitig_6585 (128 bp) and unitig_6856 (116 bp), each with the 63-base
amplicon, primers facing and the probe between them. The run took 29 minutes on 40 threads. The two
short ones carry `inclusion=1182/1376` and `inclusion=1176/1376` in their headers, the fraction `-p` below
1 makes the tool record.

**Where they sit in the file is the practical catch.** With `-p` below 1 the regions are ordered by how
many inclusion genomes hold them, so a repeat-targeted region — held by fewer of them than a core one —
lands at the bottom: these ranked 6,748, **7,013 and 7,014 of 7,014**. `design` defaults to
`--max-regions 50`, which on this run would never reach them. A run after a repeat needs
`--max-regions 0`, or to pull the regions out by name, and the ordering is worth knowing before reading
the top of the file and concluding nothing was found.

## One thing to decide: "present in every inclusion genome" counts any hit

`blast.presence_in_genomes` calls a region present in a genome when blast returns **any** hit for it below
the evalue cut-off, with no requirement on how much of the region that hit covers. The draft group above
is the measured case: unitig_1828 is 338 bases, the only hit of it in GCF_964187125.1 under the cut-off is
**41 identical bases at the start of a contig**, 12% of the region, and the published amplicon is not in
that assembly at all. The run reported the region as present in all thirty inclusion genomes.

Nothing downstream is fooled -- `design` reports `inclusion_copies_min 0` for every assay on that region,
which is what says an inclusion genome lacks the amplicon -- so the pipeline's answer is right and it is
the find stage's own log line, "N contig(s) are present in all inclusion genomes", that claims more than
was checked. This is the repeat case making a general property visible: a region partly present in a
genome is commonest where the assembly broke, which is where repeats are.

Two ways to settle it:

- **require a fraction of the region to align** for "present", the analogue of `MIN_AMPLICON_COVERAGE` on
  the design side. It changes what every run reports, so it needs measuring against all the answer keys
  before it could be a default;
- **report it**: carry the worst inclusion genome's coverage into the presence table and the log, and
  leave the decision to the reader. Non-breaking.

**Settled the same day in [`2026-10-16_v1.9.1_presence/`](../2026-10-16_v1.9.1_presence/SUMMARY.md), by
measuring the first option against every answer key: it is the second.** A rule asking for 0.9 of the
region would drop only 6 of the 7,014 regions the full panel reports at `-p 0.85` — and all three regions
carrying this assay are among those 6, since the genomes holding 0.9 of them number 1,115, 1,100 and 1,119
against the 1,170 that run requires. On the complete-genome panels such a rule would gain nothing either:
nothing is partly held there and all eight recovered published assays sit in regions every genome holds
whole. So it would leave a run's regions almost untouched and remove the repeat-targeted answer key. The
fraction is now reported, in the presence table and in a log line, and nothing is kept or dropped
differently.

## If you are after a repeat-targeted marker

What this run says to do, in order:

1. **count the copies first, and count the right thing.** `-d` is a ceiling on a kmer's *pooled* count
   across the group, so what it must reach is the **highest pooled count of a kmer of the target**,
   divided by the number of inclusion genomes. The mean copy number is a convenient lower bound and it
   is **not sufficient** — corrected after review, with this record's own draft group as the
   counterexample. Counted exactly over the element's 744 99-mers:

   | group | highest pooled count of one | the smallest `-d` that keeps them all | the record's `-d` | kmers lost |
   |---|---|---|---|---|
   | 30 complete genomes | 331 | 11.03, so 12 | 12 (the mean, 11.63) | none |
   | 15 complete + 15 contig | 162 | 5.40, so 6 | 5 (the mean, 4.93) | **18**, at element positions 433-450 |

   The draft group's run recovered the assay anyway, through the kmers below its ceiling — and the
   regions it reported begin at element position 451, exactly where the surviving kmers begin. Note also
   that the amplicon is 63 bases against `-k 99`, so no kmer lies *inside* it: the amplicon copy count is
   a proxy for the pooled count of a kmer spanning it, and on the whole panel it overstates (3.78 against
   a measured 3.11 for the highest-counting element kmer) while on the draft group the mean understates;
2. **use `-a unitig`.** With the same kmers, SKESA kept nothing of the element — its longest hit to it is
   19 bases — while `unitig` kept 40 contigs holding 100 bases or more of it, the longest three being 289,
   225 and 197. Set `--min-region` to the kmer size as well, since the pieces are short;
3. **check the assemblies.** If the inclusion group is draft-heavy, the repeat may be missing from most of
   it, and then `-p` has to come down — to 0.85 on this panel — or the panel has to be rebuilt from
   assemblies that resolved the repeat. Long reads resolve IS711; short-read contigs break inside it;
4. **read the bottom of the file**, and pass `design --max-regions 0`. With `-p` below 1 a repeat-targeted
   region is held by fewer genomes than a core one and is therefore ranked last;
5. **read `inclusion_copies_min`.** It is the number that says whether every inclusion genome really
   carries the amplicon. A `0` there is the honest version of a region the find stage called present.

## The runs

Every run of this record is in [`runs.tsv`](runs.tsv), and the per-genome census in
[`is711_census.tsv`](is711_census.tsv): one row per panel genome, with its assembly level, its whole
element copies, its amplicon copies, the longest piece of the element it holds and whether a piece of it
ends at a contig boundary.

```bash
# the panel, its blast database and the census
bash validation/brucella/run.sh /work 48 128        # panel.py and census.py are what is needed
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=M94960&rettype=fasta" \
    -o /work/is711.fasta
python validation/repeats/is711_census.py /work -t 32

# the two small groups, the 2 x 2, the draft group and the designs
bash validation/repeats/run_is711.sh /work /is711 16 32

# the whole panel, reusing the genus target's folders from the Brucella record
R=/work/runs/genus_as_labelled_p1.0
for opts in "-d 1 -p 1.0" "-d 6 -p 1.0" "-d 6 -p 0.85"; do
  primer-finder find -i $R/inclusion -e $R/exclusion -o full_${opts// /_} \
    -t 40 -m 128 -k 99 --min-region 99 -a unitig $opts
  python validation/brucella/check_assay.py full_${opts// /_} --target genus \
    --assays validation/repeats/is711_assay.tsv
done
```

## What this does and does not say

Every number here is in silico. The three oligos being where they are says the published target is
reachable on this panel at these settings; it says nothing about how the assay behaves in a reaction,
where specificity and sensitivity depend on conditions no genome file records. "Specific on this panel"
means no exclusion genome of this panel holds the amplicon, which is a statement about 315 assemblies.

What the record does establish:

- a real published repeat-targeted assay, on 1,691 real genomes, is recovered by primer-finder **only**
  with `-d` at the group's mean copy number and `-a unitig`, and each of those two alone is not enough;
- the copy numbers `design` reports for a real repeat are right: 4 to 40 per inclusion genome on the
  complete-genome group, against 4 to 40 counted by blast, and 0 as the minimum on a group where some
  assemblies have lost the element;
- the reason a repeat-targeted assay gives false negatives against public draft assemblies is in the
  assemblies: 557 of 873 contig-level assemblies of this panel hold no whole copy of IS711, and 556 of
  those 557 hold pieces of it ending at a contig boundary. No setting recovers sequence that is not in
  the file, which is why `-p` has to drop on the full panel and why the right fix for a repeat-targeted
  target is a panel of assemblies that resolved it.
