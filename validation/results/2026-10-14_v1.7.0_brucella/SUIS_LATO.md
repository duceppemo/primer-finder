# *B. suis* with *B. canis* in the inclusion group

2026-10-14, primer-finder 1.7.0. [`ALTERNATIVES.md`](ALTERNATIVES.md) found no *B. suis* region at any
setting and said the target needed its groups redrawn, because *B. canis* is a lineage **inside**
*B. suis*: no stretch is present in every *B. suis* and absent from every *B. canis*. So a `suis_lato`
target was added, with both species in the inclusion group, and the question asked again.

**It still finds nothing, and now the reason is measured rather than inferred.** There are three separate
obstacles, and the deepest is the one that also stops the published species assays.

## First: this asks a different question from the published BS assay

The published BS carries all three oligos exactly in **115 of the 194** *B. suis* and *B. canis* genomes —
120 of 120 *B. suis*, but only **1 of 74** *B. canis* (`bs_lato_scorecard.tsv`):

| | genomes | BS covers |
|---|---|---|
| `suis` (the paper's target) | 120 | 114 = 95.00% |
| `suis_lato` (*suis* + *canis*) | 194 | 115 = **59.28%** |

So BS is a *B. suis* sensu stricto assay. Nothing from `suis_lato` would replace it: an assay from this
target detects *B. canis* as well, which is a different diagnostic claim.

## Second: at the default kmer size there is no specific kmer at all

| panel | `-p` | required | result |
|---|---|---|---|
| corrected | 1.0 | 193/193 | no inclusion-specific kmer |
| corrected | 0.99 | 192 | no inclusion-specific kmer |
| corrected | 0.95 | 184 | no inclusion-specific kmer |
| as labelled | 1.0 / 0.99 / 0.95 | 194 / 193 / 185 | no inclusion-specific kmer |

Every 99-mer shared by the *B. suis* and *B. canis* genomes is also somewhere in the exclusion group —
*B. abortus*, *B. melitensis*, *B. ceti*, *B. ovis*, *B. neotomae*, *B. microti* and the 315 relatives.

A kmer survives the subtraction only if it spans a difference, so a **longer** kmer spans one more often
and is the lever to try, not a shorter one (`suis_lato_ksweep.tsv`):

| `-k` | `-p 1.0` | `-p 0.95` |
|---|---|---|
| 201 | no specific kmer | **15 specific kmers**, assembler refused them |
| 151 | no specific kmer | no specific kmer |
| 127 | no specific kmer | no specific kmer |
| 71 | no specific kmer | no specific kmer |

## Third: the 15 kmers are one locus, and the find rule rejects it

The 15 specific 201-mers overlap heavily — 55 of their ordered pairs share 150 bases or more — and merge
into **one 215 bp stretch** (`suis_lato_stretches.fasta`). SKESA will not assemble 15 reads, which is why
the run ends at "could not assemble", but the stretch is real. It is present at 99% identity or better
over at least 90% of its length in:

| | | | |
|---|---|---|---|
| ***B. suis* 120/120** | ***B. canis* 74/74** | *B. abortus* 392/398 | *B. melitensis* 553/734 |
| *B. ceti* 13/21 | *B. ovis* 0/18 | *B. neotomae* 5/6 | *B. microti* 5/5 |
| non-classical relatives **0/315** | | | |

So the locus is not specific by absence — it is present in most of the exclusion group too. What makes its
201-mers specific is a handful of bases. Aligned to the nearest homologue in each:

| against | differences in the 215 bp | spacing | the find rule |
|---|---|---|---|
| *B. abortus* 2308 | **one**, at position 179 | — | cannot keep it: one difference |
| *B. melitensis* 16M | **two**, at 145 and 179 | **34 bases apart** | cannot keep it: needs ≤ 25 |

```
position 145:  ...TGAAGGCCGCGCCGCACTGC[G]GGAAAGGGAGTGAGGATGGG...
position 179:  ...GGATGGGACAGAAACCACTT[T]GTTTCACAAAGGCAGATATT...
```

Nine bases. Had the two differences fallen 25 rather than 34 apart, the rule would have kept the region
and an assay could have been designed on it.

## What this adds up to

Redrawing the groups was the right move — it turned "nothing at all" into one identified candidate locus —
but it does not get a *B. suis* sensu lato assay out of primer-finder 1.7.0. The three obstacles in order:

1. **the default `-k 99` is too short** for this group: a 99-mer spanning the locus is not unique to it.
   `-k 201` with `-p 0.95` is what surfaces the signal, and neither of those is a setting anyone would
   reach for without this record;
2. **the assembler will not build a contig from 15 kmers**, so the run ends in an error even though the
   kmer step succeeded. The same limitation as the `-k 31` failure in [`SUMMARY.md`](SUMMARY.md), from the
   other end of the kmer range. **Since fixed**: `-a unitig` builds the 215-base stretch from those 15
   kmers and the run reaches step 3, where it is turned down for the reason below. `unitig_demo.tsv`
   records both assemblers on this run and on two others;
3. **the find rule rejects the locus anyway**: one difference from *B. abortus*, two 34 bases apart from
   *B. melitensis*. This is the same wall as the three published species assays, and it is not about
   settings.

The locus is in `suis_lato_stretches.fasta`, and the run now reports it on its own: `-a unitig` builds
the stretch from those 15 kmers and the last-resort pass keeps it on its single difference, as

```
>unitig_1 [178] inclusion=193/193 last_resort        215 bases
```

Position 178 counted from 0 is the base measured by hand above. Two of the three obstacles are therefore
gone, and the third stands: the region is reported **only** because the last-resort rule accepted one
difference, which is what the mark says, and an assay on it would rest on that single base. See
[`../2026-10-14_v1.7.0_unitig/SUMMARY.md`](../2026-10-14_v1.7.0_unitig/SUMMARY.md). An assay there
would discriminate *B. suis* and *B. canis* from *B. abortus* on **a single base at position 179** — the
same class of assay as the published BA, and the same reason primer-finder declines to propose it.

## The wider conclusion for this panel

Across every target tried, the classical *Brucella* clade is uniform enough that species-level
discrimination rests on **single scattered SNPs**, not on blocks of unique sequence. That is why the
published species assays exist and why primer-finder does not propose them: its rule asks for two
differences close enough to sit under one primer, which is a deliberate judgement about what makes an
assay robust, and this clade rarely offers it. The genus level is the opposite case — whole stretches are
present in every classical *Brucella* and absent from every relative, and there primer-finder does well
(see [`WHOLE_CLADE.md`](WHOLE_CLADE.md)).

## Files

| File | What it holds |
|---|---|
| `suis_lato_runs.tsv` | the six runs at the default kmer size, both panels, three `-p` values |
| `suis_lato_ksweep.tsv` | the eight runs sweeping `-k` to 201, with the specific-kmer count of each |
| `suis_lato_stretches.fasta` | the 215 bp locus the 15 specific 201-mers merge into |
| `bs_lato_scorecard.tsv` | the published BS assay against *B. suis* alone and against *suis* + *canis* |
