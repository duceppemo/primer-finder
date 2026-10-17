# What "present in every inclusion genome" was checking, and what it now says

> **Erratum, same day, after a review of the measurements.** The cost of a coverage rule was computed
> wrongly below and the figure has been corrected in place, with the original kept in the open: this record
> first said such a rule "would drop 31.4% of what the full *Brucella* panel reports at `-p 0.85`". That
> 31.4% is the share of regions **some** genome holds only part of, which is not what applying the rule
> would do — a coverage rule changes the presence verdict, so a region falls only when fewer than the
> `-p` quorum of genomes hold enough of it. Applied that way the cost on that run is **6 regions of
> 7,014**, and the subsidiary claim that "even a rule asking for a quarter of the region drops one of the
> three" is **false**: at 0.25 all three survive. The two readings coincide at `-p 1.0`, so every other row
> of the table was unaffected. **The decision does not change, and the corrected numbers make the case
> more sharply**: the 0.9 rule drops only 6 of 7,014 regions, and all three regions carrying the published
> IS711 assay are among those 6. `partial_presence.py` now costs a rule the way the run applies presence
> and reports both readings separately.

2026-10-16, primer-finder 1.9.1, measured with the reporting change that follows (unreleased at the time
of the runs). The question came out of the IS711 record of the same day
([`2026-10-16_v1.9.1_is711/`](../2026-10-16_v1.9.1_is711/SUMMARY.md)), which found a case of it and set
out two ways to settle it. This record measures both and takes one.

## The rule, as it was

`blast.presence_in_genomes` called a region present in a genome when blast returned **any** hit for it
under `-evalue 1e-10`, with no requirement on how much of the region that hit covered. That verdict does
two things: it decides whether a region is reported at all, and — with `-p` below 1 — it is the count the
regions are ordered by, which is what `design --max-regions` then takes from the top.

The case that raised it: on a draft-heavy *Brucella* group a 338-base region was reported as present in
all thirty inclusion genomes while one of them, a contig-level *B. canis* assembly, held **41 identical
bases of it at the start of a contig** and not the published amplicon at all.

## What the rule lets through, on every answer key

`validation/partial_presence.py` re-runs the presence search on a finished run's reported regions and
records, for each, the fraction the stingiest genome holds. "Partly" below means under 0.9, counting only
the genomes that hold any of the region at all — the genomes that hold none are the `-p` story, and are
already in the region's header.

"Partly held" is descriptive, not a cost: it counts regions some genome holds under 0.9 of. What a rule
would actually drop is the column after it, applied as the run applies presence.

| run | `-p` | regions | partly held | a 0.9 rule drops | lowest | the published assay's region |
|---|---|---|---|---|---|---|
| *Xylella* fastidiosa, 6 complete | 1.0 | 782 | 2 (0.3%) | 2 | 0.579 | 1.000 |
| *Xylella* morus, 5 complete | 1.0 | 514 | 1 (0.2%) | 1 | 0.778 | 1.000 |
| *Xylella* multiplex, 8 complete | 1.0 | 1,223 | 1 (0.1%) | 1 | 0.843 | 1.000 |
| *Xylella* pauca, 4 complete | 1.0 | 3,471 | 5 (0.1%) | 5 | 0.463 | 1.000 |
| *Listeria* CC5, 22 complete | 1.0 | 84 | 0 | 0 | 1.000 | 1.000 |
| *Listeria* CC6, 36 complete | 1.0 | 67 | 0 | 0 | 1.000 | 1.000 |
| *Listeria* CC7, 47 complete | 1.0 | 85 | 0 | 0 | 1.000 | 1.000 |
| *Listeria* CC8, 52 complete | 1.0 | 83 | 0 | 0 | 1.000 | 1.000 |
| *Brucella* genus, 1,376 draft-heavy | 1.0 | 398 | 18 (4.5%) | 18 | 0.301 | not recovered here |
| IS711, 30 complete | 1.0 | 10,414 | 14 (0.1%) | 14 | 0.201 | 1.000, 1.000, 1.000 |
| IS711, 15 complete + 15 contig | 1.0 | 10,131 | 31 (0.3%) | 31 | 0.121 | **0.121** |
| IS711, whole panel | 0.85 | 7,014 | 2,204 (31.4%) | **6** | 0.021 | **0.112, 0.289, 0.310** |

Every run but the last is at `-p 1.0`, where a region has to be in every genome and the two columns
therefore coincide: one genome holding too little is one genome short of the quorum. Only the `-p 0.85`
run separates them, and it separates them by a factor of 370.

Every run of that table is in [`partial_presence.tsv`](partial_presence.tsv), with how many regions a
rule at 0.9, 0.75 or 0.5 would drop and what each published assay's region is held to. The per-region
figures stay beside each run, as `partial_presence.tsv` in its own folder.

Two things follow, and they point the same way.

**On complete-genome panels the question is moot.** Nothing is partly held on the *Listeria* panels, and
one to five regions per run on *Xylella*. All eight published assays recovered on those panels sit in
regions every inclusion genome holds in full, so a coverage rule would have changed nothing there — it
would also have gained nothing.

**On draft-heavy panels, and on repeats, a coverage rule takes the few regions that matter most.**
Costed properly — a genome holds a region when it aligns at least the cut-off, and the region survives
when the `-p` quorum of genomes do — a rule at 0.9 drops **6 of the 7,014** regions the full *Brucella*
panel reports at `-p 0.85`. **All three regions carrying the published IS711 assay are among those 6**:
against the 1,170 genomes that run requires, the genomes holding 0.9 of them number 1,115, 1,100 and
1,119. At 0.75 the same three go (1,140, 1,151, 1,126); at 0.5 two of the three; at 0.25 none of them, and
nothing else in the run is dropped either. On the draft-heavy 30-genome group, which ran at `-p 1.0`, the
one recovered region goes at **every** cut-off down to 0.25, because 29 of its 30 genomes reach them and
30 are required.

So the cost is tiny in count and total in effect: a coverage rule would leave 99.9% of a draft-heavy run's
regions untouched and remove precisely the repeat-targeted answer key. That is the shape of the hazard this
project has spent four releases on — a rule that looks harmless in aggregate and takes the hard case.

## The decision: report the fraction, keep the rule

Presence stays "any hit at all". The reasons, in order of weight:

1. a coverage rule costs the repeat recovery outright, as above — it is the same class of false negative
   this project has spent four releases removing, reintroduced one step earlier in the pipeline;
2. the question it would answer is already answered, better, one step later. `design` checks each amplicon
   it proposes against every inclusion genome and reports `inclusion_copies_min`; a 0 there is the precise
   version of "an inclusion genome does not carry this", about the amplicon that would actually be
   ordered rather than about the region that contains it;
3. a region whose shared part carries a usable assay is worth reporting. Dropping it because its edges
   are patchy throws away the part the genomes agree on.

What was wrong was not the rule but the silence. So:

- **`presence_in_genomes` returns the covered fraction** rather than a yes or no, merging overlapping
  pieces (a repeat is hit once per copy over the same columns, and adding those up would read as a genome
  holding several times the region);
- **the presence table carries it**: `4_blast/inclusion_blast_hits.tsv` holds the fraction to three
  decimals where it held 1 or 0. Presence is still anything above 0, so a reader comparing with an older
  run reads the same rows;
- **the presence search asks for every subject** (`-max_target_seqs 500`) instead of the best one. With
  one subject a region split across two contigs of a draft assembly is counted from the longer piece
  only, which understates the fraction exactly where assemblies break — which is where repeats are;
- **the run says so**: after the "present in all inclusion genomes" line, a run that kept a region some
  genome holds only part of reports how many and the worst figure, and points at the table.

Nothing is kept, dropped or ordered differently. The same regions come out in the same order.

**One cell of one run does not reproduce, and it is not this change.** Of the 398 regions of the
*Brucella* genus run, 397 re-measure as the run recorded them. For `Contig_321_54.9391` the run's table
says one genome holds it and a re-blast today finds that genome's best hit at **1.81e-10**, just the wrong
side of the 1e-10 cut-off, at both subject caps — so under the re-measurement that region would not be
reported at `-p 1.0`. A presence verdict resting on 33 of 135 bases at the edge of the cut-off is not
stable across a blast version or a database rebuild, which is the record's own point in miniature: the
table now says the region was held to 24% by that genome, where before it said `1`.

## The ordering, which was the other half of the question

Whether the `-p` ordering should weigh coverage rather than count genomes was the one thing the IS711
record left open, since `design --max-regions 50` takes its regions from the top. Measured, it should not:
on the two panels where recoveries happen nothing is partly held, so a coverage-weighted order would be
identical; on the panels where it is not, the repeat-targeted regions are already last because fewer
genomes hold them at all, and weighing coverage would put them no higher. The honest answer for a repeat
hunt remains `design --max-regions 0`, which the IS711 record gives.

## What it costs, and what asking for every subject buys

Asking blast for 500 subjects instead of the best one is the only part of this that is not free. Measured
on the *Brucella* genus run's own inputs — its 721 candidate regions against its 1,376 inclusion genomes,
40 threads, the two caps back to back on the same machine:

| `-max_target_seqs` | the step | regions in every genome | regions a genome holds only part of |
|---|---|---|---|
| 1 | 25.0 s | 426 | 45 |
| 500 | 30.7 s | 426 | **31** |

Two things in that table. The cost is **23%** of one step, which on that run is seconds against minutes of
wall time. And the fourth column is the reason to pay it: with one subject, 45 regions looked partly held
and 14 of those were an artifact of counting only the longest piece — a region that a draft assembly had
split between two of its contigs was being credited with one of them. The third column is the same either
way, which is the other thing worth stating: the regions reported do not change.

## What this does and does not say

These are fractions of a region that blast aligns, over the genomes a run was given. A high fraction is
not a statement that an assay designed there will amplify, and a low one is not a statement that it will
not: what an amplicon does in a reaction depends on conditions no genome file records. The figure's only
claim is about sequence in assemblies.

It also says nothing about *why* a genome holds part of a region. A broken assembly, a genuine deletion
and a diverged copy all look the same here. `is711_census.py` in the IS711 record separates the first from
the others for one element by asking whether the pieces end at contig boundaries; nothing in the pipeline
does that in general.
