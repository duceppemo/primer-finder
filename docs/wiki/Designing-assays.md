# Designing assays

`primer-finder find` reports **regions**. `primer-finder design` turns them into **assays**: it runs Primer3
on the regions of a finished run, works out what would make each assay selective, and writes the primer
files that [insilicoPCR](https://github.com/duceppemo/insilicoPCR) reads — optionally running it as well, to
check the assays against both groups of genomes.

```bash
primer-finder design results/ -o assays/ -t 16
primer-finder design results/ -o assays/ -t 16 --insilico-pcr /path/to/insilicoPCR-linux-x64
```

The genomes come from what the run recorded in `run_info.json`, so the two folders need not be given again;
`-i` and `-e` override that if they have moved.

## What makes an assay selective

A region is inclusion-specific as a whole: every kmer of it is absent from the exclusion genomes. **An assay
designed inside it is not specific by itself.** Only two things make it so, and the command tells them apart:

1. **Absence.** The amplicon does not exist in any exclusion genome, so there is nothing to amplify. This is
   the strongest case and it is how the published assays this tool is
   [validated](Validation) against work — their oligos cover no differing base at all; the stretch they sit
   in is simply missing from the other subspecies.
2. **Differences.** The amplicon does exist there, but an oligo sits on bases that differ — the lower-case
   bases of `final_kmers.fasta`. The assay then depends on a mismatch stopping the reaction.

An assay that is neither — amplicon present, no oligo on a difference — would amplify both groups. Those are
listed in `assays.tsv` with `specific_by=nothing` and go no further.

## How Primer3 is asked

Three kinds of request go in for each region, and the results are pooled and de-duplicated:

| Request | What it asks for | Why |
|---|---|---|
| plain | the best assay anywhere in the region | the region may be absent from the exclusion genomes, in which case any assay in it works |
| `SEQUENCE_FORCE_LEFT_END` / `SEQUENCE_FORCE_RIGHT_END` | a primer whose 3' end sits on the last (or first) base of a run of differences | allele-specific PCR: a mismatch at the 3' end hinders extension |
| `SEQUENCE_INTERNAL_OVERLAP_JUNCTION_LIST` | a probe straddling the middle of a run of differences | a probe that cannot bind the exclusion template |

The forced requests often come back empty, which is not an error: no oligo of the required length and
melting temperature may fit there. **A forced request cannot produce a bad oligo** — the constraints below
are part of every request, forced or not, and Primer3 returns nothing rather than something outside them.
The runs are taken longest-first, up to four per region.

## What an oligo has to satisfy

These are **hard constraints**, not score terms: Primer3 returns nothing rather than an oligo outside them,
so a doomed assay is never proposed and never reaches the ranking. The ranking below only ever orders assays
that already satisfy all of this.

| | Primer | Probe |
|---|---|---|
| Length | 18–25 nt, 20 optimal | 18–27 nt, 22 optimal |
| Melting temperature | 58–63 °C, 60 optimal | 62–72 °C, 68 optimal |
| GC | 30–70% | 30–80% |
| Ambiguous bases | none | none |
| G or C at the 3' end | `--gc-clamp`, 1 by default | — |
| Hairpin | melts below `--max-hairpin-tm`, 47 °C by default | same |
| Self-dimer, whole oligo and 3' end | below `--max-dimer-tm`, 47 °C by default | same |
| Dimer with the other primer, whole oligo and 3' end | below `--max-dimer-tm` | — |

### The 3' end: `--gc-clamp`

A G or a C at the 3'-most base holds the primer down where extension starts, so one is asked for by default
(`PRIMER_GC_CLAMP=1`); `--gc-clamp 2` asks for two, `--gc-clamp 0` for none.

**The clamp is not applied to a primer whose 3' end is pinned on a differing base.** Those are the
`SEQUENCE_FORCE_LEFT_END` / `SEQUENCE_FORCE_RIGHT_END` requests, where the 3' base is whatever the genomes
made it — asking for a G or a C there is asking for a difference that may not exist. On an A/T-rich target
the two requirements are almost never satisfiable at once: with the clamp left on, seven of eight forced-end
requests on the *Xylella* set returned nothing. So the clamp applies to every request that lets Primer3
choose the end, and is dropped for the ones that pin it.

### Hairpins and dimers

An oligo that folds back on itself, pairs with a copy of itself, or pairs with its partner is spent before it
ever reaches the template — and a 3'-end dimer is worse than an internal one, because the polymerase can
extend it. Primer3 models all of these thermodynamically (`PRIMER_MAX_HAIRPIN_TH`, `PRIMER_MAX_SELF_ANY_TH`
and `_SELF_END_TH`, `PRIMER_PAIR_MAX_COMPL_ANY_TH` and `_COMPL_END_TH`, and the `PRIMER_INTERNAL_*`
equivalents for the probe), and the design step sets every one of them, for the probe as well as the primers.
Nothing extra has to be installed: this is Primer3's own thermodynamic alignment, so there is one melting
temperature model for the oligo and for the structures it might form.

The default limit of 47 °C is Primer3's own, which is roughly 10 °C below the annealing temperature these
oligos are designed for: a structure that has melted by the time the reaction anneals does not compete with
the template. Lower it to be stricter (`--max-hairpin-tm 40`), raise it on a target where nothing else will
fit — at the cost of candidates that may fold.

What Primer3 predicted is **also reported**, so a reaction run under other conditions can be judged:
`pair_dimer_tm` and `pair_dimer_end_tm` per assay, and `forward_hairpin_tm`, `forward_self_dimer_tm` and the
`reverse_` and `probe_` equivalents per oligo. A `0.0` means no structure was predicted, not that none was
looked for.

### The reaction those temperatures are predicted for

A melting temperature is only meaningful for a given reaction: salt stabilises a duplex, magnesium more so
per mole, dNTPs chelate magnesium away, and oligo concentration shifts the equilibrium. Primer3 is therefore
told what reaction to predict for, and the defaults are an ordinary TaqMan qPCR:

| Option | Default | Primer3 tag |
|---|---|---|
| `--monovalent` | 50 mM | `PRIMER_SALT_MONOVALENT` |
| `--divalent` | 3 mM Mg²⁺ | `PRIMER_SALT_DIVALENT` |
| `--dntp` | 0.8 mM total (0.2 mM each) | `PRIMER_DNTP_CONC` |
| `--primer-conc` | 250 nM | `PRIMER_DNA_CONC` |
| `--probe-conc` | 200 nM | `PRIMER_INTERNAL_DNA_CONC` |

The probe gets its own concentration, since it is normally used below the primers. These are not cosmetic:
they change both the temperatures reported and which oligos come back at all — on one test region the same
probe was 57.9 °C under Primer3's bare defaults and 60.0 °C under these, which is the difference between
failing and passing a 58 °C floor.

**Set them to your own master mix** if it differs; the values used are recorded in `design_info.json` under
`parameters.conditions`, so a table of assays can always be traced back to the reaction it was designed for.
Changing them does not change the oligos' behaviour at the bench, only which ones this step proposes.

## The scoring scheme

The order is a **heuristic for which assays to look at first**, not a prediction of what will work at the
bench. Nothing here has been tested in a laboratory. Read it as "these are the ones worth trying", and let
in silico PCR and then the bench decide.

Assays are ordered by this key, each step breaking the ties of the one before:

1. **Absence beats everything.** An assay whose amplicon no exclusion genome holds has nothing to amplify
   there. This is the one step that does not depend on how a reaction behaves.
2. **How many differences the primers cover, in total.**
3. **What those differences replace**, weighted: a position where the exclusion genomes have a G or a C
   counts 1, an A or a T counts 0.5, and one whose base is not known counts 0.75. G and C pair with three
   hydrogen bonds against two, so a mismatch there costs more to form. This only ever separates assays that
   cover the same number of differences.
4. **The longest run of differences ending at a primer's 3' end.**
5. **Differences within five bases of a 3' end**, run or not.
6. **Primer3's pair penalty, by band** (<1, <2, <4, worse), so that a probe covering more differences never
   wins over an assay that is clearly better made.
7. **Copies of the amplicon in the inclusion genomes**, since a repeated target usually gives a better limit
   of detection. See [below](#repeated-targets-and--d).
8. **Differences under the probe**, then the exact penalty, then the names, so that a run is reproducible.

### Why that order, and how much to trust it

Steps 2 to 4 follow common practice in allele-specific design rather than any measurement made here:

- a single mismatch, even at the 3' end, is often not enough on its own. Allele-specific PCR is known for
  "low discriminating power", and a common remedy is to introduce a **second, artificial mismatch** in the
  primer — the basis of double-mismatch allele-specific qPCR
  ([Lefever et al. 2019](https://doi.org/10.1038/s41598-019-38581-z)). More differences under a primer is
  the same idea, arrived at without having to engineer one;
- *where* a mismatch sits, and **which bases are involved**, changes its effect
  ([Sharma et al. 2022](https://doi.org/10.1016/j.jmoldx.2022.08.005) summarise this for PCR while
  characterising it for RPA). Hence two tie-breaks: the 3' end is the position that matters most for
  extension, and a difference replacing a G or a C is worth more than one replacing an A or a T, since
  G:C holds with three hydrogen bonds and A:T with two.

The bases are read from the exclusion genomes' own alignments of each amplicon, so what is weighed is what
an oligo would actually have to mismatch, genome by genome, rather than what the region happens to carry.
`assays.tsv` reports the weight (`primer_variant_weight`) and the count of G/C positions per primer
(`forward_strong_variants`, `reverse_strong_variants`), so the order can be checked.

**How much a given mismatch costs depends on the assay, not only on the sequence**: annealing temperature,
polymerase, magnesium, cycling and template concentration all change it. A design that discriminates under
one set of conditions may not under another, and optimisation at the bench can recover an assay this
ranking puts low — or ruin one it puts high.

### What the in silico PCR numbers do and do not say

Running the designed assays through in silico PCR gives numbers that look like evidence for the order above.
They are weaker than they look, and are reported here only so that nobody mistakes them for more. On 287
difference-based assays from *Xylella fastidiosa* subsp. *multiplex* regions, against 17 exclusion genomes,
the share that amplified only the inclusion group rises with the number of differences under the primers:

| differences under the primers | 0 (probe only) | 1 | 2 | 3 | 4 or more |
|---|---|---|---|---|---|
| qPCR mode, selective | 43% | 36% | 60% | 64% | 94% |
| PCR mode, selective | 0% | 7% | 39% | 40% | 91% |

Three reasons not to read that as a result about PCR:

1. **It is close to circular.** in silico PCR decides whether a primer binds by counting mismatches against
   its own tolerance (`--mismatches`). More differences therefore means fewer reported hits almost by
   construction. The numbers describe the model's rule as much as the biology.
2. **The model cannot see a difference at the very 3' end** — the position the ranking cares most about.
   See the three zones [below](#what---mismatches-actually-governs). Assays whose only differences are one
   or two bases at a 3' end are reported selective in none of the 65 cases seen, at **any** tolerance,
   because those bases are not counted. That is a property of the alignment, and it is the one number from
   an earlier version of this page that has had to be withdrawn as evidence: it says nothing either way
   about whether such an assay would discriminate at the bench.
3. **It is one organism, one dataset, one set of parameters.** 25 genomes, 20 regions, one in silico tool.

So: the ranking is a sensible order to work through, supported by what others have published about
mismatches, and the in silico numbers are a consistency check on one dataset. They are not a measurement of
how these assays behave in a tube. The full sweep is in
[`validation/results/2026-10-08_v1.3.0/SUMMARY.md`](https://github.com/duceppemo/primer-finder/blob/master/validation/results/2026-10-08_v1.3.0/SUMMARY.md).

One thing in that record is **not** circular, and is worth acting on: of the 178 assays that are specific by
absence — no designed mismatch anywhere, their amplicon simply missing from the exclusion genomes — 30
amplify *something else* in an exclusion genome at the default tolerance, and 50 do at `--mismatches 2`.
Off-target amplification is a real specificity risk that has nothing to do with the ranking, and it is the
main reason to run the check at all.

## Checking the assays with in silico PCR

insilicoPCR is a separate program and stays one. `primer-finder design` writes what it reads, and will run
it for you:

```bash
primer-finder design results/ -o assays/ --insilico-pcr /path/to/insilicoPCR-linux-x64
```

`--insilico-pcr` takes the folder of an extracted portable release (which brings its own Java), its jar, or
a launcher script. Without it, the same thing can be run later with the `run_insilico_pcr.sh` the command
writes.

**Two files, two runs.** insilicoPCR cannot report both kinds of assay at once: one probe in a primer file
puts the whole report in qPCR mode, where an assay without a probe is never positive. So every usable assay
is written twice — `assays_qpcr.fasta` with its probe, `assays_pcr.fasta` with the primers alone — and each
file is run against both groups. That gives two answers per assay, and they are worth having separately:

- **`pcr_selective`** says whether the primer pair discriminates on its own;
- **`qpcr_selective`** says whether the whole assay reports positive, with the probe having to bind too.

An assay is `selective` in a mode when it amplifies **every** inclusion genome and **no** exclusion genome.

### What `--mismatches` actually governs

insilicoPCR takes a mismatch tolerance, and primer-finder passes `--mismatches 1` by default. What that
governs was measured rather than assumed — a primer pair in a synthetic template, the template mutated one
base at a time at a known distance from the primer's 3' end, run at tolerances from 0 to 10
([`mismatch_zones.py`](https://github.com/duceppemo/primer-finder/blob/master/validation/results/2026-10-08_v1.3.0/mismatch_zones.py)).
A primer turns out to have three zones, and the tolerance only controls one:

| Where the mismatch is | What insilicoPCR does |
|---|---|
| the last 2 bases | **free.** Not counted as a mismatch at any tolerance, including 0. blast trims an unmatched base off the end of its alignment, and insilicoPCR reports the trim as a negative `EndMismatch` offset and calls the primer bound. |
| 3 or 4 bases from the 3' end | **never amplifies**, at any tolerance — 10 was tested. A trim that long is rejected, and keeping the mismatch scores worse for blast than trimming. |
| 5 or more bases in | **what `--mismatches` decides**, counted per primer: at `--mismatches 1` each primer may carry one, so an assay may carry two. |

This is a property of the alignment, not of PCR, and there is no setting that changes it. Two consequences:
a difference at the very 3' end — what the forced-end requests aim for — cannot be rewarded by this check,
and a difference 3 or 4 bases in is treated as fatal when it may not be.

**Why 1 and not 0, 2 or 3.** `--mismatches 0` assumes any mismatch five or more bases from the 3' end stops
a primer, which is the optimistic end of what is known: a single internal mismatch frequently does not stop
amplification, which is why double-mismatch allele-specific designs exist
([Lefever et al. 2019](https://doi.org/10.1038/s41598-019-38581-z)). It also hid 12 of the 30 off-target
amplifications above. `--mismatches 2` assumes two mismatches *per primer* still bind, which is more than
that literature supports — two deliberate mismatches are what an allele-specific design uses to
discriminate — so it belongs in a second, stricter pass rather than in the default. `--mismatches 3`
measured nothing that 2 did not: on the *Xylella* set the two give identical verdicts in PCR mode, one assay
apart in qPCR mode.

```bash
primer-finder design results/ -o assays/ --insilico-pcr <dir>                   # the default, -m 1
primer-finder design results/ -o assays/ --insilico-pcr <dir> --mismatches 2    # the stricter check
```

An assay that still passes at `--mismatches 2` rests on something the model cannot explain away. Raising the
tolerance only ever finds more exclusion amplification, never less, so it can only take assays away.

### When a `no` is the model's blind spot

Because the last two bases are free, an assay can be reported non-selective for a reason that is about the
alignment rather than about the assay. `assays.tsv` therefore carries **`qpcr_exclusion_terminal_only`** and
**`pcr_exclusion_terminal_only`**: of the exclusion genomes an assay amplified, how many did so *only* where
a difference sits in the last two bases of a primer. A genome counts only when every amplicon reported for
it needed the trim — one clean amplicon is a real amplification.

When that number equals `*_exclusion_amplified`, **every** exclusion amplification behind the `no` is one
insilicoPCR could not have refused, and the design is one in silico PCR cannot judge in either direction.
The command says so as it runs, and on the *Xylella* set it is most of them:

| PCR mode, 465 usable assays | |
|---|---|
| selective | 220 |
| `no`, but every exclusion amplification rests on an uncounted 3'-end difference | **155** |
| `no`, partly | 7 |
| `no`, on amplifications the model can defend | 83 |

Of those 155, 127 rest on a difference the ranking deliberately put at a 3' end (96 of one base, 22 of two,
3 of three) and 28 are specific by absence — an off-target amplicon somewhere else in an exclusion genome
that itself only binds through an ignored terminal mismatch, which is a weak off-target rather than a clean
one. Either way: the column tells you which `no` to argue with, and 83 rather than 245 is the number of
assays this check actually rejects on its own evidence.

The inclusion side needs no such column. Of 3,720 assay-and-genome amplifications in the inclusion group,
**none** depended on a trimmed primer end: every inclusion genome that counted had at least one clean
amplicon, which is what should happen when the assays were designed on sequence all of them share.

### When the regions were found with `-p` below 1

`-p/--min-inclusion` lets `find` keep a region that some inclusion genomes lack. An assay designed on such a
region cannot amplify genomes that do not hold it, so judging it against all of them would be unfair. The
design step reads the `-p` of the run out of its `run_info.json` and applies the same bar:

| `qpcr_selective` / `pcr_selective` | What it means |
|---|---|
| `yes` | amplifies every inclusion genome and no exclusion genome |
| `partial (88%)` | amplifies 88% of the inclusion genomes — at or above the run's `-p` — and no exclusion genome |
| `no` | below the threshold, or it amplifies an exclusion genome |

The percentage is also in `qpcr_inclusion_percent` and `pcr_inclusion_percent`, and is rounded down so that
it never overstates the coverage. Amplifying an exclusion genome is never excused, whatever `-p` was. Assays
that amplify every inclusion genome are listed before those that merely meet the threshold.

## Repeated targets and `-d`

A target present several times in each genome gives more template per cell, which usually improves the limit
of detection. `assays.tsv` reports `inclusion_copies_min` and `inclusion_copies_max`, counted by blasting
each amplicon against the inclusion genomes, and the ranking prefers more copies.

**The default `-d 1` of the `find` command discards repeated regions before any of this.** `-d` caps how
often a kmer may occur per inclusion genome, and KMC counts over the whole group: with `N` genomes and
`-d 1`, a kmer occurring twice in each genome has a count of `2N`, above the `-cx N` limit, and is dropped.
Verified with two genomes holding a region twice: `-d 1` returns only the single-copy kmers, `-d 2` returns
both. So if a multi-copy target is what you want:

```bash
primer-finder -i inclusion/ -e exclusion/ -o results/ -d 2   # or more
primer-finder design results/ -o assays/
```

and look at `inclusion_copies_min` in the table. The trade-off is that `-d` above 1 also lets through kmers
that occur several times inside a single genome, which is why it is not the default.

## What comes out

| File | What it holds |
|---|---|
| `assays.tsv` | every assay, the usable ones first, with the numbers behind the order (see [Outputs](Outputs)) |
| `assays_qpcr.fasta` | the assays that have a probe, for insilicoPCR in qPCR mode |
| `assays_pcr.fasta` | every assay, primers only, for insilicoPCR in standard PCR mode |
| `run_insilico_pcr.sh` | the in silico PCR of both files against both groups, to run later |
| `design_info.json` | the parameters, the programs and their versions, and the counts |
| `primer3/` | exactly what was sent to Primer3 and what came back |
| `amplicons/` | the blast databases and hits of the amplicon checks |
| `insilico_pcr/` | what insilicoPCR wrote, when it was run |

## Limits

- **The ranking is a prediction, not a result.** It orders assays by what makes them likely to work; in
  silico PCR is the check, and the wet lab is the answer.
- **The numbers behind the ranking come from one dataset.** They are from the *Xylella* validation set, with
  one organism, 25 genomes and 20 regions. They are consistent with what is known about 3' mismatches, but
  they are not a universal law — and in silico PCR cannot see a difference in the last two bases of a
  primer at all, so it can never confirm the step of the ranking that cares most about it.
- **Only the first 50 regions are designed on by default** (`--max-regions`), because the regions are already
  ordered most-promising-first and in silico PCR of thousands of assays takes a while. `--max-regions 0`
  does all of them.
- **The structure predictions are predictions.** Primer3's thermodynamic model is the same one behind its
  melting temperatures; it is good enough to throw out the obvious failures, not a substitute for looking at
  a candidate before ordering it.
- **Nothing here checks the assay against anything but these genomes**: no cross-reaction with what is not in
  the two folders, and no check against a transcriptome or a host.
