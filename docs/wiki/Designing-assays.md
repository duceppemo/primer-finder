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
melting temperature may fit there. **A forced request cannot produce a bad oligo** — the melting
temperature, length and GC limits are part of every request, and Primer3 returns nothing rather than
something outside them. The runs are taken longest-first, up to four per region.

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
They are weaker than they look, and are reported here only so that nobody mistakes them for more:

On 199 difference-based assays from *Xylella fastidiosa* subsp. *multiplex* regions, against 17 exclusion
genomes in PCR mode, the share that amplified only the inclusion group rose with the number of differences
under the primers (1: 20%, 2: 72%, 3: 91%, 4 or more: 100%), and assays with a run of two differences at a
3' end **and nothing else** were selective in none of the eight cases seen.

Two reasons not to read that as a result about PCR:

1. **It is close to circular.** in silico PCR decides whether a primer binds by counting mismatches against
   its own tolerance (`--mismatches`, 0 by default). More mismatches therefore means fewer reported hits
   almost by construction. The numbers describe the model's rule as much as the biology.
2. **It is one organism, one dataset, one set of parameters.** 25 genomes, 20 regions, one in silico tool.

So: the ranking is a sensible order to work through, supported by what others have published about
mismatches, and the in silico numbers are a consistency check on one dataset. They are not a measurement of
how these assays behave in a tube.

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
  they are not a universal law.
- **Only the first 50 regions are designed on by default** (`--max-regions`), because the regions are already
  ordered most-promising-first and in silico PCR of thousands of assays takes a while. `--max-regions 0`
  does all of them.
- **Nothing here checks the assay against anything but these genomes**: no secondary structure beyond what
  Primer3 scores, no cross-reaction with what is not in the two folders.
