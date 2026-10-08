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

Assays are ordered by this key, each step breaking the ties of the one before:

1. **Absence beats everything.** An assay whose amplicon no exclusion genome holds cannot amplify the wrong
   group, and no number of mismatches is as good as nothing to amplify.
2. **How many differences the primers cover, in total.** This is the strongest predictor of whether an assay
   really is selective, which is not the obvious answer — see the numbers below.
3. **The longest run of differences ending at a primer's 3' end.** At equal totals this is the better place
   for them: a mismatch at the 3' end hinders extension more than one in the middle.
4. **Differences within five bases of a 3' end**, run or not.
5. **Primer3's pair penalty, by band** (<1, <2, <4, worse). Chemistry comes before the remaining signals, so
   a probe covering more differences never wins over an assay that is clearly better made.
6. **Copies of the amplicon in the inclusion genomes.** A target present several times per genome usually
   gives a better limit of detection. See [below](#repeated-targets-and--d) — this is only ever above 1 for
   a run made with `-d 2` or more.
7. **Differences under the probe**, then the exact penalty, then the region and assay names so that a run is
   reproducible.

### Why total differences outrank a run at the 3' end

A single mismatch at the 3' end does not always stop amplification. That is visible in the numbers, which
come from 199 difference-based assays designed on *Xylella fastidiosa* subsp. *multiplex* regions and run
through in silico PCR against 17 exclusion genomes (PCR mode, primers only):

| Differences under the primers | 1 | 2 | 3 | 4 or more |
|---|---|---|---|---|
| Amplified only the inclusion group | 20% | 72% | 91% | 100% |

| Longest run at a 3' end | 0 | 1 | 2 | 3 or more |
|---|---|---|---|---|
| Amplified only the inclusion group | 96% | 41% | 65% | 100% |

The run on its own is not the whole story: assays with a run of two and **only** two differences in total
were selective in **none** of the eight cases seen, while assays with three or more differences were
selective in 55 of 57 whatever their run.

The reason is in how the differences are found. A lower-case base in `final_kmers.fasta` differs from **at
least 90%** of the exclusion genomes that hold the region, not from every one of them (see
[Methods](Methods)). One or two such bases can still match a particular genome; four rarely do. More
differences therefore means the assay holds up across the whole exclusion group, which is what selectivity
means here.

This is why the order is: total first, then where they sit.

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
On the *multiplex* regions, 254 of 315 assays were selective in qPCR mode and 224 in PCR mode: the 30 that
differ are the ones the probe rescues.

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
