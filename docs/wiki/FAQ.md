# FAQ

## "No kmer at all is shared by the inclusion genomes"

The inclusion group shared nothing before any subtraction, so the exclusion group is not the reason: either
the inclusion genomes are too diverse for this kmer size, or one of them does not belong in the group. A
smaller `-k`, or `-p` below 1 so that a kmer need not be in every one of them, gives it room; remedies 1, 2,
5 and 6 below apply, and 3 is worth trying here in particular.

## "No inclusion-specific kmer was found"

Every kmer the inclusion genomes **did** share is also somewhere in the exclusion group — a different
problem from the one above, and the run says which of the two it hit. Things to try, roughly in order:

1. **Check the groups.** One exclusion genome that actually belongs to the inclusion group is enough to
   remove every candidate. The reverse is just as bad.
2. **Check the assemblies.** A kmer must be present without a single mismatch in *all* inclusion genomes, so
   a draft assembly with a gap in the right place, a contaminated one, or one of much lower quality than the
   others silently removes candidates. Compare the genomes first
   ([genome_comparator](https://github.com/duceppemo/genome_comparator)).
3. **Lower `-p`** (`-p 0.9`): accept the regions that most, rather than all, of the inclusion genomes share.
   The assay will not amplify the genomes that lack the region, and each record says how many hold it, so
   this is also the quickest way to find out whether one bad assembly is behind an empty run.
4. **Reduce the inclusion group** to the genomes you are most confident about, if a few of them are of
   poorer quality.
5. **Use shorter kmers** (`-k 51`, `-k 31`) — but only when the inclusion group is too diverse to share
   anything, which the "no kmer at all" message above says outright. A shorter kmer is less likely to be
   broken by a sequencing error or a small variant, but it also spans a difference less often, so it finds
   **fewer** specific regions, not more: 70 over the nine *Listeria* runs at `-k 31` against 306 at 99, two
   of them empty. See [Methods](Methods#1-the-kmers-only-the-inclusion-group-has-kmc).
6. **Raise `-d`** (`-d 2`, `-d 3`): a kmer that occurs twice in some inclusion genome is kept. This helps
   only where the region of interest is repeated; it cannot make a kmer the exclusion group holds
   specific.

## "Nothing usable came out of skesa"

Not an error. SKESA and SPAdes are built for sequencing reads and can get nothing at all out of a kmer
set — most often at a kmer size away from the default, where the specific stretches are short and few.
The run assembles the kmers directly instead, as `-a unitig` does, says so in a warning, and records
`assembler_used` in `run_info.json`. Nothing else about the run changes.

It fires when the assembler exits with an error, or returns less than 0.1 assembled bases per kmer
against 0.7 to 1.4 for a run that worked. It cannot be turned off, and it applies to SPAdes exactly as it does to SKESA. The numbers behind the
threshold are in the [validation record](https://github.com/duceppemo/primer-finder/blob/main/validation/results/2026-10-14_v1.7.0_unitig/SUMMARY.md).

## "No contig differs enough from the exclusion genome"

Inclusion-specific kmers were found and assembled, but none of the contigs has two differences within 25
bases of each other when mapped to the exclusion genome. The region is specific, but by a single mismatch
here and there, which does not make a selective assay.

**Since the release after 1.7.0** (the Unreleased section of the
[changelog](https://github.com/duceppemo/primer-finder/blob/main/CHANGELOG.md)) **a run does not stop here.** It tries again accepting a single difference and reports what it
finds marked `last_resort`, which is a lead rather than a result: an assay resting on one mismatch is the
hardest kind to make work. So this message now means one of two things — the run was given
`--no-last-resort`, or **no contig differs from the exclusion genome at all**, which is a different
problem: the regions are specific to the inclusion group as a set of kmers but identical to the one genome
they were mapped to. Try another reference (`-r`), since the rule is judged against one genome only.

## "last_resort" in the headers, and in a warning

The run found nothing under its usual rule and reported what a single difference would give instead. Treat
those regions as leads: check them against the whole exclusion group yourself, and expect
`--min-primer-differences` to set aside any assay designed on them unless its amplicon is absent from the
exclusion genomes altogether. `--no-last-resort` turns the pass off.

## "No candidate contig is present in all the inclusion genomes"

The candidates came from the kmers of only some of the inclusion genomes. `4_blast/inclusion_blast_hits.tsv`
shows exactly which genome is missing which candidate; it is usually one or two genomes with incomplete
assemblies. This can also happen when `-d` is higher than 1, since a kmer occurring several times in one
genome then passes the counting step. `-p 0.9` keeps those candidates, marked with the number of genomes
that hold them.

The last-resort pass answers this message too, not only the one above: a contig that differs by a single
base may be shared by every inclusion genome where the ones with two differences are not. So seeing this
error means the weaker rule was tried as well and found nothing either — the message says so — or that
`--no-last-resort` was given.

## "No contig passed every filter"

The run finished, but `final_kmers.fasta` is empty: the differences found against the single reference
genome did not hold up against the whole exclusion group. The 90% rule means a position has to differ in
nearly every exclusion genome that matches the contig. Comparing `4_blast/all_inclusion_contigs.fasta` —
what went into this step — with the empty `final_kmers.fasta` shows what was dropped; neither file needs
`--keep-intermediate`.

With the defaults this message means the last-resort pass found nothing either — it logs a line of its own
saying so, and `run_info.json` then reads `min_differences: 1` with `last_resort: false`. With
`--no-last-resort` it means only that the strict rule found nothing.

## "The output folder path contains a space" / "the output folder and the inclusion folder overlap"

Two things the run refuses before doing any work:

- BLAST cannot open a database whose path holds a space, and it only finds out after writing it, so the
  output folder must have a space-free path. The **input** folders may contain spaces: each genome is linked
  into the output folder before blast sees it.
- The output folder must be outside both input folders. They are searched recursively, so results written
  inside one of them would be picked up as input genomes by the next run.

The space rule applies to `primer-finder design -o` as well, for the same reason: that command builds BLAST
databases of its own under its output folder.

## I have more than 255 inclusion genomes

That works. KMC stores a count in one byte by default, and it applies `-ci` to that saturated value, so
"present in all 300 genomes" would match nothing; primer-finder raises the ceiling (`-cs`) whenever the
inclusion count times `-d` passes 255 — which with `-d 2` is 128 genomes, not 256. Before 1.0.0 such a run
reported that no inclusion-specific kmer existed.

## Which kmer size should I use?

99 is a good default for bacteria: long enough to be specific, short enough to be found. Shorter kmers
(31, 51) give **fewer** candidates, not more, and at 51 and below whole runs start coming back empty — the
measurement is in [Methods](Methods#1-the-kmers-only-the-inclusion-group-has-kmc). Longer ones are stricter
and more sensitive to assembly error. The kmer size
also sets the minimum length of a specific region that can be found at all: a region shorter than `k` cannot
yield a kmer of its own.

## How many genomes do I need?

Enough inclusion genomes to cover the diversity of the group you want to amplify (a handful is already
useful, a dozen is better), and as many exclusion genomes as you can, especially the near neighbours — those
are what the assay has to discriminate against. Adding distant exclusion genomes costs little and changes
little.

## Can I use reads instead of assemblies?

No. primer-finder expects assembled genomes. Reads would add sequencing errors to the kmer counts, which the
exact-match rules cannot tolerate.

## Does it design the primers?

Yes, since 1.2.0: `primer-finder design results/ -o assays/` runs Primer3 on the regions, keeps the assays
that could tell the two groups apart, and writes the primer files
[insilicoPCR](https://github.com/duceppemo/insilicoPCR) reads — and runs it for you with `--insilico-pcr`.
Every oligo it proposes has already passed Primer3's length, melting temperature and GC limits, and
hairpin and dimer limits for the oligo on its own and paired with its partner; and unless the request
pinned a 3' end on a differing base, each primer carries a G or a C at that end — so an assay that cannot
work is not offered. See [Designing assays](Designing-assays), which lists all of those,
the reaction they are predicted for, and explains how the assays are ordered and how much that order is
worth.

`find` on its own still just reports the regions, with the specific bases in lower case, if you would rather
design them yourself.

## Why is the result not exactly the same as in the 2022 version?

Several rules were fixed; they are listed in
[CHANGELOG.md](https://github.com/duceppemo/primer-finder/blob/main/CHANGELOG.md). The ones that change
which contigs come out:

- the exclusion kmers are now really all subtracted. The old command line passed `-cx1e9` to KMC, which KMC
  reads as `1`, so only the kmers seen **exactly once** in the whole exclusion group were subtracted;
  candidates in repeated exclusion regions survived that they no longer do;
- the reference exclusion genome is the first one alphabetically instead of a random one, so a run can be
  repeated;
- the variant positions in the headers are positions in the contig, not in the blast alignment;
- contigs that map to the minus strand have the right bases marked;
- the two-differences rule covers insertions and deletions, differences that touch, and clipped ends. It
  was "within 21 bases" in 1.0.0; since 1.7.0 the pair may span 25 bases, the widest primer, while the
  clipped end is still 21;
- a position has to differ in every copy of the region a given exclusion genome holds, and genomes that do
  not hold the region at all no longer count towards the 90%;
- two shared variant positions are enough to keep a contig; the old code needed four because of an
  off-by-three in its loop.

## Can I run it on a cluster?

Yes; it is a single process that runs command-line programs. `-t` and `-m` are capped at what the control
group a scheduler puts the job in actually allows, so asking for more than the allocation only logs a
warning. There is no resuming: a run that is killed starts again from the beginning.

## Where do I report a problem?

[Open an issue](https://github.com/duceppemo/primer-finder/issues/new/choose) with `primer_finder.log`
(rerun with `--debug`) and `run_info.json` from the output folder.
