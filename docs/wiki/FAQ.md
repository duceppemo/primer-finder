# FAQ

## "No inclusion-specific kmer was found"

Every kmer the inclusion genomes share is also somewhere in the exclusion group. Things to try, roughly in
order:

1. **Check the groups.** One exclusion genome that actually belongs to the inclusion group is enough to
   remove every candidate. The reverse is just as bad.
2. **Check the assemblies.** A kmer must be present without a single mismatch in *all* inclusion genomes, so
   a draft assembly with a gap in the right place, a contaminated one, or one of much lower quality than the
   others silently removes candidates. Compare the genomes first
   ([genome_comparator](https://github.com/duceppemo/genome_comparator)).
3. **Use shorter kmers** (`-k 51`, `-k 31`): a shorter kmer is less likely to be broken by a sequencing
   error or a small variant, and finds shorter specific regions.
4. **Raise `-d`** (`-d 2`, `-d 3`): a kmer that occurs twice in some inclusion genome is kept.
5. **Reduce the inclusion group** to the genomes you are most confident about, if a few of them are of
   poorer quality.

## "No contig differs enough from the exclusion genome"

Inclusion-specific kmers were found and assembled, but none of the contigs has two differences within 21
bases of each other when mapped to the exclusion genome. The region is specific, but by a single mismatch
here and there, which does not make a selective assay. Try another exclusion genome as the reference (`-r`),
or shorter kmers so that other regions come through.

## "No candidate contig is present in all the inclusion genomes"

The candidates came from the kmers of only some of the inclusion genomes. `4_blast/inclusion_blast_hits.tsv`
shows exactly which genome is missing which candidate; it is usually one or two genomes with incomplete
assemblies. This can also happen when `-d` is higher than 1, since a kmer occurring several times in one
genome then passes the counting step.

## "No contig passed every filter"

The run finished, but `final_kmers.fasta` is empty: the differences found against the single reference
genome did not hold up against the whole exclusion group. The 90% rule means a position has to differ in
nearly every exclusion genome that matches the contig. Looking at `3_candidates/best_kmers.fasta` with
`--keep-intermediate` shows what was dropped and why.

## "The output folder path contains a space" / "the output folder and the inclusion folder overlap"

Two things the run refuses before doing any work:

- BLAST cannot open a database whose path holds a space, and it only finds out after writing it, so the
  output folder must have a space-free path. The **input** folders may contain spaces: each genome is linked
  into the output folder before blast sees it.
- The output folder must be outside both input folders. They are searched recursively, so results written
  inside one of them would be picked up as input genomes by the next run.

## I have more than 255 inclusion genomes

That works. KMC counts to 255 by default, which would make "present in all 300 genomes" impossible to ask
for, so primer-finder raises the ceiling (`-cs`) when the group is that large. Before 1.0.0 such a run
reported that no inclusion-specific kmer existed.

## Which kmer size should I use?

99 is a good default for bacteria: long enough to be specific, short enough to be found. Shorter kmers (31,
51) give more and shorter candidates and tolerate variation better; longer ones are stricter. The kmer size
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

No. It reports the regions, with the specific bases in lower case, and the design is up to you (Primer3 and
an in-silico PCR against the exclusion genomes). Put the specific bases in the primers, at their 3' end
rather than in the probe.

## Why is the result not exactly the same as in the 2022 version?

Several rules were fixed; they are listed in
[CHANGELOG.md](https://github.com/duceppemo/primer-finder/blob/master/CHANGELOG.md). The ones that change
which contigs come out:

- the exclusion kmers are now really all subtracted. The old command line passed `-cx1e9` to KMC, which KMC
  reads as `1`, so only the kmers seen **exactly once** in the whole exclusion group were subtracted;
  candidates in repeated exclusion regions survived that they no longer do;
- the reference exclusion genome is the first one alphabetically instead of a random one, so a run can be
  repeated;
- the variant positions in the headers are positions in the contig, not in the blast alignment;
- contigs that map to the minus strand have the right bases marked;
- "two differences within 21 bases" now covers insertions and deletions, differences that touch, and
  clipped ends of at least 21 bases;
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
