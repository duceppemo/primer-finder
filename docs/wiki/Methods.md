# Methods

What each step does, and the exact rule it applies. The step numbers match the subfolders of the output.

## 1. The kmers only the inclusion group has (KMC)

```
kmc -k99 -t<threads> -m<memory> -fm -ci<N> -cx<N × d>  @inclusion_list.txt  inclusion  work/
kmc -k99 -t<threads> -m<memory> -fm -ci1 -cx1000000000 @exclusion_list.txt  exclusion  work/
kmc_tools -t<threads> simple inclusion exclusion kmers_subtract inclusion_specific
kmc_tools -t<threads> transform inclusion_specific dump kmers.txt
```

KMC counts the kmers of all the files of a group together, so with `N` inclusion genomes a kmer that occurs
at least `N` times (`-ci N`) is, in practice, a kmer every inclusion genome carries. `-p/--min-inclusion`
lowers that threshold to a fraction of `N`, rounded up, for the regions that most but not all of the
inclusion genomes share. With more than 255
inclusion genomes, `-cs` is added as well: KMC's counters stop at 255 by default, so without it `-ci 300`
would match nothing at all. `-cx N × d` (`-d`,
default 1) drops the kmers that occur more often than that, which are the repeated regions: they would give
an assay that amplifies several places at once. The trade-off is that a kmer which happens to occur `N`
times inside fewer genomes also passes this step — which is why step 4 checks the presence of every
candidate in every inclusion genome with blast.

For the exclusion group, every kmer counts, however rare: `kmers_subtract` then removes from the inclusion
set every kmer seen anywhere in the exclusion group. What is left is the inclusion-specific kmers, written as
a fasta file (`kmer_0`, `kmer_1`, ...).

A shorter `-k` finds shorter specific regions and more of them; a longer one is more specific and more
sensitive to assembly errors. 99 is a good starting point for bacteria.

## 2. Assembling the kmers (SKESA or SPAdes)

```
skesa --cores <threads> --mem <memory> --fasta kmers.fasta --contigs_out assembly.fasta
# or
spades.py --s 1 kmers.fasta --isolate --only-assembler --threads <threads> --memory <memory> -o spades/
```

The kmers overlap one another by `k-1` bases, so assembling them turns thousands of kmers into a handful of
contigs: the specific regions themselves. This is only a way to have fewer and longer sequences to examine —
the sequences are not new information. SKESA is the default and is faster here; SPAdes (`-a spades`)
sometimes joins regions that SKESA leaves apart. `--only-assembler` is used because kmers carry no quality
values to correct.

## 3. The differences worth an assay (minimap2)

```
minimap2 -t <threads> -a --eqx <one exclusion genome> assembly.fasta > mapping.sam
```

`--eqx` makes minimap2 write `=` and `X` in the cigar string instead of `M`, which is what lets primer-finder
tell a match from a mismatch. The contigs are mapped to **one** exclusion genome: the first of the exclusion
folder in alphabetical order, or the one given with `-r`. Only the primary alignment of each contig is read.

A contig is kept when, judged on that cigar string:

- it does not map at all — the whole contig is missing from that exclusion genome, the best case; or
- a mismatch run, an insertion or a deletion is longer than one base; or
- two differences (mismatch, insertion or deletion) are fewer than 21 matching bases apart, including two
  that touch, and including two separated only by a deletion (a deletion takes up no base of the contig); or
- an end of the contig is clipped over at least 21 bases, which means the aligner could not place a whole
  primer's worth of sequence anywhere in that genome.

21 bases is about the length of a PCR primer: two differences that close can sit in the same primer or probe,
which is what makes an assay selective. A single mismatch somewhere in the contig is not enough, because a
primer carrying one mismatch still amplifies the exclusion template often enough to be useless.

The kept contigs are written to `best_kmers.fasta`, the ones with the most differing bases first — which puts
a contig the exclusion genome does not hold at all at the top — with those bases in lower case. Mismatched,
inserted and clipped bases are marked; a deletion is in the cigar string only, since those bases are not in
the contig. When a contig maps to the minus strand, the marks are put back on the forward sequence, which is
the one written out.

Because only one exclusion genome is used here, the differences found may be particular to that genome.
Step 5 is what checks them against the whole exclusion group.

## 4. Present in every inclusion genome (blast)

```
makeblastdb -in <each inclusion genome> -dbtype nucl -out <output folder>/...
blastn -db <genome> -query best_kmers.fasta -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid evalue"
```

Each inclusion genome is searched for each candidate; a hit with an e-value of 1e-10 or better counts as
present. The presence of every candidate in every genome is written to `inclusion_blast_hits.tsv`, and a
candidate missing from even one inclusion genome is dropped: an assay must amplify all of them. This is
where the kmer counting of step 1 is checked properly, since "occurs at least `N` times" is not quite "is in
`N` genomes". With `-p/--min-inclusion`, a candidate has to reach the same number of genomes as the kmers
did; the ones that survive say how many hold them (`inclusion=11/12`) and the widest are listed first. The genomes
are searched in parallel (`-t`), and the blast databases are built inside the output folder, never next to
the input genomes.

Gzipped genomes are decompressed into the output folder for this step: blast cannot read gzip.

## 5. Different from every exclusion genome (blast)

```
blastn -db <each exclusion genome> -query all_inclusion_contigs.fasta -evalue 1e-10 \
       -max_target_seqs 10 -outfmt "6 qseqid qstart qend evalue qseq sseq"
```

Every hit of a contig in an exclusion genome is read, not just the best one, and each gives the positions of
the contig that differ from that copy: mismatches, bases the contig has and the genome does not, and the
position where bases the genome has are missing from the contig. The aligned sequences blast returns are
walked base by base, so the positions are positions **in the contig**, counted from 0.

A position counts as different in a genome only when **every copy** of that region in the genome differs
there: one matching copy is enough for the assay to amplify that genome. A position is then kept when it
differs in at least **90%** of the exclusion genomes that align it — a few exclusion genomes carrying the
inclusion allele do not spoil an assay, and the threshold leaves room for them. A genome whose hits do not
reach that position is left out of the count rather than counted as identical: it does not hold that region,
which only makes the assay more selective. A contig is kept when:

- no exclusion genome hits it at all (it is specific on its own), in which case it is written as it came out
  of step 3; or
- more than one kept position remains, and two of them are fewer than 21 bases apart.

The survivors are `final_kmers.fasta`, with those positions listed in the header and in lower case in the
sequence.

## Designing the assay

Design the primers and the probe so that the lower-case positions fall in the **primers**, ideally at their
3' end, rather than in the probe: a mismatch there costs the exclusion template more amplification than one
under the probe. Primer3 and an in-silico PCR against the exclusion genomes are the natural next steps;
primer-finder does not do them (yet).

## Limits

- **Perfect matches only.** A kmer must be in every inclusion genome without a single mismatch, and in no
  exclusion genome. One bad assembly in either group can remove every candidate. Draft assemblies with
  missing regions are the usual reason a run finds nothing. `-p/--min-inclusion` relaxes the inclusion side
  of this, at the price of an assay that does not amplify every inclusion genome; nothing relaxes the
  exclusion side, where a single matching genome makes an assay useless.
- **Assembled genomes only.** Reads are not an input; assemble them first.
- **One exclusion genome in step 3.** The candidate list depends on which genome that is (`-r`). The
  differences are then re-checked against the whole group, so this affects what gets looked at, not what
  survives.
- **No assay is designed or tested.** The output is the region to design on, not a primer pair.
- **Plasmids and repeats.** `-d 1` discards anything repeated inside a genome; a kmer carried by a plasmid
  that moves between the two groups will not be specific.
- **Where the output folder may be.** It must be outside both input folders — they are searched recursively,
  so results written inside one would come back as input genomes — and its path must not contain a space,
  which BLAST cannot handle in a database path. The input folders themselves may contain spaces: every
  genome is linked into the output folder before blast sees it.
