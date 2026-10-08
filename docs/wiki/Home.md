# primer-finder

primer-finder finds the sequences that tell one group of assembled genomes from another, so that a selective
(q)PCR assay can be designed on them.

Give it two folders of genomes — the ones the assay should amplify (**inclusion**) and the ones it must not
(**exclusion**) — and it:

1. counts the kmers **shared by every inclusion genome** and the kmers of **any** exclusion genome (KMC),
   and subtracts the second set from the first;
2. **assembles** the inclusion-specific kmers into contigs (SKESA, or SPAdes), so that there are fewer and
   longer sequences to examine;
3. **maps** the contigs to one exclusion genome (minimap2) and keeps those whose differences — mismatches,
   insertions, deletions — are close enough together to fit in one primer or probe;
4. **checks** every candidate against all the genomes with blast: it must be present in every inclusion
   genome, and its differences must hold up against every exclusion genome, not just the one it was mapped
   to.

The answer is `final_kmers.fasta`: one record per candidate region, with the specific bases in lower case and
their positions in the header, the most promising first. Design the assay on those bases — ideally with the
variable positions in the primers, at their 3' end, rather than in the probe.

primer-finder only keeps perfect matches: a kmer must be in **all** the inclusion genomes with no mismatch,
and in **none** of the exclusion genomes. It is therefore very sensitive to the quality of the assemblies and
to how the genomes were assigned to the two groups. Curate the input genomes;
[genome_comparator](https://github.com/duceppemo/genome_comparator) helps with that.

| Page | Contents |
|---|---|
| [Installation](Installation) | conda, bioconda, from source, checking the installation |
| [Usage](Usage) | inputs, every option, choosing the two groups, performance |
| [Methods](Methods) | what each step does, the filtering rules, limits |
| [Outputs](Outputs) | every file and field |
| [Example](Example) | the simulated dataset and its expected result |
| [Validation](Validation) | finding three published qPCR assays in public genomes |
| [FAQ](FAQ) | troubleshooting, "no contig passed" |
| [Development](Development) | tests, continuous integration, releases |
