# Outputs

```
results/
├── final_kmers.fasta          the answer: the candidate regions
├── run_info.json              parameters, genomes, program versions, counts
├── primer_finder.log          the whole log of the run
├── 1_kmers/
│   ├── inclusion_list.txt     the inclusion genomes, as given to KMC
│   ├── exclusion_list.txt     the exclusion genomes, as given to KMC
│   ├── inclusion_specific_99mers.fasta   the kmers only the inclusion group has
│   ├── kmers.txt              the KMC dump (--keep-intermediate)
│   ├── inclusion.kmc_pre/suf, exclusion.kmc_pre/suf, inclusion_specific.kmc_pre/suf
│   │                          the KMC databases (--keep-intermediate)
│   └── kmc_work/              KMC's temporary folder (--keep-intermediate)
├── 2_assembly/
│   └── assembly.fasta         those kmers assembled into contigs
├── 3_candidates/
│   ├── best_kmers.fasta       the contigs whose differences could carry an assay
│   └── mapping.sam            their alignment to the exclusion genome (--keep-intermediate)
└── 4_blast/
    ├── inclusion_blast_hits.tsv   each candidate against each inclusion genome
    ├── all_inclusion_contigs.fasta  the candidates present in every inclusion genome
    ├── inclusion_db/, exclusion_db/  blast databases and raw hits (--keep-intermediate)
```

Without `--keep-intermediate`, the files marked above are removed when the run finishes; everything else is
kept. A run that stops on an error keeps them all, whatever the option says, so that the step that failed can
be examined; `run_info.json` is only written by a run that reaches the end.

## final_kmers.fasta

One record per candidate region, the most promising first (the ones with the most differences):

```
>Contig_1_67.5167 [72, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 172, 179, 187]
CACCAAAATCCTCCTGGCGAGAGGTTAGATCATCAGAGTCCGTGTACAGAAACGGCCATTGCGAACTACTGTcCCATGTGAACT...
```

- The **name** is the contig name the assembler gave (SKESA puts its length and coverage in it).
- The **header** lists the positions that differ from the exclusion genomes, counted from 0 in this
  sequence.
- Those positions are **lower case** in the sequence; everything else is upper case. Design the primers so
  that the lower-case bases fall in them, ideally at their 3' end.

A contig that no exclusion genome matched at all keeps the description it had in `best_kmers.fasta` (see
below) instead of a list of positions: all of it is specific.

With `-p/--min-inclusion` below 1, the description also holds `inclusion=11/12`: how many of the inclusion
genomes hold this region. The records are then ordered by that count first, so the regions that cover the
most inclusion genomes come first. Without the option every region is in every inclusion genome and no tag
is written.

## best_kmers.fasta

The candidates after the mapping step, before blast, the ones with the most differing bases first. Their
description is the cigar string of their alignment to the one exclusion genome they were mapped to, which
reads as a summary of the differences:

```
>Contig_1_67.5167 72=1X7=1X6=1X34=15I50=1X72=
```

`=` match, `X` mismatch, `I` bases the contig has and the exclusion genome does not, `D` the other way
round, `S` soft-clipped ends. A contig that did not map at all is described as `<length>I`.

Mismatched, inserted and clipped bases are in lower case here too. Deleted bases cannot be marked — they are
not in the contig — so they only show in the cigar string.

## inclusion_blast_hits.tsv

Which inclusion genome holds which candidate: `1` present (a blast hit at 1e-10 or better), `0` absent. The
candidates with a `0` are the ones dropped at this step; the table is the place to look when a region you
expected disappeared. A column is named after the genome's path inside the inclusion folder, so genomes with
the same file name in different subfolders (`sampleA/contigs.fasta`, `sampleB/contigs.fasta`) stay apart.

```
contig	inclusion_1.fasta	inclusion_2.fasta	inclusion_3.fasta
Contig_1_67.5167	1	1	1
Contig_2_41.2	1	0	1
```

## run_info.json

Everything needed to repeat or report the run:

| Field | Contents |
|---|---|
| `version` | the version of primer-finder |
| `command_line` | the command as it was typed |
| `parameters` | the input folders, output folder, `kmer_size`, `duplication`, `assembler`, `threads`, `memory_gb` after capping |
| `inclusion_genomes`, `exclusion_genomes` | every genome file used, in the order they were given to KMC |
| `reference` | the exclusion genome the contigs were mapped to |
| `programs` | what each program reported as its version |
| `counts` | `kmers` (inclusion-specific kmers), `contigs` (assembled), `candidates` (after mapping), `in_all_inclusion`, `final` |
| `seconds` | how long the run took |

## primer_finder.log

The same lines the run printed, with timestamps. With `--debug`, also every command line and the output of
every program — this is what to attach to a bug report, together with `run_info.json`.


## The `design` command

`primer-finder design` writes its own folder:

```
assays/
├── assays.tsv              every assay, the usable ones first
├── assays_qpcr.fasta       the assays that have a probe, for insilicoPCR in qPCR mode
├── assays_pcr.fasta        every assay, primers only, for insilicoPCR in standard PCR mode
├── run_insilico_pcr.sh     the in silico PCR of both files against both groups
├── design_info.json        parameters, programs and versions, counts
├── primer_finder_design.log
├── primer3/                what was sent to Primer3 and what came back
├── amplicons/              the blast checks of each amplicon against both groups
└── insilico_pcr/           what insilicoPCR wrote, when it was run
```

### assays.tsv

One row per assay. The usable ones come first, in the order [Designing assays](Designing-assays) explains;
the ones that would amplify both groups are kept at the end for the record.

| Column | What it is |
|---|---|
| `assay` | its name, which is also the prefix of its oligos in the fasta files (`<assay>-F`, `-R`, `-P`) |
| `region` | the candidate region it was designed in |
| `specific_by` | `absence` (no exclusion genome holds the amplicon), `difference` (an oligo sits on differing bases), or `nothing` (it would amplify both groups) |
| `exclusion_genomes_with_amplicon` | how many exclusion genomes hold the amplicon in one piece |
| `best_terminal_run` | the longest run of differences ending at a primer's 3' end |
| `inclusion_copies_min`, `inclusion_copies_max` | copies of the amplicon per inclusion genome; above 1 only for a run made with `-d 2` or more |
| `product_size`, `penalty`, `penalty_band` | the amplicon length and Primer3's pair penalty, with the band the ranking uses |
| `forward`, `reverse`, `probe` | the sequences, with `_start`, `_tm`, `_gc` and `_variants` (differences covered) for each |
| `forward_hairpin_tm`, `forward_self_dimer_tm` | and the same for `reverse` and `probe`: where Primer3 predicts the oligo folds on itself, and where it pairs with a copy of itself, in °C. `0.0` means no structure was predicted; both are below `--max-dimer-tm` / `--max-hairpin-tm` or the oligo was never proposed |
| `pair_dimer_tm`, `pair_dimer_end_tm` | where the two primers pair with each other, over the whole oligo and at their 3' ends, in °C |
| `primer_variant_weight` | the differences under the primers, weighted by the base each replaces: 1 for a G or C, 0.5 for an A or T, 0.75 when unknown |
| `forward_strong_variants`, `forward_terminal_run`, `forward_near_3prime` | and the same for `reverse`: how many of its differences replace a G or a C, and where they sit relative to the 3' end |
| `inclusion_total`, `exclusion_total` | how many genomes were tested, when in silico PCR was run |
| `qpcr_inclusion_amplified`, `qpcr_inclusion_percent`, `qpcr_exclusion_amplified`, `qpcr_selective` | what insilicoPCR said in qPCR mode, with the probe |
| `pcr_inclusion_amplified`, `pcr_inclusion_percent`, `pcr_exclusion_amplified`, `pcr_selective` | and in standard PCR mode, primers only |

The structure temperatures are predicted for the reaction the run was given (`--monovalent`, `--divalent`,
`--dntp`, `--primer-conc`, `--probe-conc`), which `design_info.json` records under `parameters.conditions`.

`selective` is `yes` when the assay amplified every inclusion genome and no exclusion genome, `partial
(88%)` when it reached the `-p/--min-inclusion` the regions were found with but not all of them, and `no`
otherwise. Amplifying an exclusion genome is `no` whatever the threshold.
