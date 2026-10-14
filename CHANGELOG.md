# Changelog

## Unreleased

### Added

- **A `suis_lato` target, and the measured reason there is no *B. suis* assay** (record in
  `validation/results/2026-10-14_v1.7.0_brucella/SUIS_LATO.md`). *B. canis* is a lineage inside
  *B. suis*, so the `suis` target reports nothing at any setting; this target puts both species in the
  inclusion group. The group definitions of all six targets now live in one place,
  `validation/brucella/groups.py`, which the four scripts that need them import — the census output is
  byte-identical to the committed one after that refactor.

  It still finds no assay, in three separate places. At the default `-k 99` there is no inclusion-specific
  kmer at all, on either panel, at `-p` 1.0, 0.99 or 0.95. At `-k 201` with `-p 0.95` there are 15, which
  overlap into **one 215 bp locus** — and SKESA will not build a contig from 15 kmers, the same limitation
  as the `-k 31` failure already recorded, from the other end of the range. The locus itself is in the
  record as a fasta file, and the find rule would reject it anyway: it differs from *B. abortus* by **one**
  base and from *B. melitensis* by **two, 34 bases apart**, where the rule asks for two within 25.

  Also measured: the published BS assay carries all three oligos in 120 of 120 *B. suis* but **1 of 74
  *B. canis***, so it is a sensu stricto assay and nothing from this target would replace it.

  The wider conclusion for this panel: across every target, species-level discrimination in the classical
  *Brucella* clade rests on single scattered SNPs rather than blocks of unique sequence, which is what the
  find rule declines by design. The genus level is the opposite case, and there primer-finder does well.

- **A `genus_wide` target, and a proposed genus assay that spans the whole classical clade**
  (record in `validation/results/2026-10-14_v1.7.0_brucella/WHOLE_CLADE.md`). The four classical-clade
  species the paper's panel does not hold — *B. pinnipedialis*, *B. amazoniensis*, *B. vulpis* and
  *B. inopinata* — are now labelled `classical_extra` by `panel.py` instead of being held out, and
  `group_genomes.py` has a target that puts all twelve species in the inclusion group. The eight-species
  `genus` target is unchanged, so the record of the published assays stays reproducible.

  At the default `-p 1.0`, **56 regions** are present in every one of the 1,393 genomes and absent from
  all 315 relatives, where the eight-species group gives 398: adding 17 genomes of four species removes
  86% of the candidates, most of it *B. inopinata*. `design` proposed 334 assays on them, and
  **`Contig_83_60.1233_assay19` carries all three oligos exactly in all 1,393 genomes of all twelve
  species, in none of the 315 relatives, with no predicted hairpin or dimer at all.** The published
  genus assay scores 1,388 of 1,393 on the same group.

  This reverses the conclusion of `ALTERNATIVES.md` for the genus assay: asked the eight-species
  question, every candidate that reached 100% gave up species outside the panel, so none was worth
  acting on; asked the twelve-species question, one gives up nothing measured. It still has no wet-lab
  data, and its primers are more AT-rich than the published ones, which is the kind of difference only a
  reaction settles.

- **Two more tools in the *Brucella* harness, and a record of what primer-finder proposes of its own**
  (`validation/brucella/score_assays.py`, `thermo.py`, record in
  `validation/results/2026-10-14_v1.7.0_brucella/ALTERNATIVES.md`). `score_assays.py` puts any set of
  assays — published or proposed — against the whole panel on one yardstick; `thermo.py` measures Tm and
  secondary structure for any oligos under the design command's own default conditions, each oligo in a
  request of its own so that one value over a limit does not blank the rest.

  Asked for alternatives to the published assays, `design` proposes 520 for the genus target and 15 and 98
  for *B. abortus* and *B. melitensis*. Five genus candidates reach 100% in-silico coverage of the panel
  species where the published assay reaches 99.85%, with nothing in the exclusion group; the species
  candidates gain one and two genomes and rest on two differences under one primer rather than one. None
  is an improvement worth acting on: the margins are a fraction of a percent, the genus candidates trade
  panel coverage for breadth across the genus, and only the published assays have wet-lab validation.

  The useful finding is about the defaults. **Three of the four published assays exceed
  `--max-hairpin-tm` 47 °C** — the *B. melitensis* forward primer is predicted at 76.4 °C — and all four
  are validated on 235 target and 192 non-target strains. `design` would never have proposed them. That is
  the only wet-lab evidence in the validation suite bearing on that limit, and it says the limit is
  conservative enough to reject assays that work.

  Also measured: no *B. suis* region exists on this panel at any `-k` or `-p` tried, because *B. canis*
  sits inside *B. suis* and 74 *B. canis* genomes are in the exclusion group. That target needs its groups
  redrawn rather than its settings changed.

- **A third validation record, and the first against assays this software itself produced**
  (`validation/brucella/`, record in `validation/results/2026-10-14_v1.7.0_brucella/`). The four published
  *Brucella* assays of Hoover et al. 2026 (doi:10.1128/spectrum.02974-25) were designed on unique core
  sequences found with `selective_primer_finder 0.2.1`, which this package is the rewrite of.

  The genus-level assay is recovered, at `-p 0.99`, in a **127 bp** region — which 1.5.0 could not have
  reported at all, since both assemblers were left at their own 200 bp minimum until 1.6.0. On the
  *Listeria* panel removing that floor changed no verdict; here it is the difference between recovering a
  published assay and not. `design`, run on that region, proposes a primer pair sharing 18 of 19 bases
  with the published forward primer and an amplicon of 92-94 bp against the published 93.

  The three species assays are **not** recovered, and the record measures why: they rest on one mismatch
  (*B. abortus*), on two mismatches 114 bases apart (*B. melitensis*), or on a 21-base stretch that does
  not align (*B. suis*), while the find rule asks for two differences within 25 bases so that one primer
  can carry both. The kmer step does find the signatures — 16 of the 25 51-mer windows of the *B. abortus*
  amplicon are among the inclusion-specific 51-mers — so the limit is that rule, which 0.2.1 did not have.
  Recovering this class of assay would need the `find` side to keep a region on a single difference; that
  is not changed here.

  Two limitations the record exposes, neither fixed yet: at `-k 31` SKESA cannot assemble the specific
  kmers at all ("Was not able to assemble anything", with `--min_count 1` and `--kmer 21` making no
  difference), so the advice the "too few kmers" error gives runs into a different wall and a region
  shorter than roughly 50 bases cannot come out of the pipeline whatever the kmer size; and one genome
  filed under the wrong species removes every kmer of the region that distinguishes it — on the panel as
  labelled, the species targets return nothing at all, where holding out eight demonstrably mislabelled
  assemblies gives 1 to 11 regions.

## 1.7.0

### Fixed

- **`design --insilico-pcr` could not check a gzipped or nested genome folder, and failed the whole run
  when it met one.** insilicoPCR reads only the top level of the folder it is given and hands each file to
  its own BLAST as it stands, so a genome in a subfolder was silently skipped and a `.fasta.gz` made that
  BLAST call fail, which took the command down with it. primer-finder accepts both everywhere else, and
  the bundled example holds a `.fasta.gz`, so the documented quick start was the case that broke. Both
  groups are now staged into `insilico_genomes/` as one flat folder of plain fasta per group — subfolders
  flattened, gzipped genomes written out plain, a name used twice numbered so two samples stay two — and
  the script written for running the check later points at the same folders.
- **A failure in the in-silico PCR check discarded the whole design run.** Everything before the check is
  the expensive part of the command, so `assays.tsv` is now written without the verdict columns and the
  command exits 1 with the reason, rather than raising and leaving nothing behind.
- **A run that stopped short left its KMC databases behind**, and an exclusion database for a few hundred
  genomes is a gigabyte. Clean-up now happens on the way out either way, not only after a run that
  finished.
- **`-m 1` was refused by KMC at run time.** Its own help says `-m` takes 1 to 1024, but 1 gives "min
  memory must be at least 2GB", so the floor is now 2 — which also keeps a tight control group or a small
  machine from computing a default of 1 and failing on the first KMC call.
- **An empty kmer result blamed the exclusion genomes even when the inclusion genomes shared nothing to
  begin with.** The two causes want opposite advice, so they now give different messages: the inclusion
  group may be too diverse for this kmer size, or the exclusion group really does hold every shared kmer.
  The count that tells them apart is only taken when the result is empty.
- **An IDT order sheet whose amplicon size reads "1,234" named its oligos after a 1 bp amplicon**, because
  only the "1" was read. The thousands separator is now part of the number.
- **`design` accepted an output path containing a space**, then failed part way through on a BLAST call
  that cannot quote one — `find` had refused it up front since 1.1.0. It also now refuses `-t` below 1 and
  caps `-t` at the number of CPUs, as `find` does.

### Changed

- **The find step keeps a contig when two of its differences could sit under one primer, which now means a
  span of 25 bases rather than 22.** The rule and the design step disagreed: `find` asked for at most 20
  matching bases between a pair, a span of 22, while the design step's widest primer is 25 and its pinned
  requests aim at 25-base windows. Pairs spanning 23, 24 or 25 bases were therefore thrown away although a
  primer could carry both. On the *Listeria* panel the wider rule gives **332 candidate regions instead of
  306**, +8.5%, with no target losing any; the four recoverable published assays are still recovered, at
  identical ranks. Measured in `validation/results/2026-10-14_v1.6.0/`.

  `mapping.PRIMER_LENGTH` is replaced by two constants, because it was doing two jobs with one number:
  `MAX_PRIMER_SPAN` = 25 is the pair rule and the window of the exclusion step, and `MIN_CLIPPED` = 21 is
  how much of a contig end the aligner must fail to place for that end alone to be worth an assay. The
  second keeps its old value: it is a different quantity, and raising it would only throw regions away.
- Four more documentation corrections, found by reading the pages against the code:
  - the wiki said `-cs` is added "with more than 255 inclusion genomes"; the trigger is `N x d` above 255,
    so 128 genomes at `-d 2`. KMC applies `-ci`/`-cx` to the one-byte saturated count, which is now stated
    from a measurement: with no `-cs`, one kmer repeated 300 times is reported by `-ci 250` with a count of
    255 and not reported at all by `-ci 256`;
  - Methods said the design step carries no further an assay with fewer than two differences "under one
    oligo". It counts primers only, as the same page says correctly two sections earlier;
  - the *Listeria* table on the Validation page was headed "primer-finder 1.6.0" but its region counts came
    from the 25-base span above, which 1.6.0 does not have;
  - the advice for an empty run put shorter kmers near the top, which the measurement on the same page
    contradicts.
- The wiki claimed a shorter `-k` finds "more" regions. It finds **fewer** — 70 regions over the nine
  *Listeria* runs at `-k 31` against 306 at 99, with two runs producing nothing at all at 31 — because a
  kmer survives only if no exclusion genome holds it, so it has to span a difference, and a long kmer spans
  one more often. The default stays at 99: it is near the top of the yield with no empty runs, and `-k 127`
  adds 4% more candidates and no recovered assay.

## 1.6.0

### Fixed

- **A candidate region shorter than 200 bases was never reported.** Neither assembler returns a contig
  below that length unless told to, and the command was not telling them, so an inclusion-specific stretch
  shorter than 200 bases was dropped before anything downstream could see it — silently, since there was
  nothing in the output to say so. A kmer is 99 bases by default, so the floor is now the kmer size: the
  shortest a region can be, and `--min-region` sets it to something else.

  **Only SKESA can be asked for less.** SPAdes has no such option, and its floor is in graph simplification
  rather than in what it writes — measured: stretches of 150, 180 and 199 bases are already absent from its
  simplified graph, while 210 and longer come through. So `-a spades` still does not report the short
  regions, and now warns instead of dropping them quietly.

  This only bites where the two groups are closely related, which is why the *Xylella* validation never
  showed it and the *Listeria* one did: on one clonal complex of a single species the assembly went from 67
  contigs to 398 and the regions reported from 25 to 66, 41 of them shorter than 200 bases. It changed which
  assays that validation recovers not at all — the four it recovered before, it still recovers, at much the
  same ranks — so this is a blind spot removed rather than a result improved. `run_info.json` records the
  value under `parameters.min_region`.
- **The requests that pin a primer's 3' end were skipped on short regions.** The guards asked for the whole
  minimum amplicon *beyond* the pinned base, but an amplicon starts where the primer starts — up to 24 bases
  earlier. On a long region that cost little; on the ~100-base regions the floor above now makes reportable
  it skipped the pins the region was kept for. On a 106-base *Listeria* region with three differences at
  66–75 the command made **no** pinned request at all; it now makes the one that fits, and the best assay
  from that region carries three differences on one primer where it carried two.
- **A difference could be left with no request aimed at it.** The windows were measured once, so a window
  overlapping a richer one was discarded whole rather than measured again over the differences left. With
  differences at 0, 5, 24, 30 and 35 only the 24–35 window was aimed at, and the pair at 0 and 5 — a target
  in its own right — was never used. Of the 432 regions of the *Listeria* runs, 26 held an uncovered pair
  while using fewer than the four targets they are allowed.
- **`--min-primer-differences 1` admitted assays with no difference under either primer.** The filter only
  ran above 1, so at 1 the older, looser rule applied and an assay whose two differences were both under its
  probe was carried forward — into the PCR-mode primer file, where it cannot discriminate at all. It now
  means the same thing at every setting.
- A region that is entirely differences gets no pinned requests. Such a region is specific by absence, so
  there is nothing to aim at, and the eight requests spent on each were the slowest kind.
- `--min-region` refuses a negative value up front instead of failing in the assembler after KMC has counted
  every genome. Dead code removed: `variant_runs`, left behind when targets replaced runs.
- **One region shorter than the smallest product ended the whole design run.** Primer3 answers such a
  request with an error rather than an empty result, and `design` treats a Primer3 error as fatal. Reaching
  it needed a region shorter than `--product-size`'s floor, which is exactly what reporting regions down to
  one kmer produced: `primer-finder design` on a *Listeria* run with `--product-size 120-200` died on the
  first of 57 regions. Regions with no room for a product are now skipped, with a line saying how many and
  why, and it is only an error when none is left.
- The room guard for a pinned reverse primer did not stop at the end of the region, so it could ask Primer3
  for a product longer than the template. Harmless in practice — the plain request for such a region
  failed first — but wrong on its own terms.
- `--min-primer-differences 0` turns the filter off, which is what the "everything was set aside" error now
  suggests when the setting is already 1. The option's help said "1 keeps them", which stopped being true
  when 1 began asking for a difference under a primer.
- **`--max-regions` was applied before the regions with no room were left out**, so a run could design on
  fewer regions than asked while eligible ones went unexamined, and `design_info.json` counted the skipped
  ones as designed on. On a *Listeria* run with `--product-size 120-200 --max-regions 50`, 8 of the first 50
  regions were too short, so it designed on 42 and left 12 eligible regions out. The length filter now runs
  first, and the count is what was designed on.
- **`--product-size` is checked.** It was passed through unvalidated, and since the design step now reads
  the smallest product itself, a reversed range such as `150-70` produced a confident wrong diagnosis —
  "every region is shorter than the 150-base smallest product … a lower --product-size would help" — for
  what is really a malformed argument. A backwards range, a non-number, a zero start and an empty value are
  each refused by name.
- **A truncated gzipped genome said which file it was.** It previously raised `EOFError` from the
  decompressor, which `cli.main` does not catch, so the run ended in a traceback that named no file — and a
  genome folder holds hundreds. It is now a plain error naming the file and saying it is truncated. (`find`
  was unaffected: KMC rejects such a file itself.)

## 1.5.0

### Changed

- **The differences a probe covers no longer qualify an assay.** `--min-oligo-differences` is now
  `--min-primer-differences` (the old name is still accepted) and looks at the primers only, however many
  differences the probe carries. The 2-difference rule was argued from work on primers, and that argument
  does not carry over: a probe here is longer than a primer (18–27 bases against 18–25) and hotter (62–72 °C
  against 58–63 °C), so in a reaction annealing near 60 °C it has far more binding energy in hand.
  Measured, rather than argued: a conventional TaqMan probe still gave a detectable signal through **five**
  mismatches, and under standard conditions neither it nor an MGB probe was sequence-specific
  ([Yao et al. 2006](https://doi.org/10.1016/j.mcp.2006.03.003)); a probe discriminates by being short, a
  12-base MGB probe melting at the same temperature as an unmodified 27-base one (Kutyavin et al. 2000,
  *Nucleic Acids Res* 28:655–661); and in a 5'-nuclease assay the discrimination sits in the primers, where
  one 3'-region mismatch moves the quantification cycle by under 1.5 to over 7 cycles depending on the base
  pair, and by up to sevenfold between master mixes
  ([Stadhouders et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2797725)). An assay whose only
  differences sit under its probe is therefore set aside like any other no primer can tell apart; the probe
  count is still reported and still breaks ties in the ranking. On the *Xylella* set this moves 73 more
  assays to the set-aside list, 40 of them with no primer difference at all. What the evidence does and does
  not settle is set out in the wiki.
- `best_oligo_variants` in `assays.tsv` is renamed **`best_primer_variants`** and now reports the most
  differences either primer covers on its own, which is what the filter reads. `probe_variants` already
  reported the probe's own count.

## 1.4.0

### Changed

- **The `design` command now aims its Primer3 requests at the differences a single oligo could cover,
  rather than at runs of consecutive ones.** This is what the `find` step selects a region for: it keeps a
  region whose differences could sit in one oligo, which means a run of them *or* two of them fewer than 21
  bases apart. A pair seven bases apart is not a run, so the old target missed it and pinned a primer's 3'
  end on a single base instead — wasting the pair the region was kept for. On the *Xylella* validation set
  that was the usual case: 19 of the 20 regions designed on hold such a pair and 8 hold nothing else.
  Aiming at any 25-base window instead raised the share of difference-based assays whose best single oligo
  covers two or more differences from 43% to 77%, everything else unchanged.

### Added

- `--min-oligo-differences`, **2 by default**: an assay that rests on differences rather than on absence
  must cover that many under one oligo, or it is listed but carried no further — not written to the
  insilicoPCR files and not ranked among the candidates. One mismatch, even at a 3' end, often does not
  stop amplification, so an assay that spends only one of a region's close pair is the weak case the
  `find` rule exists to avoid. Assays specific by absence are exempt. `--min-oligo-differences 1` keeps
  them, and when the filter leaves nothing the error says so.
- `best_oligo_variants` in `assays.tsv`: the most differences any one of an assay's oligos covers, which is
  the number the filter and the `find` rule are both about.

## 1.3.0

### Added

- The `design` command rejects the oligos that cannot work, rather than proposing them: a hairpin or a dimer
  — with itself or with the other primer, over the whole oligo and at the 3' end where a polymerase could
  extend it — that melts at or above 47 °C now disqualifies a candidate, for the probe as well as the
  primers. The limits are Primer3's own thermodynamic model, so nothing extra has to be installed, and they
  are adjustable with `--max-hairpin-tm` and `--max-dimer-tm`. What was predicted is reported in
  `assays.tsv` as `pair_dimer_tm`, `pair_dimer_end_tm` and the per-oligo `*_hairpin_tm` and
  `*_self_dimer_tm`.
- The reaction those temperatures are predicted for is now stated instead of left at Primer3's bare
  defaults: 50 mM monovalent cations, 3 mM Mg²⁺, 0.8 mM dNTPs, 250 nM primers and 200 nM probe, which is an
  ordinary TaqMan qPCR. Each is an option (`--monovalent`, `--divalent`, `--dntp`, `--primer-conc`,
  `--probe-conc`) and all of them are recorded in `design_info.json`. They change both the temperatures
  reported and which oligos come back, so a run can be matched to the master mix it is for.
- `--gc-clamp`: how many G or C bases a primer must carry at its 3' end, where extension starts. One by
  default. Primer3 takes one clamp per request, so it is dropped from the requests that pin a 3' end on a
  differing base — that base is whatever the genomes made it, and with the clamp left on there, seven of
  eight forced-end requests on the *Xylella* set returned nothing. Those requests therefore ignore this
  option, and no clamp applies to the probe, which Primer3 has no tag for.
- A verdict of its own, `undecided (3' end)`, for an assay that covers the inclusion group but amplifies
  exclusion genomes only where a difference sits in the last two bases of a primer. `no` would claim the
  check had refused it on evidence worth believing, and it has not: those bases are not counted at any
  tolerance, so there was nothing it could refuse. Those assays rank above the ones refused on evidence the
  check can defend and below the ones it cleared, and `design_info.json` counts them.
- `qpcr_exclusion_terminal_only` and `pcr_exclusion_terminal_only` in `assays.tsv`: of the exclusion genomes
  an assay amplified, how many did so only where a difference sits in the last two bases of a primer — the
  bases insilicoPCR does not count. When it equals `*_exclusion_amplified`, every amplification behind the
  verdict is one the tool could not have refused, so it is the model's blind spot rather than a
  cross-reaction. Only the amplicons that needed no counted mismatch are weighed, which makes the count the
  same at any `--mismatches`, and a genome counts only when none of those is an exact match. On the
  *Xylella* set this separates 96 unjudgeable assays from the 149 the check refuses on its own evidence, and
  the run says so as it goes.

### Changed

- In silico PCR now lets each primer bind through **one** mismatch by default, instead of none
  (`--mismatches`). What that tolerance governs was measured rather than assumed: insilicoPCR does not count
  the last two bases of a primer at all — blast trims an unmatched base off its alignment — never amplifies
  through a mismatch 3 or 4 bases from the 3' end whatever the tolerance, and only counts the mismatches 5
  or more bases in, per primer. A single internal mismatch frequently does not stop a real reaction, so
  assuming it does flattered every assay resting on one difference, and it hid 12 of the 30 absence-based
  assays in the *Xylella* set that amplify something off-target in an exclusion genome. `--mismatches 2` is
  a stricter second pass; 3 measured nothing 2 did not.
- The claim that assays resting on a run of differences at a 3' end were never selective in silico is
  **withdrawn as evidence**. It is still what happens, at every tolerance, but the cause is that the
  alignment cannot see those bases — not anything about PCR. The sweep, the controlled test behind it and
  the numbers that replace it are in `validation/results/2026-10-08_v1.3.0/`.

- A difference is weighed by the base it replaces as well as by where it is: a position where the exclusion
  genomes have a G or a C counts for more than one where they have an A or a T, since G:C holds with three
  hydrogen bonds and A:T with two. The bases come from the exclusion genomes' own alignments of each
  amplicon, and the weight separates only assays covering the same number of differences. `assays.tsv`
  reports it (`primer_variant_weight`, `forward_strong_variants`, `reverse_strong_variants`).

## 1.2.0

### Added

- `primer-finder design`: Primer3 on the regions of a finished run, and the assays it proposes ordered by a
  heuristic for which to try first — not a prediction of bench behaviour, and nothing here is wet-lab
  tested. An assay works either because no exclusion genome holds its amplicon or
  because an oligo sits on bases that differ; one that is neither would amplify both groups and is set
  aside. Primer3 is also asked for primers whose 3' end is forced onto a run of differing bases, and for
  probes straddling one, so that allele-specific assays are proposed as well as convenient ones.
- The assays are written as the primer files
  [insilicoPCR](https://github.com/duceppemo/insilicoPCR) reads, twice over — with the probe and without,
  since it cannot report both kinds at once — together with a script that runs it over both groups.
  `--insilico-pcr` runs it directly and records, per assay and per mode, whether it amplifies every
  inclusion genome and no exclusion genome. When the regions were found with `-p` below 1, an assay that
  reaches that share counts too and is reported as `partial (88%)`: the regions themselves were allowed to
  be missing from some inclusion genomes, so the assays on them are judged on the same terms.
- The copies of each amplicon in the inclusion genomes are counted and preferred, since a repeated target
  usually improves the limit of detection. Note that `find -d 1`, the default, discards repeated regions
  before the design step can see them.
- A logo, in `docs/images/`.

### Fixed

- A Primer3 request whose answer is hard to find no longer holds up the command. Each request is run on its
  own with a time limit, and pinning the probe over a position — the expensive constraint, and the one that
  ran for 45 minutes on a real region — gets a shorter one. A request that runs out of time is dropped and
  counted; its region is still designed on from the others.
- Nothing is forced where the rest of the assay would not fit, and a forced request only searches the
  neighbourhood of what it aims at rather than the whole region.
- `primer3_core --version` prints nothing, so `run_info.json` recorded a copyright line as its version;
  `--about` is used instead.

## 1.1.0

### Added

- `-p/--min-inclusion`: the fraction of the inclusion genomes a region has to be in, instead of all of them.
  The fraction is rounded up to a whole number of genomes, and both the kmer counting and the blast check
  use it. Below 1, each record in `all_inclusion_contigs.fasta` and `final_kmers.fasta` carries
  `inclusion=<holding>/<total>`, and the regions that cover the most inclusion genomes come first.
  The exclusion side is unchanged: a region must still be absent from every exclusion genome. This was the
  first of the two things the 2022 README listed as wanted.
- The package is published to PyPI as well as bioconda (`pip install primer-finder` gives the command, not
  the programs it runs).
- A validation suite (`validation/xylella/`) that recovers four published *Xylella fastidiosa* subspecies
  qPCR assays from public genomes, with a dated record of the run.

### Changed

- A file that is not a fasta is reported with the count and the first few names, instead of every path: a
  `-i` pointing at the wrong folder could print thousands of them.

## 1.0.0

First packaged release. The two scripts (`primer_finder.py` and `IDT_results_converter.py`) became the
`primer_finder` package with a `primer-finder` command, a test suite and documentation. The pipeline itself
is the same one: shared inclusion kmers (KMC), minus the exclusion kmers, assembled (SKESA or SPAdes),
mapped to one exclusion genome (minimap2) and checked against every genome (blast).

### Added

- A `primer-finder` command (`primer-finder -i inclusion/ -e exclusion/ -o results/`), installable from
  bioconda, plus `primer-finder idt` for IDT order sheets. `primer_finder.py` and
  `IDT_results_converter.py` still work.
- `run_info.json` with the parameters, the genomes, the version of every program used and the counts of
  each step, and `primer_finder.log` with the whole log, in the output folder.
- `--keep-intermediate` keeps the KMC databases, the SAM file and the blast databases; the run removes them
  by default.
- Gzipped genomes (`.fa.gz`, `.fasta.gz`, `.fna.gz`) work throughout: they are decompressed where blast
  needs plain fasta.
- A bundled example with simulated genomes and a known answer (`example/run_example.sh`), run by the CI.

### Changed

- No Python dependency at all: pysam, Biopython, pandas, psutil, pyahocorasick and xlrd are gone.
  Biopython's `NcbiblastnCommandline` wrapper, which the previous version used, no longer exists.
- bowtie2 and samtools are no longer needed. The contigs are mapped once with minimap2, whose SAM output is
  read directly, instead of being split into one file per contig and mapped one by one.
- The exclusion genome the contigs are mapped to is the first one in alphabetical order, not a random one,
  so that a run can be repeated. `-r` still forces a particular genome.
- The candidates are checked against each exclusion genome separately. A variant position has to differ in
  at least 90% of the exclusion genomes that the contig hits, as before, but the genomes are now counted
  instead of the blast alignments, and the inclusion and exclusion genomes are searched in parallel.
- The results are in subfolders (`1_kmers/`, `2_assembly/`, `3_candidates/`, `4_blast/`), with
  `final_kmers.fasta` at the top. `inclusion_blast_hits.tsv` holds 1 and 0 instead of `True` and `False`.
- Every hit of a contig in an exclusion genome is examined, not only the best one, and a position counts as
  different in that genome only if every copy of the region differs there: one matching copy is enough for
  an assay to amplify it. A genome whose hits do not reach a position no longer counts towards the 90%
  either, since it does not hold that region at all.
- SPAdes is called with `--only-assembler` (kmers carry no quality values to correct) and `--memory`, which
  the previous version did not pass.
- The blast databases are built in a folder of their own per genome, named after the genome's position in
  the list, and blast runs in that folder with relative paths. Two genomes with the same file name in
  different subfolders (`sampleA/contigs.fasta`, `sampleB/contigs.fasta`) therefore no longer collide, and
  the columns of `inclusion_blast_hits.tsv` are named after the path inside the input folder.
- An output folder that overlaps an input folder, or whose path contains a space (which BLAST cannot use in
  a database path), is refused before any work is done. The input folders may contain spaces.
- Python 3.10 or later (was 3.6).

### Fixed

- The variant positions in the headers of `final_kmers.fasta` are positions in the contig. They used to be
  positions in the blast alignment, which differ as soon as the hit does not start at the first base of the
  contig or holds a gap, so the lower-case bases could be the wrong ones.
- A contig that maps to the minus strand of the exclusion genome has its differences marked on the right
  bases: the cigar string describes the reverse complement, which the previous version applied to the
  forward sequence.
- Deletions no longer shift the lower-case marks: a deletion consumes no base of the contig.
- The rule "two differences within 21 bases" is applied to insertions and deletions as well as mismatches,
  and no longer misses the differences in the last three operations of the cigar string.
- `-k`/`--kmer_size` is checked (1-256): the previous check could not fail.
- `makeblastdb` no longer writes its index files into the inclusion and exclusion folders, which fails when
  they are read only.
- A genome listed in both groups, an output folder holding the input genomes, a file that is not a fasta and
  an unreadable kmer size are reported as errors instead of running anyway.
- A missing program is reported once, with the conda package that provides it, instead of failing midway.
- All the exclusion kmers are now subtracted. The old `-cx1e9` was read by KMC as `1`, so only kmers seen
  exactly once in the whole exclusion group were removed, and candidates in repeated exclusion regions
  survived. This changes the result of every run.
- With more than 255 inclusion genomes, `-ci <N>` could not match anything, because KMC counts to 255
  unless `-cs` says otherwise: such a run reported that no inclusion-specific kmer existed. The ceiling is
  now raised when the group needs it.
- Two differences that touch each other (a mismatch next to an inserted base, or two mismatches on either
  side of a deletion) count as two differences; the rule used to need a run of matching bases between them.
  A clipped end of at least 21 bases — a whole primer's worth of sequence the exclusion genome does not hold
  — now keeps the contig too.
- Two shared variant positions are enough to keep a contig, as documented; an off-by-three in the previous
  loop meant four were needed.
- `best_kmers.fasta` lists the contigs with the most differing bases first, so a contig the exclusion genome
  does not hold at all comes first instead of last.
- Fasta files whose headers are not UTF-8 (an accented strain name) are read instead of ending the run with
  a `UnicodeDecodeError` traceback; duplicate record names in an assembly are reported.
- An unreadable input, a read-only output folder or an `-o` that names an existing file is reported as an
  error instead of a traceback.
- The memory limit of a control group (Slurm, systemd, a container) is found wherever the job's cgroup sits,
  not only at the root of the cgroup filesystem, so the default `-m` is 85% of what the job may use rather
  than of the whole machine.
- `primer-finder idt` reads every format Excel saves: "CSV UTF-8" (with its byte-order mark), the legacy
  Windows CSV, "Unicode Text" (UTF-16) and the semicolon-separated CSV of comma-decimal locales. It no
  longer takes a styled but empty first row as the header, nor glues a phonetic run onto a cell value, and
  it says which sheet it read when a workbook holds several.
