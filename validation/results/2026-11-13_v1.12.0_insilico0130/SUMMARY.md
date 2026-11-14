# 2026-11-13 — insilicoPCR 0.12.0 to 0.14.0: one mismatch near an oligo's end counts as one

insilicoPCR 0.12.0 (2026-11-12) changes how a mismatch near the end of a primer or a probe is counted,
which is what primer-finder's in silico verdicts rest on; 0.12.1 and 0.13.0 (2026-11-13) followed within
a day. This record asks what the three releases change for primer-finder, answers it by reading the
diffs, measuring the hit filters one base at a time, and designing the *Xylella* panel and two smaller
sets under 0.11.0 and the newer releases, and says what primer-finder changed in response.

**Answer: primer-finder now requires insilicoPCR 0.12.0**, and counts the 3'-end check on its original
terms. The strict primer rules are unchanged and no PCR verdict moves. 23 qPCR verdicts move at the
defaults, all on a probe mismatched near one of its ends, in the direction the release intends. The
3'-end check now sees a mismatch two to four bases from a primer's 3' end, which no earlier release could
readmit. 0.12.1 and 0.13.0 give the same results as 0.12.0, byte for byte, on everything here, and so
does 0.14.0 (2026-11-13, section 6).

primer-finder 1.12.0 (`v1.12.0`, `3d427d3`) for every comparison between releases, Primer3 2.6.1, blast
2.17.0+ (bundled by every release), insilicoPCR 0.11.0, 0.12.0, 0.13.0 and 0.14.0 (the linux-x64 releases,
sha256-checked), every arm in the same `primer-finder_dev` environment. The end of section 5 runs the
code this record was committed with.

## 1. What the releases changed

From insilicoPCR's changelog and `git diff v0.11.0 v0.13.0`:

- **0.12.0, hit filters.** blast (`blastn-short`, +1/−3) trims a mismatch followed by three or fewer
  matching bases off the end of its local alignment, so one mismatch near an oligo's end reads as several
  unaligned bases. Until 0.12.0 insilicoPCR judged the hit by what was left. It now reads the trimmed
  bases back from the genome and counts the ones that do not pair:
  - *primers, 3' end*: one mismatch 1-3 bases from the terminus counted once per trimmed base, in the
    `-m` budget and in the weighted 3' penalty; now once. At `--3p-penalty 0` (every strict run
    primer-finder makes) this changes nothing, since the last five bases have to be perfect either way.
  - *probes, both ends*: a mismatch in the first or last two bases was trimmed and not counted (`0`
    probe mismatches); one three or four bases in left 3-4 bases unaligned, more than the 2 a probe may
    have, and was refused at any `--probe-mismatches`. Now every probe base the genome has is compared.
  - *primers, 5' end*: unchanged (`--5p-truncation` still counts the bases blast left out, paired or not).
- **0.12.0, report.** Two columns, `Forward5PrimeUnaligned` and `Reverse5PrimeUnaligned`, before
  `SizeCheck`: the 5' bases of each primer that do not pair. primer-finder does not read them.
- **0.12.0, command line.** An input mistake ends the run with an `Error:` line instead of a stack trace;
  exit code 1 as before.
- **0.12.1.** Gzipped assemblies are decompressed before blast (primer-finder hands insilicoPCR plain
  fasta); failing tools are reported as messages.
- **0.13.0.** FASTQ mode extends reads to find products longer than a read (primer-finder does not use
  FASTQ mode). Every input file is checked first, and a sample with a faulty file is left out with a
  warning, exit 0 (until now an empty file was left out without a word). A memory fix for very large
  contigs, introduced by 0.12.0. The hit filters and the report columns are not touched: the diff of
  `BlastHitParser`, `TrimmedEnds` and `ConsolidatedReport` is comments, a type and a row order used in
  FASTQ mode.

What that means for what primer-finder relies on:

| what primer-finder relies on | 0.11.0 → 0.13.0 |
|---|---|
| `-i -p -o -t -m -l --probe-mismatches --3p-penalty`, `-v` | unchanged; the new options are FASTQ mode and the input check |
| `report.tsv`: `Sample`, `Gene`, the mismatch and 3'-window columns, `SizeCheck` last | unchanged; two columns added before `SizeCheck` |
| which primer hits the strict runs accept | unchanged (sections 2-4) |
| which probe hits are accepted | **changed**: every probe base counted |
| what `--3p-penalty 4` readmits in the lenient run | **changed**: a mismatch 2-4 bases from the 3' end, and two in one window |

## 2. The hit filters, measured

`mismatch_zones.py` is the 0.8.0 record's synthetic-template test with three sections added: a 22-base
primer pair (and a 24-base probe) in a template mutated at known distances from an oligo's end, one
genome per case. Sections 1-4 under 0.11.0 and 0.12.0 differ only here (`yes` = amplified):

| one mismatch, bases from the 3' end, at `-m 4` | `--3p-penalty` 0 / 1 / 2 / 3 / 4 under 0.11.0 | under 0.12.0 |
|---|---|---|
| 1 (terminal) | no / no / no / no / yes | no / no / no / no / yes |
| 2 | no at every penalty | no / no / yes / yes / yes |
| 3 | no at every penalty | no / no / no / yes / yes |
| 4 | no at every penalty | no / yes / yes / yes / yes |
| 5 | no / yes / yes / yes / yes | no / yes / yes / yes / yes |

0.12.0 prices each position as insilicoPCR's own table does: 2 for the second or third base (3 when the
primer's base is a G or a C, as the third is here), 1 for the fourth. Runs of two or more mismatches
ending at the terminus are refused by both releases at every setting. At `--3p-penalty 0`, and for the
5' end, the two releases give identical tables.

Section 5 is the 3'-end check of a design, at `--mismatches 1`: the strict pass (`-m 1`), the lenient
pass (`-m 2 --probe-mismatches 1 --3p-penalty 4`), and what `read_three_prime_held` of primer-finder 1.12.0
counts as held:

| mismatches in the template, bases from the 3' end | held under 0.11.0 | held under 0.12.0 |
|---|---|---|
| 1, or 5 | yes | yes |
| 2, 3 or 4 | no | **yes** |
| 6, 7, 8 | no (the strict pass amplifies: a real cross-reaction) | no |
| 1+5, 2+4, 2+5, 3+4, 3+5, 4+5 (two in one window) | no | **yes** |
| 3+10 (one in the window, one internal) | no | **yes** |
| 1+10, 5+10 | yes | yes |
| 1+2, 3+4+5 | no | no |
| 8+12 (both internal) | no | no |

Sections 6 and 7 are the probe, at `-m 1`:

| probe mismatches | `--probe-mismatches` 0 / 1 / 2 under 0.11.0 | under 0.12.0 |
|---|---|---|
| one, 1st or 2nd base from either end | yes / yes / yes, reported as `0` | no / yes / yes, reported as `1` |
| one, 3rd or 4th base from either end | no / no / no | no / yes / yes |
| one, 5th base in or further | no / yes / yes | no / yes / yes |
| two, one of them 1-2 bases from an end | no / yes / yes, reported as `1` | no / no / yes |
| two, one of them 3 bases from an end | no / no / no | no / no / yes |
| two, one at each end's first base | yes / yes / yes, reported as `0` | no / no / yes |
| two, both inside | no / no / yes | no / no / yes |

So under 0.11.0 a probe's budget meant little near its ends: up to two mismatches were free there, and a
single one three or four bases in could not be tolerated at all. 0.12.0's tables are what
`--probe-mismatches` says. 0.13.0's output is identical to 0.12.0's, all seven sections.

## 3. The example and the contig-edge set

The example (17 assays, eight synthetic genomes): designed under 0.11.0 and 0.12.0
(`../2026-11-10_v1.11.0_insilico0110/compare_example.sh`), `assays.tsv` and both primer fastas are
byte-identical and every report row is the same once the two new columns are set aside; it holds no
probe or 3'-window case that the new counting decides differently. Under 0.13.0 every file is
byte-identical to 0.12.0's, through `design` and through the `run_insilico_pcr.sh` it writes.

The contig-edge set of `../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py` (300 assays on 600 kb of
*L. monocytogenes*, broken inside 100 primer sites, a copy with one base in 45 substituted, and *L.
innocua*), run by `compare_edges.py` the way a design runs it, with sizes declared:

| | 0.11.0 | 0.12.0 | |
|---|---|---|---|
| PCR, strict: amplified (assay, sample) pairs | 560 | 560 | the same pairs; 586 products, the same size and mismatch counts |
| PCR, held only by the 3' end | 37 | 99 | all 37 kept; 81 lenient products more, each with a 3'-window imperfection |
| qPCR, strict: amplified pairs | 516 | 529 | 512 in both; 17 gained, each with a probe mismatch the strict budget of 1 allows; 4 lost, whose probes 0.11.0 had aligned over 22-23 of 24 bases |
| qPCR, held only by the 3' end | 30 | 73 | 29 kept |

On 12 qPCR products both releases report, the probe's place and count move from `0` mismatches over a
shortened alignment to `1` over the whole probe. Nothing moved that the release does not announce (the
script checks it). 0.13.0's four reports are byte-identical to 0.12.0's.

## 4. The *Xylella* panel

`run_panel.sh` rebuilds the panel (the 25 RefSeq genomes, the four subspecies' find results) and designs
it under 0.11.0 and 0.13.0 at three settings; `compare_panel.py` joins the arms by oligo set. 2,430
assays are shared at every setting (Primer3 proposed 27 assays on pauca's `Contig_2759_66.3578` in one
run and 24 in the other, the wobble the earlier records describe). The default round was also run under
0.12.0: every file is byte-identical to 0.13.0's outside that region.

At the defaults (`--mismatches 1`, the probe budget following it):

| | PCR | qPCR |
|---|---|---|
| verdicts moved | none | 16 `yes` → `no`, 7 `no` → `yes` |
| `exclusion_amplified` | the same on every assay | up on 17, down on 7 |
| `yes` assays flagged `held_by_3prime`, 0.11.0 | 332 of 2,047 (16%) | 229 of 2,133 (11%) |
| the same, 0.13.0 read by primer-finder 1.12.0 | 570 of 2,047 (28%) | 362 of 2,124 (17%) |

Every qPCR verdict that moved, moved on its probe. The 16 `yes` → `no` gain an exclusion genome whose
probe now binds through one mismatch, refused by 0.11.0 because it sat three or four bases from an end.
The 7 `no` → `yes` lose 26 (assay, exclusion genome) pairs that 0.11.0 had counted with one mismatch on a
probe aligned over 19 or 21 bases: the trimmed end held a second mismatch, over the default budget.

The other two settings, as the wiki argues from them:

| against the default | 0.11.0 | 0.13.0 |
|---|---|---|
| `--probe-mismatches 2`, qPCR verdicts | 13 `yes` → `no` | 28 `yes` → `no` |
| `--probe-mismatches 2`, PCR verdicts | none | none |
| `--mismatches 0`, PCR verdicts | 151 `no` → `yes`, 5 `yes` → `no` | the same |
| `--mismatches 0`, qPCR verdicts | 123 `no` → `yes`, 6 `yes` → `no` | 141 `no` → `yes`, 7 `yes` → `no` |

The 0.11.0 column reproduces the 0.8.0 record's 13 exactly. `--probe-mismatches 2` is still a stricter
screen and never a rescue, now twice as strict; the argument for `--mismatches 1` stands as written.

## 5. The 3'-end check: what it counts, and what primer-finder changed

The lenient run passes `--3p-penalty 4`, the price of a mismatched G/C terminal base, so that one
imperfection anywhere in a primer's last five bases is readmitted. Under 0.12.0 it also readmits two in
one window when they come to 4 or less (the fourth and fifth bases score 1 each, the second or third 2),
which 0.11.0's trimming never let through. The held (assay, exclusion genome) pairs of the default round,
by the 3'-window imperfections of the lightest product:

| | one | one in each primer | two or more in one primer | total |
|---|---|---|---|---|
| PCR, 0.11.0 | 1,868 | 140 | 0 | 2,008 |
| PCR, 0.13.0 | 2,651 | 720 | 818 | 4,189 |
| qPCR, 0.11.0 | 1,164 | 63 | 0 | 1,227 |
| qPCR, 0.13.0 | 1,380 | 209 | 298 | 1,887 |

The column documents and ranks by "an imperfection at the far 3' end, which a polymerase occasionally
reads through". Two mismatches in the last five bases of one primer are not that, and counting them
would have changed what the column means without anything saying so. So `read_three_prime_held` now
leaves out a row with more than one imperfection in a primer's window (`MOST_IN_WINDOW`); one in each
primer is counted, as it always was.

The code this record was committed with, designing the panel's default round under 0.13.0, against
primer-finder 1.12.0 under 0.13.0 on the same find results (2,433 assays, the same oligos throughout):
every one of the 24 reports is byte-identical, no verdict and no column other than
`*_exclusion_held_by_3prime` differs, and that one only goes down (115 assays in PCR mode, 49 in qPCR).

| | PCR | qPCR |
|---|---|---|
| held (assay, exclusion genome) pairs, 1.12.0 → this code | 4,189 → 3,303 | 1,887 → 1,584 |
| `yes` assays flagged, 0.11.0 | 332 of 2,047 (16%) | 229 of 2,133 (11%) |
| `yes` assays flagged, 0.13.0 with this code | **501 of 2,047 (24%)** | **330 of 2,124 (16%)** |

So the check now flags one selective assay in four in PCR mode and one in six in qPCR mode, where 0.11.0
flagged one in six and one in nine: the difference is the mismatches two to four bases from a 3' end
that 0.11.0's trimming hid. The ranking breaks ties on the count, so the order of the table moves (543 of
675 multiplex assays sit at another rank than under 1.12.0, morus 257 of 471), while the first assay of
every target and 79 of the first 80 stay.

## 6. insilicoPCR 0.14.0

0.14.0 followed the same day. Its changelog and `git diff v0.13.0 v0.14.0`:

- The application now runs the same pipeline as the command line instead of a copy of it.
- In FASTQ mode, the samples processed at once share the threads and the memory. A sample that runs
  out of memory is run again, alone, at the end.
- A failed BLAST+ or BBTools command is reported as `<program> failed with exit code <n>`.
- The gel code was split into classes, and tests were added.

`BlastHitParser`, `TrimmedEnds`, `ConsolidatedReport` and `CliConfig` are not touched. `Dispatcher`
only returns its exit code instead of calling `System.exit` itself, and `-v` still prints
`insilicoPCR 0.14.0`. The installed release is the published zip (sha256-checked, the same files).
Measured against 0.13.0, both read by the code this record was committed with:

- **Mismatch zones:** `mismatch_zones.py` gives the same output, and all 303 reports are byte-identical.
- **Example:** every file is byte-identical, through `design` and through `run_insilico_pcr.sh`.
- **Contig-edge set:** the same 560 PCR and 529 qPCR pairs, the same 93 and 70 held only by the 3' end,
  and the same products, sizes, mismatch counts and probe places.
- **Panel:** at `--probe-mismatches 2` and `--mismatches 0`, the 36 reports, `assays.tsv` and primer
  fastas of each round are byte-identical. At the default, the 0.13.0 arm is the round of section 5. Its
  files are byte-identical except pauca's, where Primer3 again proposed 27 assays on one contig against
  24. Joined by oligo set (2,430 assays at every setting), no verdict, amplified count, held count or
  rank moves.

So 0.14.0 changes nothing primer-finder reads.

## What this does and does not show

- It shows primer-finder needs a change to run on 0.12.0 and newer: the qPCR verdicts and the 3'-end
  check count other things than under 0.11.0, and mixing reports of the two would not show anywhere. So
  `design` refuses an older release before Primer3, `run_insilico_pcr.sh` asks for 0.12.0, and a report
  without the two columns 0.12.0 added is refused. primer-finder 1.12.0 remains the release for 0.11.0.
- It shows the strict primer rules, and so every PCR verdict, are unchanged on 2,430 assays at three
  settings, and on the contig-edge set.
- It shows every qPCR verdict that moved, moved on the probe rule the release announces, in its
  direction.
- It does not show which counting a polymerase agrees with. 0.12.0 makes the model say what its options
  say; whether a probe mismatched three bases from its end still reports at the bench is a question for
  the bench.
- It does not cover FASTQ mode, which primer-finder does not use.

## Reproducing

```bash
# kmc, skesa, minimap2, blast and primer3_core on PATH (the primer-finder_dev environment), the NCBI
# datasets tool for the downloads, and the releases extracted. Every script needs primer-finder 1.12.0
# (git worktree add <dir> v1.12.0, with this folder copied in), which still runs 0.11.0.
python mismatch_zones.py /path/to/insilicoPCR-0.11.0 work_0.11.0 > zones_0.11.0.txt
python mismatch_zones.py /path/to/insilicoPCR-0.13.0 work_0.13.0 > zones_0.13.0.txt

bash ../2026-11-10_v1.11.0_insilico0110/compare_example.sh <out>/example /path/to/0.11.0 /path/to/0.13.0 8

python ../2026-11-10_v1.11.0_insilico0110/make_contig_edges.py lmono_EGDe.fasta linnocua_Clip11262.fasta <out>/edges
python compare_edges.py <out>/edges /path/to/insilicoPCR-0.11.0 /path/to/insilicoPCR-0.13.0

bash run_panel.sh <work> <out>/panel /path/to/insilicoPCR-0.11.0 /path/to/insilicoPCR-0.13.0 8
python compare_panel.py <out>/panel

# Section 6, from this commit (which refuses 0.11.0): the same scripts, with 0.13.0 and 0.14.0.
python mismatch_zones.py /path/to/insilicoPCR-0.14.0 work_0.14.0 > zones_0.14.0.txt
bash run_panel.sh <work> <out>/panel14 /path/to/insilicoPCR-0.13.0 /path/to/insilicoPCR-0.14.0 8
```

A panel round takes about 18 GB, mostly `amplicons/`, so a full `run_panel.sh` needs about 110 GB.

The contig-edge set took about 6 minutes for both releases, and each round of the panel (eight designs
side by side, 8 threads each) about 10 minutes on 64 cores.
