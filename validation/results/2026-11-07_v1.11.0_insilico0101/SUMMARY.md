# 2026-11-07 — primer-finder 1.11.0 + insilicoPCR 0.10.1: nothing moves

insilicoPCR went from 0.8.0 to 0.10.1 in three releases (0.9.0 on 2026-11-03, 0.10.0 and 0.10.1 on
2026-11-06/07). This record asks whether primer-finder 1.11.0, written against 0.8.0, needs anything
changed to run on 0.10.1, and answers it twice: by reading what the releases changed, and by designing
the *Xylella* panel's **same find results** under both releases and joining the outputs.

**Answer: no change is needed.** On the panel, every file the in silico check writes or reads is
byte-identical between the two releases, apart from one region where Primer3 itself proposed a different
number of assays. Joined by oligo set, no in silico count or verdict moved on any of 4,084 assays.

primer-finder 1.11.0 (`a3e602a`), Primer3 2.6.1, blast 2.17.0+, insilicoPCR 0.10.1 (the `v0.10.1`
linux-x64 release, sha256-checked) against the 0.8.0 arm of
[`../2026-11-03_3prime_dependence/`](../2026-11-03_3prime_dependence/) (`assays_held`, primer-finder
`fdaf1c9`). The design code of 1.11.0 differs from `fdaf1c9` only by one rewrapped docstring line, so
the insilicoPCR release is the only thing that differs between the arms. Both releases bundle the same
blast 2.17.0+ and JDK.

## 1. What the releases changed

From insilicoPCR's changelog and `git diff v0.8.0 v0.10.1`:

- **0.9.0**: bundled BBTools 39.94 → 40.04 (a threaded-reader fix for multi-member `.fastq.gz`); the
  second FASTQ baiting pass runs with `hdist=0`.
- **0.10.0**: FASTQ mode calls amplicons the Tadpole mini-assembly missed from the reads themselves
  (rows marked `FASTQ read-level call`), writes `report_fastq_diagnostics.tsv`, and exposes the
  Tadpole settings; the GUI bait/assemble commands are shared with the command line.
- **0.10.1**: FASTQ-mode and GUI fixes, `--read-min-fraction`, longer `--evalue` help text.

All of it is FASTQ mode, the GUI or BBTools, which only FASTQ mode runs. primer-finder hands insilicoPCR
plain fasta genomes only (`insilico_genomes/`), so none of it runs on its behalf. Of what primer-finder
relies on, the diff leaves unchanged:

| what primer-finder relies on | 0.8.0 → 0.10.1 |
|---|---|
| `-i -p -o -t -m -l --probe-mismatches --3p-penalty` | unchanged; the new options are all FASTQ-mode |
| `-v` prints `insilicoPCR <x.y.z>` | unchanged; `check_version` compares the numbers, so 0.10.1 > 0.8.0 |
| `consolidated_report/report.tsv` columns and row rules | unchanged; read-level rows are appended only for FASTQ samples |
| `report(n).tsv` beside an existing report | unchanged |
| bundled blast 2.17.0+ | unchanged; only `runtime/common/bbmap` differs in the runtime |

## 2. The panel, designed under both releases

`run_panel.sh` designs the four targets with 1.11.0 + 0.10.1 on the find results that
`../2026-11-03_v1.10.0_insilico080/run_panel.sh` built; `compare_versions.py` joins them against
`assays_held`. Each design runs six insilicoPCR passes: PCR and qPCR files against the inclusion
group, the exclusion group, and the lenient exclusion pass behind `*_exclusion_held_by_3prime`.

| target | assays (0.8.0 / 0.10.1) | `assays.tsv` and both fastas | the six `report.tsv` | in silico columns moved |
|---|---|---|---|---|
| multiplex | 1149 / 1149 | byte-identical | byte-identical | 0 |
| fastidiosa | 1045 / 1045 | byte-identical | byte-identical | 0 |
| pauca | 901 / 904 | identical outside `Contig_2759_66.3578` | identical outside `Contig_2759_66.3578` | 0 |
| morus | 989 / 989 | byte-identical | byte-identical | 0 |

27 of the 36 files are byte-identical. The other nine are pauca's, and every difference in them sits on
`Contig_2759_66.3578`, where Primer3 proposed 27 assays this time against 24 in the baseline. That is
the three-proposal wobble the 3'-end record already noted on the same region (there it went 27 → 24),
a looseness of the candidate step that renumbers the region's assays; it is not insilicoPCR's. The 24
baseline assays on that region all reappear among the 27, with the same in silico results.

Joined by oligo set (forward, reverse, probe), 4,084 assays are shared, none is missing from the 0.10.1
arm, 3 are new to it (the wobble), and not one differs in any of the ten in silico columns: amplified
counts and percentages, `*_exclusion_held_by_3prime`, `pcr_selective`, `qpcr_selective`.

The four designs took 7 minutes side by side on 8 threads each, as under 0.8.0.

## What this does and does not show

- It shows primer-finder 1.11.0 runs on insilicoPCR 0.10.1 unchanged, and that the verdicts and the
  3'-end flag of [`../2026-11-03_3prime_dependence/`](../2026-11-03_3prime_dependence/) and
  [`../2026-11-03_v1.10.0_insilico080/`](../2026-11-03_v1.10.0_insilico080/) hold as measured under
  0.10.1. Neither record needs an erratum.
- It does **not** cover insilicoPCR's FASTQ mode, which primer-finder does not use, and says nothing of
  releases after 0.10.1. A later release that touches assembly-mode calling, the report or the options
  above should be run through `run_panel.sh` and `compare_versions.py` again.

## Reproducing

```bash
# the work folder of ../2026-11-03_v1.10.0_insilico080/ (its run_panel.sh builds it), then:
bash run_panel.sh                       # edit W, O and I0101 at the top
python compare_versions.py <work> <out>  # exits 1 if anything moved outside a Primer3 wobble
```
