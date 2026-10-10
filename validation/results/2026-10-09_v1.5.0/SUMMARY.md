# 2026-10-09 — primer-finder 1.5.0, why a probe's mismatches do not qualify an assay

1.4.0 added `--min-oligo-differences 2`: an assay resting on differences had to cover two of them under one
oligo, which is what [`find`](../../../docs/wiki/Methods.md) selects a region for. The argument behind the
number came from work on **primers**. This record is why it does not carry over to the probe, and what
restricting the filter to the primers costs.

No new sequencing or in silico run: the measurement is on the 1.4.0 run of
[`../2026-10-08_v1.4.0/SUMMARY.md`](../2026-10-08_v1.4.0/SUMMARY.md) — 25 *Xylella fastidiosa* genomes,
subsp. *multiplex* against the other four subspecies, 20 regions designed on, 511 assays proposed.

## The asymmetry

primer-finder asks Primer3 for:

| | primer | probe |
|---|---|---|
| length | 18–25 nt, 20 optimal | **18–27 nt, 22 optimal** |
| melting temperature | 58–63 °C | **62–72 °C** |

The probe reaches its higher melting temperature mostly by being longer, and the reaction anneals near
60 °C. It therefore starts with considerably more binding energy in hand than a primer does, and two
mismatches spread over a longer duplex take away a smaller share of it.

## What the literature says

- **A conventional TaqMan probe still gave a detectable signal through five mismatches**, and under standard
  conditions neither it nor a 3'-MGB probe was sequence-specific. The authors recommended MGB chemistry for
  discrimination and said reaction conditions still had to be optimised and specificity verified even then.
  Yao Y, Nellåker C, Karlsson H (2006) *Mol Cell Probes* 20(5):311–316,
  [doi:10.1016/j.mcp.2006.03.003](https://doi.org/10.1016/j.mcp.2006.03.003).
- **A probe discriminates by being short, not by carrying more mismatches.** A 12-base MGB probe melts at
  the same temperature (65 °C) as an unmodified 27-base one, and shortening the duplex is what makes a
  single mismatch a large fraction of its stability — the basis of short MGB/LNA probes with the difference
  near the centre for SNP genotyping. Kutyavin IV et al. (2000) *Nucleic Acids Res* 28(2):655–661,
  [PMC102528](https://pmc.ncbi.nlm.nih.gov/articles/PMC102528).
- **In a 5'-nuclease assay the discrimination sits in the primers.** A single mismatch in a primer's 3'
  region shifts the quantification cycle by under 1.5 to over 7 cycles depending on which base pair it is,
  with up to sevenfold variation between master mixes. Stadhouders R et al. (2010) *J Mol Diagn*
  12(1):109–117, [PMC2797725](https://pmc.ncbi.nlm.nih.gov/articles/PMC2797725). Lefever S et al. (2013)
  *Clin Chem* 59(10):1470–1480 is the companion survey of 3'-end mismatches.
- The rule for primers stands on its own evidence: a single mismatch is often not enough, which is the
  stated reason double-mismatch allele-specific designs exist. Lefever S et al. (2019) *Sci Rep* 9:2150,
  [doi:10.1038/s41598-019-38581-z](https://doi.org/10.1038/s41598-019-38581-z).

## What it costs

Of the 504 assays that could tell the groups apart, 309 rest on a difference:

| | assays resting on a difference |
|---|---|
| two or more under **any one oligo**, probe included (1.4.0) | 241 |
| two or more under **one primer** (1.5.0) | **168** |
| moved to the set-aside list | 73 |
| …of which the primers carry **no** difference at all | 40 |
| …of which the primers carry exactly one | 33 |

Carried forward in all: 511 proposed → 363 (195 by absence, 168 by difference), 141 set aside, 7 that would
amplify both groups. **All 20 regions are still represented.**

## What this does and does not settle

- It does **not** show that a probe with two mismatches always binds. Yao et al. is one study, one target,
  2006, under standard conditions, and "a detectable signal" is not "indistinguishable": a probe mismatch
  that costs several cycles may still discriminate against a validated threshold.
- Everything here is condition-dependent. Stadhouders' sevenfold spread between master mixes is the finding,
  not a footnote.
- **In silico PCR cannot arbitrate**, in either direction: it decides binding by counting mismatches, so it
  will always report a two-mismatch probe as not binding — the very assumption at issue. In the 1.3.0 sweep
  all 58 probe-only assays were reported selective in qPCR mode and none in PCR mode; that gap was the model
  trusting the probe, not evidence about it.
- So the change is a refusal to treat a probe's mismatches as qualification, not a claim that such a probe
  never works. `--min-primer-differences 1` accepts those assays for anyone who would rather judge them.
- Nothing here has been tested at a bench.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/assays -t 32 --max-regions 20
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/loose  -t 32 \
    --max-regions 20 --min-primer-differences 1      # what 1.3.0 carried forward
```

The two counts of the cost table are `best_primer_variants` against the maximum over all three oligos, read
off `assays.tsv` of the second run.
