# 2026-11-03 — how much of a `yes` rests on the far 3' end, measured on the panel

The 3'-end dependence check (unreleased, `fdaf1c9`, post-1.10.0) reruns the exclusion group tolerating
an imperfection at the far 3' end and counts, per assay, the exclusion genomes that amplify only there —
`*_exclusion_held_by_3prime`. [`../2026-11-03_v1.10.0_insilico080/`](../2026-11-03_v1.10.0_insilico080/)
measured what the lenient settings readmit; this record asks what the column actually finds on the same
*Xylella* panel: the four targets' designs rerun with the check on the **same find results** as that
record, so everything joins against its three arms. `analyze_held.py` in this folder prints the numbers.

primer-finder `fdaf1c9`, insilicoPCR 0.8.0, Primer3 2.6.1, blast 2.17.0+, defaults throughout.

## The counts

Of the panel's selective assays, the share held off at least one exclusion genome **only** by the far
3' end of a primer:

| mode | flagged / selective | held 1 genome | 2 | 3–9 | 10 or more |
|---|---|---|---|---|---|
| PCR | **332 / 2047** (16%) | 114 | 46 | 112 | 60 |
| qPCR | **229 / 2133** (11%) | 85 | 32 | 86 | 26 |

Per target (PCR / qPCR): multiplex 98 / 78, fastidiosa 107 / 81, pauca 31 / 12, morus 96 / 58.

The verdicts themselves did not move: joined by oligo sequence, every assay shared with the 1.10.0
record's `assays_new` arm carries the same `pcr_selective` and `qpcr_selective` — the annotation is an
annotation. (Joined by *name* ten verdicts appear to move, all on `pauca/Contig_2759_66.3578`, where
this run's Primer3 proposed 24 assays against the earlier 904-run's 27 and the numbering shifted; the
24 shared oligo sets agree exactly. The three-proposal wobble on one region of 199 is a pre-existing
looseness of the candidate step, not of this check, and is noted here so the name-join is not trusted
on that region.)

## Against the retired `undecided (3' end)` verdict

The old verdict parked 34 PCR and 12 qPCR assays on this panel, and the 1.10.0 record showed all of
them resolving to `yes`. The flag catches **17 of the 34** and **8 of the 12** — and the ones it does
not catch are not missed, they are refused on the model's own terms: every one of their old exclusion
bindings was a **two-base** 3' trim (`EndMismatch = -2` in the 0.6.1 reports — e.g.
`multiplex/Contig_551_74.3836_assay6`, `Contig_2971_75.9495_assay18`), a run the lenient pass still
refuses because its window prices it 3+2, past the penalty of 4. A single readable base flags; a run of
two does not. That is the distinction the old model could not draw and the reason the verdict was
retired rather than translated.

## What the column adds beyond the retired verdict

The flag marks **315 PCR and 221 qPCR `yes` assays the old model never parked** — ten times the old
verdict's reach. Of the 315, the old pair called 142 `yes` with nothing attached and 173 `no`. Neither
was informative about this risk, and for the same reason: the held bindings sit at sites the old engine
did not report at all — refused by its fixed two-base 5' allowance, or never paired as a plausible
product — so their one-readable-base-from-amplifying 3' ends were invisible however the verdict fell.
The old verdict could only flag an assay when the terminal binding was the *only* exclusion evidence
left after the engine's own blind spots; the measured rerun sees the sites first and asks the question
second.

## What this does and does not show

- It shows the column finds real weight: one selective assay in six (PCR mode) on this panel is one
  readable base away from amplifying an exclusion genome, and the ranking now sinks those below equally
  backed assays — including 142 that previously ranked as unqualified `yes`.
- It shows the flag is strictly narrower than the verdict it replaces, and why: runs of two terminal
  mismatches stay refused at the lenient settings, which is the model taking the orthodox position
  twice rather than abstaining once.
- It does **not** show a flagged assay would cross-react at a bench, or an unflagged one would not. The
  flag is the model's own lenient reading of its own refusals; the polymerase was not consulted.

## Reproducing

```bash
# the same work folder as ../2026-11-03_v1.10.0_insilico080/ (its run_panel.sh builds it), then:
for t in multiplex fastidiosa pauca morus; do
  primer-finder design <work>/$t/results -o <work>/$t/assays_held -t 8 \
      --insilico-pcr /path/to/insilicoPCR-linux-x64
done
python analyze_held.py <work>
```

The four designs took 7 minutes side by side on 8 threads each.
