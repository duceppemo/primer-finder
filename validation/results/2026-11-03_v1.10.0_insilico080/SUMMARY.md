# 2026-11-03 — primer-finder 1.10.0 + insilicoPCR 0.8.0, the verdicts remeasured

Every in silico PCR number in the earlier records was produced by insilicoPCR 0.6.1, whose last-two-bases
blind spot shaped primer-finder's `undecided (3' end)` verdict and whose raw row pairings shaped 1.9.2's
re-checking. insilicoPCR 0.8.0 replaced both with a deliberate model, and 1.10.0 trusts it. This record
measures what that trust buys, twice over: the new engine's hit-filtering zones on a controlled synthetic
template, and the whole *Xylella* panel scored by both tool pairs **on the same find results and the same
proposed assays**, so every verdict joins by assay name.

primer-finder 1.10.0 (and 1.9.3 for the paired design), insilicoPCR 0.8.0 and 0.6.1 (portable releases),
libprimer3 2.6.1, blast 2.17.0+. Panel: the 25 complete RefSeq genomes of `validation/xylella/`, four
targets (multiplex, fastidiosa, pauca, morus), defaults throughout (`-m 1`). [`run_panel.sh`](run_panel.sh)
is the whole run; one `find` per target, then three designs on its results: 1.10.0 + 0.8.0, 1.9.3 + 0.6.1,
and 1.10.0 + 0.8.0 with `--probe-mismatches 2`.

## 1. The zones of 0.8.0, measured ([`mismatch_zones.py`](mismatch_zones.py))

A 22-mer pair in a synthetic template, the template mutated one base at a time at a known offset from a
primer end; each case its own genome. The 0.6.1 predecessor of this script
([`../2026-10-08_v1.3.0/mismatch_zones.py`](../2026-10-08_v1.3.0/mismatch_zones.py)) found three accidental
zones; 0.8.0's are deliberate, and two of them sit somewhere other than where its own documentation's
penalty table suggests:

| where the difference is | at the defaults (`--3p-penalty 0`) | rescuable? |
|---|---|---|
| the 3'-terminal base | refused at every `-m` (0, 1, 2, 3, 10 tested) | yes: `--3p-penalty` 3 for an A/T-ending primer, 4 for G/C-ending — measured both — plus one unit of `-m`, which the trimmed base consumes |
| 2–4 bases from the 3' end | refused at every `-m` | **no, not even at `--3p-penalty 4` with `-m 4`.** blast trims the whole tail rather than aligning through, so each unaligned base scores in the window and the penalties stack past 4 |
| 5 bases from the 3' end | refused at every `-m` | yes: `--3p-penalty 1` |
| 6 or more bases in | what `-m` decides, one unit each | — |
| a run of 2–4 ending at the 3' end | never, at any tested `-m` | not at `--3p-penalty` ≤ 4 |
| 1–4 bases from the **5' end** | **free**: amplifies at `-m 0` | `--5p-truncation 2` narrows it to the last 2, measured |
| 5 bases from the 5' end | costs one unit of `-m` (blast keeps it internal rather than trimming 5) | — |

So the 3'-end zone the ranking aims at — a difference on the very last base — is now judged and judged
against the exclusion binding, which is what retired the `undecided (3' end)` verdict; and the 5' end is
*more* forgiving than 0.6.1's two-base trim, which matters below.

## 2. The panel: 2,433 assays scored by both tool pairs ([`compare_verdicts.py`](compare_verdicts.py))

Identical assays, identical genomes, defaults on both sides. The verdict cross-tables, `partial` folded
into its neighbours (none occurred at `-p 1`):

```
pcr_selective        1.10.0+0.8.0:  yes       no    total
1.9.3+0.6.1  yes               1801       56     1857
             undecided (3'end)   34        0       34
             no                 212      330      542
             total             2047      386     2433

qpcr_selective       1.10.0+0.8.0:  yes       no    total
1.9.3+0.6.1  yes               1974       32     2006
             undecided (3'end)   12        0       12
             no                 147      268      415
             total             2133      300     2433
```

- **Every formerly-`undecided (3' end)` assay — all 34 in PCR mode, all 12 in qPCR mode — is now a clean
  `yes`**: full inclusion coverage, zero exclusion amplification. The question that label existed to park
  ("the model cannot see the difference this assay rests on") is answered in the assays' favour, which is
  what the ranking always bet by placing them above the refused.
- **212 (PCR) and 147 (qPCR) former `no`s are now `yes`.** The old `no`s rested on exclusion bindings the
  new engine refuses: a primer whose last five bases do not match has nothing for the polymerase to extend.
- **56 (PCR) and 32 (qPCR) went the other way.** Of the 56, **53 gained an exclusion amplification**: the
  more forgiving 5' end (four free bases where 0.6.1's alignment reported nothing past a two-base trim)
  and indel-tolerant binding report sites the old engine never saw — including distant off-target pairings
  within `-l 10000`, e.g. an 8,400-base product off a one-mismatch reverse-primer site in
  `GCF_004016415.2` (sandyi) that flips `multiplex/Contig_183_67.6843_assay12`. A selectivity that rested
  on a 5'-proximal difference was never selectivity a bench could trust, so this is the cautious direction.
  The other **3 lost an inclusion genome**: an inclusion variant sits in a primer's 3' window, which the
  old engine could not see either.
- Per target, `pcr yes` / `qpcr yes` of the usable assays:

  | target | usable | 1.9.3 + 0.6.1 | 1.10.0 + 0.8.0 | design wall-clock, old → new |
  |---|---|---|---|---|
  | multiplex | 675 | 505 / 554 (undecided 9 / 2) | 565 / 598 | 498 s → 205 s |
  | fastidiosa | 529 | 351 / 390 (undecided 25 / 9) | 422 / 441 | 656 s → 162 s |
  | pauca | 758 | 655 / 672 (undecided 0 / 0) | 660 / 676 | 391 s → 199 s |
  | morus | 471 | 346 / 390 (undecided 0 / 1) | 400 / 418 | 625 s → 252 s |

  The 2–3× speedup is the retired re-checking: 1.9.3 read the detailed reports row by row to tell
  amplicons from phantom pairings; 0.8.0 never writes the phantoms.

- **The counts of the 1.3.0-era records are not reproduced here and were not expected to be**: those were
  465 assays from 20 regions under a 1.3.0 pipeline; this is 2,433 from the same panel under 1.10.0
  defaults. The 96-`undecided` figure quoted in older pages belongs to that run; this panel's answer to
  the same question is the 34 + 12 above, and all of them resolved to `yes`.

## 3. `--probe-mismatches 2` against the default

The suggestion in insilicoPCR's documentation — a probe hybridizes hotter, so give it `-m`+1 — reads as an
inclusivity rescue. On this panel it is the opposite: inclusion coverage was already complete (the assays
are designed on sequence every inclusion genome shares), so the looser budget had nothing to rescue and
could only find more exclusion probe bindings. **13 qPCR verdicts moved, every one `yes` → `no`**, and
none moved the other way. So the default (`--probe-mismatches` unset, the `-m` value) stands as
primer-finder's default, and `--probe-mismatches 2` is useful as a *stricter* second screen: an assay
whose `yes` survives it does not depend on the probe refusing a near-match.

## What this changes

- `docs/wiki/Designing-assays.md`'s zone table now carries the measured zones, including the two that
  differ from insilicoPCR's own penalty arithmetic (2–4 from the 3' end unrescuable; 5 free bases at the
  5' end measured as 4 plus one `-m`).
- The `undecided (3' end)` retirement in 1.10.0 is vindicated by outcome, not only by argument: on this
  panel, nothing the label would have parked needed parking.
- Still not evidence about PCR: nothing here was near a bench. What is measured is the model's own rule,
  both engines' on the same assays.
