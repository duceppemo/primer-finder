# 2026-11-03 — primer-finder 1.10.0 + insilicoPCR 0.8.0, what the new amplicon calling changes

insilicoPCR 0.8.0 replaced the accidental zones that
[`../2026-10-08_v1.3.0/SUMMARY.md`](../2026-10-08_v1.3.0/SUMMARY.md) measured in 0.6.1 — the two free
terminal bases, the 3–4-base dead zone — with a deliberate model: `-m` counts mismatches, indels **and
3' bases blast left unaligned**; a weighted check guards the last `--3p-window` bases (default 5, maximum
penalty `--3p-penalty`, default 0); and up to `--5p-truncation` unaligned 5' bases are free (default 5,
was a fixed 2). It also calls amplicons itself — convergent orientation, product-size cap — where 1.9.2+
had to reconstruct that from raw hit pairs. primer-finder 1.10.0 speaks that interface.

This record asks two questions:

1. does the new model do what it says? Each claim is measured the way the 1.3.0 record measured 0.6.1,
   rather than quoted (section 1);
2. what do the verdicts of a real panel do when the **same find results** are designed under
   1.9.3 + 0.6.1 and under 1.10.0 + 0.8.0 (section 2)? The question of the release is where the old
   `undecided (3' end)` verdicts land, since the blind spot that forced that verdict is gone. A third
   design measures `--probe-mismatches 2` against the 1.10.0 defaults.

The 0.6.1 measurements in the 1.3.0 record remain true of 0.6.1; nothing there is contradicted. This
record is about what replaced it.

primer-finder 1.10.0 (`v1.10.0`) against a `v1.9.3` checkout, Primer3 2.6.1, blast 2.17.0+,
insilicoPCR 0.8.0 and 0.6.1, each with its bundled blast 2.17.0+. Both design arms ran inside the same
`primer-finder_dev` environment, so only the two versions differ.

## 1. What the 0.8.0 knobs govern, measured

Same instrument as the 1.3.0 record, re-aimed: a 22-mer primer pair in a synthetic template, the
*template* mutated — one mismatch at a known distance from the left primer's 3' end (`offN`, 1 = the
terminal base), or a run ending at the 3' end (`runN`) — one genome per case, so every report row names
its case. `mismatch_zones.py` in this folder prints these tables.

### The 3' side at the defaults (`--3p-penalty 0`), across `-m`

| case | `-m 0` | `-m 1` | `-m 2` | `-m 3` | `-m 10` | reported |
|---|---|---|---|---|---|---|
| perfect | yes | yes | yes | yes | yes | `m=0 3'window=0` |
| off1 … off5 | no | no | no | no | no | — |
| off6, off7, off8 | no | yes | yes | yes | yes | `m=1 3'window=0` |
| run2, run3, run4 | no | no | no | no | no | — |

The 0.6.1 free zone is gone: a terminal mismatch, and the run of two that the free zone used to wave
through, are refused at **every** `-m`, 10 included. At the default a hit must be perfect over the last
5 bases, and `-m` only governs mismatches further in — which is now the documented behavior rather than
an accident of blast.

### The weighted window (`-m 4`, across `--3p-penalty`)

| case | `p=0` | `p=1` | `p=2` | `p=3` | `p=4` | reported |
|---|---|---|---|---|---|---|
| perfect | yes | yes | yes | yes | yes | `m=0 3'window=0` |
| off1 | no | no | no | no | yes | `m=1 3'window=1` |
| off2, off3, off4 | no | no | no | no | no | — |
| off5 | no | yes | yes | yes | yes | `m=1 3'window=1` |
| off6, off7, off8 | yes | yes | yes | yes | yes | `m=1 3'window=0` |

And the terminal base's G/C surcharge, isolated — one primer ending in A, one in G:

| case | `p=0` | `p=1` | `p=2` | `p=3` | `p=4` |
|---|---|---|---|---|---|
| terminal mismatch, primer ends A | no | no | no | **yes** | yes |
| terminal mismatch, primer ends G | no | no | no | no | **yes** |

That is the documented scale, reproduced exactly: a terminal mismatch costs 3, +1 when it breaks a G/C
pair (the first table's primer ends in G, hence its threshold of 4), and a mismatch five bases in costs
1. The README's own example — "a single A/T terminal mismatch requires `--3p-penalty 3`; a G/C terminal
mismatch requires 4" — is what the measurement found.

**One subtlety the documentation does not spell out:** a single mismatch 2–4 bases from the 3' end never
passes, even at `--3p-penalty 4`. blast will not carry an alignment across a mismatch that close to the
end — it trims the tail instead — and every trimmed base is scored in the window, so the run of unaligned
bases costs more than the one underlying mismatch would (off2 ≈ 3+2, off3 ≈ 3+2+2, …). Raising
`--3p-penalty` to 4 readmits exactly the terminal-mismatch case and nothing deeper. This is the same
blast behavior that made 0.6.1's dead zone; 0.8.0 inherits it as strictness rather than as a blind spot.

### The 5' side (`-m 0` unless stated)

| case | `-m 0` | `-m 1` | `--5p-truncation 2`, `-m 0` | reported |
|---|---|---|---|---|
| perfect | yes | yes | yes | `m=0 3'window=0` |
| off1 … off4 | yes | yes | off1–off2 yes, off3–off4 no | `m=0 3'window=0` |
| off5 … off8 | no | yes | no | `m=1 3'window=0` |

A mismatch 1–4 bases from the 5' end makes blast trim the tail, and up to 5 trimmed 5' bases are free —
so those bind at `-m 0`, as the default `--5p-truncation 5` says. From 5 bases in, blast keeps the
alignment and the mismatch counts against `-m` like any internal one. Tightening `--5p-truncation` to 2
refuses the 3- and 4-base trims. (0.6.1's allowance was a fixed 2, which section 2 shows was costing it
real hits.)

## 2. The same find results, designed under both versions

The four *Xylella* targets of `validation/xylella/`, 25 genomes. `find` ran once per target under
1.10.0; all three design arms consumed the same results folder, so the assays are identical and verdicts
join by name. Designs at the defaults (`-m 1`), 50 regions, up to 3 assays each:

| target | inclusion/exclusion | usable assays | qPCR selective, old → new | PCR selective, old → new |
|---|---|---|---|---|
| multiplex | 8 / 17 | 675 | 554 → 598 | 505 → 565 |
| fastidiosa | 6 / 19 | 529 | 390 → 441 | 351 → 422 |
| pauca | 4 / 21 | 758 | 672 → 676 | 655 → 660 |
| morus | 5 / 20 | 471 | 390 → 418 | 346 → 400 |

`compare_verdicts.py` in this folder joins the verdicts per assay. 2,433 assays were judged by both
versions:

**PCR mode** (old verdict in rows, new in columns):

| 1.9.3 + 0.6.1 | yes | no | total |
|---|---|---|---|
| yes | 1801 | 56 | 1857 |
| undecided (3' end) | 34 | 0 | 34 |
| no | 212 | 330 | 542 |
| **total** | **2047** | **386** | **2433** |

**qPCR mode:**

| 1.9.3 + 0.6.1 | yes | no | total |
|---|---|---|---|
| yes | 1974 | 32 | 2006 |
| undecided (3' end) | 12 | 0 | 12 |
| no | 147 | 268 | 415 |
| **total** | **2133** | **300** | **2433** |

88% of PCR verdicts and 92% of qPCR verdicts did not move. Both arms ran 1.9.3-or-later *readers*, so
none of this re-measures the phantom-pairing fix of 1.9.2 — the movement is the 0.8.0 hit model itself.

### Where the `undecided (3' end)` verdicts landed

All of them — 34 in PCR mode (25 fastidiosa, 9 multiplex) and 12 in qPCR mode (9 fastidiosa,
2 multiplex, 1 morus) — are now **yes**,
each with every inclusion genome amplified and zero exclusion amplifications. That verdict existed
because 0.6.1 could not count a difference under a primer's last two bases, so the check had amplicons it
could neither accept nor refuse. 0.8.0 refuses those bindings, so the exclusion amplifications are gone
and the check can finally answer. Note what kind of answer it is: the model now **takes a position** — a
primer whose 3' end sits on a mismatch does not extend, at the default — where 0.6.1 abstained. The
position is the orthodox one (it is why allele-specific primers put the discriminating base at the 3'
end), but it is a modeling choice, not bench evidence.

### What drives the rest of the movement

Classified from each arm's own report columns:

- **`no → yes` (212 PCR, 147 qPCR): 211 and 147 of them are old exclusion amplifications that 0.8.0
  refuses** — bindings that got through 0.6.1's uncounted terminal bases and now fail the 3'-end check.
  The one remaining PCR case is `pauca/Contig_928_73.7927_assay1`, which 0.6.1 scored 3/4 inclusion
  genomes: in ALM4 its site carries a 5'-proximal SNP, blast trims more of the tail than 0.6.1's fixed
  2-base allowance, and the hit was dropped. 0.8.0's `--5p-truncation 5` accepts it — a clean 95 bp
  product, counted mismatches 0 — which is precisely the false-negative fix its changelog claims.
- **`yes → no` (56 PCR, 32 qPCR): 53 and 29 are new exclusion amplifications**, the same 5'-truncation
  lever pointed the other way — exclusion-genome variants whose 5'-trimmed bindings 0.6.1 was refusing on
  an allowance that section 1 shows real sites exceed. The remaining 3 (same in both modes) are one morus
  region, `Contig_1123_56.952` assays 9–11, which falls from 5/5 to 3/5 inclusion genomes: in Riv19 and
  Riv25 every binding 0.6.1 reported carried `ReverseEndMismatch = -2` — the reverse primer's last two
  bases do not match those genomes' copies of the site. The free zone was hiding an inclusion-side
  defect: the region (a repeat, with copies 622 kb apart in two genomes) spells a reverse-primer 3' end
  that two of the five morus genomes do not hold. The old `yes` was wrong on the model's own terms.

### `--probe-mismatches 2` against the default

13 verdicts moved out of 2,433, qPCR mode, and all 13 moved `yes → no` (7 multiplex, 3 fastidiosa,
3 morus; pauca unchanged). The direction is forced: a looser probe budget can only *add* probe bindings,
so it can only surface exclusion amplifications the default probe refused to confirm — and on this panel
that is all it did. Nothing was rescued on the inclusion side. The knob is a stress test here, not a
sensitivity gain.

## What this does and does not show

- It shows 0.8.0's hit model does what its documentation says, knob by knob, including the G/C surcharge
  — and one undocumented consequence of blast's trimming, that `--3p-penalty` cannot readmit a lone
  mismatch 2–4 bases from the 3' end.
- It shows the old `undecided (3' end)` verdict dissolves entirely on this panel, and in the direction
  the ranking already leaned: those assays were ranked above refused ones, and all 46 are now clean
  passes.
- It shows the movement cuts both ways on defensible grounds: the model stopped refusing real 5'-variant
  bindings (one inclusion genome recovered, 53 exclusion amplifications surfaced) and stopped hiding
  3'-end defects (three false `yes` verdicts withdrawn on one repeat region).
- It does **not** show that any assay works. Nothing here was near a bench; the 46 resolved verdicts
  rest on the model's position that a 3'-terminal mismatch stops extension, which is a choice the bench
  can disagree with.

## Reproducing

```bash
# section 1
python mismatch_zones.py /path/to/insilicoPCR-linux-x64   # 0.8.0

# section 2: find once per target under 1.10.0, then three designs on the same results
bash validation/xylella/run.sh <work> 8 16                # or group_genomes.py + find per target
for t in multiplex fastidiosa pauca morus; do
  primer-finder design <work>/$t/results -o <work>/$t/assays_new    -t 8 --insilico-pcr <0.8.0-dir>
  <1.9.3>/primer_finder.py design <work>/$t/results -o <work>/$t/assays_old -t 8 --insilico-pcr <0.6.1-dir>
  primer-finder design <work>/$t/results -o <work>/$t/assays_probe2 -t 8 --insilico-pcr <0.8.0-dir> \
      --probe-mismatches 2
done
python compare_verdicts.py <work>
```

The twelve design runs took 66 minutes on 8 threads, each including its four in silico PCR passes over
the 25 genomes.
