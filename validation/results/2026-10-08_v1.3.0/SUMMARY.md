# 2026-10-08 — primer-finder 1.3.0, how many mismatches in silico PCR should allow

`primer-finder design` runs the assays it proposes through
[insilicoPCR](https://github.com/duceppemo/insilicoPCR), which takes a `-m/--mismatches` tolerance. Until
this version primer-finder passed `0`. This record asks what that choice does, on the *Xylella* set of
[`../2026-10-08_v1.1.0/SUMMARY.md`](../2026-10-08_v1.1.0/SUMMARY.md), and settles the default at **1**.

Two things were run: a controlled test of what insilicoPCR's mismatch tolerance actually governs, and the
same 465 assays scored at `-m` 0, 1, 2 and 3.

**Still not evidence about PCR.** Nothing here was near a bench. What changed is that the model's own rule
is now measured rather than assumed, which is also why one sentence of the 1.2.0 documentation was wrong.

primer-finder 1.3.0 (commit after `bfc444f`), Primer3 2.6.1, blast 2.17.0+, insilicoPCR 0.6.1 with its
bundled blast 2.17.0+.

## 1. What `-m` governs: three zones along the primer

A 22-mer primer pair was placed in a synthetic 1.3 kb template, and the template — not the primer — was
mutated: one mismatch at a known distance from the left primer's 3' end, or a run of mismatches ending at
it. Each case is its own genome in the input folder, so every row of insilicoPCR's report names the case it
came from. Run at `-m` 0, 1, 2, 3, 5 and 10.

| Mismatch position | Amplified at `-m 0` | at `-m 1` | at `-m 10` | What insilicoPCR reports |
|---|---|---|---|---|
| none (control) | yes | yes | yes | `ForwardMismatches=0`, `ForwardEndMismatch=0` |
| last base (1 from the 3' end) | **yes** | yes | yes | `0`, `-1` |
| 2 from the 3' end | **yes** | yes | yes | `0`, `-2` |
| run of 2 ending at the 3' end | **yes** | yes | yes | `0`, `-2` |
| 3 from the 3' end | **no** | **no** | **no** | — |
| 4 from the 3' end | **no** | **no** | **no** | — |
| run of 3 or 4 ending at the 3' end | no | no | no | — |
| 5, 6, 7 or 8 from the 3' end | no | yes | yes | `1`, `0` |

So the primer has three zones, and `-m` only controls one of them:

1. **the last two bases are free.** A mismatch there is not counted as a mismatch at all. blast trims an
   unmatched base off the end of its alignment, and insilicoPCR reports how much was trimmed as a negative
   `EndMismatch` offset and calls the primer bound. This happens at every `-m`, including 0.
2. **3 and 4 bases from the 3' end is a dead zone.** Those never amplify, at any `-m` — up to 10 was
   tested. Keeping such a mismatch scores worse for blast than trimming three or four bases off, and a
   trim that long is rejected.
3. **5 or more bases in, `-m` decides.** Those mismatches are counted, one per primer: at `-m 1` the
   report's highest `ForwardMismatches` is 1, at `-m 2` it is 2. The tolerance is per primer, so an assay
   may carry twice as many as the number given.

**This is a property of the model, not of PCR.** It is what a blast alignment does at its ends, and it is
why the numbers below should not be read as measurements of how a mismatch behaves in a tube.

### What that corrects

The 1.2.0 documentation reported that assays whose only differences were a run of two at a primer's 3' end
were selective in none of the eight cases seen, and offered it as a weak consistency check on the ranking.
The mechanism is now clear and it is not about PCR: insilicoPCR does not count those two bases. The
observation holds at every `-m` — 0 of 60 assays with a single terminal difference and 0 of 5 with a run of
two are selective in PCR mode, whatever the tolerance — because the model cannot see a terminal difference
at all. It says nothing either way about whether such an assay would discriminate at the bench.

## 2. The same 465 assays at four tolerances

```bash
bash validation/xylella/run.sh <work> 16 32          # 25 genomes, multiplex: 8 inclusion, 17 exclusion
primer-finder design <work>/multiplex/results -o <work>/multiplex/assays -t 16 --max-regions 20
# then insilicoPCR over both groups and both primer files, at -m 0, 1, 2, 3
```

711 regions, 20 designed on, 472 assays proposed, **465 usable**: 178 specific by absence, 287 by a
difference under an oligo. 7 would amplify both groups and were carried no further.

Selective = amplifies all 8 inclusion genomes and none of the 17 exclusion genomes.

| `-m` | qPCR mode (with probe) | PCR mode (primers only) |
|---|---|---|
| 0 | 408 / 465 | 295 / 465 |
| 1 | 311 / 465 | 220 / 465 |
| 2 | 299 / 465 | 187 / 465 |
| 3 | 298 / 465 | 187 / 465 |

**All 465 amplified all 8 inclusion genomes at every tolerance.** The inclusion group is the positive
control here — these assays were designed on sequence shared by all 8 — and it never failed, so the
differences above are entirely about the exclusion group.

**`-m 3` is indistinguishable from `-m 2`.** Identical in PCR mode, one assay apart in qPCR mode. Of the
1287 exclusion hits at `-m 3`, only 21 are ones `-m 2` did not already have. There is nothing above 2 to
buy on this dataset.

### By what makes the assay selective

| `-m` | absence, qPCR | absence, PCR | difference, qPCR | difference, PCR |
|---|---|---|---|---|
| 0 | 169 / 178 | 160 / 178 | 239 / 287 | 135 / 287 |
| 1 | 162 / 178 | 148 / 178 | 149 / 287 | 72 / 287 |
| 2 | 157 / 178 | 128 / 178 | 142 / 287 | 59 / 287 |
| 3 | 156 / 178 | 128 / 178 | 142 / 287 | 59 / 287 |

The **absence** column is the interesting one, and the circularity that weakens the rest does not apply to
it: those assays carry no designed mismatch, their amplicon is simply missing from the exclusion genomes, so
anything insilicoPCR reports there is the primers finding **somewhere else in the genome to amplify**. At
`-m 0`, 18 of 178 did. At `-m 1`, 30. At `-m 2`, 50. Off-target amplification is a real specificity risk,
it has nothing to do with the designed difference, and `-m 0` hides two thirds of what `-m 2` finds.

### Where the exclusion hits come from

Every hit in the exclusion group, in PCR mode, classified by what insilicoPCR itself reported:

| `-m` | perfect match | terminal mismatch only (the free zone) | with a counted mismatch | total |
|---|---|---|---|---|
| 0 | 163 | 274 | 0 | 437 |
| 1 | 163 | 274 | 309 | 746 |
| 2 | 163 | 274 | 829 | 1266 |
| 3 | 163 | 274 | 850 | 1287 |

The first two columns do not move, by construction: 163 hits are exact matches — off-target amplification —
and 274 are the allele-specific primers whose 3' end sits on a difference the model does not count.

### By how many differences sit under the primers

The 287 difference-based assays, selective / total:

| `-m` | 0 (probe only) | 1 | 2 | 3 | 4 or more |
|---|---|---|---|---|---|
| **qPCR mode** | | | | | |
| 0 | 58/58 (100%) | 65/105 (62%) | 60/67 (90%) | 24/25 (96%) | 32/32 (100%) |
| 1 | 25/58 (43%) | 38/105 (36%) | 40/67 (60%) | 16/25 (64%) | 30/32 (94%) |
| 2 | 25/58 (43%) | 37/105 (35%) | 37/67 (55%) | 14/25 (56%) | 29/32 (91%) |
| **PCR mode** | | | | | |
| 0 | 0/58 (0%) | 32/105 (30%) | 53/67 (79%) | 19/25 (76%) | 31/32 (97%) |
| 1 | 0/58 (0%) | 7/105 (7%) | 26/67 (39%) | 10/25 (40%) | 29/32 (91%) |
| 2 | 0/58 (0%) | 4/105 (4%) | 20/67 (30%) | 7/25 (28%) | 28/32 (88%) |

More differences still means a larger selective share at every tolerance, and the gradient steepens: at
`-m 1`, an assay resting on one difference passes in PCR mode 7% of the time against 91% for one resting on
four or more. Read that as the model's rule, not as biology — it is close to a restatement of the
tolerance. The useful part is that the ordering survives the tolerance changing, which is not guaranteed.

### Does the ranking survive?

Selective among the first N usable assays, in the order `assays.tsv` writes them:

| `-m` | top 10 | top 25 | top 50 | top 100 |
|---|---|---|---|---|
| 0, PCR | 10/10 | 25/25 | 50/50 | 100/100 |
| 1, PCR | 10/10 | 25/25 | 50/50 | 97/100 |
| 2, PCR | 10/10 | 24/25 | 47/50 | 90/100 |
| 0–2, qPCR | 10/10 | 25/25 | 50/50 | 99–100/100 |

The top of the ranking holds at every tolerance, so the order is not an artefact of `-m 0`.

## 3. Why the default is 1

- `-m 0` assumes **any** mismatch five or more bases from the 3' end stops a primer. That is the optimistic
  end of what is known: a single internal mismatch frequently does not stop amplification, which is the
  stated reason double-mismatch allele-specific designs exist
  ([Lefever et al. 2019](https://doi.org/10.1038/s41598-019-38581-z)). It also hid 12 of the 30 absence-based
  assays that amplify something off-target.
- `-m 2` assumes **two** mismatches per primer still bind. That is beyond what the same literature
  supports: two mismatches are what a deliberate allele-specific design uses to discriminate. It is a
  sensible stress test, not a default.
- `-m 3` measures nothing `-m 2` did not.

So **1**, with 2 available for a stricter pass:

```bash
primer-finder design results/ -o assays/ --insilico-pcr <dir>                   # -m 1, the default
primer-finder design results/ -o assays/ --insilico-pcr <dir> --mismatches 2    # the stricter check
```

An assay that passes at `-m 2` rests on something the model cannot explain away. 187 of these 465 did, in
PCR mode, with primers alone.

## The command at its new default

Run end to end, with the design step driving insilicoPCR itself rather than a script doing it afterwards:

```bash
primer-finder design <work>/multiplex/results -o <work>/multiplex/assays -t 32 --max-regions 20 \
    --insilico-pcr /path/to/insilicoPCR-linux-x64
```

2 min 25 s for 20 regions on 32 threads, including all four in silico PCR runs over the 25 genomes. It
reported 311 of 465 selective in qPCR mode and 220 of 465 in PCR mode — the same two numbers the sweep
above computed independently from the reports, which is the check that the default is wired through.

### What the `no` verdicts are made of

The same run, split by what insilicoPCR itself reported for each exclusion amplification. A genome counts as
terminal when every amplicon reported for it needed a base trimmed off a primer's 3' end:

| | qPCR mode | PCR mode |
|---|---|---|
| selective | 311 | 220 |
| `no`, every exclusion amplification terminal | 105 | 155 |
| `no`, some of them terminal | 2 | 7 |
| `no`, on amplifications the model can defend | 47 | 83 |

So in PCR mode the check refuses 90 of 465 assays, and only 83 of those on amplifications with no terminal
component at all; the 155 fully terminal ones it cannot judge, because every amplification it reported for
them rests on bases it does not count. Of those 155, **121** carry a run of differences at a 3' end — 96 of
one base, 22 of two, 3 of three, which is the shape the ranking prefers — and 34 carry none. Cut by what
makes them selective instead, 127 are difference-based and 28 specific by absence, where the off-target
amplicon elsewhere in the exclusion genome itself only binds through an ignored terminal mismatch. The two
cuts cross: of the 127 difference-based, 102 have a terminal run; of the 28 absence-based, 19 do.

On the inclusion side the question does not arise: of 3,720 assay-and-genome amplifications, none depended
on a trimmed primer end. Every inclusion genome that counted had at least one clean amplicon.

`assays.tsv` reports this per assay as `qpcr_exclusion_terminal_only` and `pcr_exclusion_terminal_only`, and
those assays get a verdict of their own — `undecided (3' end)` rather than `no`, since there was nothing the
check could refuse. They are ranked above the assays it refused on evidence it can defend. The counts, from
the run's own `assays.tsv` (`design_info.json` records the first two rows; the third is the remainder):

| `*_selective` | qPCR mode | PCR mode |
|---|---|---|
| `yes` | 311 | 220 |
| `undecided (3' end)` | 105 | 155 |
| `no` | 49 | 90 |

## What this does and does not show

- It shows what insilicoPCR's `-m` governs, exactly, including that it does not govern the last two bases
  of a primer and cannot be made to.
- It shows the ranking's order is stable across tolerances, and that the absence-based assays — the ones
  the ranking puts first — are the ones that hold up.
- It shows that most of what this check does not clear, it does not clear for a reason it cannot defend:
  of 245 non-selective assays in PCR mode, 155 rest entirely on bases it does not count and only 90 are
  refusals it can stand behind.
- It does **not** show that any of these assays works. No assay here has been near a bench, and the
  gradient with the number of differences remains close to a restatement of the model's own rule.

## Reproducing

```bash
bash validation/xylella/run.sh /path/to/work 16 32
primer-finder design /path/to/work/multiplex/results -o /path/to/work/multiplex/assays \
    -t 16 --max-regions 20
for m in 0 1 2 3; do
  for kind in qpcr pcr; do
    for group in inclusion exclusion; do
      <insilicoPCR> -i /path/to/work/multiplex/$group \
        -p /path/to/work/multiplex/assays/assays_${kind}.fasta \
        -o /path/to/work/sweep/m$m/${kind}_${group} -t 32 -m $m
    done
  done
done
```

The three-zone test is `mismatch_zones.py` in this folder: it writes the synthetic genomes, runs
insilicoPCR at each tolerance, and prints the table of section 1.
