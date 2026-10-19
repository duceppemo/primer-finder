# 2026-10-19 — what blast's low-complexity masking was hiding, and a region held in pieces

Two findings of the twelve-part review of 2026-10-18/19, both of the kind where every component behaves
correctly and the composition is wrong, measured here before they were fixed and after.

**Nothing here was near a bench.** The runs are of the tool against sequence, and the synthetic panels are
built to have one property each. Where a number is about blast or about primer-finder, it is a measurement
of the program; it says nothing about how an assay behaves in a tube.

## 1. `blastn` cannot find what it masks

`blastn` masks low-complexity stretches of its **query** by default (SDUST), and a masked stretch is not
found however exactly the subject holds it. Nine 200-base shapes, each planted in the middle of the same
6 kb random genome and so held by it exactly, searched for in that genome — `low_complexity.py`:

| shape | hits, masking on | hits, `-dust no` | dustmasker | distinct 11-mers | triplet score |
|---|---|---|---|---|---|
| (AAT) x 67 | **0** | 7 | all 200 | 3 | 32.7 |
| random AT | **0** | 3 | 196 of 200 | 186 | 12.1 |
| random AC | **0** | 1 | all 200 | 179 | 12.4 |
| A x 200 | **0** | 7 | all 200 | 1 | 99.0 |
| random ACGT | 1 | 1 | nothing | 190 | 1.4 |
| 80% AT, random | 1 | 1 | nothing | 190 | 3.6 |
| (ACGTT) x 40 | **0** | 7 | all 200 | 5 | 19.4 |
| (ACGTTGGCAT) x 20 | **0** | 7 | all 200 | 10 | 9.5 |
| half ACGT, half (AAT) | 1 | 1 | the second half | 103 | 9.5 |

Six of the nine are invisible to a default search of a genome that holds them exactly. `dustmasker` —
blast's own masker, from the same package — agrees with `blastn` on all nine, including which half of the
mixed one is masked, which is why primer-finder asks it rather than guessing.

**It rules out guessing.** A random two-letter sequence holds 186 of its 190 possible 11-base seeds and is
masked all the same, so seed diversity does not predict this; the whole-sequence triplet score separates
these nine but with no margin to speak of (3.6 unmasked against 9.5 masked).

## 2. What that did to a run, end to end

Three inclusion and three exclusion genomes, 6 kb of random sequence either side of a 200-base `AAT`
element that **both** groups carry, with one base of the element differing between them. That difference
makes the 51-mers over it inclusion-specific, so `find` assembles a 101-base region of pure low-complexity
sequence that the exclusion genomes hold 100 of the 101 bases of. `one_difference.py`:

| the run | what it reported |
|---|---|
| no second look at all (every version up to 1.9.2) | **stopped**: "No candidate contig is present in all the inclusion genomes" |
| a second look in the presence search only | 1 region of 101 bases, with all 101 recorded as absent from every exclusion genome |
| as it now stands | nothing to report |

Each line is a different failure or the absence of one:

- **the region was lost from its own inclusion group.** All three inclusion genomes carry it exactly, and
  the presence search could not see it, so the run ended with nothing to design on;
- **and then claimed absent from the exclusion group.** With the presence search fixed and the exclusion
  search still masking, the whole 101-base region is recorded as missing from every exclusion genome —
  the strongest claim this tool makes, which ranks above every difference-based assay and is exempt from
  `--min-primer-differences`. minimap2 leaves such a contig unmapped too, so nothing else in the run
  contradicts it;
- **as it stands there is nothing to report, for the right reason.** The unmasked search finds the region
  in each exclusion genome 38 times over — every frame of a three-base period aligns — and while a given
  alignment carries the planted mismatch, others cover that base and match it. `shared_variants` counts a
  position as a difference only where **every** copy in the genome differs there, since one matching copy
  is enough for an assay to amplify that genome. So the single difference in a microsatellite is not one,
  which is a fact about periodic sequence rather than about this tool.

## 3. On a real panel: no change at all

The *Xylella fastidiosa* set of [`../2026-10-08_v1.1.0/SUMMARY.md`](../2026-10-08_v1.1.0/SUMMARY.md), 25
complete RefSeq genomes, all four published subspecies targets, run twice — as the code stands and with
the unmasked second look disabled (`absence_before_after.py`):

| target | regions | bases | regions with an absent stretch | absent bases |
|---|---|---|---|---|
| multiplex | 949 | 482,042 | 14 | 7,202 |
| fastidiosa | 642 | 232,083 | 6 | 2,201 |
| pauca | 2,047 | 816,514 | 61 | 11,018 |
| morus | 434 | 167,485 | 3 | 2,029 |

**Identical in both runs, on every target**: the same regions, the same absent stretches, the same bases.
No region lost its absence and none was added. Complete genomes of a GC-balanced organism, and the
inclusion-specific contigs of these panels hold no low-complexity stretch for the masker to take out.

So this fix changes nothing in the records of this panel, and none of them needs an erratum for it. What
it closes is a defect of the shape in section 2, which did not arise here — a run on draft assemblies, or
on a target in or beside a microsatellite, is where it would.

**It costs nothing measurable.** One `dustmasker` call per search step, and a second `blastn` only for the
stretches the masker took out — none, on this panel. Wall clock for the whole `find` run, as it stands
against with the second look disabled: multiplex 36 s / 33 s, fastidiosa 33 s / 32 s, pauca 41 s / 40 s,
morus 23 s / 25 s, on 16 threads. The search is asked about the stretches a genome did **not** align only
where the masker masked them, which is what keeps it from being a second search of everything a genome
lacks — most of the query, on most genomes, since that is the point of the search.

## 4. A region held in pieces

The second finding, measured on the shape that produces it: three genomes carrying a 25-base unit six
times over between 6 kb flanks, and an exclusion group with no copy of it.

At `-d 1` the kmers inside the repeat occur six times each and are dropped as repeated, which is what `-d`
is for and is documented. The walk then splices the kmers either side of the repeat into
flank + repeat + flank, and the region that comes out carries **two** copies of the unit where the genome
has six: 150 bases that are in no genome at all, on either strand. Every step is right about its own
question — the kmers are inclusion-specific, the splice is what the kmers determine, and both 100-base
halves of the result are in every genome — and until now the run said the region was present in all three
inclusion genomes and nothing more. `design` then priced products from it: on the run the review found,
140, 150 and 142 bases where insilicoPCR measured 284, 294 and 286.

The presence fraction cannot see this, by design: it merges intervals of the **region** and asks nothing
about where in the genome they came from, which is exactly what makes it the right number for a draft
assembly that broke a region between two of its own contigs. So two things are now measured instead of one:

- `find` asks each genome how much of the region it holds as **one stretch** — the best chain of hits
  running the same way through both — and names a region no genome holds in one piece, in a warning and in
  `run_info.json` as `held_in_pieces`. On the panel above: 0.667 for that 150-base region, in all three
  genomes;
- `design` measures every amplicon it counts **in the genome**, from the first base of one primer site to
  the last of the other, and reports it as `inclusion_product_span` beside Primer3's `product_size`. Where
  they disagree the run says so once, naming the worst assay.

**Nothing is dropped for either.** The region is specific, an assay inside either half may be perfectly
good, and the product size is what cannot be trusted. That was the decision asked for and taken on
2026-10-19: measure and report, change no verdict.

## What this does and does not show

- It shows what blast's masking hides, exactly, and that blast's own masker is the only thing that
  predicts it among those tried.
- It shows one run where masking turned a region every inclusion genome holds into first a lost region and
  then an unearned absence, and that neither happens now.
- It shows the fix changes nothing on the *Xylella* panel, at no measurable cost.
- It does **not** show what either fix does to the *Brucella* or *Listeria* panels, which were not re-run
  here: both are draft-heavy, which is where a low-complexity contig is likeliest.
- It does **not** show that any assay works. No assay here has been near a bench.

## Reproducing

```bash
python validation/results/2026-10-19_v1.9.3_masking/low_complexity.py <work>
python validation/results/2026-10-19_v1.9.3_masking/one_difference.py <work> 8 8
bash validation/xylella/run.sh <xylella_work> 16 24      # the panel, if it is not already there
python validation/results/2026-10-19_v1.9.3_masking/absence_before_after.py <xylella_work> 16 24
```

The repeat panel of section 4 is built by
`tests/test_pipeline.py::test_a_region_spliced_across_a_tandem_repeat_with_the_real_programs`, which
asserts the shape and the figure.
