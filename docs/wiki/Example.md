# Example

The repository holds a small simulated dataset whose answer is known, so that an installation can be checked
in a few seconds.

```bash
bash example/run_example.sh [output_folder] [threads] [memory_GB]
```

It writes the dataset, runs primer-finder on it, and checks the results. The default output folder is
`example/example_output`, with 4 threads and 8 GB.

## The dataset

`example/make_example.py` writes eight genomes: four inclusion and four exclusion, one of each group
gzipped, so that both file types are exercised.

Every genome is the same 20 kb random backbone (the generator is seeded, so the dataset is identical every
time) with its own 20 scattered single-base changes, none of them within 500 bases of the planted region.
What the four inclusion genomes all share, and the exclusion genomes do not, is one planted region around
position 9,000:

| Planted variant | Position |
|---|---|
| mismatch | 9,000 |
| mismatch | 9,008 |
| mismatch | 9,015 |
| 15-base insertion `GGCATTACGTTAACC` | 9,050 |
| mismatch | 9,100 |

Three mismatches within 16 bases, an insertion, and one isolated mismatch that could not carry an assay on
its own. `example/data/truth.json` records all of it, with the private changes of each genome.

## What primer-finder finds

With the default 99-mers, the region yields 214 inclusion-specific kmers, which assemble into one contig of
260 bases. It maps to the exclusion genome with
`72=1X7=1X6=1X34=15I50=1X72=` — in this dataset the contig happens to be assembled in the opposite
orientation — is present in all four inclusion genomes, and keeps 19 variant positions shared by all four
exclusion genomes:

```
>Contig_1_67.5167 [72, 123, 124, ..., 137, 172, 179, 187]
CACC...TGTcCCATG...TAGACggttaacgtaatgccTGTG...TTCCCaGCGCATcTCGTTAGaGGAT...
```

The three clustered mismatches (172, 179, 187 — seven and eight bases apart) and the insertion (123-137) are
what an assay would be designed on.

## What the checks verify

`example/check_example.py` is run at the end and prints one line per check:

- exactly one contig in `final_kmers.fasta`;
- it carries the planted insertion;
- the planted variants are in lower case, and at least three of them are within 21 bases;
- the header lists the variant positions;
- the contig is a substring of all four inclusion genomes, in either orientation, and of none of the
  exclusion genomes;
- `run_info.json` counts one final contig, and the first exclusion genome was used as the reference.

It ends with `All checks passed`. The CI runs the same script on every push, in a conda environment built
from `environment.yml`, so the example also checks that the pipeline still works with current versions of
KMC, SKESA, minimap2 and blast.

## Running it by hand

```bash
python example/make_example.py /tmp/example/data
primer-finder -i /tmp/example/data/inclusion -e /tmp/example/data/exclusion -o /tmp/example/results \
    -t 4 -m 8 --keep-intermediate
python example/check_example.py /tmp/example/results /tmp/example/data
```

`--keep-intermediate` is a good way to see what each step produced; see [Outputs](Outputs).
