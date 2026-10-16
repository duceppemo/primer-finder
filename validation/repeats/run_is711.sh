#!/usr/bin/env bash
# Does primer-finder find a published repeat-targeted assay? The genus-level IS711 real-time PCR assay of
# Bounaadja et al. (2009), against the Brucella panel.
#
#   bash validation/repeats/run_is711.sh /path/to/work /path/to/is711 32 64
#
# `work` is the Brucella harness's folder: it needs panel.tsv (panel.py), the panel blast database
# (census.py) and is711_census.tsv (is711_census.py). The second folder holds the groups and the runs.
#
# On the complete-genome group it runs the 2 x 2 that separates the two barriers in the way of a
# repeat-targeted marker: the duplication ceiling (-d) and the assembler. On the draft-heavy group it runs
# the setting that works, and designs on the region: the per-genome copy counts are what a draft panel
# costs, and they are what the design step reports.
set -euo pipefail

WORK=${1:?usage: run_is711.sh /work /is711 [threads] [memory_gb]}
OUT=${2:?usage: run_is711.sh /work /is711 [threads] [memory_gb]}
THREADS=${3:-16}
MEMORY=${4:-32}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ASSAY=$HERE/is711_assay.tsv
CHECK=$HERE/../brucella/check_assay.py

# -d has to reach the mean copy number of the group, which is a property of the group and not of the
# element: is711_panels.py prints both means, and these are what it prints.
D_COMPLETE=12
D_DRAFT=5

python "$HERE/is711_panels.py" "$WORK" -o "$OUT"

for assembler in skesa unitig; do
  for duplication in 1 $D_COMPLETE; do
    name=complete_d${duplication}_${assembler}
    echo "=== $name"
    rm -rf "$OUT/$name"
    primer-finder find -i "$OUT/complete/inclusion" -e "$OUT/complete/exclusion" -o "$OUT/$name" \
      -t "$THREADS" -m "$MEMORY" -k 99 --min-region 99 -d "$duplication" -a "$assembler" \
      > "$OUT/$name.log" 2>&1
    echo "    $(grep -c '>' "$OUT/$name/final_kmers.fasta") region(s) reported"
    python "$CHECK" "$OUT/$name" --target genus --assays "$ASSAY" | tail -4 || true
  done
done

name=draft_d${D_DRAFT}_unitig
echo "=== $name"
rm -rf "$OUT/$name"
primer-finder find -i "$OUT/draft/inclusion" -e "$OUT/draft/exclusion" -o "$OUT/$name" \
  -t "$THREADS" -m "$MEMORY" -k 99 --min-region 99 -d "$D_DRAFT" -a unitig > "$OUT/$name.log" 2>&1
echo "    $(grep -c '>' "$OUT/$name/final_kmers.fasta") region(s) reported"
python "$CHECK" "$OUT/$name" --target genus --assays "$ASSAY" --json "$OUT/$name.json" | tail -5 || true

for run in complete_d${D_COMPLETE}_unitig draft_d${D_DRAFT}_unitig; do
  echo "=== design on $run"
  rm -rf "$OUT/${run}_design" "$OUT/${run}_design_regions"
  python "$HERE/is711_design.py" "$OUT/$run" -o "$OUT/${run}_design" -t "$THREADS" \
      --assays "$ASSAY" || true
done
echo "=== done"
