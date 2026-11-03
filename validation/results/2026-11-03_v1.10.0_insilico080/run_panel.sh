#!/usr/bin/env bash
# Rerun the Xylella validation panel for primer-finder 1.10.0 + insilicoPCR 0.8.0, with a paired
# 1.9.3 + 0.6.1 design on the SAME find results so verdicts can be joined per assay, and a third
# design measuring --probe-mismatches 2.
set -uo pipefail

S=/tmp/claude-1000/-home-bioinfo-prog-primer-finder/1df804a1-6989-4610-bef1-ce42c6758ed4/scratchpad
W=$S/xylella-v1.10.0
ENVBIN=/home/bioinfo/miniconda3/envs/primer-finder_dev/bin
PY=$ENVBIN/python
export PATH=$ENVBIN:$PATH
NEW=/home/bioinfo/prog/primer-finder/primer_finder.py
OLD=$S/pf-1.9.3/primer_finder.py
XV=/home/bioinfo/prog/primer-finder/validation/xylella
I080=/home/bioinfo/prog/insilicoPCR-linux-x64
I061=$S/insilico-061/insilicoPCR-linux-x64

status=0
for target in multiplex fastidiosa pauca morus; do
    echo "=============== ${target} $(date +%H:%M:%S)"
    $PY "$XV/group_genomes.py" "$W/genomes/genomes" "$W/$target" --target "$target" || { status=1; continue; }
    if [ ! -f "$W/$target/results/final_kmers.fasta" ]; then
        $PY "$NEW" -i "$W/$target/inclusion" -e "$W/$target/exclusion" \
            -o "$W/$target/results" -t 8 -m 16 || { echo "find failed for $target"; status=1; continue; }
    fi
    echo "--- design new (1.10.0 + insilicoPCR 0.8.0) $(date +%H:%M:%S)"
    $PY "$NEW" design "$W/$target/results" -o "$W/$target/assays_new" -t 8 \
        --insilico-pcr "$I080" || { echo "design_new failed for $target"; status=1; }
    echo "--- design old (1.9.3 + insilicoPCR 0.6.1) $(date +%H:%M:%S)"
    $PY "$OLD" design "$W/$target/results" -o "$W/$target/assays_old" -t 8 \
        --insilico-pcr "$I061" || { echo "design_old failed for $target"; status=1; }
    echo "--- design probe2 (1.10.0, --probe-mismatches 2) $(date +%H:%M:%S)"
    $PY "$NEW" design "$W/$target/results" -o "$W/$target/assays_probe2" -t 8 \
        --insilico-pcr "$I080" --probe-mismatches 2 || { echo "design_probe2 failed for $target"; status=1; }
done
echo "=============== all targets done, status=$status $(date +%H:%M:%S)"
exit $status
