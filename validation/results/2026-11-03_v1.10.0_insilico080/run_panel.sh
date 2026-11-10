#!/usr/bin/env bash
# Rerun the Xylella validation panel for primer-finder 1.10.0 + insilicoPCR 0.8.0, with a paired
# 1.9.3 + 0.6.1 design on the SAME find results so verdicts can be joined per assay, and a third
# design measuring --probe-mismatches 2.
#
#   bash run_panel.sh <work_folder> <insilicoPCR-0.8.0> <insilicoPCR-0.6.1> <primer-finder-1.9.3> [threads] [memory_GB]
#
# <work_folder>         gets genomes/ (downloaded with ../../xylella/get_genomes.py when missing), and per
#                       target inclusion/, exclusion/, results/ (the find, skipped when already there) and
#                       assays_new/, assays_old/, assays_probe2/
# <insilicoPCR-x.y.z>   an extracted portable release of that version, its jar, or a launcher script
# <primer-finder-1.9.3> a checkout of the v1.9.3 tag, e.g. `git worktree add <dir> v1.9.3`
#
# The "new" arm runs the primer-finder this script sits in; check out v1.10.0 to reproduce the record
# exactly. kmc, skesa, minimap2, blast and primer3_core must be on PATH; $PYTHON picks the interpreter.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage="usage: run_panel.sh <work_folder> <insilicoPCR-0.8.0> <insilicoPCR-0.6.1> <primer-finder-1.9.3> [threads] [memory_GB]"
W="${1:?${usage}}"
I080="${2:?${usage}}"
I061="${3:?${usage}}"
OLD="${4:?${usage}}/primer_finder.py"
threads="${5:-8}"
memory="${6:-16}"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
NEW="${here}/../../../primer_finder.py"
XV="${here}/../../xylella"

if [ -z "${PY}" ]; then
    echo "No python interpreter found on PATH" >&2
    exit 1
fi
if [ ! -f "${OLD}" ]; then
    echo "No primer_finder.py in the 1.9.3 checkout: ${OLD}" >&2
    exit 1
fi
if [ ! -d "$W/genomes" ]; then
    "$PY" "$XV/get_genomes.py" "$W" || exit 1
fi

status=0
for target in multiplex fastidiosa pauca morus; do
    echo "=============== ${target} $(date +%H:%M:%S)"
    "$PY" "$XV/group_genomes.py" "$W/genomes" "$W/$target" --target "$target" || { status=1; continue; }
    if [ ! -f "$W/$target/results/final_kmers.fasta" ]; then
        "$PY" "$NEW" -i "$W/$target/inclusion" -e "$W/$target/exclusion" \
            -o "$W/$target/results" -t "$threads" -m "$memory" || { echo "find failed for $target"; status=1; continue; }
    fi
    echo "--- design new (1.10.0 + insilicoPCR 0.8.0) $(date +%H:%M:%S)"
    "$PY" "$NEW" design "$W/$target/results" -o "$W/$target/assays_new" -t "$threads" \
        --insilico-pcr "$I080" || { echo "design_new failed for $target"; status=1; }
    echo "--- design old (1.9.3 + insilicoPCR 0.6.1) $(date +%H:%M:%S)"
    "$PY" "$OLD" design "$W/$target/results" -o "$W/$target/assays_old" -t "$threads" \
        --insilico-pcr "$I061" || { echo "design_old failed for $target"; status=1; }
    echo "--- design probe2 (1.10.0, --probe-mismatches 2) $(date +%H:%M:%S)"
    "$PY" "$NEW" design "$W/$target/results" -o "$W/$target/assays_probe2" -t "$threads" \
        --insilico-pcr "$I080" --probe-mismatches 2 || { echo "design_probe2 failed for $target"; status=1; }
done
echo "=============== all targets done, status=$status $(date +%H:%M:%S)"
exit $status
