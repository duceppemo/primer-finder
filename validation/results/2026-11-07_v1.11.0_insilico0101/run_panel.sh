#!/usr/bin/env bash
# Design the Xylella panel with primer-finder 1.11.0 + insilicoPCR 0.10.1, on the SAME find results the
# 0.8.0 arms of ../2026-11-03_v1.10.0_insilico080/ and ../2026-11-03_3prime_dependence/ used, so the two
# insilicoPCR releases can be joined per assay. The four targets run side by side.
#
#   bash run_panel.sh <work_folder> <output_folder> <insilicoPCR-0.10.1> [threads_per_target]
#
# <work_folder>         the one ../2026-11-03_v1.10.0_insilico080/run_panel.sh builds, with the 3'-end
#                       record's assays_held/ arm designed in it (see that record's SUMMARY)
# <output_folder>       gets one design folder and one log per target
# <insilicoPCR-0.10.1>  an extracted portable release, its jar, or a launcher script
#
# Runs the primer-finder this script sits in; check out v1.11.0 to reproduce the record exactly. blast
# and primer3_core must be on PATH; $PYTHON picks the interpreter. Then:
#   python compare_versions.py <work_folder> <output_folder>
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage="usage: run_panel.sh <work_folder> <output_folder> <insilicoPCR-0.10.1> [threads_per_target]"
W="${1:?${usage}}"
O="${2:?${usage}}"
I0101="${3:?${usage}}"
threads="${4:-8}"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
PF="${here}/../../../primer_finder.py"

if [ -z "${PY}" ]; then
    echo "No python interpreter found on PATH" >&2
    exit 1
fi

mkdir -p "$O"
status=0
pids=()
for t in multiplex fastidiosa pauca morus; do
    ( "$PY" "$PF" design "$W/$t/results" -o "$O/$t" -t "$threads" \
          --insilico-pcr "$I0101" > "$O/$t.log" 2>&1
      code=$?
      echo "$t exit $code $(date +%T)"
      exit $code ) &
    pids+=($!)
done
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
exit $status
