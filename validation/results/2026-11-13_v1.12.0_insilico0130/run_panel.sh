#!/usr/bin/env bash
# Design the Xylella panel with the same primer-finder under two insilicoPCR releases -- 0.11.0 and 0.13.0
# for the record -- on the SAME find results, so the two can be joined per assay. Three settings, each
# under both releases:
#
#   default   primer-finder's own (--mismatches 1, the probe budget following it)
#   probe2    --probe-mismatches 2, the stricter probe screen
#   m0        --mismatches 0, where the probe budget is 0 too
#
#   bash run_panel.sh <work_folder> <output_folder> <insilicoPCR-old> <insilicoPCR-new> [threads_per_design]
#
# <work_folder>        gets genomes/ (downloaded with ../../xylella/get_genomes.py when missing, which
#                      needs the NCBI datasets tool) and per target inclusion/, exclusion/ and results/
#                      (the find, skipped when already there)
# <output_folder>      gets <setting>_<old|new>/<target>/ and a log beside each
# <insilicoPCR-x.y.z>  an extracted portable release of that version, its jar, or a launcher script
#
# Runs the primer-finder this script sits in, and needs the commit this record was added in, which is
# primer-finder 1.12.0: a later primer-finder requires insilicoPCR 0.12.0 and refuses to design with
# 0.11.0. kmc, skesa, minimap2, blast and primer3_core must be on PATH; $PYTHON picks the interpreter. A
# design that is already there, verdicts and all, is not run again. The eight designs of a setting run
# side by side. Then:
#   python compare_panel.py <output_folder>
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage="usage: run_panel.sh <work_folder> <output_folder> <insilicoPCR-old> <insilicoPCR-new> [threads_per_design]"
W="${1:?${usage}}"
O="${2:?${usage}}"
OLD="${3:?${usage}}"
NEW="${4:?${usage}}"
threads="${5:-8}"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
PF="${here}/../../../primer_finder.py"
XV="${here}/../../xylella"
TARGETS="multiplex fastidiosa pauca morus"

if [ -z "${PY}" ]; then
    echo "No python interpreter found on PATH" >&2
    exit 1
fi
if [ ! -d "$W/genomes" ]; then
    "$PY" "$XV/get_genomes.py" "$W" || exit 1
fi

status=0
for target in $TARGETS; do
    "$PY" "$XV/group_genomes.py" "$W/genomes" "$W/$target" --target "$target" || exit 1
    if [ ! -f "$W/$target/results/final_kmers.fasta" ]; then
        "$PY" "$PF" -i "$W/$target/inclusion" -e "$W/$target/exclusion" -o "$W/$target/results" \
            -t "$threads" -m 16 > "$W/$target.find.log" 2>&1 || { echo "find failed for $target" >&2; exit 1; }
    fi
done

design() {  # design <setting> <options...>: the four targets under both releases, side by side
    local setting="$1"; shift
    local pids=()
    echo "=============== ${setting} $(date +%T)"
    for arm in old new; do
        release="$OLD"; [ "$arm" = new ] && release="$NEW"
        mkdir -p "$O/${setting}_${arm}"
        for target in $TARGETS; do
            head -n 1 "$O/${setting}_${arm}/$target/assays.tsv" 2> /dev/null | grep -q pcr_selective && continue
            ( "$PY" "$PF" design "$W/$target/results" -o "$O/${setting}_${arm}/$target" -t "$threads" \
                  --insilico-pcr "$release" "$@" > "$O/${setting}_${arm}/$target.log" 2>&1
              code=$?
              [ $code -eq 0 ] || echo "${setting}_${arm}/$target exit $code" >&2
              exit $code ) &
            pids+=($!)
        done
    done
    for pid in "${pids[@]}"; do
        wait "$pid" || status=1
    done
}

design default
design probe2 --probe-mismatches 2
design m0 --mismatches 0
echo "=============== done, status=$status $(date +%T)"
exit $status
