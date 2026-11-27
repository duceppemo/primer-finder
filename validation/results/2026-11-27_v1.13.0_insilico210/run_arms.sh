#!/usr/bin/env bash
# Design the example and the Xylella panel twice, on the SAME find results, at the default settings:
#
#   old  primer-finder 1.13.0 with insilicoPCR 0.14.0 (the last pair checked before the Tm rule; 0.14.0
#        gives byte-identical reports to 0.12.0 and 0.13.0, ../2026-11-13_v1.12.0_insilico0130/)
#   new  the primer-finder this script sits in, with insilicoPCR 2.1.0
#
#   bash run_arms.sh <primer-finder-1.13.0> <insilicoPCR-0.14.0> <insilicoPCR-2.1.0> <panel_work> <out> [threads]
#
# <primer-finder-1.13.0>  a checkout of the v1.13.0 tag (git worktree add <dir> v1.13.0)
# <insilicoPCR-x.y.z>     an extracted portable release of that version
# <panel_work>            the Xylella work folder of ../2026-11-13_v1.12.0_insilico0130/run_panel.sh: genomes/
#                         and per target inclusion/, exclusion/ and results/ (built there when missing)
# <out>                   gets example/ (the example's find), then <old|new>/<target>/ and a log beside each
#
# kmc, skesa, minimap2, blast and primer3_core must be on PATH; $PYTHON picks the interpreter. A design
# that is already there, verdicts and all, is not run again. The ten designs run side by side. The old
# arm's panel designs take about 18 GB in all, nearly all of it insilicoPCR 0.14.0's detailed_report/
# folders; the new arm's about 0.4 GB. Then:
#   python compare_arms.py <out>
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage="usage: run_arms.sh <primer-finder-1.13.0> <insilicoPCR-0.14.0> <insilicoPCR-2.1.0> <panel_work> <out> [threads]"
OLD_PF="${1:?${usage}}/primer_finder.py"
OLD="${2:?${usage}}"
NEW="${3:?${usage}}"
W="${4:?${usage}}"
O="${5:?${usage}}"
threads="${6:-8}"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
root="${here}/../../.."
NEW_PF="${root}/primer_finder.py"
XV="${root}/validation/xylella"
TARGETS="multiplex fastidiosa pauca morus"

if [ -z "${PY}" ]; then
    echo "No python interpreter found on PATH" >&2
    exit 1
fi
mkdir -p "$O"
if [ ! -f "$O/example/results/final_kmers.fasta" ]; then
    PYTHON="$PY" bash "${root}/example/run_example.sh" "$O/example" "$threads" 8 > "$O/example.log" 2>&1 \
        || { echo "the example search failed, see $O/example.log" >&2; exit 1; }
fi
if [ ! -d "$W/genomes" ]; then
    "$PY" "$XV/get_genomes.py" "$W" || exit 1
fi
for target in $TARGETS; do
    "$PY" "$XV/group_genomes.py" "$W/genomes" "$W/$target" --target "$target" || exit 1
    if [ ! -f "$W/$target/results/final_kmers.fasta" ]; then
        "$PY" "$NEW_PF" -i "$W/$target/inclusion" -e "$W/$target/exclusion" -o "$W/$target/results" \
            -t "$threads" -m 16 > "$W/$target.find.log" 2>&1 || { echo "find failed for $target" >&2; exit 1; }
    fi
done

status=0
pids=()
echo "=============== designs $(date +%T)"
for arm in old new; do
    pf="$OLD_PF"; release="$OLD"
    [ "$arm" = new ] && { pf="$NEW_PF"; release="$NEW"; }
    mkdir -p "$O/$arm"
    for target in example $TARGETS; do
        results="$W/$target/results"
        [ "$target" = example ] && results="$O/example/results"
        head -n 1 "$O/$arm/$target/assays.tsv" 2> /dev/null | grep -q pcr_selective && continue
        ( "$PY" "$pf" design "$results" -o "$O/$arm/$target" -t "$threads" --insilico-pcr "$release" \
              > "$O/$arm/$target.log" 2>&1
          code=$?
          [ $code -eq 0 ] || echo "$arm/$target exit $code" >&2
          exit $code ) &
        pids+=($!)
    done
done
for pid in "${pids[@]}"; do
    wait "$pid" || status=1
done
echo "=============== done, status=$status $(date +%T)"
exit $status
