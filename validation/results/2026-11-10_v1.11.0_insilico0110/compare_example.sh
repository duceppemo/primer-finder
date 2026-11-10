#!/usr/bin/env bash
# Design the example dataset with the same primer-finder under two insilicoPCR releases, and say whether
# anything the in silico check writes or reads differs. Then, under the newer release only, design again
# into the same folder and run the run_insilico_pcr.sh the design wrote: 0.11.0 writes a report.html
# beside report.tsv, and neither a rerun nor the script may be disturbed by it.
#
#   bash compare_example.sh <output_folder> <insilicoPCR-0.10.1> <insilicoPCR-0.11.0> [threads]
#
# <insilicoPCR-...>  an extracted portable release, its jar, or a launcher script
#
# Runs the primer-finder this script sits in; check out v1.11.0 to reproduce the record exactly. kmc,
# skesa, minimap2, blast and primer3_core must be on PATH; $PYTHON picks the interpreter. Exits 1 if
# anything differs.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
usage="usage: compare_example.sh <output_folder> <insilicoPCR-0.10.1> <insilicoPCR-0.11.0> [threads]"
O="${1:?${usage}}"
OLD="${2:?${usage}}"
NEW="${3:?${usage}}"
threads="${4:-8}"
PY="${PYTHON:-$(command -v python || command -v python3 || true)}"
root="${here}/../../.."
PF="${root}/primer_finder.py"
report="consolidated_report/report.tsv"

if [ -z "${PY}" ]; then
    echo "No python interpreter found on PATH" >&2
    exit 1
fi

mkdir -p "$O"
PYTHON="$PY" bash "${root}/example/run_example.sh" "$O/example" "$threads" 8 > "$O/example.log" 2>&1 \
    || { echo "the example search failed, see $O/example.log" >&2; exit 1; }

status=0
same() {  # same <label> <file> <file>
    if cmp -s "$2" "$3"; then
        echo "   byte-identical  $1"
    else
        echo "   DIFFERS         $1"
        status=1
    fi
}

for arm in old new; do
    release="$OLD"; [ "$arm" = new ] && release="$NEW"
    "$PY" "$PF" design "$O/example/results" -o "$O/$arm" -t "$threads" --insilico-pcr "$release" \
        > "$O/$arm.log" 2>&1 || { echo "design ($arm) failed, see $O/$arm.log" >&2; exit 1; }
    echo "$arm: $(grep -o 'insilicoPCR: .*' "$O/$arm.log")"
done

echo "== the two releases"
for file in assays.tsv assays_pcr.fasta assays_qpcr.fasta; do
    same "$file" "$O/old/$file" "$O/new/$file"
done
for dir in "$O"/old/insilico_pcr/*/; do
    pass="$(basename "$dir")"
    rows=$(($(wc -l < "$dir/$report") - 1))
    same "$pass report.tsv ($rows rows)" "$dir/$report" "$O/new/insilico_pcr/$pass/$report"
done

echo "== the newer release, designed again into the same folder"
cp -r "$O/new" "$O/new_first"
"$PY" "$PF" design "$O/example/results" -o "$O/new" -t "$threads" --insilico-pcr "$NEW" \
    > "$O/new_again.log" 2>&1 || { echo "the second design failed, see $O/new_again.log" >&2; exit 1; }
same "assays.tsv" "$O/new_first/assays.tsv" "$O/new/assays.tsv"
for dir in "$O"/new_first/insilico_pcr/*/; do
    pass="$(basename "$dir")"
    same "$pass report.tsv" "$dir/$report" "$O/new/insilico_pcr/$pass/$report"
done
strays=$(find "$O/new/insilico_pcr" -name 'report(*' | wc -l)
echo "   numbered report(n) files left behind: $strays"
[ "$strays" -eq 0 ] || status=1

echo "== the newer release, through the run_insilico_pcr.sh the design wrote"
bash "$O/new/run_insilico_pcr.sh" "$O/script" > "$O/script.log" 2>&1 \
    || { echo "run_insilico_pcr.sh failed, see $O/script.log" >&2; exit 1; }
for dir in "$O"/new/insilico_pcr/*/; do
    pass="$(basename "$dir")"
    same "assays_$pass report.tsv" "$dir/$report" "$O/script/assays_$pass/$report"
done

exit $status
