#!/usr/bin/env bash
# Validate primer-finder against the published Xylella fastidiosa subspecies qPCR assays of
# Dupas et al. 2019 (doi:10.3389/fpls.2019.01732).
#
#   bash validation/xylella/run.sh <work_folder> [threads] [memory_GB]
#
# It downloads 25 complete RefSeq genomes (NCBI datasets), then for each of the four subspecies assays
# runs primer-finder with that subspecies as the inclusion group and the others as the exclusion group,
# and checks that the published primers and probe fall inside one reported region.
#
# Needs: primer-finder and its programs, and the NCBI datasets command line tool.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work="${1:?usage: run.sh <work_folder> [threads] [memory_GB]}"
threads="${2:-8}"
memory="${3:-16}"
python_bin="${PYTHON:-$(command -v python || command -v python3 || true)}"

"${python_bin}" "${here}/get_genomes.py" "${work}"

status=0
for target in multiplex fastidiosa pauca morus; do
    echo
    echo "=============== ${target}"
    "${python_bin}" "${here}/group_genomes.py" "${work}/genomes" "${work}/${target}" --target "${target}"
    primer-finder \
        -i "${work}/${target}/inclusion" \
        -e "${work}/${target}/exclusion" \
        -o "${work}/${target}/results" \
        -t "${threads}" -m "${memory}"
    "${python_bin}" "${here}/check_assay.py" "${work}/${target}/results" --target "${target}" \
        --json "${work}/${target}/check.json" || status=1
done

echo
if [ "${status}" -eq 0 ]; then
    echo "All four published assays were recovered"
else
    echo "At least one assay was not recovered" >&2
fi
exit "${status}"
