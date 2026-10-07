#!/usr/bin/env bash
# Run primer-finder on the example dataset and check the results.
#
#   bash example/run_example.sh [output_folder] [threads]
#
# It needs primer-finder and its programs (kmc, skesa, minimap2, blast) on PATH; see the wiki's
# Installation page. It takes a few seconds.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${1:-${here}/example_output}"
threads="${2:-4}"
memory="${3:-8}"

python "${here}/make_example.py" "${out}/data"

if command -v primer-finder > /dev/null; then
    finder=(primer-finder)
else  # Running from a clone, without installing
    finder=(python "${here}/../primer_finder.py")
fi

"${finder[@]}" \
    -i "${out}/data/inclusion" \
    -e "${out}/data/exclusion" \
    -o "${out}/results" \
    -t "${threads}" \
    -m "${memory}"

python "${here}/check_example.py" "${out}/results" "${out}/data"
