#!/usr/bin/env bash
# Validate primer-finder against the published Listeria monocytogenes clonal-complex qPCR assays of
# Felix et al. 2023 (doi:10.1128/spectrum.03954-22).
#
#   bash validation/listeria/run.sh <work_folder> [threads] [memory_GB]
#
# It downloads every complete RefSeq L. monocytogenes genome, gives each one a clonal complex from its MLST
# profile, and then for each clonal complex with a panel worth running:
#
#   - builds the inclusion group (that clonal complex) and the exclusion group (all the others),
#   - runs primer-finder,
#   - looks for the published primers and probe of that clonal complex in the regions it reports.
#
# `census.py` is run first, and is the thing to read beside the results: it says, independently of
# primer-finder, which genomes actually hold each published assay. An assay that is missing from a genome of
# its own clonal complex, or present in a genome of another, cannot be reported by a tool whose rule is
# "in every inclusion genome, different from every exclusion genome" -- and should not be.
#
# Needs: primer-finder and its programs, plus mlst and the NCBI datasets tool:
#   conda create -n primer-finder_listeria -c conda-forge -c bioconda mlst ncbi-datasets-cli "perl>=5.32"
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work="${1:?usage: run.sh <work_folder> [threads] [memory_GB]}"
threads="${2:-16}"
memory="${3:-32}"
python_bin="${PYTHON:-$(command -v python || command -v python3 || true)}"

# The clonal complexes whose assay targets the whole complex, and which have at least 15 complete genomes.
# CC1, CC14, CC37 and CC121 are left out on purpose: the paper gives them two assays each, which resolve
# subdivisions *within* the complex, so a complex-wide inclusion group is not what those assays detect.
targets=(CC2 CC3 CC4 CC5 CC6 CC7 CC8 CC9 CC224)

"${python_bin}" "${here}/type_genomes.py" "${work}" -t "${threads}"
"${python_bin}" "${here}/census.py" "${work}" -t "${threads}"

status=0
for target in "${targets[@]}"; do
    echo
    echo "=============== ${target}"
    "${python_bin}" "${here}/group_genomes.py" "${work}/genomes.tsv" "${work}/${target}" --target "${target}"
    primer-finder \
        -i "${work}/${target}/inclusion" \
        -e "${work}/${target}/exclusion" \
        -o "${work}/${target}/results" \
        -t "${threads}" -m "${memory}"
    "${python_bin}" "${here}/check_assay.py" "${work}/${target}/results" --target "${target}" \
        --json "${work}/${target}/check.json" || status=1
done

echo
echo "Done. Read census.tsv beside the per-target check.json: a target whose assay the census shows is"
echo "missing from its own group, or present in another, is one primer-finder is right to leave out."
exit "${status}"
