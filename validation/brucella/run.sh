#!/usr/bin/env bash
# Validate primer-finder against the published Brucella qPCR assays of Hoover et al. 2026
# (Microbiol Spectr 14:e02974-25, doi:10.1128/spectrum.02974-25): one genus-level singleplex and a
# three-species multiplex, whose regions were found with selective_primer_finder 0.2.1 -- the direct
# ancestor of this package.
#
#   bash validation/brucella/run.sh <work_folder> [threads] [memory_GB]
#
# It rebuilds the panel from RefSeq (the paper's accessions are not published; see panel.py), verifies the
# published amplicon lengths against a real genome, takes an independent blast census of where each assay
# actually occurs, and then for each of the four targets:
#
#   - builds the inclusion group (the species; the paper's eight for `genus`; all twelve classical-clade
#     species for `genus_wide`) and the exclusion group (the other classical species of the paper's eight
#     for a species target, and the non-classical relatives in every case),
#   - runs primer-finder,
#   - looks for the published primers and probe in the regions it reports.
#
# `census.tsv` is the thing to read beside the results: it says, independently of primer-finder, which
# genomes hold each published assay and what `-p` a run would need to report it at all. An assay missing
# from a genome of its own group, or present in one of the exclusion group, cannot be reported by a tool
# whose rule is "in every inclusion genome, different from every exclusion genome" -- and should not be.
#
# Needs: primer-finder and its programs, plus the NCBI datasets tool:
#   conda install -c conda-forge ncbi-datasets-cli
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
work="${1:?usage: run.sh <work_folder> [threads] [memory_GB]}"
threads="${2:-16}"
memory="${3:-32}"
python_bin="${PYTHON:-$(command -v python || command -v python3 || true)}"

# Every target `group_genomes.py` can build. Only those with a row in assays.tsv can be checked against a
# published assay; `genus_wide` and `suis_lato` are this project's own, so they are run and reported but
# not graded.
targets=(genus genus_wide abortus melitensis suis suis_lato)
published=(genus abortus melitensis suis)

"${python_bin}" "${here}/panel.py" "${work}"
"${python_bin}" "${here}/verify_amplicon.py" "${work}"
"${python_bin}" "${here}/census.py" "${work}" -t "${threads}"

status=0
for target in "${targets[@]}"; do
    echo
    echo "=============== ${target}"
    "${python_bin}" "${here}/group_genomes.py" "${work}/panel.tsv" "${work}/${target}" --target "${target}"
    # Guarded: `set -e` would otherwise abandon every target after one primer-finder cannot finish, which
    # is a thing that happens -- the record's own runs.tsv has a target that failed at -p 1.0.
    if ! primer-finder \
        -i "${work}/${target}/inclusion" \
        -e "${work}/${target}/exclusion" \
        -o "${work}/${target}/results" \
        -t "${threads}" -m "${memory}"; then
        echo "primer-finder found nothing for ${target}; carrying on with the rest"
        status=1
        continue
    fi
    if [[ " ${published[*]} " == *" ${target} "* ]]; then
        "${python_bin}" "${here}/check_assay.py" "${work}/${target}/results" --target "${target}" \
            --json "${work}/${target}/check.json" || status=1
    else
        echo "${target} has no published assay to check against; see its results folder"
    fi
done

echo
echo "Done. Read census.tsv beside the per-target check.json: a target whose assay the census shows is"
echo "missing from its own group, or present in the exclusion group, is one primer-finder is right to"
echo "leave out. The p_needed column is the share of the group whose genomes hold all three oligos, so it"
echo "is an upper bound on what -p could achieve -- not a setting that would be enough, since a kmer-sized"
echo "window has to be intact as well, which is a stricter thing."
exit "${status}"
