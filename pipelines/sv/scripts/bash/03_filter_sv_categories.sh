#!/usr/bin/env bash
# =============================================================================
# 03_filter_sv_categories.sh
# Applies QC filters (config: filtering.*) to the flattened SV TSV, then
# splits PASSing records into per-category files under results/<category>/.
#
# Category mapping comes from config: sv_categories.* (SVTYPE lists).
# gene_fusions is NOT populated here — see scripts/R/05_annotate_variants.R.
# This script does tag BND records with a non-NA GENE_SYMBOL as
# candidate_fusion=TRUE in an extra column, for the R step to pair up.
#
# Usage:
#   bash scripts/bash/03_filter_sv_categories.sh [path/to/pipeline_config.yaml]
#
# Input:  data/processed/<sample_id>.sv_flat.tsv   (from 02_vcf_to_tsv.sh)
# Output: results/<category>/<sample_id>.<category>.tsv
#         results/qc_summary/<sample_id>.filtering_summary.tsv
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/03_filter_sv_categories_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [03_filter_sv_categories.sh] Filtering + categorizing SVs for: ${SAMPLE_ID} ==="

FLAT_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.sv_flat.tsv"
if [[ ! -f "${FLAT_TSV}" ]]; then
    echo "[FATAL] Flattened TSV not found: ${FLAT_TSV}"
    echo "        Run scripts/bash/02_vcf_to_tsv.sh first."
    exit 1
fi

echo "Filters: min_sv_length=${MIN_SV_LENGTH}  min_read_support=${MIN_READ_SUPPORT}  min_qual=${MIN_QUAL}  min_vaf=${MIN_VAF}  pass_only=true"

FILTERED_TSV="${WORK_DIR}/${SAMPLE_ID}.sv_filtered.tsv"
mkdir -p "${WORK_DIR}"

awk -F'\t' -v OFS='\t' \
    -v min_len="${MIN_SV_LENGTH}" \
    -v min_support="${MIN_READ_SUPPORT}" \
    -v min_qual="${MIN_QUAL}" \
    -v min_vaf="${MIN_VAF}" \
'
NR==1 { print $0, "candidate_fusion"; next }
{
    filter=$5; svlen=$7; support=$10; qual=$4; vaf=$12; svtype=$6
    chrom=$1; chr2=$9; gene=$16

    abs_len = (svlen < 0 ? -svlen : svlen)

    if (filter != "PASS") next
    if (qual < min_qual) next
    if (support < min_support) next
    if (vaf != "." && vaf+0 < min_vaf) next
    if (svlen != "." && abs_len < min_len) next

    candidate_fusion = "FALSE"
    if (svtype == "BND" && gene != "NA" && chr2 != "." && chr2 != chrom) {
        candidate_fusion = "TRUE"
    }
    print $0, candidate_fusion
}
' "${FLAT_TSV}" > "${FILTERED_TSV}"

N_INPUT=$(( $(wc -l < "${FLAT_TSV}") - 1 ))
N_PASS=$(( $(wc -l < "${FILTERED_TSV}") - 1 ))
echo "Input records: ${N_INPUT}  ->  Passing filters: ${N_PASS}"

# A list plus a case function, NOT an associative array.
#
# Associative arrays are a bash 4 feature, and macOS ships bash 3.2 as
# /bin/bash. On 3.2 `declare -A` is not supported, so `[deletions]="DEL"` is
# parsed as an INDEXED subscript and evaluated arithmetically; bash then tries
# to resolve `deletions` as a variable and `set -u` aborts the script with
# "deletions: unbound variable".
#
# The failure was near-silent: this stage exited 0 having written no category
# files at all, and stage 05 failed later on the missing translocations file.
# A stage that reports success while producing nothing is the exact failure
# class the testing philosophy in CONTRIBUTING.md is written against.
#
# Iterating a fixed list also makes the output order deterministic, which
# associative-array iteration was not.
SV_CATEGORY_LIST="deletions insertions duplications inversions translocations complex_rearrangements"

svtypes_for_category() {
    case "$1" in
        deletions)              echo "DEL" ;;
        insertions)             echo "INS" ;;
        duplications)           echo "DUP" ;;
        inversions)             echo "INV" ;;
        translocations)         echo "BND TRA" ;;
        complex_rearrangements) echo "CPX BND_CLUSTER" ;;
        *)                      echo "" ;;
    esac
}

HEADER=$(head -1 "${FILTERED_TSV}")
QC_SUMMARY="${OUTPUT_DIR}/qc_summary/${SAMPLE_ID}.filtering_summary.tsv"
mkdir -p "${OUTPUT_DIR}/qc_summary"
echo -e "category\tsvtypes_included\tn_variants" > "${QC_SUMMARY}"

for category in ${SV_CATEGORY_LIST}; do
    svtypes="$(svtypes_for_category "${category}")"
    out_file="${OUTPUT_DIR}/${category}/${SAMPLE_ID}.${category}.tsv"
    mkdir -p "${OUTPUT_DIR}/${category}"

    pattern=$(echo "${svtypes}" | tr ' ' '|')

    echo "${HEADER}" > "${out_file}"
    awk -F'\t' -v OFS='\t' -v pat="^(${pattern})\$" '$6 ~ pat' "${FILTERED_TSV}" >> "${out_file}"

    n=$(( $(wc -l < "${out_file}") - 1 ))
    echo "  ${category} (${svtypes}): ${n} variants -> ${out_file}"
    echo -e "${category}\t${svtypes}\t${n}" >> "${QC_SUMMARY}"
done

# --- The loop must actually have produced its outputs ------------------------
# Without this, a failure inside the loop leaves the stage exiting 0 having
# written nothing, and the problem only surfaces two stages later as
# "translocations.tsv not found" — which points at the wrong script. That is
# exactly what happened when `declare -A` failed silently on bash 3.2.
missing_categories=""
for category in ${SV_CATEGORY_LIST}; do
    [[ -f "${OUTPUT_DIR}/${category}/${SAMPLE_ID}.${category}.tsv" ]] \
        || missing_categories="${missing_categories} ${category}"
done
if [[ -n "${missing_categories}" ]]; then
    echo "" >&2
    echo "[FATAL] No output file written for:${missing_categories}" >&2
    echo "        The categorization loop did not run to completion. This stage" >&2
    echo "        cannot succeed without these files, so it fails here rather" >&2
    echo "        than letting a later stage report a confusing missing input." >&2
    exit 1
fi

FUSION_CANDIDATES="${OUTPUT_DIR}/gene_fusions/${SAMPLE_ID}.gene_fusion_candidates.tsv"
mkdir -p "${OUTPUT_DIR}/gene_fusions"
echo "${HEADER}" > "${FUSION_CANDIDATES}"
awk -F'\t' -v OFS='\t' '$NF == "TRUE"' "${FILTERED_TSV}" >> "${FUSION_CANDIDATES}"
n_fusion=$(( $(wc -l < "${FUSION_CANDIDATES}") - 1 ))
echo "  gene_fusion candidates (unpaired BND breakpoints hitting a gene): ${n_fusion} -> ${FUSION_CANDIDATES}"
echo "  NOTE: these are single-breakpoint candidates only. True fusion pairs"
echo "        (both breakends resolved to gene partners) are built in"
echo "        scripts/R/05_annotate_variants.R — this file is its input."
echo -e "gene_fusion_candidates\tBND (unpaired, pre-R-pairing)\t${n_fusion}" >> "${QC_SUMMARY}"

echo ""
echo "=== [03_filter_sv_categories.sh] Done. QC summary: ${QC_SUMMARY} ==="
