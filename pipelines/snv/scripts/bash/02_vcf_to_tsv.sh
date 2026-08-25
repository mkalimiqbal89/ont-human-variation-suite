#!/usr/bin/env bash
# =============================================================================
# 02_vcf_to_tsv.sh (pipelines/snv)
# Flattens Clair3 SNV/Indel VCF records into a single analysis-ready TSV file
# using streaming bcftools query extraction.
#
# Usage:
#   bash scripts/bash/02_vcf_to_tsv.sh [path/to/pipeline_config.yaml]
#
# Output:
#   data/processed/<sample_id>.snv_flat.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${1:-}"

LOG_FILE="${LOG_DIR}/02_vcf_to_tsv_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [02_vcf_to_tsv.sh] Flattening Clair3 SNV VCF for sample: ${SAMPLE_ID} ==="

SNV_VCF="${INPUT_DIR}/${SNV_VCF_NAME}"
if [[ ! -f "${SNV_VCF}" ]]; then
    echo "[FATAL] SNV VCF not found: ${SNV_VCF}"
    echo "        Run scripts/bash/01_validate_inputs.sh first."
    exit 1
fi

mkdir -p "${REPO_DIR}/data/processed" "${WORK_DIR}"
RAW_TMP="${WORK_DIR}/$(basename "${SNV_VCF}" .vcf.gz).raw.tsv"
FINAL_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.snv_flat.tsv"

# --- Extract raw VCF fields via bcftools query ---
echo "Running bcftools query on SNV VCF..."
bcftools query \
    -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT]\t[%DP]\t[%AF]\t[%AD]\t%TYPE\n' \
    "${SNV_VCF}" > "${RAW_TMP}" || {
        echo "[FATAL] bcftools query failed on ${SNV_VCF}"
        rm -f "${RAW_TMP}"
        exit 1
    }

N_RAW=$(wc -l < "${RAW_TMP}" | tr -d ' ')
echo "  Extracted ${N_RAW} raw records to temp file."

# --- Format output with header and calculate VAF if needed ---
{
    echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tVAF\tAD\tTYPE"
    awk -F'\t' 'BEGIN{OFS="\t"}
    {
        chrom = $1
        pos = $2
        id = $3
        ref = $4
        alt = $5
        qual = $6
        filter = $7
        gt = $8
        dp = $9
        af = $10
        ad = $11
        var_type = $12

        if (id == "" || id == ".") id = "NA"
        if (qual == "" || qual == ".") qual = "0"
        if (filter == "" || filter == ".") filter = "PASS"
        if (gt == "" || gt == ".") gt = "NA"
        if (dp == "" || dp == ".") dp = "0"

        # Calculate VAF from AD (alt_count / total_depth) if AF missing
        vaf = "NA"
        if (af != "" && af != ".") {
            vaf = af
        } else if (ad != "" && ad != "." && dp + 0 > 0) {
            n_ad = split(ad, ad_arr, ",")
            if (n_ad >= 2) {
                vaf = sprintf("%.4f", ad_arr[2] / dp)
            }
        }
        if (vaf == "NA" || vaf == "") vaf = "0.0000"

        # Infer TYPE if missing
        if (var_type == "" || var_type == ".") {
            if (length(ref) == 1 && length(alt) == 1) {
                var_type = "SNP"
            } else {
                var_type = "INDEL"
            }
        }

        print chrom, pos, id, ref, alt, qual, filter, gt, dp, vaf, ad, var_type
    }' "${RAW_TMP}"
} > "${FINAL_TSV}"

rm -f "${RAW_TMP}"

N_FINAL=$(( $(wc -l < "${FINAL_TSV}" | tr -d ' ') - 1 ))
echo ""
echo "=== [02_vcf_to_tsv.sh] Done. ${N_FINAL} SNV/Indel records written to: ${FINAL_TSV} ==="
