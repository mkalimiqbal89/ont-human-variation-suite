#!/usr/bin/env bash
# =============================================================================
# 02_vcf_to_tsv.sh (pipelines/cnv)
# Flattens the Spectre CNV VCF into a single analysis-ready TSV file. Extracts
# positional coordinates, SV type, length, copy number, fold change, GT,
# and quality filter flags.
#
# Usage:
#   bash scripts/bash/02_vcf_to_tsv.sh [path/to/pipeline_config.yaml]
#
# Output:
#   data/processed/<sample_id>.cnv_flat.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${1:-}"

LOG_FILE="${LOG_DIR}/02_vcf_to_tsv_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [02_vcf_to_tsv.sh] Flattening CNV VCF for sample: ${SAMPLE_ID} ==="

CNV_VCF="${INPUT_DIR}/${CNV_VCF_NAME}"
if [[ ! -f "${CNV_VCF}" ]]; then
    echo "[FATAL] CNV VCF not found: ${CNV_VCF}"
    echo "        Run scripts/bash/01_validate_inputs.sh first."
    exit 1
fi

mkdir -p "${REPO_DIR}/data/processed" "${WORK_DIR}"
RAW_TMP="${WORK_DIR}/$(basename "${CNV_VCF}" .vcf.gz).raw.tsv"
FINAL_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.cnv_flat.tsv"

# --- Extract raw VCF fields via bcftools query ---
echo "Running bcftools query on CNV VCF..."
bcftools query \
    -f '%CHROM\t%POS\t%ID\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/CN\t[%CN]\t[%FC]\t[%GT]\n' \
    "${CNV_VCF}" > "${RAW_TMP}" || {
        echo "[FATAL] bcftools query failed on ${CNV_VCF}"
        rm -f "${RAW_TMP}"
        exit 1
    }

N_RAW=$(wc -l < "${RAW_TMP}" | tr -d ' ')
echo "  Extracted ${N_RAW} raw records to temp file."

# --- Clean and calculate length / copy number fields using awk ---
{
    echo -e "CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tCN\tFC\tGT\tQUAL\tFILTER\tID"
    awk -F'\t' 'BEGIN{OFS="\t"}
    {
        chrom = $1
        pos = $2
        id = $3
        qual = $4
        filter = $5
        svtype = $6
        end = $7
        svlen = $8
        info_cn = $9
        fmt_cn = $10
        fc = $11
        gt = $12

        if (svtype == "." || svtype == "") {
            svtype = "CNV"
        }

        if (end == "." || end == "") {
            end = pos
        }

        if (svlen == "." || svlen == "" || svlen == "0") {
            if (end >= pos) {
                svlen = end - pos + 1
            } else {
                svlen = pos - end + 1
            }
        }
        if (svlen < 0) svlen = -svlen

        cn = "NA"
        if (fmt_cn != "." && fmt_cn != "") {
            cn = fmt_cn
        } else if (info_cn != "." && info_cn != "") {
            cn = info_cn
        }

        if (fc == "" || fc == ".") fc = "NA"
        if (gt == "" || gt == ".") gt = "NA"
        if (qual == "" || qual == ".") qual = "0"
        if (filter == "" || filter == ".") filter = "PASS"

        print chrom, pos, end, svlen, svtype, cn, fc, gt, qual, filter, id
    }' "${RAW_TMP}"
} > "${FINAL_TSV}"

rm -f "${RAW_TMP}"

N_FINAL=$(( $(wc -l < "${FINAL_TSV}" | tr -d ' ') - 1 ))
echo ""
echo "=== [02_vcf_to_tsv.sh] Done. ${N_FINAL} CNV records written to: ${FINAL_TSV} ==="
