#!/usr/bin/env bash
# =============================================================================
# 02_vcf_to_tsv.sh
# Flattens the Sniffles2 SV VCF into a single analysis-ready TSV. Extracts
# only scalar/summary fields (never raw ALT sequence, RNAMES, or full ANN
# blocks — those can be enormous for long insertions and are not needed
# downstream). The first SnpEff Gene_Name per record is pulled out of ANN
# into its own column since it is used for gene-fusion/disruption tagging
# in scripts/R/05_annotate_variants.R.
#
# Usage:
#   bash scripts/bash/02_vcf_to_tsv.sh [path/to/pipeline_config.yaml]
#
# Output:
#   data/processed/<sample_id>.sv_flat.tsv
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/02_vcf_to_tsv_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [02_vcf_to_tsv.sh] Flattening SV VCF for sample: ${SAMPLE_ID} ==="

SV_VCF="${INPUT_DIR}/${SV_VCF_NAME}"
if [[ ! -f "${SV_VCF}" ]]; then
    echo "[FATAL] SV VCF not found: ${SV_VCF}"
    echo "        Run scripts/bash/01_validate_inputs.sh first."
    exit 1
fi

mkdir -p "${REPO_DIR}/data/processed"
RAW_TMP="${WORK_DIR}/$(basename "${SV_VCF}" .vcf.gz).raw_with_ann.tsv"
FINAL_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.sv_flat.tsv"
mkdir -p "${WORK_DIR}"

# --- Step 1: bcftools query -> scalar fields + full ANN (temp, discarded) ---
echo "Running bcftools query..."
bcftools query \
    -f '%CHROM\t%POS\t%ID\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END\t%INFO/CHR2\t%INFO/SUPPORT\t%INFO/STRAND\t%INFO/VAF\t[%GT]\t[%DR]\t[%DV]\t%INFO/ANN\n' \
    "${SV_VCF}" > "${RAW_TMP}"

N_RAW=$(wc -l < "${RAW_TMP}")
echo "  Extracted ${N_RAW} raw records to temp file."

# --- Step 2: derive GENE_SYMBOL from first ANN block, drop full ANN --------
echo "Deriving GENE_SYMBOL and ANNOTATION_IMPACT from ANN (discarding raw ANN)..."
{
    echo -e "CHROM\tPOS\tID\tQUAL\tFILTER\tSVTYPE\tSVLEN\tEND\tCHR2\tSUPPORT\tSTRAND\tVAF\tGT\tDR\tDV\tGENE_SYMBOL\tANNOTATION_IMPACT"
    awk -F'\t' 'BEGIN{OFS="\t"}
    {
        ann = $16
        gene = "NA"; impact = "NA"
        if (ann != "." && ann != "") {
            n = split(ann, blocks, ",")
            m = split(blocks[1], f, "|")
            if (m >= 4) {
                if (f[4] != "") gene = f[4]
                if (f[3] != "") impact = f[3]
            }
        }
        print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,gene,impact
    }' "${RAW_TMP}"
} > "${FINAL_TSV}"

rm -f "${RAW_TMP}"

N_FINAL=$(( $(wc -l < "${FINAL_TSV}") - 1 ))
echo ""
echo "=== [02_vcf_to_tsv.sh] Done. ${N_FINAL} SV records written to: ${FINAL_TSV} ==="
echo "Preview:"
head -3 "${FINAL_TSV}" | cut -c1-220
