#!/usr/bin/env bash
# =============================================================================
# 03_filter_snv_categories.sh (pipelines/snv)
# Filters flat SNV/Indel records by minimum depth (DP), quality (QUAL), variant
# allele frequency (VAF), and FILTER status. Categorizes calls into SNVs and
# Indels. Extracts ClinVar Pathogenic/Likely Pathogenic variants if available.
#
# Usage:
#   bash scripts/bash/03_filter_snv_categories.sh [path/to/pipeline_config.yaml]
#
# Outputs:
#   results/snvs/<sample_id>.snvs.tsv
#   results/indels/<sample_id>.indels.tsv
#   results/clinvar/<sample_id>.clinvar_pathogenic.tsv
#   results/qc_summary/<sample_id>.filtering_summary.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${1:-}"

LOG_FILE="${LOG_DIR}/03_filter_snv_categories_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [03_filter_snv_categories.sh] Filtering & categorizing small variants for: ${SAMPLE_ID} ==="
echo "Filters: min_dp=${MIN_DP}  min_qual=${MIN_QUAL}  min_vaf=${MIN_VAF}  pass_only=${PASS_ONLY}"

FLAT_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.snv_flat.tsv"
if [[ ! -f "${FLAT_TSV}" ]]; then
    echo "[FATAL] Input flat TSV not found: ${FLAT_TSV}"
    echo "        Run scripts/bash/02_vcf_to_tsv.sh first."
    exit 1
fi

SNV_FILE="${OUTPUT_DIR}/snvs/${SAMPLE_ID}.snvs.tsv"
INDEL_FILE="${OUTPUT_DIR}/indels/${SAMPLE_ID}.indels.tsv"
QC_FILE="${OUTPUT_DIR}/qc_summary/${SAMPLE_ID}.filtering_summary.tsv"
CLINVAR_FILE="${OUTPUT_DIR}/clinvar/${SAMPLE_ID}.clinvar_pathogenic.tsv"

HEADER="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tGT\tDP\tVAF\tAD\tTYPE"
echo -e "${HEADER}" > "${SNV_FILE}"
echo -e "${HEADER}" > "${INDEL_FILE}"

# AWK streaming filter
awk -v min_dp="${MIN_DP}" \
    -v min_qual="${MIN_QUAL}" \
    -v min_vaf="${MIN_VAF}" \
    -v pass_only="${PASS_ONLY}" \
    -v snv_out="${SNV_FILE}" \
    -v indel_out="${INDEL_FILE}" \
    -v qc_out="${QC_FILE}" \
    -F'\t' 'BEGIN {
    total = 0
    passed = 0
    fail_dp = 0
    fail_qual = 0
    fail_vaf = 0
    fail_filter = 0
    n_snv = 0
    n_indel = 0
}
NR == 1 { next }
{
    total++
    ref = $4
    alt = $5
    qual = $6 + 0
    filter = $7
    dp = $9 + 0
    vaf = $10 + 0
    var_type = toupper($12)

    # Filter checks
    if (dp < min_dp) {
        fail_dp++
        next
    }
    if (qual < min_qual) {
        fail_qual++
        next
    }
    if (vaf < min_vaf) {
        fail_vaf++
        next
    }
    if (pass_only == "true" && filter != "PASS" && filter != ".") {
        fail_filter++
        next
    }

    passed++

    is_snv = 0
    if (var_type ~ /SNP/ || var_type ~ /SNV/) {
        is_snv = 1
    } else if (length(ref) == 1 && length(alt) == 1) {
        is_snv = 1
    }

    if (is_snv) {
        n_snv++
        print $0 >> snv_out
    } else {
        n_indel++
        print $0 >> indel_out
    }
}
END {
    print "metric\tcount" > qc_out
    print "total_raw_records\t" total >> qc_out
    print "passed_filter_records\t" passed >> qc_out
    print "failed_dp_filter\t" fail_dp >> qc_out
    print "failed_qual_filter\t" fail_qual >> qc_out
    print "failed_vaf_filter\t" fail_vaf >> qc_out
    print "failed_status_filter\t" fail_filter >> qc_out
    print "snvs_count\t" n_snv >> qc_out
    print "indels_count\t" n_indel >> qc_out

    print "  Total raw input variants : " total
    print "  Passing filters          : " passed
    print "    SNVs (substitutions)   : " n_snv
    print "    Indels (ins/del)       : " n_indel
}' "${FLAT_TSV}"

# Optional ClinVar pathogenic variant extraction
CLINVAR_VCF="${INPUT_DIR}/${SNV_VCF_CLINVAR_NAME}"
echo -e "CHROM\tPOS\tID\tREF\tALT\tCLNSIG\tCLNDN" > "${CLINVAR_FILE}"
if [[ -f "${CLINVAR_VCF}" && -s "${CLINVAR_VCF}" ]]; then
    echo "Extracting ClinVar Pathogenic variants from ${CLINVAR_VCF}..."
    bcftools query \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\n' \
        "${CLINVAR_VCF}" 2>/dev/null | \
        awk -F'\t' 'BEGIN{OFS="\t"} {
            sig = toupper($6)
            if (sig ~ /PATHOGENIC/ || sig ~ /LIKELY_PATHOGENIC/) {
                print $0
            }
        }' >> "${CLINVAR_FILE}" || true
    N_CLINVAR=$(( $(wc -l < "${CLINVAR_FILE}" | tr -d ' ') - 1 ))
    echo "  Extracted ${N_CLINVAR} ClinVar pathogenic variants."
else
    echo "  ClinVar VCF not provided or empty; skipping ClinVar extraction."
fi

echo ""
echo "=== [03_filter_snv_categories.sh] Done. Summary: ${QC_FILE} ==="
