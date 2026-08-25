#!/usr/bin/env bash
# =============================================================================
# 03_filter_cnv_categories.sh (pipelines/cnv)
# Filters raw flattened CNV records by minimum length, quality, FILTER status,
# and copy number thresholds. Splits QC-passing calls into deletions and
# duplications/gains tables, and outputs a filtering summary TSV.
#
# Usage:
#   bash scripts/bash/03_filter_cnv_categories.sh [path/to/pipeline_config.yaml]
#
# Outputs:
#   results/deletions/<sample_id>.deletions.tsv
#   results/duplications/<sample_id>.duplications.tsv
#   results/qc_summary/<sample_id>.filtering_summary.tsv
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${1:-}"

LOG_FILE="${LOG_DIR}/03_filter_cnv_categories_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [03_filter_cnv_categories.sh] Filtering & categorizing CNVs for: ${SAMPLE_ID} ==="
echo "Filters: min_cnv_length=${MIN_CNV_LENGTH}  min_qual=${MIN_QUAL}  pass_only=${PASS_ONLY}  max_cn_loss=${MAX_CN_LOSS}  min_cn_gain=${MIN_CN_GAIN}"

FLAT_TSV="${REPO_DIR}/data/processed/${SAMPLE_ID}.cnv_flat.tsv"
if [[ ! -f "${FLAT_TSV}" ]]; then
    echo "[FATAL] Input TSV not found: ${FLAT_TSV}"
    echo "        Run scripts/bash/02_vcf_to_tsv.sh first."
    exit 1
fi

DEL_FILE="${OUTPUT_DIR}/deletions/${SAMPLE_ID}.deletions.tsv"
DUP_FILE="${OUTPUT_DIR}/duplications/${SAMPLE_ID}.duplications.tsv"
QC_FILE="${OUTPUT_DIR}/qc_summary/${SAMPLE_ID}.filtering_summary.tsv"

# Prepare headers
HEADER="CHROM\tPOS\tEND\tSVLEN\tSVTYPE\tCN\tFC\tGT\tQUAL\tFILTER\tID"
echo -e "${HEADER}" > "${DEL_FILE}"
echo -e "${HEADER}" > "${DUP_FILE}"

# AWK filtering script
awk -v min_len="${MIN_CNV_LENGTH}" \
    -v min_qual="${MIN_QUAL}" \
    -v pass_only="${PASS_ONLY}" \
    -v max_loss="${MAX_CN_LOSS}" \
    -v min_gain="${MIN_CN_GAIN}" \
    -v del_out="${DEL_FILE}" \
    -v dup_out="${DUP_FILE}" \
    -v qc_out="${QC_FILE}" \
    -F'\t' 'BEGIN {
    total = 0
    passed = 0
    fail_len = 0
    fail_qual = 0
    fail_filter = 0
    n_del = 0
    n_dup = 0
}
NR == 1 { next }
{
    total++
    chrom = $1
    pos = $2 + 0
    end = $3 + 0
    svlen = $4 + 0
    svtype = toupper($5)
    cn_str = $6
    qual = $9 + 0
    filter = $10

    # Filter checks
    if (svlen < min_len) {
        fail_len++
        next
    }
    if (qual < min_qual) {
        fail_qual++
        next
    }
    if (pass_only == "true" && filter != "PASS" && filter != ".") {
        fail_filter++
        next
    }

    passed++

    cn_val = -1
    if (cn_str != "NA" && cn_str != "." && cn_str != "") {
        cn_val = cn_str + 0
    }

    # Categorize as DEL or DUP
    is_del = 0
    is_dup = 0

    if (svtype ~ /DEL/ || svtype ~ /LOSS/) {
        is_del = 1
    } else if (svtype ~ /DUP/ || svtype ~ /GAIN/ || svtype ~ /AMP/) {
        is_dup = 1
    } else if (cn_val >= 0) {
        if (cn_val <= max_loss) is_del = 1
        else if (cn_val >= min_gain) is_dup = 1
    }

    if (is_del) {
        n_del++
        print $0 >> del_out
    } else if (is_dup) {
        n_dup++
        print $0 >> dup_out
    } else {
        # Default to DEL if CN < 2, DUP if CN > 2, otherwise DEL
        if (cn_val >= 0 && cn_val < 2) {
            n_del++
            print $0 >> del_out
        } else {
            n_dup++
            print $0 >> dup_out
        }
    }
}
END {
    print "metric\tcount" > qc_out
    print "total_raw_records\t" total >> qc_out
    print "passed_filter_records\t" passed >> qc_out
    print "failed_length_filter\t" fail_len >> qc_out
    print "failed_qual_filter\t" fail_qual >> qc_out
    print "failed_status_filter\t" fail_filter >> qc_out
    print "deletions_count\t" n_del >> qc_out
    print "duplications_count\t" n_dup >> qc_out

    print "  Total input records: " total
    print "  Passing filters    : " passed
    print "    Deletions (loss) : " n_del
    print "    Duplications (gain): " n_dup
}' "${FLAT_TSV}"

echo ""
echo "=== [03_filter_cnv_categories.sh] Done. Summary: ${QC_FILE} ==="
