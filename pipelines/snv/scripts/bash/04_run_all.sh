#!/usr/bin/env bash
# =============================================================================
# 04_run_all.sh (pipelines/snv)
# Orchestrates the full Small Variant (SNV/Indel) analysis pipeline end-to-end:
#   01_validate_inputs.sh -> 02_vcf_to_tsv.sh -> 03_filter_snv_categories.sh
#   -> 05_annotate_snvs.R -> 06_summary_stats.R -> 07_generate_report.R
#
# Stops immediately on the first failing step (fail-fast), preserving log
# integrity and preventing partial/stale downstream outputs.
#
# Usage:
#   bash scripts/bash/04_run_all.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

RUN_LOG="${LOG_DIR}/04_run_all_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

RUN_START=$(date +%s)
echo "=== [04_run_all.sh] Starting SNV pipeline run for sample: ${SAMPLE_ID} ==="
echo "Config: ${CONFIG_FILE}"
echo "Run log: ${RUN_LOG}"
echo ""

STEP_RESULTS=()

run_step() {
    local step_name="$1"
    shift
    local step_start
    step_start=$(date +%s)

    echo "-----------------------------------------------------------------"
    echo ">>> STEP: ${step_name}"
    echo "-----------------------------------------------------------------"

    if "$@"; then
        local elapsed=$(( $(date +%s) - step_start ))
        echo ">>> STEP OK: ${step_name} (${elapsed}s)"
        STEP_RESULTS+=("OK   ${step_name} (${elapsed}s)")
    else
        local rc=$?
        local elapsed=$(( $(date +%s) - step_start ))
        echo ">>> STEP FAILED: ${step_name} (exit ${rc}, after ${elapsed}s)"
        STEP_RESULTS+=("FAIL ${step_name} (exit ${rc})")
        print_summary
        echo ""
        echo "=== [04_run_all.sh] ABORTED — fix the failure above and re-run."
        exit 1
    fi
    echo ""
}

print_summary() {
    echo ""
    echo "=== Run summary ==="
    for r in "${STEP_RESULTS[@]}"; do
        echo "  ${r}"
    done
}

run_step "01_validate_inputs"       bash "${SCRIPT_DIR}/01_validate_inputs.sh" "${CONFIG_FILE}"
run_step "02_vcf_to_tsv"            bash "${SCRIPT_DIR}/02_vcf_to_tsv.sh" "${CONFIG_FILE}"
run_step "03_filter_snv_categories" bash "${SCRIPT_DIR}/03_filter_snv_categories.sh" "${CONFIG_FILE}"
run_step "05_annotate_snvs"         Rscript "${REPO_DIR}/scripts/R/05_annotate_snvs.R" "${CONFIG_FILE}"
run_step "06_summary_stats"         Rscript "${REPO_DIR}/scripts/R/06_summary_stats.R" "${CONFIG_FILE}"
run_step "07_generate_report"       Rscript "${REPO_DIR}/scripts/R/07_generate_report.R" "${CONFIG_FILE}"

TOTAL_ELAPSED=$(( $(date +%s) - RUN_START ))
print_summary
echo ""
echo "=== [04_run_all.sh] Pipeline completed successfully in ${TOTAL_ELAPSED}s ==="
echo "Results: ${OUTPUT_DIR}"
echo "Run log: ${RUN_LOG}"
