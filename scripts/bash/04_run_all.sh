#!/usr/bin/env bash
# =============================================================================
# 04_run_all.sh
# Orchestrates the full SV categorization pipeline end-to-end:
#   01_validate_inputs.sh -> 02_vcf_to_tsv.sh -> 03_filter_sv_categories.sh
#   -> 05_annotate_variants.R
#
# Stops immediately on the first failing step (set -e), so a broken run
# never silently produces partial/stale downstream outputs. Each step's own
# log file (already written by that script into logs/) is preserved; this
# wrapper additionally writes a single top-level run log summarizing what
# ran, how long it took, and pass/fail per step.
#
# Usage:
#   bash scripts/bash/04_run_all.sh [path/to/pipeline_config.yaml]
#
# Exit codes:
#   0 = all steps completed
#   1 = a step failed (see the printed step name and its own log for detail)
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

# Source once up front purely to resolve LOG_DIR/SAMPLE_ID for this wrapper's
# own log file; each sub-script re-sources 00_setup_env.sh independently, so
# this pipeline is safe to re-run any single step manually later without
# going through this wrapper.
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

RUN_LOG="${LOG_DIR}/04_run_all_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

RUN_START=$(date +%s)
echo "=== [04_run_all.sh] Starting full pipeline run for sample: ${SAMPLE_ID} ==="
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
        echo "=== [04_run_all.sh] ABORTED — fix the failure above, then either"
        echo "    re-run this wrapper (safe, re-runs everything) or re-run"
        echo "    just the failed step's script directly (see docs/TROUBLESHOOTING.md)."
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

# NOTE: archiving to internal storage (08_archive_results.sh) is deliberately
# NOT included in this automated chain. Archiving copies real results
# (including real sample identifiers) to institutional storage and should be
# a conscious, separate action, not something that runs silently as part of
# every pipeline execution -- run it explicitly when you're ready to keep a
# permanent record of this run:
#   bash scripts/bash/08_archive_results.sh
run_step "01_validate_inputs"       bash "${SCRIPT_DIR}/01_validate_inputs.sh" "${CONFIG_FILE}"
run_step "02_vcf_to_tsv"            bash "${SCRIPT_DIR}/02_vcf_to_tsv.sh" "${CONFIG_FILE}"
run_step "03_filter_sv_categories"  bash "${SCRIPT_DIR}/03_filter_sv_categories.sh" "${CONFIG_FILE}"
run_step "05_annotate_variants"     Rscript "${REPO_DIR}/scripts/R/05_annotate_variants.R" "${CONFIG_FILE}"
run_step "06_summary_stats"         Rscript "${REPO_DIR}/scripts/R/06_summary_stats.R" "${CONFIG_FILE}"
run_step "07_generate_report"       Rscript "${REPO_DIR}/scripts/R/07_generate_report.R" "${CONFIG_FILE}"

TOTAL_ELAPSED=$(( $(date +%s) - RUN_START ))
print_summary
echo ""
echo "=== [04_run_all.sh] Pipeline completed successfully in ${TOTAL_ELAPSED}s ==="
echo "Results: ${OUTPUT_DIR}"
echo "Run log: ${RUN_LOG}"
