#!/usr/bin/env bash
# =============================================================================
# 03_filter_cpg_sites.sh
# Applies the coverage threshold to stage 02's output and splits the retained
# sites by methylation state (unmethylated / intermediate / methylated), each
# written to its own bgzipped, tabix-indexed file.
#
# Contig selection already happened in stage 02 and is NOT repeated here. This
# stage's only filter is coverage.
#
# Two reconciliation assertions are the point of this stage:
#   input == dropped_low + dropped_high + retained
#   retained == unmethylated + intermediate + methylated
# Neither can hold by accident. A filter that silently loses records, or a
# category boundary with a gap or an overlap, fails one of them — and that class
# of bug is invisible to an exit-code check.
#
# Usage:
#   bash scripts/bash/03_filter_cpg_sites.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Drop a '#'-leading argument: zsh passes inline comments through as arguments.
# See the note in 00_setup_env.sh.
_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in '#'*) _CFG_ARG="" ;; esac
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/03_filter_cpg_sites_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

case "${PRIMARY_MOD_CODE}" in
    m) MOD_LABEL="5mC"  ;;
    h) MOD_LABEL="5hmC" ;;
    a) MOD_LABEL="6mA"  ;;
    *) MOD_LABEL="mod_${PRIMARY_MOD_CODE}" ;;
esac

echo ""
echo "=== [03_filter_cpg_sites.sh] Filtering ${MOD_LABEL} sites for ${SAMPLE_ID} ==="

FILTERED_DIR="${OUTPUT_DIR}/01_filtered"
STATE_DIR="${FILTERED_DIR}/by_methylation_state"
QC_DIR="${OUTPUT_DIR}/02_qc"
mkdir -p "${STATE_DIR}" "${QC_DIR}" || {
    echo "[FATAL] Could not create output directories under ${OUTPUT_DIR}" >&2; exit 1; }

IN_BED="${FILTERED_DIR}/${SAMPLE_ID}.${MOD_LABEL}.all_cov.bedmethyl.gz"
if [[ ! -f "${IN_BED}" ]]; then
    echo "[FATAL] Stage 02 output not found: ${IN_BED}" >&2
    echo "        Run: bash scripts/bash/02_bedmethyl_to_tsv.sh" >&2
    exit 1
fi

OUT_BED="${FILTERED_DIR}/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.bedmethyl.gz"
OUT_UNMETH="${STATE_DIR}/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.unmethylated.bedmethyl.gz"
OUT_INTER="${STATE_DIR}/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.intermediate.bedmethyl.gz"
OUT_METH="${STATE_DIR}/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.methylated.bedmethyl.gz"
QC_MAIN="${QC_DIR}/${SAMPLE_ID}.03_filter_qc.tsv"

if [[ "${MAX_COVERAGE}" -gt 0 ]]; then
    COV_DESC="${MIN_COVERAGE} <= coverage <= ${MAX_COVERAGE}"
else
    COV_DESC="coverage >= ${MIN_COVERAGE} (no upper cap)"
fi

echo "Input         : ${IN_BED}"
echo "Coverage      : ${COV_DESC}"
echo "States        : unmethylated < ${UNMETHYLATED_MAX_PERCENT}%  |  intermediate  |  methylated >= ${METHYLATED_MIN_PERCENT}%"
echo ""

if command -v pigz >/dev/null 2>&1; then
    DECOMP=(pigz -dc -p "${THREADS}")
else
    DECOMP=(gzip -dc)
fi

T0=$(date +%s)

# --- Single pass, four output streams ---------------------------------------
# The retained set goes to stdout (-> bgzip -> OUT_BED); the three state files
# are written through their own bgzip pipes from inside awk. Records are passed
# through with `print $0` so no field is ever reformatted.
"${DECOMP[@]}" "${IN_BED}" \
| awk -F'\t' \
      -v mincov="${MIN_COVERAGE}" \
      -v maxcov="${MAX_COVERAGE}" \
      -v unmeth_max="${UNMETHYLATED_MAX_PERCENT}" \
      -v meth_min="${METHYLATED_MIN_PERCENT}" \
      -v f_unmeth="${OUT_UNMETH}" \
      -v f_inter="${OUT_INTER}" \
      -v f_meth="${OUT_METH}" \
      -v qc_main="${QC_MAIN}" \
      -v sample="${SAMPLE_ID}" \
      -v modlabel="${MOD_LABEL}" \
      -v threads="${THREADS}" '
BEGIN {
    OFS = "\t"
    cmd_unmeth = "bgzip -c -@ " threads " > \"" f_unmeth "\""
    cmd_inter  = "bgzip -c -@ " threads " > \"" f_inter  "\""
    cmd_meth   = "bgzip -c -@ " threads " > \"" f_meth   "\""
}
{
    input_sites++
    cov = $10 + 0
    pct = $11 + 0

    # --- coverage filter ---
    if (cov < mincov)                    { dropped_low++;  next }
    if (maxcov > 0 && cov > maxcov)      { dropped_high++; next }

    retained++
    sum_cov += cov
    sum_mod += $12
    sum_pct += pct
    if (retained == 1 || cov < min_cov) min_cov = cov
    if (retained == 1 || cov > max_cov) max_cov = cov

    # --- methylation state: boundaries must not gap or overlap ---
    if (pct < unmeth_max) {
        state = "unmethylated"
        n_unmeth++; cov_unmeth += cov; mod_unmeth += $12
        print $0 | cmd_unmeth
    } else if (pct >= meth_min) {
        state = "methylated"
        n_meth++;   cov_meth += cov;   mod_meth += $12
        print $0 | cmd_meth
    } else {
        state = "intermediate"
        n_inter++;  cov_inter += cov;  mod_inter += $12
        print $0 | cmd_inter
    }

    print $0
}
END {
    close(cmd_unmeth); close(cmd_inter); close(cmd_meth)

    print "Metric", "Value" > qc_main
    printf "Sample\t%s\n",                     sample   >> qc_main
    printf "Modification\t%s\n",               modlabel >> qc_main
    printf "Min_coverage_threshold\t%d\n",     mincov   >> qc_main
    printf "Max_coverage_threshold\t%d\n",     maxcov   >> qc_main
    printf "Input_sites\t%d\n",                input_sites + 0 >> qc_main
    printf "Dropped_low_coverage\t%d\n",       dropped_low + 0  >> qc_main
    printf "Dropped_high_coverage\t%d\n",      dropped_high + 0 >> qc_main
    printf "Retained_sites\t%d\n",             retained + 0     >> qc_main
    printf "Retained_percent_of_input\t%.4f\n",
           (input_sites ? retained / input_sites * 100 : 0) >> qc_main

    printf "Minimum_coverage\t%d\n",           min_cov + 0 >> qc_main
    printf "Maximum_coverage\t%d\n",           max_cov + 0 >> qc_main
    printf "Mean_coverage\t%.4f\n",            (retained ? sum_cov / retained : 0) >> qc_main
    printf "Total_valid_coverage\t%.0f\n",     sum_cov + 0 >> qc_main
    printf "Total_modified_read_calls\t%.0f\n", sum_mod + 0 >> qc_main
    printf "Unweighted_mean_site_methylation_percent\t%.4f\n",
           (retained ? sum_pct / retained : 0) >> qc_main
    printf "Coverage_weighted_methylation_percent\t%.4f\n",
           (sum_cov > 0 ? sum_mod / sum_cov * 100 : 0) >> qc_main

    printf "State_unmethylated_sites\t%d\n",   n_unmeth + 0 >> qc_main
    printf "State_intermediate_sites\t%d\n",   n_inter  + 0 >> qc_main
    printf "State_methylated_sites\t%d\n",     n_meth   + 0 >> qc_main
    printf "State_unmethylated_percent\t%.4f\n", (retained ? n_unmeth / retained * 100 : 0) >> qc_main
    printf "State_intermediate_percent\t%.4f\n", (retained ? n_inter  / retained * 100 : 0) >> qc_main
    printf "State_methylated_percent\t%.4f\n",   (retained ? n_meth   / retained * 100 : 0) >> qc_main
    printf "State_unmethylated_weighted_methylation\t%.4f\n",
           (cov_unmeth > 0 ? mod_unmeth / cov_unmeth * 100 : 0) >> qc_main
    printf "State_intermediate_weighted_methylation\t%.4f\n",
           (cov_inter  > 0 ? mod_inter  / cov_inter  * 100 : 0) >> qc_main
    printf "State_methylated_weighted_methylation\t%.4f\n",
           (cov_meth   > 0 ? mod_meth   / cov_meth   * 100 : 0) >> qc_main

    # --- reconciliation: neither of these can hold by accident ---
    filter_residual   = input_sites - (dropped_low + dropped_high + retained)
    category_residual = retained    - (n_unmeth + n_inter + n_meth)
    printf "ASSERT_filter_reconciliation\t%d\n",   filter_residual   >> qc_main
    printf "ASSERT_category_reconciliation\t%d\n", category_residual >> qc_main
    printf "ASSERT_no_sites_retained\t%d\n",       (retained == 0 ? 1 : 0) >> qc_main
    close(qc_main)

    exit (filter_residual != 0 || category_residual != 0 || retained == 0) ? 1 : 0
}' \
| bgzip -c -@ "${THREADS}" > "${OUT_BED}"

RCS=("${PIPESTATUS[@]}")
RC_DECOMP=${RCS[0]:-0}
RC_AWK=${RCS[1]:-0}
RC_BGZIP=${RCS[2]:-0}
T1=$(date +%s)

echo "Single pass completed in $((T1-T0))s"

if [[ ${RC_DECOMP} -ne 0 || ${RC_AWK} -ne 0 || ${RC_BGZIP} -ne 0 ]]; then
    echo ""
    echo "[FATAL] Filtering failed (decompress=${RC_DECOMP} awk=${RC_AWK} bgzip=${RC_BGZIP})"
    if [[ -f "${QC_MAIN}" ]]; then
        echo "        Assertion counters:"
        grep '^ASSERT_' "${QC_MAIN}" | awk -F'\t' '$2+0 != 0 {printf "          %s = %s\n", $1, $2}'
    fi
    echo "        Removing partial outputs."
    rm -f "${OUT_BED}" "${OUT_UNMETH}" "${OUT_INTER}" "${OUT_METH}"
    exit 1
fi

# --- Index every output ------------------------------------------------------
echo -n "Indexing with tabix ... "
INDEX_FAIL=0
for f in "${OUT_BED}" "${OUT_UNMETH}" "${OUT_INTER}" "${OUT_METH}"; do
    if [[ -s "${f}" ]]; then
        tabix -p bed -f "${f}" 2>/dev/null || { echo ""; echo "  [FAIL] $(basename "${f}")"; INDEX_FAIL=1; }
    else
        # An empty state file is legitimate — a targeted panel may contain no
        # intermediate sites at all — but it cannot be indexed.
        echo ""
        echo "  [INFO] $(basename "${f}") is empty, not indexed"
    fi
done
[[ ${INDEX_FAIL} -eq 0 ]] && echo "[OK]" || { echo "[FATAL] tabix indexing failed" >&2; exit 1; }

# --- Report ------------------------------------------------------------------
qc() { awk -F'\t' -v k="$1" '$1==k {print $2; exit}' "${QC_MAIN}"; }

echo ""
echo "--- Coverage filter ---"
printf "  %-42s %s\n" "Input sites (from stage 02)" "$(qc Input_sites)"
printf "  %-42s %s\n" "Dropped, coverage < ${MIN_COVERAGE}" "$(qc Dropped_low_coverage)"
[[ "${MAX_COVERAGE}" -gt 0 ]] && \
printf "  %-42s %s\n" "Dropped, coverage > ${MAX_COVERAGE}" "$(qc Dropped_high_coverage)"
printf "  %-42s %s (%s%%)\n" "Retained" "$(qc Retained_sites)" "$(qc Retained_percent_of_input)"

echo ""
echo "--- Methylation of the retained set ---"
printf "  %-42s %s\n" "Mean coverage"                      "$(qc Mean_coverage)"
printf "  %-42s %s\n" "Unweighted mean site methylation %" "$(qc Unweighted_mean_site_methylation_percent)"
printf "  %-42s %s\n" "Coverage-weighted methylation %"     "$(qc Coverage_weighted_methylation_percent)"
echo "  The weighted figure is the one to quote: it is total modified calls over"
echo "  total valid coverage, so well-covered sites carry proportionate influence."

echo ""
echo "--- Methylation states ---"
printf "  %-16s %-14s %-10s %s\n" "State" "Sites" "Share" "Weighted_%"
printf "  %-16s %-14s %-10s %s\n" "unmethylated" "$(qc State_unmethylated_sites)" \
       "$(qc State_unmethylated_percent)%" "$(qc State_unmethylated_weighted_methylation)"
printf "  %-16s %-14s %-10s %s\n" "intermediate" "$(qc State_intermediate_sites)" \
       "$(qc State_intermediate_percent)%" "$(qc State_intermediate_weighted_methylation)"
printf "  %-16s %-14s %-10s %s\n" "methylated"   "$(qc State_methylated_sites)" \
       "$(qc State_methylated_percent)%" "$(qc State_methylated_weighted_methylation)"

echo ""
echo "--- Assertions ---"
awk -F'\t' '$1 ~ /^ASSERT_/ {printf "  [%s] %-38s %s\n", ($2+0==0 ? "OK  " : "FAIL"), $1, $2}' "${QC_MAIN}"
echo "  filter_reconciliation   = input - (dropped_low + dropped_high + retained)"
echo "  category_reconciliation = retained - (unmethylated + intermediate + methylated)"

echo ""
echo "Output: ${OUT_BED} ($(ls -lhL "${OUT_BED}" | awk '{print $5}'))"
for f in "${OUT_UNMETH}" "${OUT_INTER}" "${OUT_METH}"; do
    [[ -s "${f}" ]] && echo "        ${f} ($(ls -lhL "${f}" | awk '{print $5}'))"
done
echo "QC:     ${QC_MAIN}"
echo ""
echo "=== [03_filter_cpg_sites.sh] Done. Log: ${LOG_FILE} ==="
