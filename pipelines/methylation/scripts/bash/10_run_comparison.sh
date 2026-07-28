#!/usr/bin/env bash
# =============================================================================
# 10_run_comparison.sh
# Validates a sample sheet, then runs the cross-sample comparison.
#
# Validation happens here rather than in R because the failure modes are about
# files and paths, and catching them before R starts gives a clearer message
# than an R error two-thirds of the way through.
#
# Sample sheet: TSV with a header. Required columns:
#     sample_id     as used in the output filenames (e.g. SAMPLE_02)
#     results_dir   the results/ directory holding that sample's outputs
# Optional:
#     label         display name in tables and figures (defaults to sample_id)
#     group         free-text grouping, carried into the global table
#
# Samples may share a results_dir: outputs are sample-prefixed, so one directory
# can hold several samples. Lines beginning with # are ignored.
#
# Usage:
#   bash scripts/bash/10_run_comparison.sh <sample_sheet.tsv> [options]
#
# Options:
#   --name NAME     comparison name, used for the output subdirectory
#   --config FILE   pipeline config (default config/pipeline_config.yaml)
#   --out DIR       output directory (default <output_dir>/11_comparison/<name>)
#   --template      write an example sample sheet and exit
#   --dry-run       validate only
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SHEET=""
NAME=""
CONFIG_ARG=""
OUT_ARG=""
DRY_RUN=0
TEMPLATE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        '#'*)       shift ;;
        --name)     NAME="$2"; shift 2 ;;
        --config)   CONFIG_ARG="$2"; shift 2 ;;
        --out)      OUT_ARG="$2"; shift 2 ;;
        --template) TEMPLATE=1; shift ;;
        --dry-run)  DRY_RUN=1; shift ;;
        -h|--help)  sed -n '2,34p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)         echo "[FATAL] Unknown option: $1" >&2; exit 2 ;;
        *)          SHEET="$1"; shift ;;
    esac
done

if [[ ${TEMPLATE} -eq 1 ]]; then
    TPL="${SHEET:-${REPO_DIR}/config/sample_sheet.example.tsv}"
    cat > "${TPL}" <<'EOF'
# Sample sheet for 10_run_comparison.sh
# Tab-separated. sample_id and results_dir are required.
# Samples may share a results_dir: outputs are sample-prefixed.
sample_id	label	group	results_dir
SAMPLE_01	HLH_S0002	case	/Volumes/Extreme_SSD/Bioinformatics_KAL/ont-methylation-pipeline/results
SAMPLE_02	HLH_S0001	case	/Volumes/Extreme_SSD/Bioinformatics_KAL/ont-methylation-pipeline/results
EOF
    echo "Wrote template: ${TPL}"
    exit 0
fi

if [[ -z "${SHEET}" ]]; then
    echo "[FATAL] No sample sheet given." >&2
    echo "        Write one with: bash scripts/bash/10_run_comparison.sh --template" >&2
    exit 2
fi
if [[ ! -f "${SHEET}" ]]; then
    echo "[FATAL] Sample sheet not found: ${SHEET}" >&2
    exit 1
fi

CONFIG_FILE="${CONFIG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/10_run_comparison_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo ""
echo "=== [10_run_comparison.sh] Cross-sample comparison ==="
echo "Sample sheet: ${SHEET}"

# --- Parse and validate the sheet -------------------------------------------
HEADER="$(grep -v '^#' "${SHEET}" | head -n1)"
col_index() {
    echo "${HEADER}" | awk -F'\t' -v want="$1" \
        '{for (i=1;i<=NF;i++) if ($i==want) {print i; exit}}'
}
C_ID="$(col_index sample_id)"
C_DIR="$(col_index results_dir)"
C_LAB="$(col_index label)"

ERRORS=0
if [[ -z "${C_ID}" || -z "${C_DIR}" ]]; then
    echo "[FATAL] Sample sheet header must contain 'sample_id' and 'results_dir'." >&2
    echo "        Found: ${HEADER}" >&2
    echo "        Columns must be TAB separated, not spaces." >&2
    exit 1
fi

echo ""
echo "--- Samples ---"
N_SAMPLES=0
SEEN_IDS=""
# Read the data rows, skipping comments and the header.
while IFS=$'\t' read -r -a row; do
    sid="${row[$((C_ID - 1))]:-}"
    rdir="${row[$((C_DIR - 1))]:-}"
    lab="${sid}"
    [[ -n "${C_LAB}" ]] && lab="${row[$((C_LAB - 1))]:-${sid}}"
    [[ -z "${sid}" ]] && continue

    N_SAMPLES=$((N_SAMPLES + 1))
    printf "  %-14s %-14s %s\n" "${sid}" "${lab}" "${rdir}"

    # Duplicate sample_id would silently compare a sample against itself.
    case " ${SEEN_IDS} " in
        *" ${sid} "*)
            echo "    [FAIL] duplicate sample_id '${sid}'"
            ERRORS=$((ERRORS + 1)) ;;
        *) SEEN_IDS="${SEEN_IDS} ${sid}" ;;
    esac

    if [[ ! -d "${rdir}" ]]; then
        echo "    [FAIL] results_dir does not exist: ${rdir}"
        ERRORS=$((ERRORS + 1))
        continue
    fi

    # Stage 05 and 06 outputs are what the comparison consumes.
    for rel in "06_gene_summary/${sid}.gene_methylation_summary.tsv" \
               "07_promoter_summary/${sid}.promoter_methylation_summary.tsv" \
               "08_cpg_island_summary/${sid}.cpg_island_methylation_summary.tsv" \
               "03_global_methylation/${sid}.global_methylation_summary.tsv"; do
        if [[ ! -f "${rdir}/${rel}" ]]; then
            echo "    [FAIL] missing ${rel}"
            echo "           run stages 05 and 06 for ${sid} first"
            ERRORS=$((ERRORS + 1))
        fi
    done
done < <(grep -v '^#' "${SHEET}" | tail -n +2)

echo ""
if [[ "${N_SAMPLES}" -lt 2 ]]; then
    echo "[FATAL] Need at least 2 samples; found ${N_SAMPLES}." >&2
    echo "        Check the sheet is TAB separated and has a header row." >&2
    exit 1
fi
if [[ ${ERRORS} -gt 0 ]]; then
    echo "[FATAL] ${ERRORS} problem(s) with the sample sheet. See above." >&2
    exit 1
fi
echo "  [OK]   ${N_SAMPLES} samples, all required inputs present"

# --- Output directory --------------------------------------------------------
if [[ -z "${NAME}" ]]; then
    NAME="$(echo "${SEEN_IDS}" | tr -s ' ' '_' | sed 's/^_//;s/_$//')"
    [[ ${#NAME} -gt 60 ]] && NAME="comparison_${N_SAMPLES}_samples"
fi
OUT_DIR="${OUT_ARG:-${OUTPUT_DIR}/11_comparison/${NAME}}"
echo ""
echo "Comparison  : ${NAME}"
echo "Output      : ${OUT_DIR}"

if [[ ${DRY_RUN} -eq 1 ]]; then
    echo ""
    echo "Dry run — sheet validated, nothing computed."
    exit 0
fi

mkdir -p "${OUT_DIR}" || { echo "[FATAL] could not create ${OUT_DIR}" >&2; exit 1; }
cp -p "${SHEET}" "${OUT_DIR}/sample_sheet.tsv" || {
    echo "[FATAL] could not copy the sample sheet into the output" >&2; exit 1; }

# --- Run ---------------------------------------------------------------------
echo ""
T0=$(date +%s)
Rscript "${REPO_DIR}/scripts/R/09_compare_samples.R" \
        "${SHEET}" "${OUT_DIR}" "${CONFIG_FILE}"
RC=$?
T1=$(date +%s)

echo ""
if [[ ${RC} -ne 0 ]]; then
    echo "=== [10_run_comparison.sh] FAILED (exit ${RC}) after $((T1-T0))s. Log: ${LOG_FILE} ==="
    exit "${RC}"
fi
echo "Completed in $((T1-T0))s"
echo "=== [10_run_comparison.sh] Done. Log: ${LOG_FILE} ==="
