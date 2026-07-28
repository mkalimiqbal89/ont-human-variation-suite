#!/usr/bin/env bash
# =============================================================================
# 04_run_all.sh
# Runs the per-sample pipeline end to end and records a provenance manifest.
#
# Stops at the first failing stage. A half-finished methylation run that keeps
# going produces outputs that look complete, and the only sign of trouble is a
# number being quietly wrong several stages later.
#
# The manifest exists for the methods section: it pins the config file and its
# checksum, the reference files and their checksums, every tool version, the git
# commit if the repo is under version control, and the wall time of each stage.
# Reproducing a figure a year from now needs all of that.
#
# Usage:
#   bash scripts/bash/04_run_all.sh [config.yaml] [options]
#
# Options:
#   --from NN      start at stage NN (01, 02, 03, 05, 06, 07)
#   --to NN        stop after stage NN
#   --only NN      run just stage NN
#   --list         show the stage list and exit
#   --dry-run      print what would run, without running it
#   --skip-validate  skip stage 01 (equivalent to --from 02)
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
R_DIR="${REPO_DIR}/scripts/R"

CONFIG_ARG=""
FROM=""
TO=""
DRY_RUN=0
LIST_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        '#'*)            shift ;;   # zsh passes inline comments through
        --from)          FROM="$2"; shift 2 ;;
        --to)            TO="$2";   shift 2 ;;
        --only)          FROM="$2"; TO="$2"; shift 2 ;;
        --skip-validate) FROM="02"; shift ;;
        --list)          LIST_ONLY=1; shift ;;
        --dry-run)       DRY_RUN=1; shift ;;
        -h|--help)       sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)              echo "[FATAL] Unknown option: $1" >&2; exit 2 ;;
        *)               CONFIG_ARG="$1"; shift ;;
    esac
done

CONFIG_FILE="${CONFIG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# --- Stage table -------------------------------------------------------------
# id : script path relative to the repo : description
STAGE_IDS=(01 02 03 05 06 07)
declare -A STAGE_SCRIPT=(
    [01]="scripts/bash/01_validate_inputs.sh"
    [02]="scripts/bash/02_bedmethyl_to_tsv.sh"
    [03]="scripts/bash/03_filter_cpg_sites.sh"
    [05]="scripts/bash/05_annotate_regions.sh"
    [06]="scripts/bash/06_summary_stats.sh"
    [07]="scripts/R/07_generate_report.R"
)
declare -A STAGE_DESC=(
    [01]="Validate inputs and reference bundle"
    [02]="Extract primary mod code, single full pass"
    [03]="Coverage filter and methylation-state splits"
    [05]="Gene / promoter / CpG-island aggregation"
    [06]="Global, distribution, chromosome, feature-class tables"
    [07]="Publication figures"
)

if [[ ${LIST_ONLY} -eq 1 ]]; then
    echo "Stages:"
    for id in "${STAGE_IDS[@]}"; do
        avail="present"; [[ -f "${REPO_DIR}/${STAGE_SCRIPT[$id]}" ]] || avail="NOT YET WRITTEN"
        printf "  %s  %-52s %-12s %s\n" "${id}" "${STAGE_DESC[$id]}" "[${avail}]" "${STAGE_SCRIPT[$id]}"
    done
    exit 0
fi

# --- Environment -------------------------------------------------------------
# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/04_run_all_${RUN_STAMP}.log"
exec > >(tee -a "${LOG_FILE}") 2>&1

MANIFEST="${OUTPUT_DIR}/${SAMPLE_ID}.run_manifest.tsv"

echo ""
echo "============================================================"
echo " ONT Methylation Pipeline — full run"
echo "============================================================"
echo "Sample  : ${SAMPLE_ID} (${RAW_SAMPLE_PREFIX})"
echo "Config  : ${CONFIG_FILE}"
echo "Started : $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# --- Portable checksum -------------------------------------------------------
# macOS ships shasum, not sha256sum; the HPC usually has the reverse.
hash_file() {
    [[ -f "$1" ]] || { echo "absent"; return; }
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "no_sha_tool"
    fi
}

# --- Which stages to run -----------------------------------------------------
in_range() {
    local id="$1"
    [[ -n "${FROM}" && "${id}" < "${FROM}" ]] && return 1
    [[ -n "${TO}"   && "${id}" > "${TO}"   ]] && return 1
    return 0
}

TO_RUN=()
for id in "${STAGE_IDS[@]}"; do
    in_range "${id}" || continue
    if [[ ! -f "${REPO_DIR}/${STAGE_SCRIPT[$id]}" ]]; then
        echo "[SKIP] stage ${id} (${STAGE_DESC[$id]}) — not yet written"
        continue
    fi
    TO_RUN+=("${id}")
done

if [[ ${#TO_RUN[@]} -eq 0 ]]; then
    echo "[FATAL] No stages selected. Try --list." >&2
    exit 1
fi

echo "Stages to run: ${TO_RUN[*]}"
echo ""

if [[ ${DRY_RUN} -eq 1 ]]; then
    for id in "${TO_RUN[@]}"; do
        printf "  would run  %s  %s\n" "${id}" "${STAGE_SCRIPT[$id]}"
    done
    echo ""
    echo "Dry run — nothing executed."
    exit 0
fi

# --- Manifest header ---------------------------------------------------------
mkdir -p "${OUTPUT_DIR}"
{
    printf 'Field\tValue\n'
    printf 'pipeline\tont-methylation-pipeline\n'
    printf 'run_started\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf 'run_stamp\t%s\n' "${RUN_STAMP}"
    printf 'sample_id\t%s\n' "${SAMPLE_ID}"
    printf 'raw_sample_prefix\t%s\n' "${RAW_SAMPLE_PREFIX}"
    printf 'host\t%s\n' "$(uname -n)"
    printf 'os\t%s\n' "$(uname -sr)"
    printf 'config_file\t%s\n' "${CONFIG_FILE}"
    printf 'config_sha256\t%s\n' "$(hash_file "${CONFIG_FILE}")"
    printf 'reference_config\t%s\n' "${REF_CONFIG}"
    printf 'reference_config_sha256\t%s\n' "$(hash_file "${REF_CONFIG}")"
    printf 'genome_build\t%s\n' "${GENOME_BUILD}"
    printf 'primary_mod_code\t%s\n' "${PRIMARY_MOD_CODE}"
    printf 'min_coverage\t%s\n' "${MIN_COVERAGE}"
    printf 'contig_selection\t%s\n' \
        "$([[ "${PRIMARY_CONTIGS_ONLY}" == "true" ]] && echo "allowlist:${INCLUDE_CONTIGS_REGEX}" || echo "blocklist:${EXCLUDE_CONTIGS_REGEX}")"
    printf 'gene_bed_sha256\t%s\n' "$(hash_file "${GENE_BED}")"
    printf 'promoter_bed_sha256\t%s\n' "$(hash_file "${PROMOTER_BED}")"
    printf 'cpg_island_bed_sha256\t%s\n' "$(hash_file "${CPG_ISLAND_BED}")"

    # Git provenance, carefully. Three things can go wrong and all of them
    # produce a manifest that looks fine:
    #   1. `rev-parse HEAD` on a repo with no commits prints the literal string
    #      "HEAD" to stdout and still exits non-zero — hence --verify.
    #   2. git searches parent directories, so it can report a commit from an
    #      enclosing repository that has nothing to do with this pipeline.
    #      Compare --show-toplevel against REPO_DIR.
    #   3. A recorded commit is misleading if the working tree is dirty.
    GIT_TOP=""
    if command -v git >/dev/null 2>&1; then
        GIT_TOP="$(git -C "${REPO_DIR}" rev-parse --show-toplevel 2>/dev/null || true)"
    fi
    if [[ -z "${GIT_TOP}" ]]; then
        printf 'git_commit\tnot_a_git_repo\n'
    elif [[ "${GIT_TOP}" != "${REPO_DIR}" ]]; then
        printf 'git_commit\tenclosing_repo_only:%s\n' "${GIT_TOP}"
    elif GIT_SHA="$(git -C "${REPO_DIR}" rev-parse --verify HEAD 2>/dev/null)"; then
        printf 'git_commit\t%s\n' "${GIT_SHA}"
        printf 'git_dirty\t%s\n' \
            "$([[ -n "$(git -C "${REPO_DIR}" status --porcelain 2>/dev/null)" ]] && echo yes || echo no)"
        printf 'git_branch\t%s\n' \
            "$(git -C "${REPO_DIR}" rev-parse --abbrev-ref HEAD 2>/dev/null)"
    else
        printf 'git_commit\trepo_initialised_but_no_commits\n'
    fi

    if [[ -f "${LOG_DIR}/tool_versions.txt" ]]; then
        while IFS= read -r line; do
            printf 'tool\t%s\n' "${line}"
        done < "${LOG_DIR}/tool_versions.txt"
    fi
} > "${MANIFEST}"

# --- Run ---------------------------------------------------------------------
declare -a RESULT_ID=() RESULT_STATUS=() RESULT_SECS=()
OVERALL_RC=0
T_ALL0=$(date +%s)

for id in "${TO_RUN[@]}"; do
    script="${REPO_DIR}/${STAGE_SCRIPT[$id]}"
    echo "------------------------------------------------------------"
    echo ">>> Stage ${id}: ${STAGE_DESC[$id]}"
    echo "------------------------------------------------------------"
    t0=$(date +%s)

    if [[ "${script}" == *.R ]]; then
        Rscript "${script}" "${CONFIG_FILE}"
        rc=$?
    else
        bash "${script}" "${CONFIG_FILE}"
        rc=$?
    fi

    t1=$(date +%s)
    secs=$((t1 - t0))
    RESULT_ID+=("${id}")
    RESULT_SECS+=("${secs}")

    if [[ ${rc} -ne 0 ]]; then
        RESULT_STATUS+=("FAILED(${rc})")
        echo ""
        echo "[FATAL] Stage ${id} failed with exit ${rc} after ${secs}s."
        echo "        Stopping here rather than running later stages on incomplete"
        echo "        output. Fix the cause and resume with:"
        echo "            bash scripts/bash/04_run_all.sh ${CONFIG_ARG} --from ${id}"
        OVERALL_RC=1
        break
    fi
    RESULT_STATUS+=("ok")
    echo ""
    echo "<<< Stage ${id} completed in ${secs}s"
    echo ""
done

T_ALL1=$(date +%s)

# --- Manifest footer + summary ----------------------------------------------
{
    for i in "${!RESULT_ID[@]}"; do
        printf 'stage\t%s\t%s\t%ss\n' "${RESULT_ID[$i]}" "${RESULT_STATUS[$i]}" "${RESULT_SECS[$i]}"
    done
    printf 'run_finished\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf 'total_seconds\t%s\n' "$((T_ALL1 - T_ALL0))"
    printf 'overall_status\t%s\n' "$([[ ${OVERALL_RC} -eq 0 ]] && echo success || echo failed)"
} >> "${MANIFEST}"

echo "============================================================"
echo " Summary"
echo "============================================================"
printf "  %-6s %-52s %-14s %s\n" "Stage" "Description" "Status" "Time"
for i in "${!RESULT_ID[@]}"; do
    id="${RESULT_ID[$i]}"
    printf "  %-6s %-52s %-14s %ss\n" "${id}" "${STAGE_DESC[$id]}" "${RESULT_STATUS[$i]}" "${RESULT_SECS[$i]}"
done
echo ""
printf "  Total: %ss\n" "$((T_ALL1 - T_ALL0))"
echo ""
echo "Manifest: ${MANIFEST}"
echo "Log     : ${LOG_FILE}"
echo ""

if [[ ${OVERALL_RC} -eq 0 ]]; then
    echo "=== [04_run_all.sh] Run completed successfully. ==="
else
    echo "=== [04_run_all.sh] Run FAILED. ==="
fi
exit ${OVERALL_RC}
