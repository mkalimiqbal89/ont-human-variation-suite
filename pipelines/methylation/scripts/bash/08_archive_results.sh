#!/usr/bin/env bash
# =============================================================================
# 08_archive_results.sh
# Archives one sample's analysis to the longitudinal store:
#     <archive_root>/<SAMPLE_ID>/<RUN_STAMP>/
# with checksums, provenance, and an entry appended to archive_index.tsv.
#
# EVERY GUARD BELOW EXISTS BECAUSE THE EQUIVALENT SV STAGE GOT IT WRONG FIRST:
#
#   * Sample filtering by EXACT filename. results/ is shared between samples, so
#     copying directories wholesale contaminates the archive with another
#     sample's data. Files are selected with -name "${SAMPLE_ID}.*" and a
#     post-copy assertion re-checks every archived basename.
#
#   * Explicit exit-on-failure after mkdir and cp. Silent continuation produces
#     an archive that looks complete and is not.
#
#   * Preflight permission and free-space checks, before copying anything.
#     Finding out the destination is read-only after 400 MB is wasteful and
#     leaves a partial archive behind.
#
#   * Checksums are verified by reading back what was written, not just
#     generated. An unverified checksum file proves nothing.
#
# Usage:
#   bash scripts/bash/08_archive_results.sh [config.yaml] [--dry-run]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CONFIG_ARG=""
DRY_RUN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        '#'*)      shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -*)        echo "[FATAL] Unknown option: $1" >&2; exit 2 ;;
        *)         CONFIG_ARG="$1"; shift ;;
    esac
done
CONFIG_FILE="${CONFIG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/08_archive_results_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"

echo ""
echo "=== [08_archive_results.sh] Archiving ${SAMPLE_ID} ==="

# --- Archiving disabled? -----------------------------------------------------
if [[ -z "${ARCHIVE_ROOT}" || "${ARCHIVE_ROOT}" == "/path/to/internal/archive/root" ]]; then
    echo "[SKIP] archive.archive_root is unset or still the template placeholder."
    echo "       Archiving disabled. Set it in ${CONFIG_FILE} to enable."
    exit 0
fi

DEST_SAMPLE="${ARCHIVE_ROOT}/${SAMPLE_ID}"
DEST="${DEST_SAMPLE}/${RUN_STAMP}"
INDEX="${ARCHIVE_ROOT}/archive_index.tsv"

echo "Archive root : ${ARCHIVE_ROOT}"
echo "Destination  : ${DEST}"
echo "Intermediates: ${ARCHIVE_INCLUDE_INTERMEDIATES} (all_cov file)"
echo "Compress     : ${ARCHIVE_COMPRESS}"
echo ""

# --- Build the file list, sample-filtered ------------------------------------
# find with an exact basename pattern, NUL-delimited. Never `cp -r results/`:
# that directory is shared across samples.
echo "--- Selecting files for ${SAMPLE_ID} ---"

FILE_LIST="$(mktemp "${TMPDIR:-/tmp}/methyl_archive_list.XXXXXX")"
trap 'rm -f "${FILE_LIST}"' EXIT

if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "[FATAL] Results directory not found: ${OUTPUT_DIR}" >&2
    echo "        Nothing to archive. Run the pipeline first." >&2
    exit 1
fi

while IFS= read -r -d '' f; do
    base="$(basename "${f}")"
    # Skip the large regenerable intermediate unless explicitly requested.
    if [[ "${ARCHIVE_INCLUDE_INTERMEDIATES}" != "true" && "${base}" == *".all_cov."* ]]; then
        continue
    fi
    printf '%s\n' "${f}" >> "${FILE_LIST}"
done < <(find "${OUTPUT_DIR}" -type f -name "${SAMPLE_ID}.*" -print0 2>/dev/null)

N_RESULT_FILES=$(wc -l < "${FILE_LIST}" | tr -d ' ')
if [[ "${N_RESULT_FILES}" -eq 0 ]]; then
    echo "[FATAL] No files named '${SAMPLE_ID}.*' under ${OUTPUT_DIR}" >&2
    echo "        Either the pipeline has not run for this sample, or" >&2
    echo "        sample.sample_id does not match the outputs on disk." >&2
    exit 1
fi

# Logs for this sample are not sample-prefixed (they are per-stage, per-run), so
# they are taken wholesale from LOG_DIR — which is per-repo, not per-sample.
# Recorded as a known caveat in provenance rather than silently mixed in.
N_LOG_FILES=$(find "${LOG_DIR}" -maxdepth 1 -type f 2>/dev/null | wc -l | tr -d ' ')

TOTAL_KB=$(awk '{print}' "${FILE_LIST}" | tr '\n' '\0' | xargs -0 du -k 2>/dev/null \
           | awk '{s+=$1} END {print s+0}')
echo "  ${N_RESULT_FILES} result files, $((TOTAL_KB / 1024)) MB"
echo "  ${N_LOG_FILES} log files"
awk '{print "    " $0}' "${FILE_LIST}" | sed "s|${OUTPUT_DIR}/||" | head -n 12
[[ "${N_RESULT_FILES}" -gt 12 ]] && echo "    ... and $((N_RESULT_FILES - 12)) more"

# --- Preflight ---------------------------------------------------------------
echo ""
echo "--- Preflight ---"

# Create the archive root if needed, and fail loudly if we cannot.
if [[ ! -d "${ARCHIVE_ROOT}" ]]; then
    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "  [DRY]  would create ${ARCHIVE_ROOT}"
    else
        mkdir -p "${ARCHIVE_ROOT}" || {
            echo "  [FAIL] cannot create archive root: ${ARCHIVE_ROOT}" >&2
            echo "         Check the volume is mounted and writable." >&2
            exit 1; }
        echo "  [OK]   created ${ARCHIVE_ROOT}"
    fi
else
    echo "  [OK]   archive root exists"
fi

# Writability, tested by actually writing. `-w` lies on some network mounts.
if [[ ${DRY_RUN} -eq 0 ]]; then
    common_check_writable "${ARCHIVE_ROOT}" || exit 1
fi

# Free space, with a 20% margin.
NEED_KB=$(( TOTAL_KB * 12 / 10 ))
common_check_free_space "${ARCHIVE_ROOT}" "${NEED_KB}" || exit 1

# Refuse to overwrite an existing run directory.
if [[ -e "${DEST}" ]]; then
    echo "  [FAIL] destination already exists: ${DEST}" >&2
    echo "         Run stamps are per-second; wait a moment and retry." >&2
    exit 1
fi
echo "  [OK]   destination is free"

if [[ ${DRY_RUN} -eq 1 ]]; then
    echo ""
    echo "Dry run — nothing copied."
    exit 0
fi

# --- Copy --------------------------------------------------------------------
echo ""
echo "--- Copying ---"
mkdir -p "${DEST}/results" "${DEST}/logs" "${DEST}/config" || {
    echo "[FATAL] could not create destination tree under ${DEST}" >&2; exit 1; }

COPIED=0
while IFS= read -r src; do
    rel="${src#"${OUTPUT_DIR}"/}"
    tgt="${DEST}/results/${rel}"
    mkdir -p "$(dirname "${tgt}")" || {
        echo "[FATAL] mkdir failed for $(dirname "${tgt}")" >&2; exit 1; }
    cp -p "${src}" "${tgt}" || {
        echo "[FATAL] cp failed: ${src} -> ${tgt}" >&2
        echo "        Archive is incomplete; not updating the index." >&2
        exit 1; }
    COPIED=$((COPIED + 1))
done < "${FILE_LIST}"
echo "  [OK]   ${COPIED} result files"

if [[ "${COPIED}" -ne "${N_RESULT_FILES}" ]]; then
    echo "[FATAL] copied ${COPIED} of ${N_RESULT_FILES} files" >&2
    exit 1
fi

# Logs and configs.
if [[ "${N_LOG_FILES}" -gt 0 ]]; then
    find "${LOG_DIR}" -maxdepth 1 -type f -print0 \
        | xargs -0 -I{} cp -p {} "${DEST}/logs/" || {
            echo "[FATAL] copying logs failed" >&2; exit 1; }
    echo "  [OK]   ${N_LOG_FILES} log files"
fi

cp -p "${CONFIG_FILE}" "${DEST}/config/" || { echo "[FATAL] copying config failed" >&2; exit 1; }
cp -p "${REF_CONFIG}"  "${DEST}/config/" || { echo "[FATAL] copying reference config failed" >&2; exit 1; }
echo "  [OK]   config files"

# --- Contamination assertion -------------------------------------------------
# The whole reason this stage is filename-filtered. Verify it actually held.
echo ""
echo "--- Contamination check ---"
FOREIGN=$(find "${DEST}/results" -type f ! -name "${SAMPLE_ID}.*" | wc -l | tr -d ' ')
if [[ "${FOREIGN}" -ne 0 ]]; then
    echo "  [FAIL] ${FOREIGN} archived file(s) are not named ${SAMPLE_ID}.*" >&2
    find "${DEST}/results" -type f ! -name "${SAMPLE_ID}.*" | head -n 5 >&2
    echo "         Removing the archive rather than keeping a contaminated one." >&2
    rm -rf "${DEST}"
    exit 1
fi
echo "  [OK]   all ${COPIED} archived result files belong to ${SAMPLE_ID}"

# --- Checksums, then verify by reading back ---------------------------------
# common_write_checksums (common/lib_common.sh): generates then immediately
# re-reads the checksum file to confirm it matches what was actually archived.
echo ""
echo "--- Checksums ---"
CHECKSUM_FILE="${DEST}/checksums.sha256"
CHECKSUM_STATUS="$(common_write_checksums "${DEST}" "checksums.sha256" results config)"
case "${CHECKSUM_STATUS}" in
    verified:*)
        N_SUMS="${CHECKSUM_STATUS#verified:}"
        echo "  [OK]   ${N_SUMS} checksums written and verified against the archived files"
        ;;
    unavailable)
        echo "  [WARN] no sha256sum or shasum available — checksums skipped"
        ;;
    *)
        echo "  [FAIL] checksum verification failed — archive is not trustworthy" >&2
        rm -rf "${DEST}"
        exit 1
        ;;
esac

# --- Provenance --------------------------------------------------------------
PROV="${DEST}/provenance.tsv"
{
    printf 'Field\tValue\n'
    printf 'sample_id\t%s\n'          "${SAMPLE_ID}"
    printf 'raw_sample_prefix\t%s\n'  "${RAW_SAMPLE_PREFIX}"
    printf 'archived_at\t%s\n'        "$(date '+%Y-%m-%d %H:%M:%S %Z')"
    printf 'run_stamp\t%s\n'          "${RUN_STAMP}"
    printf 'archived_by\t%s@%s\n'     "$(id -un)" "$(uname -n)"
    printf 'os\t%s\n'                 "$(uname -sr)"
    printf 'source_results_dir\t%s\n' "${OUTPUT_DIR}"
    printf 'source_input_dir\t%s\n'   "${INPUT_DIR}"
    printf 'genome_build\t%s\n'       "${GENOME_BUILD}"
    printf 'primary_mod_code\t%s\n'   "${PRIMARY_MOD_CODE}"
    printf 'min_coverage\t%s\n'       "${MIN_COVERAGE}"
    printf 'result_files\t%s\n'       "${COPIED}"
    printf 'log_files\t%s\n'          "${N_LOG_FILES}"
    printf 'size_kb\t%s\n'            "${TOTAL_KB}"
    printf 'checksums\t%s\n'          "${CHECKSUM_STATUS}"
    printf 'intermediates_included\t%s\n' "${ARCHIVE_INCLUDE_INTERMEDIATES}"
    # Honest caveat: LOG_DIR is per-repository, not per-sample, so if several
    # samples were processed against the same checkout the logs directory may
    # contain other samples' stage logs. Results are strictly filtered; logs are not.
    printf 'log_caveat\t%s\n' "logs are per-repo not per-sample; may include other samples' stage logs"
    if [[ -f "${OUTPUT_DIR}/${SAMPLE_ID}.run_manifest.tsv" ]]; then
        printf 'run_manifest\tarchived\n'
    else
        printf 'run_manifest\tabsent (sample not run via 04_run_all.sh)\n'
    fi
} > "${PROV}" || { echo "[FATAL] could not write provenance" >&2; exit 1; }
echo ""
echo "  [OK]   provenance -> ${PROV##*/}"

# --- Optional compression ----------------------------------------------------
FINAL_PATH="${DEST}"
if [[ "${ARCHIVE_COMPRESS}" == "true" ]]; then
    echo ""
    echo "--- Compressing ---"
    TARBALL="${DEST_SAMPLE}/${RUN_STAMP}.tar.gz"
    if tar -czf "${TARBALL}" -C "${DEST_SAMPLE}" "${RUN_STAMP}"; then
        # Only remove the directory once the tarball is proven readable.
        if tar -tzf "${TARBALL}" >/dev/null 2>&1; then
            rm -rf "${DEST}"
            FINAL_PATH="${TARBALL}"
            echo "  [OK]   ${TARBALL##*/} ($(ls -lhL "${TARBALL}" | awk '{print $5}')), directory removed"
        else
            echo "  [FAIL] tarball is not readable; keeping the directory" >&2
            rm -f "${TARBALL}"
        fi
    else
        echo "  [FAIL] tar failed; keeping the uncompressed directory" >&2
        rm -f "${TARBALL}"
    fi
fi

# --- Index -------------------------------------------------------------------
# Appended only after everything above succeeded: the index must never point at
# an archive that failed to complete.
if [[ ! -f "${INDEX}" ]]; then
    printf 'Sample_id\tRaw_prefix\tRun_stamp\tArchived_at\tGenome_build\tMod\tMin_coverage\tResult_files\tSize_kb\tChecksums\tPath\n' \
        > "${INDEX}" || { echo "[FATAL] could not create ${INDEX}" >&2; exit 1; }
fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${SAMPLE_ID}" "${RAW_SAMPLE_PREFIX}" "${RUN_STAMP}" \
    "$(date '+%Y-%m-%d %H:%M:%S')" "${GENOME_BUILD}" "${PRIMARY_MOD_CODE}" \
    "${MIN_COVERAGE}" "${COPIED}" "${TOTAL_KB}" "${CHECKSUM_STATUS}" "${FINAL_PATH}" \
    >> "${INDEX}" || { echo "[FATAL] could not append to ${INDEX}" >&2; exit 1; }

echo ""
echo "--- Archive index ---"
echo "  ${INDEX}"
echo "  $(( $(wc -l < "${INDEX}" | tr -d ' ') - 1 )) archived run(s) recorded:"
awk -F'\t' 'NR>1 {printf "    %-12s %-18s %-10s %s files\n", $1, $3, $6, $8}' "${INDEX}" | tail -n 5

echo ""
echo "Archived to: ${FINAL_PATH}"
echo ""
echo "=== [08_archive_results.sh] Done. Log: ${LOG_FILE} ==="
