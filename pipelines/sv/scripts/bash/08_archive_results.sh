#!/usr/bin/env bash
# =============================================================================
# 08_archive_results.sh
#
# Copies this run's full results (real sample identifiers, all categories,
# QC summary, figures) plus provenance metadata (the exact config used, tool
# versions, pipeline git commit) into a timestamped, permanent record on
# internal institutional storage. This is entirely separate from the public
# GitHub repo -- nothing this script touches is ever git-tracked, and the
# destination (config: archive.archive_root) is expected to be institutional
# storage (e.g. /media/HCGR/NGS_RAW_DATA/), never a path inside the repo.
#
# Each run appends a new timestamped snapshot rather than overwriting a
# previous archive for the same sample, so re-running with different
# thresholds later preserves a full history rather than clobbering it. A
# running master index (archive_index.tsv) is appended to on every call,
# giving a single queryable log across every sample ever archived.
#
# Usage:
#   bash scripts/bash/08_archive_results.sh [path/to/pipeline_config.yaml]
#
# Requires config: archive.archive_root to be set. If unset, exits cleanly
# with a message rather than failing -- archiving is opt-in, not assumed.
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/08_archive_results_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "=== [08_archive_results.sh] Archiving results for sample: ${SAMPLE_ID} ==="

# --- Read archive config ------------------------------------------------------
ARCHIVE_ROOT="$(grep -E '^[[:space:]]*archive_root:' "${CONFIG_FILE}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//')"
COMPRESS="$(grep -E '^[[:space:]]*compress:' "${CONFIG_FILE}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//')"

if [[ -z "${ARCHIVE_ROOT}" || "${ARCHIVE_ROOT}" == *"/path/to/"* ]]; then
    echo "[INFO] archive.archive_root not set (or still a placeholder) in ${CONFIG_FILE}."
    echo "       Archiving is opt-in and disabled. Set archive.archive_root to enable it."
    exit 0
fi

if [[ ! -d "${ARCHIVE_ROOT}" ]]; then
    echo "[FATAL] archive_root does not exist or is not accessible: ${ARCHIVE_ROOT}"
    echo "        Confirm the internal storage mount is available from this host."
    exit 1
fi

# --- Preflight: confirm this user can actually write here before doing --
# anything. Catches ownership/permission mismatches immediately with a
# clear message, instead of failing halfway through and printing a false
# "Done" at the end.
if ! common_check_writable "${ARCHIVE_ROOT}"; then
    echo "[FATAL] No write permission in archive_root: ${ARCHIVE_ROOT}"
    echo "        Running as: $(whoami)"
    echo "        Check ownership/permissions (e.g. 'ls -ld ${ARCHIVE_ROOT}')"
    echo "        and fix as root: chown -R \$(whoami):\$(whoami) ${ARCHIVE_ROOT}"
    echo "        or chmod -R g+ws ${ARCHIVE_ROOT} if multiple users archive here."
    exit 1
fi

# --- Build destination path ---------------------------------------------------
# Grouped by RAW identifier (matches how the raw Epi2ME output for this
# sample is already organized on internal storage), timestamped so repeat
# runs (e.g. after threshold changes) accumulate rather than overwrite.
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
DEST_DIR="${ARCHIVE_ROOT}/${RAW_SAMPLE_PREFIX}/pipeline_run_${TIMESTAMP}"

if [[ -d "${OUTPUT_DIR}" ]] && [[ -z "$(find "${OUTPUT_DIR}" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
    echo "[FATAL] ${OUTPUT_DIR} exists but is empty. Run the pipeline (04_run_all.sh) before archiving."
    exit 1
fi
if [[ ! -d "${OUTPUT_DIR}" ]]; then
    echo "[FATAL] Results directory not found: ${OUTPUT_DIR}. Run the pipeline first."
    exit 1
fi

echo "Source : ${OUTPUT_DIR}"
echo "Dest   : ${DEST_DIR}"
mkdir -p "${DEST_DIR}" || { echo "[FATAL] Could not create ${DEST_DIR}"; exit 1; }

# --- Copy results --------------------------------------------------------------
echo ""
echo "Copying results..."
# Copy ONLY this sample's files, preserving the relative directory structure
# (results/deletions/, results/qc_summary/figures/, etc.) but never pulling
# in another sample's data -- results/ is a shared directory across every
# sample this pipeline has ever processed, and a naive whole-directory copy
# would silently mix other samples' data into this sample's archive folder.
mkdir -p "${DEST_DIR}/results" || { echo "[FATAL] Could not create ${DEST_DIR}/results"; exit 1; }
N_COPIED=0
while IFS= read -r -d '' f; do
    rel="${f#${OUTPUT_DIR}/}"
    dest_path="${DEST_DIR}/results/${rel}"
    mkdir -p "$(dirname "${dest_path}")" || { echo "[FATAL] Could not create $(dirname "${dest_path}")"; exit 1; }
    cp "${f}" "${dest_path}" || { echo "[FATAL] Could not copy ${f}"; exit 1; }
    N_COPIED=$((N_COPIED+1))
done < <(find "${OUTPUT_DIR}" -type f -name "${SAMPLE_ID}.*" -print0)

if [[ ${N_COPIED} -eq 0 ]]; then
    echo "[FATAL] No files matching '${SAMPLE_ID}.*' found under ${OUTPUT_DIR}"
    echo "        Nothing to archive -- did the pipeline actually run for this sample_id?"
    rm -rf "${DEST_DIR}"
    exit 1
fi
echo "  ${N_COPIED} files copied (filtered to ${SAMPLE_ID}.* only)."

# --- Copy provenance: exact config used, tool versions, git commit ----------
echo "Copying provenance metadata..."
mkdir -p "${DEST_DIR}/provenance" || { echo "[FATAL] Could not create provenance dir"; exit 1; }
cp "${CONFIG_FILE}" "${DEST_DIR}/provenance/pipeline_config_used.yaml" || { echo "[FATAL] Could not copy config"; exit 1; }

REF_CONFIG_RELATIVE="$(grep -E '^[[:space:]]*config_file:' "${CONFIG_FILE}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//')"
if [[ -n "${REF_CONFIG_RELATIVE}" ]]; then
    if [[ "${REF_CONFIG_RELATIVE}" = /* ]]; then REF_CONFIG="${REF_CONFIG_RELATIVE}"; else REF_CONFIG="${REPO_DIR}/${REF_CONFIG_RELATIVE}"; fi
    [[ -f "${REF_CONFIG}" ]] && cp "${REF_CONFIG}" "${DEST_DIR}/provenance/reference_paths_used.yaml"
fi

[[ -f "${LOG_DIR}/tool_versions.txt" ]] && cp "${LOG_DIR}/tool_versions.txt" "${DEST_DIR}/provenance/tool_versions.txt"

GIT_COMMIT="unknown"
if git -C "${REPO_DIR}" rev-parse --short HEAD > /dev/null 2>&1; then
    GIT_COMMIT="$(git -C "${REPO_DIR}" rev-parse --short HEAD)"
    GIT_DIRTY=""
    git -C "${REPO_DIR}" diff --quiet 2>/dev/null || GIT_DIRTY=" (uncommitted changes present at archive time)"
    echo "${GIT_COMMIT}${GIT_DIRTY}" > "${DEST_DIR}/provenance/pipeline_git_commit.txt"
fi

cat > "${DEST_DIR}/provenance/archive_manifest.txt" << EOF
sample_id: ${SAMPLE_ID}
raw_sample_prefix: ${RAW_SAMPLE_PREFIX}
archived_at: ${TIMESTAMP}
archived_from_host: $(hostname)
pipeline_git_commit: ${GIT_COMMIT}
source_output_dir: ${OUTPUT_DIR}
EOF
[[ -f "${DEST_DIR}/provenance/archive_manifest.txt" ]] || { echo "[FATAL] Could not write archive_manifest.txt"; exit 1; }

# --- Checksums for integrity verification ------------------------------------
# SHA-256, generated then immediately verified by reading the archived files
# back (common_write_checksums, common/lib_common.sh) -- a checksum file that
# was only ever generated proves nothing about what actually landed on disk.
# This replaces a previous md5sum-only pass that never verified itself.
echo "Computing checksums..."
CHECKSUM_STATUS="$(common_write_checksums "${DEST_DIR}" "checksums.sha256" results provenance)"
case "${CHECKSUM_STATUS}" in
    verified:*)
        N_FILES="${CHECKSUM_STATUS#verified:}"
        echo "  ${N_FILES} files checksummed and verified."
        ;;
    unavailable)
        echo "  [WARN] no sha256sum or shasum available -- checksums skipped."
        ;;
    *)
        echo "[FATAL] Checksum verification failed -- archive is not trustworthy"
        rm -rf "${DEST_DIR}"
        exit 1
        ;;
esac

# --- Optional compression -----------------------------------------------------
if [[ "${COMPRESS}" == "true" ]]; then
    echo "Compressing archive (archive.compress=true)..."
    ( cd "$(dirname "${DEST_DIR}")" && tar -czf "$(basename "${DEST_DIR}").tar.gz" "$(basename "${DEST_DIR}")" && rm -rf "$(basename "${DEST_DIR}")" )
    echo "  Archived (compressed): ${DEST_DIR}.tar.gz"
    FINAL_LOCATION="${DEST_DIR}.tar.gz"
else
    FINAL_LOCATION="${DEST_DIR}"
fi

# --- Append to master index across all samples ever archived ----------------
INDEX_FILE="${ARCHIVE_ROOT}/archive_index.tsv"
if [[ ! -f "${INDEX_FILE}" ]]; then
    echo -e "timestamp\tsample_id\traw_sample_prefix\tgit_commit\tn_deletions\tn_insertions\tn_duplications\tn_inversions\tn_translocations\tn_gene_fusions\tn_gene_disruptions\tarchive_path" > "${INDEX_FILE}" \
        || { echo "[FATAL] Could not create index file: ${INDEX_FILE}"; exit 1; }
fi

get_count() {
    local category="$1"
    local tsv="${OUTPUT_DIR}/${category}/${SAMPLE_ID}.${category}.tsv"
    [[ -f "${tsv}" ]] && echo $(( $(wc -l < "${tsv}") - 1 )) || echo "NA"
}
N_DEL=$(get_count "deletions")
N_INS=$(get_count "insertions")
N_DUP=$(get_count "duplications")
N_INV=$(get_count "inversions")
N_TRA=$(get_count "translocations")
GF_FILE="${OUTPUT_DIR}/gene_fusions/${SAMPLE_ID}.gene_fusions.tsv"
GD_FILE="${OUTPUT_DIR}/gene_fusions/${SAMPLE_ID}.gene_disruptions.tsv"
N_GF=$([[ -f "${GF_FILE}" ]] && echo $(( $(wc -l < "${GF_FILE}") - 1 )) || echo "NA")
N_GD=$([[ -f "${GD_FILE}" ]] && echo $(( $(wc -l < "${GD_FILE}") - 1 )) || echo "NA")

echo -e "${TIMESTAMP}\t${SAMPLE_ID}\t${RAW_SAMPLE_PREFIX}\t${GIT_COMMIT}\t${N_DEL}\t${N_INS}\t${N_DUP}\t${N_INV}\t${N_TRA}\t${N_GF}\t${N_GD}\t${FINAL_LOCATION}" >> "${INDEX_FILE}" \
    || { echo "[FATAL] Could not append to index file: ${INDEX_FILE}"; exit 1; }

echo ""
echo "=== [08_archive_results.sh] Done ==="
echo "  Archived to : ${FINAL_LOCATION}"
echo "  Checksums   : ${DEST_DIR}/checksums.sha256 (or inside the .tar.gz if compressed)"
echo "  Master index: ${INDEX_FILE} (now $(($(wc -l < "${INDEX_FILE}") - 1)) sample runs recorded)"
