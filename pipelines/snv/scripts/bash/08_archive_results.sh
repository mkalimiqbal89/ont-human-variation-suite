#!/usr/bin/env bash
# =============================================================================
# 08_archive_results.sh (pipelines/snv)
# Archives SNV analysis outputs to institutional longitudinal storage, computes
# SHA-256 checksums, writes run provenance, and updates archive_index_snv.tsv.
#
# Usage:
#   bash scripts/bash/08_archive_results.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}"

if [[ -z "${ARCHIVE_ROOT:-}" || "${ARCHIVE_ROOT}" == "/path/to/internal/archive/root" ]]; then
    echo "=== [08_archive_results.sh] Archiving is disabled (archive_root unset or template default) ==="
    exit 0
fi

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
ARCHIVE_DIR="${ARCHIVE_ROOT}/${SAMPLE_ID}/snv/${TIMESTAMP}"

echo "=== [08_archive_results.sh] Archiving SNV results for ${SAMPLE_ID} ==="
echo "Destination: ${ARCHIVE_DIR}"

mkdir -p "${ARCHIVE_DIR}/results" "${ARCHIVE_DIR}/logs"

# Copy results & logs
if [[ -d "${OUTPUT_DIR}" ]]; then
    cp -r "${OUTPUT_DIR}/"* "${ARCHIVE_DIR}/results/" 2>/dev/null || true
fi
if [[ -d "${LOG_DIR}" ]]; then
    cp -r "${LOG_DIR}/"* "${ARCHIVE_DIR}/logs/" 2>/dev/null || true
fi

# Compute SHA-256 checksums
CHECKSUM_FILE="${ARCHIVE_DIR}/checksums.sha256"
echo "Computing SHA-256 checksums..."
(
    cd "${ARCHIVE_DIR}"
    if command -v sha256sum >/dev/null 2>&1; then
        find results logs -type f -exec sha256sum {} + > "${CHECKSUM_FILE}"
    else
        find results logs -type f -exec shasum -a 256 {} + > "${CHECKSUM_FILE}"
    fi
)

# Write provenance record
PROV_FILE="${ARCHIVE_DIR}/provenance.txt"
{
    echo "sample_id: ${SAMPLE_ID}"
    echo "pipeline: snv"
    echo "timestamp: ${TIMESTAMP}"
    echo "host: $(hostname)"
    echo "os: $(uname -s) $(uname -m)"
    echo "git_commit: $(git -C "${REPO_DIR}" rev-parse HEAD 2>/dev/null || echo 'unknown')"
    echo "checksum_count: $(wc -l < "${CHECKSUM_FILE}" | tr -d ' ')"
} > "${PROV_FILE}"

# Update archive index file
INDEX_FILE="${ARCHIVE_ROOT}/archive_index_snv.tsv"
if [[ ! -f "${INDEX_FILE}" ]]; then
    echo -e "sample_id\tpipeline\ttimestamp\tarchive_path" > "${INDEX_FILE}"
fi
echo -e "${SAMPLE_ID}\tsnv\t${TIMESTAMP}\t${ARCHIVE_DIR}" >> "${INDEX_FILE}"

echo "=== [08_archive_results.sh] Archiving COMPLETE ==="
