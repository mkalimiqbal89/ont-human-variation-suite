#!/usr/bin/env bash
# =============================================================================
# 00_setup_env.sh
# Verifies required tools are on PATH, records their versions for the
# manuscript methods section, creates working directories, and exports
# config values as shell variables for downstream scripts to source.
#
# Usage:
#   source scripts/bash/00_setup_env.sh [path/to/pipeline_config.yaml]
#
# NOTE: this script must be SOURCED (not executed) by other scripts in this
# pipeline, since it exports variables into the calling shell. It can also be
# run standalone just to check the environment.
# =============================================================================

set -uo pipefail  # no -e here: we want to report ALL missing tools, not exit on first

# --- Resolve repo root regardless of where this is called from -------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "[FATAL] Config file not found: ${CONFIG_FILE}" >&2
    return 1 2>/dev/null || exit 1
fi

echo "=== [00_setup_env.sh] Using config: ${CONFIG_FILE} ==="

# --- Minimal YAML value getter (avoids requiring yq on the HPC) ------------
# Reads a simple "key: value" line under a given top-level/nested key path.
# For this pipeline's flat-ish config, a targeted grep/sed is sufficient and
# avoids adding a yq dependency requirement to the HPC environment.
yaml_get() {
    local key="$1"
    grep -E "^[[:space:]]*${key}:" "${CONFIG_FILE}" | head -n1 \
        | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//'
}

# --- Export config values as shell variables --------------------------------
export SAMPLE_ID="$(yaml_get 'sample_id')"
export RAW_SAMPLE_PREFIX="$(yaml_get 'raw_sample_prefix')"
export INPUT_DIR="$(yaml_get 'input_dir')"
export OUTPUT_DIR="$(yaml_get 'output_dir')"
export WORK_DIR="$(yaml_get 'work_dir')"
export REPO_DIR_CFG="$(yaml_get 'repo_dir')"
export MIN_SV_LENGTH="$(yaml_get 'min_sv_length')"
export MIN_READ_SUPPORT="$(yaml_get 'min_read_support')"
export MIN_QUAL="$(yaml_get 'min_qual')"
export MIN_VAF="$(yaml_get 'min_vaf')"
export THREADS="$(yaml_get 'threads')"
LOG_DIR_RAW="$(yaml_get 'log_dir')"
if [[ "${LOG_DIR_RAW}" = /* ]]; then
    export LOG_DIR="${LOG_DIR_RAW}"
else
    export LOG_DIR="${REPO_DIR}/${LOG_DIR_RAW}"
fi

# Substitute ${sample.raw_sample_prefix} placeholders used in input_files
# section. Input filenames use the RAW prefix (actual Epi2ME output naming),
# never the anonymized SAMPLE_ID, since that's what's really on disk.
resolve_filename() {
    yaml_get "$1" | sed "s|\${sample.raw_sample_prefix}|${RAW_SAMPLE_PREFIX}|g"
}
export SV_VCF_NAME="$(resolve_filename 'sv_vcf')"
export CNV_VCF_NAME="$(resolve_filename 'cnv_vcf')"
export SNV_VCF_NAME="$(resolve_filename 'snv_vcf')"
export STR_VCF_NAME="$(resolve_filename 'str_vcf')"

echo "Sample ID     : ${SAMPLE_ID}"
echo "Input dir     : ${INPUT_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Work dir      : ${WORK_DIR}"

# --- Check required tools ---------------------------------------------------
REQUIRED_TOOLS=(bcftools bedtools tabix Rscript)
MISSING_TOOLS=()

echo ""
echo "=== Tool check ==="
mkdir -p "${LOG_DIR}"
TOOL_VERSION_LOG="${LOG_DIR}/tool_versions.txt"
: > "${TOOL_VERSION_LOG}"

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "${tool}" >/dev/null 2>&1; then
        VER=$("${tool}" --version 2>&1 | head -n1)
        if [[ $? -ne 0 || "${VER}" == *"error"* || "${VER}" == *"cannot open shared object"* ]]; then
            echo "  [FAIL] ${tool} found on PATH but failed to run: ${VER}"
            MISSING_TOOLS+=("${tool}")
        else
            echo "  [OK]   ${tool} -> ${VER}"
            echo "${tool}: ${VER}" >> "${TOOL_VERSION_LOG}"
        fi
    else
        echo "  [MISS] ${tool} not found on PATH"
        MISSING_TOOLS+=("${tool}")
    fi
done

if [[ ${#MISSING_TOOLS[@]} -gt 0 ]]; then
    echo ""
    echo "[FATAL] Missing required tools: ${MISSING_TOOLS[*]}"
    echo "        On the HPC, try: module avail | grep -iE 'bcftools|bedtools|htslib|R/'"
    echo "        then: module load <name> for each, and re-source this script."
    return 1 2>/dev/null || exit 1
fi

echo ""
echo "Tool versions recorded to: ${TOOL_VERSION_LOG}"

# --- Create working directories ---------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"/{gene_fusions,translocations,inversions,deletions,insertions,duplications,complex_rearrangements,qc_summary}

echo ""
echo "=== [00_setup_env.sh] Environment ready. ==="
