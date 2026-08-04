#!/usr/bin/env bash
# =============================================================================
# 00_setup_env.sh (pipelines/snv)
#
# SCAFFOLD ONLY. Resolves the generic config fields and does the generic tool
# check every pipeline in the suite needs. There is no stage 01+ here yet, so
# this deliberately does NOT check SNV-specific tools (beyond the common
# VCF-handling baseline) or read SNV-specific filtering config — there is
# nothing downstream to consume them yet. Extend the
# common_guard_duplicate_keys key list and the common_check_tools tool list
# below as real stages are written, following pipelines/sv/scripts/bash/
# 00_setup_env.sh or pipelines/methylation/scripts/bash/00_setup_env.sh as a
# model.
#
# Usage:
#   source scripts/bash/00_setup_env.sh [path/to/pipeline_config.yaml]
#
# NOTE: this script must be SOURCED (not executed), since it exports
# variables into the calling shell. Safe to source from bash or zsh.
# =============================================================================

set -uo pipefail  # no -e here: we want to report ALL missing tools, not exit on first

# --- Bootstrap: find ourselves, then the suite root --------------------------
if [ -n "${BASH_SOURCE:-}" ]; then
    _SELF="${BASH_SOURCE[0]}"
else
    _SELF="$0"
fi
if [ ! -e "${_SELF}" ]; then
    echo "[FATAL] Cannot resolve the path to 00_setup_env.sh (got: '${_SELF}')." >&2
    echo "        Source it by path from the pipeline root, e.g.:" >&2
    echo "          source scripts/bash/00_setup_env.sh" >&2
    return 1 2>/dev/null || exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${_SELF}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"     # the pipeline root
SUITE_DIR="$(cd "${REPO_DIR}/../.." && pwd)"      # the suite root

LIB="${SUITE_DIR}/common/lib_common.sh"
if [ ! -f "${LIB}" ]; then
    echo "[FATAL] Shared library not found: ${LIB}" >&2
    echo "        This pipeline expects to live at <suite>/pipelines/<name>/." >&2
    return 1 2>/dev/null || exit 1
fi
# shellcheck source=../../../../common/lib_common.sh
. "${LIB}"

# --- Config ------------------------------------------------------------------
_CFG_ARG="$(common_config_arg "${1:-}")"
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "[FATAL] Config file not found: ${CONFIG_FILE}" >&2
    return 1 2>/dev/null || exit 1
fi

echo "=== [00_setup_env.sh] Using config: ${CONFIG_FILE} ==="

common_guard_duplicate_keys "${CONFIG_FILE}" \
    sample_id raw_sample_prefix \
    input_dir output_dir work_dir repo_dir \
    config_file snv_vcf snv_vcf_clinvar \
    threads log_dir archive_root \
    || { return 1 2>/dev/null || exit 1; }

# --- Export config values ----------------------------------------------------
SAMPLE_ID="$(common_yaml_get "${CONFIG_FILE}" 'sample_id')"
RAW_SAMPLE_PREFIX="$(common_yaml_get "${CONFIG_FILE}" 'raw_sample_prefix')"
INPUT_DIR="$(common_yaml_get "${CONFIG_FILE}" 'input_dir')"
OUTPUT_DIR="$(common_yaml_get "${CONFIG_FILE}" 'output_dir')"
WORK_DIR="$(common_yaml_get "${CONFIG_FILE}" 'work_dir')"
REPO_DIR_CFG="$(common_yaml_get "${CONFIG_FILE}" 'repo_dir')"
THREADS="$(common_yaml_get "${CONFIG_FILE}" 'threads')"
ARCHIVE_ROOT="$(common_yaml_get "${CONFIG_FILE}" 'archive_root')"
LOG_DIR="$(common_resolve_dir "$(common_yaml_get "${CONFIG_FILE}" 'log_dir')" "${REPO_DIR}")"

export SAMPLE_ID RAW_SAMPLE_PREFIX INPUT_DIR OUTPUT_DIR WORK_DIR REPO_DIR_CFG
export THREADS ARCHIVE_ROOT LOG_DIR

SNV_VCF_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'snv_vcf' "${RAW_SAMPLE_PREFIX}")"
SNV_VCF_CLINVAR_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'snv_vcf_clinvar' "${RAW_SAMPLE_PREFIX}")"
export SNV_VCF_NAME SNV_VCF_CLINVAR_NAME

echo "Sample ID     : ${SAMPLE_ID}"
echo "Input dir     : ${INPUT_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Work dir      : ${WORK_DIR}"

# --- Tools -------------------------------------------------------------------
# Baseline VCF-handling tools only. Add whatever this pipeline's real stages
# turn out to need.
echo ""
echo "=== Tool check ==="
mkdir -p "${LOG_DIR}"
TOOL_VERSION_LOG="${LOG_DIR}/tool_versions.txt"
export TOOL_VERSION_LOG

common_check_tools "${TOOL_VERSION_LOG}" bcftools tabix \
    || { return 1 2>/dev/null || exit 1; }

echo ""
echo "Tool versions recorded to: ${TOOL_VERSION_LOG}"

# --- Working directories -----------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${LOG_DIR}"

echo ""
echo "=== [00_setup_env.sh] Environment ready. (scaffold pipeline — no stage 01+ yet) ==="
