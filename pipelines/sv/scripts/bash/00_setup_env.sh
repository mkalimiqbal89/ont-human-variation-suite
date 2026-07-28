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
#
# Safe to source from bash OR zsh (the macOS default login shell). Stages 01+
# are bash scripts and should be invoked with `bash scripts/bash/NN_*.sh`.
#
# The generic config-reading and tool-checking logic lives in
# common/lib_common.sh, shared with every other pipeline in the suite. Only
# what is specific to structural variants stays here.
# =============================================================================

set -uo pipefail  # no -e here: we want to report ALL missing tools, not exit on first

# --- Bootstrap: find ourselves, then the suite root --------------------------
# This block cannot live in the shared library, because locating the library is
# exactly what it does.
#
# BASH_SOURCE does not exist in zsh, so sourcing this from an interactive macOS
# zsh prompt would otherwise resolve SCRIPT_DIR to the caller's cwd and look for
# the config in the wrong place. zsh sets $0 to the script path when sourcing
# (FUNCTION_ARGZERO, on by default), which covers the gap without any
# zsh-specific syntax that bash could not parse.
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
    config_file \
    sv_vcf cnv_vcf snv_vcf str_vcf \
    min_sv_length min_read_support min_qual min_vaf \
    threads log_dir archive_root \
    || { return 1 2>/dev/null || exit 1; }

# --- Export config values ----------------------------------------------------
SAMPLE_ID="$(common_yaml_get "${CONFIG_FILE}" 'sample_id')"
RAW_SAMPLE_PREFIX="$(common_yaml_get "${CONFIG_FILE}" 'raw_sample_prefix')"
INPUT_DIR="$(common_yaml_get "${CONFIG_FILE}" 'input_dir')"
OUTPUT_DIR="$(common_yaml_get "${CONFIG_FILE}" 'output_dir')"
WORK_DIR="$(common_yaml_get "${CONFIG_FILE}" 'work_dir')"
REPO_DIR_CFG="$(common_yaml_get "${CONFIG_FILE}" 'repo_dir')"
MIN_SV_LENGTH="$(common_yaml_get "${CONFIG_FILE}" 'min_sv_length')"
MIN_READ_SUPPORT="$(common_yaml_get "${CONFIG_FILE}" 'min_read_support')"
MIN_QUAL="$(common_yaml_get "${CONFIG_FILE}" 'min_qual')"
MIN_VAF="$(common_yaml_get "${CONFIG_FILE}" 'min_vaf')"
THREADS="$(common_yaml_get "${CONFIG_FILE}" 'threads')"
ARCHIVE_ROOT="$(common_yaml_get "${CONFIG_FILE}" 'archive_root')"
LOG_DIR="$(common_resolve_dir "$(common_yaml_get "${CONFIG_FILE}" 'log_dir')" "${REPO_DIR}")"

export SAMPLE_ID RAW_SAMPLE_PREFIX INPUT_DIR OUTPUT_DIR WORK_DIR REPO_DIR_CFG
export MIN_SV_LENGTH MIN_READ_SUPPORT MIN_QUAL MIN_VAF THREADS ARCHIVE_ROOT LOG_DIR

SV_VCF_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'sv_vcf'  "${RAW_SAMPLE_PREFIX}")"
CNV_VCF_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'cnv_vcf' "${RAW_SAMPLE_PREFIX}")"
SNV_VCF_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'snv_vcf' "${RAW_SAMPLE_PREFIX}")"
STR_VCF_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'str_vcf' "${RAW_SAMPLE_PREFIX}")"
export SV_VCF_NAME CNV_VCF_NAME SNV_VCF_NAME STR_VCF_NAME

echo "Sample ID     : ${SAMPLE_ID}"
echo "Input dir     : ${INPUT_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Work dir      : ${WORK_DIR}"

# --- Tools -------------------------------------------------------------------
echo ""
echo "=== Tool check ==="
mkdir -p "${LOG_DIR}"
TOOL_VERSION_LOG="${LOG_DIR}/tool_versions.txt"
export TOOL_VERSION_LOG

common_check_tools "${TOOL_VERSION_LOG}" bcftools bedtools tabix Rscript \
    || { return 1 2>/dev/null || exit 1; }

echo ""
echo "Tool versions recorded to: ${TOOL_VERSION_LOG}"

# --- Working directories -----------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}/gene_fusions" \
         "${OUTPUT_DIR}/translocations" \
         "${OUTPUT_DIR}/inversions" \
         "${OUTPUT_DIR}/deletions" \
         "${OUTPUT_DIR}/insertions" \
         "${OUTPUT_DIR}/duplications" \
         "${OUTPUT_DIR}/complex_rearrangements" \
         "${OUTPUT_DIR}/qc_summary"

echo ""
echo "=== [00_setup_env.sh] Environment ready. ==="
