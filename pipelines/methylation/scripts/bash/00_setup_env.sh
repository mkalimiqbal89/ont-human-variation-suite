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
# what is specific to methylation stays here.
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
    mod_bedmethyl mod_bedmethyl_hap1 mod_bedmethyl_hap2 mod_bedmethyl_ungrouped \
    phased primary_mod_code expected_columns \
    min_coverage max_coverage \
    primary_contigs_only include_contigs_regex expected_primary_contigs \
    exclude_contigs_regex \
    unmethylated_max_percent methylated_min_percent \
    min_cpgs_per_feature gene_key keep_annotated_cpgs \
    threads log_dir archive_root compress include_intermediates \
    || { return 1 2>/dev/null || exit 1; }

# cfg <key> [default]
# Reads a config value, falling back to a default. Defaults exist so that a
# config predating a newly added key keeps working rather than silently
# resolving to an empty string — which, for an allowlist regex, would discard
# the entire genome.
cfg() {
    local v
    v="$(common_yaml_get "${CONFIG_FILE}" "$1")"
    [ -z "${v}" ] && v="${2:-}"
    echo "${v}"
}

SAMPLE_ID="$(cfg sample_id)"
RAW_SAMPLE_PREFIX="$(cfg raw_sample_prefix)"
INPUT_DIR="$(cfg input_dir)"
OUTPUT_DIR="$(cfg output_dir)"
WORK_DIR="$(cfg work_dir)"
REPO_DIR_CFG="$(cfg repo_dir)"
PHASED="$(cfg phased false)"
PRIMARY_MOD_CODE="$(cfg primary_mod_code m)"
EXPECTED_COLUMNS="$(cfg expected_columns 18)"
MIN_COVERAGE="$(cfg min_coverage 10)"
MAX_COVERAGE="$(cfg max_coverage 0)"
EXCLUDE_CONTIGS_REGEX="$(cfg exclude_contigs_regex)"
PRIMARY_CONTIGS_ONLY="$(cfg primary_contigs_only true)"
EXPECTED_PRIMARY_CONTIGS="$(cfg expected_primary_contigs 24)"
UNMETHYLATED_MAX_PERCENT="$(cfg unmethylated_max_percent 20)"
METHYLATED_MIN_PERCENT="$(cfg methylated_min_percent 80)"
MIN_CPGS_PER_FEATURE="$(cfg min_cpgs_per_feature 5)"
GENE_KEY="$(cfg gene_key gene_id)"
KEEP_ANNOTATED_CPGS="$(cfg keep_annotated_cpgs false)"
THREADS="$(cfg threads 4)"
ARCHIVE_ROOT="$(cfg archive_root)"
ARCHIVE_COMPRESS="$(cfg compress false)"
ARCHIVE_INCLUDE_INTERMEDIATES="$(cfg include_intermediates false)"

INCLUDE_CONTIGS_REGEX="$(common_yaml_get "${CONFIG_FILE}" 'include_contigs_regex')"
if [ -z "${INCLUDE_CONTIGS_REGEX}" ]; then
    INCLUDE_CONTIGS_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'
    echo "[INFO] filtering.include_contigs_regex not set; defaulting to chr1-22,X,Y"
fi

LOG_DIR="$(common_resolve_dir "$(cfg log_dir logs)" "${REPO_DIR}")"

export SAMPLE_ID RAW_SAMPLE_PREFIX INPUT_DIR OUTPUT_DIR WORK_DIR REPO_DIR_CFG
export PHASED PRIMARY_MOD_CODE EXPECTED_COLUMNS MIN_COVERAGE MAX_COVERAGE
export EXCLUDE_CONTIGS_REGEX PRIMARY_CONTIGS_ONLY INCLUDE_CONTIGS_REGEX
export EXPECTED_PRIMARY_CONTIGS UNMETHYLATED_MAX_PERCENT METHYLATED_MIN_PERCENT
export MIN_CPGS_PER_FEATURE GENE_KEY KEEP_ANNOTATED_CPGS THREADS
export ARCHIVE_ROOT ARCHIVE_COMPRESS ARCHIVE_INCLUDE_INTERMEDIATES LOG_DIR

MOD_BEDMETHYL_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'mod_bedmethyl' "${RAW_SAMPLE_PREFIX}")"
MOD_BEDMETHYL_HAP1_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'mod_bedmethyl_hap1' "${RAW_SAMPLE_PREFIX}")"
MOD_BEDMETHYL_HAP2_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'mod_bedmethyl_hap2' "${RAW_SAMPLE_PREFIX}")"
MOD_BEDMETHYL_UNGROUPED_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'mod_bedmethyl_ungrouped' "${RAW_SAMPLE_PREFIX}")"
MOSDEPTH_SUMMARY_NAME="$(common_resolve_filename "${CONFIG_FILE}" 'mosdepth_summary' "${RAW_SAMPLE_PREFIX}")"
export MOD_BEDMETHYL_NAME MOD_BEDMETHYL_HAP1_NAME MOD_BEDMETHYL_HAP2_NAME
export MOD_BEDMETHYL_UNGROUPED_NAME MOSDEPTH_SUMMARY_NAME

# --- Reference config --------------------------------------------------------
REF_CONFIG="$(common_ref_config "${CONFIG_FILE}" "${REPO_DIR}")"
if [ -z "${REF_CONFIG}" ]; then
    echo "[FATAL] reference.config_file not set in ${CONFIG_FILE}" >&2
    return 1 2>/dev/null || exit 1
fi
if [ ! -f "${REF_CONFIG}" ]; then
    echo "[FATAL] Reference config not found: ${REF_CONFIG}" >&2
    return 1 2>/dev/null || exit 1
fi
export REF_CONFIG

GENOME_BUILD="$(common_yaml_get "${REF_CONFIG}" 'build')"
CHROM_SIZES="$(common_yaml_get "${REF_CONFIG}" 'chrom_sizes')"
GENE_BED="$(common_yaml_get "${REF_CONFIG}" 'gene_bed')"
PROMOTER_BED="$(common_yaml_get "${REF_CONFIG}" 'promoter_bed')"
CPG_ISLAND_BED="$(common_yaml_get "${REF_CONFIG}" 'cpg_island_bed')"
GENCODE_GTF="$(common_yaml_get "${REF_CONFIG}" 'gencode_gtf')"
export GENOME_BUILD CHROM_SIZES GENE_BED PROMOTER_BED CPG_ISLAND_BED GENCODE_GTF

echo "Sample ID     : ${SAMPLE_ID}"
echo "Raw prefix    : ${RAW_SAMPLE_PREFIX}"
echo "Input dir     : ${INPUT_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Work dir      : ${WORK_DIR}"
echo "Genome build  : ${GENOME_BUILD}"
echo "Phased mods   : ${PHASED}"
echo "Primary mod   : ${PRIMARY_MOD_CODE}   (min coverage ${MIN_COVERAGE})"

# --- Tools -------------------------------------------------------------------
# bedtools : all annotation intersects
# tabix    : region queries on the bgzipped bedMethyl
# bgzip    : recompressing filtered output so it can be tabix-indexed
# Rscript  : figure and comparison stages
# awk/sort : the workhorses of stages 02-06
echo ""
echo "=== Tool check ==="
mkdir -p "${LOG_DIR}"
TOOL_VERSION_LOG="${LOG_DIR}/tool_versions.txt"
export TOOL_VERSION_LOG

common_check_tools "${TOOL_VERSION_LOG}" bedtools tabix bgzip Rscript awk sort \
    || { return 1 2>/dev/null || exit 1; }

echo ""
echo "Tool versions recorded to: ${TOOL_VERSION_LOG}"

# --- Working directories -----------------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}/01_filtered" \
         "${OUTPUT_DIR}/02_qc" \
         "${OUTPUT_DIR}/03_global_methylation" \
         "${OUTPUT_DIR}/04_chromosome_summary" \
         "${OUTPUT_DIR}/05_annotated_cpgs" \
         "${OUTPUT_DIR}/06_gene_summary" \
         "${OUTPUT_DIR}/07_promoter_summary" \
         "${OUTPUT_DIR}/08_cpg_island_summary" \
         "${OUTPUT_DIR}/09_figures" \
         "${OUTPUT_DIR}/10_igv_tracks"

echo ""
echo "=== [00_setup_env.sh] Environment ready. ==="
