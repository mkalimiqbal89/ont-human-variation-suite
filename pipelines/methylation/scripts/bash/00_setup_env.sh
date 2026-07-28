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
# Safe to source from bash OR zsh (macOS default login shell). Stages 01+ are
# bash scripts and should be invoked with `bash scripts/bash/NN_*.sh`.
# =============================================================================

set -uo pipefail  # no -e here: we want to report ALL missing tools, not exit on first

# --- Resolve repo root regardless of where this is called from -------------
# BASH_SOURCE does not exist in zsh, so sourcing this from an interactive macOS
# zsh prompt would otherwise resolve SCRIPT_DIR to the caller's cwd and look for
# the config in the wrong place. zsh sets $0 to the script path when sourcing
# (FUNCTION_ARGZERO, on by default), which covers the gap without needing any
# zsh-specific syntax that bash could not parse.
if [ -n "${BASH_SOURCE:-}" ]; then
    _SELF="${BASH_SOURCE[0]}"
else
    _SELF="$0"
fi

if [ ! -e "${_SELF}" ]; then
    echo "[FATAL] Cannot resolve the path to 00_setup_env.sh (got: '${_SELF}')." >&2
    echo "        Source it by path from the repo root, e.g.:" >&2
    echo "          source scripts/bash/00_setup_env.sh" >&2
    return 1 2>/dev/null || exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${_SELF}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# --- Argument sanity ---------------------------------------------------------
# zsh does NOT strip inline `#` comments in interactive shells, unlike bash. So
# copying a documented command such as
#     source scripts/bash/00_setup_env.sh   # check the environment
# from a README into a zsh prompt delivers "#" here as $1, and the only symptom
# is a baffling "config file not found: #". No legitimate config path begins
# with '#', so drop it and say why.
_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in
    '#'*)
        echo "[WARN] Ignoring argument '${_CFG_ARG}' — that looks like a shell comment."
        echo "       zsh does not strip inline '#' comments in interactive shells."
        echo "       Run 'setopt interactive_comments' for bash-like behaviour, or"
        echo "       just paste the command without the trailing comment."
        _CFG_ARG=""
        ;;
esac

CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "[FATAL] Config file not found: ${CONFIG_FILE}" >&2
    return 1 2>/dev/null || exit 1
fi

echo "=== [00_setup_env.sh] Using config: ${CONFIG_FILE} ==="

# --- Minimal YAML value getter (avoids requiring yq on the HPC) ------------
# Reads a simple "key: value" line under a given top-level/nested key path.
# For this pipeline's flat-ish config, a targeted grep/sed is sufficient and
# avoids adding a yq dependency requirement to the HPC environment.
#
# CONSEQUENCE: this returns the FIRST match anywhere in the file, so every key
# read here must be unique across the whole config. guard_duplicate_keys()
# below enforces that rather than trusting it — a duplicate key silently
# returning the wrong value is exactly the kind of bug that is invisible until
# a result is wrong.
yaml_get() {
    local key="$1"
    grep -E "^[[:space:]]*${key}:" "${CONFIG_FILE}" | head -n1 \
        | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//'
}

guard_duplicate_keys() {
    local key dupes=()
    for key in "$@"; do
        local n
        n=$(grep -cE "^[[:space:]]*${key}:" "${CONFIG_FILE}")
        if [[ "${n}" -gt 1 ]]; then
            dupes+=("${key} (x${n})")
        fi
    done
    if [[ ${#dupes[@]} -gt 0 ]]; then
        echo "[FATAL] Duplicate config keys in ${CONFIG_FILE}: ${dupes[*]}" >&2
        echo "        yaml_get() returns the first match, so one of these is" >&2
        echo "        being silently ignored. Fix the config before running." >&2
        return 1
    fi
    return 0
}

CONFIG_KEYS=(
    sample_id raw_sample_prefix
    input_dir output_dir work_dir repo_dir
    config_file
    mod_bedmethyl mod_bedmethyl_hap1 mod_bedmethyl_hap2 mod_bedmethyl_ungrouped
    phased primary_mod_code expected_columns
    min_coverage max_coverage
    primary_contigs_only include_contigs_regex expected_primary_contigs
    exclude_contigs_regex
    unmethylated_max_percent methylated_min_percent
    min_cpgs_per_feature gene_key keep_annotated_cpgs
    threads log_dir archive_root compress include_intermediates
)
guard_duplicate_keys "${CONFIG_KEYS[@]}" || { return 1 2>/dev/null || exit 1; }

# --- Export config values as shell variables --------------------------------
export SAMPLE_ID="$(yaml_get 'sample_id')"
export RAW_SAMPLE_PREFIX="$(yaml_get 'raw_sample_prefix')"
export INPUT_DIR="$(yaml_get 'input_dir')"
export OUTPUT_DIR="$(yaml_get 'output_dir')"
export WORK_DIR="$(yaml_get 'work_dir')"
export REPO_DIR_CFG="$(yaml_get 'repo_dir')"
export PHASED="$(yaml_get 'phased')"
export PRIMARY_MOD_CODE="$(yaml_get 'primary_mod_code')"
export EXPECTED_COLUMNS="$(yaml_get 'expected_columns')"
export MIN_COVERAGE="$(yaml_get 'min_coverage')"
export MAX_COVERAGE="$(yaml_get 'max_coverage')"
export EXCLUDE_CONTIGS_REGEX="$(yaml_get 'exclude_contigs_regex')"

# Contig selection. Defaults are applied for any key a pre-existing config
# predates, so an older pipeline_config.yaml keeps working instead of silently
# resolving to an empty regex — which, for an allowlist, would discard the whole
# genome.
PRIMARY_CONTIGS_ONLY="$(yaml_get 'primary_contigs_only')"
[[ -z "${PRIMARY_CONTIGS_ONLY}" ]] && PRIMARY_CONTIGS_ONLY="true"
export PRIMARY_CONTIGS_ONLY

INCLUDE_CONTIGS_REGEX="$(yaml_get 'include_contigs_regex')"
if [[ -z "${INCLUDE_CONTIGS_REGEX}" ]]; then
    INCLUDE_CONTIGS_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'
    echo "[INFO] filtering.include_contigs_regex not set; defaulting to chr1-22,X,Y"
fi
export INCLUDE_CONTIGS_REGEX

EXPECTED_PRIMARY_CONTIGS="$(yaml_get 'expected_primary_contigs')"
[[ -z "${EXPECTED_PRIMARY_CONTIGS}" ]] && EXPECTED_PRIMARY_CONTIGS="24"
export EXPECTED_PRIMARY_CONTIGS
UNMETHYLATED_MAX_PERCENT="$(yaml_get 'unmethylated_max_percent')"
[[ -z "${UNMETHYLATED_MAX_PERCENT}" ]] && UNMETHYLATED_MAX_PERCENT="20"
export UNMETHYLATED_MAX_PERCENT

METHYLATED_MIN_PERCENT="$(yaml_get 'methylated_min_percent')"
[[ -z "${METHYLATED_MIN_PERCENT}" ]] && METHYLATED_MIN_PERCENT="80"
export METHYLATED_MIN_PERCENT

MIN_CPGS_PER_FEATURE="$(yaml_get 'min_cpgs_per_feature')"
[[ -z "${MIN_CPGS_PER_FEATURE}" ]] && MIN_CPGS_PER_FEATURE="5"
export MIN_CPGS_PER_FEATURE

export GENE_KEY="$(yaml_get 'gene_key')"

KEEP_ANNOTATED_CPGS="$(yaml_get 'keep_annotated_cpgs')"
[[ -z "${KEEP_ANNOTATED_CPGS}" ]] && KEEP_ANNOTATED_CPGS="false"
export KEEP_ANNOTATED_CPGS
export THREADS="$(yaml_get 'threads')"
export ARCHIVE_ROOT="$(yaml_get 'archive_root')"

ARCHIVE_COMPRESS="$(yaml_get 'compress')"
[[ -z "${ARCHIVE_COMPRESS}" ]] && ARCHIVE_COMPRESS="false"
export ARCHIVE_COMPRESS

ARCHIVE_INCLUDE_INTERMEDIATES="$(yaml_get 'include_intermediates')"
[[ -z "${ARCHIVE_INCLUDE_INTERMEDIATES}" ]] && ARCHIVE_INCLUDE_INTERMEDIATES="false"
export ARCHIVE_INCLUDE_INTERMEDIATES

# LOG_DIR: only prepend REPO_DIR when the configured value is relative.
# Prepending unconditionally leaks test scratch dirs into the real repo.
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
export MOD_BEDMETHYL_NAME="$(resolve_filename 'mod_bedmethyl')"
export MOD_BEDMETHYL_HAP1_NAME="$(resolve_filename 'mod_bedmethyl_hap1')"
export MOD_BEDMETHYL_HAP2_NAME="$(resolve_filename 'mod_bedmethyl_hap2')"
export MOD_BEDMETHYL_UNGROUPED_NAME="$(resolve_filename 'mod_bedmethyl_ungrouped')"
export MOSDEPTH_SUMMARY_NAME="$(resolve_filename 'mosdepth_summary')"

# --- Resolve the reference config ------------------------------------------
# Read the path from the pipeline config rather than hardcoding it, so a test
# config can point at test references.
REF_CONFIG_RAW="$(yaml_get 'config_file')"
if [[ -z "${REF_CONFIG_RAW}" ]]; then
    echo "[FATAL] reference.config_file not set in ${CONFIG_FILE}" >&2
    return 1 2>/dev/null || exit 1
elif [[ "${REF_CONFIG_RAW}" = /* ]]; then
    export REF_CONFIG="${REF_CONFIG_RAW}"
else
    export REF_CONFIG="${REPO_DIR}/${REF_CONFIG_RAW}"
fi

if [[ ! -f "${REF_CONFIG}" ]]; then
    echo "[FATAL] Reference config not found: ${REF_CONFIG}" >&2
    return 1 2>/dev/null || exit 1
fi

ref_get() {
    grep -E "^[[:space:]]*${1}:" "${REF_CONFIG}" | head -n1 \
        | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//'
}
export GENOME_BUILD="$(ref_get 'build')"
export CHROM_SIZES="$(ref_get 'chrom_sizes')"
export GENE_BED="$(ref_get 'gene_bed')"
export PROMOTER_BED="$(ref_get 'promoter_bed')"
export CPG_ISLAND_BED="$(ref_get 'cpg_island_bed')"
export GENCODE_GTF="$(ref_get 'gencode_gtf')"

echo "Sample ID     : ${SAMPLE_ID}"
echo "Raw prefix    : ${RAW_SAMPLE_PREFIX}"
echo "Input dir     : ${INPUT_DIR}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Work dir      : ${WORK_DIR}"
echo "Genome build  : ${GENOME_BUILD}"
echo "Phased mods   : ${PHASED}"
echo "Primary mod   : ${PRIMARY_MOD_CODE}   (min coverage ${MIN_COVERAGE})"

# --- Check required tools ---------------------------------------------------
# bedtools : all annotation intersects
# tabix    : region queries on the bgzipped bedMethyl
# bgzip    : recompressing filtered output as bgzip so it can be tabix-indexed
# Rscript  : summary/figure stages
# awk/sort : the workhorses of stages 02-06
REQUIRED_TOOLS=(bedtools tabix bgzip Rscript awk sort)
MISSING_TOOLS=()

echo ""
echo "=== Tool check ==="
mkdir -p "${LOG_DIR}"
TOOL_VERSION_LOG="${LOG_DIR}/tool_versions.txt"
: > "${TOOL_VERSION_LOG}"

for tool in "${REQUIRED_TOOLS[@]}"; do
    if command -v "${tool}" >/dev/null 2>&1; then
        # Both the exit status AND the output must be checked: a binary can be
        # present on PATH and still be unusable (e.g. a broken shared library),
        # which reports [OK] if only `command -v` is trusted.
        VER=$("${tool}" --version 2>&1 | head -n1)
        RC=$?
        if [[ ${RC} -ne 0 \
              || "${VER}" == *"cannot open shared object"* \
              || "${VER}" == *"error while loading shared libraries"* \
              || "${VER}" == *"command not found"* ]]; then
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
    echo "        On macOS: brew install bedtools htslib"
    echo "        On the HPC: module avail | grep -iE 'bedtools|htslib|R/'"
    echo "        then: module load <name> for each, and re-source this script."
    return 1 2>/dev/null || exit 1
fi

echo ""
echo "Tool versions recorded to: ${TOOL_VERSION_LOG}"

# --- Create working directories ---------------------------------------------
mkdir -p "${OUTPUT_DIR}" "${WORK_DIR}" "${LOG_DIR}"
mkdir -p "${OUTPUT_DIR}"/{01_filtered,02_qc,03_global_methylation,04_chromosome_summary,05_annotated_cpgs,06_gene_summary,07_promoter_summary,08_cpg_island_summary,09_figures,10_igv_tracks}

echo ""
echo "=== [00_setup_env.sh] Environment ready. ==="
