#!/usr/bin/env bash
# =============================================================================
# 01_validate_inputs.sh
# Confirms the expected Epi2ME wf-human-variation methylation outputs exist,
# are readable, have the expected bedMethyl column layout, contain the
# configured modification codes, and use contig names that actually match the
# reference annotation. Fails loudly and specifically so problems are caught
# before hours are spent downstream.
#
# Usage:
#   bash scripts/bash/01_validate_inputs.sh [path/to/pipeline_config.yaml]
#
# Environment overrides:
#   SKIP_INTEGRITY=true   skip the full-file gzip integrity test (which must
#                         decompress the entire 600 MB+ file). Use only when
#                         re-running validation repeatedly during development.
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# Drop a '#'-leading argument: zsh passes inline comments through as arguments.
# See the note in 00_setup_env.sh.
_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in '#'*) _CFG_ARG="" ;; esac
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"
SKIP_INTEGRITY="${SKIP_INTEGRITY:-false}"

# Lines sampled from the head of the file for the cheap column/mod-code checks.
SAMPLE_LINES=200000

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/01_validate_inputs_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo ""
echo "=== [01_validate_inputs.sh] Validating inputs for sample: ${SAMPLE_ID} ==="

ERRORS=0

# --- Helpers -----------------------------------------------------------------

# Reports whether a contig name is UCSC-style ("chr1") or Ensembl-style ("1").
# A mismatch between the bedMethyl and the annotation BEDs makes every
# `bedtools intersect` return zero rows — which looks like "this sample has no
# methylation in genes" rather than like an error. Cheap to check, so check it.
contig_style() {
    local name="$1"
    if [[ "${name}" == chr* ]]; then echo "ucsc"; else echo "ensembl"; fi
}

first_contig_plain() {
    head -n1 "$1" 2>/dev/null | cut -f1
}

first_contig_gz() {
    gzip -cd "$1" 2>/dev/null | head -n1 | cut -f1
}

# --- bedMethyl presence ------------------------------------------------------
echo ""
echo "--- bedMethyl inputs (modkit pileup via wf-human-variation) ---"
echo "  config says phased=${PHASED}"

BEDMETHYL_FILES=()
if [[ "${PHASED}" == "true" ]]; then
    BEDMETHYL_FILES=(
        "${INPUT_DIR}/${MOD_BEDMETHYL_HAP1_NAME}"
        "${INPUT_DIR}/${MOD_BEDMETHYL_HAP2_NAME}"
        "${INPUT_DIR}/${MOD_BEDMETHYL_UNGROUPED_NAME}"
    )
else
    BEDMETHYL_FILES=("${INPUT_DIR}/${MOD_BEDMETHYL_NAME}")
fi

# Catch the case where the config and the actual output disagree, in either
# direction. Getting this wrong silently analyses the wrong file.
UNEXPECTED=()
if [[ "${PHASED}" == "true" ]]; then
    [[ -f "${INPUT_DIR}/${MOD_BEDMETHYL_NAME}" ]] && \
        UNEXPECTED+=("${MOD_BEDMETHYL_NAME} (combined file present but phased=true)")
else
    for hap in "${MOD_BEDMETHYL_HAP1_NAME}" "${MOD_BEDMETHYL_HAP2_NAME}" "${MOD_BEDMETHYL_UNGROUPED_NAME}"; do
        [[ -f "${INPUT_DIR}/${hap}" ]] && \
            UNEXPECTED+=("${hap} (haplotype file present but phased=false)")
    done
fi
if [[ ${#UNEXPECTED[@]} -gt 0 ]]; then
    echo "  [WARN] Config/output mismatch — these files exist but are not being used:"
    for u in "${UNEXPECTED[@]}"; do echo "         ${u}"; done
    echo "         Check modifications.phased in ${CONFIG_FILE}"
fi

PRIMARY_BEDMETHYL=""

for BM in "${BEDMETHYL_FILES[@]}"; do
    LABEL="$(basename "${BM}")"
    echo ""
    echo "  ${LABEL}"

    if [[ ! -f "${BM}" ]]; then
        echo "    [MISS] not found: ${BM}"
        ERRORS=$((ERRORS+1))
        continue
    fi

    SIZE=$(ls -lhL "${BM}" | awk '{print $5}')
    echo "    [OK]   present (${SIZE})"
    [[ -z "${PRIMARY_BEDMETHYL}" ]] && PRIMARY_BEDMETHYL="${BM}"

    # BGZF detection: a bgzipped file can be tabix-indexed for region queries,
    # a plain gzip file cannot. Informational, not an error.
    MAGIC=$(od -An -tx1 -N4 "${BM}" 2>/dev/null | tr -d ' \n')
    if [[ "${MAGIC}" == "1f8b0804" ]]; then
        echo "    [INFO] BGZF-compressed — can be tabix-indexed for region queries"
    else
        echo "    [INFO] plain gzip (magic ${MAGIC}) — recompress with bgzip if"
        echo "           region queries are needed later"
    fi

    # Full integrity test. This decompresses the entire file, so it is the
    # single most expensive check in this script — hence the opt-out.
    if [[ "${SKIP_INTEGRITY}" == "true" ]]; then
        echo "    [SKIP] gzip integrity test (SKIP_INTEGRITY=true)"
    else
        echo -n "    ...... gzip integrity test (decompresses whole file, please wait) "
        T0=$(date +%s)
        if gzip -t "${BM}" 2>/dev/null; then
            T1=$(date +%s)
            echo "-> [OK] ($((T1-T0))s)"
        else
            echo "-> [FAIL]"
            echo "    [FAIL] corrupt or truncated gzip stream: ${BM}"
            ERRORS=$((ERRORS+1))
            continue
        fi
    fi

    # Column layout and modification codes, from a head sample.
    #
    # This is deliberately a SAMPLE, not a full-file scan: a full NF check on a
    # whole-genome bedMethyl costs minutes, and 02_bedmethyl_to_tsv.sh already
    # reads every record in one pass. The full-file column assertion belongs
    # there, where it is free, rather than being paid for twice.
    echo "    ...... column + mod-code check on first ${SAMPLE_LINES} records"
    gzip -cd "${BM}" 2>/dev/null | head -n "${SAMPLE_LINES}" | awk -F'\t' \
        -v expected="${EXPECTED_COLUMNS}" -v primary="${PRIMARY_MOD_CODE}" '
    {
        n++
        if (NF != expected) bad_nf++
        code[$4]++
        if ($10 !~ /^[0-9]+$/)        bad_cov++
        if ($11 !~ /^[0-9.]+$/)       bad_pct++
        if ($11 + 0 < 0 || $11 + 0 > 100) bad_pct_range++
    }
    END {
        printf "    [%s]   records sampled=%d  wrong_column_count=%d (expected %d cols)\n",
               (bad_nf ? "FAIL" : "OK"), n, bad_nf+0, expected
        printf "    [%s]   non-numeric coverage=%d  non-numeric percent=%d  percent_out_of_range=%d\n",
               ((bad_cov + bad_pct + bad_pct_range) ? "FAIL" : "OK"),
               bad_cov+0, bad_pct+0, bad_pct_range+0
        printf "    [INFO] modification codes in sample:"
        for (k in code) printf " %s=%d", k, code[k]
        printf "\n"
        if (!(primary in code))
            printf "    [FAIL] primary_mod_code \"%s\" not present in sample\n", primary
        exit ((bad_nf + bad_cov + bad_pct + bad_pct_range) || !(primary in code)) ? 1 : 0
    }'
    if [[ ${PIPESTATUS[2]} -ne 0 ]]; then
        ERRORS=$((ERRORS+1))
    fi
done

# --- Optional context file ---------------------------------------------------
echo ""
echo "--- Alignment coverage context (optional) ---"
MOSDEPTH="${INPUT_DIR}/${MOSDEPTH_SUMMARY_NAME}"
if [[ -f "${MOSDEPTH}" ]]; then
    TOTAL_MEAN=$(awk -F'\t' '$1=="total" {print $4}' "${MOSDEPTH}" | head -n1)
    echo "  [OK]   ${MOSDEPTH_SUMMARY_NAME}: genome mean depth = ${TOTAL_MEAN}x"
    echo "         (sanity check: filtering.min_coverage=${MIN_COVERAGE} should be"
    echo "          comfortably below this, or most CpGs will be discarded)"
else
    echo "  [WARN] ${MOSDEPTH_SUMMARY_NAME} not found — skipping depth cross-check"
fi

# --- Reference bundle --------------------------------------------------------
echo ""
echo "--- Reference files (per reference.config_file: ${REF_CONFIG}) ---"

# label -> path -> expected column count
check_bed() {
    local label="$1" path="$2" expected_cols="$3"
    if [[ -z "${path}" ]]; then
        echo "  [FAIL] ${label}: not set in ${REF_CONFIG}"
        ERRORS=$((ERRORS+1)); return
    fi
    if [[ ! -f "${path}" ]]; then
        echo "  [MISS] ${label}: ${path} not found"
        ERRORS=$((ERRORS+1)); return
    fi
    local n_rows bad
    read -r n_rows bad < <(awk -F'\t' -v e="${expected_cols}" '
        NF != e { bad++ } END { print NR, bad+0 }' "${path}")
    if [[ "${bad}" -gt 0 ]]; then
        echo "  [FAIL] ${label}: ${bad}/${n_rows} rows do not have ${expected_cols} tab-separated columns"
        ERRORS=$((ERRORS+1))
    else
        echo "  [OK]   ${label}: ${n_rows} rows, all ${expected_cols} columns  ($(basename "${path}"))"
    fi
}

check_bed "gene_bed"        "${GENE_BED}"        7
check_bed "promoter_bed"    "${PROMOTER_BED}"    7
check_bed "cpg_island_bed"  "${CPG_ISLAND_BED}"  10
check_bed "chrom_sizes"     "${CHROM_SIZES}"     2

# gene_bed and promoter_bed are generated from the same gene set, so a row
# count mismatch means the promoter file is stale relative to the gene file.
if [[ -f "${GENE_BED}" && -f "${PROMOTER_BED}" ]]; then
    # tr -d ' ': BSD wc pads its output with leading spaces
    NG=$(wc -l < "${GENE_BED}" | tr -d ' ')
    NP=$(wc -l < "${PROMOTER_BED}" | tr -d ' ')
    if [[ "${NG}" -ne "${NP}" ]]; then
        echo "  [WARN] gene_bed has ${NG} rows but promoter_bed has ${NP}"
        echo "         promoters_2kb.bed is normally 1:1 with genes.bed. If these"
        echo "         differ, one was regenerated without the other."
    else
        echo "  [OK]   gene_bed and promoter_bed are 1:1 (${NG} rows each)"
    fi
fi

if [[ -n "${GENCODE_GTF}" && ! -f "${GENCODE_GTF}" ]]; then
    echo "  [WARN] gencode_gtf not found: ${GENCODE_GTF} (provenance only, not required)"
fi

# --- Contig naming concordance ----------------------------------------------
echo ""
echo "--- Contig naming concordance ---"
if [[ -n "${PRIMARY_BEDMETHYL}" ]]; then
    BM_CONTIG="$(first_contig_gz "${PRIMARY_BEDMETHYL}")"
    declare -a STYLES=()
    declare -a SOURCES=()

    STYLES+=("$(contig_style "${BM_CONTIG}")");  SOURCES+=("bedMethyl (${BM_CONTIG})")
    for pair in "gene_bed:${GENE_BED}" "promoter_bed:${PROMOTER_BED}" \
                "cpg_island_bed:${CPG_ISLAND_BED}" "chrom_sizes:${CHROM_SIZES}"; do
        nm="${pair%%:*}"; pth="${pair#*:}"
        if [[ -f "${pth}" ]]; then
            c="$(first_contig_plain "${pth}")"
            STYLES+=("$(contig_style "${c}")"); SOURCES+=("${nm} (${c})")
        fi
    done

    for i in "${!STYLES[@]}"; do
        printf "  %-28s -> %s\n" "${SOURCES[$i]}" "${STYLES[$i]}"
    done

    UNIQUE_STYLES=$(printf '%s\n' "${STYLES[@]}" | sort -u | wc -l | tr -d ' ')
    if [[ "${UNIQUE_STYLES}" -ne 1 ]]; then
        echo "  [FAIL] contig naming styles disagree."
        echo "         Every bedtools intersect downstream would return ZERO rows,"
        echo "         which reads as 'no methylation in genes' rather than as an"
        echo "         error. Harmonise the naming before continuing."
        ERRORS=$((ERRORS+1))
    else
        echo "  [OK]   all sources use consistent contig naming"
    fi
else
    echo "  [SKIP] no readable bedMethyl file to compare against"
fi

# --- Summary ------------------------------------------------------------------
echo ""
if [[ ${ERRORS} -gt 0 ]]; then
    echo "=== [01_validate_inputs.sh] FAILED with ${ERRORS} error(s). See log: ${LOG_FILE} ==="
    exit 1
else
    echo "=== [01_validate_inputs.sh] All required inputs validated OK. ==="
    echo "Log saved to: ${LOG_FILE}"
fi
