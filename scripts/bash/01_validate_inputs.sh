#!/usr/bin/env bash
# =============================================================================
# 01_validate_inputs.sh
# Confirms the expected Epi2ME wf-human-variation outputs exist, are valid
# bgzipped/tabixed VCFs, and reports basic sanity stats (variant counts,
# contig list) before any filtering happens. Fails loudly and specifically
# so problems are caught before hours are spent downstream.
#
# Usage:
#   bash scripts/bash/01_validate_inputs.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONFIG_FILE="${1:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/01_validate_inputs_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo ""
echo "=== [01_validate_inputs.sh] Validating inputs for sample: ${SAMPLE_ID} ==="

ERRORS=0

check_vcf() {
    local label="$1"
    local path="$2"
    local skip_count="${3:-false}"

    if [[ ! -f "${path}" ]]; then
        echo "  [MISS] ${label}: ${path} not found"
        ERRORS=$((ERRORS+1))
        return
    fi

    # Confirm it's actually bgzipped (Epi2ME VCFs are bgzip, not plain gzip)
    if ! file "${path}" | grep -qi "gzip"; then
        echo "  [WARN] ${label}: does not look gzipped — check format: ${path}"
    fi

    # bcftools view -h will fail cleanly on a corrupt/non-VCF file
    if ! bcftools view -h "${path}" > /dev/null 2>&1; then
        echo "  [FAIL] ${label}: not a readable VCF (bcftools could not parse header): ${path}"
        ERRORS=$((ERRORS+1))
        return
    fi

    if [[ "${skip_count}" == "true" ]]; then
        echo "  [OK]   ${label}: ${path}"
        echo "         (variant count skipped — this pipeline doesn't process this VCF downstream; counting ~5M+ records is slow and not worth it here)"
    else
        local n_variants
        n_variants=$(bcftools view -H "${path}" 2>/dev/null | wc -l)
        local n_contigs
        n_contigs=$(bcftools view -H "${path}" 2>/dev/null | cut -f1 | sort -u | wc -l)

        echo "  [OK]   ${label}: ${path}"
        echo "         variants=${n_variants}  contigs_with_calls=${n_contigs}"
    fi

    # Tabix index check (needed for region queries later)
    if [[ ! -f "${path}.tbi" && ! -f "${path}.csi" ]]; then
        echo "  [WARN] ${label}: no .tbi/.csi index found — indexing now"
        tabix -p vcf "${path}" 2>/dev/null || echo "  [FAIL] ${label}: tabix indexing failed"
    fi
}

echo ""
echo "--- SV VCF (Sniffles2) ---"
check_vcf "SV VCF" "${INPUT_DIR}/${SV_VCF_NAME}"

echo ""
echo "--- CNV VCF (Spectre) ---"
check_vcf "CNV VCF" "${INPUT_DIR}/${CNV_VCF_NAME}"

echo ""
echo "--- SNV VCF (Clair3) ---"
check_vcf "SNV VCF" "${INPUT_DIR}/${SNV_VCF_NAME}" true

# --- Reference bundle check --------------------------------------------------
echo ""
echo "--- Reference files (per reference.config_file in pipeline_config.yaml) ---"
REF_CONFIG_RELATIVE="$(grep -E '^[[:space:]]*config_file:' "${CONFIG_FILE}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//')"
if [[ -z "${REF_CONFIG_RELATIVE}" ]]; then
    echo "  [FAIL] reference.config_file not set in ${CONFIG_FILE}"
    ERRORS=$((ERRORS+1))
    REF_CONFIG=""
elif [[ "${REF_CONFIG_RELATIVE}" = /* ]]; then
    REF_CONFIG="${REF_CONFIG_RELATIVE}"
else
    REF_CONFIG="${REPO_DIR}/${REF_CONFIG_RELATIVE}"
fi

if [[ -n "${REF_CONFIG}" && ! -f "${REF_CONFIG}" ]]; then
    echo "  [FAIL] Reference config file not found: ${REF_CONFIG}"
    ERRORS=$((ERRORS+1))
    REF_CONFIG=""
fi
if [[ -n "${REF_CONFIG}" ]]; then
    ref_get() {
        grep -E "^[[:space:]]*${1}:" "${REF_CONFIG}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//'
    }
    GENOME_FASTA="$(ref_get 'fasta')"
    GENE_BED="$(ref_get 'gene_bed')"

    for f in "${GENOME_FASTA}" "${GENE_BED}"; do
        if [[ -f "${f}" ]]; then
            echo "  [OK]   $(basename "${f}") found"
        else
            echo "  [MISS] ${f} not found (required before scripts/R/05_annotate_variants.R)"
            ERRORS=$((ERRORS+1))
        fi
    done
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
