#!/usr/bin/env bash
# =============================================================================
# 01_validate_inputs.sh (pipelines/cnv)
# Validates input Spectre CNV VCF files and reference resource files before
# executing analysis stages. Refuses to proceed if required files are missing
# or corrupted.
#
# Usage:
#   bash scripts/bash/01_validate_inputs.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Source environment
# shellcheck source=00_setup_env.sh
. "${SCRIPT_DIR}/00_setup_env.sh" "${1:-}"

echo ""
echo "=== [01_validate_inputs.sh] Validating CNV inputs ==="

ERRORS=0

check_vcf() {
    local label="$1"
    local path="$2"
    local optional="${3:-false}"

    if [[ ! -f "${path}" ]]; then
        if [[ "${optional}" == "true" ]]; then
            echo "  [INFO] ${label}: not found (optional file, skipping)"
        else
            echo "  [FAIL] ${label}: file not found: ${path}"
            ERRORS=$((ERRORS+1))
        fi
        return
    fi

    local size
    size="$(wc -c < "${path}" | tr -d ' ')"
    if [[ "${size}" -eq 0 ]]; then
        echo "  [FAIL] ${label}: file is 0 bytes: ${path}"
        ERRORS=$((ERRORS+1))
        return
    fi

    echo "  [OK]   ${label}: ${path} (${size} bytes)"

    # Tabix index check
    if [[ ! -f "${path}.tbi" && ! -f "${path}.csi" ]]; then
        echo "  [WARN] ${label}: no .tbi/.csi index found — indexing now"
        tabix -p vcf "${path}" 2>/dev/null || echo "  [FAIL] ${label}: tabix indexing failed"
    fi
}

echo ""
echo "--- CNV VCF (Spectre) ---"
check_vcf "CNV VCF" "${INPUT_DIR}/${CNV_VCF_NAME}"

# --- Reference bundle check --------------------------------------------------
echo ""
echo "--- Reference files ---"
REF_CONFIG="$(common_ref_config "${CONFIG_FILE}" "${REPO_DIR}" 2>/dev/null || true)"
if [[ -z "${REF_CONFIG}" ]]; then
    echo "  [FAIL] reference.config_file not set in ${CONFIG_FILE}"
    ERRORS=$((ERRORS+1))
elif [[ ! -f "${REF_CONFIG}" ]]; then
    echo "  [FAIL] Reference config file not found: ${REF_CONFIG}"
    ERRORS=$((ERRORS+1))
else
    GENE_BED="$(grep -E "^[[:space:]]*gene_bed:" "${REF_CONFIG}" | head -n1 | sed -E 's/^[^:]+:[[:space:]]*"?//; s/"?[[:space:]]*$//' || true)"
    if [[ -n "${GENE_BED}" ]]; then
        REF_DIR="$(cd "$(dirname "${REF_CONFIG}")" && pwd)"
        GENE_BED_ABS="$(common_resolve_dir "${GENE_BED}" "${REF_DIR}")"
        if [[ ! -f "${GENE_BED_ABS}" ]]; then
            GENE_BED_ABS="$(common_resolve_dir "${GENE_BED}" "${REPO_DIR}")"
        fi
        if [[ -f "${GENE_BED_ABS}" ]]; then
            echo "  [OK]   GENCODE gene BED: ${GENE_BED_ABS}"
        else
            echo "  [WARN] GENCODE gene BED specified but not found: ${GENE_BED_ABS} (annotation step will handle missing file gracefully)"
        fi
    fi
fi


echo ""
if [[ "${ERRORS}" -gt 0 ]]; then
    echo "=== [01_validate_inputs.sh] Validation FAILED with ${ERRORS} error(s) ==="
    exit 1
fi

echo "=== [01_validate_inputs.sh] Validation PASSED ==="
