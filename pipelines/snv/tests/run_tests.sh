#!/usr/bin/env bash
# =============================================================================
# run_tests.sh (pipelines/snv/tests)
# Regression test suite for SNV pipeline (pipelines/snv). Runs 04_run_all.sh
# against synthetic fixture test_sample.wf_snp.vcf.gz and asserts exact
# content-level record counts, Ti/Tv ratio, ClinVar pathogenic extraction, and
# figure generation.
#
# Usage:
#   bash pipelines/snv/tests/run_tests.sh
# =============================================================================

set -uo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${TEST_DIR}/.." && pwd)"
SUITE_DIR="$(cd "${REPO_DIR}/.." && pwd)"
FIXTURES_DIR="${TEST_DIR}/fixtures"
RUN_DIR="$(mktemp -d)"

echo "=================================================="
echo "  SNV PIPELINE REGRESSION TEST SUITE"
echo "=================================================="
echo "Repo dir:     ${REPO_DIR}"
echo "Fixtures dir: ${FIXTURES_DIR}"
echo "Scratch dir:  ${RUN_DIR}"
echo ""

mkdir -p "${RUN_DIR}/results" "${RUN_DIR}/work" "${RUN_DIR}/logs"

# Ensure test fixtures exist
if [[ ! -f "${FIXTURES_DIR}/test_sample.wf_snp.vcf.gz" ]]; then
    echo "Fixtures missing. Building now via make_fixtures.sh..."
    bash "${TEST_DIR}/make_fixtures.sh"
fi

# Build config for scratch run
TEST_CONFIG="${RUN_DIR}/pipeline_config.yaml"
sed \
  -e "s|input_dir:.*|input_dir: \"${FIXTURES_DIR}\"|" \
  -e "s|output_dir:.*|output_dir: \"${RUN_DIR}/results\"|" \
  -e "s|work_dir:.*|work_dir: \"${RUN_DIR}/work\"|" \
  -e "s|repo_dir:.*|repo_dir: \"${REPO_DIR}\"|" \
  -e "s|config_file:.*|config_file: \"${FIXTURES_DIR}/test_reference_paths.yaml\"|" \
  "${FIXTURES_DIR}/test_pipeline_config.yaml" > "${TEST_CONFIG}"

sed -i.bak "s|log_dir:.*|log_dir: \"${RUN_DIR}/logs\"|" "${TEST_CONFIG}"
rm -f "${TEST_CONFIG}.bak"

# Run pipeline orchestrator
if ! bash "${REPO_DIR}/scripts/bash/04_run_all.sh" "${TEST_CONFIG}"; then
    echo "=== [run_tests.sh] Pipeline execution FAILED ==="
    exit 1
fi

PASS=0
FAIL=0
RESULTS_DIR="${RUN_DIR}/results"

assert_count() {
    local label="$1"
    local file="$2"
    local expected="$3"

    if [[ ! -f "${file}" ]]; then
        echo "  [FAIL] ${label} — file not found: ${file}"
        FAIL=$((FAIL+1))
        return
    fi
    local actual
    actual=$(( $(wc -l < "${file}" | tr -d ' ') - 1 ))
    if [[ "${actual}" -eq "${expected}" ]]; then
        echo "  [PASS] ${label}: expected ${expected}, got ${actual}"
        PASS=$((PASS+1))
    else
        echo "  [FAIL] ${label}: expected ${expected}, got ${actual} — file: ${file}"
        FAIL=$((FAIL+1))
    fi
}

assert_file_exists() {
    local label="$1"
    local file="$2"

    if [[ -f "${file}" ]]; then
        local size
        size=$(wc -c < "${file}" | tr -d ' ')
        if [[ "${size}" -gt 0 ]]; then
            echo "  [PASS] ${label}: exists and non-empty (${size} bytes)"
            PASS=$((PASS+1))
        else
            echo "  [FAIL] ${label}: file is 0 bytes: ${file}"
            FAIL=$((FAIL+1))
        fi
    else
        echo "  [FAIL] ${label}: missing file: ${file}"
        FAIL=$((FAIL+1))
    fi
}

echo ""
echo "=== Assertions ==="
assert_count "snvs table count"               "${RESULTS_DIR}/snvs/TEST_SAMPLE.snvs.tsv" 3
assert_count "indels table count"             "${RESULTS_DIR}/indels/TEST_SAMPLE.indels.tsv" 2
assert_count "clinvar pathogenic count"      "${RESULTS_DIR}/clinvar/TEST_SAMPLE.clinvar_pathogenic.tsv" 1
assert_count "annotated SNVs count"           "${RESULTS_DIR}/gene_snvs/TEST_SAMPLE.annotated_snvs.tsv" 5
assert_count "combined variants count"        "${RESULTS_DIR}/qc_summary/TEST_SAMPLE.all_snvs_combined.tsv" 5

# Verify figures
assert_file_exists "figure 1 (snv indel counts)"         "${RESULTS_DIR}/qc_summary/figures/TEST_SAMPLE.01_snv_indel_counts.png"
assert_file_exists "figure 2 (vaf depth distribution)"   "${RESULTS_DIR}/qc_summary/figures/TEST_SAMPLE.02_vaf_and_depth_distribution.png"
assert_file_exists "figure 3 (titv mutation spectrum)"  "${RESULTS_DIR}/qc_summary/figures/TEST_SAMPLE.03_titv_and_mutation_spectrum.png"
assert_file_exists "figure 4 (chromosome snv density)"  "${RESULTS_DIR}/qc_summary/figures/TEST_SAMPLE.04_chromosome_snv_density.png"

# Check summary statistics values
QC_SUMMARY="${RESULTS_DIR}/qc_summary/TEST_SAMPLE.filtering_summary.tsv"
STATS_FILE="${RESULTS_DIR}/qc_summary/TEST_SAMPLE.summary_statistics.tsv"
assert_file_exists "QC filtering summary TSV" "${QC_SUMMARY}"
assert_file_exists "Summary statistics TSV"   "${STATS_FILE}"

RAW_RECS=$(grep "^total_raw_records" "${QC_SUMMARY}" | cut -f2 || echo "0")
PASS_RECS=$(grep "^passed_filter_records" "${QC_SUMMARY}" | cut -f2 || echo "0")
TITV_RATIO=$(grep "^titv_ratio" "${STATS_FILE}" | cut -f2 || echo "0")

if [[ "${RAW_RECS}" -eq 6 ]]; then
    echo "  [PASS] QC summary: total raw records == 6"
    PASS=$((PASS+1))
else
    echo "  [FAIL] QC summary: expected 6 raw records, got ${RAW_RECS}"
    FAIL=$((FAIL+1))
fi

if [[ "${PASS_RECS}" -eq 5 ]]; then
    echo "  [PASS] QC summary: passed filter records == 5"
    PASS=$((PASS+1))
else
    echo "  [FAIL] QC summary: expected 5 passed records, got ${PASS_RECS}"
    FAIL=$((FAIL+1))
fi

if [[ "${TITV_RATIO}" == "2" || "${TITV_RATIO}" == "2.000" || "${TITV_RATIO}" == "2.0" ]]; then
    echo "  [PASS] Summary stats: Ti/Tv ratio == 2.0 (2 transitions / 1 transversion)"
    PASS=$((PASS+1))
else
    echo "  [FAIL] Summary stats: expected Ti/Tv ratio 2.0, got ${TITV_RATIO}"
    FAIL=$((FAIL+1))
fi

# Cleanup scratch directory
rm -rf "${RUN_DIR}"

echo ""
echo "=================================================="
echo "  SUMMARY OF SNV TEST RESULTS"
echo "=================================================="
echo "  Total tests: $((PASS + FAIL))"
echo "  Passed:      ${PASS}"
echo "  Failed:      ${FAIL}"
echo "=================================================="

if [[ "${FAIL}" -eq 0 ]]; then
    echo "  ALL SNV TESTS PASSED!"
    exit 0
else
    echo "  SOME SNV TESTS FAILED!"
    exit 1
fi
