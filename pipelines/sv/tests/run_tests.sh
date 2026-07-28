#!/usr/bin/env bash
# =============================================================================
# run_tests.sh
# Regression test: runs the full pipeline (04_run_all.sh) against a small,
# hand-crafted synthetic VCF (tests/fixtures/test_sample.wf_sv.vcf.gz) with
# known expected outcomes, then asserts the actual output counts match.
#
# The fixture includes at least one case per SV category, a sub-QC-threshold
# record per filter (to confirm filtering actually excludes it), and the
# real BND/ANN records validated by hand during pipeline development
# (transcript_ablation, feature_fusion with the gene pair in the SECOND
# annotation block, bidirectional_gene_fusion appearing twice as a
# duplicate-event case).
#
# Usage:
#   bash tests/run_tests.sh
#
# Exit code 0 = all assertions passed. Non-zero = at least one failed
# (see printed FAIL lines for which).
# =============================================================================

set -uo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${TEST_DIR}/.." && pwd)"
FIXTURES_DIR="${TEST_DIR}/fixtures"
RUN_DIR="$(mktemp -d)"

echo "=== ont-sv-pipeline regression test ==="
echo "Repo dir:      ${REPO_DIR}"
echo "Fixtures dir:  ${FIXTURES_DIR}"
echo "Scratch dir:   ${RUN_DIR}"
echo ""

mkdir -p "${RUN_DIR}/results" "${RUN_DIR}/work" "${RUN_DIR}/logs"

# --- Generate a real config for this test run, with real absolute paths ----
TEST_CONFIG="${RUN_DIR}/pipeline_config.yaml"
sed \
  -e "s|input_dir:.*|input_dir: \"${FIXTURES_DIR}\"|" \
  -e "s|output_dir:.*|output_dir: \"${RUN_DIR}/results\"|" \
  -e "s|work_dir:.*|work_dir: \"${RUN_DIR}/work\"|" \
  -e "s|repo_dir:.*|repo_dir: \"${REPO_DIR}\"|" \
  -e "s|config_file:.*|config_file: \"${FIXTURES_DIR}/test_reference_paths.yaml\"|" \
  "${FIXTURES_DIR}/test_pipeline_config.yaml" > "${TEST_CONFIG}"

# logs/ referenced relative to repo_dir by 00_setup_env.sh; point it at scratch
# -i.bak, not bare -i: BSD sed on macOS requires an argument to -i and would
# otherwise consume the script expression as the backup suffix, then fail. GNU
# sed accepts both. The .bak file is removed immediately.
sed -i.bak "s|log_dir:.*|log_dir: \"${RUN_DIR}/logs\"|" "${TEST_CONFIG}"
rm -f "${TEST_CONFIG}.bak"

echo "Generated test config: ${TEST_CONFIG}"
echo ""

# --- Run the full pipeline against the fixture ------------------------------
if ! bash "${REPO_DIR}/scripts/bash/04_run_all.sh" "${TEST_CONFIG}"; then
    echo ""
    echo "=== [run_tests.sh] Pipeline itself failed to complete — see output above ==="
    exit 1
fi

# --- Assertions --------------------------------------------------------------
PASS=0
FAIL=0
RESULTS_DIR="${RUN_DIR}/results"

assert_count() {
    local label="$1"
    local file="$2"
    local expected="$3"
    if [[ ! -f "${file}" ]]; then
        echo "FAIL: ${label} — file not found: ${file}"
        FAIL=$((FAIL+1))
        return
    fi
    local actual=$(( $(wc -l < "${file}") - 1 ))  # minus header
    if [[ "${actual}" -eq "${expected}" ]]; then
        echo "PASS: ${label} (expected ${expected}, got ${actual})"
        PASS=$((PASS+1))
    else
        echo "FAIL: ${label} (expected ${expected}, got ${actual}) — ${file}"
        FAIL=$((FAIL+1))
    fi
}

echo ""
echo "=== Assertions ==="
assert_count "deletions"              "${RESULTS_DIR}/deletions/TEST_SAMPLE.deletions.tsv"                       2
assert_count "insertions"             "${RESULTS_DIR}/insertions/TEST_SAMPLE.insertions.tsv"                     1
assert_count "duplications"           "${RESULTS_DIR}/duplications/TEST_SAMPLE.duplications.tsv"                 1
assert_count "inversions"             "${RESULTS_DIR}/inversions/TEST_SAMPLE.inversions.tsv"                     1
assert_count "translocations"         "${RESULTS_DIR}/translocations/TEST_SAMPLE.translocations.tsv"             4
assert_count "complex_rearrangements" "${RESULTS_DIR}/complex_rearrangements/TEST_SAMPLE.complex_rearrangements.tsv" 0
assert_count "gene_fusions"           "${RESULTS_DIR}/gene_fusions/TEST_SAMPLE.gene_fusions.tsv"                 3
assert_count "gene_disruptions"       "${RESULTS_DIR}/gene_fusions/TEST_SAMPLE.gene_disruptions.tsv"             1

# Content-level checks, not just counts — confirm the ANN multi-block parsing
# edge case (gene pair in the SECOND annotation block) still resolves correctly.
echo ""
if grep -q "ZSWIM5" "${RESULTS_DIR}/gene_fusions/TEST_SAMPLE.gene_fusions.tsv" 2>/dev/null; then
    echo "PASS: multi-block ANN parsing recovers ZSWIM5 gene pair from second block"
    PASS=$((PASS+1))
else
    echo "FAIL: ZSWIM5 not found in gene_fusions.tsv — multi-block ANN parsing regressed"
    FAIL=$((FAIL+1))
fi

n_dup_flagged=$(awk -F'\t' 'NR>1 && $NF=="TRUE"' "${RESULTS_DIR}/gene_fusions/TEST_SAMPLE.gene_fusions.tsv" 2>/dev/null | wc -l)
if [[ "${n_dup_flagged}" -eq 2 ]]; then
    echo "PASS: possible_duplicate_event correctly flags 2 records (the LINC00707/LOC105378800 pair called twice)"
    PASS=$((PASS+1))
else
    echo "FAIL: expected 2 records flagged possible_duplicate_event, got ${n_dup_flagged}"
    FAIL=$((FAIL+1))
fi

# --- Regression check: archiving must report THIS sample's counts, not a --
# neighboring sample's, when results/ is shared across multiple samples
# (a real bug found in development: get_count() previously grabbed
# whichever *.tsv file the filesystem listed first in a shared directory,
# rather than the sample-specific filename).
echo ""
echo "=== Archiving regression check (shared results/ directory across samples) ==="

DECOY_HEADER="CHROM	POS	ID	QUAL	FILTER	SVTYPE	SVLEN	END	CHR2	SUPPORT	STRAND	VAF	GT	DR	DV	GENE_SYMBOL	ANNOTATION_IMPACT"
{ echo "${DECOY_HEADER}"; for i in 1 2 3 4 5 6 7 8 9; do echo -e "chr9\t${i}000\tdecoy${i}\t50\tPASS\tDEL\t-999\t1999\t.\t10\t+-\t0.5\t0/1\t10\t10\tNA\tMODIFIER"; done; } \
    > "${RESULTS_DIR}/deletions/DECOY_SAMPLE.deletions.tsv"

sed -e "s|archive_root:.*|archive_root: \"${RUN_DIR}/archive\"|" -e "s|compress:.*|compress: false|" "${TEST_CONFIG}" > "${TEST_CONFIG}.archive_test" 2>/dev/null
if ! grep -q "^archive:" "${TEST_CONFIG}"; then
    printf '\narchive:\n  archive_root: "%s/archive"\n  compress: false\n' "${RUN_DIR}" >> "${TEST_CONFIG}"
fi
mkdir -p "${RUN_DIR}/archive"

bash "${REPO_DIR}/scripts/bash/08_archive_results.sh" "${TEST_CONFIG}" > /dev/null 2>&1

ARCHIVE_INDEX="${RUN_DIR}/archive/archive_index.tsv"
if [[ -f "${ARCHIVE_INDEX}" ]]; then
    n_del_archived=$(awk -F'\t' 'NR>1 && $2=="TEST_SAMPLE" {print $5}' "${ARCHIVE_INDEX}" | tail -1)
    if [[ "${n_del_archived}" -eq 2 ]]; then
        echo "PASS: archived deletion count for TEST_SAMPLE is correct (2) despite a 9-record DECOY_SAMPLE sharing the same results/deletions/ directory"
        PASS=$((PASS+1))
    else
        echo "FAIL: archived deletion count is ${n_del_archived}, expected 2 (decoy sample's count may have leaked in — check get_count() in 08_archive_results.sh)"
        FAIL=$((FAIL+1))
    fi

    n_decoy_files_in_archive=$(find "${RUN_DIR}/archive" -type f -name "DECOY_SAMPLE.*" 2>/dev/null | wc -l)
    if [[ "${n_decoy_files_in_archive}" -eq 0 ]]; then
        echo "PASS: TEST_SAMPLE's archive contains zero DECOY_SAMPLE files (no cross-sample contamination)"
        PASS=$((PASS+1))
    else
        echo "FAIL: found ${n_decoy_files_in_archive} DECOY_SAMPLE file(s) inside TEST_SAMPLE's archive — cross-sample contamination regressed"
        FAIL=$((FAIL+1))
    fi
else
    echo "FAIL: archive_index.tsv was not created — archiving step itself failed"
    FAIL=$((FAIL+1))
fi

echo ""
echo "=== Test run complete: ${PASS} passed, ${FAIL} failed ==="
echo "(Scratch dir left at ${RUN_DIR} for inspection; safe to delete.)"

if [[ ${FAIL} -gt 0 ]]; then
    exit 1
fi
exit 0
