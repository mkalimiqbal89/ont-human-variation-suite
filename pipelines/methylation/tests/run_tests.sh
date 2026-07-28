#!/usr/bin/env bash
# =============================================================================
# tests/run_tests.sh
# Regression suite for the ONT methylation pipeline.
#
# Assertions are CONTENT-LEVEL, not just exit codes. Exit codes alone miss the
# bugs that matter: a stage can succeed while silently dropping records,
# double-counting, or writing a truncated output that looks plausible.
#
# All scratch output goes to a mktemp directory OUTSIDE the repo. Nothing here
# may write into results/ or logs/ — a test run must never be mistakable for a
# real run, and must never leave artefacts behind in the repo.
#
# Usage:
#   bash tests/run_tests.sh
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
FIXTURE_DIR="${SCRIPT_DIR}/fixtures"
PREFIX="TEST_SAMPLE"

SCRATCH="$(mktemp -d "${TMPDIR:-/tmp}/ont-methyl-tests.XXXXXX")"
trap 'rm -rf "${SCRATCH}"' EXIT

PASS=0
FAIL=0

ok()   { PASS=$((PASS+1)); printf '    [PASS] %s\n' "$1"; }
bad()  { FAIL=$((FAIL+1)); printf '    [FAIL] %s\n' "$1"; }

# assert_eq <label> <expected> <actual>
assert_eq() {
    if [[ "$2" == "$3" ]]; then ok "$1 ($3)"; else bad "$1: expected '$2', got '$3'"; fi
}

echo "============================================================"
echo " ONT Methylation Pipeline — regression suite"
echo "============================================================"
echo "Repo    : ${REPO_DIR}"
echo "Scratch : ${SCRATCH}"
echo ""

# --- Fixtures ----------------------------------------------------------------
echo "--- Building fixtures ---"
bash "${SCRIPT_DIR}/make_fixtures.sh" "${FIXTURE_DIR}" >/dev/null || {
    echo "[FATAL] fixture generation failed"; exit 1; }
echo "    fixtures ready"
echo ""

# Writes a config for one fixture variant into the scratch dir.
# make_config <variant> [extra_sed_expr]
#
# Every call gets its OWN directory. Keying it on the variant name alone let one
# test see the previous test's outputs — which is how "stage 03 must refuse to
# run without stage 02's output" passed spuriously.
#
# mktemp rather than an incrementing counter: this function is always called as
# DIR="$(make_config ...)", i.e. inside a subshell, so any counter variable it
# increments is discarded when that subshell exits. mktemp cannot be defeated
# that way.
make_config() {
    local variant="$1" extra="${2:-}"
    local dir
    dir="$(mktemp -d "${SCRATCH}/${variant}.XXXXXX")"
    mkdir -p "${dir}/input" "${dir}/results" "${dir}/work" "${dir}/logs"
    cp "${FIXTURE_DIR}/${PREFIX}.${variant}.wf_mods.bedmethyl.gz" \
       "${dir}/input/${PREFIX}.wf_mods.bedmethyl.gz"
    # Derived from config/pipeline_config.example.yaml rather than hand-written
    # here. A heredoc copy drifts: every config key added to the pipeline was
    # missing from the test config until something failed, and a sed override
    # targeting an absent key silently does nothing, so tests passed against
    # defaults while appearing to test the setting. Deriving from the template
    # means a new key is present in tests the moment it is added.
    sed -e "s|^  sample_id: .*|  sample_id: \"TEST_01\"|" \
        -e "s|^  raw_sample_prefix: .*|  raw_sample_prefix: \"${PREFIX}\"|" \
        -e "s|^  run_description: .*|  run_description: \"synthetic fixture: ${variant}\"|" \
        -e "s|^  input_dir: .*|  input_dir: \"${dir}/input\"|" \
        -e "s|^  output_dir: .*|  output_dir: \"${dir}/results\"|" \
        -e "s|^  work_dir: .*|  work_dir: \"${dir}/work\"|" \
        -e "s|^  repo_dir: .*|  repo_dir: \"${REPO_DIR}\"|" \
        -e "s|^  config_file: .*|  config_file: \"${FIXTURE_DIR}/test_reference_paths.yaml\"|" \
        -e "s|^  log_dir: .*|  log_dir: \"${dir}/logs\"|" \
        -e "s|^  archive_root: .*|  archive_root: \"${dir}/archive\"|" \
        -e "s|^  threads: .*|  threads: 2|" \
        -e "s|^  expected_primary_contigs: .*|  expected_primary_contigs: 1|" \
        "${REPO_DIR}/config/pipeline_config.example.yaml" > "${dir}/config.yaml"
    [[ -n "${extra}" ]] && sed -i.bak "${extra}" "${dir}/config.yaml" && rm -f "${dir}/config.yaml.bak"
    echo "${dir}"
}

qc_val() {  # qc_val <results_dir> <metric>   — stage 02 QC
    awk -F'\t' -v k="$2" '$1==k {print $2; exit}' \
        "$1/02_qc/TEST_01.02_extraction_qc.tsv" 2>/dev/null
}

qc3_val() { # qc3_val <results_dir> <metric>  — stage 03 QC
    awk -F'\t' -v k="$2" '$1==k {print $2; exit}' \
        "$1/02_qc/TEST_01.03_filter_qc.tsv" 2>/dev/null
}

# Stage 03 consumes stage 02's output, so the chain must be run in order.
run_02_then_03() {
    bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "$1/config.yaml" >/dev/null 2>&1 || return 1
    bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "$1/config.yaml" >/dev/null 2>&1
}

# =============================================================================
# TEST GROUP 1 — stage 00: config integrity
# =============================================================================
echo "--- [1] stage 00: config integrity ---"

DIR="$(make_config valid)"
bash "${REPO_DIR}/scripts/bash/00_setup_env.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "clean config sources OK" "0" "$?"

# A duplicate key must be a hard error: yaml_get() returns the first match, so
# the second value would be silently ignored.
cp "${DIR}/config.yaml" "${SCRATCH}/dup.yaml"
printf '\nextra_section:\n  min_coverage: 99\n' >> "${SCRATCH}/dup.yaml"
bash "${REPO_DIR}/scripts/bash/00_setup_env.sh" "${SCRATCH}/dup.yaml" >/dev/null 2>&1
assert_eq "duplicate config key rejected" "1" "$?"

bash "${REPO_DIR}/scripts/bash/00_setup_env.sh" "${SCRATCH}/does_not_exist.yaml" >/dev/null 2>&1
assert_eq "missing config rejected" "1" "$?"

# LOG_DIR must NOT be rewritten into the repo when the config gives an absolute
# path — otherwise test logs leak into the real repo's logs/ directory.
LOGDIR_SEEN="$(bash -c "source '${REPO_DIR}/scripts/bash/00_setup_env.sh' '${DIR}/config.yaml' >/dev/null 2>&1; echo \$LOG_DIR")"
assert_eq "absolute log_dir not prefixed with repo" "${DIR}/logs" "${LOGDIR_SEEN}"
echo ""

# =============================================================================
# TEST GROUP 2 — stage 02: structural assertions must FAIL loudly
# =============================================================================
echo "--- [2] stage 02: malformed input must fail ---"

# variant -> expected assertion counter that should be non-zero
declare -a NEG=(
    "bad_columns:ASSERT_wrong_column_count"
    "unsorted:ASSERT_unsorted_records"
    "interleaved_contigs:ASSERT_interleaved_contigs"
    "percent_out_of_range:ASSERT_percent_out_of_range"
)

for entry in "${NEG[@]}"; do
    variant="${entry%%:*}"; counter="${entry#*:}"
    DIR="$(make_config "${variant}")"
    bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
    rc=$?
    [[ ${rc} -eq 1 ]] && ok "${variant}: exits 1" || bad "${variant}: expected exit 1, got ${rc}"

    n="$(qc_val "${DIR}/results" "${counter}")"
    if [[ "${n:-0}" -gt 0 ]]; then ok "${variant}: ${counter}=${n}"
    else bad "${variant}: ${counter} was ${n:-unset}, expected > 0"; fi

    # A truncated output that looks plausible is worse than none: downstream
    # stages would analyse an incomplete genome without complaining.
    OUTF="${DIR}/results/01_filtered/TEST_01.5mC.all_cov.bedmethyl.gz"
    [[ -f "${OUTF}" ]] && bad "${variant}: partial output left behind" \
                       || ok "${variant}: partial output removed"
done

# An empty result must be an error, not an empty file.
DIR="$(make_config no_primary_code)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "no_primary_code: exits 1" "1" "$?"
echo ""

# =============================================================================
# TEST GROUP 3 — stage 02: content correctness on valid input
# =============================================================================
echo "--- [3] stage 02: content correctness ---"

DIR="$(make_config valid)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "valid: exits 0" "0" "$?"

R3="${DIR}/results"
assert_eq "valid: total records read"        "6" "$(qc_val "${R3}" Total_records_read)"
assert_eq "valid: both mod codes seen"       "2" "$(qc_val "${R3}" Distinct_modification_codes)"
# The h records must be counted but NOT retained. Getting this wrong doubles
# every site count downstream.
assert_eq "valid: h records counted"         "3" "$(qc_val "${R3}" Records_mod_code_h)"
assert_eq "valid: only m retained"           "3" "$(qc_val "${R3}" Retained_CpG_sites)"

OUTF="${R3}/01_filtered/TEST_01.5mC.all_cov.bedmethyl.gz"
assert_eq "valid: output row count"          "3" "$(gzip -cd "${OUTF}" | wc -l | tr -d ' ')"
assert_eq "valid: output has 18 columns"     "18" "$(gzip -cd "${OUTF}" | head -n1 | awk -F'\t' '{print NF}')"
assert_eq "valid: no h leaked into output"   "0" "$(gzip -cd "${OUTF}" | awk -F'\t' '$4!="m"' | wc -l | tr -d ' ')"
[[ -f "${OUTF}.tbi" ]] && ok "valid: tabix index written" || bad "valid: tabix index missing"

# Contig exclusion: the four junk contigs must go, the three chr1 sites stay.
DIR="$(make_config excluded_contigs)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "excluded_contigs: exits 0" "0" "$?"
R3="${DIR}/results"
assert_eq "excluded_contigs: 4 records dropped" "4" "$(qc_val "${R3}" Excluded_by_contig_filter)"
assert_eq "excluded_contigs: 4 contigs dropped" "4" "$(qc_val "${R3}" Excluded_contigs)"
assert_eq "excluded_contigs: 3 sites retained" "3" "$(qc_val "${R3}" Retained_CpG_sites)"
assert_eq "excluded_contigs: 1 contig retained" "1" "$(qc_val "${R3}" Contigs_retained)"
OUTF="${R3}/01_filtered/TEST_01.5mC.all_cov.bedmethyl.gz"
assert_eq "excluded_contigs: no chrM in output" "0" \
          "$(gzip -cd "${OUTF}" | awk -F'\t' '$1=="chrM"' | wc -l | tr -d ' ')"

# =============================================================================
# TEST GROUP 3b — contig selection: allowlist vs blocklist
# =============================================================================
echo ""
echo "--- [3b] contig selection ---"

# The allowlist must reject scaffolds the blocklist regex never matches.
# Under the blocklist, GL000191.1 / HLA-A* / *_fix all survive; under the
# allowlist none of them do. This is the reason the allowlist is the default.
DIR="$(make_config nonstandard_contigs)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "allowlist: exits 0" "0" "$?"
R3B="${DIR}/results"
assert_eq "allowlist: mode recorded"          "allowlist" "$(qc_val "${R3B}" Contig_selection_mode)"
assert_eq "allowlist: 4 scaffolds dropped"    "4" "$(qc_val "${R3B}" Excluded_contigs)"
assert_eq "allowlist: only chr1 retained"     "1" "$(qc_val "${R3B}" Contigs_retained)"
assert_eq "allowlist: 3 sites retained"       "3" "$(qc_val "${R3B}" Retained_CpG_sites)"
OUTF="${R3B}/01_filtered/TEST_01.5mC.all_cov.bedmethyl.gz"
assert_eq "allowlist: no non-chr contigs leaked" "0" \
          "$(gzip -cd "${OUTF}" | awk -F'\t' '$1!="chr1"' | wc -l | tr -d ' ')"
[[ -s "${R3B}/02_qc/TEST_01.02_excluded_contigs.tsv" ]] \
    && ok "allowlist: excluded-contig QC file written" \
    || bad "allowlist: excluded-contig QC file missing"

# Same fixture under the BLOCKLIST: three of the four scaffolds now survive,
# which is precisely the failure mode the allowlist prevents. Asserting the bad
# behaviour explicitly means the comparison stays honest if the regex changes.
DIR="$(make_config nonstandard_contigs "s|^  primary_contigs_only: true|  primary_contigs_only: false|")"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
R3C="${DIR}/results"
assert_eq "blocklist: mode recorded"            "blocklist" "$(qc_val "${R3C}" Contig_selection_mode)"
assert_eq "blocklist: only 1 scaffold dropped"  "1" "$(qc_val "${R3C}" Excluded_contigs)"
assert_eq "blocklist: 4 contigs leak through"   "4" "$(qc_val "${R3C}" Contigs_retained)"

# Missing chrY must warn, not fail — the female-sample case.
DIR="$(make_config no_chry "s|^  expected_primary_contigs: 1|  expected_primary_contigs: 3|")"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "unexpected contig count: still exits 0" "0" "$?"
assert_eq "unexpected contig count: warned"        "1" \
          "$(qc_val "${DIR}/results" WARN_unexpected_contig_count)"
DIR="$(make_config no_chry "s|^  expected_primary_contigs: 1|  expected_primary_contigs: 2|")"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "matching contig count: no warning"      "0" \
          "$(qc_val "${DIR}/results" WARN_unexpected_contig_count)"
echo ""

# counts_exceed_coverage is a WARNING counter, not fatal: modkit can legitimately
# produce edge cases here, so it is surfaced rather than used to abort.
DIR="$(make_config counts_exceed_coverage)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "counts_exceed_coverage: still exits 0" "0" "$?"
assert_eq "counts_exceed_coverage: flagged" "1" \
          "$(qc_val "${DIR}/results" ASSERT_counts_exceed_coverage)"
echo ""

# =============================================================================
# TEST GROUP 4 — stage 02: arithmetic on known values
# =============================================================================
echo "--- [4] stage 02: arithmetic ---"

DIR="$(make_config known_values)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
R4="${DIR}/results"
# coverage 20/40/60, modified 10/20/30 -> weighted = 60/120 = 50.0000%
assert_eq "known: total valid coverage"   "120"      "$(qc_val "${R4}" Total_valid_coverage)"
assert_eq "known: total modified calls"   "60"       "$(qc_val "${R4}" Total_modified_read_calls)"
assert_eq "known: weighted methylation"   "50.0000"  "$(qc_val "${R4}" Coverage_weighted_methylation_percent)"
assert_eq "known: unweighted methylation" "50.0000"  "$(qc_val "${R4}" Unweighted_mean_site_methylation_percent)"
assert_eq "known: mean coverage"          "40.0000"  "$(qc_val "${R4}" Mean_coverage)"
assert_eq "known: min coverage"           "20"       "$(qc_val "${R4}" Minimum_coverage)"
assert_eq "known: max coverage"           "60"       "$(qc_val "${R4}" Maximum_coverage)"

# The coverage histogram must account for every retained site, exactly once.
HIST="${R4}/02_qc/TEST_01.02_coverage_histogram.tsv"
assert_eq "known: histogram sums to site count" "3" \
          "$(awk -F'\t' 'NR>1 {s+=$2} END {print s+0}' "${HIST}")"
echo ""

# =============================================================================
# TEST GROUP 5 — stage 03: coverage filter and methylation states
# =============================================================================
echo "--- [5] stage 03: filtering and state assignment ---"

DIR="$(make_config boundaries)"
run_02_then_03 "${DIR}"
assert_eq "boundaries: exits 0" "0" "$?"
R5="${DIR}/results"

# Coverage bound is inclusive: cov 10 stays, cov 9 goes.
assert_eq "boundaries: input sites"      "8" "$(qc3_val "${R5}" Input_sites)"
assert_eq "boundaries: dropped low cov"  "1" "$(qc3_val "${R5}" Dropped_low_coverage)"
assert_eq "boundaries: retained"         "7" "$(qc3_val "${R5}" Retained_sites)"

# State boundaries: unmethylated is strictly < 20, methylated is >= 80. A site at
# exactly 20.00 belongs to intermediate; one at exactly 80.00 to methylated.
assert_eq "boundaries: unmethylated (0, 19.99)"          "2" "$(qc3_val "${R5}" State_unmethylated_sites)"
assert_eq "boundaries: intermediate (20.00, 50, 79.99)"  "3" "$(qc3_val "${R5}" State_intermediate_sites)"
assert_eq "boundaries: methylated (80.00, 100)"          "2" "$(qc3_val "${R5}" State_methylated_sites)"

# The two reconciliations. Neither can hold by accident.
assert_eq "boundaries: filter reconciliation"   "0" "$(qc3_val "${R5}" ASSERT_filter_reconciliation)"
assert_eq "boundaries: category reconciliation" "0" "$(qc3_val "${R5}" ASSERT_category_reconciliation)"

# Same reconciliation, but against the files on disk rather than awk's counters.
# A counter can agree with itself while the file writing is wrong.
F5="${R5}/01_filtered"
S5="${F5}/by_methylation_state"
MAIN_N=$(gzip -cd "${F5}/TEST_01.5mC.cov10.bedmethyl.gz" | wc -l | tr -d ' ')
SUM_N=0
for st in unmethylated intermediate methylated; do
    f="${S5}/TEST_01.5mC.cov10.${st}.bedmethyl.gz"
    n=$([[ -s "${f}" ]] && gzip -cd "${f}" | wc -l | tr -d ' ' || echo 0)
    SUM_N=$((SUM_N + n))
done
assert_eq "boundaries: main file row count"        "7" "${MAIN_N}"
assert_eq "boundaries: state files sum to main"    "${MAIN_N}" "${SUM_N}"

# Every record in a state file must actually satisfy that state's condition.
BAD_U=$(gzip -cd "${S5}/TEST_01.5mC.cov10.unmethylated.bedmethyl.gz" | awk -F'\t' '$11+0 >= 20' | wc -l | tr -d ' ')
BAD_M=$(gzip -cd "${S5}/TEST_01.5mC.cov10.methylated.bedmethyl.gz"   | awk -F'\t' '$11+0 <  80' | wc -l | tr -d ' ')
BAD_I=$(gzip -cd "${S5}/TEST_01.5mC.cov10.intermediate.bedmethyl.gz" | awk -F'\t' '$11+0 < 20 || $11+0 >= 80' | wc -l | tr -d ' ')
assert_eq "boundaries: no misfiled unmethylated" "0" "${BAD_U}"
assert_eq "boundaries: no misfiled methylated"   "0" "${BAD_M}"
assert_eq "boundaries: no misfiled intermediate" "0" "${BAD_I}"

# No site below the coverage threshold may survive into the main output.
assert_eq "boundaries: no sub-threshold coverage in output" "0" \
          "$(gzip -cd "${F5}/TEST_01.5mC.cov10.bedmethyl.gz" | awk -F'\t' '$10+0 < 10' | wc -l | tr -d ' ')"

# Upper coverage cap. With max_coverage=15 the six cov-20 sites are dropped too,
# leaving only the cov-10 site.
DIR="$(make_config boundaries "s|^  max_coverage: 0|  max_coverage: 15|")"
run_02_then_03 "${DIR}"
assert_eq "max_coverage: exits 0" "0" "$?"
R5B="${DIR}/results"
assert_eq "max_coverage: dropped high cov" "6" "$(qc3_val "${R5B}" Dropped_high_coverage)"
assert_eq "max_coverage: dropped low cov"  "1" "$(qc3_val "${R5B}" Dropped_low_coverage)"
assert_eq "max_coverage: retained"         "1" "$(qc3_val "${R5B}" Retained_sites)"
assert_eq "max_coverage: reconciliation"   "0" "$(qc3_val "${R5B}" ASSERT_filter_reconciliation)"

# Threshold above every site's coverage must fail, not write an empty file.
DIR="$(make_config boundaries "s|^  min_coverage: 10|  min_coverage: 9999|")"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "impossible threshold: exits 1" "1" "$?"
[[ -f "${DIR}/results/01_filtered/TEST_01.5mC.cov9999.bedmethyl.gz" ]] \
    && bad "impossible threshold: partial output left behind" \
    || ok "impossible threshold: partial output removed"

# Stage 03 must refuse to invent its input if stage 02 has not run.
DIR="$(make_config boundaries)"
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "missing stage 02 output: exits 1" "1" "$?"
echo ""

# =============================================================================
# TEST GROUP 6 — stage 05: annotation and aggregation keys
# =============================================================================
echo "--- [6] stage 05: annotation ---"

if ! command -v bedtools >/dev/null 2>&1; then
    echo "    [SKIP] bedtools not on PATH — stage 05 tests skipped"
else
qc5_val() { # qc5_val <results_dir> <metric>
    awk -F'\t' -v k="$2" '$1==k {print $2; exit}' \
        "$1/02_qc/TEST_01.05_annotation_qc.tsv" 2>/dev/null
}

DIR="$(make_config annotation)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/05_annotate_regions.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "annotation: exits 0" "0" "$?"
R6="${DIR}/results"

# chr1:1000 overlaps two genes, chr1:5100 one, chr1:20000 none.
assert_eq "gene: intersect rows"        "3" "$(qc5_val "${R6}" gene_intersect_rows)"
assert_eq "gene: unique CpGs annotated" "2" "$(qc5_val "${R6}" gene_unique_cpgs_annotated)"

# THE key assertion. ENSGA.1 and ENSGB.1 share the gene_name "SHARED". Keyed on
# gene_id this is 3 features; keyed on gene_name it would collapse to 2. In the
# real annotation gene_name repeats 1,615 times, Y_RNA alone 756 times.
assert_eq "gene: gene_name collision not merged (3 features)" "3" \
          "$(qc5_val "${R6}" gene_features_with_data)"

# THE other key assertion. Both islands are named "CpG:_21" at different loci.
# Keyed on coordinates this is 2 features; keyed on name it would be 1. The real
# annotation has 32,038 islands but only 485 distinct names.
assert_eq "cpg_island: name collision not merged (2 features)" "2" \
          "$(qc5_val "${R6}" cpg_island_features_with_data)"
assert_eq "cpg_island: intersect rows" "2" "$(qc5_val "${R6}" cpg_island_intersect_rows)"

# Reconciliations for all three feature types.
for l in gene promoter cpg_island; do
    assert_eq "${l}: cpg sum vs rows"      "0" "$(qc5_val "${R6}" "ASSERT_${l}_cpg_sum_vs_rows")"
    assert_eq "${l}: written vs features"  "0" "$(qc5_val "${R6}" "ASSERT_${l}_written_vs_features")"
    assert_eq "${l}: weighted in range"    "0" "$(qc5_val "${R6}" "ASSERT_${l}_weighted_out_of_range")"
done

# Output files: shape, header, and that both distinct islands really are present.
GS="${R6}/06_gene_summary/TEST_01.gene_methylation_summary.tsv"
CS="${R6}/08_cpg_island_summary/TEST_01.cpg_island_methylation_summary.tsv"
assert_eq "gene summary: data rows"     "3"  "$(awk -F'\t' 'NR>1' "${GS}" | wc -l | tr -d ' ')"
assert_eq "gene summary: 13 columns"    "13" "$(head -n1 "${GS}" | awk -F'\t' '{print NF}')"
assert_eq "island summary: data rows"   "2"  "$(awk -F'\t' 'NR>1' "${CS}" | wc -l | tr -d ' ')"
assert_eq "island summary: 15 columns"  "15" "$(head -n1 "${CS}" | awk -F'\t' '{print NF}')"
assert_eq "island summary: distinct island ids" "2" \
          "$(awk -F'\t' 'NR>1 {print $1}' "${CS}" | sort -u | wc -l | tr -d ' ')"

# Min-CpG flag: every feature here has 1-2 CpGs, below the threshold of 5.
assert_eq "gene: none meet min CpGs" "0" "$(qc5_val "${R6}" gene_features_min_cpgs_met)"
assert_eq "gene summary: all flagged no" "3" \
          "$(awk -F'\t' 'NR>1 && $13=="no"' "${GS}" | wc -l | tr -d ' ')"

# Per-feature arithmetic. chr1:1000 has coverage 20 and 15 modified calls, so
# both ENSGA.1 and ENSGB.1 must report 75.0000% weighted.
assert_eq "gene ENSGA.1 weighted methylation" "75.0000" \
          "$(awk -F'\t' '$1=="ENSGA.1" {print $12}' "${GS}")"
assert_eq "gene ENSGB.1 weighted methylation" "75.0000" \
          "$(awk -F'\t' '$1=="ENSGB.1" {print $12}' "${GS}")"
assert_eq "gene ENSGC.1 weighted methylation" "40.0000" \
          "$(awk -F'\t' '$1=="ENSGC.1" {print $12}' "${GS}")"

# Per-CpG intermediates suppressed by default; produced on request.
[[ -f "${R6}/05_annotated_cpgs/TEST_01.gene_annotated_CpGs.tsv.gz" ]] \
    && bad "per-CpG intermediate written despite keep_annotated_cpgs=false" \
    || ok "per-CpG intermediate suppressed by default"

DIR="$(make_config annotation "s|^  keep_annotated_cpgs: false|  keep_annotated_cpgs: true|")"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/05_annotate_regions.sh" "${DIR}/config.yaml" >/dev/null 2>&1
A6="${DIR}/results/05_annotated_cpgs/TEST_01.gene_annotated_CpGs.tsv.gz"
[[ -f "${A6}" ]] && ok "per-CpG intermediate written on request" \
                 || bad "per-CpG intermediate missing when requested"
if [[ -f "${A6}" ]]; then
    assert_eq "per-CpG intermediate: 25 columns" "25" \
              "$(gzip -cd "${A6}" | head -n1 | awk -F'\t' '{print NF}')"
    assert_eq "per-CpG intermediate: 3 rows" "3" \
              "$(gzip -cd "${A6}" | wc -l | tr -d ' ')"
fi

# Stage 05 must refuse to run before stage 03.
DIR="$(make_config annotation)"
bash "${REPO_DIR}/scripts/bash/05_annotate_regions.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "missing stage 03 output: exits 1" "1" "$?"
fi
echo ""

# =============================================================================
# TEST GROUP 7 — stage 06: summary statistics
# =============================================================================
echo "--- [7] stage 06: summary statistics ---"

qc6_val() { # qc6_val <results_dir> <metric>
    awk -F'\t' -v k="$2" '$1==k {print $2; exit}' \
        "$1/02_qc/TEST_01.06_summary_qc.tsv" 2>/dev/null
}
g6_val() {  # g6_val <results_dir> <metric>  — global methylation table
    awk -F'\t' -v k="$2" '$1==k {print $2; exit}' \
        "$1/03_global_methylation/TEST_01.global_methylation_summary.tsv" 2>/dev/null
}

# known_values: coverage 20/40/60, modified 10/20/30 -> weighted exactly 50%.
DIR="$(make_config known_values)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/06_summary_stats.sh"    "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "summary: exits 0" "0" "$?"
R7="${DIR}/results"

assert_eq "summary: retained sites"       "3"       "$(g6_val "${R7}" Retained_CpG_sites)"
assert_eq "summary: total coverage"       "120"     "$(g6_val "${R7}" Total_valid_coverage)"
assert_eq "summary: total modified"       "60"      "$(g6_val "${R7}" Total_modified_read_calls)"
assert_eq "summary: weighted methylation" "50.0000" "$(g6_val "${R7}" Coverage_weighted_methylation_percent)"
assert_eq "summary: mean coverage"        "40.0000" "$(g6_val "${R7}" Mean_coverage)"

# Both totals must account for every site, exactly once.
assert_eq "summary: distribution sums to sites" "0" "$(qc6_val "${R7}" ASSERT_distribution_sums_to_sites)"
assert_eq "summary: chromosomes sum to sites"   "0" "$(qc6_val "${R7}" ASSERT_chromosome_sums_to_sites)"
assert_eq "summary: agrees with stage 03"       "0" "$(qc6_val "${R7}" ASSERT_matches_stage03)"

# Distribution table shape: 10 bins, top bin labelled 90-100 because 100% folds
# into it rather than forming an eleventh bin.
D7="${R7}/03_global_methylation/TEST_01.methylation_distribution.tsv"
assert_eq "distribution: 10 bins"        "10"     "$(awk 'NR>1' "${D7}" | wc -l | tr -d ' ')"
assert_eq "distribution: top bin label"  "90-100" "$(awk -F'\t' 'END {print $1}' "${D7}")"
# All three sites are at 50% -> they must all land in the 50-59 bin.
assert_eq "distribution: all sites in 50-59" "3" \
          "$(awk -F'\t' '$1=="50-59" {print $2}' "${D7}")"

# Boundary case for binning: 100% must fall in the top bin, not overflow.
DIR="$(make_config boundaries)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/06_summary_stats.sh"    "${DIR}/config.yaml" >/dev/null 2>&1
D7B="${DIR}/results/03_global_methylation/TEST_01.methylation_distribution.tsv"
# Retained percentages are 50, 0, 19.99, 20.00, 79.99, 80.00, 100
assert_eq "binning: 0% in bin 0-9"        "1" "$(awk -F'\t' '$1=="0-9"    {print $2}' "${D7B}")"
assert_eq "binning: 19.99% in bin 10-19"  "1" "$(awk -F'\t' '$1=="10-19"  {print $2}' "${D7B}")"
assert_eq "binning: 20.00% in bin 20-29"  "1" "$(awk -F'\t' '$1=="20-29"  {print $2}' "${D7B}")"
assert_eq "binning: 79.99% in bin 70-79"  "1" "$(awk -F'\t' '$1=="70-79"  {print $2}' "${D7B}")"
assert_eq "binning: 80.00% in bin 80-89"  "1" "$(awk -F'\t' '$1=="80-89"  {print $2}' "${D7B}")"
assert_eq "binning: 100% folded into 90-100" "1" "$(awk -F'\t' '$1=="90-100" {print $2}' "${D7B}")"
assert_eq "binning: bins sum to retained"    "0" \
          "$(qc6_val "${DIR}/results" ASSERT_distribution_sums_to_sites)"

# Cross-stage disagreement must be fatal. Corrupt stage 03's recorded value and
# confirm stage 06 refuses to proceed rather than quietly reporting its own.
DIR="$(make_config known_values)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
sed -i.bak 's/^Coverage_weighted_methylation_percent\t.*/Coverage_weighted_methylation_percent\t99.9999/' \
    "${DIR}/results/02_qc/TEST_01.03_filter_qc.tsv"
bash "${REPO_DIR}/scripts/bash/06_summary_stats.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "cross-stage disagreement: exits 1" "1" "$?"

# Stage 06 must refuse to run before stage 03.
DIR="$(make_config known_values)"
bash "${REPO_DIR}/scripts/bash/06_summary_stats.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "missing stage 03 output: exits 1" "1" "$?"

# gene_type breakdown must partition each feature class exactly: the per-type
# feature counts have to sum to the pooled row, or a type is being dropped or
# double-counted.
if command -v bedtools >/dev/null 2>&1; then
    # min_cpgs_per_feature lowered to 1: the fixture's features carry 1-2 CpGs,
    # so at the default of 5 every row would be excluded and the comparison
    # would be 0 == 0 — passing without testing anything.
    DIR="$(make_config annotation "s|^  min_cpgs_per_feature: 5|  min_cpgs_per_feature: 1|")"
    SKIP_INTEGRITY=true bash "${REPO_DIR}/scripts/bash/04_run_all.sh" "${DIR}/config.yaml" >/dev/null 2>&1
    CLS="${DIR}/results/03_global_methylation/TEST_01.feature_class_summary.tsv"
    TYP="${DIR}/results/03_global_methylation/TEST_01.methylation_by_gene_type.tsv"
    if [[ -f "${CLS}" && -f "${TYP}" ]]; then
        for cls in gene_body promoter; do
            pooled="$(awk -F'\t' -v c="${cls}" '$1==c {print $2; exit}' "${CLS}")"
            bytype="$(awk -F'\t' -v c="${cls}" 'NR>1 && $1==c {s+=$3} END {print s+0}' "${TYP}")"
            assert_eq "${cls}: gene_type counts sum to pooled" "${pooled}" "${bytype}"
        done
        # protein_coding must appear as its own promoted row in the main table.
        pc="$(awk -F'\t' '$1=="promoter_protein_coding" {print $2; exit}' "${CLS}")"
        pct="$(awk -F'\t' '$1=="promoter" && $2=="protein_coding" {print $3; exit}' "${TYP}")"
        assert_eq "promoter_protein_coding row matches by-type row" "${pc}" "${pct}"
    else
        bad "gene_type breakdown files not produced"
    fi
fi

# Per-chromosome table must split contigs correctly.
DIR="$(make_config no_chry)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/06_summary_stats.sh"    "${DIR}/config.yaml" >/dev/null 2>&1
C7="${DIR}/results/04_chromosome_summary/TEST_01.chromosome_methylation_summary.tsv"
assert_eq "chromosome table: 2 rows"      "2" "$(awk 'NR>1' "${C7}" | wc -l | tr -d ' ')"
assert_eq "chromosome table: chr1 and chrX" "chr1 chrX" \
          "$(awk -F'\t' 'NR>1 {printf "%s ", $1}' "${C7}" | sed 's/ $//')"
echo ""

# =============================================================================
# TEST GROUP 8 — stage 04: orchestration
# =============================================================================
echo "--- [8] stage 04: orchestration ---"

RUNALL="${REPO_DIR}/scripts/bash/04_run_all.sh"

bash "${RUNALL}" --list >/dev/null 2>&1
assert_eq "--list exits 0" "0" "$?"

DIR="$(make_config known_values)"
assert_eq "--dry-run exits 0" "0" \
          "$(bash "${RUNALL}" "${DIR}/config.yaml" --dry-run >/dev/null 2>&1; echo $?)"
# A dry run must not create any results.
[[ -f "${DIR}/results/01_filtered/TEST_01.5mC.all_cov.bedmethyl.gz" ]] \
    && bad "--dry-run produced output" || ok "--dry-run produced no output"

# Full chain in one command. Uses the `annotation` fixture, not `known_values`:
# the latter's CpGs lie outside every test feature, and stage 05 treats zero
# overlaps as fatal (for a human genome it means a contig-naming fault, not
# biology). That is the correct behaviour, so the fixture has to overlap.
DIR="$(make_config annotation)"
SKIP_INTEGRITY=true bash "${RUNALL}" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "full chain exits 0" "0" "$?"
M8="${DIR}/results/TEST_01.run_manifest.tsv"
[[ -f "${M8}" ]] && ok "manifest written" || bad "manifest missing"
if [[ -f "${M8}" ]]; then
    assert_eq "manifest: overall status" "success" \
              "$(awk -F'\t' '$1=="overall_status" {print $2}' "${M8}")"
    # Assert the stages themselves, not a count: the count changes whenever a
    # new stage script is added, which made this test stale the moment stage 07
    # was written.
    for st in 01 02 03 05 06; do
        assert_eq "manifest: stage ${st} ok" "ok" \
                  "$(awk -F'\t' -v s="${st}" '$1=="stage" && $2==s {print $3; exit}' "${M8}")"
    done
    # Provenance the methods section depends on.
    for field in config_sha256 gene_bed_sha256 genome_build min_coverage contig_selection; do
        v="$(awk -F'\t' -v k="${field}" '$1==k {print $2; exit}' "${M8}")"
        [[ -n "${v}" ]] && ok "manifest records ${field}" || bad "manifest missing ${field}"
    done
    # git_commit must never be the literal "HEAD": that is what rev-parse prints
    # for a repo with no commits, and recording it would look like a real SHA.
    GITV="$(awk -F'\t' '$1=="git_commit" {print $2; exit}' "${M8}")"
    [[ "${GITV}" == "HEAD" ]] && bad "manifest git_commit is the literal 'HEAD'" \
                              || ok "manifest git_commit sane (${GITV})"
fi
# The full chain must produce the same figures as running stages individually.
# annotation fixture: three sites at coverage 20 with 15, 8 and 2 modified calls
# -> 25/60 = 41.6667%.
assert_eq "full chain: weighted methylation" "41.6667" \
          "$(awk -F'\t' '$1=="Coverage_weighted_methylation_percent" {print $2; exit}' \
             "${DIR}/results/03_global_methylation/TEST_01.global_methylation_summary.tsv")"

# Zero overlaps must fail loudly, with its own named counter so the report is
# not silent. known_values has no CpG inside any test feature.
DIR="$(make_config known_values)"
bash "${REPO_DIR}/scripts/bash/02_bedmethyl_to_tsv.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/03_filter_cpg_sites.sh" "${DIR}/config.yaml" >/dev/null 2>&1
bash "${REPO_DIR}/scripts/bash/05_annotate_regions.sh" "${DIR}/config.yaml" >"${SCRATCH}/nooverlap.txt" 2>&1
assert_eq "zero overlaps: exits 1" "1" "$?"
assert_eq "zero overlaps: named counter set" "1" \
          "$(awk -F'\t' '$1=="ASSERT_gene_no_overlaps" {print $2; exit}' \
             "${DIR}/results/02_qc/TEST_01.05_annotation_qc.tsv")"
grep -q "contig naming" "${SCRATCH}/nooverlap.txt" \
    && ok "zero overlaps: diagnostic names contig naming as the likely cause" \
    || bad "zero overlaps: diagnostic unhelpful"

# Fail-fast: a malformed input must stop the chain, not carry on.
DIR="$(make_config bad_columns)"
SKIP_INTEGRITY=true bash "${RUNALL}" "${DIR}/config.yaml" >"${SCRATCH}/runall_fail.txt" 2>&1
assert_eq "fail-fast: exits 1" "1" "$?"
assert_eq "fail-fast: no stage 06 output" "0" \
          "$([[ -f "${DIR}/results/03_global_methylation/TEST_01.global_methylation_summary.tsv" ]] && echo 1 || echo 0)"
M8F="${DIR}/results/TEST_01.run_manifest.tsv"
[[ -f "${M8F}" ]] && assert_eq "fail-fast: manifest records failure" "failed" \
          "$(awk -F'\t' '$1=="overall_status" {print $2}' "${M8F}")"

# --only must run exactly one stage.
DIR="$(make_config known_values)"
SKIP_INTEGRITY=true bash "${RUNALL}" "${DIR}/config.yaml" --only 01 >/dev/null 2>&1
assert_eq "--only 01 exits 0" "0" "$?"
assert_eq "--only 01: one stage recorded" "1" \
          "$(awk -F'\t' '$1=="stage"' "${DIR}/results/TEST_01.run_manifest.tsv" | wc -l | tr -d ' ')"
echo ""

# =============================================================================
# TEST GROUP 9 — stage 08: archiving
# =============================================================================
echo "--- [9] stage 08: archiving ---"

ARCHIVER="${REPO_DIR}/scripts/bash/08_archive_results.sh"

# Build a results tree holding TWO samples. results/ is shared in real use, so
# cross-sample contamination is the failure mode that matters most here.
setup_archive() { # setup_archive [extra_sed]
    local extra="${1:-}"
    local d
    d="$(mktemp -d "${SCRATCH}/archive.XXXXXX")"
    mkdir -p "${d}/results/01_filtered" "${d}/results/02_qc" \
             "${d}/results/06_gene_summary" "${d}/logs" "${d}/archive" "${d}/work"
    local s
    for s in TEST_01 OTHER_SAMPLE; do
        echo "cov10-${s}"  > "${d}/results/01_filtered/${s}.5mC.cov10.bedmethyl.gz"
        echo "allcov-${s}" > "${d}/results/01_filtered/${s}.5mC.all_cov.bedmethyl.gz"
        echo "qc-${s}"     > "${d}/results/02_qc/${s}.03_filter_qc.tsv"
        echo "gene-${s}"   > "${d}/results/06_gene_summary/${s}.gene_methylation_summary.tsv"
    done
    echo "shared-not-sample-prefixed" > "${d}/results/02_qc/shared_notes.txt"
    echo "stage log" > "${d}/logs/06_summary_stats_prior.log"
    sed -e "s|^  output_dir: .*|  output_dir: \"${d}/results\"|" \
        -e "s|^  work_dir: .*|  work_dir: \"${d}/work\"|" \
        -e "s|^  log_dir: .*|  log_dir: \"${d}/logs\"|" \
        -e "s|^  archive_root: .*|  archive_root: \"${d}/archive\"|" \
        "$(make_config known_values)/config.yaml" > "${d}/config.yaml"
    # Applied as a separate pass, not as ${extra:+-e "$extra"}: unquoted
    # expansion word-splits the sed expression on its spaces, so the edit
    # silently does nothing and the test passes against unchanged config.
    if [[ -n "${extra}" ]]; then
        sed -i.bak "${extra}" "${d}/config.yaml" && rm -f "${d}/config.yaml.bak"
    fi
    echo "${d}"
}

D9="$(setup_archive)"
bash "${ARCHIVER}" "${D9}/config.yaml" --dry-run >/dev/null 2>&1
assert_eq "archive --dry-run exits 0" "0" "$?"
assert_eq "archive --dry-run copies nothing" "0" \
          "$(find "${D9}/archive" -type f 2>/dev/null | wc -l | tr -d ' ')"

D9="$(setup_archive)"
bash "${ARCHIVER}" "${D9}/config.yaml" >/dev/null 2>&1
assert_eq "archive exits 0" "0" "$?"

# THE assertion this stage exists for: not one byte of the other sample.
assert_eq "archive: no OTHER_SAMPLE files" "0" \
          "$(find "${D9}/archive" -name 'OTHER_SAMPLE*' | wc -l | tr -d ' ')"
assert_eq "archive: no unprefixed shared file" "0" \
          "$(find "${D9}/archive/"*/*/results -type f ! -name 'TEST_01.*' 2>/dev/null | wc -l | tr -d ' ')"
# Large regenerable intermediate excluded by default.
assert_eq "archive: all_cov excluded by default" "0" \
          "$(find "${D9}/archive" -name '*all_cov*' | wc -l | tr -d ' ')"
assert_eq "archive: 3 result files" "3" \
          "$(find "${D9}/archive/"*/*/results -type f | wc -l | tr -d ' ')"

# Checksums must exist AND verify against what was written.
CK="$(find "${D9}/archive" -name checksums.sha256 | head -n1)"
[[ -s "${CK}" ]] && ok "archive: checksums written" || bad "archive: checksums missing"
if [[ -s "${CK}" ]]; then
    RUNDIR="$(dirname "${CK}")"
    if command -v sha256sum >/dev/null 2>&1; then
        ( cd "${RUNDIR}" && sha256sum -c --quiet checksums.sha256 ) >/dev/null 2>&1
    else
        ( cd "${RUNDIR}" && shasum -a 256 -c --status checksums.sha256 ) >/dev/null 2>&1
    fi
    assert_eq "archive: checksums verify" "0" "$?"
    # Corrupting an archived file must make verification fail — otherwise the
    # checksum file is decorative.
    echo "tampered" >> "${RUNDIR}/results/01_filtered/TEST_01.5mC.cov10.bedmethyl.gz"
    if command -v sha256sum >/dev/null 2>&1; then
        ( cd "${RUNDIR}" && sha256sum -c --quiet checksums.sha256 ) >/dev/null 2>&1
    else
        ( cd "${RUNDIR}" && shasum -a 256 -c --status checksums.sha256 ) >/dev/null 2>&1
    fi
    assert_eq "archive: checksums detect tampering" "1" "$?"
fi

# Provenance and index.
PV="$(find "${D9}/archive" -name provenance.tsv | head -n1)"
[[ -s "${PV}" ]] && ok "archive: provenance written" || bad "archive: provenance missing"
if [[ -s "${PV}" ]]; then
    assert_eq "provenance: records sample_id" "TEST_01" \
              "$(awk -F'\t' '$1=="sample_id" {print $2}' "${PV}")"
    assert_eq "provenance: records checksum status" "verified:5" \
              "$(awk -F'\t' '$1=="checksums" {print $2}' "${PV}")"
fi
IDX="${D9}/archive/archive_index.tsv"
[[ -s "${IDX}" ]] && ok "archive: index created" || bad "archive: index missing"
assert_eq "index: one run recorded" "1" "$(awk 'NR>1' "${IDX}" | wc -l | tr -d ' ')"

# A second archive of the same sample must APPEND, not overwrite: the point of
# the store is a longitudinal record.
sleep 1
bash "${ARCHIVER}" "${D9}/config.yaml" >/dev/null 2>&1
assert_eq "second archive exits 0" "0" "$?"
assert_eq "index: two runs recorded" "2" "$(awk 'NR>1' "${IDX}" | wc -l | tr -d ' ')"
assert_eq "archive: two run directories kept" "2" \
          "$(find "${D9}/archive/TEST_01" -mindepth 1 -maxdepth 1 -type d | wc -l | tr -d ' ')"

# include_intermediates=true must pull the all_cov file in.
D9B="$(setup_archive 's|^  include_intermediates: false|  include_intermediates: true|')"
bash "${ARCHIVER}" "${D9B}/config.yaml" >/dev/null 2>&1
assert_eq "include_intermediates: all_cov archived" "1" \
          "$(find "${D9B}/archive" -name '*all_cov*' | wc -l | tr -d ' ')"

# Nothing to archive must fail loudly, not create an empty archive entry.
D9C="$(setup_archive)"
rm -f "${D9C}/results/01_filtered/TEST_01."* "${D9C}/results/02_qc/TEST_01."* \
      "${D9C}/results/06_gene_summary/TEST_01."*
bash "${ARCHIVER}" "${D9C}/config.yaml" >/dev/null 2>&1
assert_eq "no matching files: exits 1" "1" "$?"
assert_eq "no matching files: no index created" "0" \
          "$([[ -f "${D9C}/archive/archive_index.tsv" ]] && echo 1 || echo 0)"

# Placeholder archive_root must be treated as "archiving disabled", not an error.
D9D="$(setup_archive 's|^  archive_root: .*|  archive_root: "/path/to/internal/archive/root"|')"
bash "${ARCHIVER}" "${D9D}/config.yaml" >/dev/null 2>&1
assert_eq "placeholder archive_root: exits 0 (disabled)" "0" "$?"
echo ""

# =============================================================================
# TEST GROUP 10 — stage 07: figures
# =============================================================================
echo "--- [10] stage 07: figures ---"

# Only meaningful with a real R. A stubbed Rscript that exits 0 would make every
# assertion here pass without drawing anything, which is worse than skipping.
R_REAL=0
if command -v Rscript >/dev/null 2>&1; then
    [[ "$(Rscript -e 'cat(2+2)' 2>/dev/null)" == "4" ]] && R_REAL=1
fi

if [[ ${R_REAL} -eq 0 ]]; then
    echo "    [SKIP] no working Rscript — stage 07 tests skipped"
elif ! command -v bedtools >/dev/null 2>&1; then
    echo "    [SKIP] no bedtools — stage 07 tests need stage 05 output"
else
REPORTER="${REPO_DIR}/scripts/R/07_generate_report.R"

# Case 1: no feature meets min_cpgs_per_feature. The fixture's features carry
# 1-2 CpGs against a threshold of 5, so every feature-class vector is empty.
# boxplot() on an empty list throws, so this must be a SKIP and the stage must
# still succeed. This is exactly the case that broke the full-chain test.
DIR="$(make_config annotation)"
SKIP_INTEGRITY=true bash "${REPO_DIR}/scripts/bash/04_run_all.sh" "${DIR}/config.yaml" \
    >"${SCRATCH}/fig_nofeat.txt" 2>&1
assert_eq "figures: chain succeeds when no feature meets min CpGs" "0" "$?"
F10="${DIR}/results/09_figures"
grep -q "feature-class distributions" "${SCRATCH}/fig_nofeat.txt" \
    && ok "figures: fig3 handled (skipped or drawn), not a crash" \
    || bad "figures: fig3 neither drawn nor skipped"
[[ -f "${F10}/TEST_01.fig3_feature_class_methylation.png" ]] \
    && bad "figures: fig3 written despite no qualifying features" \
    || ok "figures: fig3 correctly absent"
# The other figures must still be produced.
for fg in fig1_methylation_distribution fig2_chromosome_methylation \
          fig4_coverage_distribution fig5_methylation_states; do
    [[ -s "${F10}/TEST_01.${fg}.png" ]] && ok "figures: ${fg} written" \
                                        || bad "figures: ${fg} missing or empty"
done
# A zero-byte PNG means a device was opened and never written to.
assert_eq "figures: no zero-byte PNGs" "0" \
          "$(find "${F10}" -name '*.png' -size 0 2>/dev/null | wc -l | tr -d ' ')"
[[ -s "${F10}/TEST_01.figure_manifest.tsv" ]] && ok "figures: manifest written" \
                                              || bad "figures: manifest missing"

# Case 2: features DO qualify, so fig3 and fig6 should both appear.
DIR="$(make_config annotation "s|^  min_cpgs_per_feature: 5|  min_cpgs_per_feature: 1|")"
SKIP_INTEGRITY=true bash "${REPO_DIR}/scripts/bash/04_run_all.sh" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "figures: chain succeeds with qualifying features" "0" "$?"
F10B="${DIR}/results/09_figures"
[[ -s "${F10B}/TEST_01.fig3_feature_class_methylation.png" ]] \
    && ok "figures: fig3 drawn when features qualify" \
    || bad "figures: fig3 missing when features qualify"
[[ -s "${F10B}/TEST_01.fig6_promoter_by_gene_type.png" ]] \
    && ok "figures: fig6 drawn when gene_type rows exist" \
    || bad "figures: fig6 missing when gene_type rows exist"
assert_eq "figures: no zero-byte PNGs (case 2)" "0" \
          "$(find "${F10B}" -name '*.png' -size 0 2>/dev/null | wc -l | tr -d ' ')"

# Case 3: nothing upstream at all. No figure is possible, so the stage must fail
# rather than report success having drawn nothing.
DIR="$(make_config annotation)"
Rscript "${REPORTER}" "${DIR}/config.yaml" >/dev/null 2>&1
assert_eq "figures: exits 1 when no inputs exist" "1" "$?"
fi
echo ""

# =============================================================================
# TEST GROUP 11 — stages 09/10: cross-sample comparison
# =============================================================================
echo "--- [11] stages 09/10: cross-sample comparison ---"

COMPARER="${REPO_DIR}/scripts/bash/10_run_comparison.sh"

# --- Sheet validation (bash only; no R needed) ---
SH_DIR="$(mktemp -d "${SCRATCH}/sheets.XXXXXX")"

# Every invocation below MUST pass --config. Without it, 10_run_comparison.sh
# falls back to config/pipeline_config.yaml — which is GITIGNORED, so it exists
# on a developer's machine and does not exist in a fresh clone. On CI the script
# then died in 00_setup_env.sh before ever reaching the sheet validation, so
# these assertions passed for entirely the wrong reason locally and failed on a
# clean checkout. A test must not depend on an untracked file.
SH_CFG="$(make_config valid)/config.yaml"

bash "${COMPARER}" --template "${SH_DIR}/tpl.tsv" >/dev/null 2>&1
assert_eq "--template exits 0" "0" "$?"
[[ -s "${SH_DIR}/tpl.tsv" ]] && ok "--template writes a sheet" || bad "--template wrote nothing"
# The template must be genuinely tab separated: a space-separated sheet is the
# most likely way a hand-edited file breaks.
assert_eq "template is tab separated" "1" \
          "$(awk -F'\t' 'NR==FNR && /^sample_id/ {print (NF>=4) ? 1 : 0; exit}' "${SH_DIR}/tpl.tsv")"

bash "${COMPARER}" "${SH_DIR}/nonexistent.tsv" --config "${SH_CFG}" >/dev/null 2>&1
assert_eq "missing sheet: exits 1" "1" "$?"

printf 'sample_id\tresults_dir\nONLY_ONE\t/tmp\n' > "${SH_DIR}/one.tsv"
bash "${COMPARER}" "${SH_DIR}/one.tsv" --config "${SH_CFG}" >/dev/null 2>&1
assert_eq "single sample: exits 1" "1" "$?"

# Space- instead of tab-separated header must be rejected with a clear message.
printf 'sample_id results_dir\nA /tmp\nB /tmp\n' > "${SH_DIR}/spaces.tsv"
bash "${COMPARER}" "${SH_DIR}/spaces.tsv" --config "${SH_CFG}" >"${SH_DIR}/spaces.out" 2>&1
assert_eq "space-separated sheet: exits 1" "1" "$?"
grep -q "TAB separated" "${SH_DIR}/spaces.out" \
    && ok "space-separated sheet: message names tabs as the problem" \
    || bad "space-separated sheet: unhelpful message"

# --- Build two real samples through the pipeline ---
if ! command -v bedtools >/dev/null 2>&1; then
    echo "    [SKIP] no bedtools — cannot build per-sample inputs for comparison"
else
CMP="$(mktemp -d "${SCRATCH}/compare.XXXXXX")"
mkdir -p "${CMP}/input_a" "${CMP}/input_b" "${CMP}/results" "${CMP}/logs" "${CMP}/work"
cp "${FIXTURE_DIR}/${PREFIX}.annotation.wf_mods.bedmethyl.gz"   "${CMP}/input_a/${PREFIX}.wf_mods.bedmethyl.gz"
cp "${FIXTURE_DIR}/${PREFIX}.annotation_b.wf_mods.bedmethyl.gz" "${CMP}/input_b/${PREFIX}.wf_mods.bedmethyl.gz"

# Both samples write into ONE shared results dir, which is how the pipeline is
# actually used: outputs are sample-prefixed.
mk_cmp_cfg() { # mk_cmp_cfg <sample_id> <input_dir> <out_file>
    sed -e "s|^  sample_id: .*|  sample_id: \"$1\"|" \
        -e "s|^  raw_sample_prefix: .*|  raw_sample_prefix: \"${PREFIX}\"|" \
        -e "s|^  input_dir: .*|  input_dir: \"$2\"|" \
        -e "s|^  output_dir: .*|  output_dir: \"${CMP}/results\"|" \
        -e "s|^  work_dir: .*|  work_dir: \"${CMP}/work\"|" \
        -e "s|^  repo_dir: .*|  repo_dir: \"${REPO_DIR}\"|" \
        -e "s|^  config_file: .*|  config_file: \"${FIXTURE_DIR}/test_reference_paths.yaml\"|" \
        -e "s|^  log_dir: .*|  log_dir: \"${CMP}/logs\"|" \
        -e "s|^  threads: .*|  threads: 2|" \
        -e "s|^  expected_primary_contigs: .*|  expected_primary_contigs: 1|" \
        -e "s|^  min_cpgs_per_feature: .*|  min_cpgs_per_feature: 1|" \
        -e "s|^  min_shared_cpgs: .*|  min_shared_cpgs: 1|" \
        "${REPO_DIR}/config/pipeline_config.example.yaml" > "$3"
}
mk_cmp_cfg CMP_A "${CMP}/input_a" "${CMP}/cfg_a.yaml"
mk_cmp_cfg CMP_B "${CMP}/input_b" "${CMP}/cfg_b.yaml"

for cf in cfg_a cfg_b; do
    for st in 02_bedmethyl_to_tsv 03_filter_cpg_sites 05_annotate_regions 06_summary_stats; do
        bash "${REPO_DIR}/scripts/bash/${st}.sh" "${CMP}/${cf}.yaml" >/dev/null 2>&1
    done
done
built=1
for s in CMP_A CMP_B; do
    [[ -f "${CMP}/results/06_gene_summary/${s}.gene_methylation_summary.tsv" ]] || built=0
done
assert_eq "two samples built into a shared results dir" "1" "${built}"

printf 'sample_id\tlabel\tgroup\tresults_dir\n' > "${CMP}/sheet.tsv"
printf 'CMP_A\tsampleA\tctrl\t%s/results\n' "${CMP}" >> "${CMP}/sheet.tsv"
printf 'CMP_B\tsampleB\tcase\t%s/results\n' "${CMP}" >> "${CMP}/sheet.tsv"

bash "${COMPARER}" "${CMP}/sheet.tsv" --config "${CMP}/cfg_a.yaml" --dry-run >/dev/null 2>&1
assert_eq "valid sheet --dry-run exits 0" "0" "$?"

# Duplicate sample_id must be rejected: it would compare a sample with itself.
printf 'sample_id\tresults_dir\nCMP_A\t%s/results\nCMP_A\t%s/results\n' "${CMP}" "${CMP}" > "${CMP}/dup.tsv"
bash "${COMPARER}" "${CMP}/dup.tsv" --config "${CMP}/cfg_a.yaml" --dry-run >/dev/null 2>&1
assert_eq "duplicate sample_id: exits 1" "1" "$?"

# A sample whose stage 05 output is absent must be caught before R starts.
printf 'sample_id\tresults_dir\nCMP_A\t%s/results\nNOT_RUN\t%s/results\n' "${CMP}" "${CMP}" > "${CMP}/absent.tsv"
bash "${COMPARER}" "${CMP}/absent.tsv" --config "${CMP}/cfg_a.yaml" --dry-run >"${CMP}/absent.out" 2>&1
assert_eq "unrun sample: exits 1" "1" "$?"
grep -q "run stages 05 and 06" "${CMP}/absent.out" \
    && ok "unrun sample: message names the stages to run" \
    || bad "unrun sample: unhelpful message"

# --- The comparison itself (needs real R) ---
R_REAL=0
if command -v Rscript >/dev/null 2>&1; then
    [[ "$(Rscript -e 'cat(2+2)' 2>/dev/null)" == "4" ]] && R_REAL=1
fi
if [[ ${R_REAL} -eq 0 ]]; then
    echo "    [SKIP] no working Rscript — stage 09 comparison not executed"
else
    bash "${COMPARER}" "${CMP}/sheet.tsv" --config "${CMP}/cfg_a.yaml" \
         --name testcmp --out "${CMP}/cmp_out" >"${CMP}/cmp.log" 2>&1
    assert_eq "comparison exits 0" "0" "$?"
    O="${CMP}/cmp_out"
    for fn in global_comparison.tsv gene_comparison.tsv promoter_comparison.tsv \
              cpg_island_comparison.tsv comparison_summary.tsv sample_sheet.tsv; do
        [[ -s "${O}/${fn}" ]] && ok "comparison: ${fn} written" \
                              || bad "comparison: ${fn} missing"
    done
    [[ -s "${O}/gene_scatter.png" ]] && ok "comparison: scatter figure written" \
                                     || bad "comparison: scatter figure missing"

    # Arithmetic: chr1:1000 sits in ENSGA.1 and ENSGB.1, and goes 75% -> 10%,
    # so delta must be -65 exactly. Sign matters: positive would mean the
    # comparison is subtracting in the wrong direction.
    GC="${O}/gene_comparison.tsv"
    assert_eq "gene ENSGA.1 delta is -65" "-65" \
              "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="delta") c=i; next}
                             $1=="ENSGA.1" {printf "%g", $c}' "${GC}")"
    assert_eq "gene ENSGA.1 direction is hypo" "hypo" \
              "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="direction") c=i; next}
                             $1=="ENSGA.1" {print $c}' "${GC}")"
    # ENSGC.1 moves 40% -> 45%, only 5 pp, so it must NOT be called differential
    # even if the p-value is small. This is the effect-size threshold working.
    assert_eq "gene ENSGC.1 not differential (5 pp)" "1" \
              "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="differential") c=i; next}
                             $1=="ENSGC.1" {print ($c!="yes") ? 1 : 0}' "${GC}")"
    # Both samples' methylation columns must be present and correct.
    assert_eq "gene table has per-sample columns" "1" \
              "$(head -n1 "${GC}" | awk -F'\t' '{a=0;b=0; for(i=1;i<=NF;i++){if($i=="methyl_sampleA")a=1; if($i=="methyl_sampleB")b=1}; print (a&&b)?1:0}')"
    assert_eq "gene ENSGA.1 sampleA methylation is 75" "75" \
              "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="methyl_sampleA") c=i; next}
                             $1=="ENSGA.1" {printf "%g", $c}' "${GC}")"
    # Balance columns must exist, and a like-for-like feature must be balanced.
    assert_eq "gene table has balance columns" "1" \
              "$(head -n1 "${GC}" | awk -F'\t' '{a=0;b=0;c=0; for(i=1;i<=NF;i++){if($i=="cpg_ratio")a=1; if($i=="coverage_ratio")b=1; if($i=="balanced")c=1}; print (a&&b&&c)?1:0}')"
    assert_eq "gene ENSGA.1 is balanced (1 CpG each)" "yes" \
              "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="balanced") c=i; next}
                             $1=="ENSGA.1" {print $c}' "${GC}")"

    # The caveat must be surfaced in the run output, not buried in a doc.
    grep -q "anti-conservative" "${CMP}/cmp.log" \
        && ok "comparison: statistical caveat printed to the run log" \
        || bad "comparison: caveat not surfaced"

    # --- Coverage-imbalance filter ---
    # 6 covered CpGs against 1 in the same gene: cpg_ratio 6 > max_cpg_ratio 3.
    # A 70 pp difference must NOT be called differential on that basis. This is
    # the mappability artefact that topped the first real two-sample run.
    CMPI="$(mktemp -d "${SCRATCH}/imbal.XXXXXX")"
    mkdir -p "${CMPI}/input_c"
    cp "${FIXTURE_DIR}/${PREFIX}.annotation_imbalanced.wf_mods.bedmethyl.gz" \
       "${CMPI}/input_c/${PREFIX}.wf_mods.bedmethyl.gz"
    sed -e "s|^  sample_id: .*|  sample_id: \"CMP_C\"|" \
        -e "s|^  input_dir: .*|  input_dir: \"${CMPI}/input_c\"|" \
        "${CMP}/cfg_a.yaml" > "${CMPI}/cfg_c.yaml"
    for st in 02_bedmethyl_to_tsv 03_filter_cpg_sites 05_annotate_regions 06_summary_stats; do
        bash "${REPO_DIR}/scripts/bash/${st}.sh" "${CMPI}/cfg_c.yaml" >/dev/null 2>&1
    done
    printf 'sample_id\tlabel\tgroup\tresults_dir\n' > "${CMPI}/sheet.tsv"
    printf 'CMP_A\tsampleA\tctrl\t%s/results\n' "${CMP}"  >> "${CMPI}/sheet.tsv"
    printf 'CMP_C\tsampleC\tcase\t%s/results\n' "${CMP}"  >> "${CMPI}/sheet.tsv"
    bash "${COMPARER}" "${CMPI}/sheet.tsv" --config "${CMPI}/cfg_c.yaml" \
         --out "${CMPI}/out" >"${CMPI}/log" 2>&1
    assert_eq "imbalance: comparison exits 0" "0" "$?"
    GCI="${CMPI}/out/gene_comparison.tsv"
    if [[ -s "${GCI}" ]]; then
        assert_eq "imbalance: ENSGA.1 cpg_ratio is 6" "6" \
                  "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="cpg_ratio") c=i; next}
                                 $1=="ENSGA.1" {printf "%g", $c}' "${GCI}")"
        assert_eq "imbalance: ENSGA.1 flagged not balanced" "no" \
                  "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="balanced") c=i; next}
                                 $1=="ENSGA.1" {print $c}' "${GCI}")"
        # 70 pp difference, but must be "imbalanced", never "yes".
        assert_eq "imbalance: ENSGA.1 labelled imbalanced not differential" "imbalanced" \
                  "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="differential") c=i; next}
                                 $1=="ENSGA.1" {print $c}' "${GCI}")"
        # ENSGC.1 has 1 CpG in both samples, so it stays balanced.
        assert_eq "imbalance: ENSGC.1 still balanced" "yes" \
                  "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="balanced") c=i; next}
                                 $1=="ENSGC.1" {print $c}' "${GCI}")"
        # Balanced features must sort above imbalanced ones.
        assert_eq "imbalance: balanced features sort first" "yes" \
                  "$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="balanced") c=i; next}
                                 NR==2 {print $c; exit}' "${GCI}")"
        grep -q "imbalanced" "${CMPI}/log" \
            && ok "imbalance: count reported in the run log" \
            || bad "imbalance: not reported in the run log"
    else
        bad "imbalance: gene_comparison.tsv not produced"
    fi
fi
fi
echo ""

# =============================================================================
echo "============================================================"
printf " PASS: %d    FAIL: %d\n" "${PASS}" "${FAIL}"
echo "============================================================"

# Nothing may have leaked into the real repo.
LEAKED=0
for d in "${REPO_DIR}/results" "${REPO_DIR}/logs"; do
    if find "${d}" -name 'TEST_01.*' -o -name 'TEST_SAMPLE.*' 2>/dev/null | grep -q .; then
        echo "[WARN] test artefacts found in ${d}"
        LEAKED=1
    fi
done
[[ ${LEAKED} -eq 0 ]] && echo "No test artefacts leaked into the repo."

[[ ${FAIL} -eq 0 ]] || exit 1
