#!/usr/bin/env bash
# =============================================================================
# 02_bedmethyl_to_tsv.sh
# Extracts the primary modification code from the raw bedMethyl in a SINGLE
# pass, drops excluded contigs, and emits a bgzipped, tabix-indexed file plus
# the QC needed to justify the downstream coverage threshold.
#
# Deliberately does NOT apply a coverage filter. That belongs in stage 03, so
# that the threshold can be varied — or a sensitivity analysis run for a
# reviewer — without re-reading the full input.
#
# The output keeps the NATIVE 18-column bedMethyl layout rather than trimming to
# the four informative columns. This is intentional: it means the column numbers
# used in the Methylation_Demo walkthrough, in this pipeline, and in the test
# assertions are all the same numbers. One source of truth for column positions
# is worth more than the disk space saved by trimming.
#
# Assertions made in the same pass (all free, since every record is read anyway):
#   - every record has exactly modifications.expected_columns columns
#   - coverage and percent fields are numeric, percent within [0, 100]
#   - modified + canonical counts never exceed valid coverage
#   - records are coordinate-sorted and contigs appear in contiguous blocks
#     (both required by tabix and by `bedtools intersect -sorted`)
#
# Usage:
#   bash scripts/bash/02_bedmethyl_to_tsv.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# Drop a '#'-leading argument: zsh passes inline comments through as arguments.
# See the note in 00_setup_env.sh.
_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in '#'*) _CFG_ARG="" ;; esac
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/02_bedmethyl_to_tsv_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo ""
echo "=== [02_bedmethyl_to_tsv.sh] Extracting mod code '${PRIMARY_MOD_CODE}' for ${SAMPLE_ID} ==="

# --- Human-readable label for the modification code -------------------------
case "${PRIMARY_MOD_CODE}" in
    m) MOD_LABEL="5mC"  ;;
    h) MOD_LABEL="5hmC" ;;
    a) MOD_LABEL="6mA"  ;;
    *) MOD_LABEL="mod_${PRIMARY_MOD_CODE}" ;;
esac
echo "Modification  : ${PRIMARY_MOD_CODE} (${MOD_LABEL})"

# --- Resolve input -----------------------------------------------------------
if [[ "${PHASED}" == "true" ]]; then
    echo "[FATAL] modifications.phased=true is not yet handled by this stage." >&2
    echo "        This sample set is unphased; haplotype-resolved extraction is" >&2
    echo "        a separate stage rather than a branch inside this one." >&2
    exit 1
fi

IN_BEDMETHYL="${INPUT_DIR}/${MOD_BEDMETHYL_NAME}"
if [[ ! -f "${IN_BEDMETHYL}" ]]; then
    echo "[FATAL] Input bedMethyl not found: ${IN_BEDMETHYL}" >&2
    echo "        Run 01_validate_inputs.sh first." >&2
    exit 1
fi

# --- Output paths ------------------------------------------------------------
FILTERED_DIR="${OUTPUT_DIR}/01_filtered"
QC_DIR="${OUTPUT_DIR}/02_qc"
mkdir -p "${FILTERED_DIR}" "${QC_DIR}" || {
    echo "[FATAL] Could not create output directories under ${OUTPUT_DIR}" >&2; exit 1; }

OUT_BED="${FILTERED_DIR}/${SAMPLE_ID}.${MOD_LABEL}.all_cov.bedmethyl.gz"
QC_MAIN="${QC_DIR}/${SAMPLE_ID}.02_extraction_qc.tsv"
QC_HIST="${QC_DIR}/${SAMPLE_ID}.02_coverage_histogram.tsv"
QC_CONTIG="${QC_DIR}/${SAMPLE_ID}.02_contig_inventory.tsv"
QC_EXCL="${QC_DIR}/${SAMPLE_ID}.02_excluded_contigs.tsv"

# --- Contig selection mode ---------------------------------------------------
# Exactly one mechanism is active. Stating which, every run, means the log alone
# is enough to know what was analysed — no need to go back to the config.
if [[ "${PRIMARY_CONTIGS_ONLY}" == "true" ]]; then
    CONTIG_MODE="allowlist"
    USE_ALLOWLIST=1
    echo "Contig filter : ALLOWLIST — keeping only ${INCLUDE_CONTIGS_REGEX}"
    echo "                (exclude_contigs_regex is ignored in this mode)"
else
    CONTIG_MODE="blocklist"
    USE_ALLOWLIST=0
    echo "Contig filter : BLOCKLIST — dropping ${EXCLUDE_CONTIGS_REGEX}"
    echo "                [WARN] a blocklist only rejects the patterns it names."
    echo "                GRCh38 analysis sets carry scaffolds such as GL000191.1"
    echo "                and HLA-A*01:01:01:01 that match none of them. Prefer"
    echo "                filtering.primary_contigs_only: true"
fi

# --- Pick the fastest available decompressor --------------------------------
# pigz is multi-threaded; on a 660 MB input that is a meaningful saving. gzip is
# the guaranteed fallback so the pipeline never depends on pigz being present.
if command -v pigz >/dev/null 2>&1; then
    DECOMP=(pigz -dc -p "${THREADS}")
    echo "Decompressor  : pigz -p ${THREADS}"
else
    DECOMP=(gzip -dc)
    echo "Decompressor  : gzip (install pigz for a multi-threaded read)"
fi

echo "Input         : ${IN_BEDMETHYL}"
echo "Output        : ${OUT_BED}"
echo ""
echo "Reading every record — this is the only full pass over the raw file."
T0=$(date +%s)

# --- Single pass -------------------------------------------------------------
# Data goes to stdout (-> bgzip); QC goes to files. Records are passed through
# with `print $0`, unmodified, so no field is ever reformatted or rounded.
set -o pipefail
"${DECOMP[@]}" "${IN_BEDMETHYL}" \
| awk -F'\t' \
      -v expected="${EXPECTED_COLUMNS}" \
      -v primary="${PRIMARY_MOD_CODE}" \
      -v excl="${EXCLUDE_CONTIGS_REGEX}" \
      -v incl="${INCLUDE_CONTIGS_REGEX}" \
      -v use_allowlist="${USE_ALLOWLIST}" \
      -v contig_mode="${CONTIG_MODE}" \
      -v expected_contigs="${EXPECTED_PRIMARY_CONTIGS}" \
      -v qc_main="${QC_MAIN}" \
      -v qc_hist="${QC_HIST}" \
      -v qc_contig="${QC_CONTIG}" \
      -v qc_excl="${QC_EXCL}" \
      -v sample="${SAMPLE_ID}" \
      -v modlabel="${MOD_LABEL}" '
BEGIN {
    OFS = "\t"
    HIST_MAX = 200          # coverage values above this go into one overflow bin
    prev_chrom = ""
    prev_start = -1
}
{
    total_records++

    # --- structural assertions, applied to EVERY record ---
    if (NF != expected) { bad_nf++; next }

    code = $4
    code_count[code]++
    if (!(code in code_seen)) { code_seen[code] = 1; n_codes++ }

    if (code != primary) next
    primary_records++

    chrom = $1

    # --- contig selection: allowlist OR blocklist, never both ---
    if (use_allowlist == 1) {
        drop = (chrom !~ incl)
    } else {
        drop = (excl != "" && chrom ~ excl)
    }
    if (drop) {
        excluded_records++
        excluded_contig[chrom]++
        if (!(chrom in excl_seen)) { excl_seen[chrom] = 1; n_excl_contigs++ }
        next
    }

    start = $2 + 0
    cov   = $10
    pct   = $11
    nmod  = $12
    ncan  = $13

    # --- value assertions ---
    if (cov !~ /^[0-9]+$/)                  { bad_cov++;       next }
    if (pct !~ /^[0-9]+(\.[0-9]+)?$/)       { bad_pct++;       next }
    if (pct + 0 < 0 || pct + 0 > 100)       { bad_pct_range++; next }
    if (nmod + ncan > cov + 0)              { bad_counts++ }

    # --- sortedness: required by tabix and by bedtools -sorted ---
    if (chrom != prev_chrom) {
        if (chrom in contig_seen) {
            interleaved++          # contig reappears after another contig
        } else {
            contig_seen[chrom] = 1
            n_contigs++
            contig_order[n_contigs] = chrom
        }
        prev_start = -1
    }
    if (start < prev_start) unsorted++
    prev_chrom = chrom
    prev_start = start

    # --- tallies ---
    kept++
    contig_sites[chrom]++
    contig_cov[chrom]  += cov
    contig_mod[chrom]  += nmod
    sum_cov            += cov
    sum_mod            += nmod
    sum_pct            += pct

    if (kept == 1 || cov + 0 < min_cov) min_cov = cov + 0
    if (kept == 1 || cov + 0 > max_cov) max_cov = cov + 0

    b = (cov + 0 > HIST_MAX) ? HIST_MAX + 1 : cov + 0
    hist[b]++

    print $0
}
END {
    # ---------------- coverage histogram + cumulative ----------------
    print "Coverage", "Sites" > qc_hist
    for (c = 0; c <= HIST_MAX; c++)
        printf "%d\t%d\n", c, hist[c] + 0 >> qc_hist
    printf "%d+\t%d\n", HIST_MAX + 1, hist[HIST_MAX + 1] + 0 >> qc_hist
    close(qc_hist)

    # ---------------- per-contig inventory ----------------
    print "Chromosome", "CpG_sites", "Total_coverage", "Modified_read_calls",
          "Coverage_weighted_methylation_percent" > qc_contig
    for (i = 1; i <= n_contigs; i++) {
        c = contig_order[i]
        printf "%s\t%d\t%.0f\t%.0f\t%.4f\n", c, contig_sites[c],
               contig_cov[c], contig_mod[c],
               (contig_cov[c] > 0 ? (contig_mod[c] / contig_cov[c]) * 100 : 0) >> qc_contig
    }
    close(qc_contig)

    # ---------------- main QC ----------------
    print "Metric", "Value" > qc_main
    printf "Sample\t%s\n",                        sample        >> qc_main
    printf "Modification\t%s\n",                  modlabel      >> qc_main
    printf "Total_records_read\t%d\n",            total_records >> qc_main
    printf "Distinct_modification_codes\t%d\n",   n_codes       >> qc_main
    for (k in code_count)
        printf "Records_mod_code_%s\t%d\n", k, code_count[k]    >> qc_main
    printf "Primary_code_records\t%d\n",          primary_records >> qc_main
    printf "Contig_selection_mode\t%s\n",         contig_mode   >> qc_main
    printf "Contig_selection_pattern\t%s\n",      (use_allowlist == 1 ? incl : excl) >> qc_main
    printf "Excluded_by_contig_filter\t%d\n",     excluded_records + 0 >> qc_main
    printf "Excluded_contigs\t%d\n",              n_excl_contigs + 0   >> qc_main
    printf "Retained_CpG_sites\t%d\n",            kept + 0      >> qc_main
    printf "Contigs_retained\t%d\n",              n_contigs + 0 >> qc_main
    printf "Minimum_coverage\t%d\n",              min_cov + 0   >> qc_main
    printf "Maximum_coverage\t%d\n",              max_cov + 0   >> qc_main
    printf "Mean_coverage\t%.4f\n",               (kept ? sum_cov / kept : 0) >> qc_main
    printf "Total_valid_coverage\t%.0f\n",        sum_cov + 0   >> qc_main
    printf "Total_modified_read_calls\t%.0f\n",   sum_mod + 0   >> qc_main
    printf "Unweighted_mean_site_methylation_percent\t%.4f\n",
           (kept ? sum_pct / kept : 0) >> qc_main
    printf "Coverage_weighted_methylation_percent\t%.4f\n",
           (sum_cov > 0 ? (sum_mod / sum_cov) * 100 : 0) >> qc_main

    # ---------------- assertion results ----------------
    printf "ASSERT_wrong_column_count\t%d\n",     bad_nf + 0        >> qc_main
    printf "ASSERT_non_numeric_coverage\t%d\n",   bad_cov + 0       >> qc_main
    printf "ASSERT_non_numeric_percent\t%d\n",    bad_pct + 0       >> qc_main
    printf "ASSERT_percent_out_of_range\t%d\n",   bad_pct_range + 0 >> qc_main
    printf "ASSERT_counts_exceed_coverage\t%d\n", bad_counts + 0    >> qc_main
    printf "ASSERT_unsorted_records\t%d\n",       unsorted + 0      >> qc_main
    printf "ASSERT_interleaved_contigs\t%d\n",    interleaved + 0   >> qc_main
    # Contig-count mismatch is a WARNING, not an error: 23 instead of 24 is the
    # expected result for a female sample with no chrY calls.
    printf "WARN_unexpected_contig_count\t%d\n",
           ((expected_contigs + 0 > 0 && n_contigs != expected_contigs + 0) ? 1 : 0) >> qc_main
    close(qc_main)

    # ---------------- excluded-contig detail ----------------
    # Written to its own QC file rather than dumped to the terminal. On GRCh38
    # this is 127 lines of scaffold names, which buries the numbers that matter.
    print "Excluded_contig", "Sites_dropped" > qc_excl
    for (c in excluded_contig)
        printf "%s\t%d\n", c, excluded_contig[c] >> qc_excl
    close(qc_excl)

    fatal = bad_nf + bad_cov + bad_pct + bad_pct_range + unsorted + interleaved
    if (kept == 0) fatal++
    exit (fatal > 0) ? 1 : 0
}' \
| bgzip -c -@ "${THREADS}" > "${OUT_BED}"

# Capture the WHOLE array in one go, as the very next command. PIPESTATUS is
# rebuilt by every command — including a plain assignment — so reading it as
# RC_A=${PIPESTATUS[0]}; RC_B=${PIPESTATUS[1]} silently loses every element
# after the first.
RCS=("${PIPESTATUS[@]}")
RC_DECOMP=${RCS[0]:-0}
RC_AWK=${RCS[1]:-0}
RC_BGZIP=${RCS[2]:-0}
T1=$(date +%s)
echo ""
echo "Single pass completed in $((T1-T0))s"

# --- Fail loudly, and never leave a partial output behind -------------------
# A truncated .bedmethyl.gz that looks plausible is worse than no file at all:
# every downstream stage would silently analyse an incomplete genome.
if [[ ${RC_DECOMP} -ne 0 || ${RC_AWK} -ne 0 || ${RC_BGZIP} -ne 0 ]]; then
    echo ""
    echo "[FATAL] Extraction failed (decompress=${RC_DECOMP} awk=${RC_AWK} bgzip=${RC_BGZIP})"
    if [[ -f "${QC_MAIN}" ]]; then
        echo "        Assertion counters:"
        grep '^ASSERT_' "${QC_MAIN}" | awk -F'\t' '$2+0 > 0 {printf "          %s = %s\n", $1, $2}'
    fi
    echo "        Removing partial output: ${OUT_BED}"
    rm -f "${OUT_BED}"
    exit 1
fi

# --- Index for region queries ------------------------------------------------
echo -n "Indexing with tabix ... "
if tabix -p bed -f "${OUT_BED}" 2>/dev/null; then
    echo "[OK] ${OUT_BED}.tbi"
else
    echo "[FAIL]"
    echo "[FATAL] tabix indexing failed. The file is coordinate-sorted per the" >&2
    echo "        in-pass assertion, so this usually means bgzip did not produce" >&2
    echo "        a valid BGZF stream." >&2
    exit 1
fi

# --- Report ------------------------------------------------------------------
echo ""
echo "--- Extraction QC (${QC_MAIN}) ---"
awk -F'\t' 'NR>1 && $1 !~ /^ASSERT_/ {printf "  %-42s %s\n", $1, $2}' "${QC_MAIN}"

echo ""
echo "--- Contigs retained (${QC_CONTIG}) ---"
awk -F'\t' 'NR>1 {printf "%s ", $1} END {printf "\n"}' "${QC_CONTIG}" | fold -s -w 76 | sed 's/^/  /'

N_KEPT=$(awk -F'\t' '$1=="Contigs_retained" {print $2}' "${QC_MAIN}")
if [[ "${N_KEPT}" != "${EXPECTED_PRIMARY_CONTIGS}" ]]; then
    echo "  [WARN] ${N_KEPT} contigs retained, expected ${EXPECTED_PRIMARY_CONTIGS}."
    if [[ $(( EXPECTED_PRIMARY_CONTIGS - N_KEPT )) -eq 1 ]]; then
        echo "         One short usually means chrY carried no calls, which is expected"
        echo "         for a female sample. Confirm against the known sex of this"
        echo "         sample rather than assuming."
    elif [[ "${N_KEPT}" -lt "${EXPECTED_PRIMARY_CONTIGS}" ]]; then
        echo "         Substantially fewer contigs than expected. Check whether this"
        echo "         input is a targeted panel rather than whole-genome, and whether"
        echo "         filtering.expected_primary_contigs matches the assay."
    else
        echo "         More contigs than expected — filtering.include_contigs_regex is"
        echo "         admitting something it should not. Review the retained list above."
    fi
else
    echo "  [OK]   ${N_KEPT} contigs retained, as expected"
fi

N_EXCL_CONTIG=$(awk -F'\t' 'NR>1' "${QC_EXCL}" 2>/dev/null | wc -l | tr -d ' ')
if [[ "${N_EXCL_CONTIG:-0}" -gt 0 ]]; then
    N_EXCL_SITE=$(awk -F'\t' 'NR>1 {s+=$2} END {print s+0}' "${QC_EXCL}")
    echo ""
    echo "--- Contigs excluded: ${N_EXCL_CONTIG} contigs, ${N_EXCL_SITE} sites ---"
    echo "  Largest contributors (full list in ${QC_EXCL}):"
    awk -F'\t' 'NR>1' "${QC_EXCL}" | sort -t$'\t' -k2,2nr | head -n 5 \
        | awk -F'\t' '{printf "    %-30s %s sites\n", $1, $2}'
fi

echo ""
echo "--- Assertions ---"
awk -F'\t' '$1 ~ /^ASSERT_/ {printf "  [%s] %-38s %s\n", ($2+0==0 ? "OK  " : "FAIL"), $1, $2}' "${QC_MAIN}"
awk -F'\t' '$1 ~ /^WARN_/   {printf "  [%s] %-38s %s\n", ($2+0==0 ? "OK  " : "WARN"), $1, $2}' "${QC_MAIN}"

echo ""
echo "--- Coverage threshold guide (from ${QC_HIST}) ---"
echo "  Use this to justify filtering.min_coverage in the methods section."
awk -F'\t' '
NR > 1 {
    v = $1; sub(/\+$/, "", v)
    cov[v+0] = $2
    total += $2
    if (v+0 > maxv) maxv = v+0
}
END {
    printf "  %-12s %-14s %s\n", "Threshold", "Sites_retained", "Percent_of_all_sites"
    n = split("1 5 10 15 20 25 30 40 50", t, " ")
    for (i = 1; i <= n; i++) {
        s = 0
        for (c = t[i]; c <= maxv; c++) s += cov[c]
        printf "  >=%-10s %-14d %.2f%%\n", t[i], s, (total ? s/total*100 : 0)
    }
}' "${QC_HIST}"

OUT_SIZE=$(ls -lhL "${OUT_BED}" | awk '{print $5}')
echo ""
echo "Output: ${OUT_BED} (${OUT_SIZE})"
echo "        ${OUT_BED}.tbi"
echo "QC:     ${QC_MAIN}"
echo "        ${QC_HIST}"
echo "        ${QC_CONTIG}"
echo "        ${QC_EXCL}"
echo ""
echo "=== [02_bedmethyl_to_tsv.sh] Done. Log: ${LOG_FILE} ==="
