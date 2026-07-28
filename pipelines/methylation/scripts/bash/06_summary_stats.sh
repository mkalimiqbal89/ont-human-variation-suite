#!/usr/bin/env bash
# =============================================================================
# 06_summary_stats.sh
# Produces the genome-wide summary tables: QC, global methylation, the
# methylation distribution, the per-chromosome breakdown, and a feature-class
# comparison drawn from stage 05's per-feature tables.
#
# Every output is a plain TSV. Stage 07 (R) reads these small tables to draw
# figures and never touches the ~28 M-site data file.
#
# Bash rather than R for the same two reasons as stage 05: the demo walkthrough's
# principle that each step yields an inspectable TSV, and the fact that awk
# streams 28 M records in constant memory.
#
# Cross-stage check: the coverage-weighted methylation recomputed here must equal
# the value stage 03 recorded. Both read the same file, so any disagreement means
# the file changed between stages.
#
# Usage:
#   bash scripts/bash/06_summary_stats.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in '#'*) _CFG_ARG="" ;; esac
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/06_summary_stats_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

case "${PRIMARY_MOD_CODE}" in
    m) MOD_LABEL="5mC"  ;;
    h) MOD_LABEL="5hmC" ;;
    a) MOD_LABEL="6mA"  ;;
    *) MOD_LABEL="mod_${PRIMARY_MOD_CODE}" ;;
esac

echo ""
echo "=== [06_summary_stats.sh] Summarising ${MOD_LABEL} for ${SAMPLE_ID} ==="

IN_BED="${OUTPUT_DIR}/01_filtered/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.bedmethyl.gz"
if [[ ! -f "${IN_BED}" ]]; then
    echo "[FATAL] Stage 03 output not found: ${IN_BED}" >&2
    echo "        Run: bash scripts/bash/03_filter_cpg_sites.sh" >&2
    exit 1
fi

QC_DIR="${OUTPUT_DIR}/02_qc"
GLOBAL_DIR="${OUTPUT_DIR}/03_global_methylation"
CHROM_DIR="${OUTPUT_DIR}/04_chromosome_summary"
mkdir -p "${QC_DIR}" "${GLOBAL_DIR}" "${CHROM_DIR}" "${WORK_DIR}" || {
    echo "[FATAL] Could not create output directories" >&2; exit 1; }

QC_MAIN="${QC_DIR}/${SAMPLE_ID}.06_summary_qc.tsv"
OUT_GLOBAL="${GLOBAL_DIR}/${SAMPLE_ID}.global_methylation_summary.tsv"
OUT_DIST="${GLOBAL_DIR}/${SAMPLE_ID}.methylation_distribution.tsv"
OUT_CHROM="${CHROM_DIR}/${SAMPLE_ID}.chromosome_methylation_summary.tsv"
OUT_CLASS="${GLOBAL_DIR}/${SAMPLE_ID}.feature_class_summary.tsv"
OUT_TYPE="${GLOBAL_DIR}/${SAMPLE_ID}.methylation_by_gene_type.tsv"

if command -v pigz >/dev/null 2>&1; then
    DECOMP=(pigz -dc -p "${THREADS}")
else
    DECOMP=(gzip -dc)
fi

echo "Input : ${IN_BED}"
echo ""

CHROM_BODY="${WORK_DIR}/06_chrom_body.$$.tsv"
T0=$(date +%s)

# --- Single pass: global, distribution, per-chromosome ----------------------
"${DECOMP[@]}" "${IN_BED}" \
| awk -F'\t' \
      -v qc="${QC_MAIN}" \
      -v f_global="${OUT_GLOBAL}" \
      -v f_dist="${OUT_DIST}" \
      -v sample="${SAMPLE_ID}" \
      -v modlabel="${MOD_LABEL}" '
BEGIN { OFS = "\t" }
{
    sites++
    chrom = $1
    cov   = $10 + 0
    pct   = $11 + 0
    nmod  = $12 + 0

    total_cov += cov
    total_mod += nmod
    pct_sum   += pct

    if (sites == 1 || cov < min_cov) min_cov = cov
    if (sites == 1 || cov > max_cov) max_cov = cov

    # Per-chromosome. Contigs are contiguous in the input (asserted in stage 02),
    # so first-seen order is genomic order.
    if (!(chrom in seen)) { seen[chrom] = 1; n_chrom++; order[n_chrom] = chrom }
    c_sites[chrom]++
    c_cov[chrom]   += cov
    c_mod[chrom]   += nmod
    c_pct[chrom]   += pct

    # Distribution in 10% bins. 100% folds into the top bin, which is therefore
    # labelled 90-100 rather than 90-99.
    b = int(pct / 10) * 10
    if (b >= 100) b = 90
    hist[b]++
}
END {
    # ---- global ----
    print "Metric", "Value" > f_global
    printf "Sample\t%s\n",                                sample   >> f_global
    printf "Modification\t%s\n",                          modlabel >> f_global
    printf "Retained_CpG_sites\t%d\n",                    sites     >> f_global
    printf "Total_valid_coverage\t%.0f\n",                total_cov >> f_global
    printf "Total_modified_read_calls\t%.0f\n",           total_mod >> f_global
    printf "Minimum_coverage\t%d\n",                      min_cov   >> f_global
    printf "Maximum_coverage\t%d\n",                      max_cov   >> f_global
    printf "Mean_coverage\t%.4f\n",                       (sites ? total_cov / sites : 0) >> f_global
    printf "Unweighted_mean_site_methylation_percent\t%.4f\n",
           (sites ? pct_sum / sites : 0) >> f_global
    printf "Coverage_weighted_methylation_percent\t%.4f\n",
           (total_cov > 0 ? total_mod / total_cov * 100 : 0) >> f_global
    close(f_global)

    # ---- distribution ----
    print "Bin", "Sites", "Percent_of_sites" > f_dist
    for (i = 0; i <= 90; i += 10) {
        lab = (i == 90) ? "90-100" : (i "-" (i + 9))
        printf "%s\t%d\t%.4f\n", lab, hist[i] + 0,
               (sites ? hist[i] / sites * 100 : 0) >> f_dist
        dist_total += hist[i] + 0
    }
    close(f_dist)

    # ---- per-chromosome body, to stdout for sorting in the shell ----
    for (i = 1; i <= n_chrom; i++) {
        c = order[i]
        printf "%s\t%d\t%.0f\t%.0f\t%.4f\t%.4f\n", c, c_sites[c], c_cov[c], c_mod[c],
               c_pct[c] / c_sites[c],
               (c_cov[c] > 0 ? c_mod[c] / c_cov[c] * 100 : 0)
        chrom_site_total += c_sites[c]
    }

    # ---- QC ----
    print "Metric", "Value" > qc
    printf "Retained_CpG_sites\t%d\n",   sites     >> qc
    printf "Chromosomes\t%d\n",          n_chrom   >> qc
    printf "Coverage_weighted_methylation_percent\t%.4f\n",
           (total_cov > 0 ? total_mod / total_cov * 100 : 0) >> qc
    printf "ASSERT_distribution_sums_to_sites\t%d\n", sites - dist_total       >> qc
    printf "ASSERT_chromosome_sums_to_sites\t%d\n",   sites - chrom_site_total >> qc
    printf "ASSERT_no_sites\t%d\n",                  (sites == 0 ? 1 : 0)      >> qc
    close(qc)

    exit ((sites - dist_total) != 0 || (sites - chrom_site_total) != 0 || sites == 0) ? 1 : 0
}' > "${CHROM_BODY}"

RCS=("${PIPESTATUS[@]}")
T1=$(date +%s)
BAD=0
for i in "${!RCS[@]}"; do [[ "${RCS[$i]}" -ne 0 ]] && BAD=1; done

if [[ ${BAD} -ne 0 ]]; then
    echo "[FATAL] Summary pass failed (statuses: ${RCS[*]})"
    [[ -f "${QC_MAIN}" ]] && grep '^ASSERT_' "${QC_MAIN}" \
        | awk -F'\t' '$2+0 != 0 {printf "          %s = %s\n", $1, $2}'
    rm -f "${CHROM_BODY}" "${OUT_GLOBAL}" "${OUT_DIST}" "${OUT_CHROM}"
    exit 1
fi
echo "Single pass completed in $((T1-T0))s"

printf 'Chromosome\tCpG_sites\tTotal_coverage\tModified_read_calls\tMean_site_methylation_percent\tCoverage_weighted_methylation_percent\n' > "${OUT_CHROM}"
sort -k1,1V "${CHROM_BODY}" >> "${OUT_CHROM}"
rm -f "${CHROM_BODY}"

# --- Cross-stage consistency -------------------------------------------------
# Stage 03 and stage 06 read the same file, so these must agree exactly.
FILTER_QC="${QC_DIR}/${SAMPLE_ID}.03_filter_qc.tsv"
W06=$(awk -F'\t' '$1=="Coverage_weighted_methylation_percent" {print $2; exit}' "${QC_MAIN}")
if [[ -f "${FILTER_QC}" ]]; then
    W03=$(awk -F'\t' '$1=="Coverage_weighted_methylation_percent" {print $2; exit}' "${FILTER_QC}")
    S03=$(awk -F'\t' '$1=="Retained_sites" {print $2; exit}' "${FILTER_QC}")
    S06=$(awk -F'\t' '$1=="Retained_CpG_sites" {print $2; exit}' "${QC_MAIN}")
    MATCH=1
    [[ "${W03}" != "${W06}" ]] && MATCH=0
    [[ "${S03}" != "${S06}" ]] && MATCH=0
    printf "ASSERT_matches_stage03\t%d\n" "$((1 - MATCH))" >> "${QC_MAIN}"
    if [[ ${MATCH} -eq 1 ]]; then
        echo "Cross-stage check: agrees with stage 03 (${S06} sites, ${W06}%)"
    else
        echo "[FATAL] Disagrees with stage 03:" >&2
        echo "        stage 03: sites=${S03} weighted=${W03}" >&2
        echo "        stage 06: sites=${S06} weighted=${W06}" >&2
        echo "        Both read the same file, so it changed between stages." >&2
        exit 1
    fi
else
    echo "[WARN] ${FILTER_QC} not found — cross-stage check skipped"
    printf "ASSERT_matches_stage03\t%s\n" "NA" >> "${QC_MAIN}"
fi

# --- Feature-class comparison ------------------------------------------------
# Drawn from stage 05's per-feature tables. This is the headline biological
# comparison: promoters and CpG islands against gene bodies and the genome mean.
#
# Only features meeting annotation.min_cpgs_per_feature are included. A locus
# backed by one or two CpGs can read 0% or 100% by chance, and including those
# would widen every distribution with noise rather than signal.
GENE_TSV="${OUTPUT_DIR}/06_gene_summary/${SAMPLE_ID}.gene_methylation_summary.tsv"
PROM_TSV="${OUTPUT_DIR}/07_promoter_summary/${SAMPLE_ID}.promoter_methylation_summary.tsv"
CGI_TSV="${OUTPUT_DIR}/08_cpg_island_summary/${SAMPLE_ID}.cpg_island_methylation_summary.tsv"

printf 'Feature_class\tFeatures\tMean_weighted_methylation_percent\tMedian_weighted_methylation_percent\tQ1\tQ3\n' > "${OUT_CLASS}"

# class_row <label> <tsv> <weighted_col> <flag_col> [type_col] [type_value]
# Extract the weighted-methylation column for qualifying features, sort it
# numerically, then take mean and quartiles by index. Sorting externally keeps
# the awk POSIX-safe: BSD awk on macOS has no asort().
#
# The optional type filter restricts to one gene_type. This matters because
# GENCODE's ~78 k "genes" are mostly lncRNAs, pseudogenes and small RNAs, and
# only ~20 k are protein-coding. Pooling them hides the promoter hypomethylation
# that is the whole point of looking at promoters.
class_row() {
    local label="$1" tsv="$2" wcol="$3" fcol="$4" tcol="${5:-}" tval="${6:-}"
    if [[ ! -f "${tsv}" ]]; then
        echo "  [WARN] ${label}: ${tsv##*/} not found — run stage 05 first; row skipped"
        return
    fi
    local tmp="${WORK_DIR}/06_class_${label//[^A-Za-z0-9_]/_}.$$.tmp"
    if [[ -n "${tcol}" ]]; then
        awk -F'\t' -v w="${wcol}" -v f="${fcol}" -v tc="${tcol}" -v tv="${tval}" \
            'NR>1 && $f=="yes" && $tc==tv {print $w}' "${tsv}" | sort -g > "${tmp}"
    else
        awk -F'\t' -v w="${wcol}" -v f="${fcol}" \
            'NR>1 && $f=="yes" {print $w}' "${tsv}" | sort -g > "${tmp}"
    fi
    awk -v label="${label}" '
        { v[NR] = $1 + 0; sum += $1 }
        END {
            if (NR == 0) { printf "%s\t0\tNA\tNA\tNA\tNA\n", label; exit }
            i1 = int((NR + 1) / 4);       if (i1 < 1)  i1 = 1
            i2 = int((NR + 1) / 2);       if (i2 < 1)  i2 = 1
            i3 = int(3 * (NR + 1) / 4);   if (i3 < 1)  i3 = 1
            if (i3 > NR) i3 = NR
            printf "%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
                   label, NR, sum / NR, v[i2], v[i1], v[i3]
        }' "${tmp}" >> "${OUT_CLASS}"
    rm -f "${tmp}"
}

# by_gene_type <class_label> <tsv> <weighted_col> <flag_col> <type_col>
# Full breakdown of one feature class across every gene_type present.
#
# Values are sorted by (type, value) first, so the groups arrive contiguous and
# already ordered. That lets a single awk pass compute count, mean and quartiles
# per group without needing asort() or holding the whole table in memory.
by_gene_type() {
    local cls="$1" tsv="$2" wcol="$3" fcol="$4" tcol="$5"
    [[ -f "${tsv}" ]] || return
    awk -F'\t' -v w="${wcol}" -v f="${fcol}" -v tc="${tcol}" \
        'NR>1 && $f=="yes" {print $tc "\t" $w}' "${tsv}" \
    | sort -t$'\t' -k1,1 -k2,2g \
    | awk -F'\t' -v cls="${cls}" '
        function emit() {
            if (n == 0) return
            i1 = int((n + 1) / 4);     if (i1 < 1) i1 = 1
            i2 = int((n + 1) / 2);     if (i2 < 1) i2 = 1
            i3 = int(3 * (n + 1) / 4); if (i3 < 1) i3 = 1
            if (i3 > n) i3 = n
            printf "%s\t%s\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n",
                   cls, cur, n, sum / n, v[i2], v[i1], v[i3]
        }
        {
            if ($1 != cur) { emit(); cur = $1; n = 0; sum = 0 }
            n++; v[n] = $2 + 0; sum += $2
        }
        END { emit() }' >> "${OUT_TYPE}"
}

# Genome-wide reference line: the coverage-weighted figure over all sites.
printf 'genome_wide_all_sites\t%s\t%s\tNA\tNA\tNA\n' \
       "$(awk -F'\t' '$1=="Retained_CpG_sites" {print $2}' "${QC_MAIN}")" "${W06}" >> "${OUT_CLASS}"
class_row "gene_body"   "${GENE_TSV}" 12 13
class_row "promoter"    "${PROM_TSV}" 12 13
class_row "cpg_island"  "${CGI_TSV}"  14 15

# Protein-coding subsets, promoted into the main table because they are the rows
# most likely to be quoted. Column 3 of the gene/promoter summaries is gene_type.
class_row "gene_body_protein_coding" "${GENE_TSV}" 12 13 3 "protein_coding"
class_row "promoter_protein_coding"  "${PROM_TSV}" 12 13 3 "protein_coding"

# --- Full gene_type breakdown ------------------------------------------------
printf 'Feature_class\tGene_type\tFeatures\tMean_weighted_methylation_percent\tMedian_weighted_methylation_percent\tQ1\tQ3\n' > "${OUT_TYPE}"
by_gene_type "gene_body" "${GENE_TSV}" 12 13 3
by_gene_type "promoter"  "${PROM_TSV}" 12 13 3

# --- Report ------------------------------------------------------------------
echo ""
echo "--- Global methylation (${OUT_GLOBAL##*/}) ---"
awk -F'\t' 'NR>1 {printf "  %-42s %s\n", $1, $2}' "${OUT_GLOBAL}"

echo ""
echo "--- Distribution (${OUT_DIST##*/}) ---"
awk -F'\t' 'NR==1 {printf "  %-10s %-14s %s\n", $1, $2, $3; next}
            {printf "  %-10s %-14s %s%%\n", $1, $2, $3}' "${OUT_DIST}"

NCHR=$(( $(wc -l < "${OUT_CHROM}" | tr -d ' ') - 1 ))
echo ""
if [[ ${NCHR} -le 14 ]]; then
    echo "--- Per-chromosome (${NCHR} chromosomes) ---"
    cut -f1,2,3,5,6 "${OUT_CHROM}" | column -t -s $'\t' | sed 's/^/  /'
else
    echo "--- Per-chromosome (${NCHR} chromosomes; first 6 and last 3 shown) ---"
    { head -n 1 "${OUT_CHROM}"
      awk -F'\t' 'NR>1' "${OUT_CHROM}" | head -n 6
      printf '...\t\t\t\t\n'
      awk -F'\t' 'NR>1' "${OUT_CHROM}" | tail -n 3
    } | cut -f1,2,3,5,6 | column -t -s $'\t' | sed 's/^/  /'
    echo "  (full table in ${OUT_CHROM##*/})"
fi

echo ""
echo "--- Feature classes (${OUT_CLASS##*/}) ---"
column -t -s $'\t' "${OUT_CLASS}" | sed 's/^/  /'
echo "  Mean/median are across FEATURES, each weighted internally by coverage."
echo "  Only features meeting min_cpgs_per_feature=${MIN_CPGS_PER_FEATURE} are included."

echo ""
echo "--- By gene_type: 8 most abundant types (${OUT_TYPE##*/}) ---"
{ head -n 1 "${OUT_TYPE}"
  awk -F'\t' 'NR>1 && $1=="promoter"' "${OUT_TYPE}" | sort -t$'\t' -k3,3nr | head -n 8
} | cut -f2,3,4,5,6,7 | column -t -s $'\t' | sed 's/^/  /'
echo "  (promoters shown; gene bodies also in the file)"
echo "  GENCODE's ~78k genes are mostly non-coding, so the pooled promoter row"
echo "  above understates protein-coding promoter hypomethylation."

echo ""
echo "--- Assertions ---"
awk -F'\t' '$1 ~ /^ASSERT_/ {printf "  [%s] %-38s %s\n",
            ($2=="NA" ? "SKIP" : ($2+0==0 ? "OK  " : "FAIL")), $1, $2}' "${QC_MAIN}"

echo ""
echo "Output: ${OUT_GLOBAL}"
echo "        ${OUT_DIST}"
echo "        ${OUT_CHROM}"
echo "        ${OUT_CLASS}"
echo "        ${OUT_TYPE}"
echo "QC:     ${QC_MAIN}"
echo ""
echo "=== [06_summary_stats.sh] Done. Log: ${LOG_FILE} ==="
