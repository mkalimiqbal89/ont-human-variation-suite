#!/usr/bin/env bash
# =============================================================================
# 05_annotate_regions.sh
# Intersects the coverage-filtered CpGs with genes, promoters and CpG islands,
# and aggregates each to per-feature methylation.
#
# WHY BASH RATHER THAN R
# The filtered set is ~28 M sites for a whole genome. awk streams it in constant
# memory; reading it into R wants several GB. R is used later, on these small
# per-feature tables, for figures and cross-sample statistics.
#
# AGGREGATION KEYS — both choices are load-bearing, and both were checked
# against the real annotation rather than assumed:
#
#   genes / promoters -> gene_id (column 22)
#       genes.bed has 78,733 rows but only 77,118 distinct gene_name values.
#       Y_RNA alone appears 756 times. Keying on name would merge unrelated loci.
#
#   CpG islands -> chrom:start-end (columns 19-21), NOT the name
#       cpg_islands_hg38.bed has 32,038 islands but only 485 distinct names,
#       because the UCSC "name" field is "CpG:_<cpgNum>". CpG:_21 appears 700
#       times. Keying on name would collapse the genome to 485 rows and every
#       downstream number would be wrong in a way that still looks plausible.
#
# A CpG overlapping two features is counted in BOTH. Genes overlap, so the sum of
# per-feature CpG counts legitimately exceeds the number of unique CpGs, and
# these totals will not reconcile with the genome-wide figures from stage 03.
# Those numbers answer different questions.
#
# Usage:
#   bash scripts/bash/05_annotate_regions.sh [path/to/pipeline_config.yaml]
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

_CFG_ARG="${1:-}"
case "${_CFG_ARG}" in '#'*) _CFG_ARG="" ;; esac
CONFIG_FILE="${_CFG_ARG:-${REPO_DIR}/config/pipeline_config.yaml}"

# shellcheck source=00_setup_env.sh
source "${SCRIPT_DIR}/00_setup_env.sh" "${CONFIG_FILE}" || { echo "[FATAL] env setup failed"; exit 1; }

LOG_FILE="${LOG_DIR}/05_annotate_regions_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

case "${PRIMARY_MOD_CODE}" in
    m) MOD_LABEL="5mC"  ;;
    h) MOD_LABEL="5hmC" ;;
    a) MOD_LABEL="6mA"  ;;
    *) MOD_LABEL="mod_${PRIMARY_MOD_CODE}" ;;
esac

echo ""
echo "=== [05_annotate_regions.sh] Annotating ${MOD_LABEL} sites for ${SAMPLE_ID} ==="

IN_BED="${OUTPUT_DIR}/01_filtered/${SAMPLE_ID}.${MOD_LABEL}.cov${MIN_COVERAGE}.bedmethyl.gz"
if [[ ! -f "${IN_BED}" ]]; then
    echo "[FATAL] Stage 03 output not found: ${IN_BED}" >&2
    echo "        Run: bash scripts/bash/03_filter_cpg_sites.sh" >&2
    exit 1
fi

ANNOT_DIR="${OUTPUT_DIR}/05_annotated_cpgs"
GENE_DIR="${OUTPUT_DIR}/06_gene_summary"
PROM_DIR="${OUTPUT_DIR}/07_promoter_summary"
CGI_DIR="${OUTPUT_DIR}/08_cpg_island_summary"
QC_DIR="${OUTPUT_DIR}/02_qc"
mkdir -p "${ANNOT_DIR}" "${GENE_DIR}" "${PROM_DIR}" "${CGI_DIR}" "${QC_DIR}" "${WORK_DIR}" || {
    echo "[FATAL] Could not create output directories" >&2; exit 1; }

QC_MAIN="${QC_DIR}/${SAMPLE_ID}.05_annotation_qc.tsv"
: > "${QC_MAIN}"
printf "Metric\tValue\n" >> "${QC_MAIN}"

# Total retained sites, for the "what fraction of CpGs fall in a feature" figure.
FILTER_QC="${QC_DIR}/${SAMPLE_ID}.03_filter_qc.tsv"
TOTAL_RETAINED=$(awk -F'\t' '$1=="Retained_sites" {print $2; exit}' "${FILTER_QC}" 2>/dev/null)
[[ -z "${TOTAL_RETAINED}" ]] && TOTAL_RETAINED=0
printf "Total_retained_sites\t%s\n" "${TOTAL_RETAINED}" >> "${QC_MAIN}"
printf "Min_cpgs_per_feature\t%s\n" "${MIN_CPGS_PER_FEATURE}" >> "${QC_MAIN}"

for f in "${GENE_BED}" "${PROMOTER_BED}" "${CPG_ISLAND_BED}"; do
    if [[ ! -f "${f}" ]]; then
        echo "[FATAL] Reference annotation not found: ${f}" >&2
        echo "        Check reference_paths.yaml, then rerun 01_validate_inputs.sh" >&2
        exit 1
    fi
done

if command -v pigz >/dev/null 2>&1; then
    DECOMP=(pigz -dc -p "${THREADS}")
else
    DECOMP=(gzip -dc)
fi

echo "Input         : ${IN_BED}"
echo "Retained sites: ${TOTAL_RETAINED}"
echo "Min CpGs/feat : ${MIN_CPGS_PER_FEATURE} (features below this are flagged, not dropped)"
echo "Keep per-CpG  : ${KEEP_ANNOTATED_CPGS}"
echo ""

# =============================================================================
# awk program for BED7 features (genes, promoters)
# Intersect output: 18 bedMethyl columns + 7 feature columns = 25
#   19 chrom  20 start  21 end  22 gene_id  23 gene_name  24 gene_type  25 strand
# =============================================================================
read -r -d '' AWK_BED7 <<'AWKEOF'
BEGIN {
    OFS = "\t"
    # Writing the per-CpG intermediate from awk, rather than with
    # `tee >(gzip …)`, makes the flush deterministic: awk blocks on close() in
    # END, whereas bash does not reliably wait for a process substitution to
    # finish before the script continues.
    if (annot != "") annot_cmd = "gzip -c > \"" annot "\""
}
{
    rows++
    if (annot != "") print $0 | annot_cmd
    # Unique-CpG counting without an array: bedtools preserves -a order and the
    # input is coordinate-sorted, so rows for one CpG are adjacent. An array
    # keyed on 28 M CpGs would need gigabytes; this needs one variable.
    cur = $1 "\t" $2
    if (cur != prev_cpg) { uniq_cpgs++; prev_cpg = cur }

    key = $22
    if (!(key in seen)) {
        seen[key] = 1; nfeat++
        fname[key]  = $23; ftype[key] = $24
        fchrom[key] = $19; fstart[key] = $20; fend[key] = $21; fstrand[key] = $25
    }
    n[key]++
    cov[key] += $10
    mod[key] += $12
    pct[key] += $11
}
END {
    if (annot != "") close(annot_cmd)
    for (k in seen) {
        sum_n += n[k]
        w = (cov[k] > 0) ? mod[k] / cov[k] * 100 : 0
        u = pct[k] / n[k]
        if (w < 0 || w > 100) bad_w++
        if (u < 0 || u > 100) bad_u++
        if (n[k] >= mincpg) { suff = "yes"; n_suff++ } else suff = "no"
        printf "%s\t%s\t%s\t%s\t%d\t%d\t%s\t%d\t%.0f\t%.0f\t%.4f\t%.4f\t%s\n",
               k, fname[k], ftype[k], fchrom[k], fstart[k], fend[k], fstrand[k],
               n[k], cov[k], mod[k], u, w, suff
        written++
    }
    printf "%s_intersect_rows\t%d\n",          label, rows      >> qc
    printf "%s_unique_cpgs_annotated\t%d\n",   label, uniq_cpgs >> qc
    printf "%s_features_with_data\t%d\n",      label, nfeat     >> qc
    printf "%s_features_written\t%d\n",        label, written    >> qc
    printf "%s_features_min_cpgs_met\t%d\n",   label, n_suff + 0 >> qc
    printf "ASSERT_%s_cpg_sum_vs_rows\t%d\n",  label, rows - sum_n            >> qc
    printf "ASSERT_%s_written_vs_features\t%d\n", label, written - nfeat      >> qc
    printf "ASSERT_%s_weighted_out_of_range\t%d\n", label, bad_w + 0          >> qc
    printf "ASSERT_%s_unweighted_out_of_range\t%d\n", label, bad_u + 0        >> qc
    # Zero overlaps is fatal, and needs its own counter: without it every other
    # assertion reads 0 and the failure report has nothing to show.
    printf "ASSERT_%s_no_overlaps\t%d\n",      label, (rows == 0 ? 1 : 0)     >> qc
    close(qc)
    exit ((rows - sum_n) != 0 || (written - nfeat) != 0 || bad_w > 0 || bad_u > 0 || rows == 0) ? 1 : 0
}
AWKEOF

# =============================================================================
# awk program for the BED10 CpG-island annotation
# Intersect output: 18 + 10 = 28 columns
#   19 chrom  20 start  21 end  22 name  23 length  24 cpgNum  25 gcNum
#   26 perCpg 27 perGc  28 obsExp
# Keyed on chrom:start-end, never on name — see the header note.
# =============================================================================
read -r -d '' AWK_BED10 <<'AWKEOF'
BEGIN {
    OFS = "\t"
    # Writing the per-CpG intermediate from awk, rather than with
    # `tee >(gzip …)`, makes the flush deterministic: awk blocks on close() in
    # END, whereas bash does not reliably wait for a process substitution to
    # finish before the script continues.
    if (annot != "") annot_cmd = "gzip -c > \"" annot "\""
}
{
    rows++
    if (annot != "") print $0 | annot_cmd
    cur = $1 "\t" $2
    if (cur != prev_cpg) { uniq_cpgs++; prev_cpg = cur }

    key = $19 ":" $20 "-" $21
    if (!(key in seen)) {
        seen[key] = 1; nfeat++
        iname[key]  = $22; ilen[key] = $23; icpgnum[key] = $24
        ipercpg[key] = $26; iobsexp[key] = $28
        ichrom[key] = $19; istart[key] = $20; iend[key] = $21
    }
    n[key]++
    cov[key] += $10
    mod[key] += $12
    pct[key] += $11
}
END {
    if (annot != "") close(annot_cmd)
    for (k in seen) {
        sum_n += n[k]
        w = (cov[k] > 0) ? mod[k] / cov[k] * 100 : 0
        u = pct[k] / n[k]
        if (w < 0 || w > 100) bad_w++
        if (u < 0 || u > 100) bad_u++
        if (n[k] >= mincpg) { suff = "yes"; n_suff++ } else suff = "no"
        printf "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%.0f\t%.0f\t%.4f\t%.4f\t%s\n",
               k, iname[k], ichrom[k], istart[k], iend[k], ilen[k],
               icpgnum[k], ipercpg[k], iobsexp[k],
               n[k], cov[k], mod[k], u, w, suff
        written++
    }
    printf "%s_intersect_rows\t%d\n",          label, rows      >> qc
    printf "%s_unique_cpgs_annotated\t%d\n",   label, uniq_cpgs >> qc
    printf "%s_features_with_data\t%d\n",      label, nfeat     >> qc
    printf "%s_features_written\t%d\n",        label, written    >> qc
    printf "%s_features_min_cpgs_met\t%d\n",   label, n_suff + 0 >> qc
    printf "ASSERT_%s_cpg_sum_vs_rows\t%d\n",  label, rows - sum_n            >> qc
    printf "ASSERT_%s_written_vs_features\t%d\n", label, written - nfeat      >> qc
    printf "ASSERT_%s_weighted_out_of_range\t%d\n", label, bad_w + 0          >> qc
    printf "ASSERT_%s_unweighted_out_of_range\t%d\n", label, bad_u + 0        >> qc
    # Zero overlaps is fatal, and needs its own counter: without it every other
    # assertion reads 0 and the failure report has nothing to show.
    printf "ASSERT_%s_no_overlaps\t%d\n",      label, (rows == 0 ? 1 : 0)     >> qc
    close(qc)
    exit ((rows - sum_n) != 0 || (written - nfeat) != 0 || bad_w > 0 || bad_u > 0 || rows == 0) ? 1 : 0
}
AWKEOF

HEADER_BED7=$'Feature_id\tFeature_name\tFeature_type\tChromosome\tFeature_start\tFeature_end\tStrand\tCpG_sites\tTotal_coverage\tModified_read_calls\tMean_site_methylation_percent\tCoverage_weighted_methylation_percent\tMin_cpgs_met'
HEADER_BED10=$'Island_id\tIsland_name\tChromosome\tIsland_start\tIsland_end\tIsland_length\tAnnotated_cpg_count\tPercent_cpg\tObserved_expected_ratio\tCpG_sites\tTotal_coverage\tModified_read_calls\tMean_site_methylation_percent\tCoverage_weighted_methylation_percent\tMin_cpgs_met'

# annotate <label> <feature_bed> <out_tsv> <awk_prog_var> <header> <chrom_col> <start_col>
annotate() {
    local label="$1" bed="$2" out="$3" prog="$4" header="$5" ccol="$6" scol="$7"
    local body="${WORK_DIR}/05_${label}_body.$$.tsv"
    local annot="${ANNOT_DIR}/${SAMPLE_ID}.${label}_annotated_CpGs.tsv.gz"

    echo "--- ${label} (${bed##*/}) ---"
    local t0 t1
    t0=$(date +%s)

    local annot_arg=""
    [[ "${KEEP_ANNOTATED_CPGS}" == "true" ]] && annot_arg="${annot}"

    "${DECOMP[@]}" "${IN_BED}" \
    | bedtools intersect -a stdin -b "${bed}" -wa -wb \
    | awk -F'\t' -v label="${label}" -v qc="${QC_MAIN}" \
          -v mincpg="${MIN_CPGS_PER_FEATURE}" -v annot="${annot_arg}" "${prog}" > "${body}"

    local rcs=("${PIPESTATUS[@]}")
    local bad=0 i
    for i in "${!rcs[@]}"; do [[ "${rcs[$i]}" -ne 0 ]] && bad=1; done
    t1=$(date +%s)

    if [[ ${bad} -ne 0 ]]; then
        echo "  [FATAL] ${label} annotation failed (statuses: ${rcs[*]})"
        grep "^ASSERT_${label}" "${QC_MAIN}" 2>/dev/null \
            | awk -F'\t' '$2+0 != 0 {printf "          %s = %s\n", $1, $2}'
        if [[ "$(awk -F'\t' -v k="ASSERT_${label}_no_overlaps" '$1==k {print $2; exit}' "${QC_MAIN}" 2>/dev/null)" == "1" ]]; then
            echo "          Not a single CpG overlapped a ${label} feature. For a human"
            echo "          genome that is a configuration fault, not biology. Check:"
            echo "            - contig naming: bedMethyl uses '$(gzip -dc "${IN_BED}" 2>/dev/null | head -n1 | cut -f1)',"
            echo "              ${bed##*/} uses '$(head -n1 "${bed}" | cut -f1)'"
            echo "            - that ${bed##*/} is for the same assembly (${GENOME_BUILD})"
            echo "          Re-run 01_validate_inputs.sh, which checks naming concordance."
        fi
        rm -f "${body}" "${out}"
        return 1
    fi

    # Header written separately, body sorted by position: keeps the header first
    # regardless of sort behaviour and makes the file diff-friendly between runs.
    printf '%s\n' "${header}" > "${out}"
    sort -k"${ccol}","${ccol}"V -k"${scol}","${scol}"n "${body}" >> "${out}" || {
        echo "  [FATAL] sort failed for ${label}"; rm -f "${body}"; return 1; }
    rm -f "${body}"

    local nfeat
    nfeat=$(awk -F'\t' 'NR>1' "${out}" | wc -l | tr -d ' ')
    echo "  [OK]   ${nfeat} features -> ${out}  ($((t1-t0))s)"
    [[ "${KEEP_ANNOTATED_CPGS}" == "true" ]] && echo "         per-CpG detail -> ${annot}"
    return 0
}

FAILED=0
annotate "gene"     "${GENE_BED}"       "${GENE_DIR}/${SAMPLE_ID}.gene_methylation_summary.tsv"       "${AWK_BED7}"  "${HEADER_BED7}"  4 5 || FAILED=1
annotate "promoter" "${PROMOTER_BED}"   "${PROM_DIR}/${SAMPLE_ID}.promoter_methylation_summary.tsv"   "${AWK_BED7}"  "${HEADER_BED7}"  4 5 || FAILED=1
annotate "cpg_island" "${CPG_ISLAND_BED}" "${CGI_DIR}/${SAMPLE_ID}.cpg_island_methylation_summary.tsv" "${AWK_BED10}" "${HEADER_BED10}" 3 4 || FAILED=1

if [[ ${FAILED} -ne 0 ]]; then
    echo ""
    echo "=== [05_annotate_regions.sh] FAILED. See log: ${LOG_FILE} ==="
    exit 1
fi

# --- Report ------------------------------------------------------------------
qc() { awk -F'\t' -v k="$1" '$1==k {print $2; exit}' "${QC_MAIN}"; }

echo ""
echo "--- Annotation coverage ---"
printf "  %-14s %-14s %-16s %-14s %s\n" "Feature" "Intersect_rows" "Unique_CpGs" "Pct_of_sites" "Features"
for l in gene promoter cpg_island; do
    ROWS=$(qc "${l}_intersect_rows")
    UNIQ=$(qc "${l}_unique_cpgs_annotated")
    FEAT=$(qc "${l}_features_with_data")
    SUFF=$(qc "${l}_features_min_cpgs_met")
    PCT=$(awk -v a="${UNIQ}" -v b="${TOTAL_RETAINED}" 'BEGIN {printf "%.2f%%", (b ? a/b*100 : 0)}')
    printf "  %-14s %-14s %-16s %-14s %s (%s meet min CpGs)\n" "${l}" "${ROWS}" "${UNIQ}" "${PCT}" "${FEAT}" "${SUFF}"
done
echo "  Intersect_rows exceeds Unique_CpGs because features overlap: a CpG inside"
echo "  two genes is counted in both. That is correct for a per-feature summary."

echo ""
echo "--- Assertions ---"
awk -F'\t' '$1 ~ /^ASSERT_/ {printf "  [%s] %-44s %s\n", ($2+0==0 ? "OK  " : "FAIL"), $1, $2}' "${QC_MAIN}"
echo "  cpg_sum_vs_rows      = intersect rows - sum of per-feature CpG counts"
echo "  written_vs_features  = rows written    - features seen in the intersect"

echo ""
echo "Output: ${GENE_DIR}/${SAMPLE_ID}.gene_methylation_summary.tsv"
echo "        ${PROM_DIR}/${SAMPLE_ID}.promoter_methylation_summary.tsv"
echo "        ${CGI_DIR}/${SAMPLE_ID}.cpg_island_methylation_summary.tsv"
echo "QC:     ${QC_MAIN}"
echo ""
echo "=== [05_annotate_regions.sh] Done. Log: ${LOG_FILE} ==="
