#!/usr/bin/env bash
# =============================================================================
# tests/make_fixtures.sh
# Generates the SYNTHETIC bedMethyl fixtures used by tests/run_tests.sh.
#
# Every fixture is hand-crafted and contains no real patient data. Each one
# isolates exactly one failure mode, so a failing test names the bug.
#
# Usage:
#   bash tests/make_fixtures.sh [output_dir]
# Defaults to tests/fixtures/ next to this script.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURE_DIR="${1:-${SCRIPT_DIR}/fixtures}"
PREFIX="TEST_SAMPLE"

mkdir -p "${FIXTURE_DIR}"

# A valid record. Args: chrom start code coverage percent n_mod n_canonical
rec() {
    printf '%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\t255,0,0\t%s\t%s\t%s\t%s\t0\t0\t0\t0\t0\n' \
        "$1" "$2" "$(($2 + 1))" "$3" "$4" "$2" "$(($2 + 1))" "$4" "$5" "$6" "$7"
}

# Three coordinate-sorted chr1 CpGs, each with both an m and an h record —
# mirroring the interleaved layout modkit actually produces.
base() {
    local i s
    for i in 0 1 2; do
        s=$((1000 + i * 10))
        rec chr1 "${s}" m 20 "75.00" 15 5
        rec chr1 "${s}" h 20 "5.00"  1  19
    done
}

write() { gzip -c > "${FIXTURE_DIR}/${PREFIX}.${1}.wf_mods.bedmethyl.gz"; }

echo "Writing fixtures to ${FIXTURE_DIR}"

# 1. Clean input — must pass.
base | write valid

# 2. A record with 6 columns instead of 18.
{ base; printf 'chr1\t9000\t9001\tm\t20\t.\n'; } | write bad_columns

# 3. A record whose start goes backwards within a contig (breaks tabix).
{ base; rec chr1 500 m 20 "50.00" 10 10; } | write unsorted

# 4. chr1 reappears after chr2 — contigs not in contiguous blocks.
{ base; rec chr2 100 m 20 "50.00" 10 10; rec chr1 9999 m 20 "50.00" 10 10; } | write interleaved_contigs

# 5. percent_modified of 150.
{ base; rec chr1 9000 m 20 "150.00" 10 10; } | write percent_out_of_range

# 6. modified + canonical (30 + 30) exceeds valid coverage (20).
{ base; rec chr1 9000 m 20 "50.00" 30 30; } | write counts_exceed_coverage

# 7. Contigs that filtering.exclude_contigs_regex must drop. The three valid
#    chr1 sites must survive and the four junk contigs must not.
{
    base
    rec chrM               200 m 99 "1.00" 1 98
    rec chr1_random        200 m 99 "1.00" 1 98
    rec chrUn_GL000195v1   200 m 99 "1.00" 1 98
    rec chr1_KI270706v1_alt 200 m 99 "1.00" 1 98
} | write excluded_contigs

# 7b. Scaffolds that the BLOCKLIST regex does not match but that are still not
#     primary chromosomes. This fixture is the entire argument for the allowlist:
#     under exclude_contigs_regex every one of these survives, while
#     include_contigs_regex rejects all of them.
{
    base
    rec GL000191.1            200 m 30 "40.00" 12 18   # no chr prefix at all
    rec "HLA-A*01:01:01:01"   200 m 30 "40.00" 12 18   # ALT haplotype, no _alt suffix
    rec chr1_KI270762v1_fix   200 m 30 "40.00" 12 18   # _fix, not _random/_alt
    rec chr22_KI270879v1_alt  200 m 30 "40.00" 12 18   # _alt: blocklist catches this one
} | write nonstandard_contigs

# 7c. chrY absent — the legitimate female-sample case. Retained contig count
#     comes out at 2 (chr1, chrX) rather than triggering a hard failure.
{
    rec chr1 100 m 20 "50.00" 10 10
    rec chrX 100 m 20 "50.00" 10 10
} | write no_chry

# 8. No records carrying the primary mod code at all — an empty result must be
#    an error, not a silently empty output file.
rec chr1 1000 h 20 "5.00" 1 19 | write no_primary_code

# 8b. Boundary fixture for stage 03. Every value sits exactly on a threshold, so
#     an off-by-one in either direction changes a category count.
#
#     With min_coverage=10, unmethylated_max=20, methylated_min=80:
#       cov 9,  50%    -> dropped (coverage below threshold)
#       cov 10, 50%    -> retained, intermediate  (coverage bound is inclusive)
#       cov 20, 0%     -> unmethylated
#       cov 20, 19.99% -> unmethylated
#       cov 20, 20.00% -> intermediate  (unmethylated is strictly < 20)
#       cov 20, 79.99% -> intermediate
#       cov 20, 80.00% -> methylated    (methylated is >= 80)
#       cov 20, 100%   -> methylated
#     Expected: input 8, dropped_low 1, retained 7, unmeth 2, inter 3, meth 2.
{
    rec chr1 100 m  9 "50.00"    4  5
    rec chr1 200 m 10 "50.00"    5  5
    rec chr1 300 m 20 "0.00"     0 20
    rec chr1 400 m 20 "19.99"    4 16
    rec chr1 500 m 20 "20.00"    4 16
    rec chr1 600 m 20 "79.99"   16  4
    rec chr1 700 m 20 "80.00"   16  4
    rec chr1 800 m 20 "100.00"  20  0
} | write boundaries

# 9. Known-value fixture for arithmetic checks. Coverage-weighted methylation
#    is (10+20+30)/(20+40+60) = 60/120 = exactly 50.0000%, while the unweighted
#    mean of the site percentages is also 50.0000%. Sites differ in coverage so
#    that a bug swapping weighted for unweighted is still caught by the
#    per-site values below.
{
    rec chr1 100 m 20 "50.00" 10 10
    rec chr1 200 m 40 "50.00" 20 20
    rec chr1 300 m 60 "50.00" 30 30
} | write known_values

ls -1 "${FIXTURE_DIR}"/${PREFIX}.*.wf_mods.bedmethyl.gz | while read -r f; do
    printf '  %-28s %s records\n' "$(basename "${f}")" "$(gzip -cd "${f}" | wc -l | tr -d ' ')"
done

# 10. Annotation fixture for stage 05, designed against the two aggregation-key
#     hazards that were measured in the real reference files:
#
#     chr1:1000  -> inside genes ENSGA and ENSGB (which SHARE the gene_name
#                   "SHARED") and inside CpG island 1
#     chr1:5100  -> inside gene ENSGC and inside CpG island 2, whose name is
#                   identical to island 1's ("CpG:_21")
#     chr1:20000 -> intergenic, overlaps nothing
#
#     Correct behaviour: 3 gene features (not 2 — gene_name must not merge
#     ENSGA/ENSGB) and 2 island features (not 1 — the island name must not merge
#     two distinct loci). In the real annotation Y_RNA repeats 756 times and
#     CpG:_21 repeats 700 times, so both collisions are the normal case.
{
    rec chr1  1000 m 20 "75.00" 15  5
    rec chr1  5100 m 20 "40.00"  8 12
    rec chr1 20000 m 20 "10.00"  2 18
} | write annotation

# 10b. Second sample for the cross-sample comparison, at the SAME coordinates as
#      `annotation` but with deliberately different methylation, so the expected
#      answer is arithmetic rather than a judgement call:
#
#        chr1:1000   75% -> 10%   delta = -65 pp   (strongly hypomethylated)
#        chr1:5100   40% -> 45%   delta =  +5 pp   (below the 20 pp threshold)
#        chr1:20000  10% -> 10%   delta =   0 pp   (unchanged; intergenic anyway)
#
#      Coverage is raised to 40 so the pooled counts are large enough for a
#      Fisher test to return something meaningful.
{
    rec chr1  1000 m 40 "10.00"  4 36
    rec chr1  5100 m 40 "45.00" 18 22
    rec chr1 20000 m 40 "10.00"  4 36
} | write annotation_b

# 10c. Third sample, deliberately COVERAGE-IMBALANCED against `annotation`.
#      Gene ENSGA.1 spans chr1:900-1100. This sample places 6 covered CpGs inside
#      it where `annotation` has 1, so cpg_ratio = 6 (above the default
#      max_cpg_ratio of 3) and coverage_ratio = 6 * 30 / 20 (above nothing on its
#      own, but the CpG ratio alone must be enough to disqualify).
#
#      Methylation is set to 5% against `annotation`'s 75%, a 70 pp difference
#      that would otherwise be called differential. It must come out labelled
#      "imbalanced" instead — this is the mappability artefact that topped the
#      first real HLH_S0001 vs HLH_S0002 comparison.
{
    rec chr1  1000 m 30 "5.00" 2 28
    rec chr1  1010 m 30 "5.00" 2 28
    rec chr1  1020 m 30 "5.00" 2 28
    rec chr1  1030 m 30 "5.00" 2 28
    rec chr1  1040 m 30 "5.00" 2 28
    rec chr1  1050 m 30 "5.00" 2 28
    rec chr1  5100 m 30 "40.00" 12 18
    rec chr1 20000 m 30 "10.00" 3 27
} | write annotation_imbalanced

# --- Reference annotation for tests -----------------------------------------
# Absolute paths: stage 05 resolves these directly, and the test harness runs
# from a scratch directory where a repo-relative path would not exist.
printf 'chr1\t250000000\nchr2\t250000000\n' > "${FIXTURE_DIR}/test_chrom.sizes"

# BED7. ENSGA and ENSGB deliberately share a gene_name and overlap chr1:1000.
{
    printf 'chr1\t900\t1100\tENSGA.1\tSHARED\tprotein_coding\t+\n'
    printf 'chr1\t950\t1200\tENSGB.1\tSHARED\tlncRNA\t-\n'
    printf 'chr1\t5000\t5200\tENSGC.1\tSOLO\tprotein_coding\t+\n'
} > "${FIXTURE_DIR}/test_genes.bed"

{
    printf 'chr1\t900\t1100\tENSGA.1\tSHARED\tprotein_coding\t+\n'
    printf 'chr1\t5000\t5200\tENSGC.1\tSOLO\tprotein_coding\t+\n'
} > "${FIXTURE_DIR}/test_promoters.bed"

# BED10. Both islands carry the name "CpG:_21" at different loci.
{
    printf 'chr1\t900\t1100\tCpG:_21\t200\t21\t120\t15.0\t60.0\t0.70\n'
    printf 'chr1\t5000\t5200\tCpG:_21\t200\t21\t120\t15.0\t60.0\t0.70\n'
} > "${FIXTURE_DIR}/test_cpg_islands.bed"

cat > "${FIXTURE_DIR}/test_reference_paths.yaml" <<EOF
# Synthetic reference config for tests. Generated by make_fixtures.sh — paths are
# absolute because the harness runs stages from scratch directories.
genome:
  fasta: ""
  build: "TEST"
  chrom_sizes: "${FIXTURE_DIR}/test_chrom.sizes"
annotation:
  gene_bed: "${FIXTURE_DIR}/test_genes.bed"
  promoter_bed: "${FIXTURE_DIR}/test_promoters.bed"
  cpg_island_bed: "${FIXTURE_DIR}/test_cpg_islands.bed"
  gencode_gtf: ""
tools:
  bedtools: "bedtools"
  tabix: "tabix"
  bgzip: "bgzip"
  R: "Rscript"
EOF

echo "Done."
