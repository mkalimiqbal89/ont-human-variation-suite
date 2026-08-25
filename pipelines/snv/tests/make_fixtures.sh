#!/usr/bin/env bash
# =============================================================================
# make_fixtures.sh (pipelines/snv/tests)
# Builds synthetic Clair3 SNV VCF fixture, ClinVar VCF fixture, and dummy BED.
# =============================================================================

set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES_DIR="${TEST_DIR}/fixtures"
mkdir -p "${FIXTURES_DIR}"

VCF_TXT="${FIXTURES_DIR}/test_sample.wf_snp.vcf"
cat << 'EOF' > "${VCF_TXT}"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic Depths">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10000	snv1	A	G	50	PASS	TYPE=SNP	GT:DP:AF:AD	0/1:30:0.5000:15,15
chr1	50000	snv2	C	T	45	PASS	TYPE=SNP	GT:DP:AF:AD	0/1:25:0.4400:14,11
chr2	100000	snv3	A	C	60	PASS	TYPE=SNP	GT:DP:AF:AD	0/1:40:0.5000:20,20
chr2	200000	indel1	A	AT	50	PASS	TYPE=INDEL	GT:DP:AF:AD	0/1:35:0.4000:21,14
chr3	1000	indel2	AT	A	30	PASS	TYPE=INDEL	GT:DP:AF:AD	0/1:20:0.3500:13,7
chr3	5000	snv4	C	G	5	LowQual	TYPE=SNP	GT:DP:AF:AD	0/1:5:0.1000:4,1
EOF

bgzip -f "${VCF_TXT}"
tabix -f -p vcf "${VCF_TXT}.gz"

# Synthetic ClinVar VCF
CLIN_TXT="${FIXTURES_DIR}/test_sample.wf_snp_clinvar.vcf"
cat << 'EOF' > "${CLIN_TXT}"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical Significance">
##INFO=<ID=CLNDN,Number=1,Type=String,Description="Disease Name">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	10000	snv1	A	G	50	PASS	CLNSIG=Pathogenic;CLNDN=Cardiomyopathy
EOF

bgzip -f "${CLIN_TXT}"
tabix -f -p vcf "${CLIN_TXT}.gz"

# Dummy gene BED
cat << 'EOF' > "${FIXTURES_DIR}/dummy_genes.bed"
chr1	8000	12000	GENE_X	.	+
chr1	45000	55000	GENE_Y	.	-
chr2	95000	105000	GENE_Z	.	+
EOF

# Test pipeline config
cat << 'EOF' > "${FIXTURES_DIR}/test_pipeline_config.yaml"
sample:
  sample_id: "TEST_SAMPLE"
  raw_sample_prefix: "test_sample"
  run_description: "Synthetic test fixture for SNV pipeline"

paths:
  input_dir: "."
  output_dir: "./results"
  work_dir: "./work"
  repo_dir: "."

reference:
  config_file: "test_reference_paths.yaml"

input_files:
  snv_vcf: "test_sample.wf_snp.vcf.gz"
  snv_vcf_clinvar: "test_sample.wf_snp_clinvar.vcf.gz"

filtering:
  min_dp: 10
  min_qual: 15
  min_vaf: 0.15
  pass_only: true

snv_categories:
  snvs: ["SNV", "SNP"]
  indels: ["INDEL", "INS", "DEL"]

clinvar:
  extract_pathogenic: true

compute:
  threads: 2

logging:
  log_dir: "logs"
  log_level: "INFO"

archive:
  archive_root: ""
  compress: false
EOF

# Test reference paths config
cat << 'EOF' > "${FIXTURES_DIR}/test_reference_paths.yaml"
genome:
  fasta: "dummy_reference.fa"
  fasta_index: "dummy_reference.fa.fai"
  build: "GRCh38"

annotation:
  gene_bed: "dummy_genes.bed"

tools:
  bcftools: "bcftools"
  bedtools: "bedtools"
  tabix: "tabix"
  R: "Rscript"
EOF

echo "Synthetic SNV fixtures created in: ${FIXTURES_DIR}"
