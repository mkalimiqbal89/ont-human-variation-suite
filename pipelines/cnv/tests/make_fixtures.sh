#!/usr/bin/env bash
# =============================================================================
# make_fixtures.sh (pipelines/cnv/tests)
# Builds synthetic Spectre CNV VCF test fixture and reference files.
# =============================================================================

set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES_DIR="${TEST_DIR}/fixtures"
mkdir -p "${FIXTURES_DIR}"

VCF_TXT="${FIXTURES_DIR}/test_sample.wf_cnv.vcf"
cat << 'EOF' > "${VCF_TXT}"
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">
##FORMAT=<ID=FC,Number=1,Type=Float,Description="Fold change">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10000	cnv1	N	<DEL>	60	PASS	SVTYPE=DEL;END=25000;SVLEN=15000;CN=1	GT:CN:FC	0/1:1:-1.0
chr1	50000	cnv2	N	<DUP>	50	PASS	SVTYPE=DUP;END=75000;SVLEN=25000;CN=3	GT:CN:FC	0/1:3:1.5
chr2	100000	cnv3	N	<DEL>	55	PASS	SVTYPE=DEL;END=120000;SVLEN=20000;CN=0	GT:CN:FC	1/1:0:-2.0
chr2	200000	cnv4	N	<DUP>	10	LowQual	SVTYPE=DUP;END=210000;SVLEN=10000;CN=4	GT:CN:FC	0/1:4:2.0
chr3	5000	cnv5	N	<DEL>	40	PASS	SVTYPE=DEL;END=6000;SVLEN=1000;CN=1	GT:CN:FC	0/1:1:-1.0
EOF

bgzip -f "${VCF_TXT}"
tabix -f -p vcf "${VCF_TXT}.gz"

# Create dummy gene bed
cat << 'EOF' > "${FIXTURES_DIR}/dummy_genes.bed"
chr1	12000	18000	GENE_A	.	+
chr1	60000	70000	GENE_B	.	+
chr2	105000	115000	GENE_C	.	-
EOF

# Create test pipeline config
cat << 'EOF' > "${FIXTURES_DIR}/test_pipeline_config.yaml"
sample:
  sample_id: "TEST_SAMPLE"
  raw_sample_prefix: "test_sample"
  run_description: "Synthetic test fixture for CNV pipeline"

paths:
  input_dir: "."
  output_dir: "./results"
  work_dir: "./work"
  repo_dir: "."

reference:
  config_file: "test_reference_paths.yaml"

input_files:
  cnv_vcf: "test_sample.wf_cnv.vcf.gz"

filtering:
  min_cnv_length: 10000
  min_qual: 20
  pass_only: true
  max_cn_loss: 1
  min_cn_gain: 3

cnv_categories:
  deletions: ["DEL", "LOSS"]
  duplications: ["DUP", "GAIN", "AMP"]

compute:
  threads: 2

logging:
  log_dir: "logs"
  log_level: "INFO"

archive:
  archive_root: ""
  compress: false
EOF

# Create test reference paths config
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

echo "Synthetic CNV fixtures created in: ${FIXTURES_DIR}"
