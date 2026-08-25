# pipelines/cnv — Copy Number Variation

**Status: Released `v1.0.0`**

Reproducible, config-driven downstream pipeline for processing Spectre Copy Number Variation VCF calls (`*.wf_cnv.vcf.gz`) produced by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation).

---

## Stages

- `00_setup_env.sh`: Resolves environment variables, config fields, and validates tool binaries (`bcftools`, `bedtools`, `tabix`, `Rscript`).
- `01_validate_inputs.sh`: Validates input `*.wf_cnv.vcf.gz` and reference BED resources.
- `02_vcf_to_tsv.sh`: Extracts and flattens Spectre VCF records into `data/processed/<sample_id>.cnv_flat.tsv`.
- `03_filter_cnv_categories.sh`: Applies length, quality, FILTER status, and copy-number filters. Splits into deletions and duplications/gains tables.
- `04_run_all.sh`: Orchestrates the complete pipeline end-to-end with fail-fast execution and runtime tracking.
- `05_annotate_cnvs.R`: Intersects CNV regions with GENCODE gene BED models using `bedtools intersect`.
- `06_summary_stats.R`: Computes summary statistics (counts, mean/median lengths, total genomic span affected).
- `07_generate_report.R`: Generates publication-grade base R figures (zero external R package dependencies).
- `08_archive_results.sh`: Archives run outputs with SHA-256 integrity checksums and longitudinal provenance records.

---

## Quick Start

```bash
cd pipelines/cnv

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit both config files with sample and reference paths, then:

bash scripts/bash/04_run_all.sh
```

Run automated tests:

```bash
bash tests/run_tests.sh
```

---

## Outputs

- `results/deletions/<sample_id>.deletions.tsv`
- `results/duplications/<sample_id>.duplications.tsv`
- `results/gene_cnvs/<sample_id>.annotated_cnvs.tsv`
- `results/gene_cnvs/<sample_id>.gene_cnv_summary.tsv`
- `results/qc_summary/<sample_id>.all_cnvs_combined.tsv`
- `results/qc_summary/<sample_id>.summary_statistics.tsv`
- `results/qc_summary/figures/<sample_id>.01_cnv_counts_by_category.png`
- `results/qc_summary/figures/<sample_id>.02_cnv_size_distribution.png`
- `results/qc_summary/figures/<sample_id>.03_chromosome_cnv_distribution.png`
- `results/qc_summary/figures/<sample_id>.04_copy_number_distribution.png`
