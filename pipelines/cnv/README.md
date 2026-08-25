# ONT Epi2ME Copy Number Variation Pipeline

**Status: Released `v1.0.0`**

A reproducible, config-driven downstream pipeline for processing Spectre Copy Number Variation VCF calls (`*.wf_cnv.vcf.gz`) produced by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation).

---

## Input Data Specifications

| Property | Specification |
|---|---|
| Variant Class | Copy Number Variants (CNVs, Gains & Losses) |
| Upstream Workflow | Epi2ME `wf-human-variation` (Spectre caller) |
| Primary Input File | `<prefix>.wf_cnv.vcf.gz` (Spectre VCF) |
| Key Information | CHROM, POS, END, SVLEN, SVTYPE, CN (Copy Number), GT, QUAL, FILTER |

---

## Pipeline Stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `scripts/bash/00_setup_env.sh` | Tool check, version capture, config export, directory setup. **Must be sourced.** |
| 01 | `scripts/bash/01_validate_inputs.sh` | Input VCF file existence, tabix index checks, reference bundle verification |
| 02 | `scripts/bash/02_vcf_to_tsv.sh` | Flattens Spectre VCF into an analysis-ready TSV using `bcftools query` |
| 03 | `scripts/bash/03_filter_cnv_categories.sh` | Applies QC thresholds (length, QUAL, copy number) and splits calls into deletions/losses & duplications/gains |
| 04 | `scripts/bash/04_run_all.sh` | Main orchestrator: runs 01→07 end-to-end, stops at first failure, records manifest |
| 05 | `scripts/R/05_annotate_cnvs.R` | Intersects CNV calls with GENCODE gene bed models using `bedtools intersect` |
| 06 | `scripts/R/06_summary_stats.R` | Computes per-category CNV counts, mean/median segment lengths, and total genomic span affected |
| 07 | `scripts/R/07_generate_report.R` | Generates 4 overview figures (PNGs) across all CNV classes (base R graphics only) |
| 08 | `scripts/bash/08_archive_results.sh` | Longitudinal archiving with SHA-256 integrity verification and index logging |

---

## Requirements

- `bcftools` >= 1.20, `bedtools` >= 2.31, `tabix`/`htslib`
- R >= 4.3 (**base graphics and stats only; zero R packages required**)
- Run dependency check:
  ```bash
  bash scripts/check_dependencies.sh --pipeline cnv
  ```

---

## Quick Start

### Option A: Using the Interactive Web UI (Recommended)
Launch the configuration dashboard to set sample IDs, file paths, and reference resources visually:
```bash
python3 scripts/launch_ui.py --pipeline cnv
```

### Option B: Using the Command Line
```bash
cd pipelines/cnv

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit both config files with sample and reference paths, then:

bash scripts/bash/04_run_all.sh
```

---

## Output File Registry

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

---

## Automated Testing

Run the pipeline test suite against synthetic test fixtures:
```bash
bash pipelines/cnv/tests/run_tests.sh
```
