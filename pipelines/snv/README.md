# ONT Epi2ME Small Variant Pipeline

**Status: Released `v1.0.0`**

A reproducible, config-driven downstream pipeline for processing Clair3 Small Variant VCF calls (`*.wf_snp.vcf.gz` and optional `*.wf_snp_clinvar.vcf.gz`) produced by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation).

---

## Input Data Specifications

| Property | Specification |
|---|---|
| Variant Class | Small Variants (SNVs & Indels <50 bp) |
| Upstream Workflow | Epi2ME `wf-human-variation` (Clair3 caller) |
| Primary Input File | `<prefix>.wf_snp.vcf.gz` & optional `<prefix>.wf_snp_clinvar.vcf.gz` |
| Key Information | CHROM, POS, ID, REF, ALT, QUAL, FILTER, GT, DP, VAF, AD, TYPE |

---

## Pipeline Stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `scripts/bash/00_setup_env.sh` | Tool check, version capture, config export, directory setup. **Must be sourced.** |
| 01 | `scripts/bash/01_validate_inputs.sh` | Input VCF file existence, tabix index checks, reference bundle verification |
| 02 | `scripts/bash/02_vcf_to_tsv.sh` | Memory-efficient streaming extraction of raw Clair3 VCF records using `bcftools query` |
| 03 | `scripts/bash/03_filter_snv_categories.sh` | Streaming `awk` filter by DP, QUAL, VAF, and FILTER status. Splits calls into `snvs.tsv` & `indels.tsv`, extracts ClinVar pathogenic variants |
| 04 | `scripts/bash/04_run_all.sh` | Main orchestrator: runs 01→07 end-to-end, stops at first failure, records manifest |
| 05 | `scripts/R/05_annotate_snvs.R` | Intersects SNV/Indel calls with GENCODE gene bed models using `bedtools intersect` |
| 06 | `scripts/R/06_summary_stats.R` | Computes Transition/Transversion (Ti/Tv) ratio, 6-class mutation spectrum, mean depth, and median VAF |
| 07 | `scripts/R/07_generate_report.R` | Generates 4 overview figures (PNGs) across all small variant classes (base R graphics only) |
| 08 | `scripts/bash/08_archive_results.sh` | Longitudinal archiving with SHA-256 integrity verification and index logging |

---

## Requirements

- `bcftools` >= 1.20, `bedtools` >= 2.31, `tabix`/`htslib`
- R >= 4.3 (**base graphics and stats only; zero R packages required**)
- Run dependency check:
  ```bash
  bash scripts/check_dependencies.sh --pipeline snv
  ```

---

## Quick Start

### Option A: Using the Interactive Web UI (Recommended)
Launch the configuration dashboard to set sample IDs, file paths, and reference resources visually:
```bash
python3 scripts/launch_ui.py --pipeline snv
```

### Option B: Using the Command Line
```bash
cd pipelines/snv

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit both config files with sample and reference paths, then:

bash scripts/bash/04_run_all.sh
```

---

## Output File Registry

- `results/snvs/<sample_id>.snvs.tsv`
- `results/indels/<sample_id>.indels.tsv`
- `results/clinvar/<sample_id>.clinvar_pathogenic.tsv`
- `results/gene_snvs/<sample_id>.annotated_snvs.tsv`
- `results/gene_snvs/<sample_id>.gene_snv_summary.tsv`
- `results/qc_summary/<sample_id>.all_snvs_combined.tsv`
- `results/qc_summary/<sample_id>.summary_statistics.tsv`
- `results/qc_summary/figures/<sample_id>.01_snv_indel_counts.png`
- `results/qc_summary/figures/<sample_id>.02_vaf_and_depth_distribution.png`
- `results/qc_summary/figures/<sample_id>.03_titv_and_mutation_spectrum.png`
- `results/qc_summary/figures/<sample_id>.04_chromosome_snv_density.png`

---

## Automated Testing

Run the pipeline test suite against synthetic test fixtures:
```bash
bash pipelines/snv/tests/run_tests.sh
```
