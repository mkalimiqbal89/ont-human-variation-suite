# ONT Epi2ME Methylation Pipeline

**Status: Feature-complete (in validation)**

A reproducible, config-driven pipeline for genome-wide CpG methylation analysis from Oxford Nanopore long-read data processed by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation).

---

## Input Data Specifications

| Property | Specification |
|---|---|
| Variant Class | CpG Modification (5mC / 5hmC) |
| Upstream Workflow | Epi2ME `wf-human-variation` (`modkit pileup`) |
| Primary Input File | `<prefix>.wf_mods.bedmethyl.gz` (Unphased) or `.1.bedmethyl.gz` (Phased) |
| Key Information | 18-column bedMethyl layout (chrom, start, end, strand, coverage, % modified, mod count, canonical count) |

---

## Pipeline Stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `scripts/bash/00_setup_env.sh` | Tool check, version capture, config export, directory setup. **Must be sourced.** |
| 01 | `scripts/bash/01_validate_inputs.sh` | Input file existence, gzip integrity, column layout, mod codes, contig concordance |
| 02 | `scripts/bash/02_bedmethyl_to_tsv.sh` | Single-pass extraction of primary mod code; full-file structural assertions; coverage histogram; bgzip + tabix |
| 03 | `scripts/bash/03_filter_cpg_sites.sh` | Coverage thresholding; splits sites by methylation state (unmethylated, intermediate, methylated) |
| 04 | `scripts/bash/04_run_all.sh` | Main orchestrator: runs 01→07 end-to-end, stops at first failure, records manifest |
| 05 | `scripts/bash/05_annotate_regions.sh` | Gene, promoter, and CpG-island region aggregation via `bedtools intersect` + POSIX awk |
| 06 | `scripts/bash/06_summary_stats.sh` | Global distribution, per-chromosome, and feature-class TSV summaries; cross-checks stage 03 |
| 07 | `scripts/R/07_generate_report.R` | Generates 6 publication figures using base R graphics only |
| 08 | `scripts/bash/08_archive_results.sh` | Sample-filtered archiving with SHA-256 checksums and provenance logging |
| 09 | `scripts/R/09_compare_samples.R` | Cross-sample differential methylation comparison at gene, promoter, and island levels |

---

## Requirements

- `bedtools` >= 2.31, `tabix`/`htslib`, `bgzip`
- R >= 4.3 (**base graphics and stats only; zero R packages required**)
- Run dependency check:
  ```bash
  bash scripts/check_dependencies.sh --pipeline methylation
  ```

---

## Quick Start

### Option A: Using the Interactive Web UI (Recommended)
Launch the configuration dashboard to set sample IDs, file paths, and reference resources visually:
```bash
python3 scripts/launch_ui.py --pipeline methylation
```

### Option B: Using the Command Line
```bash
cd pipelines/methylation

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit both config files with sample and reference paths, then:

bash scripts/bash/04_run_all.sh
```

---

## Output File Registry

- `results/qc_summary/<sample_id>.cpg_site_summary.tsv`
- `results/qc_summary/<sample_id>.gene_methylation_summary.tsv`
- `results/qc_summary/<sample_id>.promoter_methylation_summary.tsv`
- `results/qc_summary/<sample_id>.cpg_island_summary.tsv`
- `results/qc_summary/figures/<sample_id>.01_methylation_distribution.png`
- `results/qc_summary/figures/<sample_id>.02_chromosome_methylation.png`
- `results/qc_summary/figures/<sample_id>.03_genomic_feature_methylation.png`
- `results/qc_summary/figures/<sample_id>.04_coverage_distribution.png`
- `results/qc_summary/figures/<sample_id>.05_methylation_states.png`
- `results/qc_summary/figures/<sample_id>.06_gene_type_methylation.png`

---

## Automated Testing

Run the pipeline test suite against synthetic test fixtures:
```bash
bash pipelines/methylation/tests/run_tests.sh
```
