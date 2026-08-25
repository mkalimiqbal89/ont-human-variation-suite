# ONT Epi2ME Structural Variant Pipeline

**Status: Released `v1.0.0`**

A reproducible, config-driven downstream pipeline for categorizing structural variants (SVs) from Oxford Nanopore long-read whole-genome sequencing processed by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation). Filters and splits Sniffles2-called SVs into deletions, insertions, duplications, inversions, translocations, and SnpEff-confirmed gene fusions/disruptions.

---

## Input Data Specifications

| Property | Specification |
|---|---|
| Variant Class | Structural Variants (SVs, >30 bp & BND translocations) |
| Upstream Workflow | Epi2ME `wf-human-variation` |
| Primary Input File | `<prefix>.wf_sv.vcf.gz` (Sniffles2 VCF) |
| Key Information | SVTYPE, SVLEN, END, QUAL, VAF, DV/DR, and SnpEff `ANN` tags |

---

## Pipeline Stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `scripts/bash/00_setup_env.sh` | Tool check, version capture, config export, directory setup. **Must be sourced.** |
| 01 | `scripts/bash/01_validate_inputs.sh` | Input VCF file existence, tabix index checks, reference bundle verification |
| 02 | `scripts/bash/02_vcf_to_tsv.sh` | Flattens Sniffles2 VCF into an analysis-ready TSV, extracting scalar fields and first ANN gene symbol |
| 03 | `scripts/bash/03_filter_sv_categories.sh` | Applies QC thresholds (length, read support, VAF, QUAL) and splits calls into per-category TSVs |
| 04 | `scripts/bash/04_run_all.sh` | Main orchestrator: runs 01→07 end-to-end, stops at first failure, records manifest |
| 05 | `scripts/R/05_annotate_variants.R` | Extracts SnpEff breakpoint ANN classifications to construct gene fusion and disruption tables |
| 06 | `scripts/R/06_summary_stats.R` | Computes per-category variant counts, mean/median length, VAF distribution, and fusion breakdown |
| 07 | `scripts/R/07_generate_report.R` | Generates 5 overview figures (PNGs) across all SV classes |
| 08 | `scripts/bash/08_archive_results.sh` | Longitudinal archiving with SHA-256 integrity verification and index logging |

---

## Requirements

- `bcftools` >= 1.20, `bedtools` >= 2.31, `tabix`/`htslib`
- R >= 4.3 (requires `yaml`, `ggplot2`, `scales` packages for SV plots)
- Run dependency check:
  ```bash
  bash scripts/check_dependencies.sh --pipeline sv
  ```

---

## Quick Start

### Option A: Using the Interactive Web UI (Recommended)
Launch the configuration dashboard to set sample IDs, file paths, and reference resources visually:
```bash
python3 scripts/launch_ui.py --pipeline sv
```

### Option B: Using the Command Line
```bash
cd pipelines/sv

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit both config files with sample and reference paths, then:

bash scripts/bash/04_run_all.sh
```

---

## Output File Registry

- `results/deletions/<sample_id>.deletions.tsv`
- `results/insertions/<sample_id>.insertions.tsv`
- `results/duplications/<sample_id>.duplications.tsv`
- `results/inversions/<sample_id>.inversions.tsv`
- `results/translocations/<sample_id>.translocations.tsv`
- `results/complex_rearrangements/<sample_id>.complex_rearrangements.tsv`
- `results/gene_fusions/<sample_id>.gene_fusions.tsv`
- `results/gene_fusions/<sample_id>.gene_disruptions.tsv`
- `results/qc_summary/<sample_id>.summary_statistics.tsv`
- `results/qc_summary/figures/<sample_id>.01_sv_counts_by_category.png`
- `results/qc_summary/figures/<sample_id>.02_sv_length_distribution.png`
- `results/qc_summary/figures/<sample_id>.03_vaf_distribution.png`
- `results/qc_summary/figures/<sample_id>.04_chromosome_distribution.png`
- `results/qc_summary/figures/<sample_id>.05_translocation_breakdown.png`

---

## Automated Testing

Run the pipeline test suite against synthetic test fixtures:
```bash
bash pipelines/sv/tests/run_tests.sh
```
