# ONT Human Variation Suite

Reproducible, config-driven downstream pipelines for Oxford Nanopore whole-genome sequencing processed by Epi2ME [`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation).

`wf-human-variation` produces the primary variant calls. This suite turns them into per-sample tables, figures, gene annotations, and archived records with the content assertions and provenance needed for defensible research.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21600726.svg)](https://doi.org/10.5281/zenodo.21600726)

---

## Pipelines Overview

| Pipeline | Variant Class | Status | Input File from `wf-human-variation` | Documentation |
|---|---|---|---|---|
| [`pipelines/sv`](pipelines/sv) | Structural Variants | Released `v1.0.0` | `<prefix>.wf_sv.vcf.gz` (Sniffles2) | [SV README](pipelines/sv/README.md) |
| [`pipelines/methylation`](pipelines/methylation) | CpG Methylation | Feature-complete | `<prefix>.wf_mods.bedmethyl.gz` (modkit) | [Methylation README](pipelines/methylation/README.md) |
| [`pipelines/cnv`](pipelines/cnv) | Copy Number | Released `v1.0.0` | `<prefix>.wf_cnv.vcf.gz` (Spectre) | [CNV README](pipelines/cnv/README.md) |
| [`pipelines/snv`](pipelines/snv) | Small Variants (SNV/Indel) | Released `v1.0.0` | `<prefix>.wf_snp.vcf.gz` (Clair3) | [SNV README](pipelines/snv/README.md) |

---

## Interactive Configuration Web UI

Instead of editing raw YAML files by hand, you can launch the **Interactive Web UI**:

```bash
python3 scripts/launch_ui.py
```

Features:
- **Visual File & Path Selector**: Set Sample ID, raw prefix, input directory, output directory, and genome references without touching YAML files.
- **Live Path Verification**: Includes a "Validate Paths" check to verify file existence on disk before saving.
- **Pipeline Runner & Real-Time Log Viewer**: Execute dependency checks and run pipelines directly from your browser with live streaming logs.

---

## Shared Architecture & Design Principles

All pipelines follow identical conventions:

- **Numbered Execution Stages**: `00_setup_env.sh` (exports env & tool checks), `01_validate_inputs.sh` (input verification), work stages (`02` to `07`), and `04_run_all.sh` (orchestration).
- **Configuration Scoping**: Parameters and reference resources live strictly in `config/pipeline_config.yaml` and `config/reference_paths.yaml`. Real config files are gitignored.
- **Content-Level Assertions**: Stage assertions verify record count reconciliations and value integrity, not just shell exit codes.
- **Fail-Loud Execution**: Failing stages remove partial outputs to prevent downstream corruption.
- **Provenance & Integrity**: Automated logging of host OS, tool versions, git commit, SHA-256 checksums, and stage runtimes.

---

## Requirements

- **OS**: macOS (Apple Silicon / Intel) or Linux HPC
- **Shell & Core**: bash 3.2+, GNU/BSD awk, POSIX coreutils
- **Bioinformatics Tools**: `bcftools` >= 1.20, `bedtools` >= 2.31, `htslib` (`bgzip`, `tabix`)
- **R Environment**: R >= 4.3 (**base graphics and `stats` only; no external R packages required for methylation, CNV, or SNV pipelines**)

Check dependencies across pipelines in one command:

```bash
bash scripts/check_dependencies.sh              # Report status
bash scripts/check_dependencies.sh --install    # Install missing packages via brew/apt/conda
bash scripts/check_dependencies.sh --pipeline sv
```

---

## Suite Quick Start

### 1. Clone the repository
```bash
git clone https://github.com/mkalimiqbal89/ont-human-variation-suite.git
cd ont-human-variation-suite
```

### 2. Configure via Interactive Web UI (Recommended)
```bash
python3 scripts/launch_ui.py
```

### 3. Or Configure via Command Line
```bash
cd pipelines/sv # or methylation, cnv, snv

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# Edit config files to specify your sample and reference paths, then:

source scripts/bash/00_setup_env.sh
bash scripts/bash/04_run_all.sh

# Archive results to institutional storage:
bash scripts/bash/08_archive_results.sh
```

---

## Automated Test Suites

Run the regression test suites across pipelines using synthetic test fixtures:

```bash
bash pipelines/sv/tests/run_tests.sh
bash pipelines/methylation/tests/run_tests.sh
bash pipelines/cnv/tests/run_tests.sh
bash pipelines/snv/tests/run_tests.sh
```

---

## Documentation Index

- [`pipelines/methylation/docs/COMPARISON_CAVEATS.md`](pipelines/methylation/docs/COMPARISON_CAVEATS.md) — Statistical considerations for cross-sample methylation comparisons.
- [`docs/AI_USAGE.md`](docs/AI_USAGE.md) — Disclosure of AI assistance in suite development.
- [`CHANGELOG.md`](CHANGELOG.md) — Release notes and version history.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — Testing philosophy, portability rules, and pull request guidelines.

---

## License & Citation

- **License**: MIT — see [`LICENSE`](LICENSE).
- **Citation**: See [`CITATION.cff`](CITATION.cff).
