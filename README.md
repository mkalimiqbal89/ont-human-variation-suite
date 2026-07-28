# ONT Human Variation Suite

Reproducible, config-driven pipelines for **downstream** analysis of Oxford
Nanopore sequencing processed by the Epi2ME
[`wf-human-variation`](https://github.com/epi2me-labs/wf-human-variation)
workflow.

`wf-human-variation` produces the calls. This suite turns them into per-sample
tables, figures and archived records, and compares samples to one another — with
the assertions and provenance needed for that work to be defensible.

---

## Pipelines

| Pipeline | Variant class | Status | Input from `wf-human-variation` |
|---|---|---|---|
| [`pipelines/sv`](pipelines/sv) | Structural variants | Released `v1.0.0` | `*.wf_sv.vcf.gz` (Sniffles2) |
| [`pipelines/methylation`](pipelines/methylation) | CpG methylation | Feature-complete, in validation | `*.wf_mods.bedmethyl.gz` (modkit) |
| `pipelines/cnv` | Copy number | Planned | `*.wf_cnv.vcf.gz` (Spectre) |
| `pipelines/snv` | Small variants | Planned | `*.wf_snp.vcf.gz` (Clair3) |

Each pipeline is self-contained: its own `config/`, `scripts/`, `tests/`,
`.gitignore` and `README`. **Start with the pipeline README, not this one.**

---

## Shared conventions

The pipelines deliberately look alike, so learning one teaches the others.

- **Numbered stages.** `00_setup_env.sh` checks tools and exports config;
  `01_validate_inputs.sh` refuses to proceed on bad input; later stages do the
  work; `04_run_all.sh` orchestrates and stops at the first failure.
- **Nothing hardcoded outside config.** One config file per sample, with a
  tracked `.example.yaml` template. Real configs are gitignored — they carry
  sample identifiers and absolute paths.
- **Shared reference resources**, built once and reused across projects:
  GENCODE-derived gene and promoter BEDs, UCSC CpG islands, chromosome sizes.
  Paths live in `reference_paths.yaml`, never in code.
- **Content-level tests.** Assertions check values and reconciliations, not just
  exit codes. Exit codes miss the failures that matter: a stage succeeding while
  silently dropping records, double-counting, or writing a truncated output that
  looks plausible to the next stage.
- **Fail loudly, leave nothing partial.** A failing stage removes its own
  incomplete output rather than leaving a file later stages will consume.
- **Provenance by default.** Each full run records config and reference
  checksums, tool versions, host, OS, git commit and per-stage timings.

---

## Requirements

- bash 4+ (stage 00 is also safe to `source` from zsh), GNU or BSD awk, coreutils
- [bedtools](https://bedtools.readthedocs.io) 2.31+
- [htslib](https://www.htslib.org) — `bgzip`, `tabix`
- [bcftools](https://samtools.github.io/bcftools/) 1.20+ *(SV pipeline only)*
- R 4.3+ — **base graphics and `stats` only; no R packages required**
- optional: `pigz`, for faster decompression of 600 MB+ inputs

Per-pipeline conda environments are provided as `environment.yml`.

Developed and tested on macOS (Apple silicon, BSD userland) and Linux HPC. Both
are supported deliberately, and the portability constraints that follow from it
are documented in `CONTRIBUTING.md`.

---

## Quick start

```bash
git clone https://github.com/mkalimiqbal89/ont-human-variation-suite.git
cd ont-human-variation-suite/pipelines/methylation

cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml
# edit both to point at your data and references, then:
source scripts/bash/00_setup_env.sh
bash scripts/bash/04_run_all.sh
```

Run the tests without any real data — they use synthetic fixtures only:

```bash
bash pipelines/methylation/tests/run_tests.sh
bash pipelines/sv/tests/run_tests.sh
```

---

## Documentation

- [`pipelines/methylation/docs/COMPARISON_CAVEATS.md`](pipelines/methylation/docs/COMPARISON_CAVEATS.md) — **read this before
  quoting a p-value from the cross-sample comparison.** Pooled-read tests on
  long-read data without biological replication are anti-conservative. The
  document explains why, and what the output is legitimately good for.
- [`docs/AI_USAGE.md`](docs/AI_USAGE.md) — disclosure of AI assistance used in
  development, including the errors it introduced and how each was caught.
- [`CHANGELOG.md`](CHANGELOG.md) — what changed, and why.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — support expectations, testing
  philosophy, portability rules.

Per-pipeline methods rationale lives in each pipeline's `docs/METHODS.md`.

---

## Status and scope

Research software developed in a cancer genomics research setting, currently
applied to haemophagocytic lymphohistiocytosis (HLH) samples. Offered in the hope
it is useful to others working downstream of `wf-human-variation`.

**This is not a validated diagnostic tool and must not be used for clinical
decision-making.**

Issues and pull requests are welcome — see `CONTRIBUTING.md` for what makes a
useful bug report, and please do not attach real patient data.

---

## Licence

MIT — see [`LICENSE`](LICENSE).

## Citation

See [`CITATION.cff`](CITATION.cff); GitHub renders a "Cite this repository"
button from it.
