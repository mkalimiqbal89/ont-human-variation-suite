# ONT Epi2ME Methylation Pipeline

A reproducible, config-driven pipeline for genome-wide CpG methylation analysis
from Oxford Nanopore data processed by the Epi2ME `wf-human-variation` workflow.

Companion to [`ont-sv-pipeline`](../ont-sv-pipeline) and built to the same
conventions, so that the two can eventually be published as one suite covering
SV, CNV, SNV and methylation.

---

## Input

The workflow's modified-base output, produced by `modkit pileup`:

| Run type | Files |
|---|---|
| Unphased | `<prefix>.wf_mods.bedmethyl.gz` |
| Phased (`--phased`) | `<prefix>.wf_mods.1.bedmethyl.gz`, `.2.`, `.ungrouped.` |

Set `modifications.phased` in `config/pipeline_config.yaml` to match; stage 01
checks both directions and warns if the config and the actual output disagree.

### bedMethyl column layout (18 columns, modkit `pileup`)

| Col | Meaning |
|----:|---|
| 1–3 | chrom, start, end (0-based, half-open) |
| 4 | modification code — `m` = 5mC, `h` = 5hmC, `a` = 6mA |
| 5 | score |
| 6 | strand |
| 7–9 | thick start, thick end, item RGB |
| **10** | **valid coverage** |
| **11** | **percent modified** |
| **12** | **modified read count** |
| **13** | **canonical (unmodified) read count** |
| 14 | other modification count |
| 15 | delete count |
| 16 | fail count |
| 17 | diff count |
| 18 | no-call count |

Columns 10–13 carry all the signal this pipeline uses.

---

## Stages

| Stage | Script | Purpose |
|---|---|---|
| 00 | `scripts/bash/00_setup_env.sh` | Tool check, version capture, config export, directory creation. **Must be sourced.** |
| 01 | `scripts/bash/01_validate_inputs.sh` | Input existence, gzip integrity, column layout, mod codes, contig-naming concordance, reference bundle |
| 02 | `scripts/bash/02_bedmethyl_to_tsv.sh` | Single-pass extraction of the primary mod code; full-file structural assertions; coverage histogram; bgzip + tabix |
| 03 | `scripts/bash/03_filter_cpg_sites.sh` | Coverage threshold; splits sites by methylation state; two reconciliation assertions |
| 04 | `scripts/bash/04_run_all.sh` | Orchestration: runs 01→07, stops at first failure, writes a provenance manifest |
| 05 | `scripts/bash/05_annotate_regions.sh` | Gene / promoter / CpG-island aggregation via `bedtools intersect` + awk |
| 06 | `scripts/bash/06_summary_stats.sh` | Global, distribution, per-chromosome and feature-class TSVs; cross-checks stage 03 |
| 07 | `scripts/R/07_generate_report.R` | Six publication figures, base R graphics only |
| 08 | `scripts/bash/08_archive_results.sh` | Sample-filtered archiving with verified checksums and provenance |
| 09 | `scripts/R/09_compare_samples.R` | Cross-sample comparison at gene, promoter and island level |
| 10 | `scripts/bash/10_run_comparison.sh` | Sample-sheet validation and comparison orchestration |

### Why stages 05 and 06 are bash, not R

Two reasons, and they point the same way.

The design principle from the demo walkthrough: every step produces a raw `.tsv`
and a human-readable summary, and R is brought in only to turn those TSVs into
figures. That separation makes each stage independently debuggable.

The practical constraint: the filtered set is ~27.9 M sites for a whole genome.
`awk` streams that in constant memory; reading it into R, even with
`data.table`, wants several GB. R is therefore used only where the input is
already small — figures (07) and cross-sample statistics on per-feature tables
(09).

---

## Quick start

> **Pasting these into zsh (the macOS default shell):** zsh does not treat `#`
> as a comment in interactive shells unless `interactive_comments` is set, so
> pasting a block that contains comment lines will produce
> `command not found: #`, and a trailing comment on a command line is passed to
> the script as an argument. Either run `setopt interactive_comments` once, or
> paste the commands without the comments. Stages 00–02 defensively ignore a
> `#`-leading argument, but the shell-level errors are outside their control.

```bash
# 1. Copy the templates and fill in real values
cp config/pipeline_config.example.yaml config/pipeline_config.yaml
cp config/reference_paths.example.yaml config/reference_paths.yaml

# 2. Check the environment (note: sourced, not executed)
source scripts/bash/00_setup_env.sh

# 3. Validate inputs before spending time downstream
bash scripts/bash/01_validate_inputs.sh
```

```bash
# 4. Extract the primary modification code (single full pass over the input)
bash scripts/bash/02_bedmethyl_to_tsv.sh

# 5. Apply the coverage threshold and split by methylation state
bash scripts/bash/03_filter_cpg_sites.sh

# 6. Aggregate to genes, promoters and CpG islands
bash scripts/bash/05_annotate_regions.sh

# 7. Genome-wide summary tables
bash scripts/bash/06_summary_stats.sh
```

Or run the whole per-sample pipeline in one command:

```bash
bash scripts/bash/04_run_all.sh
```

It stops at the first failing stage rather than running later ones on incomplete
output, and prints the resume command. Useful flags: `--list`, `--dry-run`,
`--from NN`, `--to NN`, `--only NN`, `--skip-validate`.

### Provenance manifest

Every `04_run_all.sh` run writes `results/<SAMPLE_ID>.run_manifest.tsv`
containing what a methods section needs: the config file and its SHA-256, the
reference config and each annotation BED's SHA-256, genome build, modification
code, coverage threshold, the active contig-selection mode and pattern, every
tool version, host and OS, the git commit and whether the tree was dirty, and
per-stage wall time and status.

Recording the git commit takes more care than it looks: `git rev-parse HEAD`
prints the literal string `HEAD` for a repository with no commits, and git
searches parent directories, so it will happily report an enclosing
repository's commit. The manifest distinguishes `not_a_git_repo`,
`repo_initialised_but_no_commits` and `enclosing_repo_only:<path>` from a real
SHA, and a test asserts the value is never the literal `HEAD`.

`SKIP_INTEGRITY=true` skips stage 01's full-file gzip test. Useful when
re-running validation repeatedly during development; leave it on for a real run.

Stage 00 is safe to `source` from bash or zsh. Stages 01 onward are bash
scripts — invoke them with `bash`, not `sh` or `zsh`.

---

## Tests

```bash
bash tests/run_tests.sh
```

194 content-level assertions across synthetic fixtures (more where a real R and
bedtools are present), each isolating one failure mode. Every test case gets its own `mktemp` directory outside the repo,
and the suite finishes by checking that nothing leaked into `results/` or
`logs/`. Per-case isolation is not cosmetic: while writing the stage 03 tests,
sharing one directory per fixture name let a later test see an earlier test's
output and pass for the wrong reason.

Assertions are deliberately content-level rather than exit-code-only. Exit codes
miss the bugs that matter — a stage can succeed while silently dropping records,
double-counting, or leaving a truncated output that looks plausible. The suite
therefore checks things like "the `h` records were counted but not retained",
"the coverage histogram sums to exactly the retained site count", and
"coverage-weighted methylation on a fixture with known values is exactly
50.0000%".

Regenerate fixtures alone with `bash tests/make_fixtures.sh`.

Test configs are derived from `config/pipeline_config.example.yaml` by
substitution, not hand-written. A hand-maintained copy drifts: a `sed` override
targeting a key the test config never had does nothing, so the test passes
against the default while appearing to exercise the setting. That happened twice
before the tests were restructured this way.

---

## Figures (stage 07)

Base R graphics, no packages. This is deliberate: the stage must run unattended
on an HPC where the R library path may differ from the login shell's, and a
pipeline that fails at the final stage because `ggplot2` is missing has wasted
the entire run. The config is read with a small regex reader mirroring
`yaml_get()` rather than the `yaml` package, for the same reason.

| Figure | Content |
|---|---|
| `fig1_methylation_distribution` | Per-site methylation in 10% bins |
| `fig2_chromosome_methylation` | Coverage-weighted methylation per chromosome, against the genome mean |
| `fig3_feature_class_methylation` | Boxplots: gene bodies, promoters, CpG islands, with protein-coding subsets |
| `fig4_coverage_distribution` | Coverage histogram with the chosen threshold marked |
| `fig5_methylation_states` | Site counts by unmethylated / intermediate / methylated |
| `fig6_promoter_by_gene_type` | Promoter methylation median and IQR by `gene_type` |

Each figure is skipped with a warning naming the stage to run if its input is
absent, so a partial pipeline still produces what it can. The palette is
Okabe-Ito, which is colour-blind safe and survives greyscale printing.

---

## Archiving (stage 08)

```
<archive_root>/
├── archive_index.tsv                one row per archived run
└── <SAMPLE_ID>/<RUN_STAMP>/
    ├── results/                     only files named <SAMPLE_ID>.*
    ├── logs/
    ├── config/                      pipeline_config.yaml + reference_paths.yaml
    ├── checksums.sha256             generated AND verified by read-back
    └── provenance.tsv
```

Runs are stamped, so re-analysis never overwrites history. Every guard in this
stage exists because the equivalent SV stage got it wrong first:

- **Files are selected by exact basename** (`-name "<SAMPLE_ID>.*"`), never by
  copying `results/` wholesale — that directory is shared between samples. A
  post-copy assertion re-checks every archived basename and deletes the archive
  if any file does not belong to this sample.
- **`mkdir` and `cp` failures abort.** Silent continuation produces an archive
  that looks complete and is not.
- **Preflight checks run before any copying**: archive root creatable, writable
  (tested by writing, since `-w` lies on some network mounts), and enough free
  space with a 20% margin.
- **Checksums are verified by reading back what was written.** A generated
  checksum file proves nothing on its own; a test confirms that tampering with
  an archived file makes verification fail.
- **The index is appended only after everything else succeeded**, so it never
  points at an incomplete archive.

## Comparing samples (stages 09/10)

Per-sample analysis first, then:

```bash
bash scripts/bash/10_run_comparison.sh --template config/sample_sheet.tsv
# edit the sheet, then
bash scripts/bash/10_run_comparison.sh config/sample_sheet.tsv --dry-run
bash scripts/bash/10_run_comparison.sh config/sample_sheet.tsv
```

The sheet is tab-separated with `sample_id` and `results_dir` required, `label`
and `group` optional. Samples may share a `results_dir` — outputs are
sample-prefixed, so one directory holds several samples. Two or more samples are
supported: with two you get deltas and a scatter plot, with more you get range,
standard deviation and a correlation matrix.

Stage 10 validates the sheet before R starts — missing files, duplicate
`sample_id`, fewer than two samples, and the very common case of a
space-separated sheet all fail with a message naming the actual problem.

Comparison works from stage 05's per-feature tables, which already carry
`Total_coverage` and `Modified_read_calls` per feature. So this stage never reads
a site-level file: inputs are tens of thousands of rows, and R is the right tool.

### The statistics need reading before use

**[`docs/COMPARISON_CAVEATS.md`](docs/COMPARISON_CAVEATS.md) — read it before
quoting a q-value.** In short: the pooled-read test treats every read call as an
independent observation, but nanopore reads span many CpGs so calls within a
feature are correlated; with one sample per group there is no biological
replication; and at this depth a one-point difference clears any significance
threshold. The p-values are anti-conservative.

The output is built to make this hard to ignore. Features are classified as
`differential` only when they clear **both** an effect-size threshold and the FDR
cutoff, and `significant_only` is reported as its own count — those are the
features where statistics say yes and biology very likely says no. Stage 09 also
prints the caveat to its own run log. `stat_test: none` reports effect sizes
without p-values at all, which is a defensible choice for a single pair.

---

`*.all_cov.*` is excluded by default — it is a pure intermediate, ~414 MB per
sample, and regenerable from the raw bedMethyl in about three minutes. Set
`archive.include_intermediates: true` to keep it.

---

## Results layout

```
results/
├── 01_filtered/              filtered primary-mod-code sites
├── 02_qc/                    QC summaries
├── 03_global_methylation/    genome-wide methylation + distribution
├── 04_chromosome_summary/    per-chromosome methylation
├── 05_annotated_cpgs/        CpGs intersected with genes / promoters / CGIs
├── 06_gene_summary/          per-gene aggregation
├── 07_promoter_summary/      per-promoter aggregation
├── 08_cpg_island_summary/    per-island aggregation
├── 09_figures/               publication figures
└── 10_igv_tracks/            IGV-loadable tracks
```

This is a renumbered, contiguous version of the layout used in
`Projects/Methylation_Demo`, where the annotation and gene-summary directories
were `11_Annotated_CpGs` and `12_Gene_Summary`.

---

## Design notes

Choices here are deliberate, and several were paid for during development of
`ont-sv-pipeline`:

- **Tool checks validate output, not just presence.** A binary can sit on
  `PATH` and still be unusable through a broken shared library. Stage 00 checks
  the exit status *and* inspects the output string.
- **Duplicate config keys are a hard error.** The config reader is a
  grep-based `yaml_get()`, which returns the first match anywhere in the file.
  A duplicated key would therefore be silently ignored, so stage 00 refuses to
  run rather than trusting the file.
- **Contig naming is verified across every input.** If the bedMethyl is
  UCSC-style (`chr1`) and an annotation BED is Ensembl-style (`1`), every
  `bedtools intersect` returns zero rows. That presents as "this sample has no
  methylation in genes", not as an error, so it is checked explicitly.
- **Expensive checks are paid for once.** Stage 01 samples the head of the file
  for its column check; the full-file column assertion lives in stage 02, which
  already reads every record. Validation should not re-read 660 MB to learn what
  the next stage learns for free.
- **Genes aggregate on `gene_id`, never `gene_name`.** GENCODE symbols are not
  unique.
- **Contig selection is an allowlist, not a blocklist.** A blocklist can only
  reject the junk naming patterns you thought to enumerate. GRCh38 analysis sets
  ship scaffolds — `GL000191.1`, `HLA-A*01:01:01:01`, `chr1_KI270762v1_fix` —
  that match none of `_random$`, `^chrUn` or `_alt$`, so a blocklist lets them
  through silently. `filtering.primary_contigs_only: true` keeps only what
  `include_contigs_regex` names. Both modes are implemented and exactly one is
  active per run; the stage log states which, and `tests/run_tests.sh` asserts
  the difference in both directions.
- **`chrM` is excluded by default.** Mitochondrial DNA is essentially
  unmethylated in mammals and carries very high coverage, so leaving it in drags
  coverage-weighted genome-wide averages down. On SAMPLE_02 chrM contributed 433
  sites; in the demo dataset, 435 sites at 71,466 total coverage and 0.17%
  methylation.
- **Retained contig count is checked but never fatal.** 23 instead of 24 is the
  expected result for a female sample with no chrY calls, so this warns and says
  what to verify rather than aborting a three-minute run.

---

## Reference resources

Shared across all projects, built once. See `docs/REFERENCE_PREPARATION.md`.

| File | Source |
|---|---|
| `genes.bed` (BED7) | GENCODE v50 GTF, `gene` features, GTF→BED coordinate conversion |
| `promoters_2kb.bed` (BED7) | strand-aware TSS ±2 kb via `bedtools slop` |
| `cpg_islands_hg38.bed` (BED10) | UCSC `cpgIslandExt`, index column dropped, whitespace in names replaced |
| `hg38.chrom.sizes` | UCSC goldenPath |

The genome FASTA is **not** required: modkit has already produced the pileup.

---

## Sample ID mapping

Kept consistent with `ont-sv-pipeline`:

| Pipeline ID | Epi2ME prefix |
|---|---|
| `SAMPLE_01` | `HLH_S0002_BL_EPI2ME_2.7.2` |
| `SAMPLE_02` | `HLH_S0001_BL_EPI2ME_2.7.2` |
