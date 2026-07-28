# Changelog

Notable changes to the ONT Human Variation Suite. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

---

## [Unreleased]

### Added
- **`common/lib_common.sh`** — shared config reader, duplicate-key guard,
  argument guard, path resolution and tool checker. Both pipelines now source it
  instead of carrying their own copies.
- **`scripts/check_dependencies.sh`** — reports required tools and R packages
  per pipeline, detects the platform and package manager, and can install what
  is missing. Runs each tool rather than only testing `PATH`.

### Fixed
- **`declare -A` removed from both pipelines.** Associative arrays are bash 4+
  and macOS ships bash 3.2, where the subscripts are evaluated *arithmetically*
  instead of rejected. In the SV pipeline `[deletions]` made bash resolve
  `deletions` as a variable and abort under `set -u`, so stage 03 exited 0 having
  written no category files and stage 05 failed later on a missing input. The
  methylation orchestrator used numeric keys, so `[01]` silently resolved to `1`
  and it worked only by accident — adding a stage `08` would have broken it,
  since `08` is not valid octal. The SV pipeline had therefore never worked on
  macOS despite `v1.0.0` being verified on Linux.
- Stage 03 (SV) now asserts it actually produced its category files rather than
  exiting 0 on an empty loop.
- The SV `00_setup_env.sh` gained the zsh-safe sourcing fix, the `#`-argument
  guard and the duplicate-key guard, which previously existed only in the
  methylation copy.

---

## [0.3.0] — 2026-07-28

### Added
- **Continuous integration** (`.github/workflows/tests.yml`) running both test
  suites on `ubuntu-latest` and `macos-latest`, plus a shell syntax pass, an
  assertion that the fixtures are byte-reproducible, and a check that a full test
  run leaves the working tree clean. Testing both platforms is the point: this is
  developed on macOS and run on a Linux HPC, and nearly every environment bug has
  come from that gap.
- A CI lint rejecting bash-4-only syntax, since a runner cannot easily provide
  bash 3.2.
- `pipelines/methylation/docs/METHODS.md` — rationale for every threshold,
  grounded in the validation runs.

### Fixed
- Test fixtures are gzipped with `-n` so they carry no timestamp and regenerate
  byte-identically; previously every test run dirtied the working tree.
- Three sheet-validation tests depended on the gitignored real config, so they
  passed on a developer machine and failed on a clean checkout — and two of them
  had been passing for the wrong reason.
- The SV test suite must run from its own directory: its reference config uses
  relative paths checked with a bare `[[ -f ]]`, which resolves against the cwd.
- `sed -i` without an argument in the SV test suite, which BSD sed rejects.

### Added
- **Methylation pipeline** (`pipelines/methylation`), stages 00–10, for CpG
  methylation from `modkit` bedMethyl output.
  - Single-pass extraction of the primary modification code with full-file
    structural assertions: column count, numeric and range checks, modified plus
    canonical counts never exceeding coverage, coordinate sortedness, and
    contiguous contig blocks.
  - Contig selection as an allowlist (`^chr([1-9]|1[0-9]|2[0-2]|X|Y)$`) rather
    than a blocklist. A blocklist only rejects the junk naming patterns you
    thought to enumerate; GRCh38 analysis sets ship scaffolds such as
    `GL000191.1`, `HLA-A*01:01:01:01` and `*_fix` that match none of the usual
    ones.
  - Coverage threshold plus per-site methylation-state splits
    (unmethylated / intermediate / methylated) with two reconciliation
    assertions that cannot hold by accident.
  - Gene, promoter and CpG-island aggregation. Genes keyed on `gene_id`, islands
    keyed on coordinates — both load-bearing, see below.
  - Genome-wide, distribution, per-chromosome, feature-class and per-`gene_type`
    summary tables, with a cross-stage consistency check against the filtering
    stage.
  - Six publication figures in base R graphics, no package dependencies.
  - Sample-filtered archiving with checksums verified by read-back, provenance
    metadata and a running archive index.
  - Cross-sample comparison for two or more samples, with effect sizes, optional
    Fisher/chi-squared tests, BH correction, and figures.
- **Provenance manifest** written by every full run: config and reference
  checksums, tool versions, host, OS, git commit and dirty state, per-stage
  timings.
- **`docs/COMPARISON_CAVEATS.md`** documenting why pooled-read p-values on
  long-read data without biological replication are anti-conservative, and what
  the comparison output is legitimately good for.
- **`docs/AI_USAGE.md`** disclosing AI assistance used in development, including
  errors introduced and how each was caught.
- Regression suite for the methylation pipeline: 200+ content-level assertions
  across synthetic fixtures, each isolating one failure mode.

### Changed
- **Repository restructured as a suite.** The SV pipeline moved from the
  repository root to `pipelines/sv/`; the repository was renamed from
  `ont_structural_variations_pipeline` to `ont-human-variation-suite`. GitHub
  redirects the old URL, so existing clones keep working.
- Each pipeline now owns its own `.gitignore`. Patterns containing a slash are
  anchored to the directory holding the file, so moving the SV pipeline into a
  subdirectory silently unignored its real sample configs. Keeping the rules
  beside the pipeline makes that class of mistake impossible; the suite-level
  file holds only slash-free defensive patterns.

### Fixed
- **Coverage imbalance in cross-sample comparison.** `min_shared_cpgs` required
  enough covered CpGs in each sample but not *comparable* numbers. On the first
  real two-sample run the top promoter hits had 10 covered CpGs in one sample
  against 173 in the other — mappability differences in pseudogene and
  subtelomeric loci, not methylation differences. Added `max_cpg_ratio` and
  `max_coverage_ratio`; such features are reported with their ratios but labelled
  `imbalanced` rather than `differential`, and sort below balanced features.
- `PIPESTATUS` was read element by element, so every element after the first was
  lost — the array is rebuilt by each command, including a plain assignment. This
  would have hidden `bgzip` failures.
- Optional per-CpG intermediates were written with `tee >(gzip …)`; bash does not
  reliably wait for a process substitution to flush, so the file could be
  truncated or absent. Now written from awk and closed deterministically.
- Figure stage crashed on `boxplot()` of an empty list when no feature met the
  minimum CpG count — the normal case for a small targeted panel. Now skipped,
  with an error handler that names the figure in progress and deletes partial
  PNGs.
- Boxplot colours were taken from the front of the palette rather than subset by
  the same mask as the data, so an empty middle group would mis-colour the
  survivors.
- Stage 00 could not be sourced from zsh, the macOS default shell: `BASH_SOURCE`
  does not exist there, so the repository root resolved to the caller's working
  directory. Falls back to `$0`.
- zsh does not strip inline `#` comments in interactive shells, so a command
  pasted from documentation delivered `#` as the config path argument. Stages now
  ignore a `#`-leading argument and explain why.
- The run manifest is written to the `results/` root, which a `results/*/*`
  pattern never matched — it would have been committed with real sample IDs,
  absolute paths and hostname.

---

## [1.0.0] — 2026-07-26

### Added
- Initial release of the **structural variant pipeline**: validation, Sniffles2
  VCF flattening, QC filtering and per-category splitting, gene fusion and
  disruption classification from SnpEff `ANN` breakpoint effects, summary
  statistics, figures, and archiving with checksums and provenance.
- Regression suite with content-level assertions, including the real `ANN`
  multi-block parsing edge case and a cross-sample archive contamination check.
- Validated against two clinical research samples.

*Released while the repository was named `ont_structural_variations_pipeline`.*
