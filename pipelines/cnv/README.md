# pipelines/cnv — Copy Number Variation

**Status: scaffold only. No analysis logic exists yet.**

This directory holds the skeleton the CNV pipeline will be built into, following
the same shape as `pipelines/sv/` and `pipelines/methylation/`: numbered stage
scripts under `scripts/bash/` (and `scripts/R/` for anything that needs
plotting or statistics), a per-sample `config/pipeline_config.yaml`, and
`results/` / `logs/` / `tests/` / `docs/` alongside.

**Expected input:** `*.wf_cnv.vcf.gz`, as produced by Spectre via Epi2ME
`wf-human-variation`. `pipelines/sv/scripts/bash/01_validate_inputs.sh`
already checks for this file's presence as an adjacent-file sanity check, but
does not parse or analyze it — that logic belongs here, not there.

**What exists right now:**
- `scripts/bash/00_setup_env.sh` — sources `common/lib_common.sh` (the same
  shared config-reading and tool-checking helpers every pipeline in the suite
  uses) and resolves the generic config fields every pipeline needs
  (`sample_id`, `raw_sample_prefix`, `input_dir`, `output_dir`, `work_dir`,
  `log_dir`, `archive_root`). It does **not** yet check CNV-specific tools or
  read CNV-specific filtering config, because there is no stage 01+ yet to
  need them. Extend `common_guard_duplicate_keys`'s key list and
  `common_check_tools`'s tool list here once real stages exist.
- `config/pipeline_config.example.yaml`, `config/reference_paths.example.yaml`
  — templates covering only the fields every pipeline shares. Add
  CNV-specific filtering/annotation fields (analogous to SV's
  `filtering:`/`sv_categories:` blocks or methylation's `filtering:`/
  `annotation:` blocks) as real stages are written.

**What does not exist:** input validation beyond the generic checks, VCF
parsing, CNV-specific filtering, annotation, summary statistics, figures,
archiving, and tests. None of `pipelines/sv/`'s or `pipelines/methylation/`'s
hard-won lessons (sample-ID-scoped file resolution, exhaustive multi-block
annotation parsing, archiving kept separate from orchestration, real
tool-version checks) have been *applied* here yet — they just haven't been
needed yet, since there's nothing to apply them to. Whoever builds this out
should read `CONTRIBUTING.md`'s "Adding a pipeline" section first.
