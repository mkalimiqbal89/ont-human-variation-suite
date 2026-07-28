# AI usage disclosure

---

## Tools used

| Tool | Version / model | Where used |
|---|---|---|
| Claude (Anthropic), via Cowork | Opus 5 | Shell and R scripting|

---

## Nature and scope of assistance

**Moderately AI-assisted:**

- Implementation of the methylation pipeline stages (`00`–`10`), including the
  awk single-pass extraction, filtering, annotation aggregation and summary
  statistics, and the R figure and comparison stages.
- The regression test suite (`tests/run_tests.sh`, `tests/make_fixtures.sh`) and
  its synthetic fixtures.
- Refactoring for portability (BSD versus GNU `awk`, macOS versus Linux
  `sha256sum`/`shasum`, zsh versus bash).

**Author-led, not AI-generated:**

- The scientific question, sample selection, and the decision to build a suite
  covering SV, CNV, SNV and methylation.
- All real-data execution. Every run against real patient data
  was performed by the author on the author's own
  hardware; the AI assistant never had access to the raw sequencing data.
- Reference resource preparation (GENCODE v50 → `genes.bed`,
  `promoters_2kb.bed`, UCSC `cpgIslandExt` → `cpg_islands_hg38.bed`), carried
  out and validated step by step by the author before any pipeline code existed.
- Verification of every stage against real output, including the decision to
  stop and investigate when results looked wrong.

**Jointly determined:**

- Architectural choices were discussed and decided interactively rather than
  accepted as generated. Examples where the author's decision changed the design:
  keeping bash/awk for stages 05–06 instead of R (motivated by the author's
  stated principle that each stage should emit an inspectable TSV, and by the
  ~28 M-site scale); one config file per pipeline rather than a shared global
  one; and the monorepo structure of this suite.

---

## Verification performed

The following were checked against real data by the author, not asserted by the
assistant:

- Stage 06 reproduces, to four decimal places, summary values the author had
  previously derived by hand in a separate step-by-step walkthrough
  (1331 retained sites, 87028 total coverage, 9362 modified calls,
  40.3699% unweighted, 10.7575% coverage-weighted, and per-chromosome values for
  chr6, chr20 and chrM).
- Stage 05's `bedtools intersect` row counts (1308 gene, 4264 promoter,
  125 island overlaps) match output the author generated independently with
  bedtools 2.31.1.
- Biological plausibility checks that would fail under a coordinate or column
  error: all mitochondrial genes report < 0.13% methylation despite very high
  coverage, and a CpG island reports 1.13% while the gene body containing it
  reports 37.25%.
- The regression suite (200+ content-level assertions) was executed by the
  author on the author's machine, where the real `bedtools` and `R` are
  installed. Assertions were written to check values and reconciliations, not
  merely exit codes.

---

## Errors introduced by AI assistance and subsequently caught

Recorded because it is the honest picture of how this code was produced, and
because it evidences that review was real rather than nominal.

| Error | How it surfaced |
|---|---|
| `PIPESTATUS` read element-by-element, so all but the first were lost | Failed loudly under `set -u`; would have hidden `bgzip` failures |
| `tee >(gzip …)` for an optional output — bash does not reliably wait for process substitution to flush | Test caught a missing file; replaced with an awk-managed pipe closed in `END` |
| `boxplot()` on an empty list when no feature met the CpG minimum | Full-chain test failed on the author's machine, where R is real |
| Test harness used a shell variable for unique directories inside a command substitution, so every test shared one directory and one assertion passed for the wrong reason | Investigated an anomalous pass |
| Test config maintained as a heredoc that drifted from the real config; `sed` overrides silently matched nothing | Two tests passed against defaults while appearing to test a setting |
| `sed -i` run against `.git/config` on a FUSE-mounted volume, deleting the file | Immediately detected and restored; no data lost, nothing pushed |
| Root `.gitignore` patterns stopped matching after the SV pipeline moved into a subdirectory, unignoring real sample configs | Caught during the monorepo migration, before any commit |
| Coverage-imbalance blind spot: `min_shared_cpgs` required enough CpGs per sample but not comparable numbers, so the top "differential" hits were mappability artefacts in pseudogene and subtelomeric loci | Found by the author and assistant reading the first real two-sample output together; filter added |


---

## Author's confirmation

> I, Kaleem Iqbal, confirm that I reviewed, edited and validated all AI-assisted
> output in this repository; that I made the core design decisions; that I
> performed all analysis of real data myself; and that I am accountable for the
> correctness of the software and the accuracy of this disclosure.
