# Contributing

Thanks for looking. This is research software maintained by one person in a
cancer genomics research setting, so expectations are set accordingly below
rather than left implicit.

---

## Support expectations

- **Issues are welcome** and will normally get a first response within two weeks.
- **This is not a funded project** with dedicated maintenance time. Response may
  be slower during research project commitments.
- **Bug reports about correctness take priority** over feature requests. If a
  pipeline produces a wrong number, that is the most important kind of issue you
  can file.
- **Not a diagnostic tool.** Nothing here is validated for clinical
  decision-making, and no support can be offered for clinical use.

---

## Reporting a bug

The most useful report includes:

1. Which pipeline and which stage (`pipelines/methylation`, stage `05`).
2. The stage's log from `logs/` — every stage writes a timestamped log.
3. The relevant QC TSV from `results/02_qc/`, which carries the assertion
   counters.
4. Tool versions from `logs/tool_versions.txt`.
5. Whether `tests/run_tests.sh` passes on your machine.

**Do not attach real patient data.** If a bug only reproduces on real data,
describe the shape of the input — record count, contig naming, column count,
coverage range — or construct a synthetic fixture that reproduces it. The test
fixtures in `tests/fixtures/` show the pattern; `tests/make_fixtures.sh`
generates them all and is a good starting point.

---

## Proposing a change

1. Open an issue first for anything beyond a typo, so effort is not wasted.
2. Fork, branch from `main`.
3. **Add or extend a test.** See the testing philosophy below — this is the part
   that matters most.
4. Run the affected pipeline's test suite and include the pass/fail line.
5. Keep the diff focused. Restructuring and behaviour changes in one commit are
   hard to review and hard to revert.

---

## Testing philosophy

Assertions must be **content-level**, not exit-code-level. An exit code tells you
a stage ran; it does not tell you the stage was right. Every failure that has
actually mattered in this project's history was invisible to an exit code:

- a filter silently dropping records
- a category boundary with an overlap, so sites counted twice
- aggregation keyed on a non-unique field, silently merging 32,038 CpG islands
  into 485 rows
- a truncated output that looked plausible to the next stage

So tests assert reconciliations and known values:

```
input == dropped_low + dropped_high + retained
retained == unmethylated + intermediate + methylated
sum(per-feature CpG counts) == rows in the intersect output
coverage-weighted methylation on a known fixture == exactly 50.0000%
```

A test that cannot fail is not a test. When adding a guard, add a fixture that
trips it, and check the guard actually fires — several tests in this repo passed
for the wrong reason until that was verified.

---

## Code conventions

**Portability is not optional.** This runs on macOS (BSD userland) and Linux HPC:

- **No bash 4+ syntax.** macOS ships bash **3.2** as `/bin/bash`, so no
  associative arrays (`declare -A`), `mapfile`/`readarray`, or `${var,,}` /
  `${var^^}`. Use a space-separated list plus a `case` function instead. This
  one bites quietly: on 3.2 `declare -A` is not rejected — the subscripts are
  evaluated *arithmetically*, so `[deletions]` made bash resolve `deletions` as
  a variable and abort under `set -u`, while numeric keys like `[01]` silently
  resolved to `1` and appeared to work. CI greps for these constructs.
- POSIX awk only. No `gensub`, `asort`, or `length(array)` — BSD awk lacks them.
- `sha256sum` on Linux, `shasum -a 256` on macOS. Check for both.
- `wc -l` pads output with spaces on BSD. Pipe through `tr -d ' '`.
- `source`d scripts must work in zsh as well as bash. `BASH_SOURCE` does not
  exist in zsh; fall back to `$0`.
- Never use `sed -i` on a path that might be on a network or FUSE-mounted volume.
  It writes a temp file and renames, which can fail destructively.

**Other conventions:**

- Nothing hardcoded outside `config/`. If a script needs a path or threshold, it
  comes from the config.
- Every stage validates its inputs and refuses to run on bad ones.
- A stage that fails removes its own partial output. A half-written file that the
  next stage happily consumes is worse than no file.
- Comment the *why*, especially where a choice looks odd. Much of this codebase's
  commentary explains a real bug that motivated the current shape; that context
  is the most valuable thing in the file.

---

## Adding a pipeline

New variant classes go in `pipelines/<name>/` following the existing layout:
numbered stages, `config/` with tracked `.example.yaml` templates, `tests/` with
synthetic fixtures, its own `.gitignore` (patterns are relative to the file's
directory — keeping them beside the pipeline is what stops a future move from
unignoring real configs), and a `README.md`.

Shared machinery lives in [`common/lib_common.sh`](common/lib_common.sh) — the
YAML config reader, duplicate-key guard, `#`-argument guard, relative/absolute
path resolution, the `${sample.raw_sample_prefix}` substitution, the tool
checker, and reference-config resolution. **Source it; do not copy it.**

Each `00_setup_env.sh` starts with a short bootstrap that locates the suite root
and sources the library — that block cannot itself be shared, since finding the
library is what it does. Copy it verbatim from an existing pipeline. Everything
after it should be only your pipeline's own config keys, tool list and output
directories.

This matters because the two original pipelines each had their own copy and the
copies diverged: the zsh-safe sourcing fix, the `#`-argument guard and the
duplicate-key guard existed in one and not the other, so the SV pipeline could
not be sourced from a macOS zsh prompt for weeks after that was fixed elsewhere.

---

## Licence

Contributions are accepted under the MIT licence of this repository.
