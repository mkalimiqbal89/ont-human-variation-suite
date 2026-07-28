# Troubleshooting

Real issues encountered during pipeline development, kept here so they don't
get rediscovered the hard way twice.

## "command not found" / tool check fails, but I definitely installed it

Check your prompt for which conda environment is active. A tool being
installed *somewhere* on the system is not the same as it being on `PATH`
in your *current shell*. `00_setup_env.sh` checks the tools reachable from
wherever it's run — if you see `(base)` in your prompt instead of your
actual tool environment, run `conda activate <your_env_name>` first.

If a tool is "found" but errors with something like `error while loading
shared libraries: libgsl.so.25: cannot open shared object file`, that's
still an environment mismatch, not a missing install — the binary on `PATH`
belongs to a different environment than the one with its matching
libraries. Same fix: activate the right environment.

## Config values resolve to empty strings or literal `${sample.X}` text

`pipeline_config.yaml` is not parsed by a real YAML library in the bash
scripts — `00_setup_env.sh` uses a lightweight `grep`/`sed` parser to avoid
requiring `yq` on the HPC. This means:

- Every field this pipeline reads must be a simple `key: "value"` line —
  nested structures beyond one level, multi-line values, or YAML anchors
  will not resolve correctly.
- `${sample.sample_id}` and `${sample.raw_sample_prefix}` placeholders in
  `input_files:` are resolved by explicit string substitution in
  `00_setup_env.sh` (`resolve_filename()`), not by a general templating
  engine. If you add a new placeholder pattern, you must add a
  corresponding substitution line, or it will be left as literal text.

If a value resolves empty, check that the key exists **inside the correct
top-level block** — a `grep -q <key>` check against the whole file can give
a false positive if the key name also appears as a placeholder reference
elsewhere in the file (this bit us once: `raw_sample_prefix` appeared inside
`${sample.raw_sample_prefix}` template strings in `input_files:` well before
the actual `raw_sample_prefix:` field existed under `sample:`, so a naive
existence check wrongly concluded the field was already present).

Diagnose with:
```bash
grep -n "your_key" config/pipeline_config.yaml | cat -A
```
`cat -A` shows tab/carriage-return characters (`^I`, `^M`) that plain
`grep`/`cat` render invisibly and can silently break the anchored regex
patterns the parser relies on.

## Terminal output "looks" wrong but the file is actually fine

Copying tab-delimited output out of an SSH terminal session can visually
collapse a tab to zero spaces if a value happens to end exactly on an
8-column tab stop, making two adjacent fields look glued together (e.g.
`5OR4F5` instead of `5` and `OR4F5`). Before assuming a real bug, check the
actual bytes in the file rather than trusting terminal-rendered/copy-pasted
text:
```bash
head -2 <file> | cat -A          # ^I between every field = real tabs, fine
awk -F'\t' '{print NF}' <file> | sort | uniq -c   # every row should report the same field count
```

## `01_validate_inputs.sh` is slow

`bcftools view -H <vcf> | wc -l` is run per input VCF to report a variant
count, which means fully decompressing and streaming the entire file —
on a ~5.4M-record SNV VCF this can take several minutes. This is a real
cost, not a bug; if it becomes a bottleneck on repeat runs, consider
skipping the SNV/CNV variant-count step (this pipeline doesn't process
those VCFs downstream anyway) or caching the count after the first
successful validation.

## Zero fusion/disruption calls, everything reports "unannotated"

Check whether the Epi2ME run was invoked with SnpEff annotation enabled
(`annotation: True`). If not, the VCF's `BND` records have no `ANN` INFO
field at all, and `05_annotate_variants.R` correctly has nothing to
classify — this is expected behavior given the input, not a script bug.
Confirm with:
```bash
bcftools view -h "$SV_VCF" | grep "^##INFO=<ID=ANN"
```
Missing = not annotated; re-run Epi2ME with annotation enabled, or use a
separate SV annotation tool.

## `05_annotate_variants.R` fails with an R package error

`library(yaml)` must succeed. Check with:
```bash
Rscript -e 'library(yaml)'
```
If missing, install via `Rscript -e 'install.packages("yaml",
repos="https://cloud.r-project.org")'` or, on an HPC with restricted
internet access, `conda install -c conda-forge r-yaml` inside your
pipeline's conda environment.

## A step produced zero variants for a category I expected to have some

Check `results/qc_summary/<sample_id>.filtering_summary.tsv` first — it
records the per-category count after filtering. Compare against the
pre-filter total in `data/processed/<sample_id>.sv_flat.tsv`
(`awk -F'\t' '$6=="DUP"' <file> | wc -l` for a specific SVTYPE, for
example) to see whether the category is genuinely empty in the source VCF,
or whether a filtering threshold in `pipeline_config.yaml` is responsible.

## Transferring files from this chat session to the HPC

Any file shared in this chat downloads to your **local machine**, not the
HPC — there is no direct connection between them. Standard flow:
```bash
# from a local terminal (not your SSH session)
scp downloaded_file.tar.gz hpc@<host>:/path/on/hpc/
# then, in your SSH session
tar -xzf downloaded_file.tar.gz
md5sum <extracted_file>   # compare against the checksum provided alongside the download
```
Always verify the checksum after transfer — a truncated or corrupted
transfer can otherwise fail silently (a partially-written script may still
be syntactically valid bash/R and run without obvious errors while
producing wrong output).
