# Methods — detailed rationale

This document records the *why* behind design choices in this pipeline, for
reuse in a manuscript methods section and for future maintainers (including
future-you) who wonder why a script does something a particular way.

## Input data

Structural variants were called from Oxford Nanopore long-read whole-genome
sequencing data using Epi2ME Labs' `wf-human-variation` workflow, which uses
Sniffles2 for SV calling. If SnpEff annotation was enabled for the Epi2ME
run (`annotation: True`), the resulting VCF (`<prefix>.wf_sv.vcf.gz`)
includes an `ANN` INFO field per record with SnpEff's functional annotation.

## Why VCF, not TSV

Epi2ME's `wf-human-variation` does not emit a native TSV of individual SV
calls — SVs come out as a VCF. This pipeline flattens that VCF into TSV
early (`02_vcf_to_tsv.sh`) so that downstream filtering/categorization can
use simple, auditable `awk`/shell logic, while keeping the VCF as the single
source of truth for anything requiring full annotation fidelity (see gene
fusion calling, below).

## Filtering thresholds (`config/pipeline_config.yaml: filtering`)

Default thresholds (`min_sv_length=30`, `min_read_support=4`,
`min_qual=20`, `min_vaf=0.10`, `pass_only=true`) are intentionally
permissive — they largely restate what Epi2ME's own Sniffles2 invocation
already enforces before writing `PASS` records, rather than imposing
additional stringency. On the validation sample used during pipeline
development, these thresholds removed 0 of 33,419 records — this is
expected, not a bug, and is left in the config as an explicit, auditable
statement of QC criteria even where redundant with upstream filtering.

For a **germline** analysis, we recommend raising `min_vaf` to ~0.20–0.25,
which cleanly separates the expected heterozygous (~0.5) and homozygous
(~1.0) allele-fraction clusters from alignment noise without cutting into
either. For a study specifically interested in **somatic/mosaic** events, a
low `min_vaf` (~0.10 or lower) should be retained deliberately, since a
higher threshold would discard exactly the low-allele-fraction signal of
interest. This choice should be stated explicitly in any resulting
manuscript's methods section along with the value used.

`SVLEN` in the Sniffles2 VCF is signed (negative for deletions); filtering
uses absolute value throughout.

## SV categorization

SV type -> output category mapping (`config/pipeline_config.yaml:
sv_categories`, mirrored in `03_filter_sv_categories.sh`'s `CATEGORY_MAP`):

| Category | Sniffles2 SVTYPE(s) |
|---|---|
| deletions | DEL |
| insertions | INS |
| duplications | DUP |
| inversions | INV |
| translocations | BND, TRA |
| complex_rearrangements | CPX, BND_CLUSTER (rarely emitted by Sniffles2; included for forward compatibility) |
| gene_fusions | derived, not a native SVTYPE — see below |

## Gene fusion / disruption calling

**Design choice:** this pipeline trusts SnpEff's own breakpoint-effect
classification rather than independently pairing `BND` mate coordinates and
intersecting them against a gene BED from scratch.

Rationale: Epi2ME's SnpEff annotation, when enabled, already classifies each
`BND` breakpoint's functional effect in the `ANN` INFO field's second
pipe-delimited subfield (`Annotation`). Observed values relevant to fusion
calling, from real pipeline validation data:

- `bidirectional_gene_fusion` — two overlapping/adjacent genes disrupted by
  the same breakpoint; gene pair given in the `Gene_Name` subfield,
  ampersand-joined (e.g. `GENE1&GENE2`).
- `feature_fusion` — breakpoint joins two features; when the affected region
  is intergenic, the gene pair instead appears in the `Feature_ID` subfield,
  hyphen-joined (e.g. `GENE1-GENE2`), with `Gene_Name` left empty. **A
  single `ANN` record can carry multiple comma-separated annotation
  blocks, and the usable gene pair is not guaranteed to be in the first
  block** — during development, a real `feature_fusion` record's first
  block had empty `Gene_Name`/`Feature_ID` fields entirely, with the actual
  gene pair (`ZSWIM5-LOC105378691`) only present in the second block. The
  parser in `05_annotate_variants.R` therefore scans all comma-separated
  blocks rather than assuming the first is authoritative.
- `transcript_ablation` — single gene disrupted, not a fusion; reported
  separately in `gene_disruptions.tsv`.

Records whose `ANN` blocks don't match either pattern are reported as
`other` (annotated but not a recognized fusion/disruption effect) or
`unannotated` (no `ANN` present — e.g. if the Epi2ME run had
`annotation: False`) in `fusion_annotation_full.tsv`, but are not written to
either final `gene_fusions.tsv` or `gene_disruptions.tsv`.

## Duplicate/reciprocal fusion calls

Sniffles2 can emit more than one `BND` record implicating the same
underlying fusion event — for example, two breakpoints a few kilobases
apart that both disrupt the same gene pair, likely representing the same
structural event called from reads supporting slightly different breakpoint
estimates, or true reciprocal breakpoints of a balanced rearrangement.
`05_annotate_variants.R` flags rows sharing an unordered gene pair as
`possible_duplicate_event=TRUE` in `gene_fusions.tsv`, but does **not**
automatically collapse or deduplicate them — this is deliberately left as a
manual review step, since collapsing without inspecting breakpoint distance
and orientation risks discarding a real second event.

## Reference genome and annotation

Reference genome: GRCh38 (`GCA_000001405.15_GRCh38_no_alt_analysis_set`),
matching Epi2ME's default alignment reference. Gene annotation: GENCODE
(version recorded in `config/reference_paths.yaml`; update this document
if the reference bundle changes).

## Tool versions

Exact tool versions used for a given pipeline run are recorded fresh in
`logs/tool_versions.txt` on every execution of `00_setup_env.sh` (i.e. at
the start of every script). Copy the relevant lines into a manuscript's
methods section rather than restating versions here, since this file is
static documentation and the log is the authoritative per-run record.
