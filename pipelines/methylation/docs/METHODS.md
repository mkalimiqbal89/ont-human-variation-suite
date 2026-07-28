# Methods — detailed rationale

Why the methylation pipeline makes the choices it does. Written to be the source
for a manuscript methods section, and to record the reasoning behind thresholds
that would otherwise look arbitrary.

Companion documents: [`COMPARISON_CAVEATS.md`](COMPARISON_CAVEATS.md) for the
statistics of cross-sample comparison, and the pipeline [`README`](../README.md)
for usage.

---

## Input data

Oxford Nanopore whole-genome sequencing processed by the Epi2ME
`wf-human-variation` workflow (v2.7.2), which calls modified bases with `modkit`
and emits a bedMethyl pileup: `<prefix>.wf_mods.bedmethyl.gz`.

This pipeline starts from that pileup. It does not call modified bases, and
deliberately never touches the BAM/CRAM — the per-read modification probabilities
have already been summarised into per-position counts, and re-deriving them would
duplicate `modkit` for no benefit.

### bedMethyl column layout

`modkit pileup` emits 18 tab-separated columns. Four carry all the signal this
pipeline uses:

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
| 14–18 | other-mod, delete, fail, diff and no-call counts |

Every stage asserts the column count on every record rather than trusting it.
A file with a different layout — a different `modkit` version, a hand-edited
file — fails immediately rather than producing plausible nonsense.

---

## Modification code selection

`wf-human-variation` reports `m` (5mC) and `h` (5hmC) interleaved at each
position: two records per CpG. This pipeline analyses one code at a time,
configured by `modifications.primary_mod_code` (default `m`).

5mC is the default because it is the dominant and best-characterised cytosine
modification, and because ONT's 5hmC calling carries greater uncertainty at the
coverage typical of a clinical WGS run. Setting `primary_mod_code: h` runs the
identical pipeline over 5hmC; nothing else needs to change.

Both samples analysed here contained exactly equal numbers of `m` and `h`
records (28,700,314 and 28,630,840 respectively, half of each file), consistent
with one record per code per covered CpG.

---

## Contig selection: an allowlist, not a blocklist

`filtering.primary_contigs_only: true` keeps only contigs matching
`include_contigs_regex` — by default `^chr([1-9]|1[0-9]|2[0-2]|X|Y)$`, i.e.
chr1–22, chrX, chrY.

**Why an allowlist.** A blocklist can only reject the naming patterns you thought
to enumerate. GRCh38 analysis sets ship scaffolds named `GL000191.1`,
`HLA-A*01:01:01:01` and `chr1_KI270762v1_fix`, none of which match the usual
`_random$` / `^chrUn` / `_alt$` patterns. A regression test asserts this
explicitly: on a fixture containing four such scaffolds, the blocklist rejects
one and the allowlist rejects all four.

**chrM is deliberately excluded.** Mammalian mitochondrial DNA is essentially
unmethylated, and mtDNA carries far higher coverage than the nuclear genome, so
including it drags any coverage-weighted genome-wide average downward. In the
demo dataset chrM contributed 435 sites at 71,466 total coverage and 0.17%
methylation. Add `|M` inside the regex group to include it.

**This choice made two samples comparable that otherwise were not.** The two
validation samples differ markedly in what their bedMethyl files contain:

| Sample | Contigs excluded | Records excluded |
|---|---:|---:|
| HLH_S0001 | 127 (126 scaffolds + chrM) | 135,719 |
| HLH_S0002 | 1 (chrM only) | 435 |

HLH_S0002's pileup contains no scaffold contigs at all — the two Epi2ME runs
evidently used different reference contig sets, or the second file was filtered
before archiving. Because both are reduced to the same 24 primary contigs, the
downstream analysis is unaffected and no special handling was required. Under a
blocklist tuned to one sample's naming, this would have been a silent
inconsistency between them.

A warning (never a failure) is raised if the retained contig count differs from
`expected_primary_contigs`. 23 rather than 24 normally means chrY carried no
calls, which is expected for a female sample and should be checked against
recorded sex rather than assumed.

---

## Coverage threshold

`filtering.min_coverage: 10`, applied to bedMethyl column 10 (valid coverage).

The threshold is justified from the data rather than convention. Stage 02 emits
a full coverage histogram before any filtering, and prints the fraction of sites
retained at a range of thresholds. For HLH_S0001:

| Threshold | Sites retained | Share |
|---|---:|---:|
| ≥1 | 28,564,595 | 100.00% |
| ≥5 | 28,317,025 | 99.13% |
| **≥10** | **27,869,369** | **97.57%** |
| ≥15 | 26,011,981 | 91.06% |
| ≥20 | 19,863,112 | 69.54% |
| ≥30 | 3,430,218 | 12.01% |

At mean coverage 23.8x, a threshold of 10 discards 2.4% of sites while removing
positions where a single read flips the estimate by 10 percentage points or
more. Moving to 20 would discard 30% of the genome for little gain in precision.

`max_coverage: 0` disables an upper cap. A cap is available for pileup artefacts
but was not needed here; the observed maximum was 7,138x, concentrated in
repetitive regions that the contig allowlist already removes.

Applied to the two samples:

| Sample | Sites before | Dropped | Retained | Share |
|---|---:|---:|---:|---:|
| HLH_S0001 | 28,564,595 | 695,226 | 27,869,369 | 97.57% |
| HLH_S0002 | 28,630,405 | 539,496 | 28,090,909 | 98.12% |

---

## Methylation state bins

Stage 03 splits retained sites into three states, each written to its own
tabix-indexed file:

| State | Rule |
|---|---|
| unmethylated | percent < `unmethylated_max_percent` (20) |
| intermediate | between the thresholds |
| methylated | percent ≥ `methylated_min_percent` (80) |

20/80 is the conventional split. The intermediate class is the analytically
interesting one: at adequate coverage a CpG near 50% is a candidate for
allele-specific methylation or imprinting rather than noise, and separating it
means those candidates are not buried among 28 million sites.

Observed distribution is near-identical between samples, as expected for two
blood-derived genomes:

| Sample | Unmethylated | Intermediate | Methylated |
|---|---:|---:|---:|
| HLH_S0001 | 8.95% | 14.08% | 76.97% |
| HLH_S0002 | 8.88% | 15.02% | 76.10% |

Boundaries are inclusive/exclusive as stated and are pinned by regression tests:
a site at exactly 20.00% is intermediate, one at exactly 80.00% is methylated.

---

## Annotation and aggregation keys

Stage 05 intersects retained CpGs with three feature sets using
`bedtools intersect -wa -wb`, then aggregates per feature.

**Both aggregation keys were chosen by measuring the reference files, not by
assumption.**

**Genes and promoters are keyed on `gene_id`, never `gene_name`.** GENCODE v50
`genes.bed` holds 78,733 rows but only 77,118 distinct gene names — `Y_RNA`
alone appears 756 times. Keying on name would silently merge unrelated loci.

**CpG islands are keyed on `chrom:start-end`, never the island name.**
`cpg_islands_hg38.bed` holds 32,038 islands but only **485 distinct names**,
because the UCSC `name` field is `CpG:_<cpgNum>` — `CpG:_21` appears 700 times.
Keying on name would collapse the genome to 485 rows, and every downstream number
would still look plausible.

Both hazards have dedicated regression tests using fixtures with deliberate name
collisions.

**A CpG overlapping two features is counted in both.** Genes overlap, so the sum
of per-feature CpG counts legitimately exceeds the number of unique CpGs, and
these totals do not reconcile with the genome-wide figures. They answer different
questions. A reconciliation assertion checks the weaker but meaningful identity:
sum of per-feature CpG counts equals the number of intersect output rows.

**Features with fewer than `min_cpgs_per_feature` (5) covered CpGs are flagged,
not dropped.** A locus backed by one or two CpGs reads 0% or 100% by chance. The
flag lets downstream stages exclude them while keeping the rows auditable.

Annotation coverage for HLH_S0001:

| Feature class | Features with data | Meeting ≥5 CpGs | Unique CpGs annotated |
|---|---:|---:|---:|
| gene body | 74,660 | 65,115 | 20,897,837 |
| promoter | 77,523 | 77,465 | 4,689,992 |
| CpG island | 27,560 | 27,540 | 2,056,367 |

---

## Summary statistics: weighted vs unweighted

Stage 06 reports both, and they differ substantially. **The coverage-weighted
figure is the one to quote:**

```
weighted   = total modified read calls / total valid coverage × 100
unweighted = mean of per-site percentages
```

The unweighted mean treats a site with 10 reads and one with 200 as equally
informative. The weighted figure gives each read one vote, which is the estimate
with the smaller variance. Both are reported so the difference is visible rather
than hidden by a choice.

Genome-wide results:

| Sample | Retained sites | Mean coverage | Unweighted | **Weighted** |
|---|---:|---:|---:|---:|
| HLH_S0001 | 27,869,369 | 23.80x | 81.18% | **81.88%** |
| HLH_S0002 | 28,090,909 | 27.46x | 80.30% | **81.77%** |

The per-site methylation distribution is strongly bimodal — 58.6% of sites above
90% and 8.0% below 10% for HLH_S0001 — which is the canonical shape for a
mammalian methylome and serves as a global sanity check on the whole pipeline.

Feature-class comparison for HLH_S0001, restricted to features meeting the CpG
minimum:

| Class | Features | Mean | Median | Q1 | Q3 |
|---|---:|---:|---:|---:|---:|
| genome-wide (all sites) | 27.9 M | 81.88 | — | — | — |
| gene body | 65,115 | 76.02 | 85.25 | 71.31 | 90.91 |
| promoter | 77,465 | 66.82 | 84.31 | 36.63 | 91.69 |
| CpG island | 27,540 | 29.74 | **2.72** | 0.21 | 81.80 |

CpG islands at a median of 2.72% against a genome background of 81.88% is the
expected biology and is not a subtle effect. The promoter distribution is
strongly bimodal (Q1 36.63, median 84.31), because GENCODE's ~78,000 genes are
mostly non-coding; stage 06 therefore also emits a per-`gene_type` breakdown, and
promotes the protein-coding subsets into the main table.

---

## Cross-sample comparison

Stages 09–10 compare two or more samples at gene, promoter and island level,
working from the per-feature tables rather than site-level data.

Two filters gate a `differential` call, and both must pass:

- **Effect size** — `min_delta_percent: 20` percentage points.
- **Significance** — BH-corrected q < `fdr_alpha: 0.05`.

Plus a **coverage-balance** requirement, added after inspecting the first real
two-sample run: `max_cpg_ratio: 3` and `max_coverage_ratio: 5`. `min_shared_cpgs`
requires enough covered CpGs in each sample but not *comparable* numbers, and the
initial top promoter hits had 10 covered CpGs in one sample against 173 in the
other — mappability differences in pseudogene and subtelomeric loci, not
methylation differences. Such features are reported with their ratios but
labelled `imbalanced`.

The filter removed 17% of gene calls, 27% of promoter calls and 2% of island
calls. That gradient is itself evidence it targets the right thing: CpG islands
are GC-rich unique sequence and map cleanly, whereas promoters frequently overlap
repeat-rich regions.

**The p-values are anti-conservative and must not be quoted without reading
[`COMPARISON_CAVEATS.md`](COMPARISON_CAVEATS.md).** Pooled-read tests treat every
read as independent; nanopore reads span many CpGs, and with one sample per group
there is no biological replication. Setting `stat_test: none` reports effect
sizes without p-values, which is a defensible choice for a single pair.

---

## Reference preparation

Built once and shared across all projects. Exact provenance:

| File | Source and derivation |
|---|---|
| `genes.bed` (BED7) | GENCODE v50 GTF, `gene` features only; GTF 1-based inclusive converted to BED 0-based half-open; fields: chrom, start, end, gene_id, gene_name, gene_type, strand |
| `promoters_2kb.bed` (BED7) | TSS taken strand-aware from `genes.bed` (gene start on `+`, gene end on `−`), extended ±2 kb with `bedtools slop` against `hg38.chrom.sizes` so coordinates cannot run off a chromosome |
| `cpg_islands_hg38.bed` (BED10) | UCSC `cpgIslandExt` table; leading internal index column removed; whitespace in island names replaced with underscores so shell tooling does not split them |
| `hg38.chrom.sizes` | UCSC goldenPath |

Row counts are validated at load: `genes.bed` and `promoters_2kb.bed` must be 1:1
(78,733 rows each), and a mismatch indicates one was regenerated without the
other.

**The genome FASTA is not required.** `modkit` has already produced the pileup,
so no sequence access is needed downstream.

Stage 01 verifies that contig naming is consistent across the bedMethyl and all
annotation files. A UCSC/Ensembl mismatch (`chr1` versus `1`) makes every
`bedtools intersect` return zero rows, which presents as "this sample has no
methylation in genes" rather than as an error.

---

## Tool versions

Recorded automatically to `logs/tool_versions.txt` on every run, and to the
provenance manifest with checksums of the config and every annotation file. The
validation runs used:

```
bedtools   v2.31.1
tabix      (htslib) 1.24
bgzip      (htslib) 1.24
Rscript    R version 4.6.0 (2026-04-24)
awk        version 20200816 (BSD awk, macOS)
sort       2.3-Apple (199)
```

R is used only for figures (stage 07) and cross-sample statistics (stage 09), and
requires **no R packages** — base graphics and `stats` only. Stages 02–06 are awk,
which streams ~28 M records in constant memory; the same work in R would require
several GB.

---

## What this pipeline does not do

- **Call modified bases.** That is `modkit`'s job, upstream.
- **Haplotype-resolved methylation.** Stage 02 refuses `phased: true` rather than
  half-supporting it. Both validation samples were unphased for modifications.
- **Differential methylation with replicates.** The comparison stage handles two
  or more samples but has no model for between-sample variance. With replicates,
  use `DSS`, `methylKit` or `modkit dmr` instead.
- **Region calling (DMR detection).** Comparison is per annotated feature, not
  per empirically-defined region.
