# Interpreting the cross-sample comparison

Stage 09 will happily produce a table of "significantly differentially methylated"
features with very small q-values. This document explains why most of those
p-values should not be quoted, and what the output *is* good for.

Read this before putting a q-value in a manuscript.

---

## What the test actually does

For each feature, stage 09 pools reads across all CpGs in that feature and builds
a 2×N table of modified versus unmodified read calls:

```
                sample A    sample B
modified          1,842       2,910
unmodified        4,108         690
```

It then applies Fisher's exact test (two samples) or a chi-squared test (more
than two), and corrects across features with Benjamini–Hochberg.

## Why the p-values are anti-conservative

Three separate problems, each of which inflates significance on its own.

**1. Reads are not independent observations.** The test assumes every read call
is an independent Bernoulli trial. A nanopore read is tens of kilobases long and
covers many CpGs in the same feature, so the calls within a feature are strongly
correlated — they came from the same few molecules. The effective sample size is
closer to the number of *reads* spanning the feature than the number of read
calls, and often closer still to the number of distinct alleles. Treating
correlated observations as independent is the classic route to a p-value that is
too small by orders of magnitude.

**2. There is no biological replication.** With one sample per group, the test
compares two individuals. Any difference it finds is confounded with everything
that differs between those two people: genotype, age, sex, cell-type composition
of the blood draw, time since collection, library prep batch. The test cannot
distinguish a disease-associated difference from ordinary between-individual
variation, because it has no estimate of the latter.

**3. Enormous read counts make trivial differences significant.** At ~24x
coverage across 27.9 M CpGs, a feature can easily carry several thousand read
calls. A difference of one or two percentage points — well inside technical
variation for nanopore modified-base calling — will clear any conventional
significance threshold. This is why stage 09 reports
`significant_only` as its own category: those are features where the statistics
say yes and the biology almost certainly says no.

## Coverage imbalance: the fourth problem, found on real data

The first genuine two-sample run (HLH_S0001 vs HLH_S0002) produced this as its
top promoter hit:

```
ENSG00000290383   HLH_S0001:  10 covered CpGs, coverage   100
                  HLH_S0002: 173 covered CpGs, coverage  2514
                  delta -87.5 pp,  q = 1.7e-86
```

Also in the top ten: `PKD1P5` (a pseudogene, identical numbers because it
overlaps the same locus), `MIR3180-3`, and `FOXD4L5` — all pseudogene,
subtelomeric or segmental-duplication loci where long-read placement is
unstable.

Nothing about that is a methylation difference. Two problems compound:

1. **It is a mappability difference.** One sample placed 25× more read coverage
   in the region than the other. That is a property of the repeat structure and
   the aligner, not of the epigenome.
2. **The two values do not summarise the same positions.** With 10 covered CpGs
   versus 173, the "methylation" of that promoter is computed over two different
   subsets of the feature. The delta is not one quantity measured twice.

`min_shared_cpgs` does not catch this: it requires *enough* CpGs in each sample,
not *comparable* numbers. Stage 09 therefore also computes, per feature:

- `cpg_ratio` — max/min covered-CpG count across samples
- `coverage_ratio` — max/min total coverage across samples

and sets `balanced = no` when either exceeds `max_cpg_ratio` (default 3) or
`max_coverage_ratio` (default 5). Such features keep their row and their ratios,
but are labelled `imbalanced` in the `differential` column rather than `yes`, are
sorted below balanced features, and are drawn in a distinct colour on the
scatter plot.

**If your `imbalanced_excluded` count is large, that is informative in itself**:
it says the two libraries differ in where they place reads, which is worth
understanding before interpreting anything regional.

## What the output is good for

**Effect sizes.** `delta` (two samples) or `range` and `sd` (more than two) are
descriptive statistics, not inferences. They mean what they say: this feature's
coverage-weighted methylation differs by this many percentage points. The
`min_delta_percent` threshold (default 20 pp) is deliberately well above
technical noise.

**Ranking candidates.** Combining a large effect with a small q-value is a
reasonable way to order features for follow-up. The `differential` column
requires both, which is stricter than either alone.

**Global and per-class comparison.** The genome-wide and feature-class figures
compare millions of sites and are robust to all three problems above. If two
samples differ in genome-wide methylation, or in the promoter-versus-gene-body
relationship, that is a real observation about those samples.

**Quality control.** The scatter plot and Pearson correlation across shared
features are an excellent check that two runs are comparable at all. A
correlation far below what you expect usually means a technical problem, not
biology.

## What would make the statistics defensible

In rough order of how much they help:

1. **Biological replicates.** Several samples per group, analysed with a method
   that models between-sample variance — `DSS`, `methylKit`, or a
   beta-binomial mixed model. This addresses problem 2 and partly problem 1.
2. **A read-level or region-level test.** `modkit dmr` operates on read-level
   data and handles within-region correlation better than pooling counts.
3. **Aggregating to regions with an explicit correlation structure**, rather than
   pooling reads as if independent.
4. **Permutation of sample labels** to calibrate the null — only meaningful once
   there are replicates to permute.

Until then, treat the per-feature p-values as a sorting key.

## Recommended wording

Defensible:

> Promoter methylation at *GENE* differed by 34 percentage points between the two
> samples (81.2% versus 47.1%, coverage-weighted, ≥5 covered CpGs in both).

Not defensible without replicates:

> *GENE* was significantly hypomethylated (q = 3.2 × 10⁻¹²).

If a reviewer asks whether the p-values account for read correlation and
between-individual variation, the honest answer is that they do not, and this
document is where that is written down.

---

## Configuration

| Key | Default | Effect |
|---|---|---|
| `min_shared_cpgs` | 5 | A feature is compared only where every sample has at least this many covered CpGs |
| `min_delta_percent` | 20 | Minimum absolute difference, in percentage points, for a `differential` call |
| `fdr_alpha` | 0.05 | BH-corrected threshold |
| `stat_test` | `fisher` | `fisher`, `chisq`, or `none` to report effect sizes only |
| `top_n_report` | 100 | Rows in the `*_top_differential.tsv` tables |

Setting `stat_test: none` is a reasonable choice for a single-pair comparison: it
produces effect sizes without inviting a p-value that cannot bear weight.
