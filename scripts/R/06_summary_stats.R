#!/usr/bin/env Rscript
# =============================================================================
# 06_summary_stats.R
#
# Aggregates all per-category SV TSVs (from 03_filter_sv_categories.sh) plus
# the gene fusion/disruption tables (from 05_annotate_variants.R) into:
#   1. A single combined long-format table (one row per variant, all
#      categories, with a `category` column) -- the input 07_generate_report.R
#      plots from.
#   2. A per-category summary statistics table (counts, SV length stats,
#      VAF stats, support stats) -- the table to drop into a manuscript.
#
# Usage:
#   Rscript scripts/R/06_summary_stats.R [path/to/pipeline_config.yaml]
#
# Input:  results/<category>/<sample_id>.<category>.tsv  (categories below)
#         results/gene_fusions/<sample_id>.gene_fusions.tsv
#         results/gene_fusions/<sample_id>.gene_disruptions.tsv
#         results/gene_fusions/<sample_id>.fusion_annotation_full.tsv
# Output: results/qc_summary/<sample_id>.all_variants_combined.tsv
#         results/qc_summary/<sample_id>.summary_statistics.tsv
# =============================================================================

suppressPackageStartupMessages(library(yaml))

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path("config", "pipeline_config.yaml")
if (!file.exists(config_path)) stop("[FATAL] Config file not found: ", config_path)

cat("=== [06_summary_stats.R] Using config:", config_path, "===\n")
cfg <- yaml.load_file(config_path)
sample_id  <- cfg$sample$sample_id
output_dir <- cfg$paths$output_dir
cat("Sample ID:", sample_id, "\n")

# --- Core SV categories (SVLEN-bearing types) --------------------------------
# translocations are handled separately below since BND records lack SVLEN.
core_categories <- c("deletions", "insertions", "duplications",
                      "inversions", "complex_rearrangements")

read_category <- function(category, output_dir, sample_id) {
  path <- file.path(output_dir, category, paste0(sample_id, ".", category, ".tsv"))
  if (!file.exists(path)) {
    warning("Missing expected file, skipping: ", path)
    return(NULL)
  }
  df <- tryCatch(read.delim(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$category <- category
  df
}

cat("Loading per-category tables...\n")
category_tables <- lapply(core_categories, read_category, output_dir = output_dir, sample_id = sample_id)
names(category_tables) <- core_categories

translocations_path <- file.path(output_dir, "translocations", paste0(sample_id, ".translocations.tsv"))
translocations <- if (file.exists(translocations_path)) {
  df <- read.delim(translocations_path, stringsAsFactors = FALSE)
  if (nrow(df) > 0) df$category <- "translocations"
  df
} else NULL

all_tables <- c(category_tables, list(translocations = translocations))
all_tables <- all_tables[!vapply(all_tables, is.null, logical(1))]

if (length(all_tables) == 0) stop("[FATAL] No category tables with data found. Run scripts/bash/03_filter_sv_categories.sh first.")

# --- Combine into one long table ---------------------------------------------
# Columns differ slightly (translocations lack SVLEN/END in a meaningful
# sense -- they're "." already from 02_vcf_to_tsv.sh). Keep the shared core
# columns needed for summary stats/plotting; anything category-specific
# (e.g. candidate_fusion) is preserved via rbind's NA-filling where absent.
shared_cols <- c("CHROM", "POS", "ID", "QUAL", "FILTER", "SVTYPE", "SVLEN",
                  "END", "CHR2", "SUPPORT", "STRAND", "VAF", "GT",
                  "GENE_SYMBOL", "ANNOTATION_IMPACT", "category")

combined <- do.call(rbind, lapply(all_tables, function(df) {
  missing_cols <- setdiff(shared_cols, names(df))
  for (mc in missing_cols) df[[mc]] <- NA
  df[, shared_cols]
}))
rownames(combined) <- NULL

cat("Combined table:", nrow(combined), "variants across", length(unique(combined$category)), "categories.\n")

# --- Gene fusion/disruption breakdown (subset of translocations) ------------
fusion_summary <- data.frame(
  fusion_type = character(0), n = integer(0)
)
full_bnd_path <- file.path(output_dir, "gene_fusions", paste0(sample_id, ".fusion_annotation_full.tsv"))
if (file.exists(full_bnd_path)) {
  bnd_full <- read.delim(full_bnd_path, stringsAsFactors = FALSE)
  if (nrow(bnd_full) > 0 && "fusion_type" %in% names(bnd_full)) {
    fusion_summary <- as.data.frame(table(bnd_full$fusion_type), stringsAsFactors = FALSE)
    names(fusion_summary) <- c("fusion_type", "n")
  }
}

# --- Per-category summary statistics -----------------------------------------
cat("Computing summary statistics per category...\n")

safe_stat <- function(x, fn) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  fn(x)
}

summary_rows <- lapply(split(combined, combined$category), function(df) {
  abs_svlen <- suppressWarnings(abs(as.numeric(df$SVLEN)))
  data.frame(
    category      = unique(df$category),
    n_variants    = nrow(df),
    mean_svlen_bp = round(mean(abs_svlen, na.rm = TRUE), 1),
    median_svlen_bp = round(median(abs_svlen, na.rm = TRUE), 1),
    min_svlen_bp  = suppressWarnings(min(abs_svlen, na.rm = TRUE)),
    max_svlen_bp  = suppressWarnings(max(abs_svlen, na.rm = TRUE)),
    mean_vaf      = round(safe_stat(df$VAF, mean), 3),
    median_vaf    = round(safe_stat(df$VAF, median), 3),
    mean_support  = round(safe_stat(df$SUPPORT, mean), 1),
    n_chromosomes = length(unique(df$CHROM)),
    stringsAsFactors = FALSE
  )
})
summary_stats <- do.call(rbind, summary_rows)
summary_stats <- summary_stats[order(-summary_stats$n_variants), ]

# Replace non-finite placeholders (e.g. Inf/-Inf from min/max on empty
# numeric vectors, which happens for translocations where SVLEN is always ".")
summary_stats[] <- lapply(summary_stats, function(col) {
  if (is.numeric(col)) { col[!is.finite(col)] <- NA }
  col
})

# --- Write outputs -------------------------------------------------------
out_dir <- file.path(output_dir, "qc_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

combined_out <- file.path(out_dir, paste0(sample_id, ".all_variants_combined.tsv"))
write.table(combined, combined_out, sep = "\t", row.names = FALSE, quote = FALSE)

stats_out <- file.path(out_dir, paste0(sample_id, ".summary_statistics.tsv"))
write.table(summary_stats, stats_out, sep = "\t", row.names = FALSE, quote = FALSE)

fusion_out <- file.path(out_dir, paste0(sample_id, ".fusion_type_breakdown.tsv"))
write.table(fusion_summary, fusion_out, sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== [06_summary_stats.R] Done ===\n")
cat("  Combined variant table :", nrow(combined), "rows ->", combined_out, "\n")
cat("  Summary statistics     :", nrow(summary_stats), "categories ->", stats_out, "\n")
cat("  Fusion type breakdown  :", nrow(fusion_summary), "types ->", fusion_out, "\n")
print(summary_stats, row.names = FALSE)
