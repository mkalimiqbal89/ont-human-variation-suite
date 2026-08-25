#!/usr/bin/env Rscript
# =============================================================================
# 07_generate_report.R (pipelines/cnv)
# Renders publication-quality figures using base R graphics (zero R package
# dependencies). Generates summary figures for CNV category counts, length
# distributions, chromosome spatial density, and copy-number state calls.
#
# Usage:
#   Rscript scripts/R/07_generate_report.R [path/to/pipeline_config.yaml]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg_arg <- if (length(args) > 0) args[1] else ""

get_yaml_val <- function(file, key) {
  cmd <- sprintf("grep -E '^[[:space:]]*%s:' '%s' | head -n1 | sed -E 's/^[^:]+:[[:space:]]*\"?//; s/\"?[[:space:]]*$//'", key, file)
  val <- suppressWarnings(system(cmd, intern = TRUE))
  if (length(val) == 0) return("")
  return(trimws(val))
}

script_dir <- dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))
if (length(script_dir) == 0) script_dir <- "scripts/R"
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = FALSE)

config_file <- if (nchar(cfg_arg) > 0) cfg_arg else file.path(repo_dir, "config/pipeline_config.yaml")
if (!file.exists(config_file)) {
  stop(sprintf("Config file not found: %s", config_file))
}

cat(sprintf("=== [07_generate_report.R] Using config: %s ===\n", config_file))

sample_id <- get_yaml_val(config_file, "sample_id")
output_dir <- get_yaml_val(config_file, "output_dir")

cat(sprintf("Sample ID: %s\n", sample_id))

qc_dir <- file.path(output_dir, "qc_summary")
fig_dir <- file.path(qc_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

combined_file <- file.path(qc_dir, sprintf("%s.all_cnvs_combined.tsv", sample_id))
stats_file <- file.path(qc_dir, sprintf("%s.summary_statistics.tsv", sample_id))

combined <- if (file.exists(combined_file) && file.info(combined_file)$size > 0) read.delim(combined_file, stringsAsFactors = FALSE) else data.frame()
stats_df <- if (file.exists(stats_file) && file.info(stats_file)$size > 0) read.delim(stats_file, stringsAsFactors = FALSE) else data.frame()

# Color palette definition (base R)
col_del <- "#d95f02" # orange/red for deletion/loss
col_dup <- "#7570b3" # purple/blue for duplication/gain
col_main <- "#1b9e77"

# Figure 1: Counts & Span by Category
fig1 <- file.path(fig_dir, sprintf("%s.01_cnv_counts_by_category.png", sample_id))
png(fig1, width = 900, height = 500, res = 120)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
if (nrow(stats_df) > 0) {
  counts <- stats_df$n_variants
  names(counts) <- stats_df$category
  cols <- ifelse(names(counts) == "deletions", col_del, col_dup)
  bp <- barplot(counts, main = sprintf("CNV Counts (%s)", sample_id), ylab = "Number of Call Segments",
                col = cols, border = NA, ylim = c(0, max(counts, 1) * 1.2))
  text(bp, counts + max(counts, 1) * 0.05, labels = counts, pos = 3)

  span_mb <- stats_df$total_span_bp / 1e6
  names(span_mb) <- stats_df$category
  bp2 <- barplot(span_mb, main = "Genomic Span Impacted (Mb)", ylab = "Total Span (Megabases)",
                 col = cols, border = NA, ylim = c(0, max(span_mb, 0.1) * 1.2))
  text(bp2, span_mb + max(span_mb, 0.1) * 0.05, labels = sprintf("%.2f Mb", span_mb), pos = 3)
} else {
  plot.new()
  text(0.5, 0.5, "No CNVs detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig1))

# Figure 2: CNV Length Distribution
fig2 <- file.path(fig_dir, sprintf("%s.02_cnv_size_distribution.png", sample_id))
png(fig2, width = 800, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
if (nrow(combined) > 0 && "SVLEN" %in% colnames(combined)) {
  lens <- as.numeric(combined$SVLEN)
  lens <- lens[!is.na(lens) & lens > 0]
  if (length(lens) > 0) {
    log_lens <- log10(lens)
    hist(log_lens, main = sprintf("CNV Length Distribution (%s)", sample_id),
         xlab = "Log10 Segment Length (bp)", ylab = "Frequency",
         col = col_main, border = "white", breaks = 15)
    abline(v = median(log_lens), col = "red", lty = 2, lwd = 2)
    legend("topright", legend = sprintf("Median: %.0f bp", median(lens)), col = "red", lty = 2, bty = "n")
  } else {
    plot.new(); text(0.5, 0.5, "No valid lengths")
  }
} else {
  plot.new(); text(0.5, 0.5, "No CNVs detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig2))

# Figure 3: Chromosome Distribution
fig3 <- file.path(fig_dir, sprintf("%s.03_chromosome_cnv_distribution.png", sample_id))
png(fig3, width = 1000, height = 500, res = 120)
par(mar = c(6, 5, 4, 2))
if (nrow(combined) > 0 && "CHROM" %in% colnames(combined)) {
  chrom_counts <- table(combined$CHROM)
  bp <- barplot(chrom_counts, main = sprintf("CNVs per Chromosome (%s)", sample_id),
                ylab = "Number of CNVs", col = "#386cb0", border = NA, las = 2, cex.names = 0.8)
} else {
  plot.new(); text(0.5, 0.5, "No CNVs detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig3))

# Figure 4: Copy Number Distribution
fig4 <- file.path(fig_dir, sprintf("%s.04_copy_number_distribution.png", sample_id))
png(fig4, width = 700, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
if (nrow(combined) > 0 && "CN" %in% colnames(combined)) {
  cn_vals <- combined$CN[combined$CN != "NA" & combined$CN != "."]
  if (length(cn_vals) > 0) {
    cn_table <- table(cn_vals)
    bp <- barplot(cn_table, main = sprintf("Copy Number State Calls (%s)", sample_id),
                  xlab = "Copy Number (CN)", ylab = "Count", col = "#f0027f", border = NA)
  } else {
    plot.new(); text(0.5, 0.5, "No Copy Number values")
  }
} else {
  plot.new(); text(0.5, 0.5, "No CNVs detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig4))

cat(sprintf("=== [07_generate_report.R] Done. Figures written to: %s ===\n", fig_dir))
