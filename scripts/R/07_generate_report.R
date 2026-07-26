#!/usr/bin/env Rscript
# =============================================================================
# 07_generate_report.R
#
# Generates a set of summary figures giving an overview of all SV types,
# from the tables produced by 06_summary_stats.R. PNGs only (no PDF/HTML
# report bundling) to keep this dependency-light and each figure directly
# reusable in a manuscript.
#
# Usage:
#   Rscript scripts/R/07_generate_report.R [path/to/pipeline_config.yaml]
#
# Input:  results/qc_summary/<sample_id>.all_variants_combined.tsv
#         results/qc_summary/<sample_id>.summary_statistics.tsv
#         results/qc_summary/<sample_id>.fusion_type_breakdown.tsv
# Output: results/qc_summary/figures/<sample_id>.<figure_name>.png
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path("config", "pipeline_config.yaml")
if (!file.exists(config_path)) stop("[FATAL] Config file not found: ", config_path)

cat("=== [07_generate_report.R] Using config:", config_path, "===\n")
cfg <- yaml.load_file(config_path)
sample_id  <- cfg$sample$sample_id
output_dir <- cfg$paths$output_dir
cat("Sample ID:", sample_id, "\n")

qc_dir <- file.path(output_dir, "qc_summary")
combined_path <- file.path(qc_dir, paste0(sample_id, ".all_variants_combined.tsv"))
if (!file.exists(combined_path)) {
  stop("[FATAL] ", combined_path, " not found. Run scripts/R/06_summary_stats.R first.")
}
combined <- read.delim(combined_path, stringsAsFactors = FALSE)
combined$SVLEN <- suppressWarnings(as.numeric(combined$SVLEN))
combined$VAF   <- suppressWarnings(as.numeric(combined$VAF))
combined$abs_svlen <- abs(combined$SVLEN)

fig_dir <- file.path(qc_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Consistent category order/colors across all figures (by descending count,
# but with translocations grouped near gene_fusions/disruptions conceptually)
cat_order <- names(sort(table(combined$category), decreasing = TRUE))
combined$category <- factor(combined$category, levels = cat_order)

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

save_fig <- function(plot, name, width = 7, height = 5) {
  path <- file.path(fig_dir, paste0(sample_id, ".", name, ".png"))
  ggsave(path, plot, width = width, height = height, dpi = 300, bg = "white")
  cat("  Saved:", path, "\n")
}

cat("\nGenerating figures...\n")

# --- Figure 1: SV counts by category (log scale — deletions/insertions ------
# typically dwarf everything else) ------------------------------------------
counts_df <- as.data.frame(table(combined$category))
names(counts_df) <- c("category", "n")
counts_df$category <- factor(counts_df$category, levels = cat_order)

p1 <- ggplot(counts_df, aes(x = category, y = n, fill = category)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
  scale_y_log10(labels = scales::comma) +
  labs(title = paste0("Structural Variant Counts by Category — ", sample_id),
       x = NULL, y = "Number of variants (log10 scale)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p1, "01_sv_counts_by_category")

# --- Figure 2: SV length distribution per category (excludes translocations,-
# which have no SVLEN) --------------------------------------------------------
len_df <- combined[!is.na(combined$abs_svlen) & combined$abs_svlen > 0, ]
if (nrow(len_df) > 0) {
  p2 <- ggplot(len_df, aes(x = category, y = abs_svlen, fill = category)) +
    geom_violin(scale = "width", alpha = 0.6) +
    geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.8) +
    scale_y_log10(labels = scales::comma) +
    labs(title = paste0("SV Length Distribution by Category — ", sample_id),
         x = NULL, y = "SV length, bp (log10 scale)") +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p2, "02_sv_length_distribution")
} else {
  cat("  Skipped SV length distribution — no records with numeric SVLEN.\n")
}

# --- Figure 3: VAF distribution, faceted by category ------------------------
vaf_df <- combined[!is.na(combined$VAF), ]
if (nrow(vaf_df) > 0) {
  p3 <- ggplot(vaf_df, aes(x = VAF, fill = category)) +
    geom_histogram(bins = 30, boundary = 0) +
    facet_wrap(~category, scales = "free_y") +
    labs(title = paste0("Variant Allele Fraction Distribution — ", sample_id),
         x = "VAF", y = "Count") +
    theme_pub +
    theme(strip.text = element_text(face = "bold"))
  save_fig(p3, "03_vaf_distribution", width = 9, height = 6)
} else {
  cat("  Skipped VAF distribution — no VAF values present.\n")
}

# --- Figure 4: Genome-wide distribution by chromosome, stacked by category --
chrom_order <- paste0("chr", c(1:22, "X", "Y", "M"))
combined$CHROM <- factor(combined$CHROM, levels = chrom_order[chrom_order %in% combined$CHROM])
chrom_df <- combined[!is.na(combined$CHROM), ]
if (nrow(chrom_df) > 0) {
  p4 <- ggplot(chrom_df, aes(x = CHROM, fill = category)) +
    geom_bar() +
    labs(title = paste0("SV Distribution Across Chromosomes — ", sample_id),
         x = NULL, y = "Number of variants") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "right",
          panel.grid.minor = element_blank())
  save_fig(p4, "04_chromosome_distribution", width = 10, height = 5)
}

# --- Figure 5: Translocation breakpoint effect breakdown --------------------
fusion_breakdown_path <- file.path(qc_dir, paste0(sample_id, ".fusion_type_breakdown.tsv"))
if (file.exists(fusion_breakdown_path)) {
  fusion_df <- read.delim(fusion_breakdown_path, stringsAsFactors = FALSE)
  if (nrow(fusion_df) > 0) {
    p5 <- ggplot(fusion_df, aes(x = reorder(fusion_type, -n), y = n, fill = fusion_type)) +
      geom_col() +
      geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
      labs(title = paste0("Translocation Breakpoint Classification — ", sample_id),
           subtitle = "SnpEff-derived effect per BND record",
           x = NULL, y = "Number of BND records") +
      theme_pub
    save_fig(p5, "05_translocation_breakdown")
  }
} else {
  cat("  Skipped translocation breakdown — file not found (run 05_annotate_variants.R first).\n")
}

cat("\n=== [07_generate_report.R] Done. Figures written to:", fig_dir, "===\n")
