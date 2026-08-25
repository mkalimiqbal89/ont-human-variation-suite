#!/usr/bin/env Rscript
# =============================================================================
# 07_generate_report.R (pipelines/snv)
# Renders publication-quality figures using base R graphics (zero R package
# dependencies). Generates summary figures for variant counts, VAF & depth
# distributions, Ti/Tv ratio, mutation spectrum, and chromosome variant density.
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

combined_file <- file.path(qc_dir, sprintf("%s.all_snvs_combined.tsv", sample_id))
stats_file <- file.path(qc_dir, sprintf("%s.summary_statistics.tsv", sample_id))

combined <- if (file.exists(combined_file) && file.info(combined_file)$size > 0) read.delim(combined_file, stringsAsFactors = FALSE) else data.frame()
stats_df <- if (file.exists(stats_file) && file.info(stats_file)$size > 0) read.delim(stats_file, stringsAsFactors = FALSE) else data.frame()

get_stat <- function(m) {
  val <- stats_df$value[stats_df$metric == m]
  if (length(val) == 0 || val == "NA") return(0)
  return(as.numeric(val))
}

# Color palette definition (base R)
col_snv <- "#4daf4a"
col_ins <- "#377eb8"
col_del <- "#e41a1c"
col_clinvar <- "#984ea3"

# Figure 1: Variant Class Breakdown
fig1 <- file.path(fig_dir, sprintf("%s.01_snv_indel_counts.png", sample_id))
png(fig1, width = 800, height = 500, res = 120)
par(mar = c(5, 5, 4, 2))
if (nrow(stats_df) > 0) {
  counts <- c(
    SNVs = get_stat("total_snvs"),
    Insertions = get_stat("indel_insertions"),
    Deletions = get_stat("indel_deletions"),
    ClinVar_Pathogenic = get_stat("clinvar_pathogenic_count")
  )
  cols <- c(col_snv, col_ins, col_del, col_clinvar)
  bp <- barplot(counts, main = sprintf("Small Variant Calls (%s)", sample_id),
                ylab = "Count", col = cols, border = NA, ylim = c(0, max(counts, 1) * 1.25))
  text(bp, counts + max(counts, 1) * 0.05, labels = counts, pos = 3)
} else {
  plot.new(); text(0.5, 0.5, "No variants detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig1))

# Figure 2: VAF & Depth Distribution
fig2 <- file.path(fig_dir, sprintf("%s.02_vaf_and_depth_distribution.png", sample_id))
png(fig2, width = 900, height = 500, res = 120)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
if (nrow(combined) > 0) {
  vaf_vals <- as.numeric(combined$VAF[!is.na(combined$VAF) & combined$VAF != "NA"])
  if (length(vaf_vals) > 0) {
    hist(vaf_vals, main = sprintf("VAF Distribution (%s)", sample_id),
         xlab = "Variant Allele Frequency (VAF)", ylab = "Frequency",
         col = "#377eb8", border = "white", breaks = 20, xlim = c(0, 1))
  } else {
    plot.new(); text(0.5, 0.5, "No VAF data")
  }

  dp_vals <- as.numeric(combined$DP[!is.na(combined$DP) & combined$DP != "NA"])
  if (length(dp_vals) > 0) {
    boxplot(dp_vals, main = sprintf("Read Depth (%s)", sample_id),
            ylab = "Depth (DP)", col = "#4daf4a", outline = FALSE)
  } else {
    plot.new(); text(0.5, 0.5, "No Depth data")
  }
} else {
  plot.new(); text(0.5, 0.5, "No variants detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig2))

# Figure 3: Ti/Tv Ratio & Mutation Spectrum
fig3 <- file.path(fig_dir, sprintf("%s.03_titv_and_mutation_spectrum.png", sample_id))
png(fig3, width = 900, height = 500, res = 120)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
if (nrow(stats_df) > 0) {
  titv_val <- stats_df$value[stats_df$metric == "titv_ratio"]
  ti_cnt <- get_stat("transitions_ti")
  tv_cnt <- get_stat("transversions_tv")

  # Panel A: Ti vs Tv counts & ratio
  bp <- barplot(c(Transitions=ti_cnt, Transversions=tv_cnt),
                main = sprintf("Ti/Tv Ratio: %s (%s)", titv_val, sample_id),
                ylab = "Count", col = c("#ff7f00", "#984ea3"), border = NA)
  text(bp, c(ti_cnt, tv_cnt) + max(ti_cnt, tv_cnt, 1)*0.05, labels = c(ti_cnt, tv_cnt), pos = 3)

  # Panel B: 6-class mutation spectrum
  spec_counts <- c(
    "C>T" = get_stat("spectrum_C_T"),
    "T>C" = get_stat("spectrum_T_C"),
    "C>A" = get_stat("spectrum_C_A"),
    "C>G" = get_stat("spectrum_C_G"),
    "T>A" = get_stat("spectrum_T_A"),
    "T>G" = get_stat("spectrum_T_G")
  )
  bp2 <- barplot(spec_counts, main = "Substitution Spectrum",
                 ylab = "Count", col = "#e41a1c", border = NA)
  text(bp2, spec_counts + max(spec_counts, 1)*0.05, labels = spec_counts, pos = 3)
} else {
  plot.new(); text(0.5, 0.5, "No variants detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig3))

# Figure 4: Chromosome Density
fig4 <- file.path(fig_dir, sprintf("%s.04_chromosome_snv_density.png", sample_id))
png(fig4, width = 1000, height = 500, res = 120)
par(mar = c(6, 5, 4, 2))
if (nrow(combined) > 0 && "CHROM" %in% colnames(combined)) {
  chrom_table <- table(combined$CHROM)
  barplot(chrom_table, main = sprintf("Small Variants per Chromosome (%s)", sample_id),
          ylab = "Variant Count", col = "#377eb8", border = NA, las = 2, cex.names = 0.8)
} else {
  plot.new(); text(0.5, 0.5, "No variants detected")
}
dev.off()
cat(sprintf("  Saved: %s\n", fig4))

cat(sprintf("=== [07_generate_report.R] Done. Figures written to: %s ===\n", fig_dir))
