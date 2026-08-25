#!/usr/bin/env Rscript
# =============================================================================
# 06_summary_stats.R (pipelines/cnv)
# Computes copy number variation metrics: variant counts per category (loss/gain),
# mean/median segment lengths, total genomic span impacted, and length distribution
# statistics.
#
# Usage:
#   Rscript scripts/R/06_summary_stats.R [path/to/pipeline_config.yaml]
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

cat(sprintf("=== [06_summary_stats.R] Using config: %s ===\n", config_file))

sample_id <- get_yaml_val(config_file, "sample_id")
output_dir <- get_yaml_val(config_file, "output_dir")

cat(sprintf("Sample ID: %s\n", sample_id))

del_file <- file.path(output_dir, "deletions", sprintf("%s.deletions.tsv", sample_id))
dup_file <- file.path(output_dir, "duplications", sprintf("%s.duplications.tsv", sample_id))
qc_dir <- file.path(output_dir, "qc_summary")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

combined_out <- file.path(qc_dir, sprintf("%s.all_cnvs_combined.tsv", sample_id))
stats_out <- file.path(qc_dir, sprintf("%s.summary_statistics.tsv", sample_id))

dels <- if (file.exists(del_file) && file.info(del_file)$size > 0) read.delim(del_file, stringsAsFactors = FALSE) else data.frame()
dups <- if (file.exists(dup_file) && file.info(dup_file)$size > 0) read.delim(dup_file, stringsAsFactors = FALSE) else data.frame()

if (nrow(dels) > 0) dels$CATEGORY <- "deletions"
if (nrow(dups) > 0) dups$CATEGORY <- "duplications"

combined <- rbind(dels, dups)

if (nrow(combined) > 0) {
  write.table(combined, combined_out, sep = "\t", quote = FALSE, row.names = FALSE)

  # Compute per-category stats
  categories <- unique(combined$CATEGORY)
  stats_list <- list()

  for (cat_name in categories) {
    sub_df <- combined[combined$CATEGORY == cat_name, ]
    lens <- as.numeric(sub_df$SVLEN)
    lens <- lens[!is.na(lens)]

    s_df <- data.frame(
      category = cat_name,
      n_variants = nrow(sub_df),
      mean_cnv_len_bp = if (length(lens) > 0) round(mean(lens), 1) else NA,
      median_cnv_len_bp = if (length(lens) > 0) round(median(lens), 1) else NA,
      min_cnv_len_bp = if (length(lens) > 0) min(lens) else NA,
      max_cnv_len_bp = if (length(lens) > 0) max(lens) else NA,
      total_span_bp = if (length(lens) > 0) sum(lens) else 0,
      stringsAsFactors = FALSE
    )
    stats_list[[cat_name]] <- s_df
  }

  stats_df <- do.call(rbind, stats_list)
  write.table(stats_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)
  print(stats_df)
} else {
  cat("No CNVs found across categories.\n")
  empty_df <- data.frame(category = character(), n_variants = integer(), mean_cnv_len_bp = numeric(),
                         median_cnv_len_bp = numeric(), min_cnv_len_bp = integer(), max_cnv_len_bp = integer(),
                         total_span_bp = numeric(), stringsAsFactors = FALSE)
  write.table(empty_df, combined_out, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(empty_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("Combined table written: %s\n", combined_out))
cat(sprintf("Summary stats written : %s\n", stats_out))
cat("=== [06_summary_stats.R] Done ===\n")
