#!/usr/bin/env Rscript
# =============================================================================
# 06_summary_stats.R (pipelines/snv)
# Computes biological small variant metrics: SNV/Indel counts, Transition/Transversion
# (Ti/Tv) ratio, 6-class mutation spectrum, Indel insertion/deletion ratio, mean depth,
# and median variant allele frequency (VAF).
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

snv_file <- file.path(output_dir, "snvs", sprintf("%s.snvs.tsv", sample_id))
indel_file <- file.path(output_dir, "indels", sprintf("%s.indels.tsv", sample_id))
clinvar_file <- file.path(output_dir, "clinvar", sprintf("%s.clinvar_pathogenic.tsv", sample_id))
qc_dir <- file.path(output_dir, "qc_summary")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

combined_out <- file.path(qc_dir, sprintf("%s.all_snvs_combined.tsv", sample_id))
stats_out <- file.path(qc_dir, sprintf("%s.summary_statistics.tsv", sample_id))

snvs <- if (file.exists(snv_file) && file.info(snv_file)$size > 0) read.delim(snv_file, stringsAsFactors = FALSE) else data.frame()
indels <- if (file.exists(indel_file) && file.info(indel_file)$size > 0) read.delim(indel_file, stringsAsFactors = FALSE) else data.frame()
clinvar <- if (file.exists(clinvar_file) && file.info(clinvar_file)$size > 0) read.delim(clinvar_file, stringsAsFactors = FALSE) else data.frame()

if (nrow(snvs) > 0) snvs$VARIANT_CLASS <- "SNV"
if (nrow(indels) > 0) indels$VARIANT_CLASS <- "INDEL"

combined <- rbind(snvs, indels)

# Ti/Tv Calculation function for SNVs
calc_titv <- function(snv_df) {
  if (nrow(snv_df) == 0) return(list(ti=0, tv=0, ratio=NA))
  ti <- 0
  tv <- 0

  for (i in seq_len(nrow(snv_df))) {
    ref <- toupper(snv_df$REF[i])
    alt <- toupper(snv_df$ALT[i])

    if (nchar(ref) != 1 || nchar(alt) != 1) next

    pair <- paste0(ref, alt)
    if (pair %in% c("AG", "GA", "CT", "TC")) {
      ti <- ti + 1
    } else if (pair %in% c("AC", "CA", "AT", "TA", "CG", "GC", "GT", "TG")) {
      tv <- tv + 1
    }
  }

  ratio <- if (tv > 0) round(ti / tv, 3) else NA
  return(list(ti=ti, tv=tv, ratio=ratio))
}

titv_res <- calc_titv(snvs)

# Compute 6-class mutation spectrum
calc_spectrum <- function(snv_df) {
  classes <- c("C>T", "T>C", "C>A", "C>G", "T>A", "T>G")
  counts <- setNames(rep(0, 6), classes)

  if (nrow(snv_df) > 0) {
    for (i in seq_len(nrow(snv_df))) {
      ref <- toupper(snv_df$REF[i])
      alt <- toupper(snv_df$ALT[i])
      if (nchar(ref) != 1 || nchar(alt) != 1) next

      # Collapse complementary strands to pyrimidines C and T
      if (ref == "G") { ref <- "C"; alt <- switch(alt, "A"="T", "C"="G", "T"="A", alt) }
      else if (ref == "A") { ref <- "T"; alt <- switch(alt, "G"="C", "C"="G", "T"="A", alt) }

      mut <- paste0(ref, ">", alt)
      if (mut %in% classes) {
        counts[mut] <- counts[mut] + 1
      }
    }
  }
  return(counts)
}

spectrum <- calc_spectrum(snvs)

# Indel insertion vs deletion breakdown
n_ins <- 0
n_del <- 0
if (nrow(indels) > 0) {
  for (i in seq_len(nrow(indels))) {
    l_ref <- nchar(indels$REF[i])
    l_alt <- nchar(indels$ALT[i])
    if (l_alt > l_ref) n_ins <- n_ins + 1
    else if (l_ref > l_alt) n_del <- n_del + 1
  }
}

mean_dp <- if (nrow(combined) > 0 && "DP" %in% colnames(combined)) round(mean(as.numeric(combined$DP), na.rm=TRUE), 1) else NA
med_vaf <- if (nrow(combined) > 0 && "VAF" %in% colnames(combined)) round(median(as.numeric(combined$VAF), na.rm=TRUE), 4) else NA

stats_df <- data.frame(
  metric = c("total_snvs", "total_indels", "transitions_ti", "transversions_tv", "titv_ratio",
             "indel_insertions", "indel_deletions", "clinvar_pathogenic_count",
             "mean_read_depth", "median_vaf",
             "spectrum_C_T", "spectrum_T_C", "spectrum_C_A", "spectrum_C_G", "spectrum_T_A", "spectrum_T_G"),
  value = c(nrow(snvs), nrow(indels), titv_res$ti, titv_res$tv, ifelse(is.na(titv_res$ratio), "NA", as.character(titv_res$ratio)),
            n_ins, n_del, nrow(clinvar),
            ifelse(is.na(mean_dp), "NA", as.character(mean_dp)), ifelse(is.na(med_vaf), "NA", as.character(med_vaf)),
            spectrum["C>T"], spectrum["T>C"], spectrum["C>A"], spectrum["C>G"], spectrum["T>A"], spectrum["T>G"]),
  stringsAsFactors = FALSE
)

write.table(combined, combined_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stats_df, stats_out, sep = "\t", quote = FALSE, row.names = FALSE)

print(stats_df)

cat(sprintf("Combined table written: %s\n", combined_out))
cat(sprintf("Summary stats written : %s\n", stats_out))
cat("=== [06_summary_stats.R] Done ===\n")
