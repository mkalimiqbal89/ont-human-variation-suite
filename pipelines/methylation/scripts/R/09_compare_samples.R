#!/usr/bin/env Rscript
# =============================================================================
# 09_compare_samples.R
# Compares two or more samples at gene, promoter and CpG-island level.
#
# Works from the per-feature tables written by stage 05, which already carry
# Total_coverage and Modified_read_calls per feature. That is exactly what a
# pooled-read test needs, and it means this stage never touches a site-level
# file: inputs are tens of thousands of rows, not tens of millions.
#
# BASE R + stats ONLY — no packages. fisher.test, chisq.test and p.adjust are
# all base installs. Same reasoning as stage 07: this must run unattended.
#
# STATISTICAL CAVEAT, stated here because it governs how the output may be used:
# a pooled-read test treats every read as an independent observation. Nanopore
# reads span many CpGs, so reads within a feature are correlated, and with one
# sample per group there is no biological replication. The p-values are
# ANTI-CONSERVATIVE. Use them to rank candidates for follow-up; use the effect
# size (delta) for anything you intend to claim. See docs/COMPARISON_CAVEATS.md.
#
# Usage:
#   Rscript scripts/R/09_compare_samples.R <sample_sheet.tsv> <out_dir> [config.yaml]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
args <- args[!grepl("^#", args)]
if (length(args) < 2) {
  stop("Usage: 09_compare_samples.R <sample_sheet.tsv> <out_dir> [config.yaml]",
       call. = FALSE)
}
sheet_file <- args[1]
out_dir    <- args[2]

script_path <- (function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  NA_character_
})()
repo_dir <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."))
} else normalizePath(".")
suite_dir <- normalizePath(file.path(repo_dir, "..", ".."))

config_file <- if (length(args) >= 3) args[3] else
  file.path(repo_dir, "config", "pipeline_config.yaml")

for (f in c(sheet_file, config_file)) {
  if (!file.exists(f)) stop("Not found: ", f, call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Could not create ", out_dir, call. = FALSE)

cat("=== [09_compare_samples.R] ===\n")
cat("Sample sheet:", sheet_file, "\n")
cat("Output       :", out_dir, "\n")

# --- Config ------------------------------------------------------------------
# Shared with 07_generate_report.R via R/lib/lib_common.R rather than each
# stage defining its own copy of this reader.
lib_common_r <- file.path(suite_dir, "R", "lib", "lib_common.R")
if (!file.exists(lib_common_r)) {
  stop("Cannot find ", lib_common_r, call. = FALSE)
}
source(lib_common_r)
cfg_lines <- readLines(config_file, warn = FALSE)
cfg_get <- cfg_get_factory(cfg_lines)
min_shared  <- as.numeric(cfg_get("min_shared_cpgs", "5"))
max_cpg_ratio <- as.numeric(cfg_get("max_cpg_ratio", "3"))
max_cov_ratio <- as.numeric(cfg_get("max_coverage_ratio", "5"))
min_delta   <- as.numeric(cfg_get("min_delta_percent", "20"))
fdr_alpha   <- as.numeric(cfg_get("fdr_alpha", "0.05"))
stat_test   <- cfg_get("stat_test", "fisher")
top_n       <- as.numeric(cfg_get("top_n_report", "100"))
mod_code    <- cfg_get("primary_mod_code", "m")
fig_w       <- as.numeric(cfg_get("fig_width", "2000"))
fig_h       <- as.numeric(cfg_get("fig_height", "1400"))
fig_res     <- as.numeric(cfg_get("fig_res", "150"))
mod_label   <- switch(mod_code, m = "5mC", h = "5hmC", a = "6mA", mod_code)

# --- Sample sheet ------------------------------------------------------------
sheet <- read.delim(sheet_file, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, comment.char = "#")
req <- c("sample_id", "results_dir")
missing_cols <- setdiff(req, names(sheet))
if (length(missing_cols)) {
  stop("Sample sheet must have columns: ", paste(req, collapse = ", "),
       " (missing: ", paste(missing_cols, collapse = ", "), ")", call. = FALSE)
}
if (!"label" %in% names(sheet)) sheet$label <- sheet$sample_id
if (!"group" %in% names(sheet)) sheet$group <- NA_character_
sheet$label[is.na(sheet$label) | sheet$label == ""] <-
  sheet$sample_id[is.na(sheet$label) | sheet$label == ""]

if (nrow(sheet) < 2) stop("Need at least 2 samples; sheet has ", nrow(sheet),
                          call. = FALSE)
if (anyDuplicated(sheet$sample_id)) {
  stop("Duplicate sample_id in sheet: ",
       paste(unique(sheet$sample_id[duplicated(sheet$sample_id)]), collapse = ", "),
       call. = FALSE)
}
n_s <- nrow(sheet)
cat("Samples      :", n_s, "->", paste(sheet$label, collapse = ", "), "\n")
cat("Test         :", stat_test,
    if (stat_test == "fisher" && n_s > 2) "(falls back to chisq for >2 samples)" else "", "\n")
cat("Thresholds   : min_shared_cpgs=", min_shared, "  min_delta=", min_delta,
    "pp  FDR alpha=", fdr_alpha, "\n\n", sep = "")

# --- Feature classes ---------------------------------------------------------
# subdir, filename suffix, id column, cpg column, coverage col, modified col,
# weighted col, min-cpgs-flag col, and a label column carried into the output.
classes <- list(
  gene = list(dir = "06_gene_summary", suffix = "gene_methylation_summary.tsv",
              id = "Feature_id", name = "Feature_name", type = "Feature_type",
              cpg = "CpG_sites", cov = "Total_coverage",
              mod = "Modified_read_calls",
              w = "Coverage_weighted_methylation_percent", flag = "Min_cpgs_met"),
  promoter = list(dir = "07_promoter_summary", suffix = "promoter_methylation_summary.tsv",
              id = "Feature_id", name = "Feature_name", type = "Feature_type",
              cpg = "CpG_sites", cov = "Total_coverage",
              mod = "Modified_read_calls",
              w = "Coverage_weighted_methylation_percent", flag = "Min_cpgs_met"),
  cpg_island = list(dir = "08_cpg_island_summary", suffix = "cpg_island_methylation_summary.tsv",
              id = "Island_id", name = "Island_name", type = NA_character_,
              cpg = "CpG_sites", cov = "Total_coverage",
              mod = "Modified_read_calls",
              w = "Coverage_weighted_methylation_percent", flag = "Min_cpgs_met")
)

read_tsv <- function(f) read.delim(f, sep = "\t", header = TRUE,
                                  stringsAsFactors = FALSE, quote = "",
                                  check.names = FALSE)

summary_rows <- list()
written <- character(0)

# =============================================================================
# Global comparison
# =============================================================================
cat("--- Global comparison ---\n")
glob <- data.frame()
for (i in seq_len(n_s)) {
  f <- file.path(sheet$results_dir[i], "03_global_methylation",
                 paste0(sheet$sample_id[i], ".global_methylation_summary.tsv"))
  if (!file.exists(f)) {
    cat("  [WARN] missing:", f, "\n"); next
  }
  d <- read_tsv(f)
  getv <- function(k) {
    h <- which(d[[1]] == k)
    if (!length(h)) NA_character_ else d[[2]][h[1]]
  }
  glob <- rbind(glob, data.frame(
    sample_id = sheet$sample_id[i], label = sheet$label[i], group = sheet$group[i],
    retained_sites = getv("Retained_CpG_sites"),
    mean_coverage = getv("Mean_coverage"),
    unweighted_percent = getv("Unweighted_mean_site_methylation_percent"),
    weighted_percent = getv("Coverage_weighted_methylation_percent"),
    stringsAsFactors = FALSE))
}
if (nrow(glob)) {
  f_glob <- file.path(out_dir, "global_comparison.tsv")
  write.table(glob, f_glob, sep = "\t", quote = FALSE, row.names = FALSE)
  written <- c(written, basename(f_glob))
  print(glob[, c("label", "retained_sites", "mean_coverage", "weighted_percent")],
        row.names = FALSE)
  # A large difference in mean coverage between samples is worth knowing before
  # interpreting any per-feature delta.
  mc <- suppressWarnings(as.numeric(glob$mean_coverage))
  if (all(!is.na(mc)) && length(mc) > 1 && max(mc) / min(mc) > 1.5) {
    cat("\n  [WARN] mean coverage differs by more than 1.5x across samples",
        sprintf("(%.1f to %.1f).", min(mc), max(mc)),
        "\n         Low-coverage features are noisier; min_shared_cpgs mitigates",
        "\n         but does not eliminate this.\n")
  }
} else {
  cat("  [WARN] no global summaries found\n")
}

# =============================================================================
# Per-feature comparison, one class at a time
# =============================================================================
compare_class <- function(cls_name, cls) {
  cat("\n--- ", cls_name, " ---\n", sep = "")

  tabs <- list()
  for (i in seq_len(n_s)) {
    f <- file.path(sheet$results_dir[i], cls$dir,
                   paste0(sheet$sample_id[i], ".", cls$suffix))
    if (!file.exists(f)) {
      cat("  [SKIP] ", sheet$label[i], ": missing ", basename(f), "\n", sep = "")
      return(invisible(NULL))
    }
    d <- read_tsv(f)
    need <- c(cls$id, cls$cpg, cls$cov, cls$mod, cls$w, cls$flag)
    if (!all(need %in% names(d))) {
      cat("  [SKIP] ", sheet$label[i], ": columns missing (",
          paste(setdiff(need, names(d)), collapse = ", "), ")\n", sep = "")
      return(invisible(NULL))
    }
    keep <- d[[cls$cpg]] >= min_shared
    d <- d[keep, , drop = FALSE]
    tabs[[sheet$sample_id[i]]] <- d
  }

  # Inner join on feature id: only features present in EVERY sample with enough
  # covered CpGs. Comparing a feature absent from one sample is not a comparison.
  ids <- Reduce(intersect, lapply(tabs, function(d) d[[cls$id]]))
  cat("  features passing min_shared_cpgs=", min_shared, " in all samples: ",
      length(ids), "\n", sep = "")
  if (!length(ids)) {
    cat("  [SKIP] no shared features\n")
    return(invisible(NULL))
  }

  res <- data.frame(feature_id = ids, stringsAsFactors = FALSE)
  # Names/types from the first sample; identical across samples by construction.
  d1 <- tabs[[1]][match(ids, tabs[[1]][[cls$id]]), , drop = FALSE]
  if (!is.na(cls$name) && cls$name %in% names(d1)) res$feature_name <- d1[[cls$name]]
  if (!is.na(cls$type) && cls$type %in% names(d1)) res$feature_type <- d1[[cls$type]]

  covm <- modm <- wm <- cpgm <- matrix(NA_real_, nrow = length(ids), ncol = n_s)
  for (i in seq_len(n_s)) {
    d <- tabs[[i]][match(ids, tabs[[i]][[cls$id]]), , drop = FALSE]
    covm[, i] <- as.numeric(d[[cls$cov]])
    modm[, i] <- as.numeric(d[[cls$mod]])
    wm[, i]   <- as.numeric(d[[cls$w]])
    cpgm[, i] <- as.numeric(d[[cls$cpg]])
  }
  for (i in seq_len(n_s)) {
    res[[paste0("methyl_", sheet$label[i])]] <- round(wm[, i], 4)
    res[[paste0("cpgs_",   sheet$label[i])]] <- cpgm[, i]
    res[[paste0("cov_",    sheet$label[i])]] <- covm[, i]
  }

  if (n_s == 2) {
    res$delta <- round(wm[, 2] - wm[, 1], 4)
    res$abs_delta <- abs(res$delta)
    res$direction <- ifelse(res$delta > 0, "hyper", ifelse(res$delta < 0, "hypo", "same"))
  } else {
    res$range <- round(apply(wm, 1, function(x) max(x) - min(x)), 4)
    res$sd <- round(apply(wm, 1, sd), 4)
    res$abs_delta <- res$range
  }

  # --- Statistics ---
  use_test <- stat_test
  if (use_test == "fisher" && n_s > 2) use_test <- "chisq"
  if (use_test == "none") {
    res$p_value <- NA_real_; res$q_value <- NA_real_
  } else {
    unmod <- covm - modm
    # Guard: modified counts can never exceed coverage. If they do, the input is
    # wrong and a test on it would be meaningless.
    bad <- which(unmod < 0)
    if (length(bad)) {
      cat("  [WARN] ", length(bad), " feature(s) have modified > coverage; ",
          "excluding them from testing\n", sep = "")
      unmod[unmod < 0] <- NA
    }
    pv <- rep(NA_real_, length(ids))
    for (r in seq_along(ids)) {
      m <- rbind(modm[r, ], unmod[r, ])
      if (any(is.na(m)) || any(colSums(m) <= 0)) next
      # A table with no variation has nothing to test.
      if (all(m[1, ] == 0) || all(m[2, ] == 0)) { pv[r] <- 1; next }
      pv[r] <- tryCatch({
        if (use_test == "fisher") fisher.test(m)$p.value
        else suppressWarnings(chisq.test(m)$p.value)
      }, error = function(e) NA_real_)
    }
    res$p_value <- pv
    res$q_value <- p.adjust(pv, method = "BH")
  }

  # --- Coverage balance ---------------------------------------------------
  # min_shared_cpgs requires enough CpGs in EACH sample, but says nothing about
  # whether the counts are COMPARABLE. On the first real two-sample run the top
  # promoter hits looked like this:
  #
  #   ENSG00000290383   10 CpGs / cov 100   vs   173 CpGs / cov 2514
  #   MIR3180-3         15 CpGs / cov 153   vs   181 CpGs / cov 2646
  #
  # A 17-fold difference in covered CpGs for the same promoter is a mappability
  # difference, not a methylation difference: these are pseudogene and
  # subtelomeric loci where read placement is unstable. Comparing 10 CpGs against
  # 173 also means the two "methylation" values summarise different sets of
  # positions within the feature, so the delta is not even measuring the same
  # thing twice.
  #
  # Such features are reported, with their ratios, but excluded from being called
  # differential.
  cpg_ratio <- apply(cpgm, 1, function(x) {
    mn <- min(x); if (mn <= 0) return(Inf); max(x) / mn
  })
  cov_ratio <- apply(covm, 1, function(x) {
    mn <- min(x); if (mn <= 0) return(Inf); max(x) / mn
  })
  res$cpg_ratio <- round(cpg_ratio, 3)
  res$coverage_ratio <- round(cov_ratio, 3)
  balanced <- is.finite(cpg_ratio) & is.finite(cov_ratio) &
              cpg_ratio <= max_cpg_ratio & cov_ratio <= max_cov_ratio
  res$balanced <- ifelse(balanced, "yes", "no")

  # --- Calls ---
  sig <- !is.na(res$q_value) & res$q_value < fdr_alpha
  big <- res$abs_delta >= min_delta
  res$differential <- ifelse(sig & big & balanced, "yes",
                      ifelse(sig & big & !balanced, "imbalanced",
                      ifelse(big, "effect_only",
                      ifelse(sig, "significant_only", "no"))))

  # Balanced features first, then by effect size. Without this the table is
  # topped by mappability artefacts, which is how the problem was found.
  ord <- order(res$balanced != "yes", -res$abs_delta, res$q_value, na.last = TRUE)
  res <- res[ord, , drop = FALSE]

  f_out <- file.path(out_dir, paste0(cls_name, "_comparison.tsv"))
  write.table(res, f_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  written <<- c(written, basename(f_out))

  f_top <- file.path(out_dir, paste0(cls_name, "_top_differential.tsv"))
  top <- res[res$differential == "yes", , drop = FALSE]
  if (!nrow(top)) top <- head(res, min(top_n, nrow(res)))
  else top <- head(top, min(top_n, nrow(top)))
  write.table(top, f_top, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  written <<- c(written, basename(f_top))

  n_diff <- sum(res$differential == "yes")
  n_imb  <- sum(res$differential == "imbalanced")
  n_eff  <- sum(res$differential == "effect_only")
  n_sigo <- sum(res$differential == "significant_only")
  n_unbal <- sum(res$balanced == "no")
  cat("  differential (|delta| >= ", min_delta, "pp, q < ", fdr_alpha,
      ", balanced coverage): ", n_diff, "\n", sep = "")
  cat("  would qualify but coverage is imbalanced: ", n_imb,
      "  <- excluded: CpG ratio > ", max_cpg_ratio, " or coverage ratio > ",
      max_cov_ratio, "\n", sep = "")
  cat("  large effect but not significant: ", n_eff, "\n", sep = "")
  cat("  significant but small effect:     ", n_sigo,
      "  <- with this many reads, tiny differences reach significance;\n",
      "                                        this is why the effect-size ",
      "threshold exists\n", sep = "")
  cat("  coverage-imbalanced features overall: ", n_unbal, " of ", length(ids),
      sprintf(" (%.1f%%)", 100 * n_unbal / length(ids)), "\n", sep = "")

  summary_rows[[cls_name]] <<- data.frame(
    feature_class = cls_name, shared_features = length(ids),
    differential = n_diff, imbalanced_excluded = n_imb,
    effect_only = n_eff, significant_only = n_sigo,
    coverage_imbalanced = n_unbal,
    stringsAsFactors = FALSE)

  # --- Figures ---
  if (n_s == 2) {
    fpng <- file.path(out_dir, paste0(cls_name, "_scatter.png"))
    png(fpng, width = fig_w, height = fig_h, res = fig_res)
    on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
    par(mar = c(5, 5.5, 4, 2))
    col_pt <- ifelse(res$differential == "yes", "#D55E00",
              ifelse(res$balanced == "no", "#E69F0055", "#99999955"))
    plot(wm[ord, 1], wm[ord, 2], pch = 16, cex = 0.4, col = col_pt,
         xlim = c(0, 100), ylim = c(0, 100),
         xlab = paste0(sheet$label[1], "  ", mod_label, " (%)"),
         ylab = paste0(sheet$label[2], "  ", mod_label, " (%)"),
         main = paste0(cls_name, ": ", mod_label, " methylation, ",
                       sheet$label[2], " vs ", sheet$label[1]))
    abline(0, 1, col = "grey40", lty = 2)
    abline(min_delta, 1, col = "#0072B2", lty = 3)
    abline(-min_delta, 1, col = "#0072B2", lty = 3)
    legend("topleft", bty = "n", cex = 0.8,
           legend = c(sprintf("differential (n=%d)", n_diff),
                      sprintf("coverage-imbalanced (n=%d)", n_unbal),
                      sprintf("other (n=%d)", nrow(res) - n_diff - n_unbal),
                      sprintf("dotted: +/- %g pp", min_delta)),
           col = c("#D55E00", "#E69F00", "#999999", "#0072B2"),
           pch = c(16, 16, 16, NA), lty = c(NA, NA, NA, 3))
    r <- suppressWarnings(cor(wm[, 1], wm[, 2], use = "complete.obs"))
    mtext(sprintf("Pearson r = %.4f over %d shared features", r, nrow(res)),
          side = 3, line = 0.2, cex = 0.8, col = "grey35")
    dev.off()
    written <<- c(written, basename(fpng))
  } else {
    # For >2 samples a correlation matrix is the readable summary.
    fpng <- file.path(out_dir, paste0(cls_name, "_correlation.png"))
    png(fpng, width = fig_w, height = fig_h, res = fig_res)
    on.exit(if (dev.cur() > 1) dev.off(), add = TRUE)
    cm <- suppressWarnings(cor(wm, use = "complete.obs"))
    par(mar = c(7, 7, 4, 3))
    image(seq_len(n_s), seq_len(n_s), cm, axes = FALSE, zlim = c(0, 1),
          col = colorRampPalette(c("#FFFFFF", "#0072B2"))(64),
          xlab = "", ylab = "",
          main = paste0(cls_name, ": correlation of ", mod_label,
                        " across samples"))
    axis(1, at = seq_len(n_s), labels = sheet$label, las = 2, cex.axis = 0.8)
    axis(2, at = seq_len(n_s), labels = sheet$label, las = 2, cex.axis = 0.8)
    for (a in seq_len(n_s)) for (b in seq_len(n_s)) {
      text(a, b, sprintf("%.3f", cm[a, b]), cex = 0.75,
           col = if (cm[a, b] > 0.6) "white" else "grey20")
    }
    dev.off()
    written <<- c(written, basename(fpng))
  }
  invisible(NULL)
}

for (nm in names(classes)) compare_class(nm, classes[[nm]])

# =============================================================================
# Summary
# =============================================================================
if (length(summary_rows)) {
  s <- do.call(rbind, summary_rows)
  f_sum <- file.path(out_dir, "comparison_summary.tsv")
  write.table(s, f_sum, sep = "\t", quote = FALSE, row.names = FALSE)
  written <- c(written, basename(f_sum))
  cat("\n--- Summary ---\n")
  print(s, row.names = FALSE)
}

cat("\n")
cat("Files written:", length(written), "\n")
for (w in written) cat("  ", w, "\n", sep = "")

cat("\n")
cat("IMPORTANT: p-values here come from pooling reads, which treats each read as\n")
cat("an independent observation. Nanopore reads span many CpGs, so reads within a\n")
cat("feature are correlated, and with one sample per group there is no biological\n")
cat("replication. The p-values are anti-conservative. Rank candidates with them;\n")
cat("quote the effect size (delta). See docs/COMPARISON_CAVEATS.md.\n")

if (!length(written)) {
  cat("\n[FATAL] Nothing could be compared. Check that stage 05 ran for every sample.\n")
  quit(status = 1)
}
cat("\n=== [09_compare_samples.R] Done. ===\n")
