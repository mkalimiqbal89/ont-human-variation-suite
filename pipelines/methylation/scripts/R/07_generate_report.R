#!/usr/bin/env Rscript
# =============================================================================
# 07_generate_report.R
# Publication figures, drawn from the small TSVs produced by stages 02, 03, 05
# and 06. Never reads the ~28 M-site data file.
#
# BASE R GRAPHICS ONLY — NO PACKAGES.
# Deliberate. This stage has to run unattended on an HPC where the R library
# path may differ from the login shell's, and a pipeline that fails at the last
# stage because ggplot2 is missing has wasted the whole run. Base R graphics are
# always present and perfectly capable of publication output. The config is read
# with a small regex reader rather than the yaml package, for the same reason and
# mirroring yaml_get() in 00_setup_env.sh.
#
# Every figure is optional: if its input is absent the figure is skipped with a
# warning and the others still render. A missing input is a reason to tell the
# user which stage to run, not to abort.
#
# Usage:
#   Rscript scripts/R/07_generate_report.R [path/to/pipeline_config.yaml]
# =============================================================================

suppressWarnings(rm(list = ls()))

args <- commandArgs(trailingOnly = TRUE)
# zsh passes inline '#' comments through as arguments; ignore them.
args <- args[!grepl("^#", args)]

script_path <- (function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  NA_character_
})()
repo_dir <- if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  normalizePath(".")
}

config_file <- if (length(args) >= 1) args[1] else
  file.path(repo_dir, "config", "pipeline_config.yaml")

if (!file.exists(config_file)) {
  stop("Config file not found: ", config_file, call. = FALSE)
}

cat("=== [07_generate_report.R] Using config:", config_file, "===\n")

# --- Minimal config reader ---------------------------------------------------
# Mirrors yaml_get() in 00_setup_env.sh: first line matching "^\s*key:".
cfg_lines <- readLines(config_file, warn = FALSE)
cfg_get <- function(key, default = NULL) {
  pat <- paste0("^[[:space:]]*", key, ":")
  hit <- grep(pat, cfg_lines, value = TRUE)
  if (!length(hit)) {
    if (is.null(default)) {
      stop("Config key not found and no default: ", key, call. = FALSE)
    }
    return(default)
  }
  v <- sub("^[^:]*:[[:space:]]*", "", hit[1])
  v <- sub("[[:space:]]*(#.*)?$", "", v)
  gsub('^"|"$', "", v)
}
cfg_num <- function(key, default) as.numeric(cfg_get(key, as.character(default)))

sample_id  <- cfg_get("sample_id")
output_dir <- cfg_get("output_dir")
mod_code   <- cfg_get("primary_mod_code", "m")
min_cov    <- cfg_get("min_coverage", "10")
unmeth_max <- cfg_num("unmethylated_max_percent", 20)
meth_min   <- cfg_num("methylated_min_percent", 80)
fig_w      <- cfg_num("fig_width", 2000)
fig_h      <- cfg_num("fig_height", 1400)
fig_res    <- cfg_num("fig_res", 150)
top_n      <- cfg_num("top_n_gene_types", 12)

mod_label <- switch(mod_code, m = "5mC", h = "5hmC", a = "6mA",
                    paste0("mod_", mod_code))

fig_dir <- file.path(output_dir, "09_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(fig_dir)) {
  stop("Could not create figure directory: ", fig_dir, call. = FALSE)
}

cat("Sample :", sample_id, "\n")
cat("Mod    :", mod_label, "\n")
cat("Figures:", fig_dir, "\n\n")

# --- Palette -----------------------------------------------------------------
# Colour-blind safe (Okabe-Ito). Figures must also survive greyscale printing,
# so ordering is by lightness where a sequence is implied.
col_meth   <- "#D55E00"   # vermillion  - methylated
col_inter  <- "#E69F00"   # orange      - intermediate
col_unmeth <- "#0072B2"   # blue        - unmethylated
col_neutral<- "#56B4E9"   # sky blue
col_accent <- "#009E73"   # green
col_grey   <- "#999999"

p <- function(...) file.path(...)
written <- character(0)
skipped <- character(0)

# Which figure is in progress, for the error handler below. An R error inside a
# plotting call would otherwise leave an open device (so a truncated PNG on disk)
# and a stack trace that does not say which figure was being drawn.
current_fig <- NA_character_
current_file <- NA_character_

close_devices <- function() {
  try(while (dev.cur() > 1) dev.off(), silent = TRUE)
}

options(error = function() {
  cat("\n[FATAL] 07_generate_report.R failed ",
      if (!is.na(current_fig)) paste0("while drawing: ", current_fig)
      else "before any figure was started", "\n", sep = "")
  close_devices()
  # Remove the partial PNG: a truncated figure that opens to garbage is worse
  # than no figure, especially if it is picked up for a manuscript later.
  if (!is.na(current_file) && file.exists(current_file)) {
    unlink(current_file)
    cat("        Removed partial output: ", basename(current_file), "\n", sep = "")
  }
  quit(status = 1)
})

open_png <- function(name) {
  f <- p(fig_dir, paste0(sample_id, ".", name, ".png"))
  current_fig <<- name
  current_file <<- f
  png(f, width = fig_w, height = fig_h, res = fig_res)
  f
}
finish <- function(f, label) {
  dev.off()
  current_fig <<- NA_character_
  current_file <<- NA_character_
  written <<- c(written, basename(f))
  cat("  [OK]   ", label, " -> ", basename(f), "\n", sep = "")
}
skip <- function(label, reason) {
  skipped <<- c(skipped, label)
  cat("  [SKIP] ", label, " (", reason, ")\n", sep = "")
}
read_tsv <- function(f) read.delim(f, sep = "\t", header = TRUE,
                                   stringsAsFactors = FALSE, quote = "",
                                   check.names = FALSE)

# Paths to the tables produced upstream.
f_dist   <- p(output_dir, "03_global_methylation", paste0(sample_id, ".methylation_distribution.tsv"))
f_global <- p(output_dir, "03_global_methylation", paste0(sample_id, ".global_methylation_summary.tsv"))
f_class  <- p(output_dir, "03_global_methylation", paste0(sample_id, ".feature_class_summary.tsv"))
f_type   <- p(output_dir, "03_global_methylation", paste0(sample_id, ".methylation_by_gene_type.tsv"))
f_chrom  <- p(output_dir, "04_chromosome_summary", paste0(sample_id, ".chromosome_methylation_summary.tsv"))
f_cov    <- p(output_dir, "02_qc", paste0(sample_id, ".02_coverage_histogram.tsv"))
f_filt   <- p(output_dir, "02_qc", paste0(sample_id, ".03_filter_qc.tsv"))
f_gene   <- p(output_dir, "06_gene_summary", paste0(sample_id, ".gene_methylation_summary.tsv"))
f_prom   <- p(output_dir, "07_promoter_summary", paste0(sample_id, ".promoter_methylation_summary.tsv"))
f_cgi    <- p(output_dir, "08_cpg_island_summary", paste0(sample_id, ".cpg_island_methylation_summary.tsv"))

kv <- function(f, key) {
  if (!file.exists(f)) return(NA_real_)
  d <- read_tsv(f)
  hit <- d[[1]] == key
  if (!any(hit)) return(NA_real_)
  suppressWarnings(as.numeric(d[[2]][which(hit)[1]]))
}
genome_weighted <- kv(f_global, "Coverage_weighted_methylation_percent")

cat("--- Figures ---\n")

# =============================================================================
# 1. Methylation distribution
# =============================================================================
if (!file.exists(f_dist)) {
  skip("methylation distribution", "run stage 06")
} else {
  d <- read_tsv(f_dist)
  f <- open_png("fig1_methylation_distribution")
  par(mar = c(5, 5.5, 4, 2))
  cols <- ifelse(seq_len(nrow(d)) <= 2, col_unmeth,
          ifelse(seq_len(nrow(d)) >= 9, col_meth, col_inter))
  bp <- barplot(d$Sites / 1e6, names.arg = d$Bin, col = cols, border = NA,
                xlab = paste0(mod_label, " methylation (%)"),
                ylab = "CpG sites (millions)",
                main = paste0(sample_id, ": distribution of per-site ", mod_label,
                              " methylation"),
                las = 1, cex.names = 0.9)
  pct <- sprintf("%.1f%%", d$Percent_of_sites)
  text(bp, d$Sites / 1e6, labels = pct, pos = 3, cex = 0.75, xpd = TRUE,
       col = "grey20")
  mtext(paste0("coverage >= ", min_cov, "x; bimodality is the expected shape ",
               "for a mammalian methylome"), side = 3, line = 0.2, cex = 0.8,
        col = "grey35")
  finish(f, "methylation distribution")
}

# =============================================================================
# 2. Per-chromosome methylation
# =============================================================================
if (!file.exists(f_chrom)) {
  skip("per-chromosome methylation", "run stage 06")
} else {
  d <- read_tsv(f_chrom)
  f <- open_png("fig2_chromosome_methylation")
  par(mar = c(6, 5.5, 4, 2))
  bp <- barplot(d$Coverage_weighted_methylation_percent, names.arg = d$Chromosome,
                col = col_neutral, border = NA, las = 2,
                ylab = paste0("Coverage-weighted ", mod_label, " (%)"),
                main = paste0(sample_id, ": ", mod_label, " methylation by chromosome"),
                ylim = c(0, 100), cex.names = 0.85)
  if (!is.na(genome_weighted)) {
    abline(h = genome_weighted, lty = 2, col = col_meth, lwd = 2)
    text(x = par("usr")[2], y = genome_weighted, pos = 2,
         labels = sprintf("genome %.2f%%", genome_weighted),
         col = col_meth, cex = 0.8, offset = 0.4)
  }
  # Sex chromosomes are worth flagging: chrX in a male sample is haploid, so its
  # coverage is halved relative to the autosomes and it interacts with the
  # coverage threshold differently.
  sex <- which(d$Chromosome %in% c("chrX", "chrY", "X", "Y"))
  if (length(sex)) {
    mtext("chrX/chrY differ in ploidy; interpret against known sample sex",
          side = 1, line = 4.6, cex = 0.75, col = "grey35")
  }
  finish(f, "per-chromosome methylation")
}

# =============================================================================
# 3. Feature-class distributions
# =============================================================================
have_feat <- file.exists(f_gene) && file.exists(f_prom) && file.exists(f_cgi)
if (!have_feat) {
  skip("feature-class distributions", "run stage 05")
} else {
  g <- read_tsv(f_gene); pr <- read_tsv(f_prom); ci <- read_tsv(f_cgi)
  # Only features meeting the CpG minimum: a locus backed by one or two CpGs
  # reads 0% or 100% by chance and would widen every box with noise.
  gv  <- g$Coverage_weighted_methylation_percent[g$Min_cpgs_met == "yes"]
  gvp <- g$Coverage_weighted_methylation_percent[g$Min_cpgs_met == "yes" &
                                                 g$Feature_type == "protein_coding"]
  pv  <- pr$Coverage_weighted_methylation_percent[pr$Min_cpgs_met == "yes"]
  pvp <- pr$Coverage_weighted_methylation_percent[pr$Min_cpgs_met == "yes" &
                                                  pr$Feature_type == "protein_coding"]
  cv  <- ci$Coverage_weighted_methylation_percent[ci$Min_cpgs_met == "yes"]

  dat <- list("Gene bodies\n(all)" = gv,
              "Gene bodies\n(protein-coding)" = gvp,
              "Promoters\n(all)" = pv,
              "Promoters\n(protein-coding)" = pvp,
              "CpG islands" = cv)
  # Colours are subset by the SAME mask as the data. Taking the first n colours
  # instead would silently mis-colour the surviving groups whenever a middle
  # group is empty — e.g. CpG islands drawn in the "methylated" colour.
  dat_cols_all <- c(col_meth, col_meth, col_inter, col_inter, col_unmeth)
  keep_i <- vapply(dat, length, integer(1)) > 0
  dat <- dat[keep_i]
  dat_cols <- dat_cols_all[keep_i]
}

# Dropping empty groups is not enough: if NO feature meets the CpG minimum then
# every group is empty, and boxplot() on an empty list throws rather than
# drawing nothing. That is the normal case for a small targeted panel, so it has
# to be a skip, not a crash.
if (!have_feat) {
  # already skipped above
} else if (length(dat) == 0) {
  skip("feature-class distributions",
       paste0("no feature meets min_cpgs_per_feature=",
              cfg_get("min_cpgs_per_feature", "5")))
} else {
  f <- open_png("fig3_feature_class_methylation")
  par(mar = c(6.5, 5.5, 4, 2))
  boxplot(dat, col = dat_cols,
          border = "grey25", outline = FALSE, notch = FALSE, las = 1,
          ylab = paste0("Coverage-weighted ", mod_label, " (%)"),
          main = paste0(sample_id, ": ", mod_label, " methylation by feature class"),
          ylim = c(0, 100), cex.axis = 0.85)
  if (!is.na(genome_weighted)) {
    abline(h = genome_weighted, lty = 2, col = col_grey, lwd = 2)
    text(x = 0.6, y = genome_weighted + 3,
         labels = sprintf("genome-wide %.2f%%", genome_weighted),
         col = "grey30", cex = 0.75, adj = 0)
  }
  ns <- vapply(dat, length, integer(1))
  mtext(paste0("n = ", paste(format(ns, big.mark = ","), collapse = ", "),
               "   (outliers hidden; features with < ",
               cfg_get("min_cpgs_per_feature", "5"), " CpGs excluded)"),
        side = 1, line = 5, cex = 0.72, col = "grey35")
  finish(f, "feature-class distributions")
}

# =============================================================================
# 4. Coverage distribution and the chosen threshold
# =============================================================================
if (!file.exists(f_cov)) {
  skip("coverage distribution", "run stage 02")
} else {
  d <- read_tsv(f_cov)
  d$cov <- suppressWarnings(as.numeric(sub("\\+$", "", d$Coverage)))
  d <- d[!is.na(d$cov), ]
  keep <- d$cov <= 100
  f <- open_png("fig4_coverage_distribution")
  par(mar = c(5, 5.5, 4, 2))
  plot(d$cov[keep], d$Sites[keep] / 1e6, type = "h", lwd = 2, col = col_neutral,
       xlab = "Valid coverage (reads)", ylab = "CpG sites (millions)",
       main = paste0(sample_id, ": coverage distribution of ", mod_label, " sites"),
       las = 1)
  abline(v = as.numeric(min_cov), lty = 2, col = col_meth, lwd = 2)
  retained <- sum(d$Sites[d$cov >= as.numeric(min_cov)])
  total <- sum(d$Sites)
  text(as.numeric(min_cov), par("usr")[4] * 0.9,
       labels = sprintf(" threshold = %sx\n retains %.2f%% of sites",
                        min_cov, retained / total * 100),
       col = col_meth, cex = 0.8, adj = 0)
  mtext("truncated at 100x for display; the full histogram is in the QC table",
        side = 3, line = 0.2, cex = 0.8, col = "grey35")
  finish(f, "coverage distribution")
}

# =============================================================================
# 5. Methylation states
# =============================================================================
if (!file.exists(f_filt)) {
  skip("methylation states", "run stage 03")
} else {
  n_un <- kv(f_filt, "State_unmethylated_sites")
  n_in <- kv(f_filt, "State_intermediate_sites")
  n_me <- kv(f_filt, "State_methylated_sites")
  if (any(is.na(c(n_un, n_in, n_me)))) {
    skip("methylation states", "state counts absent from stage 03 QC")
  } else {
    vals <- c(n_un, n_in, n_me) / 1e6
    labs <- c(sprintf("Unmethylated\n(< %g%%)", unmeth_max),
              sprintf("Intermediate\n(%g-%g%%)", unmeth_max, meth_min),
              sprintf("Methylated\n(>= %g%%)", meth_min))
    f <- open_png("fig5_methylation_states")
    par(mar = c(5.5, 5.5, 4, 2))
    bp <- barplot(vals, names.arg = labs, border = NA,
                  col = c(col_unmeth, col_inter, col_meth),
                  ylab = "CpG sites (millions)",
                  main = paste0(sample_id, ": CpG sites by methylation state"),
                  las = 1, cex.names = 0.9)
    tot <- sum(vals)
    text(bp, vals, labels = sprintf("%.1f%%", vals / tot * 100), pos = 3,
         cex = 0.85, xpd = TRUE, col = "grey20")
    mtext("intermediate sites at adequate coverage are candidates for allele-specific methylation",
          side = 1, line = 4.2, cex = 0.72, col = "grey35")
    finish(f, "methylation states")
  }
}

# =============================================================================
# 6. Promoter methylation by gene_type
# =============================================================================
if (!file.exists(f_type)) {
  skip("promoter methylation by gene_type", "run stage 06")
} else {
  d <- read_tsv(f_type)
  d <- d[d$Feature_class == "promoter", ]
  if (!nrow(d)) {
    skip("promoter methylation by gene_type", "no promoter rows")
  } else {
    d <- d[order(-d$Features), ]
    d <- head(d, top_n)
    d <- d[order(d$Median_weighted_methylation_percent), ]
    f <- open_png("fig6_promoter_by_gene_type")
    par(mar = c(5, 16, 4, 2))
    bp <- barplot(d$Median_weighted_methylation_percent, horiz = TRUE,
                  names.arg = paste0(d$Gene_type, "  (n=",
                                     format(d$Features, big.mark = ","), ")"),
                  col = col_accent, border = NA, las = 1, xlim = c(0, 100),
                  xlab = paste0("Median coverage-weighted ", mod_label, " (%)"),
                  main = paste0(sample_id, ": promoter ", mod_label,
                                " methylation by gene type"),
                  cex.names = 0.8)
    # Interquartile range, because the medians alone hide how bimodal promoters
    # are. Only where Q3 > Q1: arrows() warns "zero-length arrow is of
    # indeterminate angle" otherwise, which happens whenever a gene_type has a
    # single feature.
    iqr <- which(d$Q3 > d$Q1)
    if (length(iqr)) {
      arrows(d$Q1[iqr], bp[iqr], d$Q3[iqr], bp[iqr], angle = 90, code = 3,
             length = 0.03, col = "grey30", lwd = 1.5, xpd = TRUE)
    }
    if (!is.na(genome_weighted)) {
      abline(v = genome_weighted, lty = 2, col = col_meth, lwd = 2)
    }
    mtext("bars are medians, whiskers the interquartile range across promoters",
          side = 1, line = 3.4, cex = 0.75, col = "grey35")
    finish(f, "promoter methylation by gene_type")
  }
}

# =============================================================================
# Manifest
# =============================================================================
man <- p(fig_dir, paste0(sample_id, ".figure_manifest.tsv"))
writeLines(c("Figure\tSource_table",
  paste0(sample_id, ".fig1_methylation_distribution.png\t", basename(f_dist)),
  paste0(sample_id, ".fig2_chromosome_methylation.png\t",   basename(f_chrom)),
  paste0(sample_id, ".fig3_feature_class_methylation.png\t", paste(basename(f_gene), basename(f_prom), basename(f_cgi), sep = ",")),
  paste0(sample_id, ".fig4_coverage_distribution.png\t",    basename(f_cov)),
  paste0(sample_id, ".fig5_methylation_states.png\t",       basename(f_filt)),
  paste0(sample_id, ".fig6_promoter_by_gene_type.png\t",    basename(f_type))), man)

cat("\n")
cat("Figures written:", length(written), "\n")
if (length(skipped)) {
  cat("Figures skipped:", length(skipped), "->", paste(skipped, collapse = "; "), "\n")
}
cat("Manifest:", man, "\n")

if (length(written) == 0) {
  cat("\n[FATAL] No figures could be produced. Run stages 02-06 first.\n")
  quit(status = 1)
}

cat("\n=== [07_generate_report.R] Done. ===\n")
