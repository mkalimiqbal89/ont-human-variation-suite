#!/usr/bin/env Rscript
# =============================================================================
# 05_annotate_snvs.R (pipelines/snv)
# Annotates filtered SNV and Indel calls with GENCODE gene bed annotations using
# bedtools intersect. Computes per-gene small variant counts and impact burden.
#
# Usage:
#   Rscript scripts/R/05_annotate_snvs.R [path/to/pipeline_config.yaml]
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

cat(sprintf("=== [05_annotate_snvs.R] Using config: %s ===\n", config_file))

sample_id <- get_yaml_val(config_file, "sample_id")
output_dir <- get_yaml_val(config_file, "output_dir")
ref_config_rel <- get_yaml_val(config_file, "config_file")

ref_config <- if (startsWith(ref_config_rel, "/")) ref_config_rel else file.path(dirname(config_file), ref_config_rel)
gene_bed_raw <- if (file.exists(ref_config)) get_yaml_val(ref_config, "gene_bed") else ""
gene_bed <- if (nchar(gene_bed_raw) > 0 && !startsWith(gene_bed_raw, "/")) file.path(dirname(ref_config), gene_bed_raw) else gene_bed_raw

cat(sprintf("Sample ID: %s\n", sample_id))

snv_file <- file.path(output_dir, "snvs", sprintf("%s.snvs.tsv", sample_id))
indel_file <- file.path(output_dir, "indels", sprintf("%s.indels.tsv", sample_id))
gene_dir <- file.path(output_dir, "gene_snvs")
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

annotated_out <- file.path(gene_dir, sprintf("%s.annotated_snvs.tsv", sample_id))
gene_summary_out <- file.path(gene_dir, sprintf("%s.gene_snv_summary.tsv", sample_id))

snvs <- if (file.exists(snv_file) && file.info(snv_file)$size > 0) read.delim(snv_file, stringsAsFactors = FALSE) else data.frame()
indels <- if (file.exists(indel_file) && file.info(indel_file)$size > 0) read.delim(indel_file, stringsAsFactors = FALSE) else data.frame()

if (nrow(snvs) > 0) snvs$CLASS <- "SNV"
if (nrow(indels) > 0) indels$CLASS <- "INDEL"

combined <- rbind(snvs, indels)

if (nrow(combined) == 0) {
  cat("No small variants to annotate.\n")
  empty_df <- data.frame(CHROM=character(), POS=integer(), ID=character(), REF=character(), ALT=character(),
                         GENES=character(), CLASS=character(), stringsAsFactors=FALSE)
  write.table(empty_df, annotated_out, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(empty_df, gene_summary_out, sep="\t", quote=FALSE, row.names=FALSE)
  cat("=== [05_annotate_snvs.R] Done ===\n")
  quit(status = 0)
}

# Perform bedtools intersect if gene_bed exists
if (file.exists(gene_bed) && file.info(gene_bed)$size > 0) {
  cat(sprintf("Annotating %d small variants using gene BED: %s\n", nrow(combined), gene_bed))
  tmp_snv <- tempfile(fileext = ".bed")
  tmp_bed_out <- tempfile(fileext = ".tsv")

  snv_bed <- data.frame(
    chrom = combined$CHROM,
    start = combined$POS - 1,
    end = combined$POS,
    id = paste(combined$CHROM, combined$POS, combined$REF, combined$ALT, seq_len(nrow(combined)), sep = "_"),
    stringsAsFactors = FALSE
  )
  write.table(snv_bed, tmp_snv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cmd <- sprintf("bedtools intersect -a '%s' -b '%s' -wa -wb > '%s' 2>/dev/null", tmp_snv, gene_bed, tmp_bed_out)
  system(cmd)

  gene_map <- list()
  if (file.exists(tmp_bed_out) && file.info(tmp_bed_out)$size > 0) {
    intersect_data <- read.delim(tmp_bed_out, header = FALSE, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(intersect_data))) {
      id <- intersect_data[i, 4]
      g_col <- if (ncol(intersect_data) >= 8) 8 else if (ncol(intersect_data) >= 7) 7 else 4
      g_sym <- intersect_data[i, g_col]
      if (!is.na(g_sym) && g_sym != ".") {
        gene_map[[id]] <- unique(c(gene_map[[id]], g_sym))
      }
    }
  }

  genes_col <- sapply(snv_bed$id, function(id) {
    if (!is.null(gene_map[[id]])) paste(gene_map[[id]], collapse = ",") else "INTERGENIC"
  })

  combined$GENES <- genes_col
  unlink(c(tmp_snv, tmp_bed_out))
} else {
  cat("GENCODE gene BED not found; skipping positional gene overlap parsing.\n")
  combined$GENES <- "UNANNOTATED"
}

write.table(combined, annotated_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Generate Gene Summary Table
gene_list <- unlist(strsplit(combined$GENES, ","))
gene_list <- gene_list[gene_list != "INTERGENIC" & gene_list != "UNANNOTATED"]

if (length(gene_list) > 0) {
  gene_counts <- as.data.frame(table(GENE = gene_list), stringsAsFactors = FALSE)
  colnames(gene_counts) <- c("GENE", "VARIANT_COUNT")
  gene_counts <- gene_counts[order(-gene_counts$VARIANT_COUNT), ]
  write.table(gene_counts, gene_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  empty_summary <- data.frame(GENE = character(), VARIANT_COUNT = integer(), stringsAsFactors = FALSE)
  write.table(empty_summary, gene_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("Annotated table saved: %s (%d records)\n", annotated_out, nrow(combined)))
cat(sprintf("Gene summary saved    : %s\n", gene_summary_out))
cat("=== [05_annotate_snvs.R] Done ===\n")
