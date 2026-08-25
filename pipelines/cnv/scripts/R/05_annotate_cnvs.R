#!/usr/bin/env Rscript
# =============================================================================
# 05_annotate_cnvs.R (pipelines/cnv)
# Annotates filtered CNV calls with GENCODE gene bed annotations using bedtools.
# Quantifies gene overlaps, affected transcripts, and generates gene-level
# copy number alteration summaries.
#
# Usage:
#   Rscript scripts/R/05_annotate_cnvs.R [path/to/pipeline_config.yaml]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
cfg_arg <- if (length(args) > 0) args[1] else ""

# Helper to execute shell commands and parse output
run_cmd <- function(cmd) {
  res <- system(cmd, intern = TRUE)
  return(res)
}

# Read YAML helper
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

cat(sprintf("=== [05_annotate_cnvs.R] Using config: %s ===\n", config_file))

sample_id <- get_yaml_val(config_file, "sample_id")
output_dir <- get_yaml_val(config_file, "output_dir")
ref_config_rel <- get_yaml_val(config_file, "config_file")

ref_config <- if (startsWith(ref_config_rel, "/")) ref_config_rel else file.path(dirname(config_file), ref_config_rel)
gene_bed_raw <- if (file.exists(ref_config)) get_yaml_val(ref_config, "gene_bed") else ""
gene_bed <- if (nchar(gene_bed_raw) > 0 && !startsWith(gene_bed_raw, "/")) file.path(dirname(ref_config), gene_bed_raw) else gene_bed_raw

cat(sprintf("Sample ID: %s\n", sample_id))


del_file <- file.path(output_dir, "deletions", sprintf("%s.deletions.tsv", sample_id))
dup_file <- file.path(output_dir, "duplications", sprintf("%s.duplications.tsv", sample_id))
gene_dir <- file.path(output_dir, "gene_cnvs")
dir.create(gene_dir, recursive = TRUE, showWarnings = FALSE)

annotated_out <- file.path(gene_dir, sprintf("%s.annotated_cnvs.tsv", sample_id))
gene_summary_out <- file.path(gene_dir, sprintf("%s.gene_cnv_summary.tsv", sample_id))

# Read deletion & duplication tables
dels <- if (file.exists(del_file) && file.info(del_file)$size > 0) read.delim(del_file, stringsAsFactors = FALSE) else data.frame()
dups <- if (file.exists(dup_file) && file.info(dup_file)$size > 0) read.delim(dup_file, stringsAsFactors = FALSE) else data.frame()

combined <- rbind(dels, dups)

if (nrow(combined) == 0) {
  cat("No CNV calls to annotate.\n")
  empty_df <- data.frame(CHROM=character(), POS=integer(), END=integer(), SVLEN=integer(),
                         SVTYPE=character(), CN=character(), GENES=character(), OVERLAP_FRACTION=character(),
                         stringsAsFactors=FALSE)
  write.table(empty_df, annotated_out, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(empty_df, gene_summary_out, sep="\t", quote=FALSE, row.names=FALSE)
  cat("=== [05_annotate_cnvs.R] Done ===\n")
  quit(status = 0)
}

# Perform bedtools intersect if gene_bed exists
if (file.exists(gene_bed) && file.info(gene_bed)$size > 0) {
  cat(sprintf("Annotating %d CNVs using gene BED: %s\n", nrow(combined), gene_bed))
  tmp_cnv <- tempfile(fileext = ".bed")
  tmp_bed_out <- tempfile(fileext = ".tsv")

  cnv_bed <- data.frame(
    chrom = combined$CHROM,
    start = combined$POS - 1,
    end = combined$END,
    id = paste(combined$CHROM, combined$POS, combined$END, combined$SVTYPE, combined$CN, sep = "_"),
    stringsAsFactors = FALSE
  )
  write.table(cnv_bed, tmp_cnv, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  cmd <- sprintf("bedtools intersect -a '%s' -b '%s' -wa -wb > '%s' 2>/dev/null", tmp_cnv, gene_bed, tmp_bed_out)
  system(cmd)

  gene_map <- list()
  if (file.exists(tmp_bed_out) && file.info(tmp_bed_out)$size > 0) {
    intersect_data <- read.delim(tmp_bed_out, header = FALSE, stringsAsFactors = FALSE)
    # Assume gene BED has gene symbol in column 7 or 4
    for (i in seq_len(nrow(intersect_data))) {
      id <- intersect_data[i, 4]
      g_col <- if (ncol(intersect_data) >= 8) 8 else if (ncol(intersect_data) >= 7) 7 else 4
      g_sym <- intersect_data[i, g_col]
      if (!is.na(g_sym) && g_sym != ".") {
        gene_map[[id]] <- unique(c(gene_map[[id]], g_sym))
      }
    }
  }

  genes_col <- sapply(cnv_bed$id, function(id) {
    if (!is.null(gene_map[[id]])) paste(gene_map[[id]], collapse = ",") else "NONE"
  })

  combined$GENES <- genes_col
  unlink(c(tmp_cnv, tmp_bed_out))
} else {
  cat("GENCODE gene BED not found; skipping positional gene overlap parsing.\n")
  combined$GENES <- "UNANNOTATED"
}

write.table(combined, annotated_out, sep = "\t", quote = FALSE, row.names = FALSE)

# Generate Gene Summary Table
gene_list <- unlist(strsplit(combined$GENES, ","))
gene_list <- gene_list[gene_list != "NONE" & gene_list != "UNANNOTATED"]

if (length(gene_list) > 0) {
  gene_counts <- as.data.frame(table(GENE = gene_list), stringsAsFactors = FALSE)
  colnames(gene_counts) <- c("GENE", "CNV_COUNT")
  write.table(gene_counts, gene_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  empty_summary <- data.frame(GENE = character(), CNV_COUNT = integer(), stringsAsFactors = FALSE)
  write.table(empty_summary, gene_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
}

cat(sprintf("Annotated table saved: %s (%d records)\n", annotated_out, nrow(combined)))
cat(sprintf("Gene summary saved    : %s\n", gene_summary_out))
cat("=== [05_annotate_cnvs.R] Done ===\n")
