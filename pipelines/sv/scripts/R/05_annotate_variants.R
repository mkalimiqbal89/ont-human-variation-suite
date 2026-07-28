#!/usr/bin/env Rscript
# =============================================================================
# 05_annotate_variants.R
#
# Resolves gene fusions/disruptions from BND-type structural variants.
#
# Rationale: Epi2ME's wf-human-variation annotates the SV VCF with SnpEff
# (config: annotation=True), which for BND records already classifies the
# breakpoint effect via the ANN INFO field's Annotation subfield -- values
# like "bidirectional_gene_fusion", "feature_fusion", and "transcript_ablation"
# are SnpEff's own calls, not something this pipeline infers from coordinates.
# This script trusts those calls rather than re-deriving fusions by pairing
# BND mate coordinates from scratch, and extracts the gene pair from either
# the Gene_Name field (ampersand-joined, e.g. "GENE1&GENE2") or the Feature_ID
# field for intergenic fusions (hyphen-joined, e.g. "GENE1-GENE2") -- SnpEff
# uses different fields depending on annotation type, and a record's multiple
# comma-separated annotation blocks must all be scanned since the useful gene
# pair is not always in the first block (see docs/METHODS.md for a real
# example of this from validation).
#
# The flattened TSV from 02_vcf_to_tsv.sh deliberately discarded the full ANN
# field (it can be many KB per record for long insertions) so this script
# re-queries the raw SV VCF via bcftools -- but ONLY for BND records that
# already passed QC in 03_filter_sv_categories.sh, keeping this fast and
# keeping a single point of truth for filtering thresholds.
#
# Usage:
#   Rscript scripts/R/05_annotate_variants.R [path/to/pipeline_config.yaml]
#
# Input:
#   results/translocations/<sample_id>.translocations.tsv  (from bash step 3)
#   raw SV VCF (re-queried for ANN text on BND records only)
#
# Output:
#   results/gene_fusions/<sample_id>.gene_fusions.tsv          (confirmed fusions)
#   results/gene_fusions/<sample_id>.gene_disruptions.tsv      (single-gene breakpoint hits, not fusions)
#   results/gene_fusions/<sample_id>.fusion_annotation_full.tsv (all BND records, every classification, for audit)
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[1] else file.path("config", "pipeline_config.yaml")

if (!file.exists(config_path)) stop("[FATAL] Config file not found: ", config_path)
cat("=== [05_annotate_variants.R] Using config:", config_path, "===\n")

cfg <- yaml.load_file(config_path)

# --- Resolve ${sample.X} placeholders the same way 00_setup_env.sh does ----
resolve_placeholders <- function(s, cfg) {
  if (is.null(s) || !is.character(s)) return(s)
  s <- gsub("\\$\\{sample\\.sample_id\\}", cfg$sample$sample_id, s, fixed = FALSE)
  s <- gsub("\\$\\{sample\\.raw_sample_prefix\\}", cfg$sample$raw_sample_prefix, s, fixed = FALSE)
  s
}

sample_id  <- cfg$sample$sample_id
input_dir  <- cfg$paths$input_dir
output_dir <- cfg$paths$output_dir
sv_vcf     <- file.path(input_dir, resolve_placeholders(cfg$input_files$sv_vcf, cfg))

cat("Sample ID:", sample_id, "\n")
cat("SV VCF   :", sv_vcf, "\n")

if (!file.exists(sv_vcf)) stop("[FATAL] SV VCF not found: ", sv_vcf)

translocations_path <- file.path(output_dir, "translocations", paste0(sample_id, ".translocations.tsv"))
if (!file.exists(translocations_path)) {
  stop("[FATAL] ", translocations_path, " not found. Run scripts/bash/03_filter_sv_categories.sh first.")
}
translocations <- read.delim(translocations_path, stringsAsFactors = FALSE)
cat("Loaded", nrow(translocations), "QC-passing BND/translocation records.\n")

if (nrow(translocations) == 0) {
  cat("No translocation records to annotate. Writing empty output files.\n")
  empty <- translocations[0, ]
  out_dir <- file.path(output_dir, "gene_fusions")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write.table(empty, file.path(out_dir, paste0(sample_id, ".gene_fusions.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  quit(save = "no", status = 0)
}

# --- Pull ANN text for exactly these IDs via bcftools -----------------------
cat("Querying raw ANN annotations for BND records via bcftools...\n")
ids_file <- tempfile()
writeLines(translocations$ID, ids_file)

cmd <- sprintf(
  "bcftools view -i 'INFO/SVTYPE=\"BND\"' %s | bcftools query -f '%%ID\\t%%INFO/ANN\\n'",
  shQuote(sv_vcf)
)
ann_lines <- system(cmd, intern = TRUE)
ann_split <- strsplit(ann_lines, "\t", fixed = TRUE)
ann_df <- data.frame(
  ID  = vapply(ann_split, function(x) x[1], character(1)),
  ANN = vapply(ann_split, function(x) if (length(x) >= 2) x[2] else ".", character(1)),
  stringsAsFactors = FALSE
)
cat("  Retrieved ANN for", nrow(ann_df), "BND records from VCF.\n")

# --- Fusion/disruption classification from ANN, validated against real data -
# (blocks scanned exhaustively; ablation and fusion both extract gene names;
#  fusion gene pairs come from Gene_Name '&' or Feature_ID '-' depending on
#  which SnpEff populated -- see header comment above.)
parse_ann_for_fusion <- function(ann_raw) {
  if (is.na(ann_raw) || ann_raw == "." || ann_raw == "") {
    return(list(fusion_type = "unannotated", geneA = NA_character_, geneB = NA_character_))
  }
  blocks <- strsplit(ann_raw, ",")[[1]]
  any_fusion <- FALSE
  any_ablation <- FALSE
  geneA <- NA_character_; geneB <- NA_character_

  for (blk in blocks) {
    f <- strsplit(blk, "\\|")[[1]]
    if (length(f) < 7) next
    annotation   <- f[2]
    gene_name    <- f[4]
    feature_type <- f[6]
    feature_id   <- f[7]

    is_fusion_block   <- grepl("fusion", annotation, ignore.case = TRUE)
    is_ablation_block <- grepl("ablation", annotation, ignore.case = TRUE)
    if (is_fusion_block) any_fusion <- TRUE
    if (is_ablation_block) any_ablation <- TRUE

    if ((is_fusion_block || is_ablation_block) && (is.na(geneA) || is.na(geneB))) {
      if (!is.na(gene_name) && gene_name != "" && grepl("&", gene_name)) {
        parts <- strsplit(gene_name, "&")[[1]]
        geneA <- parts[1]; geneB <- parts[2]
      } else if (!is.na(feature_type) && feature_type == "intergenic_region" &&
                 !is.na(feature_id) && feature_id != "" && grepl("-", feature_id)) {
        parts <- strsplit(feature_id, "-")[[1]]
        geneA <- parts[1]; geneB <- parts[2]
      } else if (!is.na(gene_name) && gene_name != "" && is.na(geneA)) {
        geneA <- gene_name
      }
    }
  }

  fusion_type <- if (any_fusion) "gene_fusion" else if (any_ablation) "gene_disruption" else "other"
  list(fusion_type = fusion_type, geneA = geneA, geneB = geneB)
}

cat("Classifying", nrow(ann_df), "BND records by breakpoint effect...\n")
classifications <- lapply(ann_df$ANN, parse_ann_for_fusion)
ann_df$fusion_type <- vapply(classifications, function(x) x$fusion_type, character(1))
ann_df$geneA       <- vapply(classifications, function(x) x$geneA, character(1))
ann_df$geneB       <- vapply(classifications, function(x) x$geneB, character(1))

# --- Join back to the QC-passing translocation table (adds CHROM/POS/etc.) --
merged <- merge(translocations, ann_df[, c("ID", "fusion_type", "geneA", "geneB")], by = "ID", all.x = TRUE)
merged$fusion_type[is.na(merged$fusion_type)] <- "unannotated"

# --- Flag likely duplicate/reciprocal calls of the same fusion event -------
# Sniffles2 can emit near-duplicate BND records for the same underlying event
# (e.g. two breakpoints a few kb apart implicating the same gene pair). This
# does not collapse them -- that is a judgment call for manual review -- it
# only flags them so they are not silently miscounted as independent fusions.
merged$gene_pair_key <- ifelse(
  !is.na(merged$geneA) & !is.na(merged$geneB),
  paste(pmin(merged$geneA, merged$geneB), pmax(merged$geneA, merged$geneB), sep = "__"),
  NA_character_
)
dup_counts <- table(merged$gene_pair_key[!is.na(merged$gene_pair_key)])
merged$possible_duplicate_event <- !is.na(merged$gene_pair_key) &
  merged$gene_pair_key %in% names(dup_counts[dup_counts > 1])

# --- Write outputs -----------------------------------------------------------
out_dir <- file.path(output_dir, "gene_fusions")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

full_out <- file.path(out_dir, paste0(sample_id, ".fusion_annotation_full.tsv"))
write.table(merged, full_out, sep = "\t", row.names = FALSE, quote = FALSE)

fusions <- merged[merged$fusion_type == "gene_fusion", ]
fusions_out <- file.path(out_dir, paste0(sample_id, ".gene_fusions.tsv"))
write.table(fusions, fusions_out, sep = "\t", row.names = FALSE, quote = FALSE)

disruptions <- merged[merged$fusion_type == "gene_disruption", ]
disruptions_out <- file.path(out_dir, paste0(sample_id, ".gene_disruptions.tsv"))
write.table(disruptions, disruptions_out, sep = "\t", row.names = FALSE, quote = FALSE)

n_dup_flagged <- sum(fusions$possible_duplicate_event)

cat("\n=== [05_annotate_variants.R] Done ===\n")
cat("  Confirmed gene fusions      :", nrow(fusions), " ->", fusions_out, "\n")
cat("    of which flagged as possible duplicate/reciprocal calls:", n_dup_flagged, "\n")
cat("  Gene disruptions (non-fusion, single-gene breakpoint hit) :", nrow(disruptions), " ->", disruptions_out, "\n")
cat("  Full annotated BND table (all classifications, for audit):", nrow(merged), " ->", full_out, "\n")

# Append to the QC summary from step 3, if present
qc_summary <- file.path(output_dir, "qc_summary", paste0(sample_id, ".filtering_summary.tsv"))
if (file.exists(qc_summary)) {
  cat(sprintf("gene_fusions_confirmed\tSnpEff-classified fusion\t%d\n", nrow(fusions)),
      file = qc_summary, append = TRUE)
  cat(sprintf("gene_disruptions\tSnpEff-classified ablation\t%d\n", nrow(disruptions)),
      file = qc_summary, append = TRUE)
}
