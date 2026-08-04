# =============================================================================
# R/lib/lib_common.R
# Shared helpers sourced by dependency-light R stages (base graphics/stats
# only, no packages) that read the pipeline YAML config themselves rather than
# via the `yaml` package.
#
# NOT a package, not installed, only ever source()'d. That mirrors
# common/lib_common.sh: these stages have to run unattended on an HPC where
# the R library path may differ from the login shell's, and a stage that fails
# because a package cannot be found has wasted the whole run. Base R plus a
# small regex reader is always available.
#
# WHY THIS EXISTS
# pipelines/methylation/scripts/R/07_generate_report.R and
# 09_compare_samples.R each defined their own copy of this exact reader before
# this file existed. That is the same divergence risk common/lib_common.sh was
# written to avoid on the bash side -- a fix landing in one copy and not the
# other.
#
# NOTE: pipelines/sv/scripts/R/*.R deliberately do NOT use this. The SV
# pipeline already declares `yaml` as a required R package (see
# scripts/check_dependencies.sh's SV_R_PACKAGES) and reads its config with
# yaml.load_file(), which is the more correct parser when a package
# dependency is acceptable. This file is only for the packageless stages.
# =============================================================================

# -----------------------------------------------------------------------------
# cfg_get_factory(cfg_lines)
# Returns a cfg_get(key, default = NULL) closure over the given config lines
# (as read by readLines()). Mirrors common_yaml_get() in common/lib_common.sh:
# first line matching "^\\s*key:", quotes and trailing comment/space stripped.
#
# CONSEQUENCE, same as the bash reader: returns the FIRST match anywhere in
# the file, so every key read this way must be unique across the file (the
# bash side enforces this with common_guard_duplicate_keys before any R stage
# runs).
# -----------------------------------------------------------------------------
cfg_get_factory <- function(cfg_lines) {
  function(key, default = NULL) {
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
}

# -----------------------------------------------------------------------------
# cfg_num_factory(cfg_get)
# Returns a cfg_num(key, default) closure that reads a key via the given
# cfg_get and coerces it to numeric, defaulting cleanly when the key is
# absent.
# -----------------------------------------------------------------------------
cfg_num_factory <- function(cfg_get) {
  function(key, default) as.numeric(cfg_get(key, as.character(default)))
}
