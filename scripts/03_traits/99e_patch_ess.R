rm(list = ls())
gc()

# ==============================================================================
# LIGHTWEIGHT PATCH: recompute correct bulk/tail ESS for a taxon whose
# full_diagnostic_summary.csv was written before the neff_ratio() bug fix
# (see 99d_full_diagnostics.R). Reloads the combined fits only to recompute
# min_bulk_ess / min_tail_ess (+ responsible parameter), and patches those
# four columns into the existing CSV in place. Does NOT redo Rhat, treedepth,
# E-BFMI, MCSE, prevalence, Pareto k, R2, or calibration -- none of those were
# affected by the bug.
#
# Usage: Rscript scripts/03_traits/99e_patch_ess.R <taxon>
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 99e_patch_ess.R <taxon>")
taxon <- tolower(args[[1]])

MODEL_DIR <- file.path("outputs", "traits", "models", taxon)
DIAG_DIR  <- file.path("outputs", "traits", "diagnostics", taxon)
summary_path <- file.path(DIAG_DIR, "full_diagnostic_summary.csv")

if (!file.exists(summary_path)) {
  stop(sprintf("No existing full_diagnostic_summary.csv for %s -- run 99d_full_diagnostics.R instead.", taxon))
}
summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)

recompute_ess <- function(model_name) {
  fit_path <- file.path(MODEL_DIR, sprintf("%s_combined.rds", model_name))
  message(sprintf("[%s / %s] Loading %s ...", taxon, model_name, fit_path))
  fit <- readRDS(fit_path)

  draws_all <- as_draws_array(fit, variable = variables(fit))
  ess_tab <- posterior::summarise_draws(draws_all, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  n_draws_all <- posterior::ndraws(fit)
  ess_bulk_ratio <- setNames(ess_tab$ess_bulk / n_draws_all, ess_tab$variable)
  ess_tail_ratio <- setNames(ess_tab$ess_tail / n_draws_all, ess_tab$variable)
  rm(draws_all, fit)
  gc()

  min_bulk <- min(ess_bulk_ratio, na.rm = TRUE)
  min_bulk_param <- names(ess_bulk_ratio)[which.min(ess_bulk_ratio)]
  min_tail <- min(ess_tail_ratio, na.rm = TRUE)
  min_tail_param <- names(ess_tail_ratio)[which.min(ess_tail_ratio)]

  cat(sprintf(
    "[%s / %s] Corrected: min bulk ESS ratio %.3f (%s) | min tail ESS ratio %.3f (%s)\n",
    taxon, model_name, min_bulk, min_bulk_param, min_tail, min_tail_param
  ))

  list(min_bulk = min_bulk, min_bulk_param = min_bulk_param,
       min_tail = min_tail, min_tail_param = min_tail_param)
}

for (model_name in unique(summary_df$model)) {
  corrected <- recompute_ess(model_name)
  idx <- which(summary_df$model == model_name)
  summary_df$min_bulk_ess[idx] <- corrected$min_bulk
  summary_df$min_bulk_ess_param[idx] <- corrected$min_bulk_param
  summary_df$min_tail_ess[idx] <- corrected$min_tail
  summary_df$min_tail_ess_param[idx] <- corrected$min_tail_param

  # passed flag depends on min ESS too -- recheck it
  summary_df$passed[idx] <- !is.na(summary_df$n_divergences[idx]) &
    summary_df$n_divergences[idx] == 0 &
    summary_df$max_rhat[idx] < 1.01 &
    corrected$min_bulk > 0.1 &
    corrected$min_tail > 0.1
}

write.csv(summary_df, summary_path, row.names = FALSE)
message(sprintf("\nPatched %s", summary_path))
