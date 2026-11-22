rm(list = ls())
gc()

# ==============================================================================
# COLLECT FULL MCMC DIAGNOSTICS (Rhat / ESS / divergences) FOR TAXA WHOSE
# diagnostic_summary.csv IS MISSING OR INCOMPLETE.
#
# Unlike scripts/03_traits/99_model_diagnostics.R (which rebuilds the combined
# fit from 4 separate per-chain .rds files), this script loads the ALREADY
# COMBINED fit objects directly:
#   outputs/traits/models/<taxon>/logging_combined.rds
#   outputs/traits/models/<taxon>/plantation_combined.rds
# and computes the same diagnostics (max Rhat, min bulk/tail ESS ratio,
# divergent transitions). This works for both the main taxon fits and the
# *_cr_en_vu_nt sensitivity fits, whose combined objects were produced from a
# single 4-chain Stan run rather than brms::combine_models().
#
# Usage: Rscript scripts/03_traits/99b_diagnostics_from_combined.R <taxon> [<taxon2> ...]
#   e.g. Rscript scripts/03_traits/99b_diagnostics_from_combined.R amphibian mammal
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Provide at least one taxon directory name, e.g. amphibian mammal amphibian_cr_en_vu_nt")
}
taxa <- tolower(args)

run_diagnostics <- function(fit, taxon, model_name) {
  cat(sprintf(
    "\n%s\n  TAXON: %s\n  MODEL: %s\n%s\n",
    strrep("=", 60), toupper(taxon), toupper(model_name), strrep("=", 60)
  ))

  rhat_vals <- rhat(fit)
  neff_bulk <- neff_ratio(fit, pars = variables(fit))
  neff_tail <- neff_ratio(fit, type = "tail")

  total_divergences <- tryCatch({
    sampler_params <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
    sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
  }, error = function(e) {
    warning(sprintf("Could not extract sampler params for %s %s: %s", taxon, model_name, e$message))
    NA_integer_
  })

  max_rhat <- max(rhat_vals, na.rm = TRUE)
  min_bulk <- min(neff_bulk, na.rm = TRUE)
  min_tail <- min(neff_tail, na.rm = TRUE)

  cat(sprintf("Max R-hat: %.4f\n", max_rhat))
  cat(sprintf("Min bulk ESS ratio: %.3f\n", min_bulk))
  cat(sprintf("Min tail ESS ratio: %.3f\n", min_tail))
  cat(sprintf("Total divergent transitions: %s\n", ifelse(is.na(total_divergences), "NA", total_divergences)))

  passed <- !is.na(total_divergences) &&
    total_divergences == 0 &&
    max_rhat < 1.01 &&
    min_bulk > 0.1 &&
    min_tail > 0.1

  data.frame(
    taxon = taxon,
    model = model_name,
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk,
    min_tail_ess = min_tail,
    n_divergences = total_divergences,
    passed = passed
  )
}

for (taxon in taxa) {
  MODEL_DIR <- file.path("outputs", "traits", "models", taxon)
  DIAG_DIR  <- file.path("outputs", "traits", "diagnostics", taxon)
  dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

  rows <- list()

  for (model_name in c("logging", "plantation")) {
    fit_path <- file.path(MODEL_DIR, sprintf("%s_combined.rds", model_name))
    if (!file.exists(fit_path)) {
      warning(sprintf("Missing combined fit, skipping: %s", fit_path))
      next
    }
    message(sprintf("\nLoading %s / %s ...", taxon, model_name))
    fit <- readRDS(fit_path)
    rows[[model_name]] <- run_diagnostics(fit, taxon, model_name)
    rm(fit)
    gc()
  }

  if (length(rows) > 0) {
    summary_df <- do.call(rbind, rows)
    out_path <- file.path(DIAG_DIR, "diagnostic_summary.csv")
    write.csv(summary_df, out_path, row.names = FALSE)
    message(sprintf("Wrote %s", out_path))
  }
}

message("\nDone.")
