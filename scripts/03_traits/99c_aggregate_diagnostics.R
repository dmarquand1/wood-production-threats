rm(list = ls())

# ==============================================================================
# AGGREGATE FULL MCMC DIAGNOSTICS ACROSS ALL MAIN TAXA
# Reads outputs/traits/diagnostics/<taxon>/full_diagnostic_summary.csv (written
# by 99d_full_diagnostics.R) for each taxon that has one, and combines them
# into one master spreadsheet. Taxa without a full_diagnostic_summary.csv yet
# (e.g. pending a re-fit) are listed separately rather than silently dropped.
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(writexl)

taxa <- c("amphibian", "bird", "fish", "mammal", "reptile")
diag_dir <- file.path("outputs", "traits", "diagnostics")

available <- character(0)
pending <- character(0)

rows <- lapply(taxa, function(taxon) {
  path <- file.path(diag_dir, taxon, "full_diagnostic_summary.csv")
  if (!file.exists(path)) {
    pending <<- c(pending, taxon)
    return(NULL)
  }
  available <<- c(available, taxon)
  read.csv(path, stringsAsFactors = FALSE)
})

combined <- do.call(rbind, rows)
combined <- combined[order(combined$taxon, combined$model), ]

round_cols <- c("max_rhat", "min_bulk_ess", "min_tail_ess", "min_ebfmi",
                "mcse_q2.5_max", "mcse_q97.5_max", "mcse_ratio_max",
                "prop_threatened", "max_pareto_k", "r2_marginal",
                "r2_conditional", "brier_score", "auc")
for (col in round_cols) {
  if (col %in% names(combined)) combined[[col]] <- round(combined[[col]], 4)
}

# Pull run_config out into its own note sheet (one row per taxon/model is
# already carried in the main table's run_config column; this just gives a
# quick "did everything use the same MCMC settings" view up top).
config_note <- combined[, c("taxon", "model", "run_config")]

notes <- data.frame(
  note = c(
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
    sprintf("Taxa included: %s", paste(available, collapse = ", ")),
    if (length(pending) > 0) sprintf("Taxa PENDING (no full_diagnostic_summary.csv yet): %s", paste(pending, collapse = ", ")) else "All 5 main taxa included.",
    "",
    "Column notes:",
    "max_rhat / min_bulk_ess / min_tail_ess: split-Rhat and ESS-ratio (ESS/total post-warmup draws), computed via posterior::summarise_draws() -- bulk and tail are genuinely separate (see *_param columns for which parameter drives each).",
    "min_ebfmi: minimum E-BFMI across chains; values below ~0.2 indicate energy-exploration problems.",
    "mcse_ratio_max: worst-case MCSE / posterior SD among fixed effects (b_) and group-level SDs (sd_); rule of thumb is <10%.",
    "max_pareto_k / n_pareto_k_above_0.7: PSIS-LOO Pareto k-hat diagnostics; k>0.7 flags an unstable/influential observation.",
    "r2_marginal / r2_conditional: Nakagawa/Schielzeth-style Bayesian R2 partition (fixed effects only vs fixed+phylogeny+spatial), median of posterior draws -- see <taxon>/<model>_r2_bayes_partition.csv for full CI.",
    "brier_score / auc: in-sample calibration from posterior_epred (not held-out/cross-validated).",
    "passed: n_divergences==0 & max_rhat<1.01 & min_bulk_ess>0.1 & min_tail_ess>0.1.",
    "run_config: exact chains/iter/warmup/draws/adapt_delta/max_treedepth/seeds/package versions, read directly off the fitted object."
  ),
  stringsAsFactors = FALSE
)

out_csv <- file.path(diag_dir, "all_taxa_model_diagnostics.csv")
out_xlsx <- file.path(diag_dir, "all_taxa_model_diagnostics.xlsx")

write.csv(combined, out_csv, row.names = FALSE)
write_xlsx(list(diagnostics = combined, run_config = config_note, notes = notes), out_xlsx)

cat("\nWrote:\n  ", out_csv, "\n  ", out_xlsx, "\n\n", sep = "")
if (length(pending) > 0) {
  cat("PENDING (not included): ", paste(pending, collapse = ", "), "\n\n", sep = "")
}
print(combined[, c("taxon", "model", "max_rhat", "min_bulk_ess", "min_tail_ess", "n_divergences", "passed")], row.names = FALSE)
