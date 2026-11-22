rm(list = ls())
gc()

library(brms)
library(tidyverse)
library(here)

args <- commandArgs(trailingOnly = TRUE)
taxon <- if (length(args) >= 1) tolower(args[[1]]) else "fish"

OUTPUT_DIR <- here("outputs", "traits", "models", taxon)
PLOT_DIR   <- here("outputs", "traits", "diagnostics", taxon)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

load_and_combine <- function(prefix, n_chains = 4) {
  message(sprintf("\n--- Loading %s %s chains ---", taxon, prefix))
  
  fits <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    path <- file.path(OUTPUT_DIR, sprintf("%s_%s_chain%d.rds", taxon, prefix, i))
    if (!file.exists(path)) {
      stop(sprintf("Missing file: %s", path))
    }
    fits[[i]] <- readRDS(path)
    message(sprintf("  Loaded chain %d", i))
  }
  
  message("  Combining chains...")
  combine_models(fits[[1]], fits[[2]], fits[[3]], fits[[4]])
}

run_diagnostics <- function(fit, model_name) {
  cat(sprintf(
    "\n%s\n  TAXON: %s\n  MODEL: %s\n%s\n",
    strrep("=", 60), toupper(taxon), toupper(model_name), strrep("=", 60)
  ))
  
  rhat_vals <- rhat(fit)
  neff_bulk <- neff_ratio(fit, pars = variables(fit))
  neff_tail <- neff_ratio(fit, type = "tail")
  sampler_params <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
  
  total_divergences <- sum(vapply(
    sampler_params,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))
  
  cat(sprintf("Max R-hat: %.4f\n", max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("Min bulk ESS ratio: %.3f\n", min(neff_bulk, na.rm = TRUE)))
  cat(sprintf("Min tail ESS ratio: %.3f\n", min(neff_tail, na.rm = TRUE)))
  cat(sprintf("Total divergent transitions: %d\n", total_divergences))
  
  passed <- total_divergences == 0 &&
    max(rhat_vals, na.rm = TRUE) < 1.01 &&
    min(neff_bulk, na.rm = TRUE) > 0.1 &&
    min(neff_tail, na.rm = TRUE) > 0.1
  
  invisible(list(
    max_rhat = max(rhat_vals, na.rm = TRUE),
    min_bulk_ess = min(neff_bulk, na.rm = TRUE),
    min_tail_ess = min(neff_tail, na.rm = TRUE),
    n_divergences = total_divergences,
    passed = passed
  ))
}

save_trace_plots <- function(fit, model_name) {
  p_fixed <- plot(fit, variable = "^b_", regex = TRUE, ask = FALSE)
  ggsave(
    file.path(PLOT_DIR, sprintf("%s_trace_fixed_effects.png", model_name)),
    plot = p_fixed[[1]],
    width = 12, height = 10, dpi = 150
  )
  
  p_sd <- plot(fit, variable = "^sd_", regex = TRUE, ask = FALSE)
  ggsave(
    file.path(PLOT_DIR, sprintf("%s_trace_random_sd.png", model_name)),
    plot = p_sd[[1]],
    width = 8, height = 6, dpi = 150
  )
}

fit_logging <- load_and_combine("logging")
fit_plantation <- load_and_combine("plantation")

diag_logging <- run_diagnostics(fit_logging, "logging")
diag_plantation <- run_diagnostics(fit_plantation, "plantation")

save_trace_plots(fit_logging, "logging")
save_trace_plots(fit_plantation, "plantation")

summary_df <- data.frame(
  taxon = taxon,
  model = c("logging", "plantation"),
  max_rhat = c(diag_logging$max_rhat, diag_plantation$max_rhat),
  min_bulk_ess = c(diag_logging$min_bulk_ess, diag_plantation$min_bulk_ess),
  min_tail_ess = c(diag_logging$min_tail_ess, diag_plantation$min_tail_ess),
  n_divergences = c(diag_logging$n_divergences, diag_plantation$n_divergences),
  passed = c(diag_logging$passed, diag_plantation$passed)
)

write.csv(
  summary_df,
  file = file.path(PLOT_DIR, "diagnostic_summary.csv"),
  row.names = FALSE
)

saveRDS(fit_logging, file.path(OUTPUT_DIR, "logging_combined.rds"))
saveRDS(fit_plantation, file.path(OUTPUT_DIR, "plantation_combined.rds"))

message("\nCombined fits saved to ", OUTPUT_DIR)
message("Trace plots and summary saved to ", PLOT_DIR)
