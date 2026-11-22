rm(list = ls())
gc()

# ==============================================================================
# FULL MCMC / MODEL-ADEQUACY DIAGNOSTICS FOR A SINGLE TAXON
#
# Extends 99b_diagnostics_from_combined.R with:
#   - which parameter attains max Rhat / min bulk ESS / min tail ESS
#   - max-treedepth saturations, min E-BFMI across chains
#   - MCSE of the 2.5% / 97.5% posterior quantiles for fixed effects (b_) and
#     group-level SDs (sd_), and the worst MCSE/SD ratio (rule of thumb: <10%)
#   - prevalence (n species, n threatened, proportion) and grouping-factor
#     level counts (phylo_name, spatial_name)
#   - Pareto k-hat (max, and count > 0.7) from PSIS-LOO
#   - in-sample calibration (Brier score, AUC) from posterior_epred
#   - the actual run configuration read off the fit object (chains, iter,
#     warmup, post-warmup draws, adapt_delta, max_treedepth, seeds, versions)
#
# Loads outputs/traits/models/<taxon>/{logging,plantation}_combined.rds
# directly (does not refit / recombine chains).
#
# Usage: Rscript scripts/03_traits/99d_full_diagnostics.R <taxon>
#   (one taxon per process invocation, by design, to bound peak memory)
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 99d_full_diagnostics.R <taxon>  (exactly one taxon)")
}
taxon <- tolower(args[[1]])

MODEL_DIR <- file.path("outputs", "traits", "models", taxon)
DIAG_DIR  <- file.path("outputs", "traits", "diagnostics", taxon)
dir.create(DIAG_DIR, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

diagnose_one <- function(taxon, model_name) {
  fit_path <- file.path(MODEL_DIR, sprintf("%s_combined.rds", model_name))
  if (!file.exists(fit_path)) {
    warning(sprintf("Missing combined fit, skipping: %s", fit_path))
    return(NULL)
  }

  message(sprintf("\n[%s / %s] Loading %s ...", taxon, model_name, fit_path))
  fit <- readRDS(fit_path)

  cat(sprintf(
    "\n%s\n  TAXON: %s\n  MODEL: %s\n%s\n",
    strrep("=", 60), toupper(taxon), toupper(model_name), strrep("=", 60)
  ))

  # --- Rhat / ESS, with the responsible parameter name ------------------------
  # NOTE: brms::neff_ratio.brmsfit() has no real `type` argument -- it always
  # computes both ess_bulk and ess_tail internally and silently collapses them
  # via rowMins() into a single number, so calling it twice (once "for bulk",
  # once with type = "tail") returns the *same* combined value both times.
  # We instead call posterior::summarise_draws() ourselves so bulk and tail
  # are genuinely separate diagnostics.
  rhat_vals <- rhat(fit)
  max_rhat <- max(rhat_vals, na.rm = TRUE)
  max_rhat_param <- names(rhat_vals)[which.max(rhat_vals)]

  draws_all <- as_draws_array(fit, variable = variables(fit))
  ess_tab <- posterior::summarise_draws(draws_all, ess_bulk = posterior::ess_bulk, ess_tail = posterior::ess_tail)
  n_draws_all <- posterior::ndraws(fit)
  ess_bulk_ratio <- setNames(ess_tab$ess_bulk / n_draws_all, ess_tab$variable)
  ess_tail_ratio <- setNames(ess_tab$ess_tail / n_draws_all, ess_tab$variable)
  rm(draws_all)
  gc()

  min_bulk <- min(ess_bulk_ratio, na.rm = TRUE)
  min_bulk_param <- names(ess_bulk_ratio)[which.min(ess_bulk_ratio)]

  min_tail <- min(ess_tail_ratio, na.rm = TRUE)
  min_tail_param <- names(ess_tail_ratio)[which.min(ess_tail_ratio)]

  # --- Run configuration (read off the fit, not the launch script) -----------
  stan_args <- fit$fit@stan_args
  chains <- length(stan_args)
  iter <- stan_args[[1]]$iter
  warmup <- stan_args[[1]]$warmup
  seeds <- paste(vapply(stan_args, function(a) a$seed, numeric(1)), collapse = ",")
  ctrl <- stan_args[[1]]$control
  adapt_delta <- ctrl$adapt_delta %||% 0.8
  max_treedepth_setting <- ctrl$max_treedepth %||% 10
  post_warmup_draws <- posterior::ndraws(fit)

  # --- Sampler diagnostics: divergences, treedepth saturations, E-BFMI -------
  sampler_params <- tryCatch(
    rstan::get_sampler_params(fit$fit, inc_warmup = FALSE),
    error = function(e) NULL
  )
  if (!is.null(sampler_params)) {
    n_divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
    n_treedepth_sat <- sum(vapply(
      sampler_params,
      function(x) sum(x[, "treedepth__"] >= max_treedepth_setting),
      numeric(1)
    ))
  } else {
    n_divergences <- NA_integer_
    n_treedepth_sat <- NA_integer_
  }

  bfmi_vals <- tryCatch(rstan::get_bfmi(fit$fit), error = function(e) NA_real_)
  min_ebfmi <- if (all(is.na(bfmi_vals))) NA_real_ else min(bfmi_vals, na.rm = TRUE)

  cat(sprintf("Max R-hat: %.4f (%s)\n", max_rhat, max_rhat_param))
  cat(sprintf("Min bulk ESS ratio: %.3f (%s)\n", min_bulk, min_bulk_param))
  cat(sprintf("Min tail ESS ratio: %.3f (%s)\n", min_tail, min_tail_param))
  cat(sprintf("Divergences: %s | Treedepth saturations: %s | Min E-BFMI: %.3f\n",
              n_divergences, n_treedepth_sat, min_ebfmi))

  passed <- !is.na(n_divergences) &&
    n_divergences == 0 &&
    max_rhat < 1.01 &&
    min_bulk > 0.1 &&
    min_tail > 0.1

  # --- MCSE of posterior quantiles for b_ / sd_ terms -------------------------
  message(sprintf("[%s / %s] Computing MCSE on quantiles for fixed/random-SD terms...", taxon, model_name))
  key_vars <- grep("^b_|^sd_", variables(fit), value = TRUE)
  draws_key <- as_draws_df(fit, variable = key_vars)
  param_cols <- setdiff(names(draws_key), c(".chain", ".iteration", ".draw"))

  mcse_tab <- do.call(rbind, lapply(param_cols, function(p) {
    x <- draws_key[[p]]
    sdx <- sd(x)
    mcse_q <- posterior::mcse_quantile(x, probs = c(0.025, 0.975))
    data.frame(
      param = p, sd = sdx,
      mcse_q2.5 = unname(mcse_q[1]), mcse_q97.5 = unname(mcse_q[2]),
      ratio_q2.5 = unname(mcse_q[1]) / sdx, ratio_q97.5 = unname(mcse_q[2]) / sdx
    )
  }))
  mcse_tab$max_ratio <- pmax(mcse_tab$ratio_q2.5, mcse_tab$ratio_q97.5)

  write.csv(
    mcse_tab,
    file.path(DIAG_DIR, sprintf("%s_mcse_detail.csv", model_name)),
    row.names = FALSE
  )

  worst_row <- mcse_tab[which.max(mcse_tab$max_ratio), ]
  mcse_q2.5_max <- max(mcse_tab$mcse_q2.5)
  mcse_q2.5_max_param <- mcse_tab$param[which.max(mcse_tab$mcse_q2.5)]
  mcse_q97.5_max <- max(mcse_tab$mcse_q97.5)
  mcse_q97.5_max_param <- mcse_tab$param[which.max(mcse_tab$mcse_q97.5)]
  mcse_ratio_max <- worst_row$max_ratio
  mcse_ratio_max_param <- worst_row$param

  rm(draws_key)
  gc()

  # --- Prevalence and grouping-factor level counts ----------------------------
  resp_var <- all.vars(fit$formula$formula)[1]
  y <- fit$data[[resp_var]]
  y_num <- if (is.logical(y)) {
    as.numeric(y)
  } else if (is.factor(y) || is.character(y)) {
    as.numeric(as.character(y))
  } else {
    as.numeric(y)
  }
  n_obs <- nrow(fit$data)
  n_threatened <- sum(y_num == 1, na.rm = TRUE)
  prop_threatened <- n_threatened / n_obs

  group_vars <- unique(fit$ranef$group)
  n_phylo_groups <- if ("phylo_name" %in% group_vars) length(unique(fit$data$phylo_name)) else NA_integer_
  n_spatial_groups <- if ("spatial_name" %in% group_vars) length(unique(fit$data$spatial_name)) else NA_integer_

  # --- Pareto k-hat from PSIS-LOO ----------------------------------------------
  message(sprintf("[%s / %s] Computing PSIS-LOO (Pareto k-hat)...", taxon, model_name))
  loo_fit <- tryCatch(loo::loo(fit), error = function(e) {
    warning(sprintf("loo() failed for %s %s: %s", taxon, model_name, e$message))
    NULL
  })
  if (!is.null(loo_fit)) {
    pareto_k <- loo_fit$diagnostics$pareto_k
    max_pareto_k <- max(pareto_k, na.rm = TRUE)
    n_pareto_k_above_0.7 <- sum(pareto_k > 0.7, na.rm = TRUE)
  } else {
    max_pareto_k <- NA_real_
    n_pareto_k_above_0.7 <- NA_integer_
  }
  rm(loo_fit)
  gc()

  # --- In-sample calibration: Brier score and AUC ------------------------------
  message(sprintf("[%s / %s] Computing calibration (Brier / AUC)...", taxon, model_name))
  calib <- tryCatch({
    epred <- posterior_epred(fit)
    mean_prob <- colMeans(epred)
    rm(epred)
    gc()
    brier <- mean((mean_prob - y_num)^2)
    pos <- mean_prob[y_num == 1]
    neg <- mean_prob[y_num == 0]
    if (length(pos) > 0 && length(neg) > 0) {
      r <- rank(c(pos, neg))
      auc <- (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
        (length(pos) * length(neg))
    } else {
      auc <- NA_real_
    }
    list(brier = brier, auc = auc)
  }, error = function(e) {
    warning(sprintf("Calibration failed for %s %s: %s", taxon, model_name, e$message))
    list(brier = NA_real_, auc = NA_real_)
  })

  # --- Bayesian R2 partition (reuse existing partition output, not recomputed) -
  r2_path <- file.path(DIAG_DIR, sprintf("%s_r2_bayes_partition.csv", model_name))
  r2_marginal <- NA_real_
  r2_conditional <- NA_real_
  if (file.exists(r2_path)) {
    r2_tab <- read.csv(r2_path, stringsAsFactors = FALSE)
    r2_marginal <- r2_tab$median[r2_tab$component == "Fixed effects"]
    r2_conditional <- r2_tab$median[r2_tab$component == "Conditional (total)"]
  }

  run_config <- sprintf(
    "chains=%d, iter=%d, warmup=%d, post_warmup_draws=%d, adapt_delta=%s, max_treedepth=%s, seeds=%s, brms=%s, rstan=%s, Stan=%s",
    chains, iter, warmup, post_warmup_draws,
    format(adapt_delta), format(max_treedepth_setting), seeds,
    fit$version$brms %||% NA, fit$version$rstan %||% NA,
    tryCatch(as.character(rstan::stan_version()), error = function(e) NA)
  )

  result <- data.frame(
    taxon = taxon,
    model = model_name,
    max_rhat = max_rhat,
    max_rhat_param = max_rhat_param,
    min_bulk_ess = min_bulk,
    min_bulk_ess_param = min_bulk_param,
    min_tail_ess = min_tail,
    min_tail_ess_param = min_tail_param,
    n_divergences = n_divergences,
    max_treedepth_saturations = n_treedepth_sat,
    min_ebfmi = min_ebfmi,
    mcse_q2.5_max = mcse_q2.5_max,
    mcse_q2.5_max_param = mcse_q2.5_max_param,
    mcse_q97.5_max = mcse_q97.5_max,
    mcse_q97.5_max_param = mcse_q97.5_max_param,
    mcse_ratio_max = mcse_ratio_max,
    mcse_ratio_max_param = mcse_ratio_max_param,
    n_species = n_obs,
    n_threatened = n_threatened,
    prop_threatened = prop_threatened,
    n_phylo_groups = n_phylo_groups,
    n_spatial_groups = n_spatial_groups,
    max_pareto_k = max_pareto_k,
    n_pareto_k_above_0.7 = n_pareto_k_above_0.7,
    r2_marginal = r2_marginal,
    r2_conditional = r2_conditional,
    brier_score = calib$brier,
    auc = calib$auc,
    passed = passed,
    run_config = run_config,
    stringsAsFactors = FALSE
  )

  rm(fit)
  gc()

  # Checkpoint: write/replace this model's row immediately so an interrupted
  # run doesn't lose an already-finished model's results.
  out_path <- file.path(DIAG_DIR, "full_diagnostic_summary.csv")
  if (file.exists(out_path)) {
    existing <- read.csv(out_path, stringsAsFactors = FALSE)
    existing <- existing[existing$model != model_name, , drop = FALSE]
    summary_df <- rbind(existing, result)
  } else {
    summary_df <- result
  }
  write.csv(summary_df, out_path, row.names = FALSE)
  message(sprintf("[%s / %s] Checkpointed to %s", taxon, model_name, out_path))

  result
}

for (m in c("logging", "plantation")) {
  diagnose_one(taxon, m)
}

message("\nDone.")
