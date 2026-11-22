rm(list = ls())
gc()

# =============================================================================
# COMPARE FULL TRAIT MODELS TO CR/EN/VU/NT-ONLY MODELS
#
# Reads existing full-model summaries and CR/EN/VU/NT model objects. It does not
# refit models. Filtered model summaries are extracted directly from the saved
# Stan draws, so this script does not need to load brms.
#
# RStudio:
#   Edit RSTUDIO_TAXA_TO_RUN / RSTUDIO_MODEL_TO_RUN, then click Source.
#
# Command line:
#   Rscript scripts/03_traits/compare_no_lc_dd_models.R all both
#   Rscript scripts/03_traits/compare_no_lc_dd_models.R amphibian logging
# =============================================================================

# ---- RSTUDIO SETTINGS --------------------------------------------------------

RSTUDIO_TAXA_TO_RUN <- c("amphibian", "fish", "mammal", "reptile")
# Options: "all", "amphibian", "fish", "mammal", "reptile", "bird",
# or a vector such as c("amphibian", "fish").

RSTUDIO_MODEL_TO_RUN <- "both"
# Options: "both", "logging", "plantation".

# -----------------------------------------------------------------------------

find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "redListAnalysis.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = TRUE)
}

project_root <- find_project_root()
if (!identical(normalizePath(getwd(), winslash = "/", mustWork = TRUE), project_root)) {
  setwd(project_root)
  cat("Working directory set to: ", project_root, "\n", sep = "")
}

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(dplyr)
library(ggplot2)

default_taxa <- c("amphibian", "fish", "mammal", "reptile")
all_taxa <- c(default_taxa, "bird")
valid_models <- c("logging", "plantation", "both")

normalise_taxa <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE))
  x <- tolower(trimws(x))
  x <- x[x != ""]
  if (length(x) == 0 || any(x == "all")) return(default_taxa)
  x
}

args <- commandArgs(trailingOnly = TRUE)
use_settings <- interactive() || length(args) == 0

taxa_to_run <- if (use_settings) normalise_taxa(RSTUDIO_TAXA_TO_RUN) else normalise_taxa(args[[1]])
model_to_run <- if (use_settings) {
  tolower(RSTUDIO_MODEL_TO_RUN)
} else if (length(args) >= 2) {
  tolower(args[[2]])
} else {
  "both"
}

unknown_taxa <- setdiff(taxa_to_run, all_taxa)
if (length(unknown_taxa) > 0) {
  stop("Unknown taxon/taxa: ", paste(unknown_taxa, collapse = ", "), call. = FALSE)
}
if (!model_to_run %in% valid_models) {
  stop("Model must be one of: ", paste(valid_models, collapse = ", "), call. = FALSE)
}

models_to_run <- if (model_to_run == "both") c("logging", "plantation") else model_to_run

metric_value <- function(summary_df, metric_name, default = NA_character_) {
  if (nrow(summary_df) == 0) return(default)
  value <- summary_df$value[summary_df$metric == metric_name]
  if (length(value) == 0) default else value[[1]]
}

first_metric_value <- function(summary_df, metric_names, default = NA_character_) {
  for (metric_name in metric_names) {
    value <- metric_value(summary_df, metric_name, default = NA_character_)
    if (!is.na(value)) return(value)
  }
  default
}

as_binary_true <- function(x) {
  x %in% c(TRUE, "TRUE", "True", "true", 1, "1")
}

find_filtered_model_file <- function(taxon, model_name) {
  model_dir <- file.path("outputs", "traits", "models", paste0(taxon, "_cr_en_vu_nt"))
  candidates <- c(
    file.path(model_dir, sprintf("%s_combined.rds", model_name)),
    file.path(model_dir, sprintf("%s_cr_en_vu_nt_%s_pc_4chains.rds", taxon, model_name))
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) NA_character_ else existing[[1]]
}

extract_stan_sim <- function(model_file) {
  fit_obj <- readRDS(model_file)
  stan_fit <- fit_obj[["fit"]]
  sim <- attr(stan_fit, "sim")
  if (is.null(sim) || is.null(sim$samples)) {
    stop("Could not find Stan samples in: ", model_file, call. = FALSE)
  }
  list(
    samples = sim$samples,
    warmup = sim$warmup2,
    n_obs = tryCatch(nrow(fit_obj$data), error = function(e) NA_integer_),
    family = tryCatch(fit_obj$family$family, error = function(e) NA_character_),
    link = tryCatch(fit_obj$family$link, error = function(e) NA_character_)
  )
}

param_draws_by_chain <- function(sim, param) {
  lapply(seq_along(sim$samples), function(chain_id) {
    draws <- sim$samples[[chain_id]][[param]]
    if (is.null(draws)) stop("Missing parameter in Stan draws: ", param, call. = FALSE)
    warmup <- sim$warmup[[chain_id]]
    if (length(draws) <= warmup) stop("No post-warmup draws for: ", param, call. = FALSE)
    draws[(warmup + 1):length(draws)]
  })
}

summarise_param <- function(sim, param) {
  draws <- unlist(param_draws_by_chain(sim, param), use.names = FALSE)
  data.frame(
    param = param,
    estimate = mean(draws),
    std_error = sd(draws),
    q2.5 = unname(quantile(draws, 0.025)),
    q97.5 = unname(quantile(draws, 0.975))
  )
}

simple_rhat <- function(sim, param) {
  chains <- param_draws_by_chain(sim, param)
  n <- min(vapply(chains, length, integer(1)))
  if (length(chains) < 2 || n < 2) return(NA_real_)
  chains <- lapply(chains, function(x) x[seq_len(n)])
  w <- mean(vapply(chains, var, numeric(1)))
  b <- n * var(vapply(chains, mean, numeric(1)))
  if (!is.finite(w) || w == 0) return(NA_real_)
  sqrt((((n - 1) / n) * w + b / n) / w)
}

param_to_term <- function(param) {
  if (identical(param, "b_Intercept")) return("Intercept")
  sub("^b_", "", param)
}

summarise_filtered_model <- function(taxon, model_name, model_file, out_dir) {
  sim <- extract_stan_sim(model_file)
  param_names <- names(sim$samples[[1]])
  fixed_params <- grep("^b_", param_names, value = TRUE)
  random_params <- grep("^sd_", param_names, value = TRUE)
  
  fixed <- bind_rows(lapply(fixed_params, function(p) summarise_param(sim, p))) %>%
    mutate(term = vapply(param, param_to_term, character(1))) %>%
    transmute(
      taxon = paste0(taxon, "_cr_en_vu_nt"),
      model = model_name,
      term,
      estimate,
      std_error,
      q2.5,
      q97.5,
      direction = case_when(
        q2.5 > 0 ~ "positive",
        q97.5 < 0 ~ "negative",
        TRUE ~ "uncertain"
      )
    )
  
  random <- if (length(random_params) > 0) {
    bind_rows(lapply(random_params, function(p) summarise_param(sim, p))) %>%
      transmute(
        taxon = paste0(taxon, "_cr_en_vu_nt"),
        model = model_name,
        term = param,
        group = case_when(
          grepl("phylo_name", param) ~ "phylo_name",
          grepl("spatial_name", param) ~ "spatial_name",
          TRUE ~ NA_character_
        ),
        effect_type = "sd",
        estimate,
        std_error,
        q2.5,
        q97.5
      )
  } else {
    data.frame()
  }
  
  diagnostics <- data.frame(
    taxon = paste0(taxon, "_cr_en_vu_nt"),
    model = model_name,
    max_rhat = suppressWarnings(max(vapply(c(fixed_params, random_params), function(p) simple_rhat(sim, p), numeric(1)), na.rm = TRUE)),
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    n_divergences = NA_integer_,
    passed = NA
  )
  if (!is.finite(diagnostics$max_rhat)) diagnostics$max_rhat <- NA_real_
  
  overall <- data.frame(
    taxon = paste0(taxon, "_cr_en_vu_nt"),
    model = model_name,
    family = sim$family,
    link = sim$link,
    n_obs = sim$n_obs
  )
  
  write.csv(fixed, file.path(out_dir, sprintf("%s_fixed_effects_summary.csv", model_name)), row.names = FALSE)
  write.csv(random, file.path(out_dir, sprintf("%s_random_effects_summary.csv", model_name)), row.names = FALSE)
  write.csv(overall, file.path(out_dir, sprintf("%s_overall_fit_summary.csv", model_name)), row.names = FALSE)
  
  list(fixed = fixed, random = random, diagnostics = diagnostics, overall = overall)
}

compare_one_model <- function(taxon, model_name) {
  cat("\n========================================\n")
  cat("Comparing ", taxon, " ", model_name, "\n", sep = "")
  cat("========================================\n")
  
  full_diag_dir <- file.path("outputs", "traits", "diagnostics", taxon)
  filtered_diag_dir <- file.path("outputs", "traits", "diagnostics", paste0(taxon, "_cr_en_vu_nt"))
  figure_dir <- file.path("outputs", "traits", "figures")
  dir.create(filtered_diag_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  
  full_fixed_file <- file.path(full_diag_dir, sprintf("%s_fixed_effects_summary.csv", model_name))
  full_random_file <- file.path(full_diag_dir, sprintf("%s_random_effects_summary.csv", model_name))
  full_diag_file <- file.path(full_diag_dir, "diagnostic_summary.csv")
  filtered_model_file <- find_filtered_model_file(taxon, model_name)
  
  if (!file.exists(full_fixed_file)) {
    stop("Missing full-model fixed-effect summary: ", full_fixed_file, call. = FALSE)
  }
  if (is.na(filtered_model_file)) {
    stop(
      "Missing filtered model for ",
      taxon,
      " ",
      model_name,
      ". Run fit_pc_no_lc_dd.R for that model first.",
      call. = FALSE
    )
  }
  
  full_fixed <- read.csv(full_fixed_file, stringsAsFactors = FALSE) %>%
    transmute(
      term,
      estimate_with_lc_dd = estimate,
      std_error_with_lc_dd = std_error,
      q2.5_with_lc_dd = q2.5,
      q97.5_with_lc_dd = q97.5,
      direction_with_lc_dd = direction
    )
  
  filtered <- summarise_filtered_model(taxon, model_name, filtered_model_file, filtered_diag_dir)
  
  comparison <- full_fixed %>%
    inner_join(
      filtered$fixed %>%
        transmute(
          term,
          estimate_without_lc_dd = estimate,
          std_error_without_lc_dd = std_error,
          q2.5_without_lc_dd = q2.5,
          q97.5_without_lc_dd = q97.5,
          direction_without_lc_dd = direction
        ),
      by = "term"
    ) %>%
    mutate(
      delta_without_minus_with = estimate_without_lc_dd - estimate_with_lc_dd,
      absolute_delta = abs(delta_without_minus_with),
      credible_intervals_overlap =
        pmax(q2.5_with_lc_dd, q2.5_without_lc_dd) <=
        pmin(q97.5_with_lc_dd, q97.5_without_lc_dd),
      direction_changed = direction_with_lc_dd != direction_without_lc_dd,
      sensitivity_flag = case_when(
        direction_changed ~ "direction_changed",
        !credible_intervals_overlap ~ "intervals_do_not_overlap",
        absolute_delta >= 0.25 ~ "estimate_shift_ge_0.25_log_odds",
        TRUE ~ "similar"
      )
    ) %>%
    arrange(desc(sensitivity_flag != "similar"), desc(absolute_delta))
  
  comparison_file <- file.path(filtered_diag_dir, sprintf("%s_with_vs_without_lc_dd_fixed_effects.csv", model_name))
  write.csv(comparison, comparison_file, row.names = FALSE)
  
  if (file.exists(full_random_file) && nrow(filtered$random) > 0) {
    random_comparison <- read.csv(full_random_file, stringsAsFactors = FALSE) %>%
      transmute(
        term,
        random_estimate_with_lc_dd = estimate,
        random_q2.5_with_lc_dd = q2.5,
        random_q97.5_with_lc_dd = q97.5
      ) %>%
      inner_join(
        filtered$random %>%
          transmute(
            term,
            random_estimate_without_lc_dd = estimate,
            random_q2.5_without_lc_dd = q2.5,
            random_q97.5_without_lc_dd = q97.5
          ),
        by = "term"
      ) %>%
      mutate(random_delta_without_minus_with = random_estimate_without_lc_dd - random_estimate_with_lc_dd)
    
    write.csv(
      random_comparison,
      file.path(filtered_diag_dir, sprintf("%s_with_vs_without_lc_dd_random_effects.csv", model_name)),
      row.names = FALSE
    )
  }
  
  full_diag <- if (file.exists(full_diag_file)) {
    read.csv(full_diag_file, stringsAsFactors = FALSE) %>%
      filter(model == model_name) %>%
      mutate(dataset = "with_lc_dd")
  } else {
    data.frame()
  }
  
  write.csv(
    bind_rows(full_diag, filtered$diagnostics %>% mutate(dataset = "without_lc_dd")),
    file.path(filtered_diag_dir, sprintf("%s_with_vs_without_lc_dd_diagnostics.csv", model_name)),
    row.names = FALSE
  )
  
  input_summary_file <- file.path(
    "data", "processed", "traits", taxon, "cr_en_vu_nt",
    sprintf("%s_cr_en_vu_nt_model_input_summary.csv", taxon)
  )
  input_summary <- if (file.exists(input_summary_file)) {
    read.csv(input_summary_file, stringsAsFactors = FALSE)
  } else {
    data.frame(metric = character(), value = character())
  }
  
  full_data_file <- file.path("data", "processed", "traits", taxon, sprintf("%s_trait_model_final_data.csv", taxon))
  response_col <- if (model_name == "logging") "Threat_Logging" else "Threat_Plantation"
  response_prefix <- if (model_name == "logging") "n_logging" else "n_plantation"
  full_counts <- if (file.exists(full_data_file)) {
    table(as_binary_true(read.csv(full_data_file, stringsAsFactors = FALSE)[[response_col]]))
  } else {
    table()
  }
  count_value <- function(tab, nm) if (nm %in% names(tab)) as.character(unname(tab[[nm]])) else NA_character_
  
  sample_sizes <- data.frame(
    dataset = c("with_lc_dd", "without_lc_dd"),
    included_categories = c("CR;EN;VU;NT;LC;DD", metric_value(input_summary, "included_red_list_categories")),
    n_rows = c(
      first_metric_value(input_summary, c("n_original_rows", "n_original_species")),
      first_metric_value(input_summary, c("n_final_rows", "n_final_species"))
    ),
    n_threatened = c(count_value(full_counts, "TRUE"), metric_value(input_summary, paste0(response_prefix, "_threatened"))),
    n_not_threatened = c(count_value(full_counts, "FALSE"), metric_value(input_summary, paste0(response_prefix, "_not_threatened")))
  )
  write.csv(sample_sizes, file.path(filtered_diag_dir, sprintf("%s_with_vs_without_lc_dd_sample_sizes.csv", model_name)), row.names = FALSE)
  
  plot_order <- comparison %>% arrange(estimate_with_lc_dd) %>% pull(term)
  plot_data <- bind_rows(
    comparison %>% transmute(term, dataset = "With LC/DD", estimate = estimate_with_lc_dd, q2.5 = q2.5_with_lc_dd, q97.5 = q97.5_with_lc_dd),
    comparison %>% transmute(term, dataset = "Without LC/DD", estimate = estimate_without_lc_dd, q2.5 = q2.5_without_lc_dd, q97.5 = q97.5_without_lc_dd)
  ) %>%
    mutate(
      term = factor(term, levels = plot_order),
      dataset = factor(dataset, levels = c("With LC/DD", "Without LC/DD"))
    )
  
  p <- ggplot(plot_data, aes(x = estimate, y = term, color = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray45") +
    geom_errorbar(aes(xmin = q2.5, xmax = q97.5), orientation = "y", position = position_dodge(width = 0.65), width = 0, linewidth = 0.55) +
    geom_point(position = position_dodge(width = 0.65), size = 1.7) +
    scale_color_manual(values = c("With LC/DD" = "#6C6C6C", "Without LC/DD" = "#0072B2")) +
    labs(
      title = sprintf("%s %s Model: Effect of Excluding LC/DD Species", tools::toTitleCase(taxon), tools::toTitleCase(model_name)),
      x = sprintf("Effect on %s threat status (log-odds)", model_name),
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_line(color = "gray92", linetype = "dotted"),
      panel.grid.minor = element_blank()
    )
  
  plot_file <- file.path(figure_dir, sprintf("%s_%s_with_vs_without_lc_dd_coefficients.png", taxon, model_name))
  ggsave(plot_file, p, width = 10, height = max(5, 2 + length(unique(plot_data$term)) * 0.32), dpi = 300)
  
  cat("Wrote:\n")
  cat("  ", comparison_file, "\n", sep = "")
  cat("  ", plot_file, "\n", sep = "")
  invisible(comparison_file)
}

for (taxon in taxa_to_run) {
  for (model_name in models_to_run) {
    compare_one_model(taxon, model_name)
  }
}

cat("\nComparison script complete.\n")
