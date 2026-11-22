rm(list = ls())
gc()

# =============================================================================
# COMPARE BIRD LOGGING MODEL RESULTS WITH AND WITHOUT LC/DD SPECIES
#
# "With LC/DD" uses the existing all-category bird logging model summaries.
# "Without LC/DD" uses the CR/EN/VU/NT-only logging fit from:
#   outputs/traits/models/bird_cr_en_vu_nt/logging_combined.rds
#
# This script does not refit models.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(dplyr)
library(ggplot2)

taxon_original <- "bird"
taxon_filtered <- "bird_cr_en_vu_nt"
model_name <- "logging"

original_diag_dir <- file.path("outputs", "traits", "diagnostics", taxon_original)
filtered_diag_dir <- file.path("outputs", "traits", "diagnostics", taxon_filtered)
filtered_model_dir <- file.path("outputs", "traits", "models", taxon_filtered)
figure_dir <- file.path("outputs", "traits", "figures")

dir.create(filtered_diag_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

original_fixed_file <- file.path(
  original_diag_dir,
  "logging_fixed_effects_summary.csv"
)
original_random_file <- file.path(
  original_diag_dir,
  "logging_random_effects_summary.csv"
)
original_diag_file <- file.path(
  original_diag_dir,
  "diagnostic_summary.csv"
)

filtered_model_file <- file.path(filtered_model_dir, "logging_combined.rds")
if (!file.exists(filtered_model_file)) {
  filtered_model_file <- file.path(
    filtered_model_dir,
    "bird_cr_en_vu_nt_logging_pc_4chains.rds"
  )
}

if (!file.exists(original_fixed_file)) {
  stop(
    "Missing original bird logging fixed-effect summary: ",
    original_fixed_file,
    call. = FALSE
  )
}
if (!file.exists(filtered_model_file)) {
  stop(
    "Missing CR/EN/VU/NT logging model: ",
    filtered_model_file,
    "\nRun scripts/03_traits/bird/fit_pc_cr_en_vu_nt.R logging first.",
    call. = FALSE
  )
}

input_summary_file <- file.path(
  "data",
  "processed",
  "traits",
  "bird",
  "cr_en_vu_nt",
  "bird_cr_en_vu_nt_model_input_summary.csv"
)

cat("\n========================================\n")
cat("Comparing bird logging model sensitivity\n")
cat("========================================\n")

cat("Loading original fixed-effect summary...\n")
original_fixed <- read.csv(original_fixed_file, stringsAsFactors = FALSE) %>%
  transmute(
    term,
    estimate_with_lc_dd = estimate,
    std_error_with_lc_dd = std_error,
    q2.5_with_lc_dd = q2.5,
    q97.5_with_lc_dd = q97.5,
    direction_with_lc_dd = direction
  )

cat("Loading CR/EN/VU/NT logging model...\n")
filtered_fit <- readRDS(filtered_model_file)

filtered_n_obs <- tryCatch(
  nrow(filtered_fit$data),
  error = function(e) NA_integer_
)

stan_fit <- filtered_fit[["fit"]]
stan_sim <- attr(stan_fit, "sim")
stan_samples <- stan_sim[["samples"]]
warmup_by_chain <- stan_sim[["warmup2"]]

draws_by_chain <- function(param) {
  lapply(seq_along(stan_samples), function(chain_id) {
    chain_draws <- stan_samples[[chain_id]][[param]]
    if (is.null(chain_draws)) {
      stop("Parameter not found in saved Stan samples: ", param, call. = FALSE)
    }
    warmup <- warmup_by_chain[[chain_id]]
    if (length(chain_draws) <= warmup) {
      stop("No post-warmup draws found for parameter: ", param, call. = FALSE)
    }
    chain_draws[(warmup + 1):length(chain_draws)]
  })
}

summarise_param <- function(param) {
  chains <- draws_by_chain(param)
  draws <- unlist(chains, use.names = FALSE)
  data.frame(
    param = param,
    estimate = mean(draws),
    std_error = sd(draws),
    q2.5 = unname(quantile(draws, 0.025)),
    q97.5 = unname(quantile(draws, 0.975))
  )
}

split_rhat <- function(param) {
  chains <- draws_by_chain(param)
  n <- min(vapply(chains, length, integer(1)))
  chains <- lapply(chains, function(x) x[seq_len(n)])
  m <- length(chains)
  if (m < 2 || n < 2) return(NA_real_)
  
  chain_means <- vapply(chains, mean, numeric(1))
  chain_vars <- vapply(chains, var, numeric(1))
  w <- mean(chain_vars)
  b <- n * var(chain_means)
  if (!is.finite(w) || w == 0) return(NA_real_)
  var_plus <- ((n - 1) / n) * w + b / n
  sqrt(var_plus / w)
}

param_to_term <- function(param) {
  if (param == "b_Intercept") return("Intercept")
  sub("^b_", "", param)
}

fixed_params <- grep("^b_", names(stan_samples[[1]]), value = TRUE)
random_sd_params <- grep("^sd_", names(stan_samples[[1]]), value = TRUE)

filtered_fixed <- bind_rows(lapply(fixed_params, summarise_param)) %>%
  mutate(term = vapply(param, param_to_term, character(1))) %>%
  select(-param)
filtered_fixed <- filtered_fixed %>%
  transmute(
    taxon = taxon_filtered,
    model = model_name,
    term,
    estimate = estimate,
    std_error = std_error,
    q2.5 = q2.5,
    q97.5 = q97.5,
    direction = case_when(
      q2.5 > 0 ~ "positive",
      q97.5 < 0 ~ "negative",
      TRUE ~ "uncertain"
    )
  )

write.csv(
  filtered_fixed,
  file.path(filtered_diag_dir, "logging_fixed_effects_summary.csv"),
  row.names = FALSE
)

if (length(random_sd_params) > 0) {
  filtered_random <- bind_rows(lapply(random_sd_params, summarise_param)) %>%
    rename(term = param)
  filtered_random <- filtered_random %>%
    transmute(
      taxon = taxon_filtered,
      model = model_name,
      term = term,
      group = case_when(
        grepl("phylo_name", term) ~ "phylo_name",
        grepl("spatial_name", term) ~ "spatial_name",
        TRUE ~ NA_character_
      ),
      effect_type = case_when(
        grepl("^sd_", term) ~ "sd",
        grepl("^cor_", term) ~ "correlation",
        TRUE ~ "other"
      ),
      estimate = estimate,
      std_error = std_error,
      q2.5 = q2.5,
      q97.5 = q97.5
    )
} else {
  filtered_random <- data.frame()
}

write.csv(
  filtered_random,
  file.path(filtered_diag_dir, "logging_random_effects_summary.csv"),
  row.names = FALSE
)

filtered_overall <- data.frame(
  taxon = taxon_filtered,
  model = model_name,
  family = filtered_fit$family$family,
  link = filtered_fit$family$link,
  n_obs = filtered_n_obs
)

write.csv(
  filtered_overall,
  file.path(filtered_diag_dir, "logging_overall_fit_summary.csv"),
  row.names = FALSE
)

cat("Computing CR/EN/VU/NT model diagnostics...\n")
diag_filtered <- tryCatch({
  rhat_vals <- vapply(c(fixed_params, random_sd_params), split_rhat, numeric(1))
  data.frame(
    taxon = taxon_filtered,
    model = model_name,
    max_rhat = max(rhat_vals, na.rm = TRUE),
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    n_divergences = NA_integer_,
    passed = NA
  )
}, error = function(e) {
  warning("Could not compute filtered diagnostics: ", e$message, call. = FALSE)
  data.frame(
    taxon = taxon_filtered,
    model = model_name,
    max_rhat = NA_real_,
    min_bulk_ess = NA_real_,
    min_tail_ess = NA_real_,
    n_divergences = NA_integer_,
    passed = NA
  )
})

write.csv(
  diag_filtered,
  file.path(filtered_diag_dir, "diagnostic_summary.csv"),
  row.names = FALSE
)

comparison <- original_fixed %>%
  inner_join(
    filtered_fixed %>%
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

comparison_file <- file.path(
  filtered_diag_dir,
  "logging_with_vs_without_lc_dd_fixed_effects.csv"
)
write.csv(comparison, comparison_file, row.names = FALSE)

if (file.exists(original_random_file) && nrow(filtered_random) > 0) {
  original_random <- read.csv(original_random_file, stringsAsFactors = FALSE) %>%
    transmute(
      term,
      random_estimate_with_lc_dd = estimate,
      random_q2.5_with_lc_dd = q2.5,
      random_q97.5_with_lc_dd = q97.5
    )
  
  random_comparison <- original_random %>%
    inner_join(
      filtered_random %>%
        transmute(
          term,
          random_estimate_without_lc_dd = estimate,
          random_q2.5_without_lc_dd = q2.5,
          random_q97.5_without_lc_dd = q97.5
        ),
      by = "term"
    ) %>%
    mutate(
      random_delta_without_minus_with =
        random_estimate_without_lc_dd - random_estimate_with_lc_dd
    )
  
  write.csv(
    random_comparison,
    file.path(filtered_diag_dir, "logging_with_vs_without_lc_dd_random_effects.csv"),
    row.names = FALSE
  )
}

original_diag <- if (file.exists(original_diag_file)) {
  read.csv(original_diag_file, stringsAsFactors = FALSE) %>%
    filter(model == model_name) %>%
    mutate(dataset = "with_lc_dd")
} else {
  data.frame()
}

diagnostic_comparison <- bind_rows(
  original_diag,
  diag_filtered %>% mutate(dataset = "without_lc_dd")
)

write.csv(
  diagnostic_comparison,
  file.path(filtered_diag_dir, "logging_with_vs_without_lc_dd_diagnostics.csv"),
  row.names = FALSE
)

input_summary <- if (file.exists(input_summary_file)) {
  read.csv(input_summary_file, stringsAsFactors = FALSE)
} else {
  data.frame(metric = character(), value = character())
}

original_data_file <- file.path(
  "data",
  "processed",
  "traits",
  "bird",
  "bird_trait_model_final_data.csv"
)
original_logging_counts <- if (file.exists(original_data_file)) {
  original_data <- read.csv(original_data_file, stringsAsFactors = FALSE)
  table(original_data$Threat_Logging)
} else {
  table()
}

get_count <- function(count_table, value_name) {
  if (!value_name %in% names(count_table)) return(NA_character_)
  as.character(unname(count_table[[value_name]]))
}

get_metric <- function(metric_name) {
  value <- input_summary$value[input_summary$metric == metric_name]
  if (length(value) == 0) NA_character_ else value[[1]]
}

sample_comparison <- data.frame(
  dataset = c("with_lc_dd", "without_lc_dd"),
  included_categories = c(
    "CR;EN;VU;NT;LC;DD",
    get_metric("included_red_list_categories")
  ),
  n_species = c(
    get_metric("n_original_species"),
    get_metric("n_final_species")
  ),
  n_logging_threatened = c(
    get_count(original_logging_counts, "TRUE"),
    get_metric("n_logging_threatened")
  ),
  n_logging_not_threatened = c(
    get_count(original_logging_counts, "FALSE"),
    get_metric("n_logging_not_threatened")
  )
)

write.csv(
  sample_comparison,
  file.path(filtered_diag_dir, "logging_with_vs_without_lc_dd_sample_sizes.csv"),
  row.names = FALSE
)

plot_labels <- c(
  "Intercept" = "Intercept",
  "Body_Mass_log" = "Body mass",
  "log_range_size" = "Range size",
  "Habitat_Breadth_log" = "Habitat breadth",
  "Clutch_Size_log" = "Clutch size",
  "HWI_log" = "Hand-wing index",
  "Forest" = "Forest",
  "Savanna" = "Savanna",
  "Shrubland" = "Shrubland",
  "Grassland" = "Grassland",
  "Wetlands" = "Wetlands",
  "Diet_Invertebrate" = "Diet: invertebrate",
  "Diet_Fruit" = "Diet: fruit",
  "Diet_Nectar" = "Diet: nectar",
  "Diet_Seed" = "Diet: seed",
  "Diet_Vertebrate" = "Diet: vertebrate",
  "Diet_Fish" = "Diet: fish",
  "Diet_Scavenger" = "Diet: scavenger",
  "Diet_PlantMaterial" = "Diet: plant material",
  "Body_Mass_log:log_range_size" = "Body mass x range size"
)

plot_data <- bind_rows(
  comparison %>%
    transmute(
      term,
      label = ifelse(term %in% names(plot_labels), plot_labels[term], term),
      dataset = "With LC/DD",
      estimate = estimate_with_lc_dd,
      q2.5 = q2.5_with_lc_dd,
      q97.5 = q97.5_with_lc_dd
    ),
  comparison %>%
    transmute(
      term,
      label = ifelse(term %in% names(plot_labels), plot_labels[term], term),
      dataset = "Without LC/DD",
      estimate = estimate_without_lc_dd,
      q2.5 = q2.5_without_lc_dd,
      q97.5 = q97.5_without_lc_dd
    )
) %>%
  mutate(
    label = factor(
      label,
      levels = comparison %>%
        arrange(estimate_with_lc_dd) %>%
        pull(term) %>%
        { ifelse(. %in% names(plot_labels), plot_labels[.], .) }
    ),
    dataset = factor(dataset, levels = c("With LC/DD", "Without LC/DD"))
  )

comparison_plot <- ggplot(
  plot_data,
  aes(x = estimate, y = label, color = dataset)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray45") +
  geom_errorbarh(
    aes(xmin = q2.5, xmax = q97.5),
    position = position_dodge(width = 0.65),
    height = 0,
    linewidth = 0.6
  ) +
  geom_point(position = position_dodge(width = 0.65), size = 1.8) +
  scale_color_manual(
    values = c("With LC/DD" = "#6C6C6C", "Without LC/DD" = "#0072B2")
  ) +
  labs(
    title = "Bird Logging Model: Effect of Excluding LC/DD Species",
    subtitle = "Posterior median coefficients with 95% credible intervals",
    x = "Effect on logging threat status (log-odds)",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray92", linetype = "dotted"),
    panel.grid.minor = element_blank()
  )

plot_file <- file.path(
  figure_dir,
  "bird_logging_with_vs_without_lc_dd_coefficients.png"
)
ggsave(plot_file, comparison_plot, width = 10, height = 8, dpi = 300)

cat("\nWrote comparison outputs:\n")
cat("  ", comparison_file, "\n", sep = "")
cat("  ", plot_file, "\n", sep = "")
cat("  ", file.path(filtered_diag_dir, "logging_fixed_effects_summary.csv"), "\n", sep = "")
cat("  ", file.path(filtered_diag_dir, "diagnostic_summary.csv"), "\n", sep = "")

cat("\nLargest fixed-effect shifts:\n")
print(
  comparison %>%
    select(
      term,
      estimate_with_lc_dd,
      estimate_without_lc_dd,
      delta_without_minus_with,
      direction_with_lc_dd,
      direction_without_lc_dd,
      sensitivity_flag
    ) %>%
    head(10)
)
