rm(list = ls())
gc()

# ==============================================================================
# POST-PROCESSING: Simple coefficient plot for HPC Fish Models
# Models: logging_combined.rds / plantation_combined.rds
#
# This script does not refit the models. It summarises the combined HPC fits and
# creates a coefficient forest plot for the fish trait models.
#
# Interpretation:
#   - Continuous predictors: slope on the log-odds scale
#   - Binary predictors: 0->1 contrast relative to absence
#   - Categorical predictors: contrast relative to the reference level
#   - Interaction term: coefficient for MaxLength_Log : log_range_size
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)

taxon_label <- "fish"
model_dir <- file.path("outputs", "traits", "models", taxon_label)
figure_dir <- file.path("outputs", "traits", "figures")
diagnostic_dir <- file.path("outputs", "traits", "diagnostics", taxon_label)

dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & MODELS
# ------------------------------------------------------------------------------

final_data <- read.csv("data/processed/traits/fish/fish_trait_model_final_data.csv") %>%
  mutate(DemersPelag_Cat = relevel(as.factor(DemersPelag_Cat), ref = "benthopelagic"))

model_paths <- c(
  logging = file.path(model_dir, "logging_combined.rds"),
  plantation = file.path(model_dir, "plantation_combined.rds")
)

missing_models <- model_paths[!file.exists(model_paths)]
if (length(missing_models) > 0) {
  stop(
    "Missing combined model file(s):\n",
    paste(missing_models, collapse = "\n"),
    "\nRun scripts/03_traits/99_model_diagnostics.R fish after the HPC chains finish."
  )
}

message("Loading combined fish models...")
model_logging <- readRDS(model_paths[["logging"]])
model_plantation <- readRDS(model_paths[["plantation"]])

# ------------------------------------------------------------------------------
# 2. WRITE SUMMARY TABLES
# ------------------------------------------------------------------------------

write_model_summary_tables <- function(model, model_name, taxon_label, n_obs) {
  overall_summary <- data.frame(
    taxon = taxon_label,
    model = model_name,
    family = model$family$family,
    link = model$family$link,
    n_obs = n_obs
  )

  fixed_summary <- as.data.frame(fixef(model, summary = TRUE))
  fixed_summary$term <- rownames(fixed_summary)
  rownames(fixed_summary) <- NULL
  fixed_summary <- fixed_summary %>%
    transmute(
      taxon = taxon_label,
      model = model_name,
      term = term,
      estimate = Estimate,
      std_error = Est.Error,
      q2.5 = Q2.5,
      q97.5 = Q97.5,
      direction = case_when(
        Q2.5 > 0 ~ "positive",
        Q97.5 < 0 ~ "negative",
        TRUE ~ "uncertain"
      )
    )

  random_summary <- as.data.frame(posterior_summary(model, pars = c("^sd_", "^cor_")))
  if (nrow(random_summary) > 0) {
    random_summary$term <- rownames(random_summary)
    rownames(random_summary) <- NULL
    random_summary <- random_summary %>%
      transmute(
        taxon = taxon_label,
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
        estimate = Estimate,
        std_error = Est.Error,
        q2.5 = Q2.5,
        q97.5 = Q97.5
      )
  }

  write.csv(
    overall_summary,
    file.path(diagnostic_dir, sprintf("%s_overall_fit_summary.csv", model_name)),
    row.names = FALSE
  )
  write.csv(
    fixed_summary,
    file.path(diagnostic_dir, sprintf("%s_fixed_effects_summary.csv", model_name)),
    row.names = FALSE
  )
  write.csv(
    random_summary,
    file.path(diagnostic_dir, sprintf("%s_random_effects_summary.csv", model_name)),
    row.names = FALSE
  )
}

write_model_summary_tables(model_logging, "logging", taxon_label, nrow(final_data))
write_model_summary_tables(model_plantation, "plantation", taxon_label, nrow(final_data))

# ------------------------------------------------------------------------------
# 3. EXTRACT POSTERIOR DRAWS OF FIXED EFFECTS
# ------------------------------------------------------------------------------

clean_demers_label <- function(term) {
  level <- sub("^b_DemersPelag_Cat", "", term)
  level <- gsub("_", " ", level)
  level <- gsub("([a-z])([A-Z])", "\\1 \\2", level)
  paste("Demersal/Pelagic:", level)
}

extract_fixed_effect_draws <- function(model, response_label) {
  draws_df <- as_draws_df(model)

  coef_map <- c(
    # Continuous / interaction
    "b_MaxLength_Log" = "Maximum Length",
    "b_log_range_size" = "Range Size",
    "b_Salinity_Breadth_sc" = "Salinity Breadth",
    "b_Trophic_Level_sc" = "Trophic Level",
    "b_PredPreyRatioMax_log_sc" = "Predator:Prey Ratio",
    "b_MaxLength_Log:log_range_size" = "Maximum Length x Range Size",

    # Binary
    "b_AirBreathing_Bin" = "Air Breathing"
  )

  demers_terms <- grep("^b_DemersPelag_Cat", names(draws_df), value = TRUE)
  if (length(demers_terms) > 0) {
    coef_map[demers_terms] <- vapply(demers_terms, clean_demers_label, character(1))
  }

  present_terms <- names(coef_map)[names(coef_map) %in% names(draws_df)]
  missing_terms <- setdiff(names(coef_map), present_terms)

  if (length(present_terms) == 0) {
    stop("No expected fixed-effect coefficients were found in the model draws.")
  }

  if (length(missing_terms) > 0) {
    warning(
      "Some expected fish coefficients were not found and will be skipped: ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
  }

  draws_df %>%
    select(all_of(present_terms)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "term",
      values_to = ".value"
    ) %>%
    mutate(
      Clean_Label = unname(coef_map[term]),
      Response = response_label
    )
}

message("Extracting posterior coefficient draws...")
draws_logging <- extract_fixed_effect_draws(model_logging, "ThreatLogging")
draws_plantation <- extract_fixed_effect_draws(model_plantation, "ThreatPlantation")

plot_data <- bind_rows(draws_logging, draws_plantation)
saveRDS(plot_data, file.path(model_dir, "hpc_coefficient_draws.rds"))

# ------------------------------------------------------------------------------
# 4. SUMMARISE & ORDER TERMS
# ------------------------------------------------------------------------------

plot_summary <- plot_data %>%
  group_by(Clean_Label, Response) %>%
  summarise(
    median = median(.value),
    q2.5 = quantile(.value, 0.025),
    q97.5 = quantile(.value, 0.975),
    .groups = "drop"
  )

plot_order <- plot_data %>%
  group_by(Clean_Label) %>%
  summarise(Sort_Value = mean(.value), .groups = "drop")

plot_data_sorted <- plot_data %>%
  left_join(plot_order, by = "Clean_Label")

plot_summary <- plot_summary %>%
  left_join(plot_order, by = "Clean_Label")

write.csv(
  plot_summary %>% arrange(Sort_Value, Response),
  file.path(diagnostic_dir, "coefficient_plot_summary.csv"),
  row.names = FALSE
)

reference_trait_rows <- 21
reference_plot_height <- 9
plot_fixed_height <- 2
plot_height_per_trait <- (reference_plot_height - plot_fixed_height) / reference_trait_rows
plot_height <- plot_fixed_height + n_distinct(plot_data_sorted$Clean_Label) * plot_height_per_trait

# ------------------------------------------------------------------------------
# 5. SIMPLE COEFFICIENT FOREST PLOT
# ------------------------------------------------------------------------------

plot_mcmc <- ggplot(
  plot_data_sorted,
  aes(
    x = .value,
    y = reorder(Clean_Label, Sort_Value),
    color = Response
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  stat_pointinterval(
    position = position_dodge(width = 0.6),
    point_interval = median_qi,
    .width = 0.95,
    size = 1.1
  ) +
  scale_color_manual(
    values = c(
      "ThreatLogging" = "#d9be00",
      "ThreatPlantation" = "#0088d9"
    ),
    labels = c("Logging", "Plantations")
  ) +
  labs(
    title = "Trait Effects on Threat Risk in Freshwater Fish",
    subtitle = "Posterior median coefficients with 95% credible intervals",
    x = "Effect on Risk (Log-Odds)",
    y = NULL,
    color = "Threat Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    panel.grid.major.y = element_line(color = "gray95", linetype = "dotted")
  )

print(plot_mcmc)

ggsave(
  file.path(figure_dir, "fish_coefficients.png"),
  plot_mcmc,
  width = 10,
  height = plot_height
)

message("Saved fish model plot to ", file.path(figure_dir, "fish_coefficients.png"))
