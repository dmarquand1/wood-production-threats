rm(list = ls())
gc()

# ==============================================================================
# POST-PROCESSING: Simple coefficient plot for HPC Mammal Models
# Models: logging_combined.rds / plantation_combined.rds
#
# This version drops emmeans/emtrends and simply plots posterior summaries of
# the fitted model coefficients (95% credible intervals), with all traits shown
# together in a single panel.
#
# Interpretation:
#   - Continuous predictors: slope on the log-odds scale
#   - Binary predictors: 0->1 contrast relative to absence
#   - Categorical predictors: contrast relative to the reference level
#   - Interaction term: coefficient for Body_Mass_Log : log_range_size
# ==============================================================================

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)

model_dir  <- file.path("outputs", "traits", "models", "mammal")
figure_dir <- file.path("outputs", "traits", "figures")
diagnostic_dir <- file.path("outputs", "traits", "diagnostics", "mammal")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & MODELS
# ------------------------------------------------------------------------------

final_data <- read.csv("data/processed/traits/mammal/mammal_trait_model_final_data.csv") %>%
  mutate(activity_cycle = relevel(as.factor(activity_cycle), ref = "Nocturnal"))

model_logging    <- readRDS(file.path(model_dir, "logging_combined.rds"))
model_plantation <- readRDS(file.path(model_dir, "plantation_combined.rds"))

write_model_summary_tables <- function(model, model_name, taxon_label) {
  overall_summary <- data.frame(
    taxon = taxon_label,
    model = model_name,
    family = model$family$family,
    link = model$family$link,
    n_obs = nobs(model)
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
  
  write.csv(overall_summary, file.path(diagnostic_dir, sprintf("%s_overall_fit_summary.csv", model_name)), row.names = FALSE)
  write.csv(fixed_summary, file.path(diagnostic_dir, sprintf("%s_fixed_effects_summary.csv", model_name)), row.names = FALSE)
  write.csv(random_summary, file.path(diagnostic_dir, sprintf("%s_random_effects_summary.csv", model_name)), row.names = FALSE)
}

write_model_summary_tables(model_logging, "logging", "mammal")
write_model_summary_tables(model_plantation, "plantation", "mammal")

# ------------------------------------------------------------------------------
# 2. FUNCTION TO EXTRACT POSTERIOR DRAWS OF FIXED EFFECTS
# ------------------------------------------------------------------------------

extract_fixed_effect_draws <- function(model, response_label) {
  
  draws_df <- as_draws_df(model)
  
  coef_map <- c(
    # Continuous / interaction
    "b_Body_Mass_Log"                = "Body Mass",
    "b_log_range_size"               = "Range Size",
    "b_Habitat_Breadth"              = "Habitat Breadth",
    "b_Litter_Size_Log"              = "Litter Size",
    "b_Body_Mass_Log:log_range_size" = "Body Mass × Range Size",
    
    # Categorical: activity cycle (relative to Nocturnal)
    "b_activity_cycleCathemeral_Crepuscular"     = "Cathemeral/Crepuscular",
    "b_activity_cycleDiurnal"        = "Diurnal",
    
    # Habitat binaries
    "b_Forest"                       = "Forest",
    "b_Savanna"                      = "Savanna",
    "b_Shrubland"                    = "Shrubland",
    "b_Grassland"                    = "Grassland",
    "b_Wetlands"                     = "Wetlands",
    
    # Diet binaries
    "b_binary_det_inv"               = "Invertebrates",
    "b_binary_det_vend"              = "Endotherm Vertebrates",
    "b_binary_det_vect"              = "Ectotherm Vertebrates",
    "b_binary_det_vunk"              = "Unknown Vertebrates",
    "b_binary_det_scav"              = "Scavenger",
    "b_binary_det_fruit"             = "Fruit",
    "b_binary_det_nect"              = "Nectar",
    "b_binary_det_seed"              = "Seed",
    "b_binary_det_plantother"        = "Other Plant Material"
  )
  
  present_terms <- names(coef_map)[names(coef_map) %in% names(draws_df)]
  
  if (length(present_terms) == 0) {
    stop("No expected fixed-effect coefficients were found in the model draws.")
  }
  
  plot_draws <- draws_df %>%
    select(all_of(present_terms)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "term",
      values_to = ".value"
    ) %>%
    mutate(
      Clean_Label = unname(coef_map[term]),
      Response    = response_label
    )
  
  return(plot_draws)
}

# ------------------------------------------------------------------------------
# 3. EXTRACT DRAWS
# ------------------------------------------------------------------------------

draws_logging    <- extract_fixed_effect_draws(model_logging,    "ThreatLogging")
draws_plantation <- extract_fixed_effect_draws(model_plantation, "ThreatPlantation")

plot_data <- bind_rows(draws_logging, draws_plantation)

saveRDS(plot_data, file.path(model_dir, "hpc_coefficient_draws.rds"))

# ------------------------------------------------------------------------------
# 4. ORDER TERMS
# ------------------------------------------------------------------------------

plot_order <- plot_data %>%
  group_by(Clean_Label) %>%
  summarise(Sort_Value = mean(.value), .groups = "drop")

plot_data_sorted <- plot_data %>%
  left_join(plot_order, by = "Clean_Label")

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
    position       = position_dodge(width = 0.6),
    point_interval = median_qi,
    .width         = 0.95,
    size           = 1.1
  ) +
  scale_color_manual(
    values = c(
      "ThreatLogging"    = "#d9be00",
      "ThreatPlantation" = "#0088d9"
    ),
    labels = c("Logging", "Plantations")
  ) +
  labs(
    title    = "Trait Effects on Threat Risk in Mammals",
    subtitle = "Posterior median coefficients with 95% credible intervals",
    x        = "Effect on Risk (Log-Odds)",
    y        = NULL,
    color    = "Threat Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position    = "top",
    axis.title.x       = element_text(size = 15),
    axis.text.x        = element_text(size = 13),
    axis.text.y        = element_text(size = 13),
    panel.grid.major.y = element_line(color = "gray95", linetype = "dotted")
  )

print(plot_mcmc)

ggsave(
  file.path(figure_dir, "mammal_coefficients.png"),
  plot_mcmc,
  width = 10,
  height = plot_height
)
