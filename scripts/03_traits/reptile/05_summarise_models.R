rm(list = ls())
gc()

# ==============================================================================
# POST-PROCESSING: Simple coefficient plot for HPC Reptile Models
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
#   - Interaction term: coefficient for Max_Body_Mass_log : log_range_size
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)

model_dir  <- file.path("outputs", "traits", "models", "reptile")
figure_dir <- file.path("outputs", "traits", "figures")
diagnostic_dir <- file.path("outputs", "traits", "diagnostics", "reptile")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & MODELS
# ------------------------------------------------------------------------------

final_data <- read.csv("data/processed/traits/reptile/reptile_trait_model_final_data.csv") %>%
  mutate(
    Active.time       = relevel(as.factor(Active.time), ref = "Diurnal"),
    Reproductive.mode = relevel(as.factor(Reproductive.mode), ref = "oviparous"),
    Fangs             = relevel(as.factor(Fangs), ref = "non-fanged")
  )

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

write_model_summary_tables(model_logging, "logging", "reptile")
write_model_summary_tables(model_plantation, "plantation", "reptile")

# ------------------------------------------------------------------------------
# 2. FUNCTION TO EXTRACT POSTERIOR DRAWS OF FIXED EFFECTS
# ------------------------------------------------------------------------------

extract_fixed_effect_draws <- function(model, response_label) {
  
  draws_df <- as_draws_df(model)
  
  coef_map <- c(
    # Continuous / interaction
    "b_Max_Body_Mass_log"                = "Body Mass",
    "b_log_range_size"                   = "Range Size",
    "b_habitat_breadth_log"              = "Habitat Breadth",
    "b_Max_Body_Mass_log:log_range_size" = "Body Mass × Range Size",
    
    # Categorical: activity time (relative to Diurnal)
    "b_Active.timeCrepuscular_Cathemeral" = "Crepuscular/Cathemeral",
    "b_Active.timeCrepuscular_Cathermal" = "Crepuscular/Cathemeral",
    "b_Active.timeDiurnal"               = "Diurnal",
    "b_Active.timeNocturnal"             = "Nocturnal",
    
    # Categorical: reproductive mode (relative to oviparous)
    "b_Reproductive.modeviviparous_ovoviviparous"      = "Viviparous/Ovoviviparous",
    
    # Categorical: fangs (relative to non-fanged)
    "b_FangsfrontMfanged"                = "Front-fanged",
    "b_FangsrearMfanged"                 = "Rear-fanged",
    "b_FangsnonMfanged"                  = "Non-fanged",
    
    # Habitat binaries
    "b_Forest"                           = "Forest",
    "b_Savanna"                          = "Savanna",
    "b_Shrubland"                        = "Shrubland",
    "b_Grassland"                        = "Grassland",
    "b_Wetlands"                         = "Wetlands"
  )
  
  present_terms <- names(coef_map)[names(coef_map) %in% names(draws_df)]
  missing_terms <- setdiff(names(coef_map), present_terms)
  
  if (length(present_terms) == 0) {
    stop("No expected fixed-effect coefficients were found in the model draws.")
  }

  if (length(missing_terms) > 0) {
    warning(
      "Some expected reptile coefficients were not found and will be skipped: ",
      paste(missing_terms, collapse = ", "),
      call. = FALSE
    )
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

harmonise_reptile_factor_draws <- function(plot_data) {
  factor_labels <- c(
    "Crepuscular/Cathemeral", "Diurnal", "Nocturnal",
    "Front-fanged", "Non-fanged", "Rear-fanged"
  )
  
  indexed <- plot_data %>%
    group_by(Response, Clean_Label) %>%
    mutate(.draw = row_number()) %>%
    ungroup()
  
  get_draws <- function(response_label, clean_label) {
    indexed %>%
      filter(Response == response_label, Clean_Label == clean_label) %>%
      select(.draw, .value)
  }
  
  make_single_draw <- function(draws, response_label, clean_label, term_name) {
    draws %>%
      transmute(
        term = term_name,
        .value = .value,
        Clean_Label = clean_label,
        Response = response_label
      )
  }
  
  logging_activity <- indexed %>%
    filter(
      Response == "ThreatLogging",
      Clean_Label %in% c("Crepuscular/Cathemeral", "Nocturnal")
    ) %>%
    select(term, .value, Clean_Label, Response)
  
  logging_fangs <- indexed %>%
    filter(
      Response == "ThreatLogging",
      Clean_Label %in% c("Front-fanged", "Rear-fanged")
    ) %>%
    select(term, .value, Clean_Label, Response)
  
  plantation_diurnal <- get_draws("ThreatPlantation", "Diurnal")
  plantation_nocturnal <- get_draws("ThreatPlantation", "Nocturnal")
  plantation_non_fanged <- get_draws("ThreatPlantation", "Non-fanged")
  plantation_rear_fanged <- get_draws("ThreatPlantation", "Rear-fanged")
  
  plantation_activity <- bind_rows(
    make_single_draw(
      plantation_diurnal %>% mutate(.value = -.value),
      "ThreatPlantation",
      "Crepuscular/Cathemeral",
      "Active.timeCrepuscular_Cathemeral_common_ref_Diurnal"
    ),
    plantation_nocturnal %>%
      inner_join(plantation_diurnal, by = ".draw", suffix = c("_nocturnal", "_diurnal")) %>%
      transmute(.draw, .value = .value_nocturnal - .value_diurnal) %>%
      make_single_draw(
        "ThreatPlantation",
        "Nocturnal",
        "Active.timeNocturnal_common_ref_Diurnal"
      )
  )
  
  plantation_fangs <- bind_rows(
    make_single_draw(
      plantation_non_fanged %>% mutate(.value = -.value),
      "ThreatPlantation",
      "Front-fanged",
      "FangsfrontMfanged_common_ref_nonMfanged"
    ),
    plantation_rear_fanged %>%
      inner_join(plantation_non_fanged, by = ".draw", suffix = c("_rear", "_non")) %>%
      transmute(.draw, .value = .value_rear - .value_non) %>%
      make_single_draw(
        "ThreatPlantation",
        "Rear-fanged",
        "FangsrearMfanged_common_ref_nonMfanged"
      )
  )
  
  plot_data %>%
    filter(!Clean_Label %in% factor_labels) %>%
    bind_rows(logging_activity, plantation_activity, logging_fangs, plantation_fangs)
}

# ------------------------------------------------------------------------------
# 3. EXTRACT DRAWS
# ------------------------------------------------------------------------------

draws_logging    <- extract_fixed_effect_draws(model_logging,    "ThreatLogging")
draws_plantation <- extract_fixed_effect_draws(model_plantation, "ThreatPlantation")

plot_data <- bind_rows(draws_logging, draws_plantation)
plot_data <- harmonise_reptile_factor_draws(plot_data)

saveRDS(plot_data, file.path(model_dir, "hpc_coefficient_draws.rds"))

# ------------------------------------------------------------------------------
# 4. ORDER TERMS
# ------------------------------------------------------------------------------

plot_order <- plot_data %>%
  group_by(Clean_Label) %>%
  summarise(Sort_Value = mean(.value), .groups = "drop")

plot_data_sorted <- plot_data %>%
  left_join(plot_order, by = "Clean_Label")

plot_summary <- plot_data_sorted %>%
  group_by(Clean_Label, Response, Sort_Value) %>%
  summarise(
    estimate = median(.value),
    lower = quantile(.value, 0.025),
    upper = quantile(.value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    Clean_Label = factor(
      Clean_Label,
      levels = plot_order$Clean_Label[order(plot_order$Sort_Value)]
    )
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
  plot_summary,
  aes(
    y = Clean_Label,
    color = Response
  )
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  geom_segment(
    aes(x = lower, xend = upper, yend = Clean_Label),
    position = position_dodge(width = 0.6),
    linewidth = 1.1,
    lineend = "round"
  ) +
  geom_point(
    aes(x = estimate),
    position = position_dodge(width = 0.6),
    size = 2.4
  ) +
  scale_color_manual(
    values = c(
      "ThreatLogging"    = "#d9be00",
      "ThreatPlantation" = "#0088d9"
    ),
    labels = c("Logging", "Plantations")
  ) +
  labs(
    title    = "Trait Effects on Threat Risk in Reptiles",
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
  file.path(figure_dir, "reptile_coefficients.png"),
  plot_mcmc,
  width = 10,
  height = plot_height
)
