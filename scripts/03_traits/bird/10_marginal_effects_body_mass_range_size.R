rm(list = ls())
gc()

# ==============================================================================
# MARGINAL EFFECTS: Body Mass x Range Size interaction on threat probability
#
# Models: logging_combined.rds / plantation_combined.rds (bird HPC fits)
#
# For each threat model, predicts P(threatened) across a grid of log_range_size
# at three fixed Body_Mass_log values (10th / 50th / 90th percentile of the
# observed bird data), holding all other predictors at their reference level
# (continuous traits at sample mean, binary traits at 0).
#
# Predictions are computed directly from the fixed-effect posterior draws
# rather than via brms::fitted()/posterior_epred(), because those functions
# are extremely slow on these models: they reconstruct the phylogenetic and
# spatial covariance (Cholesky) structures for the gr(..., cov = ...) terms
# even when re_formula = NA drops them from the linear predictor. Since every
# dropped covariate here is either a binary held at 0 or continuous held at
# its mean, the population-level linear predictor can be computed directly:
#
#   eta = b_Intercept
#       + b_Body_Mass_log * bm + b_log_range_size * rs
#       + b_Body_Mass_log:log_range_size * bm * rs
#       + b_Habitat_Breadth_log * mean(Habitat_Breadth_log)
#       + b_Clutch_Size_log * mean(Clutch_Size_log)
#       + b_HWI_log * mean(HWI_log)
# ==============================================================================

project_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib, winslash = "/", mustWork = FALSE), .libPaths()))
}

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

model_dir  <- file.path("outputs", "traits", "models", "bird")
figure_dir <- file.path("outputs", "traits", "figures")
diagnostic_dir <- file.path("outputs", "traits", "diagnostics", "bird")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. LOAD DATA & MODELS
# ------------------------------------------------------------------------------

final_data <- read.csv("data/processed/traits/bird/bird_trait_model_final_data.csv")

# ------------------------------------------------------------------------------
# 2. BUILD PREDICTION GRID
# ------------------------------------------------------------------------------

body_mass_q <- quantile(final_data$Body_Mass_log, probs = c(0.10, 0.50, 0.90), na.rm = TRUE)
body_mass_levels <- tibble(
  percentile = factor(c("10th", "50th", "90th"), levels = c("10th", "50th", "90th")),
  Body_Mass_log = as.numeric(body_mass_q)
)

range_size_seq <- seq(
  min(final_data$log_range_size, na.rm = TRUE),
  max(final_data$log_range_size, na.rm = TRUE),
  length.out = 100
)

reference_row <- final_data %>%
  summarise(
    Habitat_Breadth_log = mean(Habitat_Breadth_log, na.rm = TRUE),
    Clutch_Size_log     = mean(Clutch_Size_log, na.rm = TRUE),
    HWI_log             = mean(HWI_log, na.rm = TRUE)
  )

newdata_grid <- expand_grid(
  body_mass_levels,
  log_range_size = range_size_seq
) %>%
  bind_cols(reference_row[rep(1, nrow(.)), ])

# ------------------------------------------------------------------------------
# 3. COMPUTE MARGINAL PREDICTIONS FROM POSTERIOR DRAWS
# ------------------------------------------------------------------------------

predict_marginal <- function(model, newdata, response_label) {
  draws <- as_draws_df(model)

  required_terms <- c(
    "b_Intercept", "b_Body_Mass_log", "b_log_range_size",
    "b_Body_Mass_log:log_range_size",
    "b_Habitat_Breadth_log", "b_Clutch_Size_log", "b_HWI_log"
  )
  missing_terms <- setdiff(required_terms, names(draws))
  if (length(missing_terms) > 0) {
    stop("Missing expected posterior draw columns: ", paste(missing_terms, collapse = ", "))
  }

  draws <- draws %>% select(all_of(required_terms))

  grid_results <- lapply(seq_len(nrow(newdata)), function(i) {
    row <- newdata[i, ]
    eta <- draws[["b_Intercept"]] +
      draws[["b_Body_Mass_log"]] * row$Body_Mass_log +
      draws[["b_log_range_size"]] * row$log_range_size +
      draws[["b_Body_Mass_log:log_range_size"]] * row$Body_Mass_log * row$log_range_size +
      draws[["b_Habitat_Breadth_log"]] * row$Habitat_Breadth_log +
      draws[["b_Clutch_Size_log"]] * row$Clutch_Size_log +
      draws[["b_HWI_log"]] * row$HWI_log

    p <- plogis(eta)
    tibble(
      estimate = median(p),
      lower    = quantile(p, 0.025),
      upper    = quantile(p, 0.975)
    )
  })

  newdata %>%
    bind_cols(bind_rows(grid_results)) %>%
    mutate(Response = response_label)
}

message("Loading logging model (large file, may take a while)...")
model_logging <- readRDS(file.path(model_dir, "logging_combined.rds"))
message("Computing marginal predictions for logging model...")
pred_logging <- predict_marginal(model_logging, newdata_grid, "ThreatLogging")
rm(model_logging)
gc()

message("Loading plantation model (large file, may take a while)...")
model_plantation <- readRDS(file.path(model_dir, "plantation_combined.rds"))
message("Computing marginal predictions for plantation model...")
pred_plantation <- predict_marginal(model_plantation, newdata_grid, "ThreatPlantation")
rm(model_plantation)
gc()

plot_data <- bind_rows(pred_logging, pred_plantation) %>%
  mutate(
    Response = factor(Response, levels = c("ThreatLogging", "ThreatPlantation"))
  )

write.csv(
  plot_data %>%
    select(Response, percentile, Body_Mass_log, log_range_size, estimate, lower, upper),
  file.path(diagnostic_dir, "marginal_effects_body_mass_range_size.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 4. PLOT
# ------------------------------------------------------------------------------

threat_labels <- c(
  "ThreatLogging"    = "Logging",
  "ThreatPlantation" = "Plantations"
)

percentile_colors <- c(
  "10th" = "#1287c9",
  "50th" = "#4a4a4a",
  "90th" = "#c7481c"
)

percentile_labels <- sprintf(
  "%s pct (Body Mass Log = %.2f)",
  body_mass_levels$percentile,
  body_mass_levels$Body_Mass_log
)
names(percentile_labels) <- body_mass_levels$percentile

marginal_plot <- ggplot(
  plot_data,
  aes(x = log_range_size, y = estimate, color = percentile, fill = percentile)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ Response, nrow = 1, labeller = labeller(Response = threat_labels)) +
  scale_color_manual(values = percentile_colors, labels = percentile_labels, name = "Body Mass") +
  scale_fill_manual(values = percentile_colors, labels = percentile_labels, name = "Body Mass") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title    = "Marginal Effect of Range Size on Threat Probability by Body Mass",
    subtitle = "Birds: predictions at 10th/50th/90th percentile Body_Mass_log, other traits held at reference values",
    x        = "Scaled log10 range size",
    y        = "Predicted probability of threat",
    caption  = "Ribbons are 95% credible intervals. Phylogenetic and spatial random effects marginalised out (re_formula = NA)."
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 15.5, hjust = 0),
    axis.title.x       = element_text(size = 13, margin = margin(t = 8)),
    axis.title.y       = element_text(size = 13, margin = margin(r = 9)),
    strip.background   = element_rect(fill = "gray92", color = "gray45"),
    strip.text         = element_text(face = "bold", size = 12.5),
    panel.grid.minor   = element_blank(),
    plot.caption       = element_text(hjust = 0, size = 9.2, color = "gray35")
  )

print(marginal_plot)

output_png <- file.path(figure_dir, "bird_marginal_effects_body_mass_range_size.png")
output_pdf <- file.path(figure_dir, "bird_marginal_effects_body_mass_range_size.pdf")

ggsave(output_png, marginal_plot, width = 11, height = 6, dpi = 300)
ggsave(output_pdf, marginal_plot, width = 11, height = 6)

message("Saved marginal effects plot to ", output_png)
message("Saved marginal effects plot to ", output_pdf)
