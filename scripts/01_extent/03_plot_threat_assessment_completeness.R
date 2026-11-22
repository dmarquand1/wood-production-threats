rm(list = ls())
gc()

# =============================================================================
# PLOT RED LIST THREAT-ASSESSMENT COMPLETENESS
#
# Uses outputs from 02_summarise_threat_assessment_completeness.R.
#
# RStudio:
#   Click Source. Plots are written to outputs/assessment_completeness/figures.
# =============================================================================

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
library(stringr)
library(tidyr)

input_dir <- file.path("outputs", "assessment_completeness")
figure_dir <- file.path(input_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

species_flags_file <- file.path(input_dir, "assessment_threat_completeness_species_flags.csv")
raw_threat_rows_file <- file.path(input_dir, "assessment_threat_completeness_raw_threat_rows.csv")

if (!file.exists(species_flags_file) || !file.exists(raw_threat_rows_file)) {
  stop(
    "Missing completeness detail files. Run scripts/01_extent/",
    "02_summarise_threat_assessment_completeness.R first."
  )
}

taxon_labels <- c(
  amphibian = "Amphibians",
  bird = "Birds",
  fish = "Freshwater fish",
  mammal = "Mammals",
  reptile = "Reptiles"
)

taxon_levels <- names(taxon_labels)

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1")
}

percent_lab <- function(x) paste0(round(100 * x), "%")

species_flags <- read.csv(species_flags_file, stringsAsFactors = FALSE) %>%
  mutate(
    taxon_key = factor(taxon_key, levels = taxon_levels, labels = taxon_labels),
    has_any_threat = to_logical(has_any_threat),
    has_scope_or_severity = to_logical(has_scope_or_severity)
  )

raw_threat_rows <- read.csv(raw_threat_rows_file, stringsAsFactors = FALSE) %>%
  mutate(
    taxon_key = factor(taxon_key, levels = taxon_levels, labels = taxon_labels),
    threat_code = as.character(threat_code),
    timing = as.character(timing),
    scope = na_if(str_squish(as.character(scope)), ""),
    severity = na_if(str_squish(as.character(severity)), ""),
    score = na_if(str_squish(as.character(score)), "")
  )

# Plot 1: assessment-detail classes by taxon.
# "Major threats only" is operationalised here as listed threat codes but no
# scope/severity fields. "Major and minor threat detail" means at least one
# listed threat row has scope or severity information.
assessment_detail <- species_flags %>%
  filter(has_any_threat) %>%
  mutate(
    detail_class = if_else(
      has_scope_or_severity,
      "Threats with scope/severity detail",
      "Threats listed only"
    ),
    detail_class = factor(
      detail_class,
      levels = c("Threats listed only", "Threats with scope/severity detail")
    )
  ) %>%
  count(taxon_key, detail_class, name = "n_species") %>%
  group_by(taxon_key) %>%
  mutate(prop_species = n_species / sum(n_species)) %>%
  ungroup()

assessment_counts_plot <- ggplot(
  assessment_detail,
  aes(x = taxon_key, y = n_species, fill = detail_class)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = scales::comma(n_species)),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    colour = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Threats listed only" = "#6B7280",
      "Threats with scope/severity detail" = "#2A9D8F"
    )
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = NULL,
    y = "Number of species with listed threats",
    fill = NULL,
    title = "Threat-assessment detail by taxon",
    subtitle = "Species with listed threats only vs species with scope/severity information"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

assessment_proportions_plot <- ggplot(
  assessment_detail,
  aes(x = taxon_key, y = prop_species, fill = detail_class)
) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  geom_text(
    aes(label = percent_lab(prop_species)),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    colour = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Threats listed only" = "#6B7280",
      "Threats with scope/severity detail" = "#2A9D8F"
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = NULL,
    y = "Proportion of species with listed threats",
    fill = NULL,
    title = "Threat-assessment detail by taxon",
    subtitle = "Proportions among species with at least one listed threat"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    plot.title.position = "plot"
  )

ggsave(
  file.path(figure_dir, "threat_assessment_detail_counts_by_taxon.png"),
  assessment_counts_plot,
  width = 9,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(figure_dir, "threat_assessment_detail_proportions_by_taxon.png"),
  assessment_proportions_plot,
  width = 9,
  height = 6,
  dpi = 300
)

# Plot 2: severity/scope/impact-score categories for ongoing wood-production
# threats, matching the modelling response as closely as possible.
wood_rows <- raw_threat_rows %>%
  filter(
    str_detect(threat_code, "^(5_3|2_2)"),
    timing == "Ongoing"
  ) %>%
  mutate(
    scope = replace_na(scope, "Missing"),
    severity = replace_na(severity, "Missing"),
    score = replace_na(score, "Missing"),
    score_group = case_when(
      str_detect(score, "^High Impact") ~ "High impact",
      str_detect(score, "^Medium Impact") ~ "Medium impact",
      str_detect(score, "^Low Impact") ~ "Low impact",
      str_detect(score, "^No/Negligible Impact") ~ "No/negligible impact",
      str_detect(score, "^Past Impact") ~ "Past impact",
      score %in% c("Unknown", "(Not specified)", "Missing") ~ score,
      TRUE ~ score
    ),
    threat_type = case_when(
      str_detect(threat_code, "^5_3") ~ "Logging",
      str_detect(threat_code, "^2_2") ~ "Plantations",
      TRUE ~ "Other wood production"
    )
  )

scope_levels <- c("Whole (>90%)", "Majority (50-90%)", "Minority (<50%)", "Unknown", "Missing")
severity_levels <- c(
  "Very Rapid Declines",
  "Rapid Declines",
  "Slow, Significant Declines",
  "Causing/Could cause fluctuations",
  "Negligible declines",
  "No decline",
  "Unknown",
  "Missing"
)
score_levels <- c(
  "High impact",
  "Medium impact",
  "Low impact",
  "No/negligible impact",
  "Past impact",
  "Unknown",
  "(Not specified)",
  "Missing"
)

wood_category_long <- bind_rows(
  wood_rows %>%
    transmute(taxon_key, threat_type, category_type = "Scope", category = factor(scope, levels = scope_levels)),
  wood_rows %>%
    transmute(taxon_key, threat_type, category_type = "Severity", category = factor(severity, levels = severity_levels)),
  wood_rows %>%
    transmute(taxon_key, threat_type, category_type = "Impact score", category = factor(score_group, levels = score_levels))
) %>%
  filter(!is.na(taxon_key), !is.na(category)) %>%
  count(taxon_key, threat_type, category_type, category, name = "n_rows") %>%
  group_by(taxon_key, threat_type, category_type) %>%
  mutate(prop_rows = n_rows / sum(n_rows)) %>%
  ungroup() %>%
  mutate(
    category_type = factor(category_type, levels = c("Scope", "Severity", "Impact score")),
    threat_type = factor(threat_type, levels = c("Logging", "Plantations"))
  )

wood_scope_severity_plot <- ggplot(
  wood_category_long,
  aes(x = category, y = prop_rows, fill = threat_type)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  facet_grid(category_type ~ taxon_key, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Logging" = "#4C78A8", "Plantations" = "#E76F51")) +
  labs(
    x = NULL,
    y = "Proportion of ongoing wood-production threat rows",
    fill = NULL,
    title = "Scope, severity and impact scores for ongoing wood-production threats",
    subtitle = "Rows are logging or timber-plantation threat codes with timing = Ongoing"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot"
  )

ggsave(
  file.path(figure_dir, "wood_production_scope_severity_score_by_taxon.png"),
  wood_scope_severity_plot,
  width = 15,
  height = 9,
  dpi = 300
)

write.csv(
  assessment_detail,
  file.path(input_dir, "plot_data_threat_assessment_detail_by_taxon.csv"),
  row.names = FALSE
)

write.csv(
  wood_category_long,
  file.path(input_dir, "plot_data_wood_production_scope_severity_score_by_taxon.csv"),
  row.names = FALSE
)

cat("\nWrote plots to: ", figure_dir, "\n", sep = "")
cat("  threat_assessment_detail_counts_by_taxon.png\n")
cat("  threat_assessment_detail_proportions_by_taxon.png\n")
cat("  wood_production_scope_severity_score_by_taxon.png\n")

# Plot 3: bivariate scope/severity heatmap for ongoing wood-production threats.
wood_scope_severity_heatmap_data <- wood_rows %>%
  mutate(
    scope = factor(replace_na(scope, "Missing"), levels = rev(scope_levels)),
    severity = factor(replace_na(severity, "Missing"), levels = severity_levels),
    threat_type = factor(threat_type, levels = c("Logging", "Plantations"))
  ) %>%
  filter(!is.na(taxon_key), !is.na(scope), !is.na(severity)) %>%
  count(taxon_key, threat_type, scope, severity, name = "n_rows") %>%
  group_by(taxon_key, threat_type) %>%
  mutate(prop_rows = n_rows / sum(n_rows)) %>%
  ungroup()

wood_scope_severity_heatmap <- ggplot(
  wood_scope_severity_heatmap_data,
  aes(x = severity, y = scope, fill = n_rows)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = n_rows), size = 2.7, colour = "white", fontface = "bold") +
  facet_grid(threat_type ~ taxon_key) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",
    labels = scales::comma,
    name = "Threat rows"
  ) +
  labs(
    x = "Severity",
    y = "Scope",
    title = "Scope and severity combinations for ongoing wood-production threats",
    subtitle = "Colour and labels show the number of logging or plantation threat rows"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot"
  )

ggsave(
  file.path(figure_dir, "wood_production_scope_by_severity_heatmap.png"),
  wood_scope_severity_heatmap,
  width = 15,
  height = 8,
  dpi = 300
)

write.csv(
  wood_scope_severity_heatmap_data,
  file.path(input_dir, "plot_data_wood_production_scope_by_severity_heatmap.csv"),
  row.names = FALSE
)

cat("  wood_production_scope_by_severity_heatmap.png\n")
