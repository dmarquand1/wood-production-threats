rm(list = ls())
gc()

project_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib, winslash = "/", mustWork = FALSE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

figure_root <- file.path("outputs", "traits", "figures")

data_paths <- c(
  amphibian = file.path("data", "processed", "traits", "amphibian", "amphibian_trait_data_matched.csv"),
  bird = file.path("data", "processed", "traits", "bird", "bird_trait_data_matched.csv"),
  fish = file.path("data", "processed", "traits", "fish", "fish_trait_data_matched.csv"),
  mammal = file.path("data", "processed", "traits", "mammal", "mammal_trait_data_matched.csv")
)

trait_specs <- data.frame(
  taxon = c(
    rep("amphibian", 10),
    rep("bird", 8),
    rep("fish", 6),
    rep("mammal", 17)
  ),
  id = c(
    "body_mass", "age_at_maturity_min", "age_at_maturity_max", "longevity_max",
    "litter_size_min", "litter_size_max", "reproductive_output",
    "offspring_size_min", "offspring_size_max", "size_at_maturity_min",
    "elevational_range", "incubation_min", "incubation_max", "fledging_min",
    "fledging_max", "breeding_season_start", "breeding_season_end", "productivity_min",
    "vulnerability", "pd50", "max_standard_length", "reported_length",
    "growth_k", "predator_prey_ratio_min",
    "adult_body_length", "adult_brain_mass", "max_longevity", "maturity",
    "age_first_reproduction", "gestation_length", "litters_per_year",
    "interbirth_interval", "neonate_mass", "weaning_age", "weaning_mass",
    "generation_length", "dispersal", "density", "home_range",
    "social_group", "altitude_breadth"
  ),
  label = c(
    "Body mass", "Minimum age at maturity", "Maximum age at maturity", "Maximum longevity",
    "Minimum litter/clutch size", "Maximum litter/clutch size", "Reproductive output",
    "Minimum offspring size", "Maximum offspring size", "Minimum size at maturity",
    "Elevational range", "Minimum incubation period", "Maximum incubation period", "Minimum fledging period",
    "Maximum fledging period", "Breeding season start", "Breeding season end", "Minimum productivity",
    "Vulnerability", "PD50", "Maximum standard length", "Reported maximum length",
    "Growth coefficient K", "Minimum predator:prey ratio",
    "Adult body length", "Adult brain mass", "Maximum longevity", "Maturity age",
    "Age at first reproduction", "Gestation length", "Litters per year",
    "Interbirth interval", "Neonate mass", "Weaning age", "Weaning mass",
    "Generation length", "Dispersal distance", "Population density", "Home range",
    "Social group size", "Elevation range"
  ),
  unit = c(
    "g", "years", "years", "years",
    "", "", "per year", "mm", "mm", "mm",
    "m", "days", "days", "days",
    "days", "", "", "",
    "", "", "cm", "cm",
    "", "", "mm", "g", "days", "days",
    "days", "days", "per year", "days", "g", "days", "g",
    "days", "km", "n/km2", "km2", "n", "m"
  ),
  column = c(
    "Body_mass_g", "Age_at_maturity_min_y", "Age_at_maturity_max_y", "Longevity_max_y",
    "Litter_size_min_n", "Litter_size_max_n", "Reproductive_output_y",
    "Offspring_size_min_mm", "Offspring_size_max_mm", "Size_at_maturity_min_mm",
    "Elevational.Range", "Incu1", "Incu2", "Fldg1",
    "Fldg2", "Brd1", "Brd2", "Prod1",
    "species_Vulnerability", "species_PD50", "estimate_MaxLengthSL", "species_Length",
    "estimate_K", "estimate_PredPreyRatioMin",
    "adult_body_length_mm", "adult_brain_mass_g", "max_longevity_d", "maturity_d",
    "age_first_reproduction_d", "gestation_length_d", "litters_per_year_n",
    "interbirth_interval_d", "neonate_mass_g", "weaning_age_d", "weaning_mass_g",
    "generation_length_d", "dispersal_km", "density_n_km2", "home_range_km2",
    "social_group_n", "altitude_breadth_m"
  ),
  transform = c(
    "log10", "none", "none", "none",
    "log10", "log10", "none", "none", "none", "none",
    "none", "none", "none", "none", "none", "none", "none", "none",
    "none", "none", "log10", "log10", "log10", "log10",
    "log10", "log10", "none", "none",
    "none", "none", "none", "none", "log10", "none", "log10",
    "none", "log10", "log10", "log10", "log10", "none"
  ),
  stringsAsFactors = FALSE
)

to_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(x)))
}

to_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

format_number <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(abs(x) >= 100, format(round(x, 0), big.mark = ","),
      ifelse(abs(x) >= 10, sprintf("%.1f", x), sprintf("%.2f", x))
    )
  )
}

plot_colors <- c(
  "Logging threatened" = "#c9a500",
  "Plantation threatened" = "#1287c9",
  "Non-threatened" = "#6e7781"
)

build_plot_data <- function(taxon_data, values) {
  logging <- to_flag(taxon_data$Threat_Logging)
  plantation <- to_flag(taxon_data$Threat_Plantation)

  bind_rows(
    data.frame(Group = "Logging threatened", value = values[logging]),
    data.frame(Group = "Plantation threatened", value = values[plantation]),
    data.frame(Group = "Non-threatened", value = values[!logging & !plantation])
  ) %>%
    filter(is.finite(value)) %>%
    mutate(
      Group = factor(
        Group,
        levels = c("Logging threatened", "Plantation threatened", "Non-threatened")
      )
    )
}

make_trait_plot <- function(plot_data, summary_data, spec) {
  axis_labels <- summary_data %>%
    mutate(axis_label = paste0(Group, "\n(n = ", n, ")")) %>%
    select(Group, axis_label) %>%
    tibble::deframe()

  y_limits <- range(plot_data$value, na.rm = TRUE)
  if (!is.finite(diff(y_limits)) || diff(y_limits) == 0) {
    y_limits <- y_limits + c(-0.5, 0.5)
  }

  if (spec$transform == "log10") {
    log_limits <- log10(y_limits)
    y_text <- 10^(log_limits[1] + 0.92 * diff(log_limits))
  } else {
    y_text <- y_limits[1] + 0.92 * diff(y_limits)
  }

  summary_data <- summary_data %>%
    mutate(
      y_text = y_text,
      label = paste0("mean = ", format_number(mean), "\nmedian = ", format_number(median))
    )

  y_label <- if (nzchar(spec$unit)) {
    paste0(spec$label, " (", spec$unit, ")")
  } else {
    spec$label
  }

  p <- ggplot(plot_data, aes(x = Group, y = value, fill = Group)) +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.75,
      outlier.alpha = 0.32,
      outlier.size = 1.35
    ) +
    geom_text(
      data = summary_data,
      aes(x = Group, y = y_text, label = label),
      inherit.aes = FALSE,
      size = 3.8,
      lineheight = 0.95,
      color = "gray15"
    ) +
    scale_fill_manual(values = plot_colors, guide = "none") +
    scale_x_discrete(labels = axis_labels) +
    labs(x = NULL, y = y_label) +
    theme_bw(base_size = 15) +
    theme(
      axis.text.x = element_text(size = 12.5, color = "gray20"),
      axis.text.y = element_text(size = 12, color = "gray20"),
      axis.title.y = element_text(size = 13.5, margin = margin(r = 10)),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 16, 12, 12)
    )

  if (spec$transform == "log10") {
    p +
      scale_y_log10(
        labels = scales::label_number(),
        expand = expansion(mult = c(0.02, 0.12))
      )
  } else {
    p +
      scale_y_continuous(
        labels = scales::label_number(),
        expand = expansion(mult = c(0.02, 0.12))
      )
  }
}

all_summaries <- list()

for (taxon in names(data_paths)) {
  taxon_data <- read.csv(data_paths[[taxon]], stringsAsFactors = FALSE, check.names = FALSE)
  taxon_specs <- trait_specs %>% filter(.data$taxon == taxon)
  output_dir <- file.path(figure_root, paste0(taxon, "_numeric_traits_by_threat"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(taxon_specs))) {
    spec <- taxon_specs[i, ]
    if (!spec$column %in% names(taxon_data)) {
      warning("Missing column for ", taxon, ": ", spec$column, call. = FALSE)
      next
    }

    values <- to_numeric(taxon_data[[spec$column]])
    plot_data <- build_plot_data(taxon_data, values)

    if (spec$transform == "log10") {
      plot_data <- plot_data %>% filter(value > 0)
    }

    if (nrow(plot_data) == 0) {
      warning("No finite values for ", taxon, " ", spec$id, call. = FALSE)
      next
    }

    summary_data <- plot_data %>%
      group_by(Group) %>%
      summarise(
        n = n(),
        mean = mean(value),
        median = median(value),
        q25 = quantile(value, 0.25),
        q75 = quantile(value, 0.75),
        min = min(value),
        max = max(value),
        .groups = "drop"
      ) %>%
      mutate(
        taxon = taxon,
        trait_id = spec$id,
        trait = spec$label
      )

    all_summaries[[paste(taxon, spec$id, sep = "_")]] <- summary_data

    trait_plot <- make_trait_plot(plot_data, summary_data, spec)

    png_path <- file.path(output_dir, paste0(taxon, "_", spec$id, "_by_threat_boxplot.png"))
    pdf_path <- file.path(output_dir, paste0(taxon, "_", spec$id, "_by_threat_boxplot.pdf"))

    ggsave(png_path, trait_plot, width = 8.8, height = 6.2, dpi = 300)
    ggsave(pdf_path, trait_plot, width = 8.8, height = 6.2)
  }

  taxon_summary <- bind_rows(all_summaries) %>%
    filter(.data$taxon == taxon) %>%
    select(taxon, trait_id, trait, Group, n, mean, median, q25, q75, min, max)

  write.csv(
    taxon_summary,
    file.path(output_dir, paste0(taxon, "_numeric_traits_by_threat_summary.csv")),
    row.names = FALSE
  )
}

combined_summary <- bind_rows(all_summaries) %>%
  select(taxon, trait_id, trait, Group, n, mean, median, q25, q75, min, max)

summary_path <- file.path(
  figure_root,
  "other_taxa_numeric_traits_by_threat_summary.csv"
)
write.csv(combined_summary, summary_path, row.names = FALSE)

message("Saved other-taxa numeric trait box plots under ", figure_root)
message("Saved combined summary table to ", summary_path)
