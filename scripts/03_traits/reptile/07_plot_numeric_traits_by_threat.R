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

data_path <- file.path(
  "data", "processed", "traits", "reptile", "reptile_trait_data_matched.csv"
)
output_dir <- file.path(
  "outputs", "traits", "figures", "reptile_numeric_traits_by_threat"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reptile_traits <- read.csv(data_path, stringsAsFactors = FALSE)

to_numeric <- function(x) {
  suppressWarnings(as.numeric(trimws(x)))
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

trait_specs <- tibble::tribble(
  ~id, ~label, ~unit, ~column, ~transform,
  "clutch_size", "Mean offspring per litter or eggs per clutch", "", "Mean.number.of.offspring.per.litter.or.number.of.eggs.per.clutch", "none",
  "max_longevity", "Maximum longevity", "years", "Maximum.Longevity..years.", "none",
  "hatchling_neonate_mass", "Hatchling/neonate mass", "g", "Hatchling.neonate.mass..g.", "log10",
  "litters_clutches_per_year", "Litters or clutches per year", "per year", "Number.of.litters.or.clutches.produced.per.year", "none",
  "egg_length", "Egg length", "mm", "Egg.length..mm.", "none",
  "egg_width", "Egg width", "mm", "Egg.width..mm.", "none",
  "max_total_length", "Maximum total length", "mm", "Maximum.total.length...TL...mm.", "log10",
  "max_svl_scl", "Maximum SVL/SCL", "mm", "Maximum.length...SVL...mm..straight.carapace.length.for.turtles...SCL...mm.", "log10",
  "elevation_range", "Elevation range", "m", NA_character_, "none",
  "mean_tb", "Mean body temperature", "degrees C", "Mean.Tb", "none",
  "minimum_tb", "Minimum body temperature", "degrees C", "Minimum.Tb", "none",
  "maximum_tb", "Maximum body temperature", "degrees C", "Maximum.Tb", "none"
)

build_trait_values <- function(spec) {
  if (spec$id == "elevation_range") {
    min_elev <- to_numeric(reptile_traits[["Minimal.elevation..m."]])
    max_elev <- to_numeric(reptile_traits[["Maximum.elevation..m."]])
    values <- max_elev - min_elev
  } else {
    values <- to_numeric(reptile_traits[[spec$column]])
  }

  values
}

build_plot_data <- function(values) {
  indexed_traits <- reptile_traits %>%
    mutate(.row_id = row_number())

  bind_rows(
    indexed_traits %>%
      filter(Threat_Logging) %>%
      transmute(Group = "Logging threatened", value = values[.row_id]),
    indexed_traits %>%
      filter(Threat_Plantation) %>%
      transmute(Group = "Plantation threatened", value = values[.row_id]),
    indexed_traits %>%
      filter(!Threat_Logging, !Threat_Plantation) %>%
      transmute(Group = "Non-threatened", value = values[.row_id])
  ) %>%
    filter(is.finite(value)) %>%
    mutate(
      Group = factor(
        Group,
        levels = c("Logging threatened", "Plantation threatened", "Non-threatened")
      )
    )
}

plot_colors <- c(
  "Logging threatened" = "#c9a500",
  "Plantation threatened" = "#1287c9",
  "Non-threatened" = "#6e7781"
)

all_summaries <- list()

for (i in seq_len(nrow(trait_specs))) {
  spec <- trait_specs[i, ]
  values <- build_trait_values(spec)
  plot_data <- build_plot_data(values)

  if (nrow(plot_data) == 0) {
    warning("No finite values for ", spec$id, call. = FALSE)
    next
  }

  if (spec$transform == "log10") {
    plot_data <- plot_data %>% filter(value > 0)
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
      trait_id = spec$id,
      trait = spec$label,
      label = paste0("mean = ", format_number(mean), "\nmedian = ", format_number(median))
    )

  all_summaries[[spec$id]] <- summary_data

  axis_labels <- summary_data %>%
    mutate(axis_label = paste0(Group, "\n(n = ", n, ")")) %>%
    select(Group, axis_label) %>%
    tibble::deframe()

  y_limits <- range(plot_data$value, na.rm = TRUE)
  if (spec$transform == "log10") {
    log_limits <- log10(y_limits)
    y_text <- 10^(log_limits[1] + 0.92 * diff(log_limits))
  } else {
    y_text <- y_limits[1] + 0.92 * diff(y_limits)
  }

  summary_data <- summary_data %>% mutate(y_text = y_text)

  y_label <- if (nzchar(spec$unit)) {
    paste0(spec$label, " (", spec$unit, ")")
  } else {
    spec$label
  }

  trait_plot <- ggplot(plot_data, aes(x = Group, y = value, fill = Group)) +
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
    trait_plot <- trait_plot +
      scale_y_log10(
        labels = scales::label_number(),
        expand = expansion(mult = c(0.02, 0.12))
      )
  } else {
    trait_plot <- trait_plot +
      scale_y_continuous(
        labels = scales::label_number(),
        expand = expansion(mult = c(0.02, 0.12))
      )
  }

  png_path <- file.path(output_dir, paste0("reptile_", spec$id, "_by_threat_boxplot.png"))
  pdf_path <- file.path(output_dir, paste0("reptile_", spec$id, "_by_threat_boxplot.pdf"))

  ggsave(png_path, trait_plot, width = 8.8, height = 6.2, dpi = 300)
  ggsave(pdf_path, trait_plot, width = 8.8, height = 6.2)
}

summary_table <- bind_rows(all_summaries) %>%
  select(trait_id, trait, Group, n, mean, median, q25, q75, min, max)

summary_path <- file.path(output_dir, "reptile_numeric_traits_by_threat_summary.csv")
write.csv(summary_table, summary_path, row.names = FALSE)

message("Saved reptile numeric trait box plots to ", output_dir)
message("Saved summary table to ", summary_path)
