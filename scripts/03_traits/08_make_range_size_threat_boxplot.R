rm(list = ls())
gc()

project_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib, winslash = "/", mustWork = FALSE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

figure_dir <- file.path("outputs", "traits", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- tibble::tribble(
  ~taxon,       ~Taxon,
  "amphibian",  "Amphibians",
  "bird",       "Birds",
  "fish",       "Fish",
  "mammal",     "Mammals",
  "reptile",    "Reptiles"
)

threat_levels <- c(
  "Logging threatened",
  "Plantation threatened",
  "Not threatened"
)

threat_colors <- c(
  "Logging threatened" = "#d9be00",
  "Plantation threatened" = "#0088d9",
  "Not threatened" = "#7a7a7a"
)

as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }

  if (is.numeric(x)) {
    return(replace_na(x != 0, FALSE))
  }

  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

pick_required_column <- function(df, candidates, path) {
  matching <- candidates[candidates %in% names(df)]
  if (length(matching) == 0) {
    stop(
      "Missing required column in ", path, ": ",
      paste(candidates, collapse = " / "),
      call. = FALSE
    )
  }
  matching[[1]]
}

read_taxon_data <- function(taxon, taxon_label) {
  path <- file.path(
    "data",
    "processed",
    "traits",
    taxon,
    sprintf("%s_trait_model_final_data.csv", taxon)
  )

  if (!file.exists(path)) {
    stop("Missing model data file: ", path, call. = FALSE)
  }

  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  range_col <- pick_required_column(df, "log_range_size", path)
  logging_col <- pick_required_column(df, c("Threat_Logging", "Threat_Logging.x"), path)
  plantation_col <- pick_required_column(df, c("Threat_Plantation", "Threat_Plantation.x"), path)
  species_col <- pick_required_column(df, c("spatial_name", "phylo_name"), path)

  tibble(
    taxon = taxon,
    Taxon = taxon_label,
    species = df[[species_col]],
    log_range_size = suppressWarnings(as.numeric(df[[range_col]])),
    Threat_Logging = as_logical_flag(df[[logging_col]]),
    Threat_Plantation = as_logical_flag(df[[plantation_col]])
  ) %>%
    filter(is.finite(log_range_size))
}

model_data <- bind_rows(Map(read_taxon_data, taxa$taxon, taxa$Taxon)) %>%
  mutate(Taxon = factor(Taxon, levels = taxa$Taxon))

plot_data <- bind_rows(
  model_data %>%
    filter(Threat_Logging) %>%
    mutate(threat_group = "Logging threatened"),
  model_data %>%
    filter(Threat_Plantation) %>%
    mutate(threat_group = "Plantation threatened"),
  model_data %>%
    filter(!Threat_Logging, !Threat_Plantation) %>%
    mutate(threat_group = "Not threatened")
) %>%
  mutate(threat_group = factor(threat_group, levels = threat_levels))

range_plot <- ggplot(plot_data, aes(x = threat_group, y = log_range_size, fill = threat_group)) +
  geom_boxplot(
    width = 0.68,
    outlier.alpha = 0.25,
    outlier.size = 0.9,
    linewidth = 0.45,
    color = "gray25"
  ) +
  facet_wrap(~ Taxon, nrow = 2) +
  scale_fill_manual(values = threat_colors, guide = "none") +
  scale_x_discrete(labels = threat_levels) +
  labs(
    title = "Range size by logging and plantation threat status",
    x = NULL,
    y = expression("Scaled log"[10] * " range size (km"^2 * ")"),
    caption = "Logging and plantation groups are non-mutually exclusive; species with both threats appear in both boxes. Not threatened = neither logging nor plantation flagged."
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1),
    strip.background = element_rect(fill = "gray92", color = "gray45"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9, color = "gray35"),
    plot.margin = margin(8, 12, 8, 8)
  )

output_png <- file.path(figure_dir, "all_taxa_range_size_by_threat_boxplot.png")
output_pdf <- file.path(figure_dir, "all_taxa_range_size_by_threat_boxplot.pdf")

ggsave(output_png, range_plot, width = 11, height = 7, dpi = 300)
ggsave(output_pdf, range_plot, width = 11, height = 7)

message("Saved range-size boxplot to ", output_png)
message("Saved range-size boxplot to ", output_pdf)
