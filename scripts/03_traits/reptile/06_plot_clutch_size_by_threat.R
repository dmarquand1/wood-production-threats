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
figure_dir <- file.path("outputs", "traits", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

clutch_col <- "Mean.number.of.offspring.per.litter.or.number.of.eggs.per.clutch"

reptile_traits <- read.csv(data_path, stringsAsFactors = FALSE)

plot_data <- bind_rows(
  reptile_traits %>%
    filter(Threat_Logging) %>%
    transmute(
      Group = "Logging threatened",
      clutch_size = suppressWarnings(as.numeric(.data[[clutch_col]]))
    ),
  reptile_traits %>%
    filter(Threat_Plantation) %>%
    transmute(
      Group = "Plantation threatened",
      clutch_size = suppressWarnings(as.numeric(.data[[clutch_col]]))
    ),
  reptile_traits %>%
    filter(!Threat_Logging, !Threat_Plantation) %>%
    transmute(
      Group = "Non-threatened",
      clutch_size = suppressWarnings(as.numeric(.data[[clutch_col]]))
    )
) %>%
  filter(is.finite(clutch_size), clutch_size > 0) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Logging threatened", "Plantation threatened", "Non-threatened")
    )
  )

group_counts <- plot_data %>%
  count(Group, name = "n") %>%
  mutate(axis_label = paste0(Group, "\n(n = ", n, ")"))

axis_labels <- setNames(group_counts$axis_label, group_counts$Group)

summary_path <- file.path(figure_dir, "reptile_clutch_size_by_threat_summary.csv")
plot_summary <- plot_data %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    median = median(clutch_size),
    mean = mean(clutch_size),
    q25 = quantile(clutch_size, 0.25),
    q75 = quantile(clutch_size, 0.75),
    min = min(clutch_size),
    max = max(clutch_size),
    .groups = "drop"
  )
write.csv(plot_summary, summary_path, row.names = FALSE)

plot_colors <- c(
  "Logging threatened" = "#c9a500",
  "Plantation threatened" = "#1287c9",
  "Non-threatened" = "#6e7781"
)

clutch_plot <- ggplot(plot_data, aes(x = Group, y = clutch_size, fill = Group)) +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.75,
    outlier.alpha = 0.35,
    outlier.size = 1.6
  ) +
  scale_fill_manual(values = plot_colors, guide = "none") +
  scale_x_discrete(labels = axis_labels) +
  scale_y_continuous(
    breaks = scales::pretty_breaks(n = 8),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Mean offspring per litter or eggs per clutch"
  ) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_text(size = 13, color = "gray20"),
    axis.text.y = element_text(size = 12, color = "gray20"),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 14, 12, 12)
  )

output_png <- file.path(figure_dir, "reptile_clutch_size_by_threat_boxplot.png")
output_pdf <- file.path(figure_dir, "reptile_clutch_size_by_threat_boxplot.pdf")

ggsave(output_png, clutch_plot, width = 8.5, height = 6, dpi = 300)
ggsave(output_pdf, clutch_plot, width = 8.5, height = 6)

message("Saved plot to ", output_png)
message("Saved plot to ", output_pdf)
message("Saved summary to ", summary_path)
