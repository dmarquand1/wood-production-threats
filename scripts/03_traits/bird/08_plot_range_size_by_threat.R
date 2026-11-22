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
  "data",
  "processed",
  "traits",
  "bird",
  "bird_trait_model_final_data.csv"
)

figure_dir <- file.path(
  "outputs",
  "traits",
  "figures",
  "bird_range_size_by_threat"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

to_flag <- function(x) {
  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }

  if (is.numeric(x)) {
    return(!is.na(x) & x != 0)
  }

  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

format_number <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(abs(x) >= 10, sprintf("%.1f", x), sprintf("%.2f", x))
  )
}

bird_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = TRUE) %>%
  transmute(
    species = spatial_name,
    body_mass_log = suppressWarnings(as.numeric(Body_Mass_log)),
    log_range_size = suppressWarnings(as.numeric(log_range_size)),
    logging = to_flag(Threat_Logging),
    timber_plantations = to_flag(Threat_Plantation)
  ) %>%
  filter(is.finite(body_mass_log), is.finite(log_range_size))

body_mass_breaks <- quantile(
  bird_data$body_mass_log,
  probs = c(0, 1 / 3, 2 / 3, 1),
  na.rm = TRUE
)

if (anyDuplicated(body_mass_breaks) > 0) {
  stop("Cannot create body-mass tertiles because Body_Mass_log quantile breaks are duplicated.", call. = FALSE)
}

body_mass_labels <- c(
  paste0("Small\n", format_number(body_mass_breaks[1]), " to ", format_number(body_mass_breaks[2])),
  paste0("Medium\n", format_number(body_mass_breaks[2]), " to ", format_number(body_mass_breaks[3])),
  paste0("Large\n", format_number(body_mass_breaks[3]), " to ", format_number(body_mass_breaks[4]))
)

binned_base <- bird_data %>%
  mutate(
    body_mass_bin = cut(
      body_mass_log,
      breaks = body_mass_breaks,
      labels = body_mass_labels,
      include.lowest = TRUE
    )
  )

plot_data_binned <- bind_rows(
  binned_base %>%
    transmute(
      species,
      body_mass_bin,
      threat = "Logging",
      status = if_else(logging, "Threatened", "Not threatened"),
      log_range_size
    ),
  binned_base %>%
    transmute(
      species,
      body_mass_bin,
      threat = "Timber plantations",
      status = if_else(timber_plantations, "Threatened", "Not threatened"),
      log_range_size
    )
) %>%
  mutate(
    threat = factor(threat, levels = c("Logging", "Timber plantations")),
    status = factor(status, levels = c("Not threatened", "Threatened")),
    fill_group = case_when(
      status == "Not threatened" ~ "Not threatened",
      threat == "Logging" ~ "Threatened by logging",
      TRUE ~ "Threatened by timber plantations"
    )
  )

binned_summary <- plot_data_binned %>%
  group_by(threat, body_mass_bin, status, fill_group) %>%
  summarise(
    n = n(),
    mean = mean(log_range_size),
    median = median(log_range_size),
    sd = sd(log_range_size),
    se = sd / sqrt(n),
    ci95 = qt(0.975, df = pmax(n - 1, 1)) * se,
    q25 = quantile(log_range_size, 0.25),
    q75 = quantile(log_range_size, 0.75),
    min = min(log_range_size),
    max = max(log_range_size),
    .groups = "drop"
  ) %>%
  mutate(
    ci95 = if_else(n > 1, ci95, 0),
    label_y = mean + ci95 + 0.08
  )

bar_colors <- c(
  "Not threatened" = "#757575",
  "Threatened by logging" = "#c7a600",
  "Threatened by timber plantations" = "#1287c9"
)

body_mass_bar_plot <- ggplot(
  binned_summary,
  aes(x = body_mass_bin, y = mean, fill = fill_group)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "gray55") +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.64,
    color = "gray24",
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(ymin = mean - ci95, ymax = mean + ci95),
    position = position_dodge(width = 0.72),
    width = 0.18,
    linewidth = 0.45,
    color = "gray24"
  ) +
  geom_text(
    aes(y = label_y, label = paste0("n=", n)),
    position = position_dodge(width = 0.72),
    size = 3.1,
    color = "gray18",
    vjust = 0
  ) +
  facet_wrap(~ threat, nrow = 1) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
  labs(
    title = "Threatened vs non-threatened range size by bird body-mass bin",
    x = expression("Body-mass bin (tertiles of " * log[10] * " body mass)"),
    y = expression("Mean scaled log"[10] * " range size"),
    caption = "Bars show means; error bars are 95% CIs. Body-mass bins use the model predictor Body_Mass_log."
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, hjust = 0),
    axis.text.x = element_text(size = 10.5, color = "gray20"),
    axis.text.y = element_text(size = 11.5, color = "gray20"),
    axis.title.x = element_text(size = 12.5, margin = margin(t = 8)),
    axis.title.y = element_text(size = 13, margin = margin(r = 9)),
    legend.position = "bottom",
    legend.text = element_text(size = 10.5),
    strip.background = element_rect(fill = "gray92", color = "gray45"),
    strip.text = element_text(face = "bold", size = 12.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9.2, color = "gray35"),
    plot.margin = margin(10, 14, 10, 10)
  )

bar_png_path <- file.path(figure_dir, "bird_log_range_size_by_threat_and_body_mass_bins_barplot.png")
bar_pdf_path <- file.path(figure_dir, "bird_log_range_size_by_threat_and_body_mass_bins_barplot.pdf")
bar_summary_path <- file.path(figure_dir, "bird_log_range_size_by_threat_and_body_mass_bins_summary.csv")

ggsave(bar_png_path, body_mass_bar_plot, width = 10.4, height = 6.2, dpi = 300)
ggsave(bar_pdf_path, body_mass_bar_plot, width = 10.4, height = 6.2)
write.csv(
  binned_summary %>%
    mutate(
      threat = as.character(threat),
      body_mass_bin = as.character(body_mass_bin),
      status = as.character(status)
    ) %>%
    select(threat, body_mass_bin, status, n, mean, median, sd, se, ci95, q25, q75, min, max),
  bar_summary_path,
  row.names = FALSE
)

scatter_data <- bind_rows(
  bird_data %>%
    transmute(
      species,
      body_mass_log,
      log_range_size,
      threat = "Logging",
      status = if_else(logging, "Threatened", "Not threatened")
    ),
  bird_data %>%
    transmute(
      species,
      body_mass_log,
      log_range_size,
      threat = "Timber plantations",
      status = if_else(timber_plantations, "Threatened", "Not threatened")
    )
) %>%
  mutate(
    threat = factor(threat, levels = c("Logging", "Timber plantations")),
    status = factor(status, levels = c("Not threatened", "Threatened")),
    color_group = case_when(
      status == "Not threatened" ~ "Not threatened",
      threat == "Logging" ~ "Threatened by logging",
      TRUE ~ "Threatened by timber plantations"
    )
  )

scatter_colors <- c(
  "Not threatened" = "#6f6f6f",
  "Threatened by logging" = "#c7a600",
  "Threatened by timber plantations" = "#1287c9"
)

interaction_scatter <- ggplot(
  scatter_data,
  aes(x = body_mass_log, y = log_range_size, color = color_group)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "gray60") +
  geom_point(alpha = 0.18, size = 0.9, stroke = 0) +
  geom_smooth(
    aes(group = color_group),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = TRUE,
    linewidth = 1.05,
    alpha = 0.16
  ) +
  facet_wrap(~ threat, nrow = 1) +
  scale_color_manual(values = scatter_colors, name = NULL) +
  labs(
    title = "Body mass and range size by bird threat status",
    x = expression(log[10] * " body mass"),
    y = expression("Scaled log"[10] * " range size"),
    caption = "Curves are quadratic fits with 95% confidence bands. Threat flags are evaluated separately for each panel."
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, hjust = 0),
    axis.text = element_text(size = 11.5, color = "gray20"),
    axis.title.x = element_text(size = 12.5, margin = margin(t = 8)),
    axis.title.y = element_text(size = 13, margin = margin(r = 9)),
    legend.position = "bottom",
    legend.text = element_text(size = 10.5),
    strip.background = element_rect(fill = "gray92", color = "gray45"),
    strip.text = element_text(face = "bold", size = 12.5),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9.2, color = "gray35"),
    plot.margin = margin(10, 14, 10, 10)
  )

scatter_png_path <- file.path(figure_dir, "bird_body_mass_range_size_by_threat_scatter.png")
scatter_pdf_path <- file.path(figure_dir, "bird_body_mass_range_size_by_threat_scatter.pdf")

ggsave(scatter_png_path, interaction_scatter, width = 10.4, height = 6.2, dpi = 300)
ggsave(scatter_pdf_path, interaction_scatter, width = 10.4, height = 6.2)

message("Saved body-mass binned bar plot to ", bar_png_path)
message("Saved body-mass binned bar plot to ", bar_pdf_path)
message("Saved body-mass binned summary table to ", bar_summary_path)
message("Saved body-mass/range-size threat scatter plot to ", scatter_png_path)
message("Saved body-mass/range-size threat scatter plot to ", scatter_pdf_path)
