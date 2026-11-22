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
  "amphibian",
  "amphibian_trait_model_final_data.csv"
)

figure_dir <- file.path(
  "outputs",
  "traits",
  "figures",
  "amphibian_range_size_by_threat"
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

amphibian_data <- read.csv(data_path, stringsAsFactors = FALSE, check.names = TRUE) %>%
  transmute(
    species = spatial_name,
    body_size_log = suppressWarnings(as.numeric(Body_Size_log)),
    log_range_size = suppressWarnings(as.numeric(log_range_size)),
    logging = to_flag(Threat_Logging),
    timber_plantations = to_flag(Threat_Plantation)
  ) %>%
  filter(is.finite(log_range_size))

plot_data <- bind_rows(
  amphibian_data %>%
    transmute(
      species,
      threat = "Logging",
      status = if_else(logging, "Threatened", "Not threatened"),
      log_range_size
    ),
  amphibian_data %>%
    transmute(
      species,
      threat = "Timber plantations",
      status = if_else(timber_plantations, "Threatened", "Not threatened"),
      log_range_size
    )
) %>%
  mutate(
    threat = factor(threat, levels = c("Logging", "Timber plantations")),
    status = factor(status, levels = c("Not threatened", "Threatened")),
    plot_group = paste(threat, status, sep = " | ")
  )

summary_data <- plot_data %>%
  group_by(threat, status) %>%
  summarise(
    n = n(),
    mean = mean(log_range_size),
    median = median(log_range_size),
    q25 = quantile(log_range_size, 0.25),
    q75 = quantile(log_range_size, 0.75),
    min = min(log_range_size),
    max = max(log_range_size),
    .groups = "drop"
  ) %>%
  mutate(
    axis_label = paste0(status, "\n(n = ", n, ")"),
    stats_label = paste0("mean = ", format_number(mean), "\nmedian = ", format_number(median))
  )

axis_labels <- summary_data %>%
  mutate(axis_key = paste(threat, status, sep = " | ")) %>%
  select(axis_key, axis_label) %>%
  tibble::deframe()

y_limits <- range(plot_data$log_range_size, na.rm = TRUE)
y_text <- y_limits[1] + 0.94 * diff(y_limits)

summary_data <- summary_data %>%
  mutate(
    y_text = y_text,
    plot_group = paste(threat, status, sep = " | ")
  )

plot_colors <- c(
  "Logging | Not threatened" = "#757575",
  "Logging | Threatened" = "#c7a600",
  "Timber plantations | Not threatened" = "#757575",
  "Timber plantations | Threatened" = "#1287c9"
)

range_plot <- ggplot(
  plot_data,
  aes(x = plot_group, y = log_range_size, fill = plot_group)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "gray55") +
  geom_boxplot(
    width = 0.58,
    linewidth = 0.7,
    outlier.alpha = 0.22,
    outlier.size = 1.1,
    color = "gray24"
  ) +
  geom_text(
    data = summary_data,
    aes(x = plot_group, y = y_text, label = stats_label),
    inherit.aes = FALSE,
    size = 3.7,
    lineheight = 0.92,
    color = "gray16"
  ) +
  facet_wrap(~ threat, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = plot_colors, guide = "none") +
  scale_x_discrete(labels = axis_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.16))) +
  labs(
    title = "Amphibian range size and wood-production threat status",
    x = NULL,
    y = expression("Scaled log"[10] * " range size"),
    caption = "Dashed line = dataset mean. Threat flags are evaluated separately; species can occur in both threatened groups."
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15.5, hjust = 0),
    axis.text.x = element_text(size = 11.5, color = "gray20"),
    axis.text.y = element_text(size = 11.5, color = "gray20"),
    axis.title.y = element_text(size = 13, margin = margin(r = 9)),
    strip.background = element_rect(fill = "gray92", color = "gray45"),
    strip.text = element_text(face = "bold", size = 12.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9.2, color = "gray35"),
    plot.margin = margin(10, 14, 10, 10)
  )

png_path <- file.path(figure_dir, "amphibian_log_range_size_by_logging_or_timber_plantation.png")
pdf_path <- file.path(figure_dir, "amphibian_log_range_size_by_logging_or_timber_plantation.pdf")
summary_path <- file.path(figure_dir, "amphibian_log_range_size_by_logging_or_timber_plantation_summary.csv")

ggsave(png_path, range_plot, width = 9.2, height = 5.8, dpi = 300)
ggsave(pdf_path, range_plot, width = 9.2, height = 5.8)
write.csv(
  summary_data %>%
    mutate(threat = as.character(threat), status = as.character(status)) %>%
    select(threat, status, n, mean, median, q25, q75, min, max),
  summary_path,
  row.names = FALSE
)

message("Saved amphibian range-size threat plot to ", png_path)
message("Saved amphibian range-size threat plot to ", pdf_path)
message("Saved summary table to ", summary_path)

body_size_breaks <- quantile(
  amphibian_data$body_size_log,
  probs = c(0, 1 / 3, 2 / 3, 1),
  na.rm = TRUE
)

if (anyDuplicated(body_size_breaks) > 0) {
  stop("Cannot create body-size tertiles because Body_Size_log quantile breaks are duplicated.", call. = FALSE)
}

body_size_labels <- c(
  paste0("Small\n", format_number(body_size_breaks[1]), " to ", format_number(body_size_breaks[2])),
  paste0("Medium\n", format_number(body_size_breaks[2]), " to ", format_number(body_size_breaks[3])),
  paste0("Large\n", format_number(body_size_breaks[3]), " to ", format_number(body_size_breaks[4]))
)

binned_base <- amphibian_data %>%
  filter(is.finite(body_size_log)) %>%
  mutate(
    body_size_bin = cut(
      body_size_log,
      breaks = body_size_breaks,
      labels = body_size_labels,
      include.lowest = TRUE
    )
  )

plot_data_binned <- bind_rows(
  binned_base %>%
    transmute(
      species,
      body_size_bin,
      threat = "Logging",
      status = if_else(logging, "Threatened", "Not threatened"),
      log_range_size
    ),
  binned_base %>%
    transmute(
      species,
      body_size_bin,
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
  group_by(threat, body_size_bin, status, fill_group) %>%
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

body_size_bar_plot <- ggplot(
  binned_summary,
  aes(x = body_size_bin, y = mean, fill = fill_group)
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
    title = "Threatened vs non-threatened range size by amphibian body-size bin",
    x = expression("Body-size bin (tertiles of " * log[10] * " body size)"),
    y = expression("Mean scaled log"[10] * " range size"),
    caption = "Bars show means; error bars are 95% CIs. Body-size bins use the model predictor Body_Size_log."
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

bar_png_path <- file.path(figure_dir, "amphibian_log_range_size_by_threat_and_body_size_bins_barplot.png")
bar_pdf_path <- file.path(figure_dir, "amphibian_log_range_size_by_threat_and_body_size_bins_barplot.pdf")
bar_summary_path <- file.path(figure_dir, "amphibian_log_range_size_by_threat_and_body_size_bins_summary.csv")

ggsave(bar_png_path, body_size_bar_plot, width = 10.4, height = 6.2, dpi = 300)
ggsave(bar_pdf_path, body_size_bar_plot, width = 10.4, height = 6.2)
write.csv(
  binned_summary %>%
    mutate(
      threat = as.character(threat),
      body_size_bin = as.character(body_size_bin),
      status = as.character(status)
    ) %>%
    select(threat, body_size_bin, status, n, mean, median, sd, se, ci95, q25, q75, min, max),
  bar_summary_path,
  row.names = FALSE
)

message("Saved body-size binned bar plot to ", bar_png_path)
message("Saved body-size binned bar plot to ", bar_pdf_path)
message("Saved body-size binned summary table to ", bar_summary_path)

scatter_data <- bind_rows(
  amphibian_data %>%
    transmute(
      species,
      body_size_log,
      log_range_size,
      threat = "Logging",
      status = if_else(logging, "Threatened", "Not threatened")
    ),
  amphibian_data %>%
    transmute(
      species,
      body_size_log,
      log_range_size,
      threat = "Timber plantations",
      status = if_else(timber_plantations, "Threatened", "Not threatened")
    )
) %>%
  filter(is.finite(body_size_log), is.finite(log_range_size)) %>%
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
  aes(x = body_size_log, y = log_range_size, color = color_group)
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
    title = "Body size and range size by amphibian threat status",
    x = expression(log[10] * " body size"),
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

scatter_png_path <- file.path(figure_dir, "amphibian_body_size_range_size_by_threat_scatter.png")
scatter_pdf_path <- file.path(figure_dir, "amphibian_body_size_range_size_by_threat_scatter.pdf")

ggsave(scatter_png_path, interaction_scatter, width = 10.4, height = 6.2, dpi = 300)
ggsave(scatter_pdf_path, interaction_scatter, width = 10.4, height = 6.2)

message("Saved body-size/range-size threat scatter plot to ", scatter_png_path)
message("Saved body-size/range-size threat scatter plot to ", scatter_pdf_path)
