# =============================================================================
#
# FIGURE S7 — PUBLICATION VERSION (Biological Conservation)
#
# Description:
# Log-odds of an ongoing natural forestry (5.3) vs plantation (2.2) threat by
# publication period of each species' current assessment, per class and
# pooled, rebuilt at print size from the table written by
# 05_threat_proportions_over_time.R (outputs/extent/tables/
# threat_proportions_over_time.csv: Haldane-Anscombe corrected log-odds and
# 95% CIs). Re-run script 05 first if the assessments change. Styling
# follows the log-odds plot in script 05; the end-of-line value labels are
# placed directly (no repel) and nudged apart only if they would touch.
#
# Output (outputs/publication_figures/Figure_S7/):
#   - Figure_S7.pdf   vector, fonts embedded (cairo_pdf)   <- submit this
#   - Figure_S7.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_S7_preview.png  300-dpi preview for checking only
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(here)
})

out_dir <- here::here("outputs", "publication_figures", "Figure_S7")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_width_mm  <- 190
fig_height_mm <- 118
font_family   <- "Arial"

summary_levels <- c("All Vertebrates", "Amphibia", "Aves", "Freshwater Fish", "Mammalia", "Reptilia")
period_levels  <- c("2010 or earlier", "2011-2015", "2016-2020", "2021-2025")
threat_colors  <- c("Natural Forestry" = "#eda100", "Plantation" = "#2a78d6")

period_summary <- read.csv(here::here("outputs", "extent", "tables", "threat_proportions_over_time.csv"))

plot_data <- bind_rows(
  period_summary %>% transmute(
    summary_group, assessment_period, threat_type = "Natural Forestry",
    log_odds = logging_log_odds, ci_lower = logging_log_odds_lower, ci_upper = logging_log_odds_upper
  ),
  period_summary %>% transmute(
    summary_group, assessment_period, threat_type = "Plantation",
    log_odds = plantation_log_odds, ci_lower = plantation_log_odds_lower, ci_upper = plantation_log_odds_upper
  )
) %>%
  mutate(
    summary_group = factor(summary_group, levels = summary_levels),
    assessment_period = factor(assessment_period, levels = period_levels),
    threat_type = factor(threat_type, levels = names(threat_colors))
  )

# End-of-line labels: latest period per panel and pathway. If the two labels
# in a panel would touch, push them apart symmetrically (in data units,
# relative to the panel's y range).
label_data <- plot_data %>%
  group_by(summary_group, threat_type) %>%
  filter(as.integer(assessment_period) == max(as.integer(assessment_period))) %>%
  ungroup() %>%
  group_by(summary_group) %>%
  mutate(
    y_range = diff(range(c(plot_data$ci_lower[plot_data$summary_group == first(summary_group)],
                           plot_data$ci_upper[plot_data$summary_group == first(summary_group)]))),
    min_gap = 0.09 * y_range,
    gap = abs(diff(range(log_odds))),
    label_y = if (n() == 2 && first(gap) < first(min_gap)) {
      mid <- mean(log_odds)
      mid + sign(log_odds - mid) * first(min_gap) / 2
    } else log_odds
  ) %>%
  ungroup()

dodge <- position_dodge(width = 0.15)

figure_s7 <- ggplot(plot_data, aes(x = assessment_period, y = log_odds, colour = threat_type, group = threat_type)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15, linewidth = 0.35, position = dodge) +
  geom_line(linewidth = 0.55, position = dodge) +
  geom_point(size = 1.3, position = dodge) +
  geom_text(
    data = label_data, aes(y = label_y, label = sprintf("%.2f", log_odds)),
    nudge_x = 0.18, hjust = 0, size = 6 / .pt, family = font_family, show.legend = FALSE
  ) +
  facet_wrap(~ summary_group, nrow = 2, scales = "free_y") +
  scale_colour_manual(name = "Production Pathway", values = threat_colors) +
  scale_y_continuous(name = "Log odds of threat, ln(p / (1 − p))", expand = expansion(mult = c(0.08, 0.12))) +
  scale_x_discrete(name = "Publication year of current assessment", expand = expansion(add = c(0.45, 0.75))) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 7, base_family = font_family) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(size = 6.5, colour = "black"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    axis.title = element_text(size = 7.5),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    strip.text = element_text(size = 7.5, face = "bold"),
    strip.background = element_rect(linewidth = 0.3),
    panel.spacing.x = unit(4, "mm"),
    panel.spacing.y = unit(3, "mm"),
    legend.position = "bottom",
    legend.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.width = unit(6, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(1, 2, 1, 1, unit = "mm")
  )

out_base <- file.path(out_dir, "Figure_S7")
ggsave(paste0(out_base, ".pdf"), figure_s7,
       width = fig_width_mm, height = fig_height_mm, units = "mm", device = cairo_pdf)
ggsave(paste0(out_base, ".tiff"), figure_s7,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 1000, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(out_base, "_preview.png"), figure_s7,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")
cat(sprintf("Saved Figure_S7 (%.0f x %.0f mm) to %s\n", fig_width_mm, fig_height_mm, out_dir))
