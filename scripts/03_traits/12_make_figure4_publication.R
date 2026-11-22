# =============================================================================
#
# FIGURE 4 — PUBLICATION VERSION (Biological Conservation)
#
# Description:
# Trait-model coefficient forest plots for the five taxa, rebuilt at final
# print size from the coefficient table written by
# 07_make_combined_taxa_coefficient_plot.R (outputs/traits/figures/
# all_taxa_coefficients.csv: posterior medians and 95% credible intervals,
# with the trait ordering used there). Re-run script 07 first if the models
# change.
#
# Layout: left column birds (A) over mammals (C); right column amphibians (B),
# fish (D), reptiles (E); shared production-pathway legend at the bottom.
# Taxon silhouettes sit in an empty corner of each panel.
#
# Output (outputs/publication_figures/Figure_4/):
#   - Figure_4.pdf   vector, fonts embedded (cairo_pdf)   <- submit this
#   - Figure_4.tiff  1000-dpi LZW-compressed bitmap (line-art requirement)
#   - Figure_4_preview.png  300-dpi preview for checking only
#
# Sizing follows the Elsevier artwork guide: full-page width = 190 mm; all
# text is 7–8 pt at final size (Arial).
#
# =============================================================================


# --- 1. SETUP ---
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(here)
})

coefficient_path <- here::here("outputs", "traits", "figures", "all_taxa_coefficients.csv")
silhouette_dir   <- here::here("archive", "reference_material", "species_silhouettes")
out_dir          <- here::here("outputs", "publication_figures", "Figure_4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fig_width_mm <- 190

font_family <- "Arial"
size_axis   <- 7    # trait labels, axis numbers, axis titles
size_legend <- 8    # legend title and labels
size_tag    <- 10   # panel tags

threat_colors <- c(ThreatLogging = "#d9be00", ThreatPlantation = "#0088d9")
threat_labels <- c(ThreatLogging = "Natural Forestry", ThreatPlantation = "Plantation")

# Panel tag, silhouette file and silhouette corner for each taxon
panels <- tibble::tribble(
  ~taxon,        ~tag, ~silhouette,          ~corner,
  "Birds",       "A",  "bird sil.png",       "top-left",
  "Amphibians",  "B",  "amphibian sil.png",  "bottom-right",
  "Mammals",     "C",  "mammal sil.png",     "bottom-right",
  "Fish",        "D",  "fish sil.png",       "bottom-right",
  "Reptiles",    "E",  "reptile sil.png",    "bottom-right"
)


# --- 2. DATA ---
# -----------------------------------------------------------------------------
plot_summary <- read.csv(coefficient_path, stringsAsFactors = FALSE) %>%
  mutate(
    Response = factor(Response, levels = names(threat_colors)),
    # Natural forestry drawn just below plantation within each trait row
    y_plot = y + ifelse(Response == "ThreatLogging", -0.18, 0.18)
  )

trait_counts <- plot_summary %>%
  distinct(Taxon, y) %>%
  count(Taxon, name = "n_traits") %>%
  tibble::deframe()


# --- 3. SILHOUETTES ---
# -----------------------------------------------------------------------------
# Recoloured solid black and cropped to their visible extent
read_silhouette <- function(file) {
  img <- png::readPNG(file.path(silhouette_dir, file))
  alpha <- img[, , dim(img)[3]]
  keep_rows <- range(which(rowSums(alpha > 0.05) > 0))
  keep_cols <- range(which(colSums(alpha > 0.05) > 0))
  alpha <- alpha[keep_rows[1]:keep_rows[2], keep_cols[1]:keep_cols[2]]
  out <- array(0, dim = c(dim(alpha), 4))
  out[, , 4] <- alpha
  out
}

# Silhouette fitted in a box (mm), pinned to a panel corner with a small inset
sil_box_mm <- c(w = 11, h = 9)
sil_inset_mm <- 1.5
silhouette_layer <- function(file, corner) {
  img <- read_silhouette(file)
  a <- ncol(img) / nrow(img)
  box_w <- sil_box_mm[["w"]]; box_h <- sil_box_mm[["h"]]
  if (a > box_w / box_h) { w <- box_w; h <- box_w / a } else { h <- box_h; w <- box_h * a }

  right <- corner == "bottom-right"
  grob <- grid::rasterGrob(
    img,
    x = if (right) unit(1, "npc") - unit(sil_inset_mm, "mm") else unit(sil_inset_mm, "mm"),
    y = if (right) unit(sil_inset_mm, "mm") else unit(1, "npc") - unit(sil_inset_mm, "mm"),
    width = unit(w, "mm"), height = unit(h, "mm"),
    just = if (right) c("right", "bottom") else c("left", "top")
  )
  annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}


# --- 4. PANELS ---
# -----------------------------------------------------------------------------
make_taxon_plot <- function(taxon_label) {
  info <- panels %>% filter(taxon == taxon_label)
  taxon_summary <- plot_summary %>% filter(Taxon == taxon_label)
  taxon_order <- taxon_summary %>% distinct(y, Clean_Label) %>% arrange(y)

  x_limits <- range(c(taxon_summary$lower, taxon_summary$upper), na.rm = TRUE)
  x_limits <- x_limits + c(-1, 1) * diff(x_limits) * 0.05

  ggplot(taxon_summary, aes(colour = Response)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_segment(
      aes(x = lower, xend = upper, y = y_plot, yend = y_plot),
      linewidth = 0.55, lineend = "round"
    ) +
    geom_point(aes(x = estimate, y = y_plot), size = 1.1) +
    silhouette_layer(info$silhouette, info$corner) +
    scale_colour_manual(values = threat_colors, labels = threat_labels, name = "Production Pathway") +
    scale_y_continuous(
      breaks = taxon_order$y, labels = taxon_order$Clean_Label,
      limits = c(0.4, max(taxon_order$y) + 0.6), expand = expansion(mult = 0)
    ) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(x = "Effect on risk (log-odds)", y = NULL, tag = info$tag) +
    theme_bw(base_size = size_axis, base_family = font_family) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(size = size_axis, colour = "black"),
      axis.title.x = element_text(size = size_axis, margin = margin(t = 1.5, unit = "mm")),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(0.8, "mm"),
      panel.border = element_rect(linewidth = 0.4, colour = "black"),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_line(colour = "grey92", linetype = "dotted", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.tag = element_text(size = size_tag, face = "bold", family = font_family),
      plot.tag.position = c(0, 1),
      plot.margin = margin(t = 1, r = 2.5, b = 1.5, l = 5, unit = "mm")
    )
}

plots <- setNames(lapply(panels$taxon, make_taxon_plot), panels$taxon)


# --- 5. ASSEMBLE ---
# -----------------------------------------------------------------------------
# Each column is its own stack with panel heights proportional to the number
# of traits, so trait rows are evenly spaced within a column and both columns
# end at the same height.
left_column <- plots$Birds / plots$Mammals +
  plot_layout(heights = trait_counts[c("Birds", "Mammals")])
right_column <- plots$Amphibians / plots$Fish / plots$Reptiles +
  plot_layout(heights = trait_counts[c("Amphibians", "Fish", "Reptiles")])

row_mm        <- 4.1   # printed height of one trait row in the left column
axis_block_mm <- 9     # x-axis text + title + margins per panel
legend_mm     <- 8
columns_mm    <- sum(trait_counts[c("Birds", "Mammals")]) * row_mm + 2 * axis_block_mm
fig_height_mm <- columns_mm + legend_mm

figure_4 <- (left_column | right_column) / guide_area() +
  plot_layout(heights = c(columns_mm, legend_mm), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = size_legend, face = "bold", family = font_family),
    legend.text = element_text(size = size_legend, family = font_family),
    legend.key.width = unit(7, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.key.spacing.x = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0)
  )


# --- 6. SAVE ---
# -----------------------------------------------------------------------------
out_base <- file.path(out_dir, "Figure_4")

ggsave(paste0(out_base, ".pdf"), figure_4,
       width = fig_width_mm, height = fig_height_mm, units = "mm", device = cairo_pdf)
ggsave(paste0(out_base, ".tiff"), figure_4,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 1000, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(out_base, "_preview.png"), figure_4,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")

cat(sprintf("Saved Figure_4 (%.0f x %.0f mm) to %s\n", fig_width_mm, fig_height_mm, out_dir))
