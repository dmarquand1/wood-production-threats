# =============================================================================
#
# FIGURES S2 & S4 — PUBLICATION VERSIONS (Biological Conservation)
#
#
# Description:
# Cross-taxon overlap of wood-production threat hotspots, rebuilt at print
# size with the shared map formatting (publication_map_helpers.R). Hotspot
# definition follows 11_combined_overlap_frequency_panel.R (and section 4 of
# 14_make_transformed_supplementary_maps.R): for each taxon, a cell is a
# hotspot if it is in the top 10% of that taxon's values on land; cells are
# coloured by how many of the five taxa they are a hotspot for.
#
#   Figure S2 (cube-root range-weighted):
#     A = weighted threat score (positive values only enter the quantile),
#         from the same matrix-derived scores as Figure S1
#     B = SES of the weighted threat score (all values enter the quantile)
#   Figure S4 (threatened species richness):
#     A = threatened richness (positive values only enter the quantile)
#     B = SES of threatened richness (all values enter the quantile)
#
# Output (outputs/publication_figures/Figure_S2/ and Figure_S4/):
#   - Figure_SN.pdf   fonts embedded (cairo_pdf)   <- submit this
#   - Figure_SN.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_SN_preview.png  300-dpi preview for checking only
#
# =============================================================================

source(here::here("scripts", "02_hotspots", "publication_map_helpers.R"))

freq_colors <- c("1" = "#FED976", "2" = "#FEB24C", "3" = "#FD8D3C", "4" = "#F03B20", "5" = "#BD0026")

# Number of taxa for which each land cell is in that taxon's top 10%
hotspot_frequency_df <- function(rasters, positive_only) {
  counts <- NULL
  for (r in rasters) {
    r_land <- terra::mask(r, world_vect)
    vals <- values(r_land, mat = FALSE)
    vals <- vals[!is.na(vals) & (!positive_only | vals > 0)]
    threshold <- quantile(vals, probs = 0.90, na.rm = TRUE)

    is_hotspot <- r_land >= threshold
    is_hotspot[is.na(is_hotspot)] <- 0
    counts <- if (is.null(counts)) is_hotspot else counts + is_hotspot
  }
  raster_df(counts, "count") %>%
    filter(count > 0) %>%
    mutate(count = factor(count, levels = 1:5))
}

plot_hotspot_frequency <- function(freq_df, tag) {
  ggplot() +
    map_base() +
    geom_raster(data = freq_df, aes(x = x, y = y, fill = count)) +
    scale_fill_manual(values = freq_colors, name = "Taxa Count", drop = FALSE) +
    labs(tag = tag) +
    base_map_theme
}

# Two maps side by side across the page (A = raw metric, B = its SES) with a
# shared horizontal legend beneath
save_overlap_figure <- function(raw_rasters, ses_rasters, name) {
  cat(name, ": hotspots\n", sep = "")
  raw_freq <- hotspot_frequency_df(raw_rasters, positive_only = TRUE)
  ses_freq <- hotspot_frequency_df(ses_rasters, positive_only = FALSE)

  map_w    <- (fig_width_mm - 8) / 2
  map_h    <- map_w / map_aspect
  legend_h <- 7

  figure <- wrap_plots(
    A = plot_hotspot_frequency(raw_freq, "A"),
    B = plot_hotspot_frequency(ses_freq, "B"),
    C = guide_area(),
    design = "AB\nCC"
  ) +
    plot_layout(
      widths = unit(rep(map_w, 2), "mm"),
      heights = unit(c(map_h, legend_h), "mm"),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = size_title, margin = margin(r = 2, unit = "mm")),
      legend.text = element_text(size = size_body, margin = margin(l = 1, unit = "mm")),
      legend.key.size = unit(3, "mm"),
      legend.key.spacing.x = unit(2.5, "mm")
    )

  save_publication_figure(figure, name, map_h + legend_h + 8)
}

# Figure S2: cube-root range-weighted threat
save_overlap_figure(
  lapply(taxa_order, weighted_threat_raster),
  lapply(taxa_order, weighted_ses_raster),
  "Figure_S2"
)

# Figure S4: threatened species richness
save_overlap_figure(
  lapply(taxa_order, function(taxon) rast(file.path(raster_dir, paste0(taxon, "_combined_richness_raster.tif")))),
  lapply(taxa_order, function(taxon) rast(file.path(ses_dir, paste0(taxon, "_has_wood_prod_ongoing_ses_raster.tif")))),
  "Figure_S4"
)
