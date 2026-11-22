# =============================================================================
#
# FIGURE S5 — PUBLICATION VERSION (Biological Conservation)
#
#
# Description:
# Per-cell comparison of species richness threatened by natural forestry
# (logging) vs plantations, one square panel per taxon (A amphibians,
# B birds, C freshwater fishes, D mammals, E reptiles) with the realm legend
# in the sixth slot. Rebuilt at print size from create_source_scatter() in
# 02_make_taxon_maps.R: same rasters, land mask, realm assignment (realm
# polygons rasterised with touches = TRUE, gaps filled by modal focal
# windows), cell filter (at least one threatened species) and 1:1 line.
#
# Also prints the percentage of cells below the 1:1 line (natural forestry
# threatens more species than plantations) per taxon and for all vertebrates
# combined, for the figure caption.
#
# Output (outputs/publication_figures/Figure_S5/):
#   - Figure_S5.pdf   vector, fonts embedded (cairo_pdf)   <- submit this
#   - Figure_S5.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_S5_preview.png  300-dpi preview for checking only
#
# =============================================================================

source(here::here("scripts", "02_hotspots", "publication_map_helpers.R"))

taxon_tags <- c(amphibian = "A", bird = "B", fish = "C", mammal = "D", reptile = "E")

realm_colors <- c(
  "Afrotropical"   = "#FF6B6B", "Australian"     = "#D4A017",
  "Madagascan"     = "#FFD54F", "Nearctic"       = "#55A630",
  "Neotropical"    = "#00C853", "Oceanina"       = "#00BFA5",
  "Oriental"       = "#00B0FF", "Palearctic"     = "#2979FF",
  "Panamanian"     = "#6200EA", "Saharo-Arabian" = "#D500F9",
  "Sino-Japanese"  = "#C51162"
)


# --- 1. REALMS ---
# -----------------------------------------------------------------------------
suppressMessages(sf::sf_use_s2(FALSE))
realms_sf <- st_read(here::here("data", "raw", "realms", "newRealms.shp"), quiet = TRUE) %>%
  st_make_valid() %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>%
  st_make_valid()
potential_cols <- names(realms_sf)[sapply(realms_sf, function(x) length(unique(x)) > 4 & length(unique(x)) < 20)]
realm_col_name <- grep("Realm", potential_cols, value = TRUE, ignore.case = TRUE)[1]
if (is.na(realm_col_name)) realm_col_name <- potential_cols[1]

realms_proj <- st_transform(realms_sf, moll_crs)
realms_proj$rid <- as.numeric(as.factor(realms_proj[[realm_col_name]]))
realm_lookup <- unique(st_drop_geometry(realms_proj)[, c("rid", realm_col_name)])
names(realm_lookup) <- c("realm_id", "Realm")

# Realm id per land cell, with unassigned land cells filled from their
# neighbours (as create_source_scatter())
realm_raster <- function(template_land) {
  realm_r <- terra::rasterize(vect(realms_proj), template_land, field = "rid", touches = TRUE)
  target <- !is.na(template_land)
  for (i in 1:30) {
    missing <- target & is.na(realm_r)
    if (global(missing, "sum", na.rm = TRUE)[1, 1] == 0) break
    win_size <- if (i <= 5) 3 else if (i <= 15) 5 else 7
    filled <- terra::focal(realm_r, w = win_size, fun = "modal", na.policy = "only", na.rm = TRUE)
    realm_r <- terra::cover(realm_r, filled)
  }
  realm_r
}


# --- 2. DATA ---
# -----------------------------------------------------------------------------
richness_raster <- function(taxon, pathway) {
  rast(file.path(raster_dir, paste0(taxon, "_", pathway, "_richness_raster.tif")))
}

source_scatter_df <- function(r_log, r_plant) {
  x_land <- terra::mask(r_log, world_vect)
  y_land <- terra::mask(r_plant, world_vect)
  stack_data <- c(x_land, y_land, realm_raster(x_land))
  names(stack_data) <- c("X", "Y", "realm_id")
  as.data.frame(stack_data, xy = FALSE) %>%
    filter(!is.na(X) & !is.na(Y) & (X > 0 | Y > 0)) %>%
    left_join(realm_lookup, by = "realm_id") %>%
    filter(!is.na(Realm))
}

scatter_data <- lapply(taxa_order, function(taxon) {
  cat("Figure S5:", taxon, "\n")
  source_scatter_df(richness_raster(taxon, "logging"), richness_raster(taxon, "plantation"))
})
names(scatter_data) <- taxa_order

# Caption statistics: share of plotted cells below the 1:1 line
all_log   <- Reduce(`+`, lapply(taxa_order, richness_raster, pathway = "logging"))
all_plant <- Reduce(`+`, lapply(taxa_order, richness_raster, pathway = "plantation"))
all_df    <- source_scatter_df(all_log, all_plant)
below_line <- c(sapply(scatter_data, function(df) mean(df$X > df$Y)), all_vertebrates = mean(all_df$X > all_df$Y))
cat("\nCells where natural forestry threatens more species than plantations:\n")
print(round(100 * below_line, 1))

caption_stats <- tibble::tibble(
  taxon = names(below_line),
  n_cells_plotted = c(sapply(scatter_data, nrow), nrow(all_df)),
  pct_natural_forestry_gt_plantation = round(100 * below_line, 1)
)
dir.create(file.path(pub_dir, "Figure_S5"), recursive = TRUE, showWarnings = FALSE)
write.csv(caption_stats, file.path(pub_dir, "Figure_S5", "Figure_S5_caption_statistics.csv"), row.names = FALSE)


# --- 3. PANELS ---
# -----------------------------------------------------------------------------
panel_mm <- 66   # square panel size

# Silhouette in the top-left corner of the panel (inside the axes)
scatter_silhouette <- function(taxon) {
  img <- silhouettes[[taxon]]
  a <- ncol(img) / nrow(img)
  box_w <- 13; box_h <- 10
  if (a > box_w / box_h) { w <- box_w; h <- box_w / a } else { h <- box_h; w <- box_h * a }
  annotation_custom(
    grid::rasterGrob(
      img, x = unit(3, "mm"), y = unit(1, "npc") - unit(1.5, "mm"),
      width = unit(w, "mm"), height = unit(h, "mm"), just = c("left", "top")
    ),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  )
}

make_scatter <- function(taxon) {
  df <- scatter_data[[taxon]]
  max_val <- max(c(df$X, df$Y))

  ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(colour = Realm), size = 0.55, stroke = 0, shape = 16) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3) +
    scatter_silhouette(taxon) +
    scale_x_continuous(limits = c(0, max_val), expand = expansion(mult = 0.04)) +
    scale_y_continuous(limits = c(0, max_val), expand = expansion(mult = 0.04)) +
    scale_colour_manual(values = realm_colors, name = "Zoogeographic Realm", drop = FALSE,
                        limits = names(realm_colors)) +
    labs(x = "Natural Forestry Threat Richness", y = "Plantation Threat Richness", tag = taxon_tags[[taxon]]) +
    theme_classic(base_size = 7, base_family = font_family) +
    theme(
      aspect.ratio = 1,
      text = element_text(colour = "black"),
      axis.text = element_text(size = 7, colour = "black"),
      axis.title = element_text(size = 7),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      plot.tag = element_text(size = size_tag, face = "bold", family = font_family),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(3.5, "mm"),
      plot.margin = margin(1.5, 3, 1.5, 1.5, unit = "mm")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.5)))
}

plots <- lapply(taxa_order, make_scatter)


# --- 4. ASSEMBLE & SAVE ---
# -----------------------------------------------------------------------------
figure_s5 <- wrap_plots(c(plots, list(guide_area())), ncol = 2) +
  plot_layout(
    widths = unit(rep(panel_mm, 2), "mm"),
    heights = unit(rep(panel_mm, 3), "mm"),
    guides = "collect"
  )

save_publication_figure(figure_s5, "Figure_S5", 3 * (panel_mm + 12))
