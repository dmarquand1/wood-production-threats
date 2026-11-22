# =============================================================================
#
# SHARED HELPERS FOR PUBLICATION MAP FIGURES
#
# Sourced by 19_make_figures_2_3_publication.R (Figures 2 & 3) and
# 20_make_figure_s1_publication.R (Figure S1) so the main-text and
# supplementary maps use identical formatting: Mollweide land background,
# magma sqrt-distributed gradient maps, symmetric blue–white–red SES maps,
# colour-bar legends, silhouettes, and the 190 mm full-page output set
# (PDF with embedded fonts, 1000-dpi TIFF, 300-dpi preview).
#
# =============================================================================

# --- 1. SETUP ---
# -----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)
library(patchwork)
library(Matrix)
library(here)

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
# Same 110 km grid as the matrices were built on (08_bivariate_SES_mapping.R)
grid_template <- rast(
  resolution = 110000, crs = moll_crs,
  extent = ext(-18040000, 18040000, -9020000, 9020000), vals = 0
)

taxa_order <- c("amphibian", "bird", "fish", "mammal", "reptile")
silhouette_files <- c(
  amphibian = "amphibian sil.png",
  bird      = "bird sil.png",
  fish      = "fish sil.png",
  mammal    = "mammal sil.png",
  reptile   = "reptile sil.png"
)

raster_dir     <- here::here("data", "processed", "hotspots", "rasters")
ses_dir        <- here::here("data", "processed", "hotspots", "null_model", "rasters")
ses_pathway_dir <- here::here("data", "processed", "hotspots", "ses_rasters")
matrix_dir     <- here::here("data", "processed", "hotspots", "null_model", "matrices")
silhouette_dir <- here::here("archive", "reference_material", "species_silhouettes")
pub_dir        <- here::here("outputs", "publication_figures")

fig_width_mm <- 190

# Typography (pt at final print size); geom_text sizes are mm, so divide by .pt
font_family <- "Arial"
pt <- function(x) x / ggplot2::.pt
size_body  <- 6    # legend tick labels
size_title <- 7    # legend titles
size_tag   <- 10   # panel tags

land_fill <- "#f0f0f0"

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>%
  st_transform(moll_crs)
world_vect <- vect(world_map)

# Map extent (land only, as coord_sf(expand = FALSE) drew it before)
map_bbox   <- st_bbox(world_map)
map_aspect <- unname((map_bbox["xmax"] - map_bbox["xmin"]) / (map_bbox["ymax"] - map_bbox["ymin"]))


# --- 2. SHARED HELPERS ---
# -----------------------------------------------------------------------------
base_map_theme <- theme_void() +
  theme(
    text = element_text(family = font_family, colour = "black"),
    plot.tag = element_text(size = size_tag, face = "bold", family = font_family, hjust = 0, vjust = 1),
    plot.tag.position = c(0, 1),
    legend.title = element_text(size = size_title, lineheight = 0.9, margin = margin(b = 1, unit = "mm")),
    legend.text = element_text(size = size_body, margin = margin(l = 0.8, unit = "mm")),
    legend.key.width = unit(2.2, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(0, 0, 0, 0)
  )


# Colour bar with an explicit length (ggplot2 otherwise makes bars 5x the
# key height, which overflows the short map rows) and no tick marks
colourbar_length_mm <- 20
min_label_gap_mm    <- 3    # 6 pt labels closer than this collide
colourbar_guide <- function() {
  guide_colourbar(theme = theme(
    legend.key.height = unit(colourbar_length_mm, "mm"),
    legend.key.width = unit(2.2, "mm"),
    legend.ticks = element_blank()
  ))
}
# Land background plus coordinate system shared by every map
map_base <- function() {
  list(
    geom_sf(data = world_map, fill = land_fill, colour = NA),
    coord_sf(
      crs = moll_crs, expand = FALSE,
      xlim = map_bbox[c("xmin", "xmax")], ylim = map_bbox[c("ymin", "ymax")]
    )
  )
}

# Raster cells overlapping land, as an x/y data frame (terra::mask as before)
raster_df <- function(r, value_name = "val") {
  df <- as.data.frame(terra::mask(r, world_vect), xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", value_name)
  mutate(df, x = round(x), y = round(y))  # remove float noise so geom_raster sees a regular grid
}

# Silhouettes: recoloured solid black, cropped to their visible extent
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
silhouettes <- lapply(silhouette_files, read_silhouette)

# Silhouette in the top-right corner of a map panel (empty ocean in the
# Mollweide projection). The box is sized in mm and converted to map units
# from the panel's printed size, so the animal keeps its aspect ratio.
sil_box_mm <- c(w = 11, h = 8)
silhouette_layer <- function(taxon, panel_w_mm, panel_h_mm) {
  img <- silhouettes[[taxon]]
  a <- ncol(img) / nrow(img)
  box_w <- sil_box_mm[["w"]]; box_h <- sil_box_mm[["h"]]
  if (a > box_w / box_h) { w <- box_w; h <- box_w / a } else { h <- box_h; w <- box_h * a }
  x_per_mm <- (map_bbox[["xmax"]] - map_bbox[["xmin"]]) / panel_w_mm
  y_per_mm <- (map_bbox[["ymax"]] - map_bbox[["ymin"]]) / panel_h_mm
  x_right <- map_bbox[["xmax"]] - (0.5 + (box_w - w) / 2) * x_per_mm
  y_top   <- map_bbox[["ymax"]] - 0.5 * y_per_mm
  annotation_raster(
    img,
    xmin = x_right - w * x_per_mm, xmax = x_right,
    ymin = y_top - h * y_per_mm, ymax = y_top,
    interpolate = TRUE
  )
}

save_publication_figure <- function(plot, name, height_mm) {
  # Ocean columns with no land cells leave gaps in x that are exact multiples
  # of the 110 km cell width; geom_raster warns about "uneven intervals" but
  # still places every cell correctly, so that one warning is muffled.
  withCallingHandlers(
    save_publication_files(plot, name, height_mm),
    warning = function(w) {
      if (grepl("uneven horizontal intervals", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

save_publication_files <- function(plot, name, height_mm) {
  out_dir <- file.path(pub_dir, name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_base <- file.path(out_dir, name)

  ggsave(paste0(out_base, ".pdf"), plot,
         width = fig_width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(paste0(out_base, ".tiff"), plot,
         width = fig_width_mm, height = height_mm, units = "mm",
         dpi = 1000, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  ggsave(paste0(out_base, "_preview.png"), plot,
         width = fig_width_mm, height = height_mm, units = "mm",
         dpi = 300, device = ragg::agg_png, bg = "white")

  cat(sprintf("Saved %s (%.0f x %.0f mm) to %s\n", name, fig_width_mm, height_mm, out_dir))
}


# --- 3. MAP ENGINES ---
# -----------------------------------------------------------------------------

# Gradient map, magma with sqrt-distributed colours and linear legend breaks
# (as plot_gradient_sqrt() in 02_make_taxon_maps.R). `r` must be NA where the
# taxon has no species (drawn as land background); cells with species but a
# value of 0 are black. The maximum is always labelled (to `label_digits`
# decimals), and the neat break below it is dropped if the two would print
# too close together on the colour bar.
plot_gradient_map <- function(r, legend_title, label_digits = 0) {
  df <- raster_df(r) %>% mutate(val = ifelse(val == 0, NA, val))
  max_val <- max(df$val, na.rm = TRUE)

  p_breaks <- pretty(c(0, max_val), n = 5)
  kept_breaks <- p_breaks[p_breaks < max_val]
  gap_mm <- (max_val - tail(kept_breaks, 1)) / max_val * colourbar_length_mm
  if (gap_mm < min_label_gap_mm) kept_breaks <- head(kept_breaks, -1)
  breaks_linear <- c(kept_breaks, max_val)
  break_labels <- c(format(kept_breaks, trim = TRUE), format(round(max_val, label_digits), nsmall = label_digits))

  n_cols <- 256
  ggplot() +
    map_base() +
    geom_raster(data = df, aes(x = x, y = y, fill = val)) +
    scale_fill_gradientn(
      colours = scales::viridis_pal(option = "magma")(n_cols),
      values = seq(0, 1, length.out = n_cols)^2,
      limits = c(0, max_val), breaks = breaks_linear, labels = break_labels,
      na.value = "black",
      name = legend_title, guide = colourbar_guide()
    ) +
    base_map_theme
}

# SES map, symmetric blue–white–red (as plot_continuous_ses_map() in
# 08_bivariate_SES_mapping.R)
plot_ses_map <- function(r_ses) {
  df <- raster_df(r_ses)
  max_abs <- max(abs(df$val), na.rm = TRUE)

  ggplot() +
    map_base() +
    geom_raster(data = df, aes(x = x, y = y, fill = val)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs, max_abs), name = "SES", guide = colourbar_guide()
    ) +
    base_map_theme
}

# Five-row figure: per taxon, a gradient map (with silhouette) on the left and
# its SES map on the right, each with its own legend; tagged A–J.
# `left_plots` / `right_plots` are named lists of ggplots keyed by taxon.
# `map_w` is the map panel width (mm); the legends sit beside the maps, so it
# must leave room for both columns' legends within the 190 mm page (reduce it
# when the left legend title is wider).
pair_row_gap <- 1.5

save_map_pair_figure <- function(left_plots, right_plots, name, map_w = 83) {
  row_h <- map_w / map_aspect
  panels <- list()
  for (taxon in taxa_order) {
    panels <- c(panels, list(
      left_plots[[taxon]] + silhouette_layer(taxon, map_w, row_h),
      right_plots[[taxon]]
    ))
  }

  figure <- wrap_plots(panels, ncol = 2, byrow = TRUE) +
    plot_layout(
      widths = unit(rep(map_w, 2), "mm"),
      heights = unit(rep(row_h, length(taxa_order)), "mm")
    ) +
    plot_annotation(tag_levels = "A") &
    theme(plot.margin = margin(t = pair_row_gap / 2, b = pair_row_gap / 2, unit = "mm"))

  save_publication_figure(figure, name, length(taxa_order) * (row_h + pair_row_gap))
}


# --- 4. BIVARIATE NATURAL FORESTRY vs PLANTATION MAPS (Figures 3, S3) ---
# -----------------------------------------------------------------------------
# Colour mixing as mix_2d_colors() in 08_bivariate_SES_mapping.R
c00 <- "#e8e8e8"  # low x, low y
c10 <- "#d9be00"  # high x (natural forestry), low y
c01 <- "#0088d9"  # low x, high y (plantation)
c11 <- "#007b00"  # high x, high y

mix_2d_colors <- function(x_vec, y_vec) {
  corners <- lapply(list(c00, c10, c01, c11), function(col) col2rgb(col)[, 1])
  mixed <- sapply(1:3, function(k) {
    bottom <- (1 - x_vec) * corners[[1]][k] + x_vec * corners[[2]][k]
    top    <- (1 - x_vec) * corners[[3]][k] + x_vec * corners[[4]][k]
    (1 - y_vec) * bottom + y_vec * top
  })
  rgb(mixed[, 1], mixed[, 2], mixed[, 3], maxColorValue = 255)
}

# Each axis is scaled to its own maximum; negatives (SES) are floored at 0,
# so cells with SES <= 0 on both axes are grey
plot_bivariate <- function(r_x, r_y) {
  s <- c(r_x, r_y)
  names(s) <- c("X", "Y")
  df <- as.data.frame(terra::mask(s, world_vect), xy = TRUE) %>%
    filter(!is.na(X) & !is.na(Y)) %>%
    mutate(x = round(x), y = round(y))
  max_x <- max(max(df$X), 1e-12)
  max_y <- max(max(df$Y), 1e-12)
  df$mix_colour <- mix_2d_colors(pmax(df$X, 0) / max_x, pmax(df$Y, 0) / max_y)

  ggplot() +
    map_base() +
    geom_raster(data = df, aes(x = x, y = y, fill = mix_colour)) +
    scale_fill_identity() +
    base_map_theme
}

# Shared square legend
bivariate_legend <- function(title) {
  grid_df <- expand.grid(x = seq(0, 1, length.out = 200), y = seq(0, 1, length.out = 200))
  grid_df$colour <- mix_2d_colors(grid_df$x, grid_df$y)

  ggplot(grid_df) +
    geom_raster(aes(x = x, y = y, fill = colour)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(title = title, x = "Natural Forestry", y = "Plantation") +
    coord_fixed() +
    theme_minimal(base_family = font_family) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      plot.title = element_text(size = size_title, face = "bold", hjust = 0.5, lineheight = 0.9,
                                margin = margin(b = 1.5, unit = "mm")),
      axis.title.x = element_text(size = size_title, face = "bold", margin = margin(t = 1, unit = "mm")),
      axis.title.y = element_text(size = size_title, face = "bold", margin = margin(r = 1, unit = "mm")),
      plot.margin = margin(0, 0, 0, 0)
    )
}

# Five-row figure: per taxon, the raw bivariate map (with silhouette) on the
# left and the SES bivariate map on the right, tagged A–J, with the shared
# square legend in a narrow third column beside the last row.
# `raw` / `ses` are named lists (by taxon) of list(x = raster, y = raster).
save_bivariate_figure <- function(raw, ses, name, legend_title) {
  legend_w <- 30
  map_w    <- (fig_width_mm - legend_w) / 2
  row_h    <- map_w / map_aspect
  gap      <- 1.5

  panels <- list()
  for (taxon in taxa_order) {
    panels <- c(panels, list(
      plot_bivariate(raw[[taxon]]$x, raw[[taxon]]$y) + silhouette_layer(taxon, map_w, row_h),
      plot_bivariate(ses[[taxon]]$x, ses[[taxon]]$y)
    ))
  }
  names(panels) <- LETTERS[1:10]

  legend_panel <- wrap_elements(
    full = bivariate_legend(legend_title) +
      theme(plot.margin = margin(l = 3, r = 1, t = 1, b = 1, unit = "mm"))
  )

  figure <- wrap_plots(
    c(panels, list(K = legend_panel)),
    design = paste("AB#", "CD#", "EF#", "GH#", "IJK", sep = "\n")
  ) +
    plot_layout(
      widths = unit(c(map_w, map_w, legend_w), "mm"),
      heights = unit(rep(row_h, length(taxa_order)), "mm")
    ) +
    plot_annotation(tag_levels = list(c(LETTERS[1:10], ""))) &
    theme(plot.margin = margin(t = gap / 2, b = gap / 2, unit = "mm"))

  save_publication_figure(figure, name, length(taxa_order) * (row_h + gap) + 4)
}

# Sum of a cell x species matrix over the species flagged TRUE, as a raster;
# cells with no species of the taxon are NA
matrix_sum_raster <- function(matrix_obj, flag) {
  idx <- which(as.logical(matrix_obj$species[[flag]]) %in% TRUE)
  values_vec <- as.numeric(Matrix::rowSums(matrix_obj$matrix[, idx, drop = FALSE]))
  out <- rast(grid_template)
  all_values <- rep(NA_real_, ncell(out))
  all_values[matrix_obj$active_cells] <- values_vec
  values(out) <- all_values
  out
}


# --- 5. CUBE-ROOT RANGE-WEIGHTED DATA (Figures S1–S3) ---
# -----------------------------------------------------------------------------
# Weighted cell x species matrices (12_build_weighted_transformed_matrices.R,
# patched by 12d_recover_missing_species.R) and their SES rasters
# (13_weighted_null_model_ses.R).
weighted_matrix_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "matrices")
weighted_ses_dir    <- here::here("data", "processed", "hotspots", "null_model_weighted", "rasters")

read_weighted_matrix <- function(taxon) {
  readRDS(file.path(weighted_matrix_dir, paste0(taxon, "_cell_species_weighted_matrix.rds")))
}

# Weighted threat score per cell, T_c = sum over threatened species of
# (a_sc / A_s)^(1/3): the matrix entries summed over species with the given
# threat flag. Computed from the matrices rather than the cached
# *_transformed_threat_raster.tif files, which were built with
# terra::rasterize(cover = TRUE) and silently omit species whose whole range
# is smaller than one grid cell.
weighted_threat_raster <- function(taxon, flag = "has_wood_prod_ongoing") {
  matrix_sum_raster(read_weighted_matrix(taxon), flag)
}

weighted_ses_raster <- function(taxon, flag = "has_wood_prod_ongoing") {
  rast(file.path(weighted_ses_dir, paste0(taxon, "_", flag, "_transformed_ses_raster.tif")))
}
