rm(list = ls())
gc()

# =============================================================================
# SUPPLEMENTARY FIGURE: CUBE-ROOT RANGE-WEIGHTED THREAT MAPS
#
# Description:
# Reproduces the main-text map formatting (gradient threat maps, continuous
# SES maps, 4x4 biscale bivariate maps, and the raw-vs-SES hotspot overlap
# panel) using the cube-root inverse-range-weighted "Transformed Threat"
# rasters instead of raw threatened-species richness. For each taxon:
#   1. Wood-production threat, continuous gradient   (plot_gradient_sqrt style)
#   2. SES wood-production threat, continuous         (plot_continuous_ses_map style)
#   3. Bivariate: Natural Forest (logging) vs Plantation threat (biscale 4x4)
#   4. SES bivariate: Natural Forest (logging) vs Plantation SES (biscale 4x4)
# Plus one cross-taxon panel: top-10% hotspot overlap frequency (1-5 taxa) for
# the weighted threatened-richness raster and its SES equivalent, side by side.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c(
  "dplyr", "ggplot2", "here", "sf", "terra", "rnaturalearth",
  "stringr", "biscale", "scales", "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

# --- 1. Global constants -----------------------------------------------------
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
map_legend_title_size <- 16
map_legend_text_size  <- 12

taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% st_transform(moll_crs)

raster_dir          <- here::here("data", "processed", "hotspots", "rasters")
ses_raster_dir       <- here::here("data", "processed", "hotspots", "null_model_weighted", "rasters")
out_dir <- here::here("outputs", "hotspots", "figures", "supplementary_transformed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

trans_raster_path <- function(taxon, threat) file.path(raster_dir, paste0(taxon, "_", threat, "_transformed_threat_raster.tif"))
ses_raster_path    <- function(taxon, flag)   file.path(ses_raster_dir, paste0(taxon, "_", flag, "_transformed_ses_raster.tif"))

# --- 2. Plotting engines (ported from the main-text map scripts) -------------

# 2a. Continuous gradient (sqrt-distributed colour, linear legend breaks)
#     -- matches plot_gradient_sqrt() in 02_make_taxon_maps.R
plot_gradient_sqrt <- function(r, legend_title = "Threat") {
  r_mask <- terra::mask(r, vect(world_map))
  df <- as.data.frame(r_mask, xy = TRUE) %>%
    setNames(c("x", "y", "val")) %>%
    filter(!is.na(val))
  if (nrow(df) == 0) return(NULL)

  df_zero <- df %>% filter(val == 0)
  df_pos  <- df %>% filter(val > 0)
  if (nrow(df_pos) == 0) return(NULL)

  max_val <- max(df_pos$val, na.rm = TRUE)
  if (!is.finite(max_val)) return(NULL)
  scale_max <- if (max_val > 0) max_val else 1

  p_breaks <- pretty(c(0, scale_max), n = 5)
  step_size <- if (length(p_breaks) > 1) p_breaks[2] - p_breaks[1] else max_val
  kept_breaks <- p_breaks[p_breaks < max_val]
  if (length(kept_breaks) > 0) {
    last_b <- tail(kept_breaks, 1)
    if ((max_val - last_b) < (0.25 * step_size)) kept_breaks <- head(kept_breaks, -1)
  }
  breaks_linear <- if (max_val > 0) c(kept_breaks, max_val) else 0

  n_cols <- 256
  colors_pal <- scales::viridis_pal(option = "magma", direction = 1)(n_cols)
  values_dist <- seq(0, 1, length.out = n_cols)^2

  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df_zero, aes(x = x, y = y), fill = "black") +
    geom_tile(data = df_pos, aes(x = x, y = y, fill = val)) +
    scale_fill_gradientn(
      colors = colors_pal, values = values_dist, breaks = breaks_linear,
      name = legend_title, limits = c(0, scale_max), labels = round(breaks_linear, 2)
    ) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.25, "cm"),
      legend.title = element_text(size = map_legend_title_size),
      legend.text = element_text(size = map_legend_text_size)
    )
}

# 2b. Continuous symmetric divergent SES map -- matches
#     plot_continuous_ses_map() in 08_bivariate_SES_mapping.R
plot_continuous_ses_map <- function(ses_raster) {
  r_masked <- terra::mask(ses_raster, terra::vect(world_map))
  df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "value")
  if (nrow(df) == 0) return(NULL)

  df_zero <- df %>% filter(value == 0)
  df_nz   <- df %>% filter(value != 0)

  max_abs <- if (nrow(df_nz) > 0) max(abs(df_nz$value), na.rm = TRUE) else 1
  if (max_abs == 0) max_abs <- 1

  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df_zero, aes(x = x, y = y), fill = "black") +
    geom_tile(data = df_nz, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs, max_abs), name = "SES"
    ) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.25, "cm"),
      legend.title = element_text(size = map_legend_title_size),
      legend.text = element_text(size = map_legend_text_size)
    )
}

# 2c. Continuous 2D bivariate color-mix -- matches plot_true_continuous_bivariate()
#     in 08_bivariate_SES_mapping.R exactly (same corner colors, no special
#     zero handling; the black-for-zero treatment is only used on the
#     single-variable gradient maps above, not the bivariate ones).
bivar_c00 <- "#e8e8e8"  # low X, low Y (neutral gray)
bivar_c10 <- "#d9be00"  # high X, low Y
bivar_c01 <- "#0088d9"  # low X, high Y
bivar_c11 <- "#007b00"  # high X, high Y

mix_2d_colors <- function(x_vec, y_vec) {
  rgb00 <- col2rgb(bivar_c00)[, 1]; rgb10 <- col2rgb(bivar_c10)[, 1]
  rgb01 <- col2rgb(bivar_c01)[, 1]; rgb11 <- col2rgb(bivar_c11)[, 1]

  r_x0 <- (1 - x_vec) * rgb00[1] + x_vec * rgb10[1]
  g_x0 <- (1 - x_vec) * rgb00[2] + x_vec * rgb10[2]
  b_x0 <- (1 - x_vec) * rgb00[3] + x_vec * rgb10[3]

  r_x1 <- (1 - x_vec) * rgb01[1] + x_vec * rgb11[1]
  g_x1 <- (1 - x_vec) * rgb01[2] + x_vec * rgb11[2]
  b_x1 <- (1 - x_vec) * rgb01[3] + x_vec * rgb11[3]

  r_final <- (1 - y_vec) * r_x0 + y_vec * r_x1
  g_final <- (1 - y_vec) * g_x0 + y_vec * g_x1
  b_final <- (1 - y_vec) * b_x0 + y_vec * b_x1

  rgb(r_final, g_final, b_final, maxColorValue = 255)
}

plot_continuous_bivariate <- function(r_x, r_y) {
  s <- c(r_x, r_y); names(s) <- c("X", "Y")
  s_masked <- terra::mask(s, vect(world_map))
  df <- as.data.frame(s_masked, xy = TRUE) %>% filter(!is.na(X) & !is.na(Y))
  if (nrow(df) == 0) return(NULL)

  max_x <- max(df$X, na.rm = TRUE); max_y <- max(df$Y, na.rm = TRUE)
  if (max_x <= 0) max_x <- 1
  if (max_y <= 0) max_y <- 1

  df <- df %>% mutate(
    X_norm = pmax(X, 0) / max_x,
    Y_norm = pmax(Y, 0) / max_y
  )
  df$mix_color <- mix_2d_colors(df$X_norm, df$Y_norm)

  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = mix_color)) +
    scale_fill_identity() +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void()
}

# 2e. Standalone shared bivariate legend (saved once, not per taxon/map) --
#     matches plot_standalone_bivariate_legend() in 08_bivariate_SES_mapping.R
plot_standalone_bivariate_legend <- function(output_path) {
  leg_grid <- expand.grid(X_norm = seq(0, 1, length.out = 300), Y_norm = seq(0, 1, length.out = 300))
  leg_grid$leg_color <- mix_2d_colors(leg_grid$X_norm, leg_grid$Y_norm)

  p_legend <- ggplot(leg_grid) +
    geom_tile(aes(x = X_norm, y = Y_norm, fill = leg_color)) +
    scale_fill_identity() +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(
      title = "Inverse Cube Root Weighted Threat",
      x = "Natural Forestry",
      y = "Plantation"
    ) +
    coord_fixed() +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 8)),
      axis.line = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(15, 15, 15, 15, "pt")
    )

  ggsave(output_path, p_legend, width = 4.5, height = 4.5, dpi = 300, bg = "white")
}

# --- 3. Per-taxon maps --------------------------------------------------------
for (taxon in taxa_list) {
  cat("\nProcessing (transformed/weighted):", toupper(taxon), "\n")

  path_comb  <- trans_raster_path(taxon, "combined")
  path_log   <- trans_raster_path(taxon, "logging")
  path_plant <- trans_raster_path(taxon, "plantation")
  if (!all(file.exists(c(path_comb, path_log, path_plant)))) {
    warning("Missing transformed threat raster(s) for ", taxon, " -- skipping.")
    next
  }
  r_comb  <- rast(path_comb)
  r_log   <- rast(path_log)
  r_plant <- rast(path_plant)

  # 1. Continuous gradient: wood production threat
  p1 <- plot_gradient_sqrt(r_comb, legend_title = "Weighted Threat\n(Cube-Root)")
  if (!is.null(p1)) {
    ggsave(file.path(out_dir, paste0(taxon, "_transformed_gradient_wood_production.png")),
           p1, width = 10, height = 5, dpi = 300, bg = "white")
  }

  # 2. SES wood production threat (continuous, symmetric divergent)
  path_ses_comb <- ses_raster_path(taxon, "has_wood_prod_ongoing")
  if (file.exists(path_ses_comb)) {
    r_ses_comb <- rast(path_ses_comb)
    p2 <- plot_continuous_ses_map(r_ses_comb)
    if (!is.null(p2)) {
      ggsave(file.path(out_dir, paste0(taxon, "_transformed_ses_wood_production.png")),
             p2, width = 10, height = 5, dpi = 300, bg = "white")
    }
  } else {
    warning("Missing weighted SES raster (combined) for ", taxon, " -- run 13_weighted_null_model_ses.R first.")
  }

  # 3. Bivariate: Natural Forest (logging) vs Plantation threat (continuous)
  p3 <- plot_continuous_bivariate(r_log, r_plant)
  if (!is.null(p3)) {
    ggsave(file.path(out_dir, paste0(taxon, "_transformed_bivariate_logging_vs_plantation.png")),
           p3, width = 10, height = 5, dpi = 300, bg = "white")
  }

  # 4. SES bivariate: Natural Forest (logging) vs Plantation SES (continuous)
  path_ses_log   <- ses_raster_path(taxon, "has_logging_ongoing")
  path_ses_plant <- ses_raster_path(taxon, "has_plantation_ongoing")
  if (file.exists(path_ses_log) && file.exists(path_ses_plant)) {
    r_ses_log   <- rast(path_ses_log)
    r_ses_plant <- rast(path_ses_plant)
    p4 <- plot_continuous_bivariate(r_ses_log, r_ses_plant)
    if (!is.null(p4)) {
      ggsave(file.path(out_dir, paste0(taxon, "_transformed_ses_bivariate_logging_vs_plantation.png")),
             p4, width = 10, height = 5, dpi = 300, bg = "white")
    }
  } else {
    warning("Missing weighted SES raster(s) (logging/plantation) for ", taxon,
            " -- run 13_weighted_null_model_ses.R first.")
  }
}

# --- 4. Cross-taxon overlap panel: raw weighted vs SES weighted --------------
# Mirrors 11_combined_overlap_frequency_panel.R (side-by-side, shared legend),
# but built from the cube-root weighted combined threat raster and its SES.
cat("\nBuilding cross-taxon hotspot overlap panel (weighted)...\n")

freq_colors <- c("1" = "#FED976", "2" = "#FEB24C", "3" = "#FD8D3C", "4" = "#F03B20", "5" = "#BD0026")

build_overlap_freq_df <- function(raster_paths, world_map, filter_positive) {
  hotspot_freq <- NULL
  for (path in raster_paths) {
    if (!file.exists(path)) stop("Missing raster: ", path)
    r <- rast(path)
    if (is.null(hotspot_freq)) hotspot_freq <- rast(r, vals = 0)

    r_mask <- terra::mask(r, vect(world_map))
    vals <- values(r_mask, mat = FALSE)
    if (filter_positive) {
      vals <- vals[vals > 0 & !is.na(vals)]
    } else {
      vals <- vals[!is.na(vals)]
    }
    if (length(vals) == 0) next

    threshold_10 <- quantile(vals, probs = 0.90, na.rm = TRUE)
    r_binary <- (r_mask >= threshold_10)
    r_binary[is.na(r_binary)] <- 0
    hotspot_freq <- hotspot_freq + r_binary
  }

  as.data.frame(terra::mask(hotspot_freq, vect(world_map)), xy = TRUE) %>%
    setNames(c("x", "y", "count")) %>% filter(count > 0) %>%
    mutate(count = factor(count, levels = 1:5))
}

raw_paths <- trans_raster_path(taxa_list, "combined")
raw_freq_df <- build_overlap_freq_df(raw_paths, world_map, filter_positive = TRUE)

ses_paths <- ses_raster_path(taxa_list, "has_wood_prod_ongoing")
ses_freq_df <- build_overlap_freq_df(ses_paths, world_map, filter_positive = FALSE)

p_raw <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = raw_freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count", drop = FALSE) +
  theme_void()

p_ses <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = ses_freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count", drop = FALSE) +
  theme_void()

combined <- (p_raw + p_ses) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

ggsave(
  file.path(out_dir, "all_vertebrates_transformed_hotspot_overlap_frequency_combined.png"),
  combined, width = 16, height = 7, dpi = 300, bg = "white"
)

# --- 5. Standalone shared bivariate legend (one asset, used by all bivariate maps) ---
plot_standalone_bivariate_legend(file.path(out_dir, "bivariate_legend.png"))

cat("\nAll cube-root weighted supplementary maps saved to:", out_dir, "\n")
