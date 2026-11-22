rm(list = ls())
gc()

# =============================================================================
# COMPARE RICHNESS VS CUBE-ROOT WEIGHTED: WHOLE-DISTRIBUTION COMPARISON
#
# Description:
# Extends 16_compare_hotspot_weighting_schemes.R (which only compared a binary
# top-10% hotspot cutoff) to the full distribution:
#   1. Spearman rank correlation across ALL nonzero cells, per taxon
#   2. Percentile-rank scatter (richness rank vs weighted rank), all cells
#   3. Continuous percentile-rank-shift maps (per taxon + cross-taxon mean) --
#      shows how every cell's relative importance changes, not just whether
#      it crosses the 90th-percentile hotspot line
#   4. Jaccard overlap curve across multiple thresholds (1/5/10/25/50%)
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c("terra", "sf", "rnaturalearth", "dplyr", "ggplot2", "here", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000, crs = moll_crs,
  extent = terra::ext(-18040000, 18040000, -9020000, 9020000), vals = 0
)

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% st_transform(moll_crs)

taxa_list <- c("amphibian", "bird", "fish", "mammal", "reptile")
raster_dir <- here::here("data", "processed", "hotspots", "rasters")
out_dir    <- here::here("outputs", "hotspots", "figures", "supplementary_transformed")
table_dir  <- here::here("outputs", "hotspots", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load and align both metrics per taxon, keep percentile ranks ------------
all_cell_data <- list()
rank_shift_stack <- list()

for (taxon in taxa_list) {
  r_rich  <- rast(file.path(raster_dir, paste0(taxon, "_combined_richness_raster.tif")))
  r_trans <- rast(file.path(raster_dir, paste0(taxon, "_combined_transformed_threat_raster.tif")))

  s <- c(r_rich, r_trans); names(s) <- c("richness", "weighted")
  s_masked <- terra::mask(s, vect(world_map))
  df <- as.data.frame(s_masked, xy = TRUE) %>% filter(richness > 0 | weighted > 0)

  df$richness_pctile <- rank(df$richness, ties.method = "average") / nrow(df)
  df$weighted_pctile <- rank(df$weighted, ties.method = "average") / nrow(df)
  df$rank_shift <- df$weighted_pctile - df$richness_pctile
  df$taxon <- taxon

  all_cell_data[[taxon]] <- df

  # rasterize the rank shift for cross-taxon aggregation
  shift_r <- rast(grid_template, vals = NA_real_)
  cells <- terra::cellFromXY(grid_template, as.matrix(df[, c("x", "y")]))
  shift_r[cells] <- df$rank_shift
  rank_shift_stack[[taxon]] <- shift_r
}

combined_df <- bind_rows(all_cell_data)

# --- 1. Full-distribution Spearman correlation --------------------------------
cat("=== Full-distribution Spearman correlation (richness vs weighted) ===\n\n")
corr_table <- combined_df %>%
  group_by(taxon) %>%
  summarise(
    n_cells = n(),
    spearman_r = cor(richness, weighted, method = "spearman"),
    .groups = "drop"
  )
print(corr_table)
write.csv(corr_table, file.path(table_dir, "hotspot_weighting_full_distribution_correlation.csv"), row.names = FALSE)

overall_r <- cor(combined_df$richness_pctile, combined_df$weighted_pctile, method = "spearman")
cat("\nPooled across all taxa (percentile ranks): Spearman r =", round(overall_r, 3), "\n\n")

# --- 2. Percentile-rank scatter, faceted by taxon -----------------------------
p_scatter <- ggplot(combined_df, aes(x = richness_pctile, y = weighted_pctile)) +
  geom_hex(bins = 60) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  scale_fill_viridis_c(option = "magma", trans = "log10", name = "Cell count") +
  facet_wrap(~ tools::toTitleCase(taxon), ncol = 3) +
  labs(
    x = "Percentile rank -- Richness (species count)",
    y = "Percentile rank -- Cube-root weighted threat"
  ) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "hotspot_weighting_rank_scatter.png"),
       p_scatter, width = 12, height = 9, dpi = 300, bg = "white")

# --- 3. Continuous rank-shift maps --------------------------------------------
cat("Building continuous rank-shift maps...\n")

shift_colors <- scale_fill_gradient2(
  low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
  name = "Percentile-rank shift\n(weighted - richness)",
  limits = c(-1, 1)
)

for (taxon in taxa_list) {
  df <- all_cell_data[[taxon]]
  p <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = rank_shift)) +
    shift_colors +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() +
    theme(legend.position = "right")
  ggsave(file.path(out_dir, paste0(taxon, "_hotspot_weighting_rank_shift.png")),
         p, width = 10, height = 5, dpi = 300, bg = "white")
}

# Cross-taxon mean rank shift (only where at least one taxon has data)
mean_shift <- mean(rast(rank_shift_stack), na.rm = TRUE)
mean_shift_df <- as.data.frame(terra::mask(mean_shift, vect(world_map)), xy = TRUE) %>%
  setNames(c("x", "y", "mean_rank_shift")) %>% filter(!is.na(mean_rank_shift))

p_mean_shift <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = mean_shift_df, aes(x = x, y = y, fill = mean_rank_shift)) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "Mean percentile-rank shift\nacross taxa (weighted - richness)"
  ) +
  coord_sf(crs = moll_crs, expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "all_taxa_mean_hotspot_weighting_rank_shift.png"),
       p_mean_shift, width = 10, height = 5, dpi = 300, bg = "white")

# --- 4. Jaccard overlap across multiple thresholds ----------------------------
cat("\n=== Jaccard overlap at multiple thresholds ===\n\n")

thresholds <- c(0.01, 0.05, 0.10, 0.25, 0.50)
jaccard_curve <- expand.grid(taxon = taxa_list, threshold = thresholds) %>%
  rowwise() %>%
  mutate(
    jaccard = {
      df <- all_cell_data[[taxon]]
      rich_top <- df$richness_pctile >= (1 - threshold)
      wt_top   <- df$weighted_pctile >= (1 - threshold)
      inter <- sum(rich_top & wt_top)
      uni   <- sum(rich_top | wt_top)
      if (uni > 0) inter / uni else NA_real_
    }
  ) %>%
  ungroup()

print(jaccard_curve %>% tidyr::pivot_wider(names_from = threshold, values_from = jaccard))
write.csv(jaccard_curve, file.path(table_dir, "hotspot_weighting_jaccard_by_threshold.csv"), row.names = FALSE)

jaccard_curve$taxon_label <- tools::toTitleCase(as.character(jaccard_curve$taxon))

p_jaccard <- ggplot(jaccard_curve, aes(x = threshold, y = jaccard, color = taxon_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(labels = scales::percent, breaks = thresholds) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Top-X%% threshold used to define \"hotspot\"", y = "Jaccard overlap (richness vs weighted)",
       color = "Taxon") +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, "hotspot_weighting_jaccard_by_threshold.png"),
       p_jaccard, width = 8, height = 6, dpi = 300, bg = "white")

cat("\nAll whole-distribution comparison outputs saved to:\n  ", out_dir, "\n  ", table_dir, "\n")
