rm(list = ls())
gc()

# =============================================================================
# COMPARE HOTSPOTS: RICHNESS (main text) VS. CUBE-ROOT RANGE-WEIGHTED
#
# Description:
# Quantifies how much the cube-root inverse-range weighting changes which
# cells are classified as top-10% wood-production-threat hotspots, relative
# to the main-text richness (simple species count) based hotspots. Compares:
#   1. Per-taxon top-10% hotspot masks (Jaccard index, % agreement, cells
#      gained/lost)
#   2. The cross-taxon hotspot-overlap-frequency maps (how many taxa agree a
#      cell is a hotspot), correlation and category-shift table
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c("terra", "sf", "rnaturalearth", "dplyr", "ggplot2", "here")
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

raster_dir     <- here::here("data", "processed", "hotspots", "rasters")
ses_raster_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "rasters")
out_dir        <- here::here("outputs", "hotspots", "figures", "supplementary_transformed")
table_dir      <- here::here("outputs", "hotspots", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

top10_mask <- function(r, world_map, filter_positive = TRUE) {
  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE)
  if (filter_positive) {
    vals <- vals[vals > 0 & !is.na(vals)]
  } else {
    vals <- vals[!is.na(vals)]
  }
  if (length(vals) == 0) return(NULL)
  threshold <- quantile(vals, probs = 0.90, na.rm = TRUE)
  binary <- (r_mask >= threshold)
  binary[is.na(binary)] <- 0
  binary
}

# --- 1. Per-taxon comparison --------------------------------------------------
cat("=== Per-taxon top-10% hotspot agreement: Richness vs Cube-Root Weighted ===\n\n")

per_taxon_stats <- data.frame()

for (taxon in taxa_list) {
  path_rich  <- file.path(raster_dir, paste0(taxon, "_combined_richness_raster.tif"))
  path_trans <- file.path(raster_dir, paste0(taxon, "_combined_transformed_threat_raster.tif"))
  if (!file.exists(path_rich) || !file.exists(path_trans)) {
    warning("Missing raster(s) for ", taxon, " -- skipping.")
    next
  }

  r_rich  <- rast(path_rich)
  r_trans <- rast(path_trans)

  mask_rich  <- top10_mask(r_rich, world_map)
  mask_trans <- top10_mask(r_trans, world_map)

  v_rich  <- values(mask_rich, mat = FALSE)
  v_trans <- values(mask_trans, mat = FALSE)
  valid <- !is.na(v_rich) & !is.na(v_trans)
  v_rich  <- v_rich[valid]
  v_trans <- v_trans[valid]

  both      <- sum(v_rich == 1 & v_trans == 1)
  rich_only <- sum(v_rich == 1 & v_trans == 0)
  trans_only <- sum(v_rich == 0 & v_trans == 1)
  union_n   <- both + rich_only + trans_only
  jaccard   <- if (union_n > 0) both / union_n else NA

  n_rich  <- sum(v_rich == 1)
  n_trans <- sum(v_trans == 1)

  cat(sprintf(
    "%-10s | richness hotspots: %4d | weighted hotspots: %4d | shared: %4d | Jaccard: %.3f | %% richness cells reclassified: %.1f%%\n",
    taxon, n_rich, n_trans, both, jaccard, 100 * rich_only / n_rich
  ))

  per_taxon_stats <- rbind(per_taxon_stats, data.frame(
    taxon = taxon, n_richness_hotspot_cells = n_rich, n_weighted_hotspot_cells = n_trans,
    n_shared = both, n_richness_only = rich_only, n_weighted_only = trans_only,
    jaccard = jaccard,
    pct_richness_hotspots_lost = 100 * rich_only / n_rich,
    pct_weighted_hotspots_new  = 100 * trans_only / n_trans
  ))
}

write.csv(per_taxon_stats, file.path(table_dir, "hotspot_weighting_comparison_per_taxon.csv"), row.names = FALSE)
cat("\nSaved:", file.path(table_dir, "hotspot_weighting_comparison_per_taxon.csv"), "\n\n")

# --- 2. Cross-taxon overlap-frequency comparison ------------------------------
cat("=== Cross-taxon hotspot-overlap-frequency: Richness vs Cube-Root Weighted ===\n\n")

freq_richness <- rast(grid_template, vals = 0)
freq_weighted <- rast(grid_template, vals = 0)

for (taxon in taxa_list) {
  r_rich  <- rast(file.path(raster_dir, paste0(taxon, "_combined_richness_raster.tif")))
  r_trans <- rast(file.path(raster_dir, paste0(taxon, "_combined_transformed_threat_raster.tif")))

  mask_rich  <- top10_mask(r_rich, world_map)
  mask_trans <- top10_mask(r_trans, world_map)

  freq_richness <- freq_richness + mask_rich
  freq_weighted <- freq_weighted + mask_trans
}

v_fr <- values(terra::mask(freq_richness, vect(world_map)), mat = FALSE)
v_fw <- values(terra::mask(freq_weighted, vect(world_map)), mat = FALSE)
valid <- !is.na(v_fr) & !is.na(v_fw)
v_fr <- v_fr[valid]; v_fw <- v_fw[valid]

corr <- cor(v_fr, v_fw, method = "spearman")
cat("Spearman correlation (richness freq vs weighted freq):", round(corr, 3), "\n\n")

# Category-shift table (0-5 taxa hotspot count, richness rows vs weighted cols)
shift_table <- table(
  Richness_taxa_count = factor(v_fr, levels = 0:5),
  Weighted_taxa_count = factor(v_fw, levels = 0:5)
)
print(shift_table)
write.csv(as.data.frame(shift_table), file.path(table_dir, "hotspot_weighting_crosstaxon_shift_table.csv"), row.names = FALSE)

# Cells that are 5-taxon hotspots under one method but not the other
five_rich_only <- sum(v_fr == 5 & v_fw < 5)
five_weighted_only <- sum(v_fw == 5 & v_fr < 5)
five_both <- sum(v_fr == 5 & v_fw == 5)
cat(sprintf(
  "\n5-taxon hotspot cells -- richness only: %d | weighted only: %d | both: %d\n",
  five_rich_only, five_weighted_only, five_both
))

# --- 3. Visual: side-by-side agreement map for the combined (all-taxa) overlap
diff_df <- data.frame(
  x = terra::xyFromCell(grid_template, which(valid))[, 1],
  y = terra::xyFromCell(grid_template, which(valid))[, 2],
  richness = v_fr, weighted = v_fw
) %>%
  mutate(category = case_when(
    richness == 0 & weighted == 0 ~ NA_character_,
    richness > 0 & weighted == 0 ~ "Hotspot lost (richness only)",
    richness == 0 & weighted > 0 ~ "Hotspot gained (weighted only)",
    weighted > richness ~ "Intensified under weighting",
    weighted < richness ~ "Weakened under weighting",
    TRUE ~ "Unchanged"
  )) %>%
  filter(!is.na(category))

cat_colors <- c(
  "Hotspot lost (richness only)"    = "#2166AC",
  "Weakened under weighting"        = "#92C5DE",
  "Unchanged"                       = "#B0B0B0",
  "Intensified under weighting"     = "#F4A582",
  "Hotspot gained (weighted only)"  = "#B2182B"
)
diff_df$category <- factor(diff_df$category, levels = names(cat_colors))

p_diff <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = diff_df, aes(x = x, y = y, fill = category)) +
  scale_fill_manual(values = cat_colors, name = "Change in taxa-hotspot count\n(cube-root weighted vs richness)", drop = FALSE) +
  coord_sf(crs = moll_crs, expand = FALSE) +
  theme_void() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "hotspot_weighting_comparison_map.png"),
       p_diff, width = 10, height = 5, dpi = 300, bg = "white")

cat("\nSaved comparison map:", file.path(out_dir, "hotspot_weighting_comparison_map.png"), "\n")
cat("\nDone.\n")
