# =============================================================================
# SCRIPT: SES HOTSPOT OVERLAP FREQUENCY
#
# Description:
# SES analogue of all_vertebrates_hotspot_overlap_frequency.png. Instead of
# counting how many taxa have raw threatened-richness in their top 10%, this
# counts how many taxa have SES(threatened richness) [null-model standardized
# effect size for has_wood_prod_ongoing richness] in their top 10%.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required_packages <- c("dplyr", "ggplot2", "here", "sf", "terra", "rnaturalearth")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% st_transform(moll_crs)

out_dir <- here::here("outputs", "hotspots", "figures", "maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ses_dir <- here::here("data", "processed", "hotspots", "null_model", "rasters")
taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

hotspot_freq_ses <- NULL

for (taxon in taxa_list) {
  cat("Processing SES overlap for:", taxon, "\n")
  path_ses <- file.path(ses_dir, paste0(taxon, "_has_wood_prod_ongoing_ses_raster.tif"))
  if (!file.exists(path_ses)) stop("Missing SES raster for ", taxon, ": ", path_ses)

  r <- rast(path_ses)
  if (is.null(hotspot_freq_ses)) hotspot_freq_ses <- rast(r, vals = 0)

  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE); vals <- vals[!is.na(vals)]
  if (length(vals) == 0) next

  threshold_10 <- quantile(vals, probs = 0.90, na.rm = TRUE)
  r_binary <- (r_mask >= threshold_10)
  r_binary[is.na(r_binary)] <- 0
  hotspot_freq_ses <- hotspot_freq_ses + r_binary
}

freq_df <- as.data.frame(terra::mask(hotspot_freq_ses, vect(world_map)), xy = TRUE) %>%
  setNames(c("x", "y", "count")) %>% filter(count > 0) %>%
  mutate(count = factor(count, levels = 1:5))

freq_colors <- c("1" = "#FED976", "2" = "#FEB24C", "3" = "#FD8D3C", "4" = "#F03B20", "5" = "#BD0026")

p_overlap_ses <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count") +
  theme_void() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "all_vertebrates_ses_hotspot_overlap_frequency.png"),
       p_overlap_ses, width = 10, height = 6, bg = "white")

cat("\nSES Hotspot Overlap Map Saved.\n")
