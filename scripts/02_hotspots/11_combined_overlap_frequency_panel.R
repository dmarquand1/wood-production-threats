# =============================================================================
# SCRIPT: COMBINED HOTSPOT OVERLAP FREQUENCY PANEL (RAW vs SES)
#
# Description:
# Side-by-side panel: raw threatened-richness top-10% overlap frequency (left)
# vs SES(threatened richness) top-10% overlap frequency (right), sharing one
# horizontal "Taxa Count" legend centred underneath both maps.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required_packages <- c("dplyr", "ggplot2", "here", "sf", "terra", "rnaturalearth", "patchwork")
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

raw_raster_dir <- here::here("data", "processed", "hotspots", "rasters")
ses_raster_dir <- here::here("data", "processed", "hotspots", "null_model", "rasters")
taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

freq_colors <- c("1" = "#FED976", "2" = "#FEB24C", "3" = "#FD8D3C", "4" = "#F03B20", "5" = "#BD0026")

# --- Helper: build a Taxa-Count overlap-frequency data frame from a list of rasters ---
build_overlap_freq_df <- function(raster_paths, world_map) {
  hotspot_freq <- NULL
  for (path in raster_paths) {
    if (!file.exists(path)) stop("Missing raster: ", path)
    r <- rast(path)
    if (is.null(hotspot_freq)) hotspot_freq <- rast(r, vals = 0)

    r_mask <- terra::mask(r, vect(world_map))
    vals <- values(r_mask, mat = FALSE); vals <- vals[!is.na(vals)]
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

# --- Raw threatened-richness overlap (values >= 0 only, matches original script) ---
cat("Building RAW overlap frequency layer...\n")
raw_paths <- file.path(raw_raster_dir, paste0(taxa_list, "_combined_richness_raster.tif"))
raw_hotspot_freq <- NULL
for (path in raw_paths) {
  r <- rast(path)
  if (is.null(raw_hotspot_freq)) raw_hotspot_freq <- rast(r, vals = 0)
  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE); vals <- vals[vals > 0 & !is.na(vals)]
  if (length(vals) == 0) next
  threshold_10 <- quantile(vals, probs = 0.90, na.rm = TRUE)
  r_binary <- (r_mask >= threshold_10)
  r_binary[is.na(r_binary)] <- 0
  raw_hotspot_freq <- raw_hotspot_freq + r_binary
}
raw_freq_df <- as.data.frame(terra::mask(raw_hotspot_freq, vect(world_map)), xy = TRUE) %>%
  setNames(c("x", "y", "count")) %>% filter(count > 0) %>%
  mutate(count = factor(count, levels = 1:5))

# --- SES overlap (all valid cells, matches 10_ses_hotspot_overlap.R) ---
cat("Building SES overlap frequency layer...\n")
ses_paths <- file.path(ses_raster_dir, paste0(taxa_list, "_has_wood_prod_ongoing_ses_raster.tif"))
ses_freq_df <- build_overlap_freq_df(ses_paths, world_map)

# --- Plots (no individual legends; collected into one shared legend) ---
p_raw <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = raw_freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count", drop = FALSE) +
  labs(title = "Threatened Richness") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

p_ses <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = ses_freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count", drop = FALSE) +
  labs(title = "SES(Threatened Richness)") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

combined <- (p_raw + p_ses) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

ggsave(
  file.path(out_dir, "all_vertebrates_hotspot_overlap_frequency_combined.png"),
  combined, width = 16, height = 7, dpi = 300, bg = "white"
)

cat("\nCombined Overlap Frequency Panel Saved.\n")
