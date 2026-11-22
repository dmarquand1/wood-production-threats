# =============================================================================
# SCRIPT: VERTEBRATE THREAT MAPPING - INDIVIDUAL MAPS GENERATOR
#
# Description:
# Generates individual high-resolution maps for Birds, Mammals, Amphibians, 
# Reptiles, Fish, and ALL VERTEBRATES combined.
#
# CONFIGURATION:
# - Bivariate Plot: 4x4 Grid, Legend Bottom-Left.
# - Scatter Plot: No 1:1 line for Richness vs Threat, Black LOESS, Aspect Ratio 1.
# - Gradient Maps: Square-root visual scaling with linear legend labels.
#   * UPDATE: Legend always shows Max Value, removes overlapping previous labels.
# - Proportion Maps: Wood-production threatened richness / total species richness.
# - Hotspot Overlap: Calculates how many taxa (0-5) are hotspots in a cell.
# =============================================================================

# --- 1. SETUP & LIBRARIES ---
local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "here", "sf", "terra", "rnaturalearth",
  "stringr", "rredlist", "biscale", "scales", "patchwork"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))
library(grid)

# --- 2. GLOBAL CONSTANTS ---
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent <- ext(-18040000, 18040000, -9020000, 9020000)
grid_template <- rast(resolution = 110000, crs = moll_crs, extent = global_extent, vals = 0)
scatter_axis_title_size <- 20
scatter_axis_text_size <- 16
map_legend_title_size <- 16
map_legend_text_size <- 12

# Realm Colors
my_colors <- c(
  "Afrotropical"   = "#FF6B6B", "Australian"     = "#D4A017",  
  "Madagascan"     = "#FFD54F", "Nearctic"       = "#55A630",
  "Neotropical"    = "#00C853", "Oceanina"       = "#00BFA5",  
  "Oriental"       = "#00B0FF", "Palearctic"     = "#2979FF",
  "Panamanian"     = "#6200EA", "Saharo-Arabian" = "#D500F9",  
  "Sino-Japanese"  = "#C51162"   
)

# --- 3. LOAD DATA ---
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% st_transform(moll_crs)

path_shp <- here::here("data", "raw", "realms", "newRealms.shp")
path_rds <- here::here("data", "realms_global_sf.rds")

if (file.exists(path_shp)) {
  suppressMessages(sf::sf_use_s2(FALSE))
  realms_sf <- st_read(path_shp, quiet = TRUE)
} else if (file.exists(path_rds)) {
  realms_sf <- readRDS(path_rds)
} else {
  realms_sf <- NULL
}

if (!is.null(realms_sf)) {
  realms_sf <- st_make_valid(realms_sf) %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>% st_make_valid()
  potential_cols <- names(realms_sf)[sapply(realms_sf, function(x) length(unique(x)) > 4 & length(unique(x)) < 20)]
  realm_col_name <- grep("Realm", potential_cols, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(realm_col_name)) realm_col_name <- potential_cols[1]
}

# Output Directory
out_dir <- here::here("outputs", "hotspots", "figures", "maps")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Canonical raster directory
raster_dir <- here::here("data", "processed", "hotspots", "rasters")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
raster_path <- function(filename) file.path(raster_dir, filename)

load_assessment_species_summary <- function(taxon_key) {
  readRDS(
    here::here(
      "data",
      "processed",
      "assessments",
      paste0(taxon_key, "_assessment_species_summary.rds")
    )
  )
}

# --- 4. HELPER FUNCTIONS ---

get_richness_raster <- function(species_names, range_sf, file_path) {
  if (file.exists(file_path)) return(rast(file_path))
  cat("  -> Calculating:", basename(file_path), "\n")
  out_rast <- rast(grid_template, vals = 0)
  subset_sf <- range_sf %>% filter(sci_name %in% species_names)
  if (nrow(subset_sf) > 0) {
    out_rast <- terra::rasterize(vect(st_transform(subset_sf, moll_crs)), 
                                 grid_template, field = 1, fun = "sum", background = 0)
  }
  writeRaster(out_rast, file_path, overwrite = TRUE)
  return(out_rast)
}

get_threatened_proportion_raster <- function(r_threatened, r_total, file_path) {
  cat("  -> Calculating:", basename(file_path), "\n")
  proportion_raster <- terra::ifel(r_total > 0, r_threatened / r_total, NA)
  writeRaster(proportion_raster, file_path, overwrite = TRUE)
  return(proportion_raster)
}

# Map: Continuous Gradient (Sqrt Scale, Smart Linear Legend)
plot_gradient_sqrt <- function(r, legend_title = "Richness", title_txt = "") {
  r_mask <- terra::mask(r, vect(world_map))
  df <- as.data.frame(r_mask, xy = TRUE) %>% setNames(c("x", "y", "val")) %>% filter(val > 0)
  if (nrow(df) == 0) return(NULL)
  
  # --- SMART BREAKS LOGIC ---
  max_val <- max(df$val, na.rm = TRUE)
  
  # 1. Generate standard neat breaks
  p_breaks <- pretty(c(0, max_val), n = 5)
  
  # 2. Determine typical step size
  if(length(p_breaks) > 1) {
    step_size <- p_breaks[2] - p_breaks[1]
  } else {
    step_size <- max_val
  }
  
  # 3. Filter breaks to those strictly less than Max
  kept_breaks <- p_breaks[p_breaks < max_val]
  
  # 4. Check proximity: If the last kept break is too close to Max (< 25% of step), remove it
  if (length(kept_breaks) > 0) {
    last_b <- tail(kept_breaks, 1)
    if ((max_val - last_b) < (0.25 * step_size)) {
      kept_breaks <- head(kept_breaks, -1) # Drop the cluttering label
    }
  }
  
  # 5. Assemble final breaks
  breaks_linear <- c(kept_breaks, max_val)
  
  # --- COLOR LOGIC (Sqrt distributed) ---
  n_cols <- 256
  colors_pal <- scales::viridis_pal(option = "magma", direction = 1)(n_cols)
  values_dist <- seq(0, 1, length.out = n_cols)^2 
  
  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x=x, y=y, fill=val)) +
    scale_fill_gradientn(
      colors = colors_pal, 
      values = values_dist, 
      breaks = breaks_linear,
      name = legend_title, 
      limits = c(0, max_val), 
      labels = round(breaks_linear, 0)
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

plot_proportion_map <- function(r, legend_title = "Threatened by\nwood production", title_txt = "") {
  r_mask <- terra::mask(r, vect(world_map))
  df <- as.data.frame(r_mask, xy = TRUE) %>%
    setNames(c("x", "y", "val")) %>%
    filter(!is.na(val))
  if (nrow(df) == 0) return(NULL)

  max_val <- max(df$val, na.rm = TRUE)
  if (max_val <= 0) return(NULL)

  breaks_linear <- pretty(c(0, max_val), n = 5)
  step_size <- if (length(breaks_linear) > 1) {
    breaks_linear[2] - breaks_linear[1]
  } else {
    max_val
  }
  kept_breaks <- breaks_linear[breaks_linear >= 0 & breaks_linear < max_val]
  if (length(kept_breaks) > 0 && (max_val - tail(kept_breaks, 1)) < (0.25 * step_size)) {
    kept_breaks <- head(kept_breaks, -1)
  }
  breaks_linear <- c(kept_breaks, max_val)
  
  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = val)) +
    scale_fill_viridis_c(
      option = "magma",
      direction = 1,
      name = legend_title,
      limits = c(0, max_val),
      breaks = breaks_linear,
      labels = scales::percent_format(accuracy = 1),
      oob = scales::squish
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.25, "cm"),
      legend.title = element_text(size = map_legend_title_size),
      legend.text = element_text(size = map_legend_text_size)
    )
}

# Map: Discrete Hotspots (Top 10, 5, 1%)
plot_hotspot_discrete <- function(r, title_suffix = "") {
  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE); vals <- vals[vals > 0 & !is.na(vals)]
  if (length(vals) == 0) return(NULL)
  
  q <- quantile(vals, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)
  df <- as.data.frame(r_mask, xy = TRUE) %>% setNames(c("x", "y", "val")) %>% filter(val > 0) %>%
    mutate(cat = case_when(val >= q[3] ~ "Top 1%", val >= q[2] ~ "Top 5%", val >= q[1] ~ "Top 10%", TRUE ~ "Present")) %>%
    mutate(cat = factor(cat, levels = c("Present", "Top 10%", "Top 5%", "Top 1%")))
  
  cols <- c("Present" = "grey80", "Top 10%" = "#FFEDA0", "Top 5%" = "#FEB24C", "Top 1%" = "#F03B20")
  
  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x=x, y=y, fill=cat)) +
    scale_fill_manual(values = cols, name = "Hotspot Level") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(0.9, "cm"),
      legend.title = element_text(size = map_legend_title_size),
      legend.text = element_text(size = map_legend_text_size)
    )
}

plot_proportion_hotspot_discrete <- function(r, title_txt = "") {
  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE)
  vals <- vals[vals > 0 & !is.na(vals)]
  if (length(vals) == 0) return(NULL)
  
  q <- quantile(vals, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)
  df <- as.data.frame(r_mask, xy = TRUE) %>%
    setNames(c("x", "y", "val")) %>%
    filter(val > 0 & !is.na(val)) %>%
    mutate(cat = case_when(
      val >= q[3] ~ "Top 1%",
      val >= q[2] ~ "Top 5%",
      val >= q[1] ~ "Top 10%",
      TRUE ~ "Present"
    )) %>%
    mutate(cat = factor(cat, levels = c("Present", "Top 10%", "Top 5%", "Top 1%")))
  
  cols <- c("Present" = "grey80", "Top 10%" = "#FFEDA0", "Top 5%" = "#FEB24C", "Top 1%" = "#F03B20")
  
  ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = cat)) +
    scale_fill_manual(values = cols, name = "Proportion\nHotspot Level") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(0.9, "cm"),
      legend.title = element_text(size = map_legend_title_size),
      legend.text = element_text(size = map_legend_text_size)
    )
}

# Scatter: Total vs Threatened (NO 1:1 LINE, BLACK LOESS)
create_richness_scatter <- function(r_total, r_threatened, title_txt = "") {
  rx_masked <- terra::mask(r_total, vect(world_map)); ry_masked <- terra::mask(r_threatened, vect(world_map))
  
  realms_proj <- st_transform(realms_sf, crs(rx_masked))
  realms_proj$rid <- as.numeric(as.factor(realms_proj[[realm_col_name]]))
  realm_r <- terra::rasterize(vect(realms_proj), rx_masked, field = "rid", touches = TRUE)
  
  target_mask <- rx_masked; target_mask[!is.na(target_mask)] <- 1 
  for(i in 1:30) {
    missing_mask <- !is.na(target_mask) & is.na(realm_r)
    if (global(missing_mask, "sum", na.rm=TRUE)[1,1] == 0) break
    win_size <- if(i <= 5) 3 else if(i <= 15) 5 else 7
    filled <- terra::focal(realm_r, w=win_size, fun="modal", na.policy="only", na.rm=TRUE)
    realm_r <- terra::cover(realm_r, filled)
  }
  
  stack_data <- c(rx_masked, ry_masked, realm_r)
  names(stack_data) <- c("Total", "Threatened", "realm_id")
  df <- as.data.frame(stack_data, xy = FALSE) %>% filter(!is.na(Total) & !is.na(Threatened) & Total > 0)
  
  lookup <- unique(st_drop_geometry(realms_proj)[, c("rid", realm_col_name)])
  colnames(lookup) <- c("realm_id", "Realm_Name")
  df <- df %>% left_join(lookup, by = "realm_id") %>% filter(!is.na(Realm_Name))
  
  ggplot(df, aes(x = Total, y = Threatened)) +
    geom_point(aes(color = Realm_Name), alpha = 1, size = 0.75) +
    geom_smooth(method = "loess", linewidth = 1, colour = "black", se = FALSE, span = 0.1) + 
    scale_color_manual(values = my_colors, name = "Zoogeographic Realm") +
    labs(x = "Total Species Richness", y = "Threatened Species Richness", title = title_txt) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = scatter_axis_title_size),
      axis.text = element_text(size = scatter_axis_text_size)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

# Scatter: Logging vs Plantation (1:1 LINE, ASPECT RATIO 1, NO GRIDS)
create_source_scatter <- function(r_x, r_y, title_txt = "") {
  rx_masked <- terra::mask(r_x, vect(world_map)); ry_masked <- terra::mask(r_y, vect(world_map))
  realms_proj <- st_transform(realms_sf, crs(rx_masked))
  realms_proj$rid <- as.numeric(as.factor(realms_proj[[realm_col_name]]))
  realm_r <- terra::rasterize(vect(realms_proj), rx_masked, field = "rid", touches = TRUE)
  target_mask <- rx_masked; target_mask[!is.na(target_mask)] <- 1 
  for(i in 1:30) {
    missing_mask <- !is.na(target_mask) & is.na(realm_r)
    if (global(missing_mask, "sum", na.rm=TRUE)[1,1] == 0) break
    win_size <- if(i <= 5) 3 else if(i <= 15) 5 else 7
    filled <- terra::focal(realm_r, w=win_size, fun="modal", na.policy="only", na.rm=TRUE)
    realm_r <- terra::cover(realm_r, filled)
  }
  stack_data <- c(rx_masked, ry_masked, realm_r)
  names(stack_data) <- c("X", "Y", "realm_id")
  df <- as.data.frame(stack_data, xy = FALSE) %>% filter(!is.na(X) & !is.na(Y) & (X > 0 | Y > 0))
  lookup <- unique(st_drop_geometry(realms_proj)[, c("rid", realm_col_name)])
  colnames(lookup) <- c("realm_id", "Realm_Name")
  df <- df %>% left_join(lookup, by = "realm_id") %>% filter(!is.na(Realm_Name))
  
  max_val <- max(c(df$X, df$Y), na.rm=TRUE)
  
  ggplot(df, aes(x = X, y = Y)) +
    geom_point(aes(color = Realm_Name), alpha = 1, size = 0.75) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
    scale_x_continuous(limits = c(0, max_val)) +
    scale_y_continuous(limits = c(0, max_val)) +
    scale_color_manual(values = my_colors, name = "Zoogeographic Realm") +
    labs(x = "Logging Threat Richness", y = "Plantation Threat Richness", title = title_txt) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = scatter_axis_title_size),
      axis.text = element_text(size = scatter_axis_text_size)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}

# Map: Bivariate (Legend Bottom Left)
plot_bivariate_map <- function(r_log, r_plant, title_txt = "") {
  s <- c(r_log, r_plant); names(s) <- c("Logging", "Plantation")
  s <- terra::mask(s, vect(world_map))
  df <- as.data.frame(s, xy = TRUE) %>% filter(Logging > 0 | Plantation > 0)
  if(nrow(df) == 0) return(NULL)
  
  df$Logging <- df$Logging + runif(nrow(df), 0, 0.000001)
  df$Plantation <- df$Plantation + runif(nrow(df), 0, 0.000001)
  
  get_cuts <- function(x) { 
    nonzero <- x[x > 0.001]
    if(length(nonzero)<10) return(NULL)
    quantile(nonzero, probs = c(0.90, 0.95, 0.99), na.rm = TRUE) 
  }
  x_cuts <- get_cuts(df$Logging); y_cuts <- get_cuts(df$Plantation)
  if(is.null(x_cuts) || is.null(y_cuts)) return(NULL)
  
  df <- df %>% mutate(
    x_cat = case_when(Logging < x_cuts[1] ~ 1, Logging < x_cuts[2] ~ 2, Logging < x_cuts[3] ~ 3, TRUE ~ 4),
    y_cat = case_when(Plantation < y_cuts[1] ~ 1, Plantation < y_cuts[2] ~ 2, Plantation < y_cuts[3] ~ 3, TRUE ~ 4),
    bi_class = paste0(x_cat, "-", y_cat)
  )
  
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = bi_class)) +
    bi_scale_fill(pal = "BlueYl", dim = 4) +
    theme_void() + 
    theme(legend.position = "none")
  
  leg <- bi_legend(pal = "BlueYl", dim = 4, size = 8, xlab = "Natural Forestry", ylab = "Plantation")
  
  final_plot <- p_map +
    patchwork::inset_element(
      leg,
      left = 0.02,
      bottom = 0.02,
      right = 0.26,
      top = 0.26,
      align_to = "panel"
    )
  
  return(final_plot)
}

plot_proportion_richness_bivariate_hotspots <- function(r_prop, r_richness, title_txt = "") {
  s <- c(r_prop, r_richness)
  names(s) <- c("Proportion", "Richness")
  s <- terra::mask(s, vect(world_map))
  df <- as.data.frame(s, xy = TRUE) %>%
    filter(!is.na(Proportion), !is.na(Richness), Proportion > 0, Richness > 0)
  if (nrow(df) == 0) return(NULL)
  
  prop_q <- quantile(df$Proportion, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)
  rich_q <- quantile(df$Richness, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)
  
  df <- df %>%
    mutate(
      prop_class = case_when(
        Proportion >= prop_q[3] ~ 4L,
        Proportion >= prop_q[2] ~ 3L,
        Proportion >= prop_q[1] ~ 2L,
        TRUE ~ 1L
      ),
      rich_class = case_when(
        Richness >= rich_q[3] ~ 4L,
        Richness >= rich_q[2] ~ 3L,
        Richness >= rich_q[1] ~ 2L,
        TRUE ~ 1L
      ),
      bi_class = paste0(prop_class, "-", rich_class)
    )
  
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = bi_class)) +
    bi_scale_fill(pal = "BlueYl", dim = 4) +
    theme_void() +
    theme(legend.position = "none")
  
  leg <- bi_legend(
    pal = "BlueYl",
    dim = 4,
    size = 8,
    xlab = "Proportion",
    ylab = "Richness"
  )
  
  p_map +
    patchwork::inset_element(
      leg,
      left = 0.02,
      bottom = 0.02,
      right = 0.28,
      top = 0.28,
      align_to = "panel"
    )
}

# --- 5. EXECUTION ---

all_pot   <- rast(grid_template, vals=0)
all_log   <- rast(grid_template, vals=0)
all_plant <- rast(grid_template, vals=0)
all_comb  <- rast(grid_template, vals=0)

taxa_rasters <- list() # For overlap analysis

taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")
inputs <- list(
  bird      = list(range = "data/processed/ranges/all_bird_ranges_dissolved.rds"),
  mammal    = list(range = "data/processed/ranges/all_mammal_ranges_dissolved.rds"),
  amphibian = list(range = "data/processed/ranges/all_amphibian_ranges_dissolved.rds"),
  reptile   = list(range = "data/processed/ranges/all_reptile_ranges_dissolved.rds"),
  fish      = list(range = "data/processed/ranges/all_fish_ranges_dissolved.rds")
)

# PART A: INDIVIDUAL TAXA LOOP
for (taxon in taxa_list) {
  cat("\nProcessing:", toupper(taxon), "\n")
  
  # 1. Load/Calc Rasters
  path_pot   <- raster_path(paste0(taxon, "_potential_richness_raster.tif"))
  path_log   <- raster_path(paste0(taxon, "_logging_richness_raster.tif"))
  path_plant <- raster_path(paste0(taxon, "_plantation_richness_raster.tif"))
  path_comb  <- raster_path(paste0(taxon, "_combined_richness_raster.tif"))
  path_prop  <- raster_path(paste0(taxon, "_wood_production_threatened_proportion_raster.tif"))
  
  if (file.exists(path_pot) & file.exists(path_log) & file.exists(path_plant) & file.exists(path_comb)) {
    r_pot <- rast(path_pot); r_log <- rast(path_log); r_plant <- rast(path_plant); r_comb <- rast(path_comb)
  } else {
    assessment_df <- load_assessment_species_summary(taxon) %>%
      filter(!category %in% c("EX", "EW"))
    ranges <- readRDS(here::here(inputs[[taxon]]$range))
    r_pot   <- get_richness_raster(unique(assessment_df$scientific_name), ranges, path_pot)
    r_log   <- get_richness_raster(unique(assessment_df %>% filter(has_logging_ongoing) %>% pull(scientific_name)), ranges, path_log)
    r_plant <- get_richness_raster(unique(assessment_df %>% filter(has_plantation_ongoing) %>% pull(scientific_name)), ranges, path_plant)
    r_comb  <- get_richness_raster(unique(assessment_df %>% filter(has_wood_prod_ongoing) %>% pull(scientific_name)), ranges, path_comb)
  }
  r_prop <- get_threatened_proportion_raster(r_comb, r_pot, path_prop)
  
  taxa_rasters[[taxon]] <- r_comb # Store for overlap
  
  # 2. Accumulate
  all_pot   <- all_pot + r_pot
  all_log   <- all_log + r_log
  all_plant <- all_plant + r_plant
  all_comb  <- all_comb + r_comb
  
  # 3. Generate & Save Individual Plots
  
  # Gradient: Total Richness
  p_grad1 <- plot_gradient_sqrt(r_pot, legend_title = "Species Richness", title_txt = paste(tools::toTitleCase(taxon), "- Total Richness"))
  if(!is.null(p_grad1)) ggsave(file.path(out_dir, paste0(taxon, "_gradient_total.png")), p_grad1, width = 10, height = 6, bg = "white")
  
  # Gradient: Threatened Richness
  p_grad2 <- plot_gradient_sqrt(r_comb, legend_title = "Threatened\nRichness", title_txt = paste(tools::toTitleCase(taxon), "- Threatened Richness"))
  if(!is.null(p_grad2)) ggsave(file.path(out_dir, paste0(taxon, "_gradient_threatened.png")), p_grad2, width = 10, height = 6, bg = "white")
  
  # Gradient: Proportion Threatened by Wood Production
  p_prop <- plot_proportion_map(r_prop, title_txt = paste(tools::toTitleCase(taxon), "- Proportion Threatened by Wood Production"))
  if(!is.null(p_prop)) ggsave(file.path(out_dir, paste0(taxon, "_proportion_wood_production_threatened.png")), p_prop, width = 10, height = 6, bg = "white")

  # Hotspots: Proportion Threatened by Wood Production
  p_prop_hot <- plot_proportion_hotspot_discrete(r_prop, title_txt = paste(tools::toTitleCase(taxon), "- Proportion Threatened Hotspots"))
  if(!is.null(p_prop_hot)) ggsave(file.path(out_dir, paste0(taxon, "_proportion_wood_production_threatened_hotspots.png")), p_prop_hot, width = 10, height = 6, bg = "white")

  # Bivariate hotspots: Proportion vs Threatened Richness
  p_prop_rich_bi <- plot_proportion_richness_bivariate_hotspots(
    r_prop,
    r_comb,
    title_txt = paste(tools::toTitleCase(taxon), "- Proportion vs Richness Hotspots")
  )
  if(!is.null(p_prop_rich_bi)) {
    ggsave(
      file.path(out_dir, paste0(taxon, "_proportion_vs_richness_bivariate_hotspots.png")),
      p_prop_rich_bi,
      width = 10,
      height = 7,
      bg = "white"
    )
  }
  
  # Bivariate
  p_bi <- plot_bivariate_map(r_log, r_plant, title_txt = paste(tools::toTitleCase(taxon), "- Logging vs Plantation"))
  if(!is.null(p_bi)) ggsave(file.path(out_dir, paste0(taxon, "_bivariate.png")), p_bi, width = 10, height = 7, bg = "white")
  
  # Hotspots
  p_hot <- plot_hotspot_discrete(r_comb, title_suffix = paste("(", taxon, ")"))
  if(!is.null(p_hot)) ggsave(file.path(out_dir, paste0(taxon, "_hotspots.png")), p_hot, width = 10, height = 6, bg = "white")
  
  # Scatter: Total vs Threatened
  p_scat1 <- create_richness_scatter(r_pot, r_comb, title_txt = paste(tools::toTitleCase(taxon), "- Total vs Threatened"))
  if(!is.null(p_scat1)) ggsave(file.path(out_dir, paste0(taxon, "_scatter_richness.png")), p_scat1, width = 8, height = 8, bg = "white")
  
  # Scatter: Sources
  p_scat2 <- create_source_scatter(r_log, r_plant, title_txt = paste(tools::toTitleCase(taxon), "- Source Comparison"))
  if(!is.null(p_scat2)) ggsave(file.path(out_dir, paste0(taxon, "_scatter_sources.png")), p_scat2, width = 8, height = 8, bg = "white")
}

# PART B: ALL VERTEBRATES (Combined)
cat("\nProcessing: ALL VERTEBRATES\n")

path_all_prop <- raster_path("all_vertebrates_wood_production_threatened_proportion_raster.tif")
all_prop <- get_threatened_proportion_raster(all_comb, all_pot, path_all_prop)

# Gradient: Total
p_all_grad1 <- plot_gradient_sqrt(all_pot, legend_title = "Total Richness", title_txt = "All Vertebrates - Total Richness")
if(!is.null(p_all_grad1)) ggsave(file.path(out_dir, "all_vertebrates_gradient_total.png"), p_all_grad1, width = 10, height = 6, bg = "white")

# Gradient: Threatened
p_all_grad2 <- plot_gradient_sqrt(all_comb, legend_title = "Threatened\nRichness", title_txt = "All Vertebrates - Threatened Richness")
if(!is.null(p_all_grad2)) ggsave(file.path(out_dir, "all_vertebrates_gradient_threatened.png"), p_all_grad2, width = 10, height = 6, bg = "white")

# Gradient: Proportion Threatened by Wood Production
p_all_prop <- plot_proportion_map(all_prop, title_txt = "All Vertebrates - Proportion Threatened by Wood Production")
if(!is.null(p_all_prop)) ggsave(file.path(out_dir, "all_vertebrates_proportion_wood_production_threatened.png"), p_all_prop, width = 10, height = 6, bg = "white")

# Hotspots: Proportion Threatened by Wood Production
p_all_prop_hot <- plot_proportion_hotspot_discrete(all_prop, title_txt = "All Vertebrates - Proportion Threatened Hotspots")
if(!is.null(p_all_prop_hot)) ggsave(file.path(out_dir, "all_vertebrates_proportion_wood_production_threatened_hotspots.png"), p_all_prop_hot, width = 10, height = 6, bg = "white")

# Bivariate hotspots: Proportion vs Threatened Richness
p_all_prop_rich_bi <- plot_proportion_richness_bivariate_hotspots(
  all_prop,
  all_comb,
  title_txt = "All Vertebrates - Proportion vs Richness Hotspots"
)
if(!is.null(p_all_prop_rich_bi)) {
  ggsave(
    file.path(out_dir, "all_vertebrates_proportion_vs_richness_bivariate_hotspots.png"),
    p_all_prop_rich_bi,
    width = 10,
    height = 7,
    bg = "white"
  )
}

# Bivariate
p_all_bi <- plot_bivariate_map(all_log, all_plant, title_txt = "All Vertebrates - Logging vs Plantation")
if(!is.null(p_all_bi)) ggsave(file.path(out_dir, "all_vertebrates_bivariate.png"), p_all_bi, width = 10, height = 7, bg = "white")

# Hotspots
p_all_hot <- plot_hotspot_discrete(all_comb, title_suffix = "(All Vertebrates)")
if(!is.null(p_all_hot)) ggsave(file.path(out_dir, "all_vertebrates_hotspots.png"), p_all_hot, width = 10, height = 6, bg = "white")

# Scatter: Total vs Threatened
p_all_scat1 <- create_richness_scatter(all_pot, all_comb, title_txt = "All Vertebrates - Total vs Threatened")
if(!is.null(p_all_scat1)) ggsave(file.path(out_dir, "all_vertebrates_scatter_richness.png"), p_all_scat1, width = 8, height = 8, bg = "white")

# Scatter: Sources
p_all_scat2 <- create_source_scatter(all_log, all_plant, title_txt = "All Vertebrates - Source Comparison")
if(!is.null(p_all_scat2)) ggsave(file.path(out_dir, "all_vertebrates_scatter_sources.png"), p_all_scat2, width = 8, height = 8, bg = "white")

# PART C: HOTSPOT OVERLAP FREQUENCY
cat("\nProcessing: HOTSPOT OVERLAP FREQUENCY\n")

hotspot_freq <- rast(grid_template, vals = 0)
for (taxon in taxa_list) {
  r <- taxa_rasters[[taxon]]
  r_mask <- terra::mask(r, vect(world_map))
  vals <- values(r_mask, mat = FALSE); vals <- vals[vals > 0 & !is.na(vals)]
  if (length(vals) > 0) {
    threshold_10 <- quantile(vals, probs = 0.90, na.rm = TRUE)
    r_binary <- (r_mask >= threshold_10)
    r_binary[is.na(r_binary)] <- 0
    hotspot_freq <- hotspot_freq + r_binary
  }
}

freq_df <- as.data.frame(terra::mask(hotspot_freq, vect(world_map)), xy = TRUE) %>% 
  setNames(c("x", "y", "count")) %>% filter(count > 0) %>%
  mutate(count = factor(count, levels = 1:5))

freq_colors <- c("1" = "#FED976", "2" = "#FEB24C", "3" = "#FD8D3C", "4" = "#F03B20", "5" = "#BD0026")

p_overlap <- ggplot() +
  geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
  geom_tile(data = freq_df, aes(x = x, y = y, fill = count)) +
  scale_fill_manual(values = freq_colors, name = "Taxa Count") +
  theme_void() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "all_vertebrates_hotspot_overlap_frequency.png"), p_overlap, width = 10, height = 6, bg = "white")

cat("\nAll Individual, Combined, and Overlap Maps Saved.\n")
