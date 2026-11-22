# =============================================================================
#
# SCRIPT: GENERATE GLOBAL THREAT & SENSITIVITY MAPS FOR VERTEBRATE BIODIVERSITY
#
# Date: 28 November 2025
#
# Description:
# This script produces global hotspot maps for different vertebrate
# classes based on the methodology of Lamb et al. (2024). It is designed to be
# run for any taxonomic group for which the required data is available.
#
# For each specified taxon, it calculates and maps:
#
#   1. "Species Threat Value" (ST) - Linear Inverse Range Weighting
#   2. "Community Sensitivity" (CS) - ST / Potential Species Threat
#   2. "Transformed Threat Value" (TT) - Cube root inverse range weighting
#   3. "Transformed CS" - TT / Potential Transformed Threat
#   4. "Species Richness" (SR) - Simple count of threatened species in each cell
#   6. "Richness CS" - SR / Total Cell Richness
#
# It does this for SEVEN threat categories:
#   1. Logging & Wood Harvesting (IUCN Code 5.3)
#   2. Wood & Pulp Plantations (IUCN Code 2.2)
#   3. A combined total timber production threat (5.3 & 2.2)
#   4. Fuelwood/Charcoal collection (derived from Red List Documnentation text)
#   5. Industrial Roundwood (derived from Red List Documentation text)
#   6. Small-scale wood production (2.2.1 & 5.3.3)
#   7. Large-scale wood production (2.2.2 & 5.3.4)
#
# =============================================================================


# --- 1. SETUP: Load necessary libraries ---
# -----------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,          # For data manipulation
  tidyr,          # For data tidying
  ggplot2,        # For creating plots
  here,           # For file path management
  sf,             # For working with spatial vector data
  terra,          # For working with spatial raster data
  rnaturalearth,  # For getting a world map
  stringr,        # For string manipulation
  rredlist,       # For rl_assessment_extract function
  cowplot,        # For arranging multiple ggplot objects
  biscale         # For creating bivariate color palettes
)


# --- 2. GLOBAL SETUP & HELPER FUNCTIONS ---
# -----------------------------------------------------------------------------
# These settings and functions are defined globally as they are common to all
# taxonomic analyses.

# --- Global Spatial Grid ---
# Using a Molleweide equal area projection
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent_moll <- ext(-18040000, 18040000, -9020000, 9020000)
grid <- rast(resolution = 110000, crs = moll_crs, extent = global_extent_moll, vals = 0)
cell_area <- res(grid)[1] * res(grid)[2]
world_map_moll_data <- ne_countries(scale = "medium", returnclass = "sf")
world_map_moll <- st_transform(
  filter(world_map_moll_data, name != "Antarctica"),
  moll_crs
)

# --- Global Plotting Palettes ---
hotspot_colors <- c("Threat Present" = "darkgrey", 
                    "Top 10% Hotspot" = "yellow",
                    "Top 5% Hotspot" = "orange", 
                    "Top 1% Hotspot" = "red")

sensitivity_colors <- c("Sensitivity Present" = "darkgrey", 
                        "Top 10% Hotspot" = "yellow", 
                        "Top 5% Hotspot" = "orange", 
                        "Top 1% Hotspot" = "red")

# Canonical raster library location
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

# Canonical figure output locations
maps_dir <- here::here("outputs", "hotspots", "figures", "maps")
realm_dir <- here::here("outputs", "hotspots", "figures", "realm")
summary_dir <- here::here("outputs", "hotspots", "figures", "summary")
dir.create(maps_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(realm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

# --- Function 1: Calculate or Load Threat Raster (Linear Range-Weighted) ---
calculate_or_load_raster <- function(species_names, all_ranges_sf, grid_template, file_path) {
  if (file.exists(file_path)) {
    cat(paste("\nLoading cached raster from:", basename(file_path), "\n"))
    return(rast(file_path))
  }
  cat(paste("\n--- Starting calculation for:", basename(file_path), "---\n"))
  
  species_sf <- all_ranges_sf %>% filter(sci_name %in% species_names)
  species_sf <- species_sf %>%
    mutate(total_area_m2 = st_area(st_transform(geom, crs = crs(grid_template))))
  
  final_raster <- rast(grid_template, vals = 0)
  total_species <- nrow(species_sf)
  start_time <- Sys.time()
  
  if (total_species == 0) {
    cat("No species in the provided list. Returning an empty raster.\n")
    return(final_raster)
  }
  
  for (i in 1:total_species) {
    current_species <- species_sf[i, ]
    cat("   -> Processing", i, "of", total_species, ":", current_species$sci_name, "\n")
    total_area <- current_species$total_area_m2
    if (as.numeric(total_area) == 0) next
    
    coverage_raster <- terra::rasterize(vect(st_transform(current_species, crs = crs(grid_template))), grid_template, cover = TRUE)
    area_in_cell <- coverage_raster * cell_area
    threat_certainty <- area_in_cell / as.numeric(total_area)
    threat_certainty[threat_certainty > 1] <- 1
    threat_certainty[is.na(threat_certainty)] <- 0
    final_raster <- final_raster + threat_certainty
    
    elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
    avg_time_per_species <- as.numeric(elapsed_time) / i
    estimated_time_remaining <- round(avg_time_per_species * (total_species - i), 1)
    cat("       Time elapsed:", round(as.numeric(elapsed_time), 1), "mins. | Est. time remaining:", estimated_time_remaining, "mins.\n")
  }
  
  cat(paste("Calculation complete. Saving raster to:", basename(file_path), "for next time...\n"))
  writeRaster(final_raster, file_path, overwrite = TRUE)
  return(final_raster)
}

# --- Function 2a: Create Hotspot Map ---
create_hotspot_map <- function(input_raster, map_title, legend_title, color_palette, 
                               add_title = TRUE, remove_axis_elements = TRUE, show_legend = TRUE) {
  cat("Generating map for:", map_title, "\n")
  raster_land <- terra::mask(input_raster, vect(world_map_moll))
  raster_values <- values(raster_land)
  raster_values <- raster_values[raster_values > 0 & !is.na(raster_values)]
  if (length(raster_values) == 0) {
    cat(" -> No data found for this category after masking. Skipping map.\n")
    return(ggplot() + theme_void() + geom_text(aes(x=0.5, y=0.5, label="No Data")) + labs(title=map_title))
  }
  
  top_10_threshold <- quantile(raster_values, 0.90, na.rm = TRUE)
  top_5_threshold <- quantile(raster_values, 0.95, na.rm = TRUE)
  top_1_threshold <- quantile(raster_values, 0.99, na.rm = TRUE)
  
  plot_df <- as.data.frame(raster_land, xy = TRUE) %>%
    setNames(c("x", "y", "score")) %>%
    filter(score > 0 & !is.na(score)) %>%
    mutate(category = case_when(
      score >= top_1_threshold ~ "Top 1% Hotspot",
      score >= top_5_threshold ~ "Top 5% Hotspot",
      score >= top_10_threshold ~ "Top 10% Hotspot",
      TRUE ~ names(color_palette)[1]
    )) %>%
    mutate(category = factor(category, levels = names(color_palette)))
  
  map_plot <- ggplot() +
    # Base land is a very light grey
    geom_sf(data = world_map_moll, fill = "#f0f0f0", color = NA) + 
    geom_tile(data = plot_df, aes(x = x, y = y, fill = category)) +
    scale_fill_manual(values = color_palette, name = legend_title) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_minimal()
  
  if (add_title) {
    map_plot <- map_plot + labs(title = map_title)
  } else {
    map_plot <- map_plot + labs(title = NULL)
  }
  
  theme_settings <- list()
  
  if (remove_axis_elements) {
    theme_settings$axis.title <- element_blank()
    theme_settings$axis.text <- element_blank()
    theme_settings$axis.ticks <- element_blank()
    theme_settings$panel.grid <- element_blank()
    
    theme_settings$plot.margin <- unit(c(0, 0, 0, 5), "mm") 
  }
  
  if (!show_legend) {
    theme_settings$legend.position <- "none"
  } else {
    theme_settings$legend.position <- "right"
  }
  
  map_plot <- map_plot + theme(!!!theme_settings) 
  
  return(map_plot)
}

create_continuous_hotspot_map <- function(input_raster, map_title, legend_title, 
                                          add_title = TRUE, remove_axis_elements = TRUE, 
                                          show_legend = TRUE, 
                                          transform_scale = "identity") { # default is linear
  
  cat("Generating map for:", map_title, "\n")
  raster_land <- terra::mask(input_raster, vect(world_map_moll))
  raster_values <- values(raster_land)
  raster_values <- raster_values[raster_values > 0 & !is.na(raster_values)]
  
  if (length(raster_values) == 0) {
    cat(" -> No data found. Skipping map.\n")
    return(ggplot() + theme_void() + geom_text(aes(x=0.5, y=0.5, label="No Data")) + labs(title=map_title))
  }
  
  plot_df <- as.data.frame(raster_land, xy = TRUE) %>%
    setNames(c("x", "y", "score")) %>%
    filter(score > 0 & !is.na(score)) 
  
  map_plot <- ggplot() +
    geom_sf(data = world_map_moll, fill = "lightgrey", color = NA) + 
    geom_tile(data = plot_df, aes(x = x, y = y, fill = score)) +
    
    # Use the 'trans' argument here
    scale_fill_viridis_c(option = "magma", name = legend_title, direction = 1, 
                         trans = transform_scale) + # <--- Applied here
    
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_minimal()
  
  if (add_title) { map_plot <- map_plot + labs(title = map_title) } 
  else { map_plot <- map_plot + labs(title = NULL) }
  
  theme_settings <- list()
  if (remove_axis_elements) {
    theme_settings$axis.title <- element_blank()
    theme_settings$axis.text <- element_blank()
    theme_settings$axis.ticks <- element_blank()
    theme_settings$panel.grid <- element_blank()
    theme_settings$plot.margin <- unit(c(0, 0, 0, 5), "mm") 
  }
  if (!show_legend) { theme_settings$legend.position <- "none" } 
  else { 
    theme_settings$legend.position <- "bottom"
  }
  
  map_plot <- map_plot + theme(!!!theme_settings) 
  
  return(map_plot)
}

# --- Function 2b: Continuous map with dark grey ----
# -----------------------------------------------------------------------------
create_continuous_hotspot_map2 <- function(input_raster, map_title, legend_title, 
                                          add_title = TRUE, remove_axis_elements = TRUE, 
                                          show_legend = TRUE, 
                                          transform_scale = "identity") { # default is linear
  
  cat("Generating map for:", map_title, "\n")
  
  # Masking and basic prep
  raster_land <- terra::mask(input_raster, terra::vect(world_map_moll))
  
  # Create the main dataframe
  plot_df <- as.data.frame(raster_land, xy = TRUE) %>%
    setNames(c("x", "y", "score")) %>%
    filter(score > 0 & !is.na(score)) 
  
  # Check for empty data
  if (nrow(plot_df) == 0) {
    cat(" -> No data found. Skipping map.\n")
    return(ggplot() + theme_void() + geom_text(aes(x=0.5, y=0.5, label="No Data")) + labs(title=map_title))
  }
  
  # --- NEW LOGIC START ---
  # 1. Calculate the median (50th percentile) of the scores
  cutoff_value <- median(plot_df$score, na.rm = TRUE)
  cat(" -> Bottom half cutoff value:", cutoff_value, "\n")
  
  # 2. Split data into 'bottom half' (low) and 'top half' (high)
  df_low  <- plot_df %>% filter(score <= cutoff_value)
  df_high <- plot_df %>% filter(score > cutoff_value)
  # --- NEW LOGIC END ---
  
  map_plot <- ggplot() +
    geom_sf(data = world_map_moll, fill = "lightgrey", color = NA) + 
    
    # Layer 1: The bottom half (Static Dark Grey)
    # Note: fill is outside aes() to make it a static color, not a mapped scale
    geom_tile(data = df_low, aes(x = x, y = y), fill = "darkgrey") +
    
    # Layer 2: The top half (Gradient)
    geom_tile(data = df_high, aes(x = x, y = y, fill = score)) +
    
    # The scale will only apply to Layer 2 because Layer 1 has no mapped fill
    scale_fill_viridis_c(option = "magma", name = legend_title, direction = 1, 
                         trans = transform_scale) + 
    
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_minimal()
  
  # Handle Titles
  if (add_title) { map_plot <- map_plot + labs(title = map_title) } 
  else { map_plot <- map_plot + labs(title = NULL) }
  
  # Handle Theme
  theme_settings <- list()
  if (remove_axis_elements) {
    theme_settings$axis.title <- element_blank()
    theme_settings$axis.text <- element_blank()
    theme_settings$axis.ticks <- element_blank()
    theme_settings$panel.grid <- element_blank()
    theme_settings$plot.margin <- unit(c(0, 0, 0, 5), "mm") 
  }
  if (!show_legend) { theme_settings$legend.position <- "none" } 
  else { 
    theme_settings$legend.position <- "right"
  }
  
  map_plot <- map_plot + theme(!!!theme_settings) 
  
  return(map_plot)
}

# --- Function 3: Calculate or Load Richness Raster ---
# -----------------------------------------------------------------------------
calculate_or_load_richness_raster <- function(species_names, all_ranges_sf, grid_template, file_path) {
  if (file.exists(file_path)) {
    cat(paste("\nLoading cached richness raster from:", basename(file_path), "\n"))
    return(rast(file_path))
  }
  cat(paste("\n--- Starting richness calculation for:", basename(file_path), "---\n"))
  
  species_sf <- all_ranges_sf %>% filter(sci_name %in% species_names)
  final_raster <- rast(grid_template, vals = 0)
  total_species <- nrow(species_sf)
  start_time <- Sys.time()
  
  if (total_species == 0) {
    cat("No species in the provided list. Returning an empty raster.\n")
    return(final_raster)
  }
  
  for (i in 1:total_species) {
    current_species <- species_sf[i, ]
    cat("   -> Processing Richness", i, "of", total_species, ":", current_species$sci_name, "\n")
    
    current_species_raster <- terra::rasterize(
      vect(st_transform(current_species, crs = crs(grid_template))), 
      grid_template, 
      field = 1,       # Burn a '1'
      background = NA, # Set background to NA
      touches = TRUE   # Include cells that are touched
    )
    
    # Convert NAs to 0s so they can be added
    current_species_raster[is.na(current_species_raster)] <- 0
    
    # Add this species' presence raster to the total sum
    final_raster <- final_raster + current_species_raster
    
    elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
    avg_time_per_species <- as.numeric(elapsed_time) / i
    estimated_time_remaining <- round(avg_time_per_species * (total_species - i), 1)
    cat("       Time elapsed:", round(as.numeric(elapsed_time), 1), "mins. | Est. time remaining:", estimated_time_remaining, "mins.\n")
  }
  
  cat(paste("Richness calculation complete. Saving raster to:", basename(file_path), "for next time...\n"))
  writeRaster(final_raster, file_path, overwrite = TRUE)
  return(final_raster)
}

# --- Function 4: Calculate or Load Transformed Threat Raster (Cube Root) ---
# NEW FUNCTION: Calculates threat certainty as (range in cell / range outside)^1/3
# -----------------------------------------------------------------------------
calculate_or_load_transformed_raster <- function(species_names, all_ranges_sf, grid_template, file_path) {
  if (file.exists(file_path)) {
    cat(paste("\nLoading cached transformed raster from:", basename(file_path), "\n"))
    return(rast(file_path))
  }
  cat(paste("\n--- Starting transformed calculation (Cube Root) for:", basename(file_path), "---\n"))
  
  species_sf <- all_ranges_sf %>% filter(sci_name %in% species_names)
  species_sf <- species_sf %>%
    mutate(total_area_m2 = st_area(st_transform(geom, crs = crs(grid_template))))
  
  final_raster <- rast(grid_template, vals = 0)
  total_species <- nrow(species_sf)
  start_time <- Sys.time()
  
  if (total_species == 0) {
    cat("No species in the provided list. Returning an empty raster.\n")
    return(final_raster)
  }
  
  for (i in 1:total_species) {
    current_species <- species_sf[i, ]
    cat("   -> Processing Transformed", i, "of", total_species, ":", current_species$sci_name, "\n")
    total_area <- current_species$total_area_m2
    if (as.numeric(total_area) == 0) next
    
    coverage_raster <- terra::rasterize(vect(st_transform(current_species, crs = crs(grid_template))), grid_template, cover = TRUE)
    area_in_cell <- coverage_raster * cell_area
    
    # Calculate base proportion (Range inside / Total Range)
    proportion_in_cell <- area_in_cell / as.numeric(total_area)
    proportion_in_cell[proportion_in_cell > 1] <- 1
    proportion_in_cell[is.na(proportion_in_cell)] <- 0
    
    # --- APPLY CUBE ROOT TRANSFORMATION ---
    threat_certainty_transformed <- proportion_in_cell^(1/3)
    
    final_raster <- final_raster + threat_certainty_transformed
    
    elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
    avg_time_per_species <- as.numeric(elapsed_time) / i
    estimated_time_remaining <- round(avg_time_per_species * (total_species - i), 1)
    cat("       Time elapsed:", round(as.numeric(elapsed_time), 1), "mins. | Est. time remaining:", estimated_time_remaining, "mins.\n")
  }
  
  cat(paste("Transformed calculation complete. Saving raster to:", basename(file_path), "for next time...\n"))
  writeRaster(final_raster, file_path, overwrite = TRUE)
  return(final_raster)
}

# --- Function 5: Save Map Helper ---
# -----------------------------------------------------------------------------
save_map_image <- function(plot_obj, taxon_name, threat_name, weighting_code, map_type = "Hotspot") {
  if (is.null(plot_obj)) return()
  
  # 1. Ensure directory exists
  out_dir <- maps_dir
  
  # 2. Construct Filename
  clean_taxon <- str_to_title(taxon_name)
  clean_taxon <- gsub(" ", "_", clean_taxon)
  filename <- paste0(clean_taxon, "_", threat_name, "_", weighting_code, "_", map_type, ".png")
  full_path <- file.path(out_dir, filename)
  
  cat(paste("   -> Saving:", filename, "\n"))
  
  # 3. Adjust Theme for Saving
  plot_to_save <- plot_obj +
    theme(
      plot.title = element_text(size = 12),       
      legend.title = element_text(size = 10),     
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # 4. Save
  ggsave(filename = full_path, 
         plot = plot_to_save, 
         width = 1500, 
         height = 750, 
         units = "px", 
         bg = "white", 
         dpi = 150)
}

# --- 3. MAIN PROCESSING FUNCTION ---
# -----------------------------------------------------------------------------
generate_threat_maps <- function(taxa_name, range_file) {
  
  taxa_name_capitalized <- str_to_title(taxa_name)
  cat(paste("\n\n--- STARTING ANALYSIS FOR:", toupper(taxa_name_capitalized), "---\n"))
  
  # --- 3a. Data Loading and Preparation ---
  cat("Loading pre-processed data...\n")
  assessment_df_consolidated <- load_assessment_species_summary(taxa_name)
  all_ranges_dissolved <- readRDS(range_file)
  
  # --- 3b. Define Species Lists ---
  extant_species_df <- assessment_df_consolidated %>% filter(!category %in% c("EX", "EW"))
  extant_species_names <- extant_species_df %>% pull(scientific_name)
  
  logging_species_names <- extant_species_df %>%
    filter(has_logging_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  plantation_species_names <- extant_species_df %>%
    filter(has_plantation_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  combined_species_names <- extant_species_df %>%
    filter(has_wood_prod_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  
  fuelwood_species_names <- extant_species_df %>%
    filter(has_fuelwood_charcoal_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  cat(paste("\nFound", length(fuelwood_species_names), "extant", taxa_name, "species threatened by fuelwood/charcoal production.\n"))
  
  industrial_roundwood_species_names <- extant_species_df %>%
    filter(has_industrial_roundwood_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  cat(paste("\nFound", length(industrial_roundwood_species_names), "extant", taxa_name, "species threatened by industrial roundwood production.\n"))
  
  small_scale_species_names <- extant_species_df %>%
    filter(has_small_scale_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  cat(paste("\nFound", length(small_scale_species_names), "extant", taxa_name, "species threatened by small-scale wood production.\n"))
  
  large_scale_species_names <- extant_species_df %>%
    filter(has_large_scale_ongoing) %>%
    pull(scientific_name) %>%
    unique()
  cat(paste("\nFound", length(large_scale_species_names), "extant", taxa_name, "species threatened by large-scale wood production.\n"))
  
  # --- 3c. Raster Generation ---
  potential_threat_raster <- calculate_or_load_raster(extant_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_potential_threat_raster.tif")))
  potential_trans_raster <- calculate_or_load_transformed_raster(extant_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_potential_transformed_threat_raster.tif")))
  sr_potential_raster <- calculate_or_load_richness_raster(extant_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_potential_richness_raster.tif")))
  
  # --- Logging ---
  logging_st_raster <- calculate_or_load_raster(logging_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_logging_threat_raster.tif")))
  logging_richness_raster <- calculate_or_load_richness_raster(logging_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_logging_richness_raster.tif")))
  logging_trans_raster <- calculate_or_load_transformed_raster(logging_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_logging_transformed_threat_raster.tif")))
  
  # --- Plantation ---
  plantation_st_raster <- calculate_or_load_raster(plantation_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_plantation_threat_raster.tif")))
  plantation_richness_raster <- calculate_or_load_richness_raster(plantation_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_plantation_richness_raster.tif")))
  plantation_trans_raster <- calculate_or_load_transformed_raster(plantation_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_plantation_transformed_threat_raster.tif")))
  
  # --- Combined ---
  combined_st_raster <- calculate_or_load_raster(combined_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_combined_threat_raster.tif")))
  combined_richness_raster <- calculate_or_load_richness_raster(combined_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_combined_richness_raster.tif")))
  combined_trans_raster <- calculate_or_load_transformed_raster(combined_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_combined_transformed_threat_raster.tif")))
  
  # --- Fuelwood ---
  fuelwood_st_raster <- calculate_or_load_raster(fuelwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_fuelwood_threat_raster.tif")))
  fuelwood_richness_raster <- calculate_or_load_richness_raster(fuelwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_fuelwood_richness_raster.tif")))
  fuelwood_trans_raster <- calculate_or_load_transformed_raster(fuelwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_fuelwood_transformed_threat_raster.tif")))
  
  # --- Industrial Roundwood ---
  industrial_roundwood_st_raster <- calculate_or_load_raster(industrial_roundwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_industrial_roundwood_threat_raster.tif")))
  industrial_roundwood_richness_raster <- calculate_or_load_richness_raster(industrial_roundwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_industrial_roundwood_richness_raster.tif")))
  industrial_roundwood_trans_raster <- calculate_or_load_transformed_raster(industrial_roundwood_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_industrial_roundwood_transformed_threat_raster.tif")))
  
  # --- Small-scale and Large-scale Rasters ---
  small_scale_st_raster <- calculate_or_load_raster(small_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_small_scale_threat_raster.tif")))
  small_scale_richness_raster <- calculate_or_load_richness_raster(small_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_small_scale_richness_raster.tif")))
  small_scale_trans_raster <- calculate_or_load_transformed_raster(small_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_small_scale_transformed_threat_raster.tif")))
  
  large_scale_st_raster <- calculate_or_load_raster(large_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_large_scale_threat_raster.tif")))
  large_scale_richness_raster <- calculate_or_load_richness_raster(large_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_large_scale_richness_raster.tif")))
  large_scale_trans_raster <- calculate_or_load_transformed_raster(large_scale_species_names, all_ranges_dissolved, grid, raster_path(paste0(taxa_name, "_large_scale_transformed_threat_raster.tif")))
  
  
  # Calculate Community Sensitivity Rasters
  cs_logging_raster <- ifel(potential_threat_raster > 0, logging_st_raster / potential_threat_raster, 0)
  cs_plantation_raster <- ifel(potential_threat_raster > 0, plantation_st_raster / potential_threat_raster, 0)
  cs_combined_raster <- ifel(potential_threat_raster > 0, combined_st_raster / potential_threat_raster, 0)
  cs_fuelwood_raster <- ifel(potential_threat_raster > 0, fuelwood_st_raster / potential_threat_raster, 0)
  cs_industrial_roundwood_raster <- ifel(potential_threat_raster > 0, industrial_roundwood_st_raster / potential_threat_raster, 0)
  cs_small_scale_raster <- ifel(potential_threat_raster > 0, small_scale_st_raster / potential_threat_raster, 0)
  cs_large_scale_raster <- ifel(potential_threat_raster > 0, large_scale_st_raster / potential_threat_raster, 0)
  
  # --- SAVE COMMUNITY SENSITIVITY RASTERS ---
  cat(paste("Saving Community Sensitivity rasters for", taxa_name, "...\n"))
  writeRaster(cs_logging_raster, raster_path(paste0(taxa_name, "_cs_logging_raster.tif")), overwrite = TRUE)
  writeRaster(cs_plantation_raster, raster_path(paste0(taxa_name, "_cs_plantation_raster.tif")), overwrite = TRUE)
  writeRaster(cs_combined_raster, raster_path(paste0(taxa_name, "_cs_combined_raster.tif")), overwrite = TRUE)
  writeRaster(cs_fuelwood_raster, raster_path(paste0(taxa_name, "_cs_fuelwood_raster.tif")), overwrite = TRUE)
  writeRaster(cs_industrial_roundwood_raster, raster_path(paste0(taxa_name, "_cs_industrial_roundwood_raster.tif")), overwrite = TRUE)
  
  # Calculate Transformed Community Sensitivity (CS Trans) for THIS Taxon ---
  cat(paste("Calculating Transformed Community Sensitivity (CS Trans) rasters for", taxa_name, "...\n"))
  
  cs_trans_logging_raster <- ifel(potential_trans_raster > 0, logging_trans_raster / potential_trans_raster, 0)
  cs_trans_plantation_raster <- ifel(potential_trans_raster > 0, plantation_trans_raster / potential_trans_raster, 0)
  cs_trans_combined_raster <- ifel(potential_trans_raster > 0, combined_trans_raster / potential_trans_raster, 0)
  cs_trans_fuelwood_raster <- ifel(potential_trans_raster > 0, fuelwood_trans_raster / potential_trans_raster, 0)
  cs_trans_industrial_roundwood_raster <- ifel(potential_trans_raster > 0, industrial_roundwood_trans_raster / potential_trans_raster, 0)
  cs_trans_small_scale_raster <- ifel(potential_trans_raster > 0, small_scale_trans_raster / potential_trans_raster, 0)
  cs_trans_large_scale_raster <- ifel(potential_trans_raster > 0, large_scale_trans_raster / potential_trans_raster, 0)
  
  # --- SAVE CS TRANS RASTERS ---
  cat(paste("Saving Transformed Community Sensitivity rasters for", taxa_name, "...\n"))
  writeRaster(cs_trans_logging_raster, raster_path(paste0(taxa_name, "_cs_trans_logging_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_plantation_raster, raster_path(paste0(taxa_name, "_cs_trans_plantation_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_combined_raster, raster_path(paste0(taxa_name, "_cs_trans_combined_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_fuelwood_raster, raster_path(paste0(taxa_name, "_cs_trans_fuelwood_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_industrial_roundwood_raster, raster_path(paste0(taxa_name, "_cs_trans_industrial_roundwood_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_small_scale_raster, raster_path(paste0(taxa_name, "_cs_trans_small_scale_raster.tif")), overwrite = TRUE)
  writeRaster(cs_trans_large_scale_raster, raster_path(paste0(taxa_name, "_cs_trans_large_scale_raster.tif")), overwrite = TRUE)
  
  # --- Create Maps for CS Trans ---
  map_cst1 <- create_hotspot_map(cs_trans_logging_raster, paste(taxa_name_capitalized, "- Transformed CS (Logging)"), "Transformed CS", sensitivity_colors)
  map_cst2 <- create_hotspot_map(cs_trans_plantation_raster, paste(taxa_name_capitalized, "- Transformed CS (Plantations)"), "Transformed CS", sensitivity_colors)
  map_cst3 <- create_hotspot_map(cs_trans_combined_raster, paste(taxa_name_capitalized, "- Transformed CS (All Wood)"), "Transformed CS", sensitivity_colors)
  map_cst4 <- create_hotspot_map(cs_trans_fuelwood_raster, paste(taxa_name_capitalized, "- Transformed CS (Fuelwood)"), "Transformed CS", sensitivity_colors)
  map_cst5 <- create_hotspot_map(cs_trans_industrial_roundwood_raster, paste(taxa_name_capitalized, "- Transformed CS (Ind. Roundwood)"), "Transformed CS", sensitivity_colors)
  
  # Calculate Species Richness Community Sensitivity (SR CS) ---
  cat(paste("Calculating SR Community Sensitivity (SR CS) rasters for", taxa_name, "...\n"))
  
  sr_cs_logging_raster <- ifel(sr_potential_raster > 0, logging_richness_raster / sr_potential_raster, 0)
  sr_cs_plantation_raster <- ifel(sr_potential_raster > 0, plantation_richness_raster / sr_potential_raster, 0)
  sr_cs_combined_raster <- ifel(sr_potential_raster > 0, combined_richness_raster / sr_potential_raster, 0)
  sr_cs_fuelwood_raster <- ifel(sr_potential_raster > 0, fuelwood_richness_raster / sr_potential_raster, 0)
  sr_cs_industrial_roundwood_raster <- ifel(sr_potential_raster > 0, industrial_roundwood_richness_raster / sr_potential_raster, 0)
  sr_cs_small_scale_raster <- ifel(sr_potential_raster > 0, small_scale_richness_raster / sr_potential_raster, 0)
  sr_cs_large_scale_raster <- ifel(sr_potential_raster > 0, large_scale_richness_raster / sr_potential_raster, 0)
  
  # --- SAVE SR CS RASTERS ---
  cat(paste("Saving SR Community Sensitivity rasters for", taxa_name, "...\n"))
  writeRaster(sr_cs_logging_raster, raster_path(paste0(taxa_name, "_sr_cs_logging_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_plantation_raster, raster_path(paste0(taxa_name, "_sr_cs_plantation_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_combined_raster, raster_path(paste0(taxa_name, "_sr_cs_combined_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_fuelwood_raster, raster_path(paste0(taxa_name, "_sr_cs_fuelwood_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_industrial_roundwood_raster, raster_path(paste0(taxa_name, "_sr_cs_industrial_roundwood_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_small_scale_raster, raster_path(paste0(taxa_name, "_sr_cs_small_scale_raster.tif")), overwrite = TRUE)
  writeRaster(sr_cs_large_scale_raster, raster_path(paste0(taxa_name, "_sr_cs_large_scale_raster.tif")), overwrite = TRUE)
  
  # --- 3d. Map Plotting ---
  cat("\n--- Generating all maps for", taxa_name_capitalized, "---\n")
  
  # Species Threat Maps (Standard)
  map1 <- create_hotspot_map(logging_st_raster, paste(taxa_name_capitalized, "- Threat from Logging (5.3)"), "Logging Threat", hotspot_colors)
  map2 <- create_hotspot_map(plantation_st_raster, paste(taxa_name_capitalized, "- Threat from Plantations (2.2)"), "Plantation Threat", hotspot_colors)
  map3 <- create_hotspot_map(combined_st_raster, paste(taxa_name_capitalized, "- Threat from All Wood Production"), "Wood Production Threat", hotspot_colors)
  map4 <- create_hotspot_map(fuelwood_st_raster, paste(taxa_name_capitalized, "- Threat from Fuelwood & Charcoal"), "Fuelwood Threat", hotspot_colors)
  map5 <- create_hotspot_map(industrial_roundwood_st_raster, paste(taxa_name_capitalized, "- Threat from Industrial Roundwood"), "Industrial Roundwood Threat", hotspot_colors)
  
  # Species Richness Hotspot Maps
  map_r1 <- create_hotspot_map(logging_richness_raster, paste(taxa_name_capitalized, "- Richness from Logging (5.3)"), "Species Richness", hotspot_colors)
  map_r2 <- create_hotspot_map(plantation_richness_raster, paste(taxa_name_capitalized, "- Richness from Plantations (2.2)"), "Species Richness", hotspot_colors)
  map_r3 <- create_hotspot_map(combined_richness_raster, paste(taxa_name_capitalized, "- Richness from All Wood Production"), "Species Richness", hotspot_colors)
  map_r4 <- create_hotspot_map(fuelwood_richness_raster, paste(taxa_name_capitalized, "- Richness from Fuelwood & Charcoal"), "Species Richness", hotspot_colors)
  map_r5 <- create_hotspot_map(industrial_roundwood_richness_raster, paste(taxa_name_capitalized, "- Richness from Industrial Roundwood"), "Species Richness", hotspot_colors)
  
  # --- Transformed Threat Maps (Cube Root) ---
  map_t1 <- create_hotspot_map(logging_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Logging)"), "Transformed Threat", hotspot_colors)
  map_t2 <- create_hotspot_map(plantation_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Plantations)"), "Transformed Threat", hotspot_colors)
  map_t3 <- create_hotspot_map(combined_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (All Wood)"), "Transformed Threat", hotspot_colors)
  map_t4 <- create_hotspot_map(fuelwood_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Fuelwood)"), "Transformed Threat", hotspot_colors)
  map_t5 <- create_hotspot_map(industrial_roundwood_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Ind. Roundwood)"), "Transformed Threat", hotspot_colors)
  
  # Small-scale and Large-scale Maps
  map_s1 <- create_hotspot_map(small_scale_st_raster, paste(taxa_name_capitalized, "- Threat from Small-Scale Wood Production"), "Small-Scale Threat", hotspot_colors)
  map_s2 <- create_hotspot_map(large_scale_st_raster, paste(taxa_name_capitalized, "- Threat from Large-Scale Wood Production"), "Large-Scale Threat", hotspot_colors)
  map_sr1 <- create_hotspot_map(small_scale_richness_raster, paste(taxa_name_capitalized, "- Richness from Small-Scale Wood Production"), "Small-Scale Richness", hotspot_colors)
  map_sr2 <- create_hotspot_map(large_scale_richness_raster, paste(taxa_name_capitalized, "- Richness from Large-Scale Wood Production"), "Large-Scale Richness", hotspot_colors)
  
  # --- Small/Large Scale Transformed Maps ---
  map_st1 <- create_hotspot_map(small_scale_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Small-Scale)"), "Transformed Threat", hotspot_colors)
  map_st2 <- create_hotspot_map(large_scale_trans_raster, paste(taxa_name_capitalized, "- Transformed Threat (Large-Scale)"), "Transformed Threat", hotspot_colors)
  
  # Community Sensitivity Maps (Standard)
  map6 <- create_hotspot_map(cs_logging_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Logging"), "Community Sensitivity", sensitivity_colors)
  map7 <- create_hotspot_map(cs_plantation_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Plantations"), "Community Sensitivity", sensitivity_colors)
  map8 <- create_hotspot_map(cs_combined_raster, paste(taxa_name_capitalized, "- Community Sensitivity to All Wood Production"), "Community Sensitivity", sensitivity_colors)
  map9 <- create_hotspot_map(cs_fuelwood_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Fuelwood & Charcoal"), "Community Sensitivity", sensitivity_colors)
  map10 <- create_hotspot_map(cs_industrial_roundwood_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Industrial Roundwood"), "Community Sensitivity", sensitivity_colors)
  map11 <- create_hotspot_map(cs_small_scale_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Small Scale Production"), "Community Sensitivity", sensitivity_colors)#
  map12 <- create_hotspot_map(cs_large_scale_raster, paste(taxa_name_capitalized, "- Community Sensitivity to Large Scale Production"), "Community Sensitivity", sensitivity_colors)
  
  # --- Transformed Community Sensitivity Maps ---
  map_cst1 <- create_hotspot_map(cs_trans_logging_raster, paste(taxa_name_capitalized, "- Transformed CS to Logging"), "Transformed CS", sensitivity_colors)
  map_cst2 <- create_hotspot_map(cs_trans_plantation_raster, paste(taxa_name_capitalized, "- Transformed CS to Plantations"), "Transformed CS", sensitivity_colors)
  map_cst3 <- create_hotspot_map(cs_trans_combined_raster, paste(taxa_name_capitalized, "- Transformed CS to All Wood"), "Transformed CS", sensitivity_colors)
  map_cst4 <- create_hotspot_map(cs_trans_fuelwood_raster, paste(taxa_name_capitalized, "- Transformed CS to Fuelwood"), "Transformed CS", sensitivity_colors)
  map_cst5 <- create_hotspot_map(cs_trans_industrial_roundwood_raster, paste(taxa_name_capitalized, "- Transformed CS to Ind. Roundwood"), "Transformed CS", sensitivity_colors)
  
  map_cst_small <- create_hotspot_map(cs_trans_small_scale_raster, paste(taxa_name_capitalized, "- Transformed CS (Small-Scale)"), "Transformed CS", sensitivity_colors)
  map_cst_large <- create_hotspot_map(cs_trans_large_scale_raster, paste(taxa_name_capitalized, "- Transformed CS (Large-Scale)"), "Transformed CS", sensitivity_colors)
  
  # --- NEW: SR Community Sensitivity Maps ---
  # (These were calculated but not plotted in the previous code)
  map_srcs1 <- create_hotspot_map(sr_cs_logging_raster, paste(taxa_name_capitalized, "- SR CS to Logging"), "SR Community Sensitivity", sensitivity_colors)
  map_srcs2 <- create_hotspot_map(sr_cs_plantation_raster, paste(taxa_name_capitalized, "- SR CS to Plantations"), "SR Community Sensitivity", sensitivity_colors)
  map_srcs3 <- create_hotspot_map(sr_cs_combined_raster, paste(taxa_name_capitalized, "- SR CS to All Wood"), "SR Community Sensitivity", sensitivity_colors)
  map_srcs4 <- create_hotspot_map(sr_cs_fuelwood_raster, paste(taxa_name_capitalized, "- SR CS to Fuelwood"), "SR Community Sensitivity", sensitivity_colors)
  map_srcs5 <- create_hotspot_map(sr_cs_industrial_roundwood_raster, paste(taxa_name_capitalized, "- SR CS to Ind. Roundwood"), "SR Community Sensitivity", sensitivity_colors)
  
  map_srcs_small <- create_hotspot_map(sr_cs_small_scale_raster, paste(taxa_name_capitalized, "- SR CS (Small-Scale)"), "SR Community Sensitivity", sensitivity_colors)
  map_srcs_large <- create_hotspot_map(sr_cs_large_scale_raster, paste(taxa_name_capitalized, "- SR CS (Large-Scale)"), "SR Community Sensitivity", sensitivity_colors)
  
  
  # --- SAVE AND PRINT MAPS ---
  cat(paste("\n--- Saving maps for:", taxa_name_capitalized, "---\n"))
  
  # Group 1: Logging (5.3)
  if (!is.null(map1)) { print(map1); save_map_image(map1, taxa_name, "Logging", "ST"); }
  if (!is.null(map_r1)) { print(map_r1); save_map_image(map_r1, taxa_name, "Logging", "SR"); }
  if (!is.null(map_t1)) { print(map_t1); save_map_image(map_t1, taxa_name, "Logging", "Trans"); }
  if (!is.null(map6)) { print(map6); save_map_image(map6, taxa_name, "Logging", "CS"); }
  if (!is.null(map_cst1)) { print(map_cst1); save_map_image(map_cst1, taxa_name, "Logging", "CS_Trans"); }
  if (!is.null(map_srcs1)) { print(map_srcs1); save_map_image(map_srcs1, taxa_name, "Logging", "SR_CS"); } # ADDED
  
  # Group 2: Plantations (2.2)
  if (!is.null(map2)) { print(map2); save_map_image(map2, taxa_name, "Plantations", "ST"); }
  if (!is.null(map_r2)) { print(map_r2); save_map_image(map_r2, taxa_name, "Plantations", "SR"); }
  if (!is.null(map_t2)) { print(map_t2); save_map_image(map_t2, taxa_name, "Plantations", "Trans"); }
  if (!is.null(map7)) { print(map7); save_map_image(map7, taxa_name, "Plantations", "CS"); }
  if (!is.null(map_cst2)) { print(map_cst2); save_map_image(map_cst2, taxa_name, "Plantations", "CS_Trans"); }
  if (!is.null(map_srcs2)) { print(map_srcs2); save_map_image(map_srcs2, taxa_name, "Plantations", "SR_CS"); } # ADDED
  
  # Group 3: Combined Wood Production
  if (!is.null(map3)) { print(map3); save_map_image(map3, taxa_name, "Wood_Production", "ST"); }
  if (!is.null(map_r3)) { print(map_r3); save_map_image(map_r3, taxa_name, "Wood_Production", "SR"); }
  if (!is.null(map_t3)) { print(map_t3); save_map_image(map_t3, taxa_name, "Wood_Production", "Trans"); }
  if (!is.null(map8)) { print(map8); save_map_image(map8, taxa_name, "Wood_Production", "CS"); }
  if (!is.null(map_cst3)) { print(map_cst3); save_map_image(map_cst3, taxa_name, "Wood_Production", "CS_Trans"); }
  if (!is.null(map_srcs3)) { print(map_srcs3); save_map_image(map_srcs3, taxa_name, "Wood_Production", "SR_CS"); } # ADDED
  
  # Group 4: Fuelwood
  if (!is.null(map4)) { print(map4); save_map_image(map4, taxa_name, "Fuelwood", "ST"); }
  if (!is.null(map_r4)) { print(map_r4); save_map_image(map_r4, taxa_name, "Fuelwood", "SR"); }
  if (!is.null(map_t4)) { print(map_t4); save_map_image(map_t4, taxa_name, "Fuelwood", "Trans"); }
  if (!is.null(map9)) { print(map9); save_map_image(map9, taxa_name, "Fuelwood", "CS"); }
  if (!is.null(map_cst4)) { print(map_cst4); save_map_image(map_cst4, taxa_name, "Fuelwood", "CS_Trans"); }
  if (!is.null(map_srcs4)) { print(map_srcs4); save_map_image(map_srcs4, taxa_name, "Fuelwood", "SR_CS"); } # ADDED
  
  # Group 5: Industrial Roundwood
  if (!is.null(map5)) { print(map5); save_map_image(map5, taxa_name, "Roundwood", "ST"); }
  if (!is.null(map_r5)) { print(map_r5); save_map_image(map_r5, taxa_name, "Roundwood", "SR"); }
  if (!is.null(map_t5)) { print(map_t5); save_map_image(map_t5, taxa_name, "Roundwood", "Trans"); }
  if (!is.null(map10)) { print(map10); save_map_image(map10, taxa_name, "Roundwood", "CS"); }
  if (!is.null(map_cst5)) { print(map_cst5); save_map_image(map_cst5, taxa_name, "Roundwood", "CS_Trans"); }
  if (!is.null(map_srcs5)) { print(map_srcs5); save_map_image(map_srcs5, taxa_name, "Roundwood", "SR_CS"); } # ADDED
  
  # Group 6: Small/Large Scale
  if (!is.null(map_s1)) { print(map_s1); save_map_image(map_s1, taxa_name, "Small_Scale", "ST"); }
  if (!is.null(map11)) { print(map11); save_map_image(map11, taxa_name, "Small_Scale", "CS"); }
  if (!is.null(map_sr1)) { print(map_sr1); save_map_image(map_sr1, taxa_name, "Small_Scale", "SR"); }
  if (!is.null(map_st1)) { print(map_st1); save_map_image(map_st1, taxa_name, "Small_Scale", "Trans"); }
  if (!is.null(map_cst_small)) { print(map_cst_small); save_map_image(map_cst_small, taxa_name, "Small_Scale", "CS_Trans"); }
  if (!is.null(map_srcs_small)) { print(map_srcs_small); save_map_image(map_srcs_small, taxa_name, "Small_Scale", "SR_CS"); } # ADDED
  
  if (!is.null(map_s2)) { print(map_s2); save_map_image(map_s2, taxa_name, "Large_Scale", "ST"); }
  if (!is.null(map12)) { print(map12); save_map_image(map12, taxa_name, "Large_Scale", "CS"); }
  if (!is.null(map_sr2)) { print(map_sr2); save_map_image(map_sr2, taxa_name, "Large_Scale", "SR"); }
  if (!is.null(map_st2)) { print(map_st2); save_map_image(map_st2, taxa_name, "Large_Scale", "Trans"); }
  if (!is.null(map_cst_large)) { print(map_cst_large); save_map_image(map_cst_large, taxa_name, "Large_Scale", "CS_Trans"); }
  if (!is.null(map_srcs_large)) { print(map_srcs_large); save_map_image(map_srcs_large, taxa_name, "Large_Scale", "SR_CS"); } # ADDED
  
  cat(paste("\n--- ANALYSIS FOR:", toupper(taxa_name_capitalized), "COMPLETE ---\n"))
}

# --- 4. EXECUTE ANALYSIS ---
# -----------------------------------------------------------------------------

# --- Run for Birds ---
generate_threat_maps(
  taxa_name = "bird",
  range_file = here::here("data", "processed", "ranges", "all_bird_ranges_dissolved.rds")
)

# --- Run for Mammals ---
generate_threat_maps(
  taxa_name = "mammal",
  range_file = here::here("data", "processed", "ranges", "all_mammal_ranges_dissolved.rds")
)

# --- Run for amphibians ---
generate_threat_maps(
  taxa_name = "amphibian",
  range_file = here::here("data", "processed", "ranges", "all_amphibian_ranges_dissolved.rds")
)

# --- Run for reptiles ---
generate_threat_maps(
  taxa_name = "reptile",
  range_file = here::here("data", "processed", "ranges", "all_reptile_ranges_dissolved.rds")
)

# --- Run for fish ---
generate_threat_maps(
  taxa_name = "fish",
  range_file = here::here("data", "processed", "ranges", "all_fish_ranges_dissolved.rds")
)


# --- 5. ALL VERTEBRATE ANALYSIS ---
# -----------------------------------------------------------------------------
# This section loads the cached rasters created in Section 4, sums them
# to create "All Vertebrate" maps, saves these new global rasters, 
# and plots the final summary maps.

cat("\n\n--- STARTING ANALYSIS FOR: ALL TERRESTRIAL VERTEBRATES ---\n")

# --- 5a. Load and sum Potential Threat Rasters ---
# This creates the denominator (Total VP) for the community sensitivity maps.
cat("Loading and summing Potential Threat (VP) rasters...\n")
vert_potential_raster <- rast(raster_path("bird_potential_threat_raster.tif")) +
  rast(raster_path("mammal_potential_threat_raster.tif")) +
  rast(raster_path("amphibian_potential_threat_raster.tif")) +
  rast(raster_path("reptile_potential_threat_raster.tif"))+
  rast(raster_path("fish_potential_threat_raster.tif"))

cat("Saving All Vertebrates Potential Threat raster...\n")
writeRaster(vert_potential_raster, raster_path("all_vertebrates_potential_threat_raster.tif"), overwrite = TRUE)

# [NEW] Load and sum Potential TRANSFORMED Threat rasters
cat("Loading and summing Potential Transformed Threat rasters...\n")
vert_potential_trans_raster <- rast(raster_path("bird_potential_transformed_threat_raster.tif")) + 
  rast(raster_path("mammal_potential_transformed_threat_raster.tif")) + 
  rast(raster_path("amphibian_potential_transformed_threat_raster.tif")) + 
  rast(raster_path("reptile_potential_transformed_threat_raster.tif")) + 
  rast(raster_path("fish_potential_transformed_threat_raster.tif"))

cat("Saving All Vertebrates Potential Transformed Threat raster...\n")
writeRaster(vert_potential_trans_raster, raster_path("all_vertebrates_potential_transformed_threat_raster.tif"), overwrite = TRUE)


# --- 5b. Load and sum Species Threat (ST) Rasters ---
# This creates the 5 summary "Species Threat" (VT) rasters.
cat("Loading and summing Species Threat (ST) rasters...\n")

vert_st_logging_raster <- rast(raster_path("bird_logging_threat_raster.tif")) +
  rast(raster_path("mammal_logging_threat_raster.tif")) +
  rast(raster_path("amphibian_logging_threat_raster.tif")) +
  rast(raster_path("reptile_logging_threat_raster.tif"))+
  rast(raster_path("fish_logging_threat_raster.tif"))

vert_st_plantation_raster <- rast(raster_path("bird_plantation_threat_raster.tif")) +
  rast(raster_path("mammal_plantation_threat_raster.tif")) +
  rast(raster_path("amphibian_plantation_threat_raster.tif")) +
  rast(raster_path("reptile_plantation_threat_raster.tif"))+
  rast(raster_path("fish_plantation_threat_raster.tif"))

vert_st_combined_raster <- rast(raster_path("bird_combined_threat_raster.tif")) +
  rast(raster_path("mammal_combined_threat_raster.tif")) +
  rast(raster_path("amphibian_combined_threat_raster.tif")) +
  rast(raster_path("reptile_combined_threat_raster.tif"))+
  rast(raster_path("fish_combined_threat_raster.tif"))

vert_st_fuelwood_raster <- rast(raster_path("bird_fuelwood_threat_raster.tif")) +
  rast(raster_path("mammal_fuelwood_threat_raster.tif")) +
  rast(raster_path("amphibian_fuelwood_threat_raster.tif")) +
  rast(raster_path("reptile_fuelwood_threat_raster.tif"))+
  rast(raster_path("fish_fuelwood_threat_raster.tif"))

vert_st_industrial_roundwood_raster <- rast(raster_path("bird_industrial_roundwood_threat_raster.tif")) +
  rast(raster_path("mammal_industrial_roundwood_threat_raster.tif")) +
  rast(raster_path("amphibian_industrial_roundwood_threat_raster.tif")) +
  rast(raster_path("reptile_industrial_roundwood_threat_raster.tif"))+
  rast(raster_path("fish_industrial_roundwood_threat_raster.tif"))

# --- Load and sum Small-scale and Large-scale ST Rasters ---

vert_st_small_scale_raster <- rast(raster_path("bird_small_scale_threat_raster.tif")) +
  rast(raster_path("mammal_small_scale_threat_raster.tif")) +
  rast(raster_path("amphibian_small_scale_threat_raster.tif")) +
  rast(raster_path("reptile_small_scale_threat_raster.tif")) +
  rast(raster_path("fish_small_scale_threat_raster.tif"))

vert_st_large_scale_raster <- rast(raster_path("bird_large_scale_threat_raster.tif")) +
  rast(raster_path("mammal_large_scale_threat_raster.tif")) +
  rast(raster_path("amphibian_large_scale_threat_raster.tif")) +
  rast(raster_path("reptile_large_scale_threat_raster.tif")) +
  rast(raster_path("fish_large_scale_threat_raster.tif"))
# --- SAVE ST RASTERS ---

cat("Saving All Vertebrates Species Threat (ST) rasters...\n")

writeRaster(vert_st_logging_raster, raster_path("all_vertebrates_logging_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_plantation_raster, raster_path("all_vertebrates_plantation_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_combined_raster, raster_path("all_vertebrates_combined_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_fuelwood_raster, raster_path("all_vertebrates_fuelwood_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_industrial_roundwood_raster, raster_path("all_vertebrates_industrial_roundwood_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_small_scale_raster, raster_path("all_vertebrates_small_scale_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_st_large_scale_raster, raster_path("all_vertebrates_large_scale_threat_raster.tif"), overwrite = TRUE)


# --- 5c. Load and sum Species Richness (SR) Rasters ---

cat("Loading and summing Species Richness (SR) rasters...\n")

vert_sr_logging_raster <- rast(raster_path("bird_logging_richness_raster.tif")) +
  rast(raster_path("mammal_logging_richness_raster.tif")) +
  rast(raster_path("amphibian_logging_richness_raster.tif")) +
  rast(raster_path("reptile_logging_richness_raster.tif")) +
  rast(raster_path("fish_logging_richness_raster.tif"))

vert_sr_plantation_raster <- rast(raster_path("bird_plantation_richness_raster.tif")) +
  rast(raster_path("mammal_plantation_richness_raster.tif")) +
  rast(raster_path("amphibian_plantation_richness_raster.tif")) +
  rast(raster_path("reptile_plantation_richness_raster.tif")) +
  rast(raster_path("fish_plantation_richness_raster.tif"))

vert_sr_combined_raster <- rast(raster_path("bird_combined_richness_raster.tif")) +
  rast(raster_path("mammal_combined_richness_raster.tif")) +
  rast(raster_path("amphibian_combined_richness_raster.tif")) +
  rast(raster_path("reptile_combined_richness_raster.tif")) +
  rast(raster_path("fish_combined_richness_raster.tif"))

vert_sr_fuelwood_raster <- rast(raster_path("bird_fuelwood_richness_raster.tif")) +
  rast(raster_path("mammal_fuelwood_richness_raster.tif")) +
  rast(raster_path("amphibian_fuelwood_richness_raster.tif")) +
  rast(raster_path("reptile_fuelwood_richness_raster.tif")) +
  rast(raster_path("fish_fuelwood_richness_raster.tif"))

vert_sr_industrial_roundwood_raster <- rast(raster_path("bird_industrial_roundwood_richness_raster.tif")) +
  rast(raster_path("mammal_industrial_roundwood_richness_raster.tif")) +
  rast(raster_path("amphibian_industrial_roundwood_richness_raster.tif")) +
  rast(raster_path("reptile_industrial_roundwood_richness_raster.tif")) +
  rast(raster_path("fish_industrial_roundwood_richness_raster.tif"))

# --- Load and sum Small-scale and Large-scale SR Rasters ---
vert_sr_small_scale_raster <- rast(raster_path("bird_small_scale_richness_raster.tif")) +
  rast(raster_path("mammal_small_scale_richness_raster.tif")) +
  rast(raster_path("amphibian_small_scale_richness_raster.tif")) +
  rast(raster_path("reptile_small_scale_richness_raster.tif")) +
  rast(raster_path("fish_small_scale_richness_raster.tif"))

vert_sr_large_scale_raster <- rast(raster_path("bird_large_scale_richness_raster.tif")) +
  rast(raster_path("mammal_large_scale_richness_raster.tif")) +
  rast(raster_path("amphibian_large_scale_richness_raster.tif")) +
  rast(raster_path("reptile_large_scale_richness_raster.tif")) +
  rast(raster_path("fish_large_scale_richness_raster.tif"))

# --- SAVE SR RASTERS ---
cat("Saving All Vertebrates Species Richness (SR) rasters...\n")

writeRaster(vert_sr_logging_raster, raster_path("all_vertebrates_logging_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_plantation_raster, raster_path("all_vertebrates_plantation_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_combined_raster, raster_path("all_vertebrates_combined_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_fuelwood_raster, raster_path("all_vertebrates_fuelwood_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_industrial_roundwood_raster, raster_path("all_vertebrates_industrial_roundwood_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_small_scale_raster, raster_path("all_vertebrates_small_scale_richness_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_large_scale_raster, raster_path("all_vertebrates_large_scale_richness_raster.tif"), overwrite = TRUE)

# --- 5c-2. Load and sum SR Potential (Total Richness) Rasters ---
cat("Loading and summing Species Richness Potential (Total Richness) rasters...\n")

vert_sr_potential_raster <- rast(raster_path("bird_potential_richness_raster.tif")) +
  rast(raster_path("mammal_potential_richness_raster.tif")) +
  rast(raster_path("amphibian_potential_richness_raster.tif")) +
  rast(raster_path("reptile_potential_richness_raster.tif")) +
  rast(raster_path("fish_potential_richness_raster.tif"))

cat("Saving All Vertebrates SR Potential raster...\n")

writeRaster(vert_sr_potential_raster, raster_path("all_vertebrates_sr_potential_raster.tif"), overwrite = TRUE)

# --- 5d. Load and sum Transformed Threat (Cube Root) Rasters ---
cat("Loading and summing Transformed Threat (Cube Root) rasters...\n")

vert_trans_logging_raster <- rast(raster_path("bird_logging_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_logging_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_logging_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_logging_transformed_threat_raster.tif")) +
  rast(raster_path("fish_logging_transformed_threat_raster.tif"))

vert_trans_plantation_raster <- rast(raster_path("bird_plantation_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_plantation_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_plantation_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_plantation_transformed_threat_raster.tif")) +
  rast(raster_path("fish_plantation_transformed_threat_raster.tif"))

vert_trans_combined_raster <- rast(raster_path("bird_combined_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_combined_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_combined_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_combined_transformed_threat_raster.tif")) +
  rast(raster_path("fish_combined_transformed_threat_raster.tif"))

vert_trans_fuelwood_raster <- rast(raster_path("bird_fuelwood_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_fuelwood_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_fuelwood_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_fuelwood_transformed_threat_raster.tif")) +
  rast(raster_path("fish_fuelwood_transformed_threat_raster.tif"))

vert_trans_industrial_roundwood_raster <- rast(raster_path("bird_industrial_roundwood_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_industrial_roundwood_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_industrial_roundwood_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_industrial_roundwood_transformed_threat_raster.tif")) +
  rast(raster_path("fish_industrial_roundwood_transformed_threat_raster.tif"))

# --- Load and sum Small-scale and Large-scale Transformed Rasters ---

vert_trans_small_scale_raster <- rast(raster_path("bird_small_scale_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_small_scale_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_small_scale_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_small_scale_transformed_threat_raster.tif")) +
  rast(raster_path("fish_small_scale_transformed_threat_raster.tif"))

vert_trans_large_scale_raster <- rast(raster_path("bird_large_scale_transformed_threat_raster.tif")) +
  rast(raster_path("mammal_large_scale_transformed_threat_raster.tif")) +
  rast(raster_path("amphibian_large_scale_transformed_threat_raster.tif")) +
  rast(raster_path("reptile_large_scale_transformed_threat_raster.tif")) +
  rast(raster_path("fish_large_scale_transformed_threat_raster.tif"))

# --- SAVE TRANSFORMED RASTERS ---

cat("Saving All Vertebrates Transformed Threat rasters...\n")

writeRaster(vert_trans_logging_raster, raster_path("all_vertebrates_logging_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_plantation_raster, raster_path("all_vertebrates_plantation_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_combined_raster, raster_path("all_vertebrates_combined_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_fuelwood_raster, raster_path("all_vertebrates_fuelwood_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_industrial_roundwood_raster, raster_path("all_vertebrates_industrial_roundwood_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_small_scale_raster, raster_path("all_vertebrates_small_scale_transformed_threat_raster.tif"), overwrite = TRUE)
writeRaster(vert_trans_large_scale_raster, raster_path("all_vertebrates_large_scale_transformed_threat_raster.tif"), overwrite = TRUE)

# --- 5e. Calculate "All Vertebrate" Community Sensitivity (CS) Rasters ---

cat("Calculating final Community Sensitivity (CS) rasters...\n")

vert_cs_logging_raster <- ifel(vert_potential_raster > 0, vert_st_logging_raster / vert_potential_raster, 0)
vert_cs_plantation_raster <- ifel(vert_potential_raster > 0, vert_st_plantation_raster / vert_potential_raster, 0)
vert_cs_combined_raster <- ifel(vert_potential_raster > 0, vert_st_combined_raster / vert_potential_raster, 0)
vert_cs_fuelwood_raster <- ifel(vert_potential_raster > 0, vert_st_fuelwood_raster / vert_potential_raster, 0)
vert_cs_industrial_roundwood_raster <- ifel(vert_potential_raster > 0, vert_st_industrial_roundwood_raster / vert_potential_raster, 0)
vert_cs_small_scale_raster <- ifel(vert_potential_raster > 0, vert_st_small_scale_raster / vert_potential_raster, 0)
vert_cs_large_scale_raster <- ifel(vert_potential_raster > 0, vert_st_large_scale_raster / vert_potential_raster, 0)

# --- SAVE CS RASTERS ---

cat("Saving All Vertebrates Community Sensitivity (CS) rasters...\n")

writeRaster(vert_cs_logging_raster, raster_path("all_vertebrates_cs_logging_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_plantation_raster, raster_path("all_vertebrates_cs_plantation_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_combined_raster, raster_path("all_vertebrates_cs_combined_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_fuelwood_raster, raster_path("all_vertebrates_cs_fuelwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_industrial_roundwood_raster, raster_path("all_vertebrates_cs_industrial_roundwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_small_scale_raster, raster_path("all_vertebrates_cs_small_scale_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_large_scale_raster, raster_path("all_vertebrates_cs_large_scale_raster.tif"), overwrite = TRUE)

# --- 5e-2. Calculate "All Vertebrate" Transformed Community Sensitivity (CS Trans) ---

cat("Calculating final All Vertebrates Transformed Community Sensitivity (CS Trans) rasters...\n")

# Use vert_potential_trans_raster as the denominator (Calculated in 5a)
vert_cs_trans_logging_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_logging_raster / vert_potential_trans_raster, 0)
vert_cs_trans_plantation_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_plantation_raster / vert_potential_trans_raster, 0)
vert_cs_trans_combined_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_combined_raster / vert_potential_trans_raster, 0)
vert_cs_trans_fuelwood_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_fuelwood_raster / vert_potential_trans_raster, 0)
vert_cs_trans_industrial_roundwood_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_industrial_roundwood_raster / vert_potential_trans_raster, 0)
vert_cs_trans_small_scale_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_small_scale_raster / vert_potential_trans_raster, 0)
vert_cs_trans_large_scale_raster <- ifel(vert_potential_trans_raster > 0, vert_trans_large_scale_raster / vert_potential_trans_raster, 0)

# --- SAVE CS TRANS RASTERS ---

cat("Saving All Vertebrates Transformed Community Sensitivity rasters...\n")

writeRaster(vert_cs_trans_logging_raster, raster_path("all_vertebrates_cs_trans_logging_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_plantation_raster, raster_path("all_vertebrates_cs_trans_plantation_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_combined_raster, raster_path("all_vertebrates_cs_trans_combined_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_fuelwood_raster, raster_path("all_vertebrates_cs_trans_fuelwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_industrial_roundwood_raster, raster_path("all_vertebrates_cs_trans_industrial_roundwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_small_scale_raster, raster_path("all_vertebrates_cs_trans_small_scale_raster.tif"), overwrite = TRUE)
writeRaster(vert_cs_trans_large_scale_raster, raster_path("all_vertebrates_cs_trans_large_scale_raster.tif"), overwrite = TRUE)

# --- 5g. Calculate "All Vertebrate" Species Richness Community Sensitivity (SR CS) ---

cat("Calculating final All Vertebrates SR Community Sensitivity (SR CS) rasters...\n")

# Formula: Total Threatened Species Count / Total Species Richness
vert_sr_cs_logging_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_logging_raster / vert_sr_potential_raster, 0)
vert_sr_cs_plantation_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_plantation_raster / vert_sr_potential_raster, 0)
vert_sr_cs_combined_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_combined_raster / vert_sr_potential_raster, 0)
vert_sr_cs_fuelwood_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_fuelwood_raster / vert_sr_potential_raster, 0)
vert_sr_cs_industrial_roundwood_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_industrial_roundwood_raster / vert_sr_potential_raster, 0)
vert_sr_cs_small_scale_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_small_scale_raster / vert_sr_potential_raster, 0)
vert_sr_cs_large_scale_raster <- ifel(vert_sr_potential_raster > 0, vert_sr_large_scale_raster / vert_sr_potential_raster, 0)

# --- SAVE SR CS RASTERS ---

cat("Saving All Vertebrates SR Community Sensitivity rasters...\n")

writeRaster(vert_sr_cs_logging_raster, raster_path("all_vertebrates_sr_cs_logging_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_plantation_raster, raster_path("all_vertebrates_sr_cs_plantation_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_combined_raster, raster_path("all_vertebrates_sr_cs_combined_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_fuelwood_raster, raster_path("all_vertebrates_sr_cs_fuelwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_industrial_roundwood_raster, raster_path("all_vertebrates_sr_cs_industrial_roundwood_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_small_scale_raster, raster_path("all_vertebrates_sr_cs_small_scale_raster.tif"), overwrite = TRUE)
writeRaster(vert_sr_cs_large_scale_raster, raster_path("all_vertebrates_sr_cs_large_scale_raster.tif"), overwrite = TRUE)

# --- 5f. Map Plotting for All Vertebrates ---
cat("\n--- Generating all final maps for All Vertebrates ---\n")

# Species Threat Maps
map_v1 <- create_hotspot_map(vert_st_logging_raster, "All Vertebrates - Threat from Logging (5.3)", "Logging Threat", hotspot_colors)
map_v2 <- create_hotspot_map(vert_st_plantation_raster, "All Vertebrates - Threat from Plantations (2.2)", "Plantation Threat", hotspot_colors)
map_v3 <- create_hotspot_map(vert_st_combined_raster, "All Vertebrates - Threat from All Wood Production", "Wood Production Threat", hotspot_colors)
map_v4 <- create_hotspot_map(vert_st_fuelwood_raster, "All Vertebrates - Threat from Fuelwood & Charcoal", "Fuelwood Threat", hotspot_colors)
map_v5 <- create_hotspot_map(vert_st_industrial_roundwood_raster, "All Vertebrates - Threat from Industrial Roundwood", "Industrial Roundwood Threat", hotspot_colors)

# Species Richness Hotspot Maps
map_vr1 <- create_hotspot_map(vert_sr_logging_raster, "All Vertebrates - Richness from Logging (5.3)", "Species Richness", hotspot_colors)
map_vr2 <- create_hotspot_map(vert_sr_plantation_raster, "All Vertebrates - Richness from Plantations (2.2)", "Species Richness", hotspot_colors)
map_vr3 <- create_hotspot_map(vert_sr_combined_raster, "All Vertebrates - Richness from All Wood Production", "Species Richness", hotspot_colors)
map_vr4 <- create_hotspot_map(vert_sr_fuelwood_raster, "All Vertebrates - Richness from Fuelwood & Charcoal", "Species Richness", hotspot_colors)
map_vr5 <- create_hotspot_map(vert_sr_industrial_roundwood_raster, "All Vertebrates - Richness from Industrial Roundwood", "Species Richness", hotspot_colors)

# Transformed Threat Maps (Cube Root)
map_vt1 <- create_hotspot_map(vert_trans_logging_raster, "All Vertebrates - Transformed Threat (Logging)", "Transformed Threat", hotspot_colors)
map_vt2 <- create_hotspot_map(vert_trans_plantation_raster, "All Vertebrates - Transformed Threat (Plantations)", "Transformed Threat", hotspot_colors)
map_vt3 <- create_hotspot_map(vert_trans_combined_raster, "All Vertebrates - Transformed Threat (All Wood)", "Transformed Threat", hotspot_colors)
map_vt4 <- create_hotspot_map(vert_trans_fuelwood_raster, "All Vertebrates - Transformed Threat (Fuelwood)", "Transformed Threat", hotspot_colors)
map_vt5 <- create_hotspot_map(vert_trans_industrial_roundwood_raster, "All Vertebrates - Transformed Threat (Ind. Roundwood)", "Transformed Threat", hotspot_colors)

# --- NEW: Transformed Community Sensitivity Maps (CS Trans) ---
map_v_cst1 <- create_hotspot_map(vert_cs_trans_logging_raster, "All Vertebrates - Transformed CS (Logging)", "Transformed CS", sensitivity_colors)
map_v_cst2 <- create_hotspot_map(vert_cs_trans_plantation_raster, "All Vertebrates - Transformed CS (Plantations)", "Transformed CS", sensitivity_colors)
map_v_cst3 <- create_hotspot_map(vert_cs_trans_combined_raster, "All Vertebrates - Transformed CS (All Wood)", "Transformed CS", sensitivity_colors)
map_v_cst4 <- create_hotspot_map(vert_cs_trans_fuelwood_raster, "All Vertebrates - Transformed CS (Fuelwood)", "Transformed CS", sensitivity_colors)
map_v_cst5 <- create_hotspot_map(vert_cs_trans_industrial_roundwood_raster, "All Vertebrates - Transformed CS (Ind. Roundwood)", "Transformed CS", sensitivity_colors)
map_v_cst_small <- create_hotspot_map(vert_cs_trans_small_scale_raster, "All Vertebrates - Transformed CS (Small-Scale)", "Transformed CS", sensitivity_colors)
map_v_cst_large <- create_hotspot_map(vert_cs_trans_large_scale_raster, "All Vertebrates - Transformed CS (Large-Scale)", "Transformed CS", sensitivity_colors)

# Small-scale and Large-scale Maps
map_vs1 <- create_hotspot_map(vert_st_small_scale_raster, "All Vertebrates - Threat from Small-Scale Wood Production", "Small-Scale Threat", hotspot_colors)
map_vs2 <- create_hotspot_map(vert_st_large_scale_raster, "All Vertebrates - Threat from Large-Scale Wood Production", "Large-Scale Threat", hotspot_colors)
map_vsr1 <- create_hotspot_map(vert_sr_small_scale_raster, "All Vertebrates - Richness from Small-Scale Wood Production", "Small-Scale Richness", hotspot_colors)
map_vsr2 <- create_hotspot_map(vert_sr_large_scale_raster, "All Vertebrates - Richness from Large-Scale Wood Production", "Large-Scale Richness", hotspot_colors)
map_v_cs_small <- create_hotspot_map(vert_cs_small_scale_raster, "All Vertebrates - Community Sensitivity (Small-Scale)", "Community Sensitivity", sensitivity_colors)
map_v_cs_large <- create_hotspot_map(vert_cs_large_scale_raster, "All Vertebrates - Community Sensitivity (Large-Scale)", "Community Sensitivity", sensitivity_colors)
map_vst1 <- create_hotspot_map(vert_trans_small_scale_raster, "All Vertebrates - Transformed Threat (Small-Scale)", "Transformed Threat", hotspot_colors)
map_vst2 <- create_hotspot_map(vert_trans_large_scale_raster, "All Vertebrates - Transformed Threat (Large-Scale)", "Transformed Threat", hotspot_colors)

# Community Sensitivity Maps
map_v6 <- create_hotspot_map(vert_cs_logging_raster, "All Vertebrates - Community Sensitivity to Logging", "Community Sensitivity", sensitivity_colors)
map_v7 <- create_hotspot_map(vert_cs_plantation_raster, "All Vertebrates - Community Sensitivity to Plantations", "Community Sensitivity", sensitivity_colors)
map_v8 <- create_hotspot_map(vert_cs_combined_raster, "All Vertebrates - Community Sensitivity to All Wood Production", "Community Sensitivity", sensitivity_colors)
map_v9 <- create_hotspot_map(vert_cs_fuelwood_raster, "All Vertebrates - Community Sensitivity to Fuelwood & Charcoal", "Community Sensitivity", sensitivity_colors)
map_v10 <- create_hotspot_map(vert_cs_industrial_roundwood_raster, "All Vertebrates - Community Sensitivity to Industrial Roundwood", "Community Sensitivity", sensitivity_colors)

# SR Community Sensitivity Maps
map_v_srcs1 <- create_hotspot_map(vert_sr_cs_logging_raster, "All Vertebrates - SR Community Sensitivity (Logging)", "SR Sensitivity", sensitivity_colors)
map_v_srcs2 <- create_hotspot_map(vert_sr_cs_plantation_raster, "All Vertebrates - SR Community Sensitivity (Plantations)", "SR Sensitivity", sensitivity_colors)
map_v_srcs3 <- create_hotspot_map(vert_sr_cs_combined_raster, "All Vertebrates - SR Community Sensitivity (All Wood)", "SR Sensitivity", sensitivity_colors)
map_v_srcs4 <- create_hotspot_map(vert_sr_cs_fuelwood_raster, "All Vertebrates - SR Community Sensitivity (Fuelwood)", "SR Sensitivity", sensitivity_colors)
map_v_srcs5 <- create_hotspot_map(vert_sr_cs_industrial_roundwood_raster, "All Vertebrates - SR Community Sensitivity (Ind. Roundwood)", "SR Sensitivity", sensitivity_colors)
map_v_srcs_small <- create_hotspot_map(vert_sr_cs_small_scale_raster, "All Vertebrates - SR Community Sensitivity (Small-Scale)", "SR Sensitivity", sensitivity_colors)
map_v_srcs_large <- create_hotspot_map(vert_sr_cs_large_scale_raster, "All Vertebrates - SR Community Sensitivity (Large-Scale)", "SR Sensitivity", sensitivity_colors)

# --- SAVE AND PRINT GLOBAL MAPS ---
cat("\n--- Saving final maps for All Vertebrates ---\n")
all_taxa_name <- "All_Vertebrates"

# Group 1: Logging/ Plantations
if (!is.null(map_v1)) { print(map_v1); save_map_image(map_v1, all_taxa_name, "Logging", "ST"); }
if (!is.null(map_v6)) { print(map_v6); save_map_image(map_v6, all_taxa_name, "Logging", "CS"); }
if (!is.null(map_vr1)) { print(map_vr1); save_map_image(map_vr1, all_taxa_name, "Logging", "SR"); }
if (!is.null(map_v_srcs1)) { print(map_v_srcs1); save_map_image(map_v_srcs1, all_taxa_name, "Logging", "SR_CS"); }
if (!is.null(map_vt1)) { print(map_vt1); save_map_image(map_vt1, all_taxa_name, "Logging", "Trans"); }
if (!is.null(map_v_cst1)) { print(map_v_cst1); save_map_image(map_v_cst1, all_taxa_name, "Logging", "CS_Trans"); }

if (!is.null(map_v2)) { print(map_v2); save_map_image(map_v2, all_taxa_name, "Plantations", "ST"); }
if (!is.null(map_vr2)) { print(map_vr2); save_map_image(map_vr2, all_taxa_name, "Plantations", "SR"); }
if (!is.null(map_v_srcs2)) { print(map_v_srcs2); save_map_image(map_v_srcs2, all_taxa_name, "Plantations", "SR_CS"); }
if (!is.null(map_vt2)) { print(map_vt2); save_map_image(map_vt2, all_taxa_name, "Plantations", "Trans"); }
if (!is.null(map_v_cst2)) { print(map_v_cst2); save_map_image(map_v_cst2, all_taxa_name, "Plantations", "CS_Trans"); }
if (!is.null(map_v7)) { print(map_v7); save_map_image(map_v7, all_taxa_name, "Plantations", "CS"); }

# Group 2: Combined Wood Production
if (!is.null(map_v3)) { print(map_v3); save_map_image(map_v3, all_taxa_name, "Wood_Production", "ST"); }
if (!is.null(map_vr3)) { print(map_vr3); save_map_image(map_vr3, all_taxa_name, "Wood_Production", "SR"); }
if (!is.null(map_v_srcs3)) { print(map_v_srcs3); save_map_image(map_v_srcs3, all_taxa_name, "Wood_Production", "SR_CS"); }
if (!is.null(map_vt3)) { print(map_vt3); save_map_image(map_vt3, all_taxa_name, "Wood_Production", "Trans"); }
if (!is.null(map_v_cst3)) { print(map_v_cst3); save_map_image(map_v_cst3, all_taxa_name, "Wood_Production", "CS_Trans"); }
if (!is.null(map_v8)) { print(map_v8); save_map_image(map_v8, all_taxa_name, "Wood_Production", "CS"); }

# Group 3: Fuelwood/Industrial Roundwood
if (!is.null(map_v4)) { print(map_v4); save_map_image(map_v4, all_taxa_name, "Fuelwood", "ST"); }
if (!is.null(map_vr4)) { print(map_vr4); save_map_image(map_vr4, all_taxa_name, "Fuelwood", "SR"); }
if (!is.null(map_v_srcs4)) { print(map_v_srcs4); save_map_image(map_v_srcs4, all_taxa_name, "Fuelwood", "SR_CS"); }
if (!is.null(map_vt4)) { print(map_vt4); save_map_image(map_vt4, all_taxa_name, "Fuelwood", "Trans"); }
if (!is.null(map_v_cst4)) { print(map_v_cst4); save_map_image(map_v_cst4, all_taxa_name, "Fuelwood", "CS_Trans"); }
if (!is.null(map_v9)) { print(map_v9); save_map_image(map_v9, all_taxa_name, "Fuelwood", "CS"); }

if (!is.null(map_v5)) { print(map_v5); save_map_image(map_v5, all_taxa_name, "Roundwood", "ST"); }
if (!is.null(map_vr5)) { print(map_vr5); save_map_image(map_vr5, all_taxa_name, "Roundwood", "SR"); }
if (!is.null(map_v_srcs5)) { print(map_v_srcs5); save_map_image(map_v_srcs5, all_taxa_name, "Roundwood", "SR_CS"); }
if (!is.null(map_vt5)) { print(map_vt5); save_map_image(map_vt5, all_taxa_name, "Roundwood", "Trans"); }
if (!is.null(map_v_cst5)) { print(map_v_cst5); save_map_image(map_v_cst5, all_taxa_name, "Roundwood", "CS_Trans"); }
if (!is.null(map_v10)) { print(map_v10); save_map_image(map_v10, all_taxa_name, "Roundwood", "CS"); }

# Group 4: Small/ Large Scale
if (!is.null(map_vs1)) { print(map_vs1); save_map_image(map_vs1, all_taxa_name, "Small_Scale", "ST"); }
if (!is.null(map_v_cs_small)) { print(map_v_cs_small); save_map_image(map_v_cs_small, all_taxa_name, "Small_Scale", "CS"); }
if (!is.null(map_vsr1)) { print(map_vsr1); save_map_image(map_vsr1, all_taxa_name, "Small_Scale", "SR"); }
if (!is.null(map_v_srcs_small)) { print(map_v_srcs_small); save_map_image(map_v_srcs_small, all_taxa_name, "Small_Scale", "SR_CS"); }
if (!is.null(map_vst1)) { print(map_vst1); save_map_image(map_vst1, all_taxa_name, "Small_Scale", "Trans"); }
if (!is.null(map_v_cst_small)) { print(map_v_cst_small); save_map_image(map_v_cst_small, all_taxa_name, "Small_Scale", "CS_Trans"); }

if (!is.null(map_vs2)) { print(map_vs2); save_map_image(map_vs2, all_taxa_name, "Large_Scale", "ST"); }
if (!is.null(map_v_cs_large)) { print(map_v_cs_large); save_map_image(map_v_cs_large, all_taxa_name, "Large_Scale", "CS"); }
if (!is.null(map_vsr2)) { print(map_vsr2); save_map_image(map_vsr2, all_taxa_name, "Large_Scale", "SR"); }
if (!is.null(map_v_srcs_large)) { print(map_v_srcs_large); save_map_image(map_v_srcs_large, all_taxa_name, "Large_Scale", "SR_CS"); }
if (!is.null(map_vst2)) { print(map_vst2); save_map_image(map_vst2, all_taxa_name, "Large_Scale", "Trans"); }
if (!is.null(map_v_cst_large)) { print(map_v_cst_large); save_map_image(map_v_cst_large, all_taxa_name, "Large_Scale", "CS_Trans"); }

cat("\n--- ALL VERTEBRATE ANALYSIS COMPLETE ---\n")

# --- 6. BIVARIATE PLOTS ---
# -----------------------------------------------------------------------------
# This section creates bivariate plots that compare hotspots of logging vs
# hotspots of planations, small scale vs large scale production, and fuelwood
# vs industrial roundwood. The biscale package is useful for creating these
# colour palettes.

# --- Define Global Bivariate Parameters for 4x4 Classification ---
# New breaks: <90%, 90-95%, 95-99%, >99%
BI_N_CLASSES <- 4
BI_QUANTILES <- c(0.90, 0.95, 0.99)

# --- Helper Function ---
create_biscale_map <- function(
    raster_x, raster_y, x_name, y_name,
    map_title, map_subtitle,
    palette_name) { # Pass the string name of the palette (e.g., "DkBlue2")
  
  cat(paste0("\n--- Generating 4x4 Biscale Map: ", palette_name, " (", x_name, " vs. ", y_name, ") ---\n"))
  
  # 1. Stack and mask
  biv_stack <- c(raster_x, raster_y)
  names(biv_stack) <- c("x_var", "y_var")
  biv_stack_land <- terra::mask(biv_stack, vect(world_map_moll))
  
  # 2. Convert to dataframe and filter for cells with *any* threat
  biv_df <- as.data.frame(biv_stack_land, xy = TRUE) %>%
    filter(x_var > 0 | y_var > 0)
  
  if (nrow(biv_df) == 0) {
    cat(" -> No overlapping data found.\n")
    return(NULL)
  }
  
  # 3. Calculate Custom Breaks (90, 95, 99)
  # Calculate these on the positive values only
  x_pos <- biv_df$x_var[biv_df$x_var > 0]
  y_pos <- biv_df$y_var[biv_df$y_var > 0]
  
  x_qt <- quantile(x_pos, BI_QUANTILES, na.rm = TRUE)
  y_qt <- quantile(y_pos, BI_QUANTILES, na.rm = TRUE)
  
  # Create break vectors for cut(): c(-Inf, 90%, 95%, 99%, Inf)
  x_breaks_full <- c(-Inf, x_qt, Inf)
  y_breaks_full <- c(-Inf, y_qt, Inf)
  
  # 4. Manually create Factors
  # As per guide: To use custom breaks, pass factors to bi_class
  biv_df <- biv_df %>%
    mutate(
      x_factor = cut(x_var, breaks = x_breaks_full, labels = FALSE, include.lowest = TRUE),
      y_factor = cut(y_var, breaks = y_breaks_full, labels = FALSE, include.lowest = TRUE),
      x_factor = factor(x_factor, levels = 1:4),
      y_factor = factor(y_factor, levels = 1:4)
    )
  
  # 5. Run bi_class
  # Pass the factors; style is ignored when x and y are factors
  biv_data <- bi_class(biv_df, x = x_factor, y = y_factor, dim = BI_N_CLASSES)
  
  # 6. create Legend Labels
  # Manually construct the list structure that bi_legend expects for 'breaks'
  # These will appear as text labels on the legend keys
  legend_labels <- c(
    paste0("<", BI_QUANTILES[1]*100, "%"),
    paste0(BI_QUANTILES[1]*100, "-", BI_QUANTILES[2]*100, "%"),
    paste0(BI_QUANTILES[2]*100, "-", BI_QUANTILES[3]*100, "%"),
    paste0(">", BI_QUANTILES[3]*100, "%")
  )
  
  custom_breaks <- list(
    bi_x = legend_labels,
    bi_y = legend_labels
  )
  
  # 7. Create Legend
  # Note: Use the "2" version of the palette for dim=4 support
  legend <- bi_legend(
    pal = palette_name,
    dim = BI_N_CLASSES,
    xlab = x_name,
    ylab = y_name,
    breaks = custom_breaks, # Pass our manual label list
    size = 8
  )
  
  # 8. Create Map
  # bi_scale_fill applies the palette color mapping
  map_plot <- ggplot() +
    geom_sf(data = world_map_moll, fill = "gainsboro", color = "white") +
    geom_tile(data = biv_data, aes(x = x, y = y, fill = bi_class)) +
    bi_scale_fill(pal = palette_name, dim = BI_N_CLASSES) +
    theme_cowplot(10)+
    
    labs(
      title = map_title,
      subtitle = map_subtitle
    ) +
    coord_sf(crs = moll_crs) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) + # Hide default legend
    annotation_custom(
      grob = ggplotGrob(legend),
      xmin = -18040000, xmax = -9040000,
      ymin = -9020000,  ymax = -2520000
    )
  
  print(map_plot)
  return(map_plot)
}

# --- Generate the 3 Maps ---

# 1. Logging vs. Plantations
# Palette: "DkBlue2" (Supports 4x4)
map_1 <- create_biscale_map(
  raster_x = rast(raster_path("all_vertebrates_logging_threat_raster.tif")),
  raster_y = rast(raster_path("all_vertebrates_plantation_threat_raster.tif")),
  x_name = "Logging",
  y_name = "Plantation",
  map_title = "Logging vs. Plantations (4x4 Biscale)",
  map_subtitle = "DkBlue2 Palette | Breaks: 90th, 95th, 99th percentiles",
  palette_name = "DkBlue2" 
)

save_map_image(map_1, all_taxa_name, "logging_vs_plantation", "linear")

# 2. Roundwood vs. Fuelwood
# Palette: "GrPink2" (Supports 4x4)
map_2 <- create_biscale_map(
  raster_x = rast(raster_path("all_vertebrates_industrial_roundwood_threat_raster.tif")),
  raster_y = rast(raster_path("all_vertebrates_fuelwood_threat_raster.tif")),
  x_name = "Ind. Roundwood",
  y_name = "Fuelwood",
  map_title = "Roundwood vs. Fuelwood (4x4 Biscale)",
  map_subtitle = "GrPink2 Palette | Breaks: 90th, 95th, 99th percentiles",
  palette_name = "GrPink2"
)

save_map_image(map_2, all_taxa_name, "roundwood_vs_fuelwood", "linear")

# 3. Small-scale vs. Large-scale
# Palette: "DkViolet2" (Supports 4x4)
map_3 <- create_biscale_map(
  raster_x = rast(raster_path("all_vertebrates_small_scale_richness_raster.tif")),
  raster_y = rast(raster_path("all_vertebrates_large_scale_richness_raster.tif")),
  x_name = "Small-Scale",
  y_name = "Large-Scale",
  map_title = "Small vs. Large Scale (4x4 Biscale)",
  map_subtitle = "DkViolet2 Palette | Breaks: 90th, 95th, 99th percentiles",
  palette_name = "DkViolet2"
)

save_map_image(map_3, all_taxa_name, "small_vs_large_scale", "richness")

# --- 7. Combining maps ---
# -----------------------------------------------------------------------------
# Use simple addition (+) to combine plots
# plot_layout() handles the grid and the shared legend
map_type="_combined_richness_raster.tif"
map_type2="_cs_sr_raster.tif"

#species threat maps
raster_name1 = paste0("all_vertebrates",map_type)
raster1=rast(raster_path(raster_name1))
p1= create_continuous_hotspot_map2(raster1, raster_name1, map_type, transform_scale = "sqrt")
print(p1)

raster_name2 = paste0("bird",map_type)
raster2=rast(raster_path(raster_name2))
p2= create_continuous_hotspot_map2(raster2, raster_name2, map_type, transform_scale = "sqrt")
print(p2)

raster_name3 = paste0("mammal",map_type)
raster3=rast(raster_path(raster_name3))
p3= create_continuous_hotspot_map2(raster3, raster_name3, map_type, transform_scale = "sqrt")
print(p3)

raster_name4 = paste0("amphibian",map_type)
raster4=rast(raster_path(raster_name4))
p4= create_continuous_hotspot_map2(raster4, raster_name4, map_type, transform_scale = "sqrt")
print(p4)

raster_name5 = paste0("reptile",map_type)
raster5=rast(raster_path(raster_name5))
p5= create_continuous_hotspot_map2(raster5, raster_name5, map_type, transform_scale = "sqrt")
print(p5)

raster_name6 = paste0("fish",map_type)
raster6=rast(raster_path(raster_name6))
p6= create_continuous_hotspot_map2(raster6, raster_name6, map_type, transform_scale = "sqrt")
print(p6)

#community threat maps
raster_name7 = paste0("all_vertebrates",map_type2)
raster7=rast(raster_path(raster_name7))
p7= create_continuous_hotspot_map(raster7, raster_name7, map_type2)
print(p7)

raster_name8 = paste0("bird",map_type2)
raster8=rast(raster_path(raster_name8))
p8= create_continuous_hotspot_map(raster8, raster_name8, map_type2)
print(p8)

raster_name9 = paste0("mammal",map_type2)
raster9=rast(raster_path(raster_name9))
p9= create_continuous_hotspot_map(raster9, raster_name9, map_type2)
print(p9)

raster_name10 = paste0("amphibian",map_type2)
raster10=rast(raster_path(raster_name10))
p10= create_continuous_hotspot_map(raster10, raster_name10, map_type2)
print(p10)

raster_name11 = paste0("reptile",map_type2)
raster11=rast(raster_path(raster_name11))
p11= create_continuous_hotspot_map(raster11, raster_name11, map_type2)
print(p11)

raster_name12 = paste0("fish",map_type2)
raster12=rast(raster_path(raster_name12))
p12= create_continuous_hotspot_map(raster12, raster_name12, map_type2)
print(p12)

library(patchwork)

combined_plot <- (p1 + p7 + p2 + p8 + p3 + p9 + p4 + p10 + p5 + p11 + p6 + p12) +
  plot_layout(
    ncol = 2, 
    nrow = 6
  ) +
  plot_annotation(
    caption = 'Top to bottom: all vertebrates, birds, mammals, amphibians, reptiles, fish
    Left to right: species threat, community threat',
  ) & 
  theme(legend.position = "right",
        legend.title=element_blank(),
        plot.title = element_blank())

# Print to screen to verify
print(combined_plot)

# Save the plot
ggsave(
  filename = file.path(summary_dir, "combined_wood_production_transformed_continuous.png"), 
  plot = combined_plot, 
  width = 9, 
  height = 12, 
  dpi = 300,
  bg = "white" # Ensures transparent backgrounds don't turn black
)

message("Plot saved successfully as 'world_maps_grid.png'")


# =============================================================================
# 8. Realm Scatter Plots for All Taxa
# =============================================================================
# This section creates scatter plots colour coded by zoogeographic realm
# -----------------------------------------------------------------------------
# 1. GLOBAL SETUP (Run Once)
# -----------------------------------------------------------------------------

  
  # Save
  clean_name <- gsub(" ", "_", display_name)
  outfile <- file.path(realm_dir, paste0("Scatter_", clean_name, "_vs_Transformed.png"))
  ggsave(outfile, plot = p, width = 10, height = 7, dpi = 300)
  
  cat(paste0("  -> Saved: ", outfile, "\n"))


# -----------------------------------------------------------------------------
# 3. EXECUTION LOOP
# -----------------------------------------------------------------------------
cat("\n--- 3. EXECUTING BATCH ANALYSIS ---\n")

# Define your groups here. 
# Left side: File prefix (e.g., 'birds' -> 'birds_sr_potential_raster.tif')
# Right side: Display name for the plot title
taxa_groups <- list(
  "all_vertebrates" = "All Vertebrates",
  "bird"           = "Birds",
  "mammal"         = "Mammals",
  "amphibian"      = "Amphibians",
  "reptile"        = "Reptiles",
  "fish"            = "Freshwater Fish" # Check if your file is 'fish_...' or 'freshwater_fish_...'
)

# Iterate through the list
for (prefix in names(taxa_groups)) {
  create_realm_scatter(group_prefix = prefix, display_name = taxa_groups[[prefix]])
}

cat("\n--- Batch Analysis Complete ---\n")

# =============================================================================
# 9. Comparative Scatter Plots (9 Plots per Taxon)
# =============================================================================
# Scatter plots to complement the bivariate plots created in section 7
# -----------------------------------------------------------------------------
# 1. GLOBAL SETUP (Run Once)
# -----------------------------------------------------------------------------
cat("\n--- 1. GLOBAL SETUP: Pre-processing Shapefile ---\n")

# Define Colors
my_colors <- c(
  "Afrotropical"   = "#FF6B6B",
  "Australian"     = "#D4A017",  
  "Madagascan"     = "#FFD54F",  
  "Nearctic"       = "#55A630",
  "Neotropical"    = "#00C853",
  "Oceanina"       = "#00BFA5",  
  "Oriental"       = "#00B0FF",  
  "Palearctic"     = "#2979FF",
  "Panamanian"     = "#6200EA",  
  "Saharo-Arabian" = "#D500F9",  
  "Sino-Japanese"  = "#C51162"   
)

# Load and Fix Shapefile
realms_path <- here::here("data", "raw", "realms", "newRealms.shp")

if (file.exists(realms_path)) {
  suppressMessages(sf::sf_use_s2(FALSE))
  realms_global_sf <- st_read(realms_path, quiet = TRUE)
  
  # Fix Date Line Wrapping (Antimeridian)
  realms_global_sf <- st_make_valid(realms_global_sf)
  realms_global_sf <- st_wrap_dateline(realms_global_sf, options = c("WRAPDATELINE=YES"))
  realms_global_sf <- st_make_valid(realms_global_sf)
  
  # Auto-detect Realm Column
  potential_cols <- names(realms_global_sf)[sapply(realms_global_sf, function(x) length(unique(x)) > 4 & length(unique(x)) < 20)]
  realm_col_name <- grep("Realm", potential_cols, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(realm_col_name)) realm_col_name <- potential_cols[1]
  
  cat(paste0("  -> Shapefile ready. Using column: '", realm_col_name, "'\n"))
} else {
  stop("CRITICAL: Realm shapefile not found.")
}

# -----------------------------------------------------------------------------
# 2. HELPER FUNCTION: PREPARE MASTER REALM RASTER (Per Taxon)
# -----------------------------------------------------------------------------
# This creates a single, gap-filled realm raster aligned to the taxon's grid
prepare_taxon_realm_raster <- function(template_raster, land_mask_vect, realms_sf, r_col) {
  
  # 1. Mask template to land (defines the grid)
  grid_masked <- terra::mask(template_raster, land_mask_vect)
  
  # 2. Project & Rasterize Realms
  realms_proj <- st_transform(realms_sf, crs(grid_masked))
  realms_proj <- st_make_valid(realms_proj)
  realms_proj$rid <- as.numeric(as.factor(realms_proj[[r_col]]))
  
  realm_r <- terra::rasterize(realms_proj, grid_masked, field = "rid", touches = TRUE)
  
  # 3. Gap Filling (Fill ALL land pixels in the mask)
  # Use the grid_masked itself as the target mask
  target_mask <- grid_masked
  target_mask[!is.na(target_mask)] <- 1 # Binary mask of land
  
  prev_missing <- Inf
  for(i in 1:30) {
    missing_mask <- !is.na(target_mask) & is.na(realm_r)
    missing_count <- global(missing_mask, "sum", na.rm=TRUE)[1,1]
    
    if (missing_count == 0) break
    
    win_size <- if(i <= 5) 3 else if(i <= 15) 5 else 7
    filled <- terra::focal(realm_r, w=win_size, fun="modal", na.policy="only", na.rm=TRUE)
    realm_r <- terra::cover(realm_r, filled)
    
    prev_missing <- missing_count
  }
  
  # Return the lookup table and the raster
  lookup <- unique(st_drop_geometry(realms_proj)[, c("rid", r_col)])
  colnames(lookup) <- c("realm_id", "Realm_Name")
  
  return(list(raster = realm_r, lookup = lookup))
}

# -----------------------------------------------------------------------------
# 3. HELPER FUNCTION: GENERATE ONE PLOT
# -----------------------------------------------------------------------------
generate_comparison_plot <- function(x_raster, y_raster, realm_raster, lookup, 
                                     x_name, y_name, measure_label, 
                                     title_prefix, filename_suffix, colors) {
  
  # Assemble Data
  stack <- c(x_raster, y_raster, realm_raster)
  names(stack) <- c("x_val", "y_val", "realm_id")
  df <- as.data.frame(stack, na.rm = FALSE)
  
  # Filter: Keep if at least one threat > 0
  df <- df %>% 
    filter(!is.na(x_val), !is.na(y_val)) %>%
    filter(x_val > 0 | y_val > 0) %>%
    left_join(lookup, by = "realm_id") %>%
    mutate(Realm_Name = ifelse(is.na(Realm_Name), "Unknown", Realm_Name)) %>%
    filter(Realm_Name != "Unknown")
  
  if (nrow(df) == 0) return(NULL)
  
  # Calculate Equal Axis Limits
  max_val <- max(c(df$x_val, df$y_val), na.rm = TRUE)
  
  # Plot
  p <- ggplot(df, aes(x = x_val, y = y_val)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60") +
    geom_point(aes(color = Realm_Name), alpha = 1, size = 1) +
    #geom_smooth(method = "gam", color = "black", size = 1.2, se = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, max_val*1.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, max_val*1.05)) +
    scale_color_manual(values = colors) +
    coord_fixed(ratio = 1) +
    labs(
      x = x_name,
      y = y_name,
      color = "Realm",
      title = paste0(title_prefix, ": ", measure_label),
      subtitle = paste0(x_name, " vs. ", y_name)
    ) +
    theme_minimal() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      plot.title = element_text(face = "bold")
    )
  
  # Save
  outfile <- file.path(realm_dir, paste0("Scatter_", filename_suffix, ".png"))
  ggsave(outfile, plot = p, width = 10, height = 7, dpi = 300)
  cat(paste0("      -> Saved: ", basename(outfile), "\n"))
}

# -----------------------------------------------------------------------------
# 4. MAIN EXECUTION LOOP
# -----------------------------------------------------------------------------
cat("\n--- 4. EXECUTING BATCH ANALYSIS (9 Maps per Taxon) ---\n")

# A. Define Lists
taxa_groups <- list(
  "all_vertebrates" = "All Vertebrates",
  "bird"           = "Birds",
  "mammal"         = "Mammals",
  "amphibian"      = "Amphibians",
  "reptile"        = "Reptiles",
  "fish"            = "Freshwater Fish"
)

# Define the 3 Comparisons (X axis file tag, Y axis file tag, X Label, Y Label)
comparisons <- list(
  list(x="logging", y="plantation", x_lab="Logging", y_lab="Plantations"),
  list(x="industrial_roundwood", y="fuelwood", x_lab="Ind. Roundwood", y_lab="Fuelwood"),
  list(x="small_scale", y="large_scale", x_lab="Small-scale", y_lab="Large-scale")
)

# Define the 3 Measures (File tag, Display Label)
measures <- list(
  list(tag="richness", label="Species Richness"),
  list(tag="transformed_threat", label="Transformed Threat"),
  list(tag="threat", label="Linear Threat")
)

# Check Land Mask
if (!exists("world_map_moll")) stop("Error: 'world_map_moll' not loaded.")
land_vect <- vect(world_map_moll)

# --- B. Taxon Loop ---
for (taxa_prefix in names(taxa_groups)) {
  
  display_name <- taxa_groups[[taxa_prefix]]
  cat(paste0("\n>>> PROCESSING TAXON: ", display_name, " <<<\n"))
  
  # 1. Load a Template Raster (Logging Richness) to build the Realm Grid
  template_path <- raster_path(paste0(taxa_prefix, "_logging_richness_raster.tif"))
  if (!file.exists(template_path)) {
    cat("  [!] Template raster missing. Skipping taxon.\n"); next
  }
  
  # 2. Build Master Realm Raster (Once per Taxon)
  cat("  -> Building Master Realm Raster (Gap Filling)...\n")
  template_r <- rast(template_path)
  realm_obj <- prepare_taxon_realm_raster(template_r, land_vect, realms_global_sf, realm_col_name)
  master_realm_r <- realm_obj$raster
  master_lookup <- realm_obj$lookup
  
  # --- C. Comparison Loop ---
  for (comp in comparisons) {
    cat(paste0("  -> Comparing: ", comp$x_lab, " vs ", comp$y_lab, "\n"))
    
    # --- D. Measure Loop ---
    for (meas in measures) {
      
      # Construct filenames
      # Format: [taxa]_[threat]_[measure]_raster.tif
      file_x <- raster_path(paste0(taxa_prefix, "_", comp$x, "_", meas$tag, "_raster.tif"))
      file_y <- raster_path(paste0(taxa_prefix, "_", comp$y, "_", meas$tag, "_raster.tif"))
      
      if (file.exists(file_x) && file.exists(file_y)) {
        
        # Load and Mask
        rx <- terra::mask(rast(file_x), land_vect)
        ry <- terra::mask(rast(file_y), land_vect)
        
        # Generate Plot
        fname_suffix <- paste0(taxa_prefix, "_", comp$x, "_vs_", comp$y, "_", meas$tag)
        
        generate_comparison_plot(
          x_raster = rx, 
          y_raster = ry, 
          realm_raster = master_realm_r, 
          lookup = master_lookup,
          x_name = paste0(comp$x_lab, " (", meas$label, ")"),
          y_name = paste0(comp$y_lab, " (", meas$label, ")"),
          measure_label = meas$label,
          title_prefix = display_name,
          filename_suffix = fname_suffix,
          colors = my_colors
        )
        
      } else {
        cat(paste0("     [!] Missing files for ", meas$tag, "\n"))
      }
    } # End Measure Loop
  } # End Comparison Loop
} # End Taxon Loop

cat("\n--- All Comparative Maps Generated Successfully ---\n")

# =============================================================================
# COMBINED PLOT: SCATTER (SQUARE & NO GRID) + HOTSPOT MAPS
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)
library(here)

# --- 1. GLOBAL SETUP: Colors & Shapefile ---

# A. Define Colors
my_colors <- c(
  "Afrotropical"   = "#FF6B6B",
  "Australian"     = "#D4A017",  
  "Madagascan"     = "#FFD54F",  
  "Nearctic"       = "#55A630",
  "Neotropical"    = "#00C853",
  "Oceanina"       = "#00BFA5",  
  "Oriental"       = "#00B0FF",  
  "Palearctic"     = "#2979FF",
  "Panamanian"     = "#6200EA",  
  "Saharo-Arabian" = "#D500F9",  
  "Sino-Japanese"  = "#C51162"   
)

# B. Pre-process the Realm Shapefile
realms_path <- here::here("data", "raw", "realms", "newRealms.shp")

if (file.exists(realms_path)) {
  cat("Loading and Fixing Realm Shapefile...\n")
  suppressMessages(sf::sf_use_s2(FALSE))
  
  realms_global_sf <- st_read(realms_path, quiet = TRUE)
  
  # FIX DATE LINE
  realms_global_sf <- st_make_valid(realms_global_sf)
  realms_global_sf <- st_wrap_dateline(realms_global_sf, options = c("WRAPDATELINE=YES"))
  realms_global_sf <- st_make_valid(realms_global_sf)
  
  # Identify Realm Column
  potential_cols <- names(realms_global_sf)[sapply(realms_global_sf, function(x) {
    len <- length(unique(x))
    return(len > 4 & len < 20) 
  })]
  realm_col_name <- grep("Realm", potential_cols, value = TRUE, ignore.case = TRUE)[1]
  if (is.na(realm_col_name)) realm_col_name <- potential_cols[1]
  
  cat(paste0("  -> Shapefile ready. Using column: '", realm_col_name, "'\n"))
  
} else {
  warning("Realm shapefile not found at: ", realms_path, ". Scatter plot may fail.")
}

# --- 2. CUSTOM FUNCTION: Realm Scatter Plot (Modified) ---

create_realm_scatter_custom <- function(r_x, r_y, display_name, realm_sf) {
  
  cat(paste0("\n>>> GENERATING SCATTER PLOT: ", display_name, " <<<\n"))
  
  # Check environment for mask
  if (!exists("world_map_moll")) stop("'world_map_moll' not found in environment. Please load it from mapping.R")
  
  # C. Mask to Land
  sr_raster_masked <- terra::mask(r_x, vect(world_map_moll))
  threat_raster_masked <- terra::mask(r_y, vect(world_map_moll))
  
  # D. Project & Rasterize Realms
  cat("  -> Rasterizing Realms...\n")
  realms_moll <- st_transform(realm_sf, crs(sr_raster_masked))
  realms_moll <- st_make_valid(realms_moll)
  
  realms_moll$rid_numeric <- as.numeric(as.factor(realms_moll[[realm_col_name]]))
  
  realm_raster <- terra::rasterize(realms_moll, sr_raster_masked, 
                                   field = "rid_numeric", touches = TRUE)
  
  # E. Aggressive Gap Filling
  cat("  -> Gap Filling...\n")
  data_mask <- sr_raster_masked
  data_mask[data_mask <= 0] <- NA 
  
  for(i in 1:30) {
    missing_mask <- !is.na(data_mask) & is.na(realm_raster)
    missing_count <- global(missing_mask, "sum", na.rm=TRUE)[1,1]
    
    if (missing_count == 0) break
    
    win_size <- if(i <= 5) 3 else if(i <= 15) 5 else 7
    filled <- terra::focal(realm_raster, w=win_size, fun="modal", na.policy="only", na.rm=TRUE)
    realm_raster <- terra::cover(realm_raster, filled)
  }
  
  # F. Assemble Dataframe
  stack_data <- c(sr_raster_masked, threat_raster_masked, realm_raster)
  names(stack_data) <- c("richness", "threatened", "realm_id")
  
  plot_data <- as.data.frame(stack_data, na.rm = FALSE)
  
  plot_data <- plot_data %>% 
    filter(!is.na(richness), richness > 0)
  
  # Join Realm Names
  lookup <- unique(st_drop_geometry(realms_moll)[, c("rid_numeric", realm_col_name)])
  colnames(lookup) <- c("realm_id", "Realm_Name")
  
  plot_data <- plot_data %>% 
    left_join(lookup, by = "realm_id") %>% 
    mutate(Realm_Name = ifelse(is.na(Realm_Name), "Unknown", Realm_Name)) %>%
    filter(Realm_Name != "Unknown")
  
  # G. Plotting
  cat("  -> Generating ggplot...\n")
  p <- ggplot(plot_data, aes(x = richness, y = threatened)) +
    geom_point(aes(color = Realm_Name), alpha = 1, size = 1) +
    geom_smooth(method = "gam", color = "black", linewidth = 1.2, se = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = my_colors) +
    labs(
      x = paste0(display_name, " Potential Richness"),
      y = "Threatened Species Count",
      color = "Zoogeographic Realm")+
    theme_minimal() +
    theme(
      # Box border
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      
      # --- NEW: Remove Grid Lines ---
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # --- NEW: Make Square ---
      aspect.ratio = 1
    )
  
  return(p)
}

# --- 3. MAIN EXECUTION ---

cat("Loading Rasters...\n")
path_threat <- raster_path("all_vertebrates_combined_richness_raster.tif")
path_total  <- raster_path("reptile_potential_richness_raster.tif")

r_wood_threat <- rast(path_threat)
r_total_rich  <- rast(path_total)

# 2. Generate Left Panel: Scatter Plot
p_scatter <- create_realm_scatter_custom(
  r_x = r_total_rich,
  r_y = r_wood_threat,
  display_name = "All Vertebrates",
  realm_sf = realms_global_sf
)

# 3. Generate Right Panels: Maps
p_map_continuous <- create_continuous_hotspot_map2(
  input_raster = r_wood_threat,
  map_title = "",
  legend_title = "Richness",
  add_title = TRUE,
  remove_axis_elements = TRUE,
  show_legend = TRUE,
  transform_scale = "identity" 
)

p_map_hotspot <- create_hotspot_map(
  input_raster = r_wood_threat,
  map_title = " ",
  legend_title = "Threat Level",
  color_palette = hotspot_colors,
  add_title = TRUE,
  remove_axis_elements = TRUE,
  show_legend = TRUE
)

# --- 4. COMBINE AND SAVE ---
cat("Combining Plots...\n")

right_col <- plot_grid(
  p_map_continuous, 
  p_map_hotspot, 
  ncol = 1, 
  rel_heights = c(1, 1),
  align = "v"
)

# Adjusted rel_widths to account for the square scatter plot
final_plot <- plot_grid(
  p_scatter, 
  right_col, 
  ncol = 2, 
  rel_widths = c(1, 2), 
  labels = c("A", "B"),
  label_size = 14
)

# Save
outfile <- file.path(summary_dir, "reptile.png")
ggsave(outfile, final_plot, width = 20, height = 10, bg = "white")

cat(paste0("Plot saved to: ", outfile, "\n"))

