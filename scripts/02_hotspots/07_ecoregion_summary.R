# ==============================================================================
# Script: Extract Average Threatened Species Richness per Ecoregion
# ==============================================================================

# Load required libraries
library(sf)
library(terra)
library(dplyr)

# 1. Define File Paths
shapefile_path <- ".../redListAnalysis/archive/reference_material/geodata/Ecoregions2017.shp"
raster_folder  <- ".../redListAnalysis/data/processed/hotspots/rasters"
output_dir     <- ".../redListAnalysis/outputs/hotspots/tables"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# Automatically find all 5 taxon rasters matching your file naming convention
raster_files <- list.files(
  path = raster_folder, 
  pattern = "_combined_richness_raster\\.tif$", 
  full.names = TRUE
)

# Verify that all 5 rasters were found
if (length(raster_files) != 5) {
  warning(paste("Expected 5 rasters, but found", length(raster_files), ". Check your folder!"))
}

# Automatically generate clean column names from the filenames
taxon_names <- paste0("mean_", gsub("_combined_richness_raster\\.tif", "", basename(raster_files)))

# ------------------------------------------------------------------------------
# 2. Load Spatial Data
# ------------------------------------------------------------------------------
cat("Loading ecoregions shapefile...\n")
ecoregions <- st_read(shapefile_path)

cat("Loading and stacking taxon rasters...\n")
raster_stack <- rast(raster_files)
names(raster_stack) <- taxon_names

# ------------------------------------------------------------------------------
# 3. Align Coordinate Reference Systems (CRS)
# ------------------------------------------------------------------------------
if (st_crs(ecoregions) != crs(raster_stack)) {
  cat("CRS mismatch detected. Reprojecting ecoregions to match rasters...\n")
  ecoregions <- st_transform(ecoregions, crs(raster_stack))
}

# ------------------------------------------------------------------------------
# 4. Extract Zonal Statistics
# ------------------------------------------------------------------------------
cat("Extracting average richness per ecoregion (110km cell resolution)...\n")

extracted_means <- terra::extract(
  raster_stack, 
  ecoregions, 
  fun = mean, 
  na.rm = TRUE, 
  touches = TRUE # Ensures small ecoregions still grab data from overlapping 110km cells
)

# ------------------------------------------------------------------------------
# 5. Combine and Clean Data
# ------------------------------------------------------------------------------
cat("Processing final table...\n")

final_data <- ecoregions %>%
  st_drop_geometry() %>%                  # Remove spatial geometry for CSV export
  mutate(ID = row_number()) %>%          # Create temporary ID matching terra's output
  left_join(extracted_means, by = "ID") %>%
  select(-ID)                             # Drop the temporary ID column

# ------------------------------------------------------------------------------
# 6. Export to CSV
# ------------------------------------------------------------------------------
output_csv <- file.path(output_dir, "ecoregions_threatened_richness_summary.csv")
write.csv(final_data, output_csv, row.names = FALSE)

cat(paste0("Success! Data exported to:\n", output_csv, "\n"))
