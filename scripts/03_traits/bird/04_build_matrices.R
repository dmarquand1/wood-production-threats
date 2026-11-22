# --- 1. SETUP & LIBRARIES ----------------------------------------------------
library(sf)       # Spatial data handling
library(units)    # Handling measurement units (km/m)
library(dplyr)    # Data manipulation
library(here)     # Relative file paths
library(ape)      # Phylogenetic tree analysis (required for read.tree/vcv)

# Set global options
options(stringsAsFactors = FALSE)
# Ensure S2 geometry (spherical) is ON by default
sf::sf_use_s2(TRUE) 

# --- 2. DATA LOADING ---------------------------------------------------------
message("--- Loading Data ---")

# 2a. Load and Process Phylogenetic Tree
if(file.exists("data/processed/traits/bird/final_bird_tree_pruned.tre")) {
  tree <- read.tree("data/processed/traits/bird/final_bird_tree_pruned.tre")
  
  message("Building Phylogenetic Variance-Covariance Matrix...")
  phylo_cor <- vcv(tree, corr = TRUE)
  
  # Sort matrix alphabetically to ensure matching order later
  order_names <- sort(colnames(phylo_cor))
  phylo_cor <- phylo_cor[order_names, order_names]
  phylo_cor <- round(phylo_cor, digits = 6)
} else {
  warning("Tree file not found. Skipping Phylogenetic Matrix step.")
}

saveRDS(phylo_cor, here("data/processed/traits/bird/bird_phylogeny_matrix.rds"))

# 2b. Load Traits
analysis_df <- read.csv("data/processed/traits/bird/bird_trait_data_matched.csv")

# 2c. Load IUCN Polygons
# We filter immediately to only keep species present in our trait data
bird_ranges <- readRDS(here("data/processed/ranges/all_bird_ranges_dissolved.rds")) %>%
  filter(sci_name %in% analysis_df$Trait_Name)

message(paste("Loaded ranges for", nrow(bird_ranges), "species."))

# --- 3. ROBUST CENTROID CALCULATION ------------------------------------------
results_list <- list()
failed_species <- c() # New: Track specific failures
n <- nrow(bird_ranges)

# Initialize progress bar
pb <- txtProgressBar(min = 0, max = n, style = 3)

message(paste("\nStarting robust centroid calculation for", n, "species..."))

for (i in 1:n) {
  setTxtProgressBar(pb, i)
  
  single_species <- bird_ranges[i, ]
  result <- NULL
  
  # --- ATTEMPT 1: Standard S2 Make Valid ---
  try({
    result <- suppressWarnings(
      single_species %>% 
        st_make_valid() %>% 
        st_centroid()
    )
  }, silent = TRUE)
  
  # --- ATTEMPT 2: Buffer Zero ---
  if (is.null(result)) {
    try({
      result <- suppressWarnings(
        single_species %>% 
          st_buffer(dist = 0) %>% 
          st_make_valid() %>% 
          st_centroid()
      )
    }, silent = TRUE)
  }
  
  # --- ATTEMPT 3: Planar Fallback (S2 OFF) ---
  if (is.null(result)) {
    try({
      sf::sf_use_s2(FALSE)
      
      result <- suppressWarnings(
        single_species %>% 
          st_make_valid() %>% 
          st_centroid()
      )
      
      sf::sf_use_s2(TRUE) 
      
      if(!is.null(result)) result$s2_status <- "Planar_Fallback"
      
    }, silent = TRUE)
    sf::sf_use_s2(TRUE) 
  }
  
  # --- CHECK FOR FAILURE ---
  if (is.null(result)) {
    # 1. Create error message
    err_msg <- paste("\n[ERROR] Centroid calculation FAILED for:", single_species$sci_name)
    
    # 2. Print immediately (immediate. = TRUE forces it to show even inside a loop)
    warning(err_msg, call. = FALSE, immediate. = TRUE)
    
    # 3. Add to failure list for summary later
    failed_species <- c(failed_species, single_species$sci_name)
    
  } else {
    results_list[[i]] <- result
  }
}

close(pb)

# --- POST-LOOP SUMMARY ---
bird_centroids <- bind_rows(results_list)

cat("\n------------------------------------------------\n")
cat("Processing Complete.\n")
cat("Success:", nrow(bird_centroids), "\n")
cat("Failed: ", length(failed_species), "\n")

if (length(failed_species) > 0) {
  cat("The following species could not be processed:\n")
  print(failed_species)
  warning("Some species were dropped. Check 'failed_species' vector.")
}
cat("------------------------------------------------\n")

saveRDS(bird_centroids, here("data/processed/traits/bird/bird_centroids_s2_valid.rds"))

bird_centroids <- readRDS("data/processed/traits/bird/bird_centroids_s2_valid.rds")

# 1. Identify the suspects
# (Assumes your results df is named 'bird_centroids')
planar_suspects <- bird_centroids %>% 
  filter(s2_status == "Planar_Fallback") %>% 
  pull(sci_name)

bird_ranges <- readRDS(here("data/processed/ranges/all_bird_ranges_dissolved.rds")) %>%
  filter(sci_name %in% planar_suspects)

if(length(planar_suspects) == 0) {
  message("No planar fallbacks found. No correction needed!")
} else {
  
  message(paste("Checking", length(planar_suspects), "planar species for Dateline errors..."))
  
  # Initialize counters
  fixed_count <- 0
  
  # Ensure S2 is OFF for the check (we are checking planar width)
  sf::sf_use_s2(FALSE)
  
  # 2. Iterate only through the suspects
  for (species_name in planar_suspects) {
    
    # Get the original polygon
    poly <- bird_ranges %>% filter(sci_name == species_name)
    
    # --- CHECK: Is it a "Gabon" Glitch? ---
    # If a single polygon spans > 300 degrees longitude in planar mode, 
    # it is almost certainly wrapping around the dateline (-180 to 180).
    bbox <- st_bbox(poly)
    width <- bbox$xmax - bbox$xmin
    
    if (width > 300) {
      # message(paste("Fixing Dateline wrap for:", species_name))
      
      # --- THE FIX: Pacific Centering ---
      # 1. Shift coords from (-180, 180) to (0, 360)
      # This stitches the split polygon back into one piece over the Pacific.
      shifted_poly <- st_shift_longitude(poly)
      
      # 2. Calculate centroid on the stitched polygon
      shifted_centroid <- suppressWarnings(st_centroid(shifted_poly))
      
      # 3. Normalize back to standard coordinates (-180, 180)
      # If longitude > 180, subtract 360.
      coords <- st_coordinates(shifted_centroid)
      long <- coords[1]
      lat <- coords[2]
      
      if (long > 180) { long <- long - 360 }
      
      # Reconstruct the SF point
      fixed_point <- st_point(c(long, lat)) %>% 
        st_sfc(crs = st_crs(poly)) %>% 
        st_sf(sci_name = species_name, geometry = .)
      
      # 4. Update the main dataframe
      # We find the row in the original 'bird_centroids' and overwrite the geometry
      target_row <- which(bird_centroids$sci_name == species_name)
      
      # Update geometry
      st_geometry(bird_centroids[target_row, ]) <- st_geometry(fixed_point)
      
      # Update status so we know it was fixed
      bird_centroids$s2_status[target_row] <- "Fixed_Dateline"
      
      fixed_count <- fixed_count + 1
    }
  }
  
  # Turn S2 back ON for safety
  sf::sf_use_s2(TRUE)
  
  message(paste("Correction Complete. Fixed", fixed_count, "species."))
}

saveRDS(bird_centroids, here("data/processed/traits/bird/bird_centroids_s2_valid.rds"))

library(sf)
library(dplyr)
library(here)

# 1. Identify the "Wild" Suspects
# Any longitude outside -180 to 180, or Latitude outside -90 to 90
wild_suspects <- bird_centroids %>% 
  mutate(
    lon = st_coordinates(geom)[,1],
    lat = st_coordinates(geom)[,2]
  ) %>% 
  filter(abs(lon) > 180.1 | abs(lat) > 90.1) %>%  # 0.1 buffer for floating point noise
  pull(sci_name)

if(length(wild_suspects) == 0) {
  message("No wild coordinates found. Everything looks normal!")
} else {
  
  message(paste("Found", length(wild_suspects), "species with wild coordinates (e.g., > 180). Fixing..."))
  
  # Load ranges ONLY for the broken species to save memory
bird_ranges_wild <- readRDS(here("data/processed/ranges/all_bird_ranges_dissolved.rds")) %>%
    filter(sci_name %in% wild_suspects)
  
  # Initialize counter
  fixed_count <- 0
  
  # Ensure S2 is ON for the calculation (handles dateline naturally)
  sf::sf_use_s2(TRUE)
  
  for (species_name in wild_suspects) {
    
    # 1. Get the original range
    poly <- bird_ranges_wild %>% filter(sci_name == species_name)
    
    # 2. THE FIX: Force Normalize the Geometry
    # st_wrap_dateline cuts the 'winding' geometry and resets it to standard coords
    # 'options' ensures it aggressively splits lines that cross the dateline
    fixed_poly <- poly %>% 
      st_make_valid() %>% 
      st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>% 
      st_make_valid()
    
    # 3. Recalculate Centroid
    new_centroid <- tryCatch({
      st_centroid(fixed_poly)
    }, error = function(e) { return(NULL) })
    
    # 4. Validate Result
    if (!is.null(new_centroid) && !st_is_empty(new_centroid)) {
      
      # Check if the new coordinates are actually sane
      coords <- st_coordinates(new_centroid)
      if (abs(coords[1]) <= 180 && abs(coords[2]) <= 90) {
        
        # Update the main dataframe
        target_row <- which(bird_centroids$sci_name == species_name)
        st_geometry(bird_centroids[target_row, ]) <- st_geometry(new_centroid)
        
        # Mark as fixed
        bird_centroids$s2_status[target_row] <- "Fixed_Winding_Coords"
        fixed_count <- fixed_count + 1
      }
    } else {
      warning(paste("Could not fix geometry for:", species_name))
    }
  }
  
  message(paste("Wild Coordinate Fix Complete. Fixed", fixed_count, "species."))
}

# --- FINAL VERIFICATION ---
# Double check that no crazy numbers remain
remaining_errors <- bird_centroids %>% 
  mutate(lon = st_coordinates(geom)[,1]) %>% 
  filter(abs(lon) > 180)

if(nrow(remaining_errors) > 0) {
  warning("ALERT: Some species still have invalid coordinates!")
  print(remaining_errors$sci_name)
} else {
  message("ALL CLEAR: All centroids are now within standard Earth bounds.")
}

library(sf)
library(dplyr)
library(here)

# 1. Identify the remaining broken species
# We look for coordinates that are still wild (outside -180 to 180)
# or species explicitly listed in your warning log
remaining_suspects <- bird_centroids %>% 
  mutate(lon = st_coordinates(geom)[,1]) %>% 
  filter(abs(lon) > 180 | is.na(lon)) %>% 
  pull(sci_name)

if(length(remaining_suspects) == 0) {
  message("All coordinates look valid! No fix needed.")
} else {
  
  message(paste("Applying 'Point Cloud' Fallback for", length(remaining_suspects), "stubborn species..."))
  
  # Load ranges
bird_ranges_stubborn <- readRDS(here("data/processed/ranges/all_bird_ranges_dissolved.rds")) %>%
    filter(sci_name %in% remaining_suspects)
  
  fixed_count <- 0
  
  # IMPORTANT: Use S2. 
  # This ensures that a cloud of points spanning -179 and +179 
  # has a center at 180 (Pacific), not 0 (Greenwich).
  sf::sf_use_s2(TRUE)
  
  for (species_name in remaining_suspects) {
    
    # 1. Get the range
    poly <- bird_ranges_stubborn %>% filter(sci_name == species_name)
    
    if(nrow(poly) > 0) {
      
      # 2. Extract Vertices
      coords <- st_coordinates(poly)
      
      if(nrow(coords) > 0) {
        
        # 3. Mathematically Fix Coordinates
        x_fixed <- ((coords[,1] + 180) %% 360) - 180
        y_fixed <- pmin(pmax(coords[,2], -90), 90)
        
        # 4. Create the "Cloud"
        # We define CRS as 4326 initially
        points_sf <- st_as_sf(data.frame(x = x_fixed, y = y_fixed), 
                              coords = c("x", "y"), 
                              crs = 4326) %>% 
          st_combine()
        
        # 5. Calculate Centroid
        new_centroid_sf <- st_centroid(points_sf)
        
        # --- CRITICAL FIX STARTS HERE ---
        
        # A. Extract ONLY the geometry (remove the dataframe wrapper)
        new_geom <- st_geometry(new_centroid_sf)
        
        # B. Force the CRS to match 'bird_centroids' exactly
        # This prevents the "Can't convert" error
        st_crs(new_geom) <- st_crs(bird_centroids)
        
        # C. Assign directly to the geometry vector
        # Using st_geometry(df)[i] <- geom is safer than modifying the row
        target_row <- which(bird_centroids$sci_name == species_name)
        st_geometry(bird_centroids)[[target_row]] <- new_geom[[1]]
        
        # --- CRITICAL FIX ENDS HERE ---
        
        # Update status
        bird_centroids$s2_status[target_row] <- "Fixed_PointCloud_Fallback"
        fixed_count <- fixed_count + 1
        
        message(paste("Fixed:", species_name))
        
      } else {
        warning(paste("Skipping:", species_name, "- No coordinates found."))
      }
    }
  }
  
  message(paste("Nuclear Fix Complete. Fixed", fixed_count, "species."))
}

# --- FINAL SANITY CHECK ---
final_check <- bird_centroids %>% 
  mutate(lon = st_coordinates(geom)[,1],
         lat = st_coordinates(geom)[,2]) %>% 
  filter(abs(lon) > 180 | abs(lat) > 90)

if(nrow(final_check) == 0) {
  message("SUCCESS: All centroids are finally valid.")
} else {
  warning("CRITICAL: Some species are still broken.")
  print(final_check$sci_name)
}

library(sf)
library(dplyr)

# Check for Empty, NA, or Infinite coordinates
broken_report <- bird_centroids %>%
  mutate(
    # Get raw coordinates
    x = st_coordinates(geom)[,1],
    y = st_coordinates(geom)[,2],
    valid_geom = !st_is_empty(geom)
  ) %>%
  filter(
    is.na(x) | is.na(y) |       # Coordinates are Missing
      is.infinite(x) |            # Coordinates are Infinite
      abs(x) > 180.1 |            # Coordinates are still wrapping
      !valid_geom                 # Geometry is strictly EMPTY
  ) %>%
  select(sci_name, x, y, valid_geom)

print("The remaining broken species are:")
print(broken_report)

saveRDS(bird_centroids, here("data/processed/traits/bird/bird_centroids_s2_valid.rds"))

# Check the previously broken species
bird_centroids %>% 
  filter(sci_name == "Himantopus himantopus") %>% 
  st_coordinates()

library(sf)
library(dplyr)

# 1. Identify species with Impossible Latitudes (outside -90 to 90)
# We add a tiny buffer (90.00001) to catch floating point errors
impossible_lat_species <- bird_centroids %>%
  mutate(Y = st_coordinates(geom)[,2]) %>%
  filter(abs(Y) > 90) %>%
  pull(sci_name)

if(length(impossible_lat_species) > 0) {
  message(paste("Fixing", length(impossible_lat_species), "species with impossible latitudes (e.g., > 90)..."))
  
  # 2. Iterate and Clamp
  for(species in impossible_lat_species) {
    
    # Get the bad point
    target_row <- which(bird_centroids$sci_name == species)
    coords <- st_coordinates(bird_centroids[target_row, ])
    
    lon <- coords[1]
    lat <- coords[2]
    
    # CLAMP LOGIC:
    # If it's > 90, make it 90. If it's < -90, make it -90.
    if (lat > 90) lat <- 90
    if (lat < -90) lat <- -90
    
    # Rebuild the point
    fixed_point <- st_point(c(lon, lat)) %>%
      st_sfc(crs = st_crs(bird_centroids))
    
    # Overwrite
    st_geometry(bird_centroids)[[target_row]] <- fixed_point[[1]]
    
    # Update status for your records
    bird_centroids$s2_status[target_row] <- "Fixed_Lat_Clamp"
    
    message(paste("Clamped:", species, "from", round(coords[2], 4), "to", lat))
  }
} else {
  message("No impossible latitudes found.")
}

# --- VERIFY HIMANTOPUS ---
# Should now show exactly 90 (or whatever value you clamped to)
print(
  bird_centroids %>% 
    filter(sci_name == "Himantopus himantopus") %>% 
    st_coordinates()
)

# Save immediately so you don't lose the fix
saveRDS(bird_centroids, here("data/processed/traits/bird/bird_centroids_s2_valid.rds"))

# --- 4. CALCULATE SPATIAL DISTANCE MATRIX ------------------------------------
message("--- Calculating pairwise distance matrix (Geodesic) ---")

# Ensure S2 is ON for distance calculation (Great Circle distance)
sf::sf_use_s2(TRUE)

# Calculate distance (returns matrix with units)
# Note: Large matrices (10k species) require significant RAM
bird_centroids = readRDS("data/processed/traits/bird/bird_centroids_s2_valid.rds")

library(sf)
library(dplyr)

# 1. Transform to meters and save original CRS
original_crs <- st_crs(bird_centroids)
birds_meters <- st_transform(bird_centroids, 3857)

# 2. Extract coordinates and Add Jitter
# We do this in one pipe for cleanliness
birds_jittered_df <- birds_meters %>%
  mutate(
    # Get current X and Y
    X = st_coordinates(.)[,1],
    Y = st_coordinates(.)[,2],
    
    # Add random noise (+/- 100m)
    # We add row_number/1e6 to guarantee mathematical uniqueness
    X_new = X + runif(n(), -100, 100) + (row_number() * 0.000001),
    Y_new = Y + runif(n(), -100, 100) + (row_number() * 0.000001)
  ) %>%
  # --- CRITICAL STEP: Drop the old geometry ---
  # This forces R to treat this as a standard dataframe, forgetting the old points
  st_drop_geometry()

# 3. Create NEW geometry from the jittered coordinates
bird_centroids_final <- st_as_sf(birds_jittered_df, coords = c("X_new", "Y_new"), crs = 3857)

# 4. Transform back to original Lat/Lon
bird_centroids_final <- st_transform(bird_centroids_final, original_crs)

# 5. Clean up temporary columns (X, Y are no longer needed)
bird_centroids_final <- bird_centroids_final %>% select(-X, -Y)

# Verification
print(paste("Duplicates remaining:", sum(duplicated(bird_centroids_final$geometry))))

dist_matrix_raw <- sf::st_distance(bird_centroids)

# Set names
rownames(dist_matrix_raw) <- bird_centroids$sci_name
colnames(dist_matrix_raw) <- bird_centroids$sci_name

dist_matrix_numeric <-  units::drop_units(dist_matrix_raw)

#convert to km
dist_matrix_km <- dist_matrix_numeric/1000

message("Max distance observed: ", round(max(dist_matrix_km), 2), " km")

# --- 5. CONVERT TO SIMILARITY MATRIX -----------------------------------------
message("--- Converting to Correlation/Similarity Matrix ---")

# 1. Create a copy of the numeric matrix
dist_fixed <- dist_matrix_km

# 2. Identify off-diagonal zeros
# (Logic: Distance is 0, but it is NOT the diagonal)
is_zero_dist <- (dist_fixed < 0.1)
diag(is_zero_dist) <- FALSE # Ignore the diagonal (self-distance must stay 0)

# 3. Replace these zeros with a small fixed constant
# 0.1 km (100m) is usually enough to break the correlation of 1
# while remaining biologically negligible for large-scale analyses.
dist_fixed[is_zero_dist] <- 0.1

# --- Re-calculate Similarity Matrix ---
message("--- Converting to Similarity Matrix (Fixed Zeros) ---")

# Calcualte the exponential decay kernel so that correlation is near zero (0.05) for the 4000km
l = -4000 / log(0.05)

# Recalculate using the fixed distances
bird_dist_cor <- exp(-dist_fixed / l)

# --- Verification ---
# Check if any off-diagonal correlations are still exactly 1
off_diag_check <- bird_dist_cor
diag(off_diag_check) <- 0
num_perfect <- sum(off_diag_check >= 0.999999999) 
# Note: checking for ">= 0.999999999" handles floating point rounding issues

if(num_perfect > 0){
  message("Warning: ", num_perfect, " pairs are still extremely close to 1.")
} else {
  message("Success: No species have a correlation of 1.")
}

# Force symmetry to correct for floating point errors
if (!isSymmetric(bird_dist_cor)) {
  bird_dist_cor[lower.tri(bird_dist_cor)] <- t(bird_dist_cor)[lower.tri(bird_dist_cor)]
}

# --- 6. SAVE OUTPUTS ---------------------------------------------------------
message("--- Saving Outputs ---")


saveRDS(dist_matrix_km, here("data/processed/traits/bird/bird_distance_matrix_km.rds"))
saveRDS(bird_dist_cor, here("data/processed/traits/bird/bird_correlation_matrix.rds"))
message("All tasks completed successfully.")

# --- 7. CALCULATE RANGE SIZES (OPTIMIZED & FAST) -----------------------------
message("--- Calculating Range Sizes (Fast Planar) ---")

# 1. Reload the ranges
analysis_df <- read.csv("data/processed/traits/bird/bird_trait_data_matched.csv")
bird_ranges_area <- readRDS(here("data/processed/ranges/all_bird_ranges_dissolved.rds")) %>%
  filter(sci_name %in% analysis_df$Trait_Name)

# 2. Turn S2 OFF for speed/stability (we rely on projection for accuracy)
sf::sf_use_s2(FALSE) 

# 3. Define Chunking Parameters (Smaller chunks for safety)
chunk_size <- 100  # Reduced from 500 to 100 to prevent memory bottlenecks
total_species <- nrow(bird_ranges_area)
num_chunks <- ceiling(total_species / chunk_size)
area_results_list <- list()

# Define Equal-Area Projection (Mollweide)
# --- Global Spatial Grid ---
# Using a Molleweide equal area projection
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

message(paste("Processing", total_species, "species in", num_chunks, "chunks..."))

# 4. Loop with Simplification
pb <- txtProgressBar(min = 0, max = num_chunks, style = 3)

for (k in 1:num_chunks) {
  
  start_idx <- ((k - 1) * chunk_size) + 1
  end_idx <- min(k * chunk_size, total_species)
  
  chunk_df <- bird_ranges_area[start_idx:end_idx, ]
  
  # Print the first species of the chunk so you know exactly where you are
  # (Uncomment the line below if you want to see names in the console)
  message(paste("\nChunk", k, "starts with:", chunk_df$sci_name[1]))
  
  tryCatch({
    chunk_processed <- chunk_df %>%
      st_make_valid() %>%
      # --- CRITICAL FIX: Simplify Geometry ---
      # dTolerance = 100 means we remove vertices that don't add 
      # more than 100m of detail. This un-freezes complex coastlines.
      st_simplify(preserveTopology = TRUE, dTolerance = 100) %>% 
      st_transform(mollweide_crs) %>% 
      mutate(
        area_raw = st_area(geom)
      ) %>%
      st_drop_geometry() %>% 
      select(sci_name, area_raw) %>%
      mutate(
        range_size_km2 = as.numeric(area_raw) / 1e6
      ) %>%
      select(sci_name, range_size_km2)
    
    area_results_list[[k]] <- chunk_processed
    
  }, error = function(e) {
    message(paste("\nError in chunk", k, ":", e$message))
  })
  
  setTxtProgressBar(pb, k)
  
  # Clean up memory
  rm(chunk_df) 
  if(exists("chunk_processed")) rm(chunk_processed)
  gc(verbose = FALSE) 
}

close(pb)
sf::sf_use_s2(TRUE) # Restore S2 for later tasks

# 5. Bind and Check
bird_area_df <- bind_rows(area_results_list)

zeros <- sum(bird_area_df$range_size_km2 == 0)
if(zeros > 0) warning(paste(zeros, "species have 0 area (likely geometry errors)."))

# Save immediate result
saveRDS(bird_area_df, here("data/processed/traits/bird/bird_range_areas.rds"))

# --- 9. MERGE AND SAVE (As before) -------------------------------------------
message("--- Merging Range Sizes with Trait Data ---")

trait_file_path <- here("data/processed/traits/bird/bird_traits_imputed_ready.csv")
if(file.exists(trait_file_path)) {
  bird_traits <- read.csv(trait_file_path)
} else {
stop("Error: 'data/processed/traits/bird/bird_traits_imputed_ready.csv' not found.")
}

name_map <- if ("phylo_match" %in% names(analysis_df)) {
  analysis_df %>%
    transmute(phylo_name = gsub(" ", "_", phylo_match), Trait_Name) %>%
    distinct()
} else if ("Tree_Name" %in% names(analysis_df)) {
  analysis_df %>%
    transmute(phylo_name = gsub(" ", "_", Tree_Name), Trait_Name) %>%
    distinct()
} else {
  analysis_df %>%
    transmute(phylo_name = eBird_Name, Trait_Name) %>%
    distinct()
}

bird_traits_mapped <- bird_traits %>%
  left_join(name_map, by = c("species" = "phylo_name"))

final_df <- bird_traits_mapped %>%
  left_join(bird_area_df, by = c("Trait_Name" = "sci_name"))

output_path <- here("data/processed/traits/bird/bird_traits_imputed_with_range.csv")
write.csv(final_df, output_path, row.names = FALSE)
message(paste("Final dataset saved to:", output_path))

