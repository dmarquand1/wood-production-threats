rm(list=ls())
gc()
# =============================================================================
# SCRIPT: CALCULATE fish SPECIES RANGE SIZES
#
# Description:
# This script calculates the total range size (in km2) for fish species.
# It uses an Equal Area projection to ensure accurate area calculation,
# based on the methodology of the inspiration script.
# =============================================================================

# --- 1. SETUP: Load necessary libraries ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,          # For data manipulation
  sf,             # For working with spatial vector data
  units,           # For handling measurement units
  here
)

# --- 2. DEFINE PROJECTION ---
# We use the Mollweide equal area projection
# This is crucial: Standard Lat/Lon (WGS84) distorts area, making high-latitude
# ranges appear larger than they are.
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# --- 3. LOAD DATA ---
analysis_df <- read.csv("data/processed/traits/fish/fish_trait_data_matched.csv")
fish_ranges <- readRDS(here("data/processed/ranges/all_fish_ranges_dissolved.rds")) %>% 
  filter(sci_name %in% analysis_df$IUCN_Name)

# Get list of unique species
species_list <- unique(fish_ranges$sci_name)
total_species <- length(species_list)

# Initialize an empty list to store results (more efficient than growing a dataframe)
results_list <- list()

# --- 4. THE LOOP ---
cat("--- Starting Range Size Calculation for", total_species, "Species ---\n")

start_time <- Sys.time()

for (i in 1:total_species) {
  
  # A. Get current species name
  current_species_name <- species_list[i]
  
  # B. Filter the spatial data for just this species
  single_species_sf <- fish_ranges %>% 
    filter(sci_name == current_species_name)
  
  # C. Transform to Equal Area Projection
  species_transformed <- st_transform(single_species_sf, crs = moll_crs)
  
  # D. Calculate Area
  # st_area returns area in square meters (because proj is in meters)
  # We sum it in case the range is split into multiple polygons
  area_m2 <- sum(st_area(species_transformed))
  
  # E. Convert to square kilometers for readability
  area_km2 <- set_units(area_m2, km^2)
  
  # F. Store in the list
  results_list[[i]] <- data.frame(
    species_name = current_species_name,
    range_size_km2 = as.numeric(area_km2)
  )
  
  # G. Progress Update
  if (i %% 10 == 0 || i == total_species) { # Update every 10 species or at the end
    elapsed_time <- difftime(Sys.time(), start_time, units = "mins")
    avg_time <- as.numeric(elapsed_time) / i
    est_remain <- round(avg_time * (total_species - i), 1)
    
    cat(paste0(" -> Processed ", i, " of ", total_species, 
               " | ", current_species_name, 
               " | Est. Remaining: ", est_remain, " mins\n"))
  }
}

# --- 5. FINALIZE OUTPUT ---
# Combine the list into a single clean Dataframe
final_range_sizes <- do.call(rbind, results_list)
final_range_sizes$log_range_size = log10(final_range_sizes$range_size_km2) %>% 
  scale()
hist(final_range_sizes$log_range_size)

# View the top of the dataframe
print(head(final_range_sizes))


cat("\nCalculation Complete.\n")

traits = read.csv("data/processed/traits/fish/fish_trait_data_matched.csv")


library(dplyr)

traits <- traits %>%
  left_join(final_range_sizes %>% select(species_name, log_range_size), 
            by = c("IUCN_Name" = "species_name"))

# Add range size to the imputed trait data
model_data = read.csv("data/processed/traits/fish/fish_traits_imputed_final.csv") %>% 
  left_join(traits %>% select(Final_Tree_Tip,log_range_size,Threat_Logging,Threat_Plantation),
            by = c("species" = "Final_Tree_Tip"))


# Make sure the trait data has the trait and phylogeny names
model_data <- model_data %>% 
  left_join(analysis_df %>% select(Final_Tree_Tip, IUCN_Name),
            by = c("species" = "Final_Tree_Tip")) %>% 
  rename(phylo_name = species) %>% 
  rename(spatial_name = IUCN_Name)

write.csv(model_data, file="data/processed/traits/fish/fish_trait_model_final_data.csv")
df=read.csv("data/processed/traits/fish/fish_trait_model_final_data.csv")
df2 = df %>% left_join(analysis_df %>% select(Species, habitat_breadth), 
                       by = c("spatial_name" = "Species"))
