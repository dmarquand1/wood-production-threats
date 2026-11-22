# =============================================================================
# SCRIPT: VERTEBRATE THREAT MAPPING - CACHE-BASED BIVARIATE SES HOTSPOTS
#
# Description:
# Natively imports your pre-built sparse cell-by-species matrices.
# Extracts both Logging and Plantation threat vectors directly from your cached 
# metadata, runs parallel null models, exports continuous SES rasters, and 
# maps them onto a 4x4 bivariate scale using biscale (BlueYl).
# =============================================================================

# --- 1. SETUP & LIBRARIES (CLEAN INITIALIZATION) ---
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(rnaturalearth)
library(ggplot2)
library(patchwork)
library(biscale)
library(Matrix)

# --- 2. GLOBAL CONSTANTS & CONFIGURATION ---
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent <- ext(-18040000, 18040000, -9020000, 9020000)
grid_template <- rast(resolution = 110000, crs = moll_crs, extent = global_extent, vals = 0)

# Set execution options
N_SIMS    <- 1000L  # Match your original simulation count
taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

# Define the exact column names used inside your cached metadata dataframes
logging_flag    <- "has_logging_ongoing"
plantation_flag <- "has_plantation_ongoing"

# --- 3. PATH SETTINGS ---
matrix_dir   <- here::here("data", "processed", "hotspots", "null_model", "matrices")
out_dir_rast <- here::here("data", "processed", "hotspots", "ses_rasters")
out_dir_maps <- here::here("outputs", "hotspots", "figures", "ses_bivariate")

dir.create(out_dir_rast, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_maps, recursive = TRUE, showWarnings = FALSE)

# Load global background map
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% 
  st_transform(moll_crs)

# --- 4. NULL MODEL MATRIX MECHANICS ---

logical_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

run_matrix_null_model <- function(sparse_mat, threat_vector, n_simulations = 1000L) {
  cat("  -> Shuffling threat attributes (", n_simulations, " steps via sparse matrix algebra)...\n", sep="")
  
  threat_vec <- logical_flag(threat_vector)
  observed_idx <- which(threat_vec)
  observed_n   <- length(observed_idx)
  
  if (observed_n == 0) stop("The threat vector has no active occurrences in this taxon.")
  
  # Calculate absolute observed richness per active cell
  obs_richness <- as.numeric(Matrix::rowSums(sparse_mat[, observed_idx, drop = FALSE]))
  
  n_cells <- nrow(sparse_mat)
  sum_null   <- numeric(n_cells)
  sumsq_null <- numeric(n_cells)
  
  # Optimized batch processing approach matching your memory framework
  draws <- replicate(
    n_simulations,
    sample.int(ncol(sparse_mat), observed_n, replace = FALSE),
    simplify = FALSE
  )
  
  drawn_species    <- unlist(draws, use.names = FALSE)
  selection_matrix <- Matrix::sparseMatrix(
    i    = drawn_species,
    j    = rep(seq_len(n_simulations), each = observed_n),
    x    = 1,
    dims = c(ncol(sparse_mat), n_simulations)
  )
  
  null_counts <- as.matrix(sparse_mat %*% selection_matrix)
  sum_null    <- rowSums(null_counts)
  sumsq_null  <- rowSums(null_counts * null_counts)
  
  # Calculate mean and standard deviation matrices
  mean_null    <- sum_null / n_simulations
  variance_null <- (sumsq_null - n_simulations * mean_null^2) / (n_simulations - 1)
  variance_null[!is.na(variance_null) & variance_null < 0 & variance_null > -1e-10] <- 0
  sd_null      <- sqrt(pmax(variance_null, 0))
  
  # Compute continuous Standardized Effect Size (SES)
  ses <- (obs_richness - mean_null) / sd_null
  ses[is.na(sd_null) | sd_null == 0] <- 0
  
  return(ses)
}

write_ses_raster_from_cache <- function(ses_values, active_cells, file_path) {
  all_values <- rep(NA_real_, ncell(grid_template))
  all_values[active_cells] <- ses_values
  
  out_rast <- rast(grid_template)
  values(out_rast) <- all_values
  names(out_rast)  <- "ses"
  
  writeRaster(out_rast, file_path, overwrite = TRUE, filetype = "GTiff", gdal = c("COMPRESS=LZW"))
  return(out_rast)
}

# --- 5. BIVARIATE SES VISUALIZATION ENGINE ---

plot_bivariate_ses_map <- function(r_log_ses, r_plant_ses, taxon_name, output_path) {
  cat("  -> Compiling 4x4 Bivariate SES Map via biscale...\n")
  
  s <- c(r_log_ses, r_plant_ses)
  names(s) <- c("Logging_SES", "Plantation_SES")
  s_masked <- terra::mask(s, vect(world_map))
  
  df <- as.data.frame(s_masked, xy = TRUE) %>% 
    filter(!is.na(Logging_SES) & !is.na(Plantation_SES))
  
  if (nrow(df) == 0) {
    cat("  [Warning] No overlapping geographic cells available to map for:", taxon_name, "\n")
    return(NULL)
  }
  
  # Map values to fixed 4-tier risk categories
  df <- df %>% mutate(
    x_cat = case_when(Logging_SES <= 0 ~ 1L, Logging_SES <= 1.96 ~ 2L, Logging_SES <= 3.0 ~ 3L, TRUE ~ 4L),
    y_cat = case_when(Plantation_SES <= 0 ~ 1L, Plantation_SES <= 1.96 ~ 2L, Plantation_SES <= 3.0 ~ 3L, TRUE ~ 4L),
    bi_class = paste0(x_cat, "-", y_cat)
  )
  
  # Render Map
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = bi_class)) +
    biscale::bi_scale_fill(pal = "BlueYl", dim = 4) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() + 
    theme(legend.position = "none")
  
  # Render Bivariate Legend Shield
  leg <- biscale::bi_legend(
    pal = "BlueYl", dim = 4, size = 7,
    xlab = "Logging SES Tiers", ylab = "Plantation SES Tiers"
  )
  
  # Composite output using patchwork
  final_plot <- p_map +
    patchwork::inset_element(
      leg, left = 0.02, bottom = 0.02, right = 0.28, top = 0.28, align_to = "panel"
    )
  
  ggsave(output_path, final_plot, width = 11, height = 7, dpi = 300, bg = "white")
}

# --- 6. HIGH-HYGIENE ISOLATED COMPUTER LOOP ---

for (taxon in taxa_list) {
  cat("\n=====================================================================\n")
  cat("PROCESSING TAXON PROFILE:", toupper(taxon), "\n")
  cat("=====================================================================\n")
  
  path_log_tif   <- file.path(out_dir_rast, paste0(taxon, "_logging_ses.tif"))
  path_plant_tif <- file.path(out_dir_rast, paste0(taxon, "_plantation_ses.tif"))
  path_map_png   <- file.path(out_dir_maps, paste0(taxon, "_ses_bivariate_4x4.png"))
  
  if (file.exists(path_log_tif) && file.exists(path_plant_tif)) {
    cat("  -> Pre-existing SES rasters discovered. Skipping null model processing.\n")
    r_log_ses   <- rast(path_log_tif)
    r_plant_ses <- rast(path_plant_tif)
  } else {
    
    # 1. Target and parse your pre-existing matrix cache file
    cache_path <- file.path(matrix_dir, paste0(taxon, "_cell_species_presence_matrix.rds"))
    if (!file.exists(cache_path)) {
      stop("Cached matrix file not found for ", taxon, " at: ", cache_path, ". Run your matrix generation framework first.")
    }
    
    cat("  -> Loading presence matrix metadata directly from cache layer...\n")
    presence <- readRDS(cache_path)
    
    # Check for metadata flag alignment safety
    if (!logging_flag %in% names(presence$species) || !plantation_flag %in% names(presence$species)) {
      stop("Threat flag names do not align with columns inside cached species metadata dataframe.")
    }
    
    # 2. Run Logging Null Model
    cat("  -> Initiating Logging Matrix Null Model...\n")
    log_scores  <- run_matrix_null_model(presence$matrix, presence$species[[logging_flag]], N_SIMS)
    r_log_ses   <- write_ses_raster_from_cache(log_scores, presence$active_cells, path_log_tif)
    rm(log_scores); tmp <- gc(verbose = FALSE)
    
    # 3. Run Plantation Null Model
    cat("  -> Initiating Plantation Matrix Null Model...\n")
    plant_scores <- run_matrix_null_model(presence$matrix, presence$species[[plantation_flag]], N_SIMS)
    r_plant_ses  <- write_ses_raster_from_cache(plant_scores, presence$active_cells, path_plant_tif)
    rm(plant_scores); tmp <- gc(verbose = FALSE)
    
    # Evict cache block from workspace
    rm(presence); tmp <- gc(verbose = FALSE)
  }
  
  # 4. Generate Bivariate Diagnostic Maps
  plot_bivariate_ses_map(r_log_ses, r_plant_ses, taxon, path_map_png)
  
  # 5. Core Memory Reset Shield
  cat("  -> Flushing allocated vector objects for next taxon cycle...\n")
  suppressWarnings(rm(r_log_ses, r_plant_ses))
  gc(verbose = FALSE)
  
  cat("SUCCESS: Finished processing cycle for:", toupper(taxon), "\n")
}

cat("\n=====================================================================\n")
cat("PROCESS COMPLETE. ALL BIVARIATE DATA SETS STORED NATIVELY.\n")
cat("=====================================================================\n")