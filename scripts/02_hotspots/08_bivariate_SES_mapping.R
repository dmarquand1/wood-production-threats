# =============================================================================
# SCRIPT: VERTEBRATE THREAT MAPPING - CONTINUOUS SES & RICHNESS BIVARIATE MATRICES
#
# Description:
# Natively imports pre-built sparse cell-by-species matrices.
# Generates four continuous maps per taxon plus a standalone universal legend.
# =============================================================================

# --- 1. SETUP & LIBRARIES (CLEAN INITIALIZATION) ---
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(rnaturalearth)
library(ggplot2)
library(patchwork)
library(Matrix)

# --- 2. GLOBAL CONSTANTS & CONFIGURATION ---
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent <- ext(-18040000, 18040000, -9020000, 9020000)
grid_template <- rast(resolution = 110000, crs = moll_crs, extent = global_extent, vals = 0)

N_SIMS    <- 1000L  
taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

logging_flag    <- "has_logging_ongoing"
plantation_flag <- "has_plantation_ongoing"

# --- 3. PATH SETTINGS ---
matrix_dir      <- here::here("data", "processed", "hotspots", "null_model", "matrices")
out_dir_rast    <- here::here("data", "processed", "hotspots", "ses_rasters")
out_dir_figures <- here::here("outputs", "hotspots", "figures", "ses_diagnostics")

dir.create(out_dir_rast, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figures, recursive = TRUE, showWarnings = FALSE)

# Load global background map
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% 
  st_transform(moll_crs)

# Define Primary Paint Blending Anchors (Global Constants for Consistency)
c00 <- "#e8e8e8" # Bottom-Left:  Low X, Low Y (Neutral Gray)
c10 <- "#d9be00" # Bottom-Right: High X, Low Y (Cyan-Blue)
c01 <- "#0088d9" # Top-Left:      Low X, High Y (Vibrant Yellow)
c11 <- "#007b00" # Top-Right:     High X, High Y (Overlapping Mix Green)

mix_2d_colors <- function(x_vec, y_vec) {
  rgb00 <- col2rgb(c00)[,1]; rgb10 <- col2rgb(c10)[,1]
  rgb01 <- col2rgb(c01)[,1]; rgb11 <- col2rgb(c11)[,1]
  
  r_x0 <- (1 - x_vec) * rgb00[1] + x_vec * rgb10[1]
  g_x0 <- (1 - x_vec) * rgb00[2] + x_vec * rgb10[2]
  b_x0 <- (1 - x_vec) * rgb00[3] + x_vec * rgb10[3]
  
  r_x1 <- (1 - x_vec) * rgb01[1] + x_vec * rgb11[1]
  g_x1 <- (1 - x_vec) * rgb01[2] + x_vec * rgb11[2]
  b_x1 <- (1 - x_vec) * rgb01[3] + x_vec * rgb11[3]
  
  r_final <- (1 - y_vec) * r_x0 + y_vec * r_x1
  g_final <- (1 - y_vec) * g_x0 + y_vec * g_x1
  b_final <- (1 - y_vec) * b_x0 + y_vec * b_x1
  
  rgb(r_final, g_final, b_final, maxColorValue = 255)
}

# --- 4. NULL MODEL & RICHNESS EXTRACTION MECHANICS ---

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
  
  obs_richness <- as.numeric(Matrix::rowSums(sparse_mat[, observed_idx, drop = FALSE]))
  
  n_cells <- nrow(sparse_mat)
  sum_null   <- numeric(n_cells)
  sumsq_null <- numeric(n_cells)
  
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
  
  mean_null    <- sum_null / n_simulations
  variance_null <- (sumsq_null - n_simulations * mean_null^2) / (n_simulations - 1)
  variance_null[!is.na(variance_null) & variance_null < 0 & variance_null > -1e-10] <- 0
  sd_null      <- sqrt(pmax(variance_null, 0))
  
  ses <- (obs_richness - mean_null) / sd_null
  ses[is.na(sd_null) | sd_null == 0] <- 0
  
  return(ses)
}

write_raster_from_values <- function(values_vec, active_cells, file_path) {
  all_values <- rep(NA_real_, ncell(grid_template))
  all_values[active_cells] <- values_vec
  
  out_rast <- rast(grid_template)
  values(out_rast) <- all_values
  names(out_rast)  <- "val"
  
  if(!is.null(file_path)) {
    writeRaster(out_rast, file_path, overwrite = TRUE, filetype = "GTiff", gdal = c("COMPRESS=LZW"))
  }
  return(out_rast)
}

# --- 5. VISUALIZATION ENGINES ---

# Engine A: Continuous Symmetrical Divergent Map
plot_continuous_ses_map <- function(ses_raster, title_text, output_path) {
  cat("  -> Generating Continuous Map (No Title):", title_text, "\n")
  
  r_masked <- terra::mask(ses_raster, terra::vect(world_map))
  df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
  names(df) <- c("x", "y", "value")
  
  if (nrow(df) == 0) return(NULL)
  
  max_abs <- max(abs(df$value), na.rm = TRUE)
  if(max_abs == 0) max_abs <- 1 
  
  p <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = value)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-max_abs, max_abs), name = "SES"
    ) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = unit(1.25, "cm")
    )
  
  ggsave(output_path, p, width = 10, height = 5, dpi = 300, bg = "white")
}

# Engine B: True Continuous 2D Map - COMPLETELY UNCLUTTERED
plot_true_continuous_bivariate <- function(r_x_metric, r_y_metric, console_log_text, output_path) {
  cat("  -> Compiling Clean Bivariate Map:", console_log_text, "\n")
  
  s <- c(r_x_metric, r_y_metric)
  names(s) <- c("X", "Y")
  s_masked <- terra::mask(s, vect(world_map))
  df <- as.data.frame(s_masked, xy = TRUE) %>% filter(!is.na(X) & !is.na(Y))
  
  if (nrow(df) == 0) return(NULL)
  
  max_x <- max(df$X, na.rm = TRUE)
  max_y <- max(df$Y, na.rm = TRUE)
  
  if(max_x <= 0) max_x <- 1
  if(max_y <= 0) max_y <- 1
  
  df <- df %>% mutate(
    X_norm = pmax(X, 0) / max_x,
    Y_norm = pmax(Y, 0) / max_y
  )
  
  df$mix_color <- mix_2d_colors(df$X_norm, df$Y_norm)
  
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "#f0f0f0", color = NA) +
    geom_tile(data = df, aes(x = x, y = y, fill = mix_color)) +
    scale_fill_identity() + 
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0, "cm")
    )
  
  ggsave(output_path, p_map, width = 11, height = 5.5, dpi = 300, bg = "white")
}

# Engine C: STANDALONE UNIVERSAL SQUARE LEGEND
plot_standalone_bivariate_legend <- function(output_path) {
  cat("  -> Generating Standalone Universal Square Legend Asset...\n")
  
  # Set up a high-res relative 0 to 1 matrix grid space
  leg_grid <- expand.grid(
    X_norm = seq(0, 1, length.out = 300),
    Y_norm = seq(0, 1, length.out = 300)
  )
  leg_grid$leg_color <- mix_2d_colors(leg_grid$X_norm, leg_grid$Y_norm)
  
  p_legend <- ggplot(leg_grid) +
    geom_tile(aes(x = X_norm, y = Y_norm, fill = leg_color)) +
    scale_fill_identity() +
    # Remove text/tick labels completely
    scale_x_continuous(breaks = NULL, expand = c(0,0)) +
    scale_y_continuous(breaks = NULL, expand = c(0,0)) +
    labs(
      title = "Threatened Richness/SES", 
      x = "Natural Forestry", 
      y = "Plantation"
    ) +
    coord_fixed() + 
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      # Big, highly legible text properties for figure panel combination layouts
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 8)),
      axis.line   = element_line(color = "black", linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15, "pt")
    )
  
  ggsave(output_path, p_legend, width = 4.5, height = 4.5, dpi = 300, bg = "white")
}


# --- 6. CORE PROCESSING CONTROL PIPELINE ---

# Step A: Run through global taxa loops (saving clean map canvases)
for (taxon in taxa_list) {
  cat("\n=====================================================================\n")
  cat("PROCESSING TAXON PROFILE:", toupper(taxon), "\n")
  cat("=====================================================================\n")
  
  path_log_tif   <- file.path(out_dir_rast, paste0(taxon, "_logging_ses.tif"))
  path_plant_tif <- file.path(out_dir_rast, paste0(taxon, "_plantation_ses.tif"))
  
  path_map_log     <- file.path(out_dir_figures, paste0(taxon, "_logging_ses_continuous.png"))
  path_map_plant   <- file.path(out_dir_figures, paste0(taxon, "_plantation_ses_continuous.png"))
  path_map_biv_ses <- file.path(out_dir_figures, paste0(taxon, "_ses_bivariate_continuous_gradient.png"))
  path_map_biv_ric <- file.path(out_dir_figures, paste0(taxon, "_richness_bivariate_continuous_gradient.png"))
  
  cache_path <- file.path(matrix_dir, paste0(taxon, "_cell_species_presence_matrix.rds"))
  if (!file.exists(cache_path)) {
    stop("Cached matrix file not found for ", taxon, " at: ", cache_path)
  }
  
  presence <- readRDS(cache_path)
  
  log_species_idx   <- which(logical_flag(presence$species[[logging_flag]]))
  plant_species_idx <- which(logical_flag(presence$species[[plantation_flag]]))
  
  log_richness_vec   <- as.numeric(Matrix::rowSums(presence$matrix[, log_species_idx, drop = FALSE]))
  plant_richness_vec <- as.numeric(Matrix::rowSums(presence$matrix[, plant_species_idx, drop = FALSE]))
  
  r_log_rich   <- write_raster_from_values(log_richness_vec, presence$active_cells, file_path = NULL)
  r_plant_rich <- write_raster_from_values(plant_richness_vec, presence$active_cells, file_path = NULL)
  
  if (file.exists(path_log_tif) && file.exists(path_plant_tif)) {
    r_log_ses   <- rast(path_log_tif)
    r_plant_ses <- rast(path_plant_tif)
  } else {
    log_scores  <- run_matrix_null_model(presence$matrix, presence$species[[logging_flag]], N_SIMS)
    r_log_ses   <- write_raster_from_values(log_scores, presence$active_cells, path_log_tif)
    rm(log_scores); tmp = gc(verbose = FALSE)
    
    plant_scores <- run_matrix_null_model(presence$matrix, presence$species[[plantation_flag]], N_SIMS)
    r_plant_ses  <- write_raster_from_values(plant_scores, presence$active_cells, path_plant_tif)
    rm(plant_scores); tmp = gc(verbose = FALSE)
  }
  
  rm(presence); tmp = gc(verbose = FALSE)
  taxon_label <- paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)))
  
  # Plot Single Axis Dashboards
  plot_continuous_ses_map(r_log_ses, paste(taxon_label, "Natural Forest Logging"), path_map_log)
  plot_continuous_ses_map(r_plant_ses, paste(taxon_label, "Industrial Plantations"), path_map_plant)
  
  # Plot Bivariate Maps (No inset legends anymore)
  plot_true_continuous_bivariate(
    r_log_ses, r_plant_ses, 
    console_log_text = paste(taxon_label, "Wood Production Hotspots (SES)"),
    output_path = path_map_biv_ses
  )
  
  plot_true_continuous_bivariate(
    r_log_rich, r_plant_rich, 
    console_log_text = paste(taxon_label, "Threatened Species Richness"),
    output_path = path_map_biv_ric
  )
  
  suppressWarnings(rm(r_log_ses, r_plant_ses, r_log_rich, r_plant_rich, log_richness_vec, plant_richness_vec))
  gc(verbose = FALSE)
}

# Step B: Export the standalone universal square legend key asset
path_universal_legend <- file.path(out_dir_figures, "universal_bivariate_square_legend.png")
plot_standalone_bivariate_legend(path_universal_legend)

cat("\n=====================================================================\n")
cat("PROCESS COMPLETE. RAW GLOBAL MAPS AND STANDALONE LEGEND SAVED.\n")
cat("=====================================================================\n")

