rm(list = ls())
gc()

# =============================================================================
# CUBE-ROOT WEIGHTED NULL MODEL: SES FOR TRANSFORMED THREAT RASTERS
#
# Description:
# Uses the weighted cell-by-species matrices built in
# 12_build_weighted_transformed_matrices.R (entries = cube-root inverse-range
# weighted per-species-per-cell contribution) to compute standardized effect
# sizes (SES) for the combined wood-production threat, logging threat, and
# plantation threat, under the same species-label-reshuffling null model used
# for the richness-based SES rasters (06_null_model_hotspots.R) -- applied to
# continuous weights instead of binary presence/absence.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c("Matrix", "dplyr", "here", "terra")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

# --- 1. Configuration --------------------------------------------------------
taxa_to_run  <- c("bird", "mammal", "amphibian", "reptile", "fish")
threat_flags <- c("has_wood_prod_ongoing", "has_logging_ongoing", "has_plantation_ongoing")
n_random  <- 1000L
null_seed <- 20260511L
rebuild   <- FALSE

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000,
  crs        = moll_crs,
  extent     = terra::ext(-18040000, 18040000, -9020000, 9020000),
  vals       = 0
)

matrix_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "matrices")
raster_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "rasters")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

logical_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

# --- 2. Weighted null model ----------------------------------------------------
run_weighted_null_model <- function(sparse_mat, threat_vector, n_simulations, seed) {
  threat_vec   <- logical_flag(threat_vector)
  observed_idx <- which(threat_vec)
  observed_n   <- length(observed_idx)
  if (observed_n == 0) stop("Threat vector has no active occurrences.")
  if (observed_n > ncol(sparse_mat)) stop("Observed species count exceeds species pool.")

  obs_weighted_sum <- as.numeric(Matrix::rowSums(sparse_mat[, observed_idx, drop = FALSE]))

  set.seed(seed)
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

  null_sums  <- as.matrix(sparse_mat %*% selection_matrix)
  sum_null   <- rowSums(null_sums)
  sumsq_null <- rowSums(null_sums^2)
  rm(null_sums, selection_matrix, draws, drawn_species)
  gc(verbose = FALSE)

  mean_null     <- sum_null / n_simulations
  variance_null <- (sumsq_null - n_simulations * mean_null^2) / (n_simulations - 1)
  variance_null[!is.na(variance_null) & variance_null < 0 & variance_null > -1e-10] <- 0
  sd_null <- sqrt(pmax(variance_null, 0))

  ses <- (obs_weighted_sum - mean_null) / sd_null
  ses[is.na(sd_null) | sd_null == 0] <- 0

  list(ses = ses, observed_n = observed_n, species_pool_n = ncol(sparse_mat))
}

write_raster_from_values <- function(values_vec, active_cells, grid, path) {
  all_values <- rep(NA_real_, terra::ncell(grid))
  all_values[active_cells] <- values_vec
  r <- terra::setValues(terra::rast(grid), all_values)
  names(r) <- "ses"
  terra::writeRaster(r, path, overwrite = TRUE, filetype = "GTiff", gdal = c("COMPRESS=LZW"))
  r
}

# --- 3. Run --------------------------------------------------------------------
for (taxon in taxa_to_run) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("SES (cube-root weighted): ", toupper(taxon), "\n", sep = "")
  cat(strrep("=", 60), "\n\n", sep = "")

  matrix_path <- file.path(matrix_dir, paste0(taxon, "_cell_species_weighted_matrix.rds"))
  if (!file.exists(matrix_path)) {
    stop("Missing weighted matrix for ", taxon, ": ", matrix_path,
         ". Run 12_build_weighted_transformed_matrices.R first.")
  }

  presence <- readRDS(matrix_path)

  for (flag in threat_flags) {
    out_path <- file.path(raster_dir, paste0(taxon, "_", flag, "_transformed_ses_raster.tif"))
    if (file.exists(out_path) && !rebuild) {
      cat("  Already computed, skipping: ", basename(out_path), "\n", sep = "")
      next
    }
    if (!flag %in% names(presence$species)) {
      stop("Flag not found in species metadata for ", taxon, ": ", flag)
    }

    cat("  Running null model for flag: ", flag, "\n", sep = "")
    result <- run_weighted_null_model(presence$matrix, presence$species[[flag]], n_random, null_seed)
    write_raster_from_values(result$ses, presence$active_cells, grid_template, out_path)
    cat("    Observed species: ", result$observed_n, " / ", result$species_pool_n, "\n", sep = "")
    cat("    Saved: ", out_path, "\n", sep = "")
    rm(result)
    gc(verbose = FALSE)
  }

  rm(presence)
  gc(verbose = FALSE)
}

cat("\nAll cube-root weighted SES rasters computed.\n")
