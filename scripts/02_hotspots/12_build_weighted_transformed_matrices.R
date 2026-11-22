rm(list = ls())
gc()

# =============================================================================
# BUILD CUBE-ROOT RANGE-WEIGHTED CELL-BY-SPECIES MATRICES
#
# Description:
# Companion to the binary presence/absence null-model matrices
# (data/processed/hotspots/null_model/matrices/*_cell_species_presence_matrix.rds).
# For each taxon, builds a sparse cell x species matrix whose entries are the
# cube-root inverse-range-weighted "Transformed Threat" contribution of that
# species to that cell: (area_in_cell / total_range_area)^(1/3), matching the
# per-species formula used in calculate_or_load_transformed_raster() in
# 01_build_richness_rasters.R. Column sums over any species subset reproduce
# the corresponding *_transformed_threat_raster.tif, which lets the same
# reshuffling null-model machinery used for richness SES be reused to build a
# cube-root-weighted SES (see 13_weighted_null_model_ses.R).
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c("Matrix", "dplyr", "sf", "terra", "here")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

# --- 1. Configuration --------------------------------------------------------
taxa_to_run <- c("bird", "mammal", "amphibian", "reptile", "fish")
rebuild <- FALSE

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000,
  crs        = moll_crs,
  extent     = terra::ext(-18040000, 18040000, -9020000, 9020000),
  vals       = 0
)
cell_area <- terra::res(grid_template)[1] * terra::res(grid_template)[2]

range_files <- c(
  amphibian = here::here("data", "processed", "ranges", "all_amphibian_ranges_dissolved.rds"),
  bird      = here::here("data", "processed", "ranges", "all_bird_ranges_dissolved.rds"),
  fish      = here::here("data", "processed", "ranges", "all_fish_ranges_dissolved.rds"),
  mammal    = here::here("data", "processed", "ranges", "all_mammal_ranges_dissolved.rds"),
  reptile   = here::here("data", "processed", "ranges", "all_reptile_ranges_dissolved.rds")
)

matrix_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "matrices")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)

logical_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

load_assessment_species <- function(taxon_key) {
  path <- here::here("data", "processed", "assessments", paste0(taxon_key, "_assessment_species_summary.rds"))
  if (!file.exists(path)) stop("Missing assessment species summary: ", path)
  dat <- readRDS(path)
  dat <- dat[!dat$category %in% c("EX", "EW"), , drop = FALSE]
  dat <- dat[!is.na(dat$scientific_name) & nzchar(dat$scientific_name), , drop = FALSE]
  dat[!duplicated(dat$scientific_name), , drop = FALSE]
}

# --- 2. Core builder -----------------------------------------------------------
build_weighted_matrix <- function(taxon_key) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Building cube-root weighted matrix: ", toupper(taxon_key), "\n", sep = "")
  cat(strrep("=", 60), "\n\n", sep = "")

  assessment <- load_assessment_species(taxon_key)
  ranges <- readRDS(range_files[[taxon_key]])
  if (!inherits(ranges, "sf")) stop("Range object is not sf: ", range_files[[taxon_key]])
  if (!"sci_name" %in% names(ranges)) stop("Range object missing sci_name: ", range_files[[taxon_key]])

  ranges <- ranges[!is.na(ranges$sci_name) & ranges$sci_name %in% assessment$scientific_name, , drop = FALSE]
  ranges <- ranges[!duplicated(ranges$sci_name), , drop = FALSE]
  if (nrow(ranges) == 0) stop("No extant assessed species with ranges for ", taxon_key)

  species_meta <- assessment[match(ranges$sci_name, assessment$scientific_name), , drop = FALSE]
  species_meta$species_id <- ranges$sci_name

  cat("Reprojecting ", nrow(ranges), " ranges to Mollweide...\n", sep = "")
  ranges_moll <- sf::st_transform(ranges, moll_crs)

  n_species <- nrow(ranges_moll)
  cell_list   <- vector("list", n_species)
  weight_list <- vector("list", n_species)
  skipped_zero_area <- 0L
  start_time <- Sys.time()

  for (i in seq_len(n_species)) {
    sp_geom <- ranges_moll[i, , drop = FALSE]
    total_area <- as.numeric(sf::st_area(sp_geom))
    if (!is.finite(total_area) || total_area <= 0) {
      skipped_zero_area <- skipped_zero_area + 1L
      next
    }

    cov <- try(terra::rasterize(terra::vect(sp_geom), grid_template, cover = TRUE), silent = TRUE)
    if (inherits(cov, "try-error")) {
      sp_geom_fix <- sf::st_make_valid(sp_geom)
      cov <- terra::rasterize(terra::vect(sp_geom_fix), grid_template, cover = TRUE)
    }

    prop <- (cov * cell_area) / total_area
    vals <- terra::values(prop, mat = FALSE)
    nz <- which(!is.na(vals) & vals > 0)
    if (length(nz) == 0) next

    vals_nz <- vals[nz]
    vals_nz[vals_nz > 1] <- 1

    cell_list[[i]]   <- nz
    weight_list[[i]] <- vals_nz^(1 / 3)

    if (i %% 250 == 0 || i == n_species) {
      elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
      avg <- elapsed / i
      remaining <- round(avg * (n_species - i), 1)
      cat("  ", i, " / ", n_species, " species (", elapsed, " min elapsed, ~",
          remaining, " min remaining)\n", sep = "")
    }
    if (i %% 500 == 0) gc(verbose = FALSE)
  }

  has_data <- !vapply(cell_list, is.null, logical(1))
  cell_list   <- cell_list[has_data]
  weight_list <- weight_list[has_data]
  species_meta_kept <- species_meta[has_data, , drop = FALSE]

  if (length(cell_list) == 0) stop("No species produced nonzero coverage for ", taxon_key)

  active_cells     <- sort(unique(unlist(cell_list, use.names = FALSE)))
  nonzero_lengths  <- lengths(cell_list)
  cell_ids         <- unlist(cell_list, use.names = FALSE)
  row_ids          <- match(cell_ids, active_cells)
  col_ids          <- rep(seq_along(cell_list), nonzero_lengths)
  x_vals           <- unlist(weight_list, use.names = FALSE)

  weighted_matrix <- Matrix::sparseMatrix(
    i = row_ids, j = col_ids, x = x_vals,
    dims = c(length(active_cells), length(cell_list)),
    dimnames = list(as.character(active_cells), species_meta_kept$species_id)
  )

  out <- list(
    matrix       = weighted_matrix,
    active_cells = active_cells,
    species      = species_meta_kept,
    taxa         = taxon_key,
    grid_spec = list(
      nrows = terra::nrow(grid_template), ncols = terra::ncol(grid_template),
      xmin  = terra::xmin(grid_template), xmax  = terra::xmax(grid_template),
      ymin  = terra::ymin(grid_template), ymax  = terra::ymax(grid_template),
      crs   = terra::crs(grid_template)
    ),
    created_at         = Sys.time(),
    skipped_zero_area  = skipped_zero_area
  )

  cache_path <- file.path(matrix_dir, paste0(taxon_key, "_cell_species_weighted_matrix.rds"))
  saveRDS(out, cache_path)
  cat("\nSaved: ", cache_path, "\n", sep = "")
  cat("  Species (nonzero contribution): ", ncol(weighted_matrix),
      " | Active cells: ", nrow(weighted_matrix),
      " | Skipped (zero-area): ", skipped_zero_area, "\n", sep = "")

  # --- Sanity check against the already-cached combined transformed raster ---
  raster_cache <- here::here(
    "data", "processed", "hotspots", "rasters",
    paste0(taxon_key, "_combined_transformed_threat_raster.tif")
  )
  if (file.exists(raster_cache)) {
    combined_idx <- which(logical_flag(species_meta_kept$has_wood_prod_ongoing))
    if (length(combined_idx) > 0) {
      row_sums <- as.numeric(Matrix::rowSums(weighted_matrix[, combined_idx, drop = FALSE]))
      full_vals <- rep(0, terra::ncell(grid_template))
      full_vals[active_cells] <- row_sums
      check_rast <- terra::setValues(terra::rast(grid_template), full_vals)

      cached_rast <- terra::rast(raster_cache)
      diff_vals <- terra::values(check_rast, mat = FALSE) - terra::values(cached_rast, mat = FALSE)
      cat("  Sanity check vs cached combined_transformed_threat_raster -> max abs diff: ",
          format(max(abs(diff_vals), na.rm = TRUE), scientific = TRUE), "\n", sep = "")
    }
  }

  invisible(out)
}

# --- 3. Run --------------------------------------------------------------------
for (taxon in taxa_to_run) {
  cache_path <- file.path(matrix_dir, paste0(taxon, "_cell_species_weighted_matrix.rds"))
  if (file.exists(cache_path) && !rebuild) {
    cat("Already built, skipping: ", taxon, "\n", sep = "")
    next
  }
  build_weighted_matrix(taxon)
  gc(verbose = FALSE)
}

cat("\nAll weighted (cube-root range-weighted) matrices built.\n")
