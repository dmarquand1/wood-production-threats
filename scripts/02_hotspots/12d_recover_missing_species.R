args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 12d_recover_missing_species.R <taxon>")
taxon_key <- args[1]

# =============================================================================
# RECOVER SMALL-RANGE SPECIES DROPPED BY rasterize(cover=TRUE)
#
# Description:
# terra::rasterize(..., cover = TRUE) silently returns zero coverage for
# polygons much smaller than a single grid cell (110km resolution -> ~12,100
# km^2 per cell); its internal sub-sampling never lands a sample point inside
# very small ranges. This drops exactly the small-range endemics that the
# cube-root weighting is meant to emphasize. This script finds every species
# present in the taxon's range file but absent from the weighted matrix built
# by 12b/12c, and recovers them using an exact vector intersection between the
# species polygon and the grid cells it *touches* (touches=TRUE reliably finds
# the right cell(s) even when cover=TRUE underflows to zero).
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(sf); library(terra); library(Matrix); library(here)
})

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000, crs = moll_crs,
  extent = terra::ext(-18040000, 18040000, -9020000, 9020000), vals = 0
)
res_xy <- terra::res(grid_template)

logical_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

cell_bbox_sf <- function(cid) {
  xy <- terra::xyFromCell(grid_template, cid)
  if (any(!is.finite(xy))) return(NULL)
  xmin <- xy[1, 1] - res_xy[1] / 2; xmax <- xy[1, 1] + res_xy[1] / 2
  ymin <- xy[1, 2] - res_xy[2] / 2; ymax <- xy[1, 2] + res_xy[2] / 2
  coords <- matrix(c(xmin, ymin, xmax, ymin, xmax, ymax, xmin, ymax, xmin, ymin), ncol = 2, byrow = TRUE)
  sf::st_sfc(sf::st_polygon(list(coords)), crs = moll_crs)
}

recover_one <- function(sp_geom, name) {
  tryCatch({
    total_area <- as.numeric(sf::st_area(sp_geom))
    if (!is.finite(total_area) || total_area <= 0) return(NULL)

    touch_r <- terra::rasterize(terra::vect(sp_geom), grid_template, field = 1, background = NA, touches = TRUE)
    touch_vals <- terra::values(touch_r, mat = FALSE)
    touched_cells <- which(!is.na(touch_vals))
    if (length(touched_cells) == 0) return(NULL)

    local_cells <- integer(0); local_weights <- numeric(0)
    for (cid in touched_cells) {
      cell_poly <- cell_bbox_sf(cid)
      if (is.null(cell_poly)) next
      inter <- suppressWarnings(sf::st_intersection(sf::st_geometry(sp_geom), cell_poly))
      inter_area <- if (length(inter) == 0) 0 else sum(as.numeric(sf::st_area(inter)))
      prop <- min(inter_area / total_area, 1)
      if (prop > 0) {
        local_cells <- c(local_cells, cid)
        local_weights <- c(local_weights, prop^(1 / 3))
      }
    }
    if (length(local_cells) == 0) return(NULL)
    list(cells = local_cells, weights = local_weights)
  }, error = function(e) {
    warning("Recovery failed for ", name, " -- ", conditionMessage(e))
    NULL
  })
}

matrix_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "matrices")
mat_path <- file.path(matrix_dir, paste0(taxon_key, "_cell_species_weighted_matrix.rds"))
prep_path <- here::here("data", "processed", "hotspots", "null_model_weighted", "_prep",
                         paste0(taxon_key, "_ranges_moll.rds"))

mat_obj <- readRDS(mat_path)
prep <- readRDS(prep_path)

in_matrix <- colnames(mat_obj$matrix)
all_ids <- prep$species_meta$species_id
missing_idx <- which(!(all_ids %in% in_matrix))
cat(taxon_key, "-- total species:", length(all_ids), "| missing from matrix:", length(missing_idx), "\n")

recovered_cell <- list(); recovered_weight <- list(); recovered_name <- character(0)
still_missing <- character(0)
start_time <- Sys.time()

for (k in seq_along(missing_idx)) {
  i <- missing_idx[k]
  name <- all_ids[i]
  sp_geom <- prep$ranges_moll[i, , drop = FALSE]
  out <- recover_one(sp_geom, name)
  if (is.null(out)) {
    still_missing <- c(still_missing, name)
  } else {
    recovered_cell[[length(recovered_cell) + 1]]     <- out$cells
    recovered_weight[[length(recovered_weight) + 1]] <- out$weights
    recovered_name <- c(recovered_name, name)
  }
  if (k %% 100 == 0 || k == length(missing_idx)) {
    elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    cat("  ", k, "/", length(missing_idx), "(", elapsed, "min )\n")
  }
}

cat("\nRecovered:", length(recovered_name), "| Still missing (genuinely no coverage):", length(still_missing), "\n")
if (length(still_missing) > 0) {
  cat("Still-missing species (likely zero-area / degenerate geometries):\n")
  print(still_missing)
}

if (length(recovered_name) > 0) {
  existing_active <- mat_obj$active_cells
  new_active_needed <- setdiff(unique(unlist(recovered_cell)), existing_active)
  full_active <- sort(c(existing_active, new_active_needed))

  old_summary <- Matrix::summary(mat_obj$matrix)
  old_row_map <- match(existing_active, full_active)
  existing_mat_expanded <- Matrix::sparseMatrix(
    i = old_row_map[old_summary$i], j = old_summary$j, x = old_summary$x,
    dims = c(length(full_active), ncol(mat_obj$matrix)),
    dimnames = list(as.character(full_active), colnames(mat_obj$matrix))
  )

  patch_row_ids <- match(unlist(recovered_cell), full_active)
  patch_col_ids <- rep(seq_along(recovered_cell), lengths(recovered_cell))
  patch_x <- unlist(recovered_weight)
  patch_mat <- Matrix::sparseMatrix(
    i = patch_row_ids, j = patch_col_ids, x = patch_x,
    dims = c(length(full_active), length(recovered_name)),
    dimnames = list(as.character(full_active), recovered_name)
  )

  combined_mat <- cbind(existing_mat_expanded, patch_mat)
  recovered_meta <- prep$species_meta[match(recovered_name, prep$species_meta$species_id), , drop = FALSE]
  combined_species_meta <- rbind(mat_obj$species, recovered_meta)

  mat_obj$matrix <- combined_mat
  mat_obj$active_cells <- full_active
  mat_obj$species <- combined_species_meta
  mat_obj$patched_recovered_species <- recovered_name
  mat_obj$patched_still_missing <- still_missing

  saveRDS(mat_obj, mat_path)
  cat("\nSaved patched matrix:", mat_path, "\n")
  cat("  Species now:", ncol(mat_obj$matrix), "/", length(all_ids),
      "| Active cells now:", nrow(mat_obj$matrix), "\n")
}

# --- Re-run sanity check ---
raster_cache <- here::here("data", "processed", "hotspots", "rasters",
                            paste0(taxon_key, "_combined_transformed_threat_raster.tif"))
if (file.exists(raster_cache)) {
  combined_idx <- which(logical_flag(mat_obj$species$has_wood_prod_ongoing))
  row_sums <- as.numeric(Matrix::rowSums(mat_obj$matrix[, combined_idx, drop = FALSE]))
  full_vals <- rep(0, terra::ncell(grid_template))
  full_vals[mat_obj$active_cells] <- row_sums
  check_rast <- terra::setValues(terra::rast(grid_template), full_vals)
  cached_rast <- terra::rast(raster_cache)
  diff_vals <- terra::values(check_rast, mat = FALSE) - terra::values(cached_rast, mat = FALSE)
  cat("\nRe-checked sanity vs cached combined_transformed_threat_raster -> max abs diff: ",
      format(max(abs(diff_vals), na.rm = TRUE), scientific = TRUE), "\n", sep = "")
}
