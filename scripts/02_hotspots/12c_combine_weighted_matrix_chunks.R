args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 12c_combine_weighted_matrix_chunks.R <taxon>")
taxon_key <- args[1]

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix); library(terra); library(here)
})

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000, crs = moll_crs,
  extent = terra::ext(-18040000, 18040000, -9020000, 9020000), vals = 0
)

logical_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

prep_path <- here::here("data", "processed", "hotspots", "null_model_weighted", "_prep",
                         paste0(taxon_key, "_ranges_moll.rds"))
if (!file.exists(prep_path)) stop("Missing prep cache: ", prep_path)
prep <- readRDS(prep_path)
n_total <- nrow(prep$ranges_moll)

chunk_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "_chunks")
chunk_files <- list.files(chunk_dir, pattern = paste0("^", taxon_key, "_chunk_\\d{6}_\\d{6}\\.rds$"), full.names = TRUE)
if (length(chunk_files) == 0) stop("No chunk files found for ", taxon_key, " in ", chunk_dir)

cat("Combining", length(chunk_files), "chunks for", taxon_key, "\n")

all_global_idx <- integer(0)
all_cell_list  <- list()
all_weight_list <- list()

# Verify full coverage of 1:n_total across chunks
covered <- rep(FALSE, n_total)

for (f in chunk_files) {
  ch <- readRDS(f)
  covered[ch$start_idx:ch$end_idx] <- TRUE
  all_global_idx  <- c(all_global_idx, ch$global_idx)
  all_cell_list   <- c(all_cell_list, ch$cell_list)
  all_weight_list <- c(all_weight_list, ch$weight_list)
}

missing_idx <- which(!covered)
if (length(missing_idx) > 0) {
  stop("Chunks do not fully cover species range for ", taxon_key,
       ". Missing indices (first 10): ", paste(head(missing_idx, 10), collapse = ", "),
       ". Total missing: ", length(missing_idx))
}

# Reorder by global index for reproducibility
ord <- order(all_global_idx)
all_global_idx  <- all_global_idx[ord]
all_cell_list   <- all_cell_list[ord]
all_weight_list <- all_weight_list[ord]

species_meta_kept <- prep$species_meta[all_global_idx, , drop = FALSE]

active_cells    <- sort(unique(unlist(all_cell_list, use.names = FALSE)))
nonzero_lengths <- lengths(all_cell_list)
cell_ids        <- unlist(all_cell_list, use.names = FALSE)
row_ids         <- match(cell_ids, active_cells)
col_ids         <- rep(seq_along(all_cell_list), nonzero_lengths)
x_vals          <- unlist(all_weight_list, use.names = FALSE)

weighted_matrix <- Matrix::sparseMatrix(
  i = row_ids, j = col_ids, x = x_vals,
  dims = c(length(active_cells), length(all_cell_list)),
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
  created_at = Sys.time(),
  n_species_total = n_total,
  n_species_kept  = length(all_cell_list)
)

matrix_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "matrices")
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
final_path <- file.path(matrix_dir, paste0(taxon_key, "_cell_species_weighted_matrix.rds"))
saveRDS(out, final_path)
cat("Saved:", final_path, "\n")
cat("  Species (nonzero contribution):", ncol(weighted_matrix),
    "/", n_total, "| Active cells:", nrow(weighted_matrix), "\n")

# --- Sanity check against the already-cached combined transformed raster ---
raster_cache <- here::here("data", "processed", "hotspots", "rasters",
                            paste0(taxon_key, "_combined_transformed_threat_raster.tif"))
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
