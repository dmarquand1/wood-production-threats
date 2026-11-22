args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript 12b_build_weighted_matrix_chunk.R <taxon> <start_idx> <end_idx>")
taxon_key <- args[1]
start_idx <- as.integer(args[2])
end_idx   <- as.integer(args[3])

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(sf); library(terra); library(here)
})

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
grid_template <- terra::rast(
  resolution = 110000, crs = moll_crs,
  extent = terra::ext(-18040000, 18040000, -9020000, 9020000), vals = 0
)
cell_area <- terra::res(grid_template)[1] * terra::res(grid_template)[2]
PER_SPECIES_TIMEOUT_SECS <- 30

prep_path <- here::here("data", "processed", "hotspots", "null_model_weighted", "_prep",
                         paste0(taxon_key, "_ranges_moll.rds"))
if (!file.exists(prep_path)) stop("Missing prep cache: ", prep_path, ". Run 12a_prep_weighted_ranges.R first.")

chunk_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "_chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)
chunk_path <- file.path(chunk_dir, sprintf("%s_chunk_%06d_%06d.rds", taxon_key, start_idx, end_idx))

if (file.exists(chunk_path)) {
  cat("Chunk already built:", chunk_path, "\n")
  quit(status = 0)
}

cat("Loading prep cache...\n")
prep <- readRDS(prep_path)
ranges_moll <- prep$ranges_moll
n_total <- nrow(ranges_moll)
end_idx <- min(end_idx, n_total)

cat(sprintf("Processing %s indices %d-%d of %d\n", taxon_key, start_idx, end_idx, n_total))

idx_seq <- start_idx:end_idx
cell_list   <- vector("list", length(idx_seq))
weight_list <- vector("list", length(idx_seq))
start_time <- Sys.time()

process_one <- function(sp_geom, name, idx) {
  tryCatch({
    setTimeLimit(elapsed = PER_SPECIES_TIMEOUT_SECS, transient = TRUE)

    total_area <- as.numeric(sf::st_area(sp_geom))
    result <- NULL
    if (is.finite(total_area) && total_area > 0) {
      cov <- try(terra::rasterize(terra::vect(sp_geom), grid_template, cover = TRUE), silent = TRUE)
      if (inherits(cov, "try-error")) {
        sp_geom_fix <- sf::st_make_valid(sp_geom)
        cov <- terra::rasterize(terra::vect(sp_geom_fix), grid_template, cover = TRUE)
      }
      prop <- (cov * cell_area) / total_area
      vals <- terra::values(prop, mat = FALSE)
      nz <- which(!is.na(vals) & vals > 0)
      if (length(nz) > 0) {
        vals_nz <- vals[nz]
        vals_nz[vals_nz > 1] <- 1
        result <- list(nz = nz, w = vals_nz^(1 / 3))
      }
    }
    setTimeLimit(elapsed = Inf, transient = TRUE)
    result
  }, error = function(e) {
    setTimeLimit(elapsed = Inf, transient = TRUE)
    warning("Skipping species (timeout/error) [", idx, "] ", name, " -- ", conditionMessage(e))
    NULL
  })
}

for (k in seq_along(idx_seq)) {
  i <- idx_seq[k]
  sp_geom <- ranges_moll[i, , drop = FALSE]
  name <- sp_geom$sci_name

  out <- process_one(sp_geom, name, i)
  if (!is.null(out)) {
    cell_list[[k]]   <- out$nz
    weight_list[[k]] <- out$w
  }

  if (k %% 200 == 0 || k == length(idx_seq)) {
    elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
    cat(sprintf("  %d / %d in chunk (%.1f min elapsed)\n", k, length(idx_seq), elapsed))
  }
}

has_data <- !vapply(cell_list, is.null, logical(1))
chunk_out <- list(
  taxon        = taxon_key,
  start_idx    = start_idx,
  end_idx      = end_idx,
  global_idx   = idx_seq[has_data],
  cell_list    = cell_list[has_data],
  weight_list  = weight_list[has_data]
)
saveRDS(chunk_out, chunk_path)
cat("Saved chunk:", chunk_path, "| species with data:", sum(has_data), "/", length(idx_seq), "\n")
