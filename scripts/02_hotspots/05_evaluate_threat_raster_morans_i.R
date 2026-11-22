# =============================================================================
# SCRIPT: EVALUATE SPATIAL AUTOCORRELATION IN THREAT RASTERS
#
# Description:
# Computes global Moran's I for the hotspot threat rasters using queen-neighbour
# cell adjacency. By default, the script evaluates rasters whose names end in
# `_threat_raster.tif` or `_transformed_threat_raster.tif`, masks them to land,
# and writes a summary table to outputs/hotspots/tables.
#
# Optional environment variables:
#   MORAN_RASTER_PATTERN   Regex used to select rasters.
#   MORAN_DIRECTIONS       Cell neighbourhood, either 4 or 8. Default: 8.
#   MORAN_PERMUTATIONS     Number of random permutations for p-values. Default: 0.
#   MORAN_SEED             Random seed for permutations. Default: 20260427.
# =============================================================================

rm(list = ls())
gc()

# --- 1. SETUP & LIBRARIES -----------------------------------------------------
local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required_packages <- c("dplyr", "here", "purrr", "sf", "terra", "rnaturalearth", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# --- 2. PATHS & OPTIONS -------------------------------------------------------
raster_dir <- here::here("data", "processed", "hotspots", "rasters")
table_dir <- here::here("outputs", "hotspots", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

output_csv <- file.path(table_dir, "threat_raster_morans_i.csv")

raster_pattern <- Sys.getenv(
  "MORAN_RASTER_PATTERN",
  "(_threat_raster|_transformed_threat_raster)\\.tif$"
)
directions <- as.integer(Sys.getenv("MORAN_DIRECTIONS", "8"))
n_permutations <- as.integer(Sys.getenv("MORAN_PERMUTATIONS", "0"))
moran_seed <- as.integer(Sys.getenv("MORAN_SEED", "20260427"))

if (!directions %in% c(4L, 8L)) {
  stop("MORAN_DIRECTIONS must be either 4 or 8.")
}
if (is.na(n_permutations) || n_permutations < 0) {
  stop("MORAN_PERMUTATIONS must be a non-negative integer.")
}

set.seed(moran_seed)

# --- 3. LAND MASK -------------------------------------------------------------
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name != "Antarctica") %>%
  sf::st_transform(moll_crs)

# --- 4. HELPER FUNCTIONS ------------------------------------------------------
parse_raster_name <- function(path) {
  filename <- basename(path)
  stem <- stringr::str_remove(filename, "\\.tif$")
  taxa <- c("all_vertebrates", "amphibian", "bird", "fish", "mammal", "reptile")
  taxon <- taxa[stringr::str_starts(stem, paste0(taxa, "_"))][1]

  metric <- dplyr::case_when(
    stringr::str_ends(stem, "_transformed_threat_raster") ~ "transformed_threat",
    stringr::str_ends(stem, "_threat_raster") ~ "threat",
    TRUE ~ "unknown"
  )

  threat_category <- stem
  if (!is.na(taxon)) {
    threat_category <- stringr::str_remove(threat_category, paste0("^", taxon, "_"))
  }
  threat_category <- stringr::str_remove(threat_category, "_transformed_threat_raster$")
  threat_category <- stringr::str_remove(threat_category, "_threat_raster$")

  data.frame(
    raster_file = filename,
    taxon = taxon,
    threat_category = threat_category,
    metric = metric,
    stringsAsFactors = FALSE
  )
}

build_land_adjacency <- function(template_raster, land_mask, directions = 8L) {
  land_raster <- terra::rasterize(
    terra::vect(land_mask),
    template_raster,
    field = 1,
    touches = TRUE
  )
  land_cells <- which(!is.na(terra::values(land_raster, mat = FALSE)))

  adjacency <- terra::adjacent(
    land_raster,
    cells = land_cells,
    directions = directions,
    pairs = TRUE,
    include = FALSE
  )

  cell_lookup <- seq_along(land_cells)
  names(cell_lookup) <- land_cells
  from_idx <- unname(cell_lookup[as.character(adjacency[, 1])])
  to_idx <- unname(cell_lookup[as.character(adjacency[, 2])])
  keep <- !is.na(from_idx) & !is.na(to_idx)

  list(
    land_cells = land_cells,
    from_idx = from_idx[keep],
    to_idx = to_idx[keep]
  )
}

calculate_morans_i <- function(raster_path, land_graph, template_raster, permutations = 0L) {
  cat("  ->", basename(raster_path), "\n")

  r <- terra::rast(raster_path)
  if (!terra::compareGeom(r, template_raster, stopOnError = FALSE)) {
    stop("Raster geometry does not match the first raster: ", basename(raster_path))
  }

  vals <- terra::values(r, mat = FALSE)
  land_values <- vals[land_graph$land_cells]
  valid <- !is.na(land_values)
  valid_values <- land_values[valid]

  if (length(valid_values) < 3) {
    return(data.frame(
      n_cells = length(valid_values),
      n_nonzero_cells = sum(valid_values != 0),
      n_adjacencies = 0L,
      mean_value = mean(valid_values),
      variance_value = stats::var(valid_values),
      morans_i = NA_real_,
      expected_i = NA_real_,
      permutation_p_two_sided = NA_real_,
      permutation_p_greater = NA_real_,
      permutation_p_less = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  z_valid <- valid_values - mean(valid_values)
  denominator <- sum(z_valid^2)

  if (denominator == 0) {
    return(data.frame(
      n_cells = length(valid_values),
      n_nonzero_cells = sum(valid_values != 0),
      n_adjacencies = 0L,
      mean_value = mean(valid_values),
      variance_value = stats::var(valid_values),
      morans_i = NA_real_,
      expected_i = -1 / (length(valid_values) - 1),
      permutation_p_two_sided = NA_real_,
      permutation_p_greater = NA_real_,
      permutation_p_less = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  valid_position <- integer(length(land_values))
  valid_position[which(valid)] <- seq_along(valid_values)
  keep <- valid[land_graph$from_idx] & valid[land_graph$to_idx]
  from_idx <- valid_position[land_graph$from_idx[keep]]
  to_idx <- valid_position[land_graph$to_idx[keep]]

  n <- length(valid_values)
  w <- length(from_idx)

  if (w == 0) {
    morans_i <- NA_real_
  } else {
    neighbour_cross_product <- sum(z_valid[from_idx] * z_valid[to_idx])
    morans_i <- (n / w) * (neighbour_cross_product / denominator)
  }

  expected_i <- -1 / (n - 1)

  permutation_p_two_sided <- NA_real_
  permutation_p_greater <- NA_real_
  permutation_p_less <- NA_real_

  if (!is.na(morans_i) && permutations > 0) {
    permuted_i <- numeric(permutations)
    for (i in seq_len(permutations)) {
      z_perm <- sample(z_valid, size = n, replace = FALSE)
      permuted_i[i] <- (n / w) * (sum(z_perm[from_idx] * z_perm[to_idx]) / denominator)
    }

    observed_departure <- abs(morans_i - expected_i)
    permutation_p_two_sided <- (sum(abs(permuted_i - expected_i) >= observed_departure) + 1) /
      (permutations + 1)
    permutation_p_greater <- (sum(permuted_i >= morans_i) + 1) / (permutations + 1)
    permutation_p_less <- (sum(permuted_i <= morans_i) + 1) / (permutations + 1)
  }

  data.frame(
    n_cells = n,
    n_nonzero_cells = sum(valid_values != 0),
    n_adjacencies = w,
    mean_value = mean(valid_values),
    variance_value = stats::var(valid_values),
    morans_i = morans_i,
    expected_i = expected_i,
    permutation_p_two_sided = permutation_p_two_sided,
    permutation_p_greater = permutation_p_greater,
    permutation_p_less = permutation_p_less,
    stringsAsFactors = FALSE
  )
}

# --- 5. RUN ANALYSIS ----------------------------------------------------------
raster_files <- list.files(
  raster_dir,
  pattern = raster_pattern,
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(raster_files) == 0) {
  stop("No rasters matched MORAN_RASTER_PATTERN in ", raster_dir)
}

cat("Evaluating Moran's I for", length(raster_files), "rasters.\n")
cat("Neighbourhood:", directions, "directions\n")
cat("Permutations:", n_permutations, "\n")

template_raster <- terra::rast(raster_files[1])
cat("Building land-cell adjacency graph once from", basename(raster_files[1]), "\n")
land_graph <- build_land_adjacency(
  template_raster = template_raster,
  land_mask = world_map,
  directions = directions
)

results <- purrr::map_dfr(raster_files, function(path) {
  dplyr::bind_cols(
    parse_raster_name(path),
    calculate_morans_i(
      raster_path = path,
      land_graph = land_graph,
      template_raster = template_raster,
      permutations = n_permutations
    )
  )
})

results <- results %>%
  dplyr::arrange(taxon, threat_category, metric)

write.csv(results, output_csv, row.names = FALSE)

cat("\nDone. Moran's I results written to:\n", output_csv, "\n")
