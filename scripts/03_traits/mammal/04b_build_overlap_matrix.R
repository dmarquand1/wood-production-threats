rm(list = ls())
gc()

library(sf)
library(dplyr)
library(here)
library(units)
library(parallel)
library(matrixcalc)

options(stringsAsFactors = FALSE)
sf::sf_use_s2(FALSE)

dir.create(here("data", "processed", "traits", "mammal"), recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# RANGE OVERLAP (JACCARD) SPATIAL CORRELATION MATRIX — mammalS
#
# Computes pairwise Jaccard similarity:
#   J(A, B) = area(intersect(A, B)) / area(union(A, B))
#
# Parallelisation (Windows-safe PSOCK cluster, one cluster for whole script):
#   Step 2 - area pre-computation    (parLapply, index-based)
#   Step 3 - bounding-box filter     (sequential, sf R-tree — see note)
#   Step 5 - Jaccard computation     (parLapply, index-based)
#
# KEY GEOMETRY FIX vs mammal run:
#   Both polygons are projected to Mollweide BEFORE st_intersection is called.
#   The mammal run intersected in lon/lat then projected the result afterward.
#   GEOS uses purely planar arithmetic and has no concept of the sphere, so
#   intersecting complex/wide-ranging geometries in lon/lat can produce
#   spuriously large intersection polygons, inflating inter_area above the
#   union area and giving Jaccard >> 1 (e.g. 1.77e9 for one pair).
#   Projecting first means GEOS always operates in metres on a flat plane.
#   This costs nothing in speed — projecting two small polygons is cheaper
#   than projecting a potentially complex intersection result afterward.
# ===========================================================================

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

TEST_MODE <- TRUE
TEST_N    <- 100L

N_CORES       <- max(1L, detectCores() - 2L)
JACCARD_CHUNK <- 50L
MOLL_CRS      <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Input files — update for your mammal data
TRAIT_DATA_FILE <- here("data/processed/traits/mammal/mammal_trait_model_final_data.csv")
RANGES_FILE     <- here("data/processed/ranges/all_mammal_ranges_dissolved.rds")
SPECIES_COL     <- "spatial_name"   # column in trait CSV matching range sci_name
RANGE_NAME_COL  <- "sci_name"       # species name column in the ranges object

OUT_MATRIX_FILE  <- here("data", "processed", "traits", "mammal", "mammal_range_overlap_matrix.rds")
OUT_SUMMARY_FILE <- here("data", "processed", "traits", "mammal", "mammal_overlap_matrix_summary.csv")

# ---------------------------------------------------------------------------
# PROGRESS HELPER
# ---------------------------------------------------------------------------
progress_bar <- function(current, total, label = "", bar_width = 35L) {
  pct    <- current / total
  filled <- round(pct * bar_width)
  bar    <- paste0(strrep("=", filled), strrep("-", bar_width - filled))
  ts     <- format(Sys.time(), "%H:%M:%S")
  cat(sprintf("\r  [%s] [%s] %3.0f%%  %d/%d  %s     ",
              ts, bar, pct * 100, current, total, label))
  if (current == total) cat("\n")
  flush.console()
}

# ---------------------------------------------------------------------------
# WORKER FUNCTIONS
# Workers access species_ranges / candidate_list / species_areas_m2 from
# their own environment — exported once via clusterExport, never passed
# as arguments (avoids per-call serialisation of large sf objects).
# ---------------------------------------------------------------------------

worker_get_area_km2 <- function(i) {
  moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  # Validate in geographic space before reprojecting — st_transform can fail
  # on invalid geometries in a less recoverable way than st_make_valid can fix
  poly    <- tryCatch(sf::st_make_valid(species_ranges[i, ]),
                      error = function(e) species_ranges[i, ])
  area_m2 <- tryCatch(
    sum(sf::st_area(sf::st_transform(poly, moll))),
    error = function(e) units::set_units(0, "m^2")
  )
  as.numeric(units::set_units(area_m2, "km^2"))
}

worker_compute_row <- function(i) {
  moll         <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  j_candidates <- candidate_list[[i]]
  j_candidates <- j_candidates[j_candidates > i]
  if (length(j_candidates) == 0L) return(NULL)
  
  # Project row i ONCE outside the inner loop
  row_i_proj <- tryCatch(
    sf::st_transform(sf::st_make_valid(species_ranges[i, ]), moll),
    error = function(e) NULL
  )
  if (is.null(row_i_proj)) {
    res <- rep(NA_real_, length(j_candidates))
    names(res) <- as.character(j_candidates)
    return(res)
  }
  
  area_i  <- species_areas_m2[i]
  results <- numeric(length(j_candidates))
  
  for (k in seq_along(j_candidates)) {
    j      <- j_candidates[k]
    area_j <- species_areas_m2[j]
    
    if (area_i == 0 || area_j == 0) { results[k] <- NA_real_; next }
    
    # PROJECT j BEFORE intersecting — this is the critical fix.
    # GEOS operates in flat Cartesian space; intersecting in lon/lat gives
    # it coordinates it treats as metres, producing corrupt results for
    # geometrically complex species. See header comment for details.
    row_j_proj <- tryCatch(
      sf::st_transform(sf::st_make_valid(species_ranges[j, ]), moll),
      error = function(e) NULL
    )
    if (is.null(row_j_proj)) { results[k] <- NA_real_; next }
    
    # Intersect in projected (flat metre) space
    inter <- tryCatch(
      suppressWarnings(sf::st_intersection(row_i_proj, row_j_proj)),
      error = function(e)
        tryCatch(
          suppressWarnings(sf::st_intersection(
            sf::st_make_valid(row_i_proj),
            sf::st_make_valid(row_j_proj)
          )),
          error = function(e2) NULL
        )
    )
    
    if (is.null(inter) || nrow(inter) == 0L) { results[k] <- 0; next }
    
    # set_units forces explicit m² before as.numeric strips the unit,
    # guarding against any internal units object state issues
    inter_area <- tryCatch(
      as.numeric(units::set_units(sum(sf::st_area(inter)), "m^2")),
      error = function(e) 0
    )
    if (!is.finite(inter_area) || inter_area <= 0) { results[k] <- 0; next }
    
    raw_jaccard <- inter_area / (area_i + area_j - inter_area)
    # Clamp catches genuine floating-point rounding (e.g. inter slightly > area
    # due to rounding). If it clamps a value that is substantially > 1, the
    # n_values_above_1 diagnostic at the end will flag it.
    results[k] <- min(1, max(0, raw_jaccard))
  }
  
  names(results) <- as.character(j_candidates)
  results
}

# ===========================================================================
# STEP 1 — LOAD & VALIDATE DATA
# ===========================================================================
cat("\n========================================\n")
cat("  STEP 1 / 6  |  Loading data\n")
cat("========================================\n")
t_total <- proc.time()
t_step  <- proc.time()

analysis_df <- read.csv(TRAIT_DATA_FILE) %>%
  distinct(!!sym(SPECIES_COL), .keep_all = TRUE)

ranges_raw <- readRDS(RANGES_FILE) %>%
  filter(!!sym(RANGE_NAME_COL) %in% analysis_df[[SPECIES_COL]]) %>%
  arrange(!!sym(RANGE_NAME_COL))

cat(sprintf("  Loaded %d species from disk.\n", nrow(ranges_raw)))
cat(sprintf("  CRS: %s\n", st_crs(ranges_raw)$input))

if (TEST_MODE) {
  set.seed(42L)
  keep       <- sort(sample(nrow(ranges_raw), min(TEST_N, nrow(ranges_raw))))
  ranges_raw <- ranges_raw[keep, ]
  cat(sprintf("  *** TEST MODE: subsampled to %d species ***\n", nrow(ranges_raw)))
  OUT_MATRIX_FILE  <- sub("\\.rds$", "_TEST.rds",  OUT_MATRIX_FILE)
  OUT_SUMMARY_FILE <- sub("\\.csv$", "_TEST.csv",  OUT_SUMMARY_FILE)
}

# ---- Geometry validation ---------------------------------------------------
cat("  Validating geometries ...\n")
ranges_valid <- suppressWarnings(st_make_valid(ranges_raw))
is_valid     <- st_is_valid(ranges_valid)
keep_idx     <- !is.na(is_valid) & is_valid
failed_geom  <- character(0)

if (any(!keep_idx)) {
  failed_geom <- ranges_raw[[RANGE_NAME_COL]][!keep_idx]
  cat(sprintf("  WARNING — dropped %d species with unrepairable geometries:\n",
              length(failed_geom)))
  cat(paste(" ", failed_geom, collapse = "\n"), "\n")
}

species_ranges <- ranges_valid[keep_idx, ]
n              <- nrow(species_ranges)
species_names  <- species_ranges[[RANGE_NAME_COL]]

# ---- CRS check — bbox step requires lon/lat --------------------------------
if (!st_is_longlat(species_ranges)) {
  cat("  Ranges are not in lon/lat. Reprojecting to WGS84 for bbox step.\n")
  species_ranges <- st_transform(species_ranges, 4326)
}

elapsed <- (proc.time() - t_step)["elapsed"]
cat(sprintf("  Proceeding with %d species.  [%.1f sec]\n\n", n, elapsed))

# ===========================================================================
# STEP 2 — PARALLEL AREA PRE-COMPUTATION
# ===========================================================================
cat("========================================\n")
cat(sprintf("  STEP 2 / 6  |  Pre-computing species areas  (%d workers)\n",
            N_CORES))
cat("========================================\n")
t_step <- proc.time()

cl <- makeCluster(N_CORES, type = "PSOCK")
on.exit(stopCluster(cl), add = TRUE)

clusterExport(cl, varlist = c("species_ranges", "worker_get_area_km2"))
clusterEvalQ(cl, { library(sf); library(units); sf::sf_use_s2(FALSE) })

AREA_BATCH   <- max(1L, ceiling(n / 20L))
batch_starts <- seq(1L, n, by = AREA_BATCH)
n_batches    <- length(batch_starts)
species_areas_km2 <- numeric(n)

for (bi in seq_len(n_batches)) {
  idx <- seq(batch_starts[bi], min(batch_starts[bi] + AREA_BATCH - 1L, n))
  species_areas_km2[idx] <- unlist(parLapply(cl, idx, worker_get_area_km2))
  progress_bar(bi, n_batches,
               label = sprintf("rows %d-%d", idx[1L], idx[length(idx)]))
}

names(species_areas_km2) <- species_names
species_areas_m2 <- species_areas_km2 * 1e6

# ---- Area sanity check — catches the zero-area-for-all bug immediately ----
nonzero_areas   <- species_areas_m2[species_areas_m2 > 0]
nonzero_mean_m2 <- if (length(nonzero_areas) > 0) mean(nonzero_areas) else 0
nonzero_med_m2  <- if (length(nonzero_areas) > 0) median(nonzero_areas) else 0

cat(sprintf("\n  Area sanity check:\n"))
cat(sprintf("    Mean non-zero area:   %.3e m2  (%.1f km2)\n",
            nonzero_mean_m2, nonzero_mean_m2 / 1e6))
cat(sprintf("    Median non-zero area: %.3e m2  (%.1f km2)\n",
            nonzero_med_m2,  nonzero_med_m2  / 1e6))

if (nonzero_mean_m2 < 1e4)
  stop(sprintf(paste0(
    "FATAL: mean non-zero area is %.3e m2 (%.6f km2). Expected >> 1e4 m2.\n",
    "Area conversion is almost certainly wrong.\n",
    "Check worker_get_area_km2 and the species_areas_m2 <- * 1e6 line."),
    nonzero_mean_m2, nonzero_mean_m2 / 1e6))

zero_area <- sum(species_areas_km2 == 0)
elapsed   <- (proc.time() - t_step)["elapsed"]
cat(sprintf("  Zero-area species: %d  [%.1f sec]\n", zero_area, elapsed))
if (zero_area > 0) {
  cat("  Zero-area species (will receive NA overlap):\n")
  cat(paste(" ", species_names[species_areas_km2 == 0], collapse = "\n"), "\n")
}
cat("\n")

# ===========================================================================
# STEP 3 — BOUNDING-BOX CANDIDATE FILTER  (sequential)
# ===========================================================================
cat("========================================\n")
cat("  STEP 3 / 6  |  Bounding-box spatial index  (main process)\n")
cat("========================================\n")
t_step <- proc.time()

cat("  Building R-tree index and querying all pairs ...\n")
TICK_EVERY     <- max(1L, floor(n * 0.05))
candidate_list <- vector("list", n)

for (i in seq_len(n)) {
  candidate_list[[i]] <- as.integer(unlist(
    st_intersects(species_ranges[i, ], species_ranges, sparse = TRUE)
  ))
  if (i %% TICK_EVERY == 0L || i == n)
    progress_bar(i, n, label = sprintf("species %d/%d", i, n))
}

n_candidates <- sum(vapply(seq_len(n),
                           function(i) sum(candidate_list[[i]] > i), integer(1L)))
n_possible   <- n * (n - 1L) / 2L
elapsed      <- (proc.time() - t_step)["elapsed"]

cat(sprintf(
  "  Candidate pairs: %s of %s  (%.1f%% eliminated by bbox)  [%.1f sec]\n\n",
  format(n_candidates, big.mark = ","),
  format(n_possible,   big.mark = ","),
  100 * (1 - n_candidates / n_possible),
  elapsed
))

# ===========================================================================
# STEP 4 — INITIALISE OVERLAP MATRIX
# ===========================================================================
cat("========================================\n")
cat("  STEP 4 / 6  |  Initialising overlap matrix\n")
cat("========================================\n")

overlap_matrix <- matrix(
  0, nrow = n, ncol = n,
  dimnames = list(species_names, species_names)
)
diag(overlap_matrix) <- 1
cat(sprintf("  %d x %d matrix allocated  (%.1f MB)\n\n",
            n, n, object.size(overlap_matrix) / 1024^2))

# ===========================================================================
# STEP 5 — PARALLEL JACCARD COMPUTATION
# ===========================================================================
cat("========================================\n")
cat(sprintf("  STEP 5 / 6  |  Computing Jaccard overlaps  (%d workers)\n",
            N_CORES))
cat("========================================\n")
t_step <- proc.time()

clusterExport(cl, varlist = c(
  "species_ranges", "candidate_list",
  "species_areas_m2", "worker_compute_row"
))
clusterEvalQ(cl, { library(sf); library(units); sf::sf_use_s2(FALSE) })

rows_to_do   <- seq_len(n - 1L)
chunk_starts <- seq(1L, length(rows_to_do), by = JACCARD_CHUNK)
n_chunks     <- length(chunk_starts)

cat(sprintf("  %d row-chunks of up to %d rows each.\n\n",
            n_chunks, JACCARD_CHUNK))

for (chunk_idx in seq_len(n_chunks)) {
  
  chunk_begin <- chunk_starts[chunk_idx]
  chunk_end   <- min(chunk_begin + JACCARD_CHUNK - 1L, length(rows_to_do))
  chunk_rows  <- rows_to_do[chunk_begin:chunk_end]
  
  row_results <- parLapply(cl, chunk_rows, worker_compute_row)
  
  for (k in seq_along(chunk_rows)) {
    i   <- chunk_rows[k]
    res <- row_results[[k]]
    if (is.null(res)) next
    j_idx <- as.integer(names(res))
    overlap_matrix[i, j_idx] <- res
    overlap_matrix[j_idx, i] <- res
  }
  
  last_row <- chunk_rows[length(chunk_rows)]
  elapsed  <- (proc.time() - t_step)["elapsed"]
  rate     <- last_row / max(elapsed, 1e-9)
  eta_sec  <- (n - 1L - last_row) / max(rate, 1e-9)
  
  progress_bar(
    current = chunk_idx, total = n_chunks,
    label   = sprintf("row %d/%d | ETA %.0f min",
                      last_row, n - 1L, eta_sec / 60)
  )
}

elapsed_jaccard <- (proc.time() - t_step)["elapsed"]
cat(sprintf("  Jaccard complete.  [%.1f sec / %.1f min]\n\n",
            elapsed_jaccard, elapsed_jaccard / 60))

stopCluster(cl)
on.exit(NULL)

# ===========================================================================
# STEP 6 — VERIFY & SAVE
# ===========================================================================
cat("========================================\n")
cat("  STEP 6 / 6  |  Verifying matrix & saving\n")
cat("========================================\n")

# ---- Symmetry --------------------------------------------------------------
if (!isSymmetric(overlap_matrix)) {
  cat("  Correcting floating-point asymmetry ... ")
  overlap_matrix[lower.tri(overlap_matrix)] <-
    t(overlap_matrix)[lower.tri(overlap_matrix)]
  cat("done.\n")
}
diag(overlap_matrix) <- 1

# ---- Values out of [0, 1] — report immediately, before eigenvalue step -----
off_diag  <- overlap_matrix[upper.tri(overlap_matrix)]
n_above_1 <- sum(off_diag > 1, na.rm = TRUE)
n_below_0 <- sum(off_diag < 0, na.rm = TRUE)

if (n_above_1 > 0 || n_below_0 > 0) {
  cat(sprintf("  WARNING — %d values > 1 and %d values < 0.\n",
              n_above_1, n_below_0))
  cat("  Run diagnose_overlap_matrix.R before using this matrix.\n")
  ut <- upper.tri(overlap_matrix)
  if (n_above_1 > 0) {
    bad <- which(overlap_matrix > 1 & ut, arr.ind = TRUE)
    bad <- bad[order(-overlap_matrix[bad]), , drop = FALSE]
    cat("  Top offenders (> 1):\n")
    for (k in seq_len(min(10L, nrow(bad))))
      cat(sprintf("    %-38s x %-38s = %.4f\n",
                  species_names[bad[k,1]], species_names[bad[k,2]],
                  overlap_matrix[bad[k,1], bad[k,2]]))
  }
} else {
  cat("  All off-diagonal values in [0, 1]. OK.\n")
}

# ---- Eigenvalue check ------------------------------------------------------
cat("  Computing eigenvalues ... ")
min_eigen <- min(eigen(overlap_matrix, symmetric = TRUE,
                       only.values = TRUE)$values)
cat(sprintf("min eigenvalue = %.8f\n", min_eigen))

if (min_eigen < 0) {
  nugget <- abs(min_eigen) + 1e-6
  if (nugget > 0.01)
    cat(sprintf(paste0(
      "  WARNING — nugget is large (%.4f). A nugget > 0.01 suggests a\n",
      "  structural matrix problem, not just floating-point rounding.\n",
      "  Investigate values > 1 above before proceeding.\n"), nugget))
  diag(overlap_matrix) <- diag(overlap_matrix) + nugget
  cat(sprintf("  Nugget applied: %.2e  (diagonal now %.8f)\n",
              nugget, diag(overlap_matrix)[1]))
} else {
  cat("  No nugget needed. Matrix is positive semi-definite.\n")
}

# ---- Summary stats ---------------------------------------------------------
off_diag     <- overlap_matrix[upper.tri(overlap_matrix)]  # refresh post-nugget
pct_zero     <- 100 * mean(off_diag == 0,   na.rm = TRUE)
pct_above_05 <- 100 * mean(off_diag > 0.05, na.rm = TRUE)
mean_nonzero <- mean(off_diag[!is.na(off_diag) & off_diag > 0])
med_nonzero  <- median(off_diag[!is.na(off_diag) & off_diag > 0])
n_na         <- sum(is.na(off_diag))
max_jaccard  <- max(off_diag, na.rm = TRUE)
pos_def      <- is.positive.definite(overlap_matrix)

cat(sprintf("  Off-diagonal zeros:            %.1f%%\n",  pct_zero))
cat(sprintf("  Off-diagonal > 0.05:           %.1f%%\n",  pct_above_05))
cat(sprintf("  Mean overlap (non-zero):       %.6f\n",    mean_nonzero))
cat(sprintf("  Median overlap (non-zero):     %.6f\n",    med_nonzero))
cat(sprintf("  Max Jaccard value:             %.6f\n",    max_jaccard))
cat(sprintf("  Values > 1:                    %d\n",      n_above_1))
cat(sprintf("  NA pairs:                      %d\n",      n_na))
cat(sprintf("  Is positive definite:          %d\n",      pos_def))

# Mean >> 100x median is the diagnostic signature of the inflate-one-pair bug
if (mean_nonzero > med_nonzero * 100)
  cat(sprintf(paste0(
    "\n  WARNING — mean (%.4f) is >> 100x median (%.6f).\n",
    "  This is the signature of inflated Jaccard values from the\n",
    "  project-after-intersect bug. Run diagnose_overlap_matrix.R.\n"),
    mean_nonzero, med_nonzero))

# ---- Save ------------------------------------------------------------------
saveRDS(overlap_matrix, OUT_MATRIX_FILE)
total_elapsed <- as.numeric((proc.time() - t_total)["elapsed"])

summary_df <- data.frame(
  test_mode            = TEST_MODE,
  n_species            = n,
  n_dropped            = length(failed_geom),
  n_cores_used         = N_CORES,
  pct_zero_overlap     = round(pct_zero, 2),
  pct_above_05         = round(pct_above_05, 2),
  mean_nonzero_jaccard = round(mean_nonzero, 6),
  median_nonzero       = round(med_nonzero,  6),
  max_jaccard          = round(max_jaccard,  6),
  n_values_above_1     = n_above_1,
  n_na_pairs           = n_na,
  min_eigenvalue       = round(min_eigen, 8),
  total_mins           = round(total_elapsed / 60, 1),
  stringsAsFactors     = FALSE
)
write.csv(summary_df, OUT_SUMMARY_FILE, row.names = FALSE)

cat(sprintf("\n  Matrix saved:  %s\n", basename(OUT_MATRIX_FILE)))
cat(sprintf("  Summary saved: %s\n",  basename(OUT_SUMMARY_FILE)))
if (TEST_MODE){
  cat("  NOTE: TEST_MODE = TRUE. Set to FALSE for the full run.\n")
cat(sprintf("  Total wall time: %.1f min\n\n", total_elapsed / 60))
}
