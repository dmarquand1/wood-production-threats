args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript 12a_prep_weighted_ranges.R <taxon>")
taxon_key <- args[1]

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(sf); library(here)
})

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

range_files <- c(
  amphibian = here::here("data", "processed", "ranges", "all_amphibian_ranges_dissolved.rds"),
  bird      = here::here("data", "processed", "ranges", "all_bird_ranges_dissolved.rds"),
  fish      = here::here("data", "processed", "ranges", "all_fish_ranges_dissolved.rds"),
  mammal    = here::here("data", "processed", "ranges", "all_mammal_ranges_dissolved.rds"),
  reptile   = here::here("data", "processed", "ranges", "all_reptile_ranges_dissolved.rds")
)

prep_dir <- here::here("data", "processed", "hotspots", "null_model_weighted", "_prep")
dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
prep_path <- file.path(prep_dir, paste0(taxon_key, "_ranges_moll.rds"))

if (file.exists(prep_path)) {
  cat("Already prepped:", prep_path, "\n")
  quit(status = 0)
}

load_assessment_species <- function(taxon_key) {
  path <- here::here("data", "processed", "assessments", paste0(taxon_key, "_assessment_species_summary.rds"))
  if (!file.exists(path)) stop("Missing assessment species summary: ", path)
  dat <- readRDS(path)
  dat <- dat[!dat$category %in% c("EX", "EW"), , drop = FALSE]
  dat <- dat[!is.na(dat$scientific_name) & nzchar(dat$scientific_name), , drop = FALSE]
  dat[!duplicated(dat$scientific_name), , drop = FALSE]
}

cat("Prepping:", taxon_key, "\n")
assessment <- load_assessment_species(taxon_key)
ranges <- readRDS(range_files[[taxon_key]])
if (!inherits(ranges, "sf")) stop("Range object is not sf: ", range_files[[taxon_key]])

ranges <- ranges[!is.na(ranges$sci_name) & ranges$sci_name %in% assessment$scientific_name, , drop = FALSE]
ranges <- ranges[!duplicated(ranges$sci_name), , drop = FALSE]
if (nrow(ranges) == 0) stop("No extant assessed species with ranges for ", taxon_key)

species_meta <- assessment[match(ranges$sci_name, assessment$scientific_name), , drop = FALSE]
species_meta$species_id <- ranges$sci_name

cat("Reprojecting", nrow(ranges), "ranges to Mollweide...\n")
ranges_moll <- sf::st_transform(ranges, moll_crs)

saveRDS(list(ranges_moll = ranges_moll, species_meta = species_meta), prep_path)
cat("Saved:", prep_path, "| n =", nrow(ranges_moll), "\n")
