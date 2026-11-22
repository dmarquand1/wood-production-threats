rm(list = ls())
gc()

# ==============================================================================
# SMALL-RANGE SPECIES HOTSPOT MAPS FOR TRAIT MODEL INTERPRETATION
#
# Purpose:
#   Map where modelled species with the smallest ranges are concentrated. This is
#   intended to help interpret why range size can have different signs in the
#   trait models across taxa.
#
# Definition:
#   "Small-ranged" species are the lowest SMALL_RANGE_QUANTILE of log_range_size
#   within each taxon. The model data use taxon-specific scaled log range size, so
#   this script also uses taxon-specific thresholds.
#
# Outputs:
#   data/processed/traits/small_range_hotspots/
#     - *_model_species_richness.tif
#     - *_small_range_q25_richness.tif
#     - *_small_range_q25_proportion.tif
#
#   outputs/traits/figures/small_range_hotspots/
#     - individual count/proportion maps per taxon
#     - combined count/proportion panels
#
#   outputs/traits/tables/small_range_hotspots/
#     - small_range_hotspot_summary.csv
#     - small_range_species_list.csv
# ==============================================================================

# --- 1. SETUP -----------------------------------------------------------------

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

required_pkgs <- c(
  "dplyr",
  "ggplot2",
  "here",
  "patchwork",
  "rnaturalearth",
  "scales",
  "sf",
  "terra",
  "tibble",
  "viridis"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Install these packages before running this script: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(rnaturalearth)
  library(sf)
  library(terra)
  library(viridis)
})

SMALL_RANGE_QUANTILE <- 0.25
REBUILD_RASTERS <- FALSE
RASTER_CHUNK_SIZE <- 250

quantile_tag <- paste0("q", round(SMALL_RANGE_QUANTILE * 100))

taxa_config <- tibble::tribble(
  ~taxon,      ~label,              ~range_file,
  "amphibian", "Amphibians",        "data/processed/ranges/all_amphibian_ranges_dissolved.rds",
  "bird",      "Birds",             "data/processed/ranges/all_bird_ranges_dissolved.rds",
  "mammal",    "Mammals",           "data/processed/ranges/all_mammal_ranges_dissolved.rds",
  "reptile",   "Reptiles",          "data/processed/ranges/all_reptile_ranges_dissolved.rds",
  "fish",      "Freshwater Fish",   "data/processed/ranges/all_fish_ranges_dissolved.rds"
)

moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent <- terra::ext(-18040000, 18040000, -9020000, 9020000)
grid_template <- terra::rast(
  resolution = 110000,
  crs = moll_crs,
  extent = global_extent,
  vals = 0
)

world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>%
  st_transform(moll_crs)

land_vect <- terra::vect(world_map)

raster_dir <- here::here("data", "processed", "traits", "small_range_hotspots")
figure_dir <- here::here("outputs", "traits", "figures", "small_range_hotspots")
table_dir <- here::here("outputs", "traits", "tables", "small_range_hotspots")

dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2. HELPERS ----------------------------------------------------------------

clean_species_names <- function(x) {
  unique(trimws(as.character(x[!is.na(x) & x != ""])))
}

read_model_data <- function(taxon) {
  path <- here::here(
    "data",
    "processed",
    "traits",
    taxon,
    paste0(taxon, "_trait_model_final_data.csv")
  )

  if (!file.exists(path)) {
    stop("Missing model data file: ", path)
  }

  model_data <- read.csv(path, stringsAsFactors = FALSE)

  required_cols <- c("spatial_name", "phylo_name", "log_range_size")
  missing_cols <- setdiff(required_cols, names(model_data))
  if (length(missing_cols) > 0) {
    stop(
      "Model data for ",
      taxon,
      " is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  model_data %>%
    mutate(
      taxon = taxon,
      spatial_name = trimws(spatial_name),
      phylo_name = trimws(phylo_name),
      log_range_size = as.numeric(log_range_size)
    ) %>%
    filter(!is.na(spatial_name), spatial_name != "", !is.na(log_range_size))
}

calculate_or_load_richness <- function(species_names, ranges_sf, output_path, label) {
  if (file.exists(output_path) && !REBUILD_RASTERS) {
    message("  Loading cached raster: ", basename(output_path))
    cached_raster <- tryCatch(
      terra::rast(output_path),
      error = function(e) {
        warning(
          "Cached raster could not be loaded and will be rebuilt: ",
          basename(output_path),
          " (",
          conditionMessage(e),
          ")",
          call. = FALSE
        )
        NULL
      }
    )

    if (!is.null(cached_raster)) {
      return(cached_raster)
    }
  }

  message("  Rasterising ", label, "...")

  species_names <- clean_species_names(species_names)
  species_sf <- ranges_sf %>%
    filter(sci_name %in% species_names)

  out_raster <- terra::rast(grid_template, vals = 0)

  if (nrow(species_sf) == 0) {
    warning("No range polygons matched for ", label, call. = FALSE)
    terra::writeRaster(out_raster, output_path, overwrite = TRUE)
    return(out_raster)
  }

  species_sf <- species_sf %>%
    mutate(.richness_count = 1)

  rasterize_chunk <- function(chunk_sf, chunk_label) {
    empty_raster <- terra::rast(grid_template, vals = 0)

    if (nrow(chunk_sf) == 0) {
      return(empty_raster)
    }

    tryCatch(
      {
        chunk_vect <- chunk_sf %>%
          st_transform(moll_crs) %>%
          terra::vect()

        chunk_raster <- terra::rasterize(
          chunk_vect,
          grid_template,
          field = ".richness_count",
          fun = "sum",
          background = 0
        )

        chunk_raster[is.na(chunk_raster)] <- 0
        chunk_raster
      },
      error = function(e) {
        if (nrow(chunk_sf) == 1) {
          warning(
            "Skipping one range polygon during rasterisation for ",
            chunk_label,
            ": ",
            conditionMessage(e),
            call. = FALSE
          )
          return(empty_raster)
        }

        message(
          "    Chunk failed for ",
          chunk_label,
          " (",
          conditionMessage(e),
          "). Splitting into smaller chunks..."
        )

        split_at <- floor(nrow(chunk_sf) / 2)
        left_raster <- rasterize_chunk(chunk_sf[seq_len(split_at), ], chunk_label)
        right_raster <- rasterize_chunk(
          chunk_sf[(split_at + 1):nrow(chunk_sf), ],
          chunk_label
        )

        left_raster + right_raster
      }
    )
  }

  total_species <- nrow(species_sf)
  chunk_starts <- seq(1, total_species, by = RASTER_CHUNK_SIZE)

  for (chunk_index in seq_along(chunk_starts)) {
    start_i <- chunk_starts[[chunk_index]]
    end_i <- min(start_i + RASTER_CHUNK_SIZE - 1, total_species)

    message(
      "    Chunk ",
      chunk_index,
      "/",
      length(chunk_starts),
      " (species ",
      start_i,
      "-",
      end_i,
      " of ",
      total_species,
      ")"
    )

    chunk_sf <- species_sf[start_i:end_i, ]
    chunk_raster <- rasterize_chunk(chunk_sf, label)
    out_raster <- out_raster + chunk_raster

    rm(chunk_sf, chunk_raster)
    gc()
  }

  out_raster[is.na(out_raster)] <- 0
  terra::writeRaster(out_raster, output_path, overwrite = TRUE)
  out_raster
}

safe_divide <- function(numerator, denominator) {
  out <- numerator / denominator
  out[denominator <= 0] <- NA
  out[is.nan(out)] <- NA
  out[is.infinite(out)] <- NA
  out
}

get_positive_land_values <- function(r) {
  r_land <- terra::mask(r, land_vect)
  vals <- terra::values(r_land, mat = FALSE)
  vals[vals > 0 & !is.na(vals)]
}

hotspot_thresholds <- function(r) {
  vals <- get_positive_land_values(r)

  if (length(vals) == 0) {
    return(c(top_10 = NA_real_, top_5 = NA_real_, top_1 = NA_real_))
  }

  unname(quantile(vals, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)) %>%
    setNames(c("top_10", "top_5", "top_1"))
}

plot_small_range_hotspots <- function(r, title_txt, legend_title = "Small-range\nrichness") {
  r_land <- terra::mask(r, land_vect)
  vals <- terra::values(r_land, mat = FALSE)
  vals <- vals[vals > 0 & !is.na(vals)]

  if (length(vals) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = title_txt) +
        annotate("text", x = 0, y = 0, label = "No matched small-range species")
    )
  }

  q <- quantile(vals, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)

  plot_df <- as.data.frame(r_land, xy = TRUE) %>%
    setNames(c("x", "y", "value")) %>%
    filter(value > 0 & !is.na(value)) %>%
    mutate(
      category = case_when(
        value >= q[[3]] ~ "Top 1%",
        value >= q[[2]] ~ "Top 5%",
        value >= q[[1]] ~ "Top 10%",
        TRUE ~ "Present"
      ),
      category = factor(category, levels = c("Present", "Top 10%", "Top 5%", "Top 1%"))
    )

  hotspot_cols <- c(
    "Present" = "grey82",
    "Top 10%" = "#FEE08B",
    "Top 5%" = "#F46D43",
    "Top 1%" = "#9E0142"
  )

  ggplot() +
    geom_sf(data = world_map, fill = "#f2f2f2", color = NA) +
    geom_tile(data = plot_df, aes(x = x, y = y, fill = category)) +
    scale_fill_manual(values = hotspot_cols, name = legend_title, drop = FALSE) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    labs(title = title_txt) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(0.7, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

plot_proportion_gradient <- function(r, title_txt, legend_title = "Small-range\nproportion") {
  r_land <- terra::mask(r, land_vect)
  plot_df <- as.data.frame(r_land, xy = TRUE) %>%
    setNames(c("x", "y", "value")) %>%
    filter(value > 0 & !is.na(value))

  if (nrow(plot_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = title_txt) +
        annotate("text", x = 0, y = 0, label = "No matched small-range species")
    )
  }

  ggplot() +
    geom_sf(data = world_map, fill = "#f2f2f2", color = NA) +
    geom_tile(data = plot_df, aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_c(
      option = "magma",
      direction = 1,
      name = legend_title,
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    coord_sf(crs = moll_crs, expand = FALSE) +
    labs(title = title_txt) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.height = grid::unit(0.9, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

save_plot <- function(plot_obj, filename, width = 10, height = 6) {
  ggsave(
    filename = file.path(figure_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300,
    bg = "white"
  )
}

# --- 3. MAIN LOOP -------------------------------------------------------------

summary_rows <- list()
species_rows <- list()
count_maps <- list()
proportion_maps <- list()

for (i in seq_len(nrow(taxa_config))) {
  taxon <- taxa_config$taxon[[i]]
  taxon_label <- taxa_config$label[[i]]
  range_file <- here::here(taxa_config$range_file[[i]])

  message("\n--- Processing ", taxon_label, " ---")

  if (!file.exists(range_file)) {
    warning("Skipping ", taxon_label, "; missing range file: ", range_file, call. = FALSE)
    next
  }

  model_data <- read_model_data(taxon)

  range_threshold <- quantile(
    model_data$log_range_size,
    probs = SMALL_RANGE_QUANTILE,
    na.rm = TRUE
  )

  model_data <- model_data %>%
    mutate(
      range_percentile = percent_rank(log_range_size),
      is_small_range = log_range_size <= range_threshold
    )

  model_species <- clean_species_names(model_data$spatial_name)
  small_species <- clean_species_names(model_data$spatial_name[model_data$is_small_range])

  species_rows[[taxon]] <- model_data %>%
    filter(is_small_range) %>%
    select(taxon, phylo_name, spatial_name, log_range_size, range_percentile) %>%
    arrange(log_range_size)

  message("  Model species: ", length(model_species))
  message("  Small-range species: ", length(small_species))
  message("  Loading ranges: ", basename(range_file))

  ranges_sf <- readRDS(range_file)

  if (!inherits(ranges_sf, "sf")) {
    stop("Range file is not an sf object: ", range_file)
  }
  if (!"sci_name" %in% names(ranges_sf)) {
    stop("Range file is missing sci_name column: ", range_file)
  }

  matched_model_species <- intersect(model_species, ranges_sf$sci_name)
  matched_small_species <- intersect(small_species, ranges_sf$sci_name)

  total_raster_path <- file.path(
    raster_dir,
    paste0(taxon, "_model_species_richness.tif")
  )
  small_raster_path <- file.path(
    raster_dir,
    paste0(taxon, "_small_range_", quantile_tag, "_richness.tif")
  )
  proportion_raster_path <- file.path(
    raster_dir,
    paste0(taxon, "_small_range_", quantile_tag, "_proportion.tif")
  )

  total_richness <- calculate_or_load_richness(
    species_names = model_species,
    ranges_sf = ranges_sf,
    output_path = total_raster_path,
    label = paste(taxon_label, "model species richness")
  )

  small_richness <- calculate_or_load_richness(
    species_names = small_species,
    ranges_sf = ranges_sf,
    output_path = small_raster_path,
    label = paste(taxon_label, quantile_tag, "small-range richness")
  )

  small_proportion <- safe_divide(small_richness, total_richness)
  names(small_proportion) <- paste0(taxon, "_small_range_proportion")
  terra::writeRaster(small_proportion, proportion_raster_path, overwrite = TRUE)

  thresholds <- hotspot_thresholds(small_richness)
  max_small_richness <- max(get_positive_land_values(small_richness), na.rm = TRUE)
  if (!is.finite(max_small_richness)) {
    max_small_richness <- 0
  }

  summary_rows[[taxon]] <- data.frame(
    taxon = taxon,
    label = taxon_label,
    small_range_quantile = SMALL_RANGE_QUANTILE,
    log_range_size_threshold = as.numeric(range_threshold),
    n_model_species = length(model_species),
    n_small_range_species = length(small_species),
    n_model_species_matched_to_ranges = length(matched_model_species),
    n_small_range_species_matched_to_ranges = length(matched_small_species),
    n_small_range_species_missing_ranges = length(setdiff(small_species, ranges_sf$sci_name)),
    max_small_range_richness = max_small_richness,
    top_10_hotspot_threshold = thresholds[["top_10"]],
    top_5_hotspot_threshold = thresholds[["top_5"]],
    top_1_hotspot_threshold = thresholds[["top_1"]]
  )

  count_plot <- plot_small_range_hotspots(
    small_richness,
    title_txt = paste0(taxon_label, ": Small-Range Species Hotspots")
  )
  proportion_plot <- plot_proportion_gradient(
    small_proportion,
    title_txt = paste0(taxon_label, ": Small-Range Share of Modelled Richness")
  )

  save_plot(
    count_plot,
    paste0(taxon, "_small_range_", quantile_tag, "_hotspots.png")
  )
  save_plot(
    proportion_plot,
    paste0(taxon, "_small_range_", quantile_tag, "_proportion.png")
  )

  count_maps[[taxon]] <- count_plot
  proportion_maps[[taxon]] <- proportion_plot

  rm(ranges_sf, total_richness, small_richness, small_proportion)
  gc()
}

# --- 4. SUMMARY OUTPUTS -------------------------------------------------------

summary_df <- bind_rows(summary_rows)
species_df <- bind_rows(species_rows)

write.csv(
  summary_df,
  file.path(table_dir, "small_range_hotspot_summary.csv"),
  row.names = FALSE
)

write.csv(
  species_df,
  file.path(table_dir, "small_range_species_list.csv"),
  row.names = FALSE
)

# --- 5. COMBINED PANELS -------------------------------------------------------

if (length(count_maps) > 0) {
  count_panel <- patchwork::wrap_plots(count_maps, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(
        "Hotspots of Small-Range Species (Lowest ",
        round(SMALL_RANGE_QUANTILE * 100),
        "% Within Each Taxon)"
      )
    ) &
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  save_plot(
    count_panel,
    paste0("all_taxa_small_range_", quantile_tag, "_hotspot_panel.png"),
    width = 13,
    height = 15
  )
}

if (length(proportion_maps) > 0) {
  proportion_panel <- patchwork::wrap_plots(proportion_maps, ncol = 2, guides = "collect") +
    patchwork::plot_annotation(
      title = paste0(
        "Local Share of Modelled Richness from Small-Range Species (Lowest ",
        round(SMALL_RANGE_QUANTILE * 100),
        "% Within Each Taxon)"
      )
    ) &
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  save_plot(
    proportion_panel,
    paste0("all_taxa_small_range_", quantile_tag, "_proportion_panel.png"),
    width = 13,
    height = 15
  )
}

message("\nSmall-range hotspot maps complete.")
message("Figures: ", figure_dir)
message("Tables: ", table_dir)
message("Rasters: ", raster_dir)
