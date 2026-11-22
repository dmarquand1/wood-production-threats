# =============================================================================
# SCRIPT: SPECIES THREATENED BY WOOD PRODUCTION PER UNIT FAO ROUNDWOOD OUTPUT
#
# Description:
# Summarises, for every country, the number of vertebrate species threatened by
# wood production (ongoing logging and/or plantation threats) whose IUCN range
# intersects that country, then joins these counts to FAO roundwood production
# (2024, m3) to derive a "threat intensity" metric:
#
#       species_per_million_m3 = n_species_wood_threatened / (roundwood_m3 / 1e6)
#
# Outputs:
#   1. A country-level table (CSV) with per-taxon counts, total threatened
#      species, FAO roundwood production, and the threat-per-production metric.
#   2. A global choropleth (Mollweide) of the metric, highlighting the worst
#      countries.
#
# Notes / caveats (see also the INTERPRETATION block near the bottom):
#   * A species "occurs" in a country if its dissolved range polygon intersects
#     the country polygon (presence, not area-weighted). This mirrors the
#     presence logic used in the richness rasterisation pipeline. Marginal range
#     clipping at borders will therefore count a species in a country; this is
#     the same spatial-misattribution caveat that affects the IUCN range data
#     generally (cf. the Arabian Peninsula artefact). An optional area-overlap
#     threshold is sketched in count_threatened_species_per_country() but left
#     off by default for speed.
#   * The metric is sensitive to the denominator: a high ratio can reflect EITHER
#     a genuinely biodiversity-costly forestry sector OR simply a small producer
#     sitting in a species-rich region. A minimum-production filter is applied
#     for the "worst country" ranking/highlight so micro-producers do not
#     dominate the headline; the full table retains every country.
# =============================================================================

# --- 1. SETUP & LIBRARIES ----------------------------------------------------
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "here",
  "sf", "rnaturalearth", "ggplot2", "scales", "viridisLite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}
invisible(lapply(required_packages, library, character.only = TRUE))

# Planar geometry: consistent with the rest of the hotspot pipeline and robust
# to the (occasionally invalid) IUCN range polygons.
sf::sf_use_s2(FALSE)

# --- 2. CONFIGURATION --------------------------------------------------------
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")

# Inputs (mirrors 02_make_taxon_maps.R)
range_inputs <- list(
  bird      = "data/processed/ranges/all_bird_ranges_dissolved.rds",
  mammal    = "data/processed/ranges/all_mammal_ranges_dissolved.rds",
  amphibian = "data/processed/ranges/all_amphibian_ranges_dissolved.rds",
  reptile   = "data/processed/ranges/all_reptile_ranges_dissolved.rds",
  fish      = "data/processed/ranges/all_fish_ranges_dissolved.rds"
)

# FAO roundwood production CSV (M49 area codes, m3, 2024).
# Adjust if you store the file elsewhere.
fao_csv <- here::here("data", "raw", "fao", "FAOSTAT_data_en_6-12-2026.csv")
if (!file.exists(fao_csv)) {
  # Fall back to a couple of common locations / the working directory.
  candidates <- c(
    here::here("FAOSTAT_data_en_6-12-2026.csv"),
    here::here("data", "raw", "FAOSTAT_data_en_6-12-2026.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) fao_csv <- hit[1]
}

# Ranking / map parameters
MIN_PRODUCTION_M3 <- 1e6   # ignore micro-producers when ranking "worst" countries
TOP_N_WORST       <- 15    # number of worst countries to highlight / label

# Optional stricter spatial attribution (slower). If > 0, a species is only
# counted in a country when the intersected area exceeds this many km2.
MIN_OVERLAP_KM2   <- 0

# Outputs
out_dir_tab  <- here::here("outputs", "hotspots", "tables")
out_dir_maps <- here::here("outputs", "hotspots", "figures", "maps")
dir.create(out_dir_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_maps, recursive = TRUE, showWarnings = FALSE)

load_assessment_species_summary <- function(taxon_key) {
  readRDS(
    here::here(
      "data", "processed", "assessments",
      paste0(taxon_key, "_assessment_species_summary.rds")
    )
  )
}

# --- 3. LOAD & CLEAN FAO ROUNDWOOD DATA --------------------------------------
cat("Loading FAO roundwood production...\n")

fao_raw <- readr::read_csv(fao_csv, show_col_types = FALSE)

fao <- fao_raw %>%
  filter(Element == "Production", Item == "Roundwood") %>%
  transmute(
    m49             = as.integer(`Area Code (M49)`),
    fao_area        = Area,
    roundwood_m3    = suppressWarnings(as.numeric(Value)),
    fao_flag        = Flag
  ) %>%
  filter(!is.na(m49)) %>%
  # FAO carries some sub/super-national aggregates that share names with members
  # (e.g. "China" 159 vs "China, mainland" 156). These aggregate M49 codes have
  # no Natural Earth polygon and so are dropped automatically by the M49 join.
  distinct(m49, .keep_all = TRUE)

cat("  -> FAO countries with a production figure:", nrow(fao), "\n")

# --- 4. COUNTRY POLYGONS + M49 JOIN KEY --------------------------------------
cat("Preparing country polygons...\n")

# Robust numeric country key: prefer ISO numeric (== M49 for almost all
# countries), fall back to the UN M49 string, with a few manual overrides for
# Natural Earth polygons that carry a "-99" code but have real FAO production.
to_code <- function(x) {
  x <- as.character(x)
  x[x %in% c("-99", "-099", "")] <- NA
  suppressWarnings(as.integer(x))
}

manual_m49 <- c(
  "France"  = 250L,
  "Norway"  = 578L
)

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  mutate(
    m49 = dplyr::coalesce(to_code(iso_n3), to_code(un_a3)),
    m49 = ifelse(is.na(m49) & admin %in% names(manual_m49),
                 manual_m49[admin], m49)
  ) %>%
  st_make_valid()

# One polygon per country code (dissolve any multi-row admin units).
world_join <- world %>%
  filter(!is.na(m49)) %>%
  group_by(m49) %>%
  summarise(
    country  = dplyr::first(name_long),
    iso_a3   = dplyr::first(iso_a3),
    .groups  = "drop"
  )

world_moll <- world_join %>%
  st_transform(moll_crs) %>%
  st_make_valid() %>%
  mutate(ctry_id = dplyr::row_number())

cat("  -> Country polygons with a usable M49 code:", nrow(world_moll), "\n")

# --- 5. COUNT THREATENED SPECIES PER COUNTRY (PER TAXON) ----------------------
# Returns a data frame: ctry_id, n  (number of distinct wood-production
# threatened species of this taxon whose range intersects the country).
count_threatened_species_per_country <- function(taxon, countries_moll) {
  cat("\nProcessing:", toupper(taxon), "\n")
  
  assessment_df <- load_assessment_species_summary(taxon) %>%
    filter(!category %in% c("EX", "EW"), has_wood_prod_ongoing)
  
  threatened_names <- unique(assessment_df$scientific_name)
  cat("  -> Wood-production threatened species:", length(threatened_names), "\n")
  if (length(threatened_names) == 0) {
    return(data.frame(ctry_id = countries_moll$ctry_id, n = 0L))
  }
  
  ranges <- readRDS(here::here(range_inputs[[taxon]]))
  ranges_sub <- ranges %>%
    filter(sci_name %in% threatened_names) %>%
    st_transform(moll_crs)
  
  # Drop empty / invalid geometries (validate only where needed -- cheaper than
  # blanket st_make_valid across tens of thousands of features).
  ranges_sub <- ranges_sub[!st_is_empty(ranges_sub), ]
  bad <- !st_is_valid(ranges_sub)
  if (any(bad, na.rm = TRUE)) {
    cat("  -> Repairing", sum(bad, na.rm = TRUE), "invalid range geometries...\n")
    ranges_sub[bad, ] <- st_make_valid(ranges_sub[bad, ])
    ranges_sub <- ranges_sub[!st_is_empty(ranges_sub), ]
  }
  cat("  -> Range features to intersect:", nrow(ranges_sub), "\n")
  
  sci_vec <- ranges_sub$sci_name
  
  if (MIN_OVERLAP_KM2 > 0) {
    # Stricter, area-weighted attribution (slower): count a species only where
    # the country x range intersection exceeds the threshold.
    inter <- suppressWarnings(st_intersection(
      countries_moll[, "ctry_id"], ranges_sub[, "sci_name"]
    ))
    if (nrow(inter) == 0) {
      counts <- data.frame(ctry_id = countries_moll$ctry_id, n = 0L)
    } else {
      inter$area_km2 <- as.numeric(st_area(inter)) / 1e6
      counts <- inter %>%
        st_drop_geometry() %>%
        filter(area_km2 >= MIN_OVERLAP_KM2) %>%
        distinct(ctry_id, sci_name) %>%
        count(ctry_id, name = "n")
    }
  } else {
    # Default: presence by intersection (fast, index-backed).
    hits <- st_intersects(countries_moll, ranges_sub)
    n_vec <- vapply(hits, function(idx) length(unique(sci_vec[idx])), integer(1))
    counts <- data.frame(ctry_id = countries_moll$ctry_id, n = n_vec)
  }
  
  # Ensure every country is present (fill non-intersected with 0).
  out <- data.frame(ctry_id = countries_moll$ctry_id) %>%
    left_join(counts, by = "ctry_id") %>%
    mutate(n = tidyr::replace_na(n, 0L))
  
  rm(ranges, ranges_sub); invisible(gc(verbose = FALSE))
  out
}

# Loop taxa and accumulate a wide country x taxon table.
taxon_counts <- data.frame(ctry_id = world_moll$ctry_id)
for (taxon in taxa_list) {
  ct <- count_threatened_species_per_country(taxon, world_moll)
  taxon_counts[[taxon]] <- ct$n[match(taxon_counts$ctry_id, ct$ctry_id)]
}

# Species sets are disjoint across taxa, so the all-vertebrate total is the row
# sum (no cross-taxon double counting).
taxon_counts <- taxon_counts %>%
  mutate(total_species = rowSums(across(all_of(taxa_list)), na.rm = TRUE))

# --- 6. JOIN TO FAO & DERIVE THE METRIC --------------------------------------
cat("\nJoining species counts to FAO production...\n")

country_meta <- world_moll %>%
  st_drop_geometry() %>%
  select(ctry_id, m49, country, iso_a3)

result <- country_meta %>%
  left_join(taxon_counts, by = "ctry_id") %>%
  left_join(fao %>% select(m49, fao_area, roundwood_m3, fao_flag), by = "m49") %>%
  mutate(
    roundwood_million_m3   = roundwood_m3 / 1e6,
    species_per_million_m3 = ifelse(
      !is.na(roundwood_m3) & roundwood_m3 > 0,
      total_species / roundwood_million_m3,
      NA_real_
    )
  ) %>%
  arrange(desc(species_per_million_m3))

# Rank (worst = highest threatened species per unit production), restricted to
# meaningful producers so micro-producers do not dominate.
result <- result %>%
  mutate(
    rank_worst = ifelse(
      !is.na(species_per_million_m3) & roundwood_m3 >= MIN_PRODUCTION_M3,
      rank(-species_per_million_m3, ties.method = "min"),
      NA_integer_
    )
  )

# --- 7. WRITE TABLE ----------------------------------------------------------
table_cols <- c(
  "country", "iso_a3", "m49", taxa_list, "total_species",
  "roundwood_m3", "roundwood_million_m3", "species_per_million_m3",
  "fao_flag", "rank_worst"
)
out_table <- result %>% select(all_of(table_cols))

table_path <- file.path(out_dir_tab, "country_wood_production_threat_intensity.csv")
readr::write_csv(out_table, table_path)
cat("  -> Table written:", table_path, "\n")

# --- 8. AUDIT: UNMATCHED RECORDS ---------------------------------------------
fao_unmatched <- fao %>% filter(!m49 %in% world_moll$m49)
poly_no_fao   <- result %>% filter(is.na(roundwood_m3))

if (nrow(fao_unmatched) > 0) {
  cat("\n[Audit] FAO areas with no Natural Earth polygon (dropped):\n")
  print(fao_unmatched %>% select(m49, fao_area, roundwood_m3), n = Inf)
}
if (nrow(poly_no_fao) > 0) {
  cat("\n[Audit] Country polygons with no FAO production figure:\n")
  print(poly_no_fao %>% select(m49, country) %>% arrange(country), n = Inf)
}

# --- 9. MAP: THREAT INTENSITY, WORST COUNTRIES HIGHLIGHTED --------------------
cat("\nBuilding map...\n")

map_df <- world_join %>%
  st_transform(moll_crs) %>%
  left_join(
    result %>% select(m49, species_per_million_m3, total_species,
                      roundwood_m3, rank_worst),
    by = "m49"
  )

worst <- map_df %>% filter(!is.na(rank_worst) & rank_worst <= TOP_N_WORST)

# Centroids for labelling the worst countries.
worst_pts <- worst %>%
  st_point_on_surface() %>%
  mutate(
    lab = paste0(country, " (#", rank_worst, ")")
  )
worst_xy <- worst_pts %>%
  mutate(
    x = sf::st_coordinates(.)[, 1],
    y = sf::st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry()

p <- ggplot() +
  geom_sf(data = map_df, aes(fill = species_per_million_m3),
          color = "white", linewidth = 0.08) +
  # Highlight the worst producers with a heavy outline.
  geom_sf(data = worst, fill = NA, color = "#111111", linewidth = 0.5) +
  ggrepel::geom_text_repel(
    data = worst_xy, aes(x = x, y = y, label = lab),
    size = 3, fontface = "bold", color = "#111111",
    box.padding = 0.4, max.overlaps = Inf, min.segment.length = 0,
    segment.color = "#555555", segment.size = 0.3
  ) +
  scale_fill_gradientn(
    colours  = viridisLite::magma(256, direction = -1),
    trans    = "log10",
    na.value = "grey92",
    name     = "Wood-production\nthreatened species\nper million m\u00b3\nroundwood",
    labels   = scales::label_number(accuracy = 0.1)
  ) +
  labs(
    title    = "Vertebrate species threatened by wood production per unit roundwood output",
    subtitle = paste0(
      "Threatened species (IUCN range intersect) \u00f7 FAO 2024 roundwood (m\u00b3). ",
      "Top ", TOP_N_WORST, " producers (\u2265 ",
      format(MIN_PRODUCTION_M3, big.mark = ",", scientific = FALSE),
      " m\u00b3) highlighted."
    )
  ) +
  coord_sf(crs = moll_crs, expand = FALSE) +
  theme_void(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey30")
  )

# ggrepel is optional; fall back to no labels if it is unavailable.
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  cat("  [Note] ggrepel not installed -- producing map without country labels.\n")
  p <- ggplot() +
    geom_sf(data = map_df, aes(fill = species_per_million_m3),
            color = "white", linewidth = 0.08) +
    geom_sf(data = worst, fill = NA, color = "#111111", linewidth = 0.5) +
    scale_fill_gradientn(
      colours  = viridisLite::magma(256, direction = -1),
      trans    = "log10", na.value = "grey92",
      name     = "Wood-production\nthreatened species\nper million m\u00b3\nroundwood",
      labels   = scales::label_number(accuracy = 0.1)
    ) +
    labs(title = "Vertebrate species threatened by wood production per unit roundwood output") +
    coord_sf(crs = moll_crs, expand = FALSE) +
    theme_void(base_size = 13) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 14))
}

map_path <- file.path(out_dir_maps, "country_wood_production_threat_intensity.png")
ggsave(map_path, p, width = 12, height = 6.5, dpi = 300, bg = "white")
cat("  -> Map written:", map_path, "\n")

# --- 10. CONSOLE SUMMARY -----------------------------------------------------
cat("\n========================  TOP", TOP_N_WORST,
    "WORST PRODUCERS  ========================\n")
cat("(threatened species per million m3 roundwood; producers >= ",
    format(MIN_PRODUCTION_M3, big.mark = ",", scientific = FALSE), " m3)\n\n", sep = "")

result %>%
  filter(!is.na(rank_worst), rank_worst <= TOP_N_WORST) %>%
  arrange(rank_worst) %>%
  transmute(
    rank = rank_worst,
    country,
    total_species,
    roundwood_million_m3 = round(roundwood_million_m3, 2),
    species_per_million_m3 = round(species_per_million_m3, 2)
  ) %>%
  print(n = TOP_N_WORST)

cat("\nDone.\n")

# =============================================================================
# INTERPRETATION
# -----------------------------------------------------------------------------
# This metric ranks countries by how many wood-production-threatened vertebrates
# are exposed PER unit of roundwood actually produced. A high value flags places
# where the forestry sector co-occurs with concentrated, vulnerable biodiversity
# relative to its output. Read it alongside two complementary views to avoid
# mis-reading the denominator:
#   * total_species (absolute exposure) -- large biodiverse producers (e.g.
#     Brazil, Indonesia, DRC) may rank low on the ratio yet carry the greatest
#     absolute number of threatened species.
#   * roundwood_m3 (the production base) -- small producers in the tropics can
#     post very high ratios simply because the denominator is tiny; the
#     MIN_PRODUCTION_M3 filter guards the headline ranking against this.
# For thesis use it is worth reporting the ratio map next to the absolute
# threatened-richness map so the two framings are visible together.
# =============================================================================