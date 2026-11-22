# =============================================================================
# SCRIPT: WOOD PRODUCTION vs SPECIES THREATENED, BY COUNTRY (IUCN LOCATIONS)
#
# Description:
# Companion to 09_country_production.R, which attributes species to countries
# by intersecting dissolved range POLYGONS with Natural Earth country borders.
# This script instead uses the country list IUCN assessors attach directly to
# each assessment (assessment$locations), which avoids all border/geometry
# artefacts (Arabian Peninsula clipping, invalid-geometry repairs, etc.) at the
# cost of relying on IUCN's own country calls, which can lag range updates for
# some species.
#
# A species counts toward a country when it has an IUCN "locations" entry for
# that country with:
#   presence == "Extant"
#   origin   %in% c("Native", "Reintroduced")
# i.e. the standard "currently present, self-sustaining, not introduced/vagrant"
# definition. Change KEEP_PRESENCE / KEEP_ORIGIN below to widen or narrow this.
#
# NOTE (gotcha): assessment$locations$description is itself a data.frame with
# an `en` column -- the country name is at locations$description$en, NOT
# locations$en. Missing this silently yields all-NA country names.
#
# Outputs:
#   outputs/hotspots/tables/country_wood_production_threat_intensity_iucn_locations.csv
#   outputs/hotspots/figures/country_wood_production_scatter_iucn_locations.png
# =============================================================================

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "here",
  "rnaturalearth", "ggplot2", "scales", "ggrepel", "purrr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

# --- CONFIGURATION -----------------------------------------------------------
KEEP_PRESENCE <- c("Extant")
KEEP_ORIGIN   <- c("Native", "Reintroduced")
FAO_YEAR      <- "2024"
MIN_PRODUCTION_M3 <- 1e6
TOP_N_LABEL       <- 15

taxon_config <- tibble::tribble(
  ~taxon_key,  ~cache_file,
  "bird",      here("data", "raw", "iucn", "iucn_assessment_list.rds"),
  "mammal",    here("data", "raw", "iucn", "mam_iucn_assessment_list.rds"),
  "amphibian", here("data", "raw", "iucn", "amphib_iucn_assessment_list.rds"),
  "reptile",   here("data", "raw", "iucn", "rept_iucn_assessment_list.rds"),
  "fish",      here("data", "raw", "iucn", "fish_iucn_assessment_list.rds")
)

out_dir_tab  <- here("outputs", "hotspots", "tables")
out_dir_figs <- here("outputs", "hotspots", "figures")
dir.create(out_dir_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figs, recursive = TRUE, showWarnings = FALSE)

clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  x %>%
    str_replace_all("�", " ") %>%
    str_replace_all(" ", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

# --- 1. EXTRACT PER-SPECIES COUNTRY LISTS FROM RAW IUCN ASSESSMENTS ---------
extract_locations_for_taxon <- function(taxon_key, cache_file) {
  cat(sprintf("--- Extracting locations: %s ---\n", taxon_key))
  assessment_list <- readRDS(cache_file)

  rows <- purrr::map(assessment_list, function(a) {
    loc <- a$locations
    if (!is.data.frame(loc) || nrow(loc) == 0) return(NULL)

    tibble(
      scientific_name = clean_names(a$taxon$scientific_name),
      category        = a$red_list_category$code,
      iso2            = loc$code,
      country_name    = clean_names(loc$description$en),   # nested field -- see header note
      presence        = loc$presence,
      origin          = loc$origin
    )
  })

  out <- bind_rows(rows) %>%
    mutate(taxon_key = taxon_key, .before = scientific_name)

  cat(sprintf("  -> %d species-country rows (%d species)\n",
              nrow(out), n_distinct(out$scientific_name)))
  out
}

all_locations <- purrr::map2_dfr(
  taxon_config$taxon_key, taxon_config$cache_file,
  extract_locations_for_taxon
)

readr::write_csv(
  all_locations,
  here("data", "processed", "assessments", "all_taxa_assessment_locations.csv")
)

locations_filtered <- all_locations %>%
  filter(presence %in% KEEP_PRESENCE, origin %in% KEEP_ORIGIN) %>%
  filter(!is.na(iso2), iso2 != "")

cat(sprintf(
  "\nAfter presence/origin filter: %d species-country rows across %d species, %d countries\n",
  nrow(locations_filtered), n_distinct(locations_filtered$scientific_name),
  n_distinct(locations_filtered$iso2)
))

# --- 2. LOAD PER-TAXON THREAT FLAGS (already computed) -----------------------
load_species_summary <- function(taxon_key) {
  readRDS(here("data", "processed", "assessments",
               paste0(taxon_key, "_assessment_species_summary.rds"))) %>%
    filter(!category %in% c("EX", "EW")) %>%
    transmute(
      taxon_key,
      scientific_name = clean_names(scientific_name),
      category,
      has_wood_prod_ongoing,
      has_logging_ongoing,
      has_plantation_ongoing
    )
}

all_threats <- purrr::map_dfr(taxon_config$taxon_key, load_species_summary)

# --- 3. JOIN THREATS TO COUNTRIES, COUNT PER COUNTRY -------------------------
species_country_threats <- locations_filtered %>%
  distinct(taxon_key, scientific_name, category, iso2, country_name) %>%
  inner_join(all_threats, by = c("taxon_key", "scientific_name", "category"))

cat(sprintf(
  "\nSpecies-country rows with a matched threat record: %d\n",
  nrow(species_country_threats)
))

country_counts <- species_country_threats %>%
  group_by(iso2, country_name) %>%
  summarise(
    total_extant_species     = n_distinct(scientific_name),
    wood_prod_species        = n_distinct(scientific_name[has_wood_prod_ongoing]),
    logging_species          = n_distinct(scientific_name[has_logging_ongoing]),
    plantation_species       = n_distinct(scientific_name[has_plantation_ongoing]),
    .groups = "drop"
  )

# --- 4. ISO2 -> M49 CROSSWALK, + TROPICAL/NON-TROPICAL CLASSIFICATION -------
# Tropical := country centroid falls between the Tropics (|latitude| <= 23.5 deg).
# A centroid-based call is a simplification for countries straddling the
# tropics (e.g. Brazil, Australia, China, the USA) -- it assigns the whole
# country to one class based on where its landmass centres, not a within-
# country split. Good enough for a first-pass colour code; flag if you want a
# stricter area-weighted tropical/non-tropical split per country.
to_code <- function(x) {
  x <- as.character(x)
  x[x %in% c("-99", "-099", "")] <- NA
  suppressWarnings(as.integer(x))
}

manual_m49 <- c("FR" = 250L, "NO" = 578L)
TROPIC_LAT <- 23.5

world_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica") %>%
  sf::st_make_valid()

world_centroid_lat <- world_sf %>%
  sf::st_point_on_surface() %>%
  mutate(centroid_lat = sf::st_coordinates(.)[, 2]) %>%
  sf::st_drop_geometry()

world_attr <- world_centroid_lat %>%
  transmute(
    iso2 = iso_a2,
    m49  = dplyr::coalesce(to_code(iso_n3), to_code(un_a3)),
    country = name_long,
    centroid_lat = centroid_lat,
    tropical = factor(
      ifelse(abs(centroid_lat) <= TROPIC_LAT, "Tropical", "Non-tropical"),
      levels = c("Tropical", "Non-tropical")
    )
  ) %>%
  mutate(m49 = ifelse(is.na(m49) & iso2 %in% names(manual_m49), manual_m49[iso2], m49)) %>%
  filter(!is.na(iso2), iso2 != "-99") %>%
  distinct(iso2, .keep_all = TRUE)

unmatched_iso2 <- setdiff(country_counts$iso2, world_attr$iso2)
if (length(unmatched_iso2) > 0) {
  cat("\n[Audit] IUCN location codes with no Natural Earth ISO2 match (dropped from FAO join):\n")
  print(country_counts %>% filter(iso2 %in% unmatched_iso2) %>% select(iso2, country_name, total_extant_species))
}

# --- 5. LOAD FAO ROUNDWOOD PRODUCTION ----------------------------------------
fao_raw <- readr::read_csv(
  here("data", "raw", "fao", "Country wood production stats.csv"),
  show_col_types = FALSE
)

fao <- fao_raw %>%
  filter(Element == "Production", Item == "Roundwood", Year == FAO_YEAR) %>%
  transmute(
    m49          = as.integer(`Area Code (M49)`),
    fao_area     = Area,
    roundwood_m3 = suppressWarnings(as.numeric(Value)),
    fao_flag     = Flag
  ) %>%
  filter(!is.na(m49)) %>%
  distinct(m49, .keep_all = TRUE)

cat(sprintf("\nFAO %s roundwood production available for %d areas\n", FAO_YEAR, nrow(fao)))

# --- 6. JOIN & DERIVE METRIC --------------------------------------------------
result <- country_counts %>%
  left_join(
    world_attr %>% select(iso2, m49, ne_country = country, centroid_lat, tropical),
    by = "iso2"
  ) %>%
  left_join(fao, by = "m49") %>%
  mutate(
    country_display = dplyr::coalesce(ne_country, country_name),
    roundwood_million_m3 = roundwood_m3 / 1e6,
    species_per_million_m3 = ifelse(
      !is.na(roundwood_m3) & roundwood_m3 > 0,
      wood_prod_species / roundwood_million_m3,
      NA_real_
    ),
    pct_wood_prod_species = ifelse(
      total_extant_species > 0,
      100 * wood_prod_species / total_extant_species,
      NA_real_
    )
  ) %>%
  select(
    iso2, m49, country = country_display, centroid_lat, tropical,
    total_extant_species, wood_prod_species, logging_species, plantation_species,
    pct_wood_prod_species,
    roundwood_m3, roundwood_million_m3, species_per_million_m3, fao_flag
  ) %>%
  arrange(desc(species_per_million_m3))

rank_lookup <- result %>%
  filter(!is.na(species_per_million_m3), roundwood_m3 >= MIN_PRODUCTION_M3) %>%
  mutate(rank_worst = rank(-species_per_million_m3, ties.method = "min")) %>%
  select(iso2, rank_worst)

result <- result %>%
  select(-any_of("rank_worst")) %>%
  left_join(rank_lookup, by = "iso2")

table_path <- file.path(out_dir_tab, "country_wood_production_threat_intensity_iucn_locations.csv")
readr::write_csv(result, table_path)
cat("\n  -> Table written:", table_path, "\n")

# --- 7. CORRELATION / SANITY CHECK -------------------------------------------
fit_data <- result %>% filter(!is.na(roundwood_m3), roundwood_m3 > 0, wood_prod_species > 0)
if (nrow(fit_data) >= 3) {
  ct <- cor.test(log10(fit_data$roundwood_m3), log10(fit_data$wood_prod_species), method = "pearson")
  cat(sprintf(
    "\nlog-log correlation (roundwood_m3 vs wood_prod_species), n=%d: r = %.3f, p = %.3g\n",
    nrow(fit_data), unname(ct$estimate), ct$p.value
  ))
}

pct_fit_data <- result %>% filter(!is.na(roundwood_m3), roundwood_m3 > 0, !is.na(pct_wood_prod_species))
if (nrow(pct_fit_data) >= 3) {
  ct_pct <- cor.test(log10(pct_fit_data$roundwood_m3), pct_fit_data$pct_wood_prod_species, method = "spearman")
  cat(sprintf(
    "Spearman correlation (roundwood_m3 vs %% species threatened), n=%d: rho = %.3f, p = %.3g\n",
    nrow(pct_fit_data), unname(ct_pct$estimate), ct_pct$p.value
  ))
}

# --- 8. SHARED COLOUR SCALE (tropical / non-tropical) -------------------------
# Okabe-Ito colourblind-safe pair.
tropical_colours <- c("Tropical" = "#E69F00", "Non-tropical" = "#0072B2")
scale_colour_tropical <- function() {
  scale_colour_manual(
    name = "Region", values = tropical_colours,
    na.value = "grey60", na.translate = TRUE,
    labels = c("Tropical", "Non-tropical", "Unclassified")
  )
}

# --- 9. SCATTER PLOT 1: PRODUCTION vs ABSOLUTE SPECIES THREATENED -------------
plot_data <- result %>%
  filter(!is.na(roundwood_m3), roundwood_m3 > 0)

label_data <- plot_data %>%
  filter(!is.na(rank_worst) & rank_worst <= TOP_N_LABEL | wood_prod_species >= quantile(plot_data$wood_prod_species, 0.95, na.rm = TRUE))

p1 <- ggplot(plot_data, aes(x = roundwood_m3, y = wood_prod_species)) +
  geom_point(aes(size = species_per_million_m3, colour = tropical), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", linewidth = 0.6, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = label_data, aes(label = country),
    size = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, colour = "grey20"
  ) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_log10() +
  scale_size_continuous(name = "species per\nmillion m³", range = c(1, 8)) +
  scale_colour_tropical() +
  labs(
    title = "Wood-production-threatened vertebrate species vs. national roundwood output",
    subtitle = paste0(
      "IUCN-assessed country presence vs FAO ", FAO_YEAR, " roundwood production (log-log axes). ",
      "Colour = country centroid within/outside the Tropics (|lat| ≤ ", TROPIC_LAT, "°)."
    ),
    x = expression("Roundwood production (m"^3*", log scale)"),
    y = "Wood-production-threatened species (log scale)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30")
  )

plot_path <- file.path(out_dir_figs, "country_wood_production_scatter_iucn_locations.png")
ggsave(plot_path, p1, width = 10, height = 7, dpi = 300, bg = "white")
cat("  -> Plot written:", plot_path, "\n")

# --- 10. SCATTER PLOT 2: PRODUCTION vs PROPORTION OF SPECIES THREATENED -------
pct_plot_data <- result %>%
  filter(!is.na(roundwood_m3), roundwood_m3 > 0, !is.na(pct_wood_prod_species))

MIN_TOTAL_SPECIES <- 20  # guard against noisy ratios from species-poor territories

pct_rank_lookup <- pct_plot_data %>%
  filter(total_extant_species >= MIN_TOTAL_SPECIES, roundwood_m3 >= MIN_PRODUCTION_M3) %>%
  mutate(pct_rank_worst = rank(-pct_wood_prod_species, ties.method = "min")) %>%
  select(iso2, pct_rank_worst)

pct_plot_data <- pct_plot_data %>%
  left_join(pct_rank_lookup, by = "iso2")

pct_label_data <- pct_plot_data %>%
  filter(!is.na(pct_rank_worst) & pct_rank_worst <= TOP_N_LABEL)

p2 <- ggplot(pct_plot_data, aes(x = roundwood_m3, y = pct_wood_prod_species)) +
  geom_point(aes(size = total_extant_species, colour = tropical), alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "grey30", linewidth = 0.6, linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = pct_label_data, aes(label = country),
    size = 3, max.overlaps = Inf, min.segment.length = 0, segment.size = 0.3, colour = "grey20"
  ) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  scale_size_continuous(name = "extant species\nassessed in country", range = c(1, 8)) +
  scale_colour_tropical() +
  labs(
    title = "Proportion of a country's species threatened by wood production",
    subtitle = paste0(
      "wood-production-threatened species ÷ all extant assessed species, by country, vs FAO ", FAO_YEAR,
      " roundwood production. Colour = country centroid within/outside the Tropics (|lat| ≤ ", TROPIC_LAT, "°). ",
      "Labelled: top ", TOP_N_LABEL, " by proportion among countries with ≥ ", MIN_TOTAL_SPECIES,
      " assessed species and ≥ ", format(MIN_PRODUCTION_M3, big.mark = ",", scientific = FALSE), " m³ production."
    ),
    x = expression("Roundwood production (m"^3*", log scale)"),
    y = "% of assessed species threatened by wood production"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 8.5, colour = "grey30")
  )

pct_plot_path <- file.path(out_dir_figs, "country_wood_production_pct_scatter_iucn_locations.png")
ggsave(pct_plot_path, p2, width = 10, height = 7, dpi = 300, bg = "white")
cat("  -> Plot written:", pct_plot_path, "\n")

# --- 9. CONSOLE SUMMARY -------------------------------------------------------
cat("\n======================== TOP", TOP_N_LABEL, "WORST (species per million m3) ========================\n")
result %>%
  filter(!is.na(rank_worst), rank_worst <= TOP_N_LABEL) %>%
  arrange(rank_worst) %>%
  transmute(
    rank = rank_worst, country, wood_prod_species,
    roundwood_million_m3 = round(roundwood_million_m3, 2),
    species_per_million_m3 = round(species_per_million_m3, 2)
  ) %>%
  print(n = TOP_N_LABEL)

cat("\n======================== TOP", TOP_N_LABEL, "WORST (% of species threatened) ========================\n")
pct_plot_data %>%
  filter(!is.na(pct_rank_worst), pct_rank_worst <= TOP_N_LABEL) %>%
  arrange(pct_rank_worst) %>%
  transmute(
    rank = pct_rank_worst, country, tropical, total_extant_species, wood_prod_species,
    pct_wood_prod_species = round(pct_wood_prod_species, 1)
  ) %>%
  print(n = TOP_N_LABEL)

cat("\nDone.\n")
