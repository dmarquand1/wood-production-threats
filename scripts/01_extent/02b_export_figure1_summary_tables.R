.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = FALSE), .libPaths()))

library(dplyr)
library(here)
library(tidyr)

output_dir <- here::here("outputs", "extent", "tables", "figure1_summary")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

iucn_levels <- c("DD", "LC", "NT", "VU", "EN", "CR")
taxa_levels <- c("Amphibia", "Aves", "Freshwater Fish", "Mammalia", "Reptilia")
summary_taxa_levels <- c("All Vertebrates", taxa_levels)
extinction_threat_levels <- c("VU", "EN", "CR")
figure1_status_levels <- c(
  "Threatened with Wood Production Threat",
  "Not threatened with Wood Production threat",
  "Data Deficient with Wood Production Threat",
  "Species without Wood Production Threat"
)
production_type_levels <- c("Plantations", "Both", "Logging")

round_pct <- function(x) {
  round(100 * x, 2)
}

safe_pct <- function(numerator, denominator) {
  denominator <- rep(denominator, length.out = length(numerator))
  ifelse(denominator > 0, round_pct(numerator / denominator), NA_real_)
}

make_subset_row <- function(data, label, keep_rows, all_extant_n) {
  subset_data <- data[keep_rows, , drop = FALSE]
  species_count <- nrow(subset_data)
  extinction_threatened_species <- sum(subset_data$is_extinction_threatened, na.rm = TRUE)
  
  tibble(
    summary_group = label,
    species_count = species_count,
    extinction_threatened_species = extinction_threatened_species,
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    ),
    denominator_label = "All extant vertebrates",
    denominator_species = all_extant_n,
    pct_of_denominator = safe_pct(species_count, all_extant_n)
  )
}

assessments <- readRDS(
  here::here("data", "processed", "assessments", "all_taxa_assessment_species_summary.rds")
) %>%
  distinct(scientific_name, .keep_all = TRUE) %>%
  mutate(
    category = factor(category, levels = iucn_levels, ordered = TRUE),
    taxa_group = factor(taxa_group, levels = taxa_levels)
  ) %>%
  filter(!category %in% c("EX", "EW"), !is.na(category)) %>%
  mutate(
    is_extinction_threatened = category %in% extinction_threat_levels,
    has_wood_prod = has_wood_prod_ongoing,
    has_logging = has_logging_ongoing,
    has_plantation = has_plantation_ongoing,
    figure1_status = case_when(
      has_wood_prod & is_extinction_threatened ~ figure1_status_levels[1],
      has_wood_prod & !is_extinction_threatened & category %in% c("LC", "NT") ~ figure1_status_levels[2],
      has_wood_prod & category == "DD" ~ figure1_status_levels[3],
      !has_wood_prod ~ figure1_status_levels[4],
      TRUE ~ NA_character_
    ),
    production_type = case_when(
      has_wood_prod & has_logging & has_plantation ~ "Both",
      has_wood_prod & has_logging & !has_plantation ~ "Logging",
      has_wood_prod & !has_logging & has_plantation ~ "Plantations",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    figure1_status = factor(figure1_status, levels = figure1_status_levels),
    production_type = factor(production_type, levels = production_type_levels)
  )

all_extant_n <- nrow(assessments)

assessments_with_combined <- bind_rows(
  assessments %>% mutate(taxa_group = as.character(taxa_group)),
  assessments %>% mutate(taxa_group = "All Vertebrates")
) %>%
  mutate(taxa_group = factor(taxa_group, levels = summary_taxa_levels))

overview_summary <- bind_rows(
  make_subset_row(assessments, "All extant vertebrates", rep(TRUE, all_extant_n), all_extant_n),
  make_subset_row(assessments, "Species with ongoing wood production threat", assessments$has_wood_prod, all_extant_n),
  make_subset_row(assessments, "Species without ongoing wood production threat", !assessments$has_wood_prod, all_extant_n),
  make_subset_row(assessments, "Species with ongoing logging threat", assessments$has_logging, all_extant_n),
  make_subset_row(assessments, "Species with ongoing plantation threat", assessments$has_plantation, all_extant_n),
  make_subset_row(
    assessments,
    "Species with both logging and plantation threats",
    assessments$has_logging & assessments$has_plantation,
    all_extant_n
  ),
  make_subset_row(
    assessments,
    "Species with logging-only threat",
    assessments$has_logging & !assessments$has_plantation,
    all_extant_n
  ),
  make_subset_row(
    assessments,
    "Species with plantation-only threat",
    !assessments$has_logging & assessments$has_plantation,
    all_extant_n
  )
)

taxon_overview <- assessments_with_combined %>%
  group_by(taxa_group) %>%
  summarise(
    total_species = n(),
    extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
    pct_extinction_threatened_overall = safe_pct(
      extinction_threatened_species,
      total_species
    ),
    wood_prod_species = sum(has_wood_prod, na.rm = TRUE),
    wood_prod_extinction_threatened_species = sum(
      has_wood_prod & is_extinction_threatened,
      na.rm = TRUE
    ),
    pct_extinction_threatened_within_wood_prod = safe_pct(
      wood_prod_extinction_threatened_species,
      wood_prod_species
    ),
    logging_species = sum(has_logging, na.rm = TRUE),
    logging_extinction_threatened_species = sum(
      has_logging & is_extinction_threatened,
      na.rm = TRUE
    ),
    pct_extinction_threatened_within_logging = safe_pct(
      logging_extinction_threatened_species,
      logging_species
    ),
    plantation_species = sum(has_plantation, na.rm = TRUE),
    plantation_extinction_threatened_species = sum(
      has_plantation & is_extinction_threatened,
      na.rm = TRUE
    ),
    pct_extinction_threatened_within_plantation = safe_pct(
      plantation_extinction_threatened_species,
      plantation_species
    ),
    .groups = "drop"
  )

panel_a_status <- assessments %>%
  count(figure1_status, name = "species_count") %>%
  complete(
    figure1_status = factor(figure1_status_levels, levels = figure1_status_levels),
    fill = list(species_count = 0)
  ) %>%
  left_join(
    assessments %>%
      group_by(figure1_status) %>%
      summarise(
        extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "figure1_status"
  ) %>%
  mutate(
    extinction_threatened_species = replace_na(extinction_threatened_species, 0L),
    denominator_label = "All extant vertebrates",
    denominator_species = all_extant_n,
    pct_of_denominator = safe_pct(species_count, denominator_species),
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    )
  )

wood_prod_species <- assessments %>%
  filter(has_wood_prod)

panel_a_production_type <- wood_prod_species %>%
  count(production_type, name = "species_count") %>%
  complete(
    production_type = factor(production_type_levels, levels = production_type_levels),
    fill = list(species_count = 0)
  ) %>%
  left_join(
    wood_prod_species %>%
      group_by(production_type) %>%
      summarise(
        extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "production_type"
  ) %>%
  mutate(
    extinction_threatened_species = replace_na(extinction_threatened_species, 0L),
    denominator_label = "Species with ongoing wood production threat",
    denominator_species = nrow(wood_prod_species),
    pct_of_denominator = safe_pct(species_count, denominator_species),
    pct_of_all_extant_species = safe_pct(species_count, all_extant_n),
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    )
  )

panel_b_status_by_taxon <- assessments_with_combined %>%
  count(taxa_group, figure1_status, name = "species_count") %>%
  complete(
    taxa_group = factor(summary_taxa_levels, levels = summary_taxa_levels),
    figure1_status = factor(figure1_status_levels, levels = figure1_status_levels),
    fill = list(species_count = 0)
  ) %>%
  left_join(
    assessments_with_combined %>%
      group_by(taxa_group, figure1_status) %>%
      summarise(
        extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("taxa_group", "figure1_status")
  ) %>%
  group_by(taxa_group) %>%
  mutate(
    extinction_threatened_species = replace_na(extinction_threatened_species, 0L),
    denominator_label = "All extant species within taxon",
    denominator_species = sum(species_count),
    pct_of_denominator = safe_pct(species_count, denominator_species),
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    )
  ) %>%
  ungroup()

wood_prod_species_with_combined <- assessments_with_combined %>%
  filter(has_wood_prod)

panel_b_production_by_taxon <- wood_prod_species_with_combined %>%
  count(taxa_group, production_type, name = "species_count") %>%
  complete(
    taxa_group = factor(summary_taxa_levels, levels = summary_taxa_levels),
    production_type = factor(production_type_levels, levels = production_type_levels),
    fill = list(species_count = 0)
  ) %>%
  left_join(
    wood_prod_species_with_combined %>%
      group_by(taxa_group, production_type) %>%
      summarise(
        extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("taxa_group", "production_type")
  ) %>%
  left_join(
    assessments_with_combined %>%
      count(taxa_group, name = "all_taxon_species"),
    by = "taxa_group"
  ) %>%
  group_by(taxa_group) %>%
  mutate(
    extinction_threatened_species = replace_na(extinction_threatened_species, 0L),
    denominator_label = "Wood-production-threatened species within taxon",
    denominator_species = sum(species_count),
    pct_of_denominator = safe_pct(species_count, denominator_species),
    pct_of_all_taxon_species = safe_pct(species_count, all_taxon_species),
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    )
  ) %>%
  ungroup()

panel_c_taxon_iucn_segments <- assessments_with_combined %>%
  count(taxa_group, category, figure1_status, name = "species_count") %>%
  complete(
    taxa_group = factor(summary_taxa_levels, levels = summary_taxa_levels),
    category = factor(iucn_levels, levels = iucn_levels, ordered = TRUE),
    figure1_status = factor(figure1_status_levels, levels = figure1_status_levels),
    fill = list(species_count = 0)
  ) %>%
  left_join(
    assessments_with_combined %>%
      group_by(taxa_group, category, figure1_status) %>%
      summarise(
        extinction_threatened_species = sum(is_extinction_threatened, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("taxa_group", "category", "figure1_status")
  ) %>%
  group_by(taxa_group, category) %>%
  mutate(
    extinction_threatened_species = replace_na(extinction_threatened_species, 0L),
    denominator_label = "Species within taxon and IUCN category",
    denominator_species = sum(species_count),
    pct_of_denominator = safe_pct(species_count, denominator_species),
    pct_extinction_threatened_within_row = safe_pct(
      extinction_threatened_species,
      species_count
    )
  ) %>%
  ungroup() %>%
  mutate(iucn_category_is_extinction_threatened = category %in% extinction_threat_levels)

master_summary <- bind_rows(
  overview_summary %>%
    transmute(
      panel = "Overview",
      summary_group,
      taxa_group = NA_character_,
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_a_status %>%
    transmute(
      panel = "Figure 1A - Overall status",
      summary_group = as.character(figure1_status),
      taxa_group = NA_character_,
      iucn_category = NA_character_,
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_a_production_type %>%
    transmute(
      panel = "Figure 1A - Production type",
      summary_group = as.character(production_type),
      taxa_group = NA_character_,
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = as.character(production_type),
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_b_status_by_taxon %>%
    transmute(
      panel = "Figure 1B - Status by taxon",
      summary_group = paste(taxa_group, figure1_status, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = NA_character_,
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_b_production_by_taxon %>%
    transmute(
      panel = "Figure 1B - Production type by taxon",
      summary_group = paste(taxa_group, production_type, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = as.character(production_type),
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_c_taxon_iucn_segments %>%
    transmute(
      panel = "Figure 1C - Taxon x IUCN x status segments",
      summary_group = paste(taxa_group, category, figure1_status, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = as.character(category),
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    )
)

all_vertebrates_only <- bind_rows(
  overview_summary %>%
    transmute(
      panel = "Overview",
      summary_group,
      taxa_group = "All Vertebrates",
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  taxon_overview %>%
    filter(taxa_group == "All Vertebrates") %>%
    transmute(
      panel = "Taxon summary",
      summary_group = "All Vertebrates",
      taxa_group = as.character(taxa_group),
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = NA_character_,
      species_count = total_species,
      denominator_label = "All extant vertebrates",
      denominator_species = total_species,
      pct_of_denominator = 100,
      extinction_threatened_species,
      pct_extinction_threatened_within_row = pct_extinction_threatened_overall
    ),
  panel_a_status %>%
    transmute(
      panel = "Figure 1A - Overall status",
      summary_group = as.character(figure1_status),
      taxa_group = "All Vertebrates",
      iucn_category = NA_character_,
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_a_production_type %>%
    transmute(
      panel = "Figure 1A - Production type",
      summary_group = as.character(production_type),
      taxa_group = "All Vertebrates",
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = as.character(production_type),
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_b_status_by_taxon %>%
    filter(taxa_group == "All Vertebrates") %>%
    transmute(
      panel = "Figure 1B - Status by taxon",
      summary_group = paste(taxa_group, figure1_status, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = NA_character_,
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_b_production_by_taxon %>%
    filter(taxa_group == "All Vertebrates") %>%
    transmute(
      panel = "Figure 1B - Production type by taxon",
      summary_group = paste(taxa_group, production_type, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = NA_character_,
      figure1_status = NA_character_,
      production_type = as.character(production_type),
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    ),
  panel_c_taxon_iucn_segments %>%
    filter(taxa_group == "All Vertebrates") %>%
    transmute(
      panel = "Figure 1C - Taxon x IUCN x status segments",
      summary_group = paste(taxa_group, category, figure1_status, sep = " | "),
      taxa_group = as.character(taxa_group),
      iucn_category = as.character(category),
      figure1_status = as.character(figure1_status),
      production_type = NA_character_,
      species_count,
      denominator_label,
      denominator_species,
      pct_of_denominator,
      extinction_threatened_species,
      pct_extinction_threatened_within_row
    )
)

write.csv(
  overview_summary,
  file = file.path(output_dir, "figure1_summary_overview.csv"),
  row.names = FALSE
)
write.csv(
  taxon_overview,
  file = file.path(output_dir, "figure1_summary_by_taxon.csv"),
  row.names = FALSE
)
write.csv(
  panel_a_status,
  file = file.path(output_dir, "figure1_panel_a_status.csv"),
  row.names = FALSE
)
write.csv(
  panel_a_production_type,
  file = file.path(output_dir, "figure1_panel_a_production_type.csv"),
  row.names = FALSE
)
write.csv(
  panel_b_status_by_taxon,
  file = file.path(output_dir, "figure1_panel_b_status_by_taxon.csv"),
  row.names = FALSE
)
write.csv(
  panel_b_production_by_taxon,
  file = file.path(output_dir, "figure1_panel_b_production_type_by_taxon.csv"),
  row.names = FALSE
)
write.csv(
  panel_c_taxon_iucn_segments,
  file = file.path(output_dir, "figure1_panel_c_taxon_iucn_segments.csv"),
  row.names = FALSE
)
write.csv(
  master_summary,
  file = file.path(output_dir, "figure1_summary_all_categories.csv"),
  row.names = FALSE
)
write.csv(
  all_vertebrates_only,
  file = file.path(output_dir, "figure1_all_vertebrates_only.csv"),
  row.names = FALSE
)

message("Figure 1 summary tables written to: ", output_dir)
