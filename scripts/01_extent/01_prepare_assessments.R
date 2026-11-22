rm(list = ls())
gc()

# ==============================================================================
# PREPARE CANONICAL IUCN ASSESSMENT PRODUCTS FOR DOWNSTREAM SCRIPTS
# ==============================================================================

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(rredlist)
library(here)

output_dir <- here("data", "processed", "assessments")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

THREAT_CODE_LOGGING <- "^5_3"
THREAT_CODE_PLANTATION <- "^2_2"
WOOD_PRODUCTION_PATTERN <- "^(5_3|2_2)"
FUELWOOD_PATTERN <- "fuelwood|fuel wood|fuel-wood|firewood|fire-wood|fire wood|charcoal"
INDUSTRIAL_ROUNDWOOD_PATTERN <- "timber|sawlogs|lumber|paper|pulp|plank"

clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  x %>%
    str_replace_all("\uFFFD", " ") %>%
    str_replace_all("\u00A0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

extract_synonym_names <- function(x) {
  if (is.null(x)) return(character(0))
  
  if (is.data.frame(x)) {
    if (all(c("genus_name", "species_name") %in% names(x))) {
      return(
        unique(clean_names(paste(x$genus_name, x$species_name)))
      )
    }
    if ("name" %in% names(x)) {
      return(unique(clean_names(x$name)))
    }
  }
  
  if (is.list(x) && !is.data.frame(x)) {
    return(unique(unlist(lapply(x, extract_synonym_names), use.names = FALSE)))
  }
  
  character(0)
}

filter_assessment_systems <- function(df, keep_systems = c("1", "3")) {
  keep_systems <- as.character(keep_systems)
  
  df %>%
    mutate(system_code = as.character(system_code)) %>%
    group_by(scientific_name, category) %>%
    filter(any(system_code %in% keep_systems, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(system_code %in% keep_systems)
}

standardize_threat_extract <- function(df) {
  # Match the older script convention:
  # code.x = Red List category, code.y = threat code
  if ("code.x" %in% names(df)) names(df)[names(df) == "code.x"] <- "category"
  if ("code.y" %in% names(df)) names(df)[names(df) == "code.y"] <- "threat"
  if ("threats" %in% names(df)) names(df)[names(df) == "threats"] <- "description"
  
  if (!"category" %in% names(df)) df$category <- NA_character_
  if (!"threat" %in% names(df)) df$threat <- NA_character_
  if (!"description" %in% names(df)) df$description <- NA_character_
  if (!"timing" %in% names(df)) df$timing <- NA_character_
  if (!"scientific_name" %in% names(df)) df$scientific_name <- NA_character_
  if (!"family_name" %in% names(df)) df$family_name <- NA_character_
  if (!"order_name" %in% names(df)) df$order_name <- NA_character_
  if (!"synonyms" %in% names(df)) df$synonyms <- vector("list", nrow(df))
  
  df %>%
    mutate(
      scientific_name = clean_names(scientific_name),
      family_name = clean_names(family_name),
      order_name = clean_names(order_name),
      description = clean_names(as.character(description)),
      threat = as.character(threat),
      timing = as.character(timing),
      category = as.character(category)
    )
}

standardize_threat_system_extract <- function(df) {
  # When systems are extracted in the same flattened table as threats and
  # documentation, the system code is kept as plain `code`.
  if ("code" %in% names(df)) names(df)[names(df) == "code"] <- "system_code"
  if (!"system_code" %in% names(df)) df$system_code <- NA_character_
  
  standardize_threat_extract(df) %>%
    mutate(system_code = as.character(system_code))
}

standardize_habitat_extract <- function(df) {
  # Match the older script convention:
  # code.x = Red List category, code.y = habitat code
  if ("code.x" %in% names(df)) names(df)[names(df) == "code.x"] <- "category"
  if ("code.y" %in% names(df)) {
    names(df)[names(df) == "code.y"] <- "habitat"
  } else if ("code" %in% names(df)) {
    names(df)[names(df) == "code"] <- "habitat"
  }
  
  if (!"category" %in% names(df)) df$category <- NA_character_
  if (!"habitat" %in% names(df)) df$habitat <- NA_character_
  if (!"suitability" %in% names(df)) df$suitability <- NA_character_
  if (!"scientific_name" %in% names(df)) df$scientific_name <- NA_character_
  if (!"family_name" %in% names(df)) df$family_name <- NA_character_
  if (!"order_name" %in% names(df)) df$order_name <- NA_character_
  
  df %>%
    mutate(
      scientific_name = clean_names(scientific_name),
      family_name = clean_names(family_name),
      order_name = clean_names(order_name),
      habitat = as.character(habitat),
      suitability = as.character(suitability),
      category = as.character(category)
    )
}

standardize_system_extract <- function(df) {
  # Match the older script convention:
  # code.x = Red List category, code.y = system code
  if ("code.x" %in% names(df)) names(df)[names(df) == "code.x"] <- "category"
  
  if ("code.y" %in% names(df)) {
    names(df)[names(df) == "code.y"] <- "system_code"
  } else if ("code" %in% names(df) && !"category" %in% names(df)) {
    names(df)[names(df) == "code"] <- "system_code"
  } else if ("code" %in% names(df) && "category" %in% names(df) && !"system_code" %in% names(df)) {
    names(df)[names(df) == "code"] <- "system_code"
  } else {
    candidate_cols <- c("systems.code", "systems_code", "system_code")
    matched_col <- intersect(candidate_cols, names(df))[1]
    if (!is.na(matched_col)) {
      names(df)[names(df) == matched_col] <- "system_code"
    } else {
      stop(
        "Could not find a system code column in extracted systems table. Columns were: ",
        paste(names(df), collapse = ", ")
      )
    }
  }
  
  if (!"category" %in% names(df)) df$category <- NA_character_
  if (!"system_code" %in% names(df)) df$system_code <- NA_character_
  if (!"scientific_name" %in% names(df)) df$scientific_name <- NA_character_
  
  df %>%
    mutate(
      scientific_name = clean_names(scientific_name),
      category = as.character(category),
      system_code = as.character(system_code)
    )
}

make_habitat_summary <- function(habitat_rows) {
  if (nrow(habitat_rows) == 0) {
    return(tibble(
      scientific_name = character(),
      category = character(),
      Forest = integer(),
      Savanna = integer(),
      Shrubland = integer(),
      Grassland = integer(),
      Wetlands = integer(),
      Rocky = integer(),
      Caves = integer(),
      Desert = integer(),
      Marine = integer(),
      Artificial = integer(),
      Other = integer(),
      Unknown = integer(),
      habitat_breadth = integer()
    ))
  }
  
  habitat_rows %>%
    mutate(
      non_marginal = is.na(suitability) | suitability != "Marginal",
      is_forest     = str_detect(habitat, "^1(\\_|$)")               & non_marginal,
      is_savanna    = str_detect(habitat, "^2(\\_|$)")               & non_marginal,
      is_shrub      = str_detect(habitat, "^3(\\_|$)")               & non_marginal,
      is_grass      = str_detect(habitat, "^4(\\_|$)")               & non_marginal,
      is_wetland    = str_detect(habitat, "^5(\\_|$)")               & non_marginal,
      is_rocky      = str_detect(habitat, "^6(\\_|$)")               & non_marginal,
      is_cave       = str_detect(habitat, "^7(\\_|$)")               & non_marginal,
      is_desert     = str_detect(habitat, "^8(\\_|$)")               & non_marginal,
      is_marine     = str_detect(habitat, "^(9|10|11|12|13)(\\_|$)") & non_marginal,
      is_artificial = str_detect(habitat, "^(14|15|16)(\\_|$)")      & non_marginal,
      is_other      = str_detect(habitat, "^17(\\_|$)")              & non_marginal,
      is_unknown    = str_detect(habitat, "^18(\\_|$)")              & non_marginal
    ) %>%
    group_by(scientific_name, category) %>%
    summarise(
      Forest     = as.integer(any(is_forest, na.rm = TRUE)),
      Savanna    = as.integer(any(is_savanna, na.rm = TRUE)),
      Shrubland  = as.integer(any(is_shrub, na.rm = TRUE)),
      Grassland  = as.integer(any(is_grass, na.rm = TRUE)),
      Wetlands   = as.integer(any(is_wetland, na.rm = TRUE)),
      Rocky      = as.integer(any(is_rocky, na.rm = TRUE)),
      Caves      = as.integer(any(is_cave, na.rm = TRUE)),
      Desert     = as.integer(any(is_desert, na.rm = TRUE)),
      Marine     = as.integer(any(is_marine, na.rm = TRUE)),
      Artificial = as.integer(any(is_artificial, na.rm = TRUE)),
      Other      = as.integer(any(is_other, na.rm = TRUE)),
      Unknown    = as.integer(any(is_unknown, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      habitat_breadth = Forest + Savanna + Shrubland + Grassland +
        Wetlands + Rocky + Caves + Desert + Marine + Artificial + Other + Unknown
    )
}

build_species_summary <- function(threat_rows, taxa_group, habitat_summary) {
  taxonomy_lookup <- threat_rows %>%
    group_by(scientific_name, category, order_name, family_name) %>%
    summarise(
      Synonyms = list(unique(unlist(lapply(synonyms, extract_synonym_names), use.names = FALSE))),
      .groups = "drop"
    )
  
  threat_summary <- threat_rows %>%
    group_by(scientific_name, category, order_name, family_name) %>%
    summarise(
      description = first(na.omit(description)),
      threats = paste(unique(na.omit(threat)), collapse = "; "),
      has_logging_ongoing = any(str_detect(threat, THREAT_CODE_LOGGING) & timing == "Ongoing", na.rm = TRUE),
      has_plantation_ongoing = any(str_detect(threat, THREAT_CODE_PLANTATION) & timing == "Ongoing", na.rm = TRUE),
      has_wood_prod_ongoing = any(str_detect(threat, WOOD_PRODUCTION_PATTERN) & timing == "Ongoing", na.rm = TRUE),
      has_fuelwood_charcoal_ongoing = any(
        !is.na(description) &
          str_detect(str_to_lower(description), FUELWOOD_PATTERN) &
          str_detect(threat, WOOD_PRODUCTION_PATTERN) &
          timing == "Ongoing",
        na.rm = TRUE
      ),
      has_industrial_roundwood_ongoing = any(
        !is.na(description) &
          str_detect(str_to_lower(description), INDUSTRIAL_ROUNDWOOD_PATTERN) &
          str_detect(threat, WOOD_PRODUCTION_PATTERN) &
          timing == "Ongoing",
        na.rm = TRUE
      ),
      has_small_scale_ongoing = any(str_detect(threat, "^(2_2_1|5_3_3)") & timing == "Ongoing", na.rm = TRUE),
      has_large_scale_ongoing = any(str_detect(threat, "^(2_2_2|5_3_4)") & timing == "Ongoing", na.rm = TRUE),
      .groups = "drop"
    )
  
  taxonomy_lookup %>%
    left_join(threat_summary, by = c("scientific_name", "category", "order_name", "family_name")) %>%
    left_join(habitat_summary, by = c("scientific_name", "category")) %>%
    mutate(
      taxa_group = taxa_group,
      description = na_if(description, ""),
      threats = na_if(threats, ""),
      has_logging_ongoing = replace_na(has_logging_ongoing, FALSE),
      has_plantation_ongoing = replace_na(has_plantation_ongoing, FALSE),
      has_wood_prod_ongoing = replace_na(has_wood_prod_ongoing, FALSE),
      has_fuelwood_charcoal_ongoing = replace_na(has_fuelwood_charcoal_ongoing, FALSE),
      has_industrial_roundwood_ongoing = replace_na(has_industrial_roundwood_ongoing, FALSE),
      has_small_scale_ongoing = replace_na(has_small_scale_ongoing, FALSE),
      has_large_scale_ongoing = replace_na(has_large_scale_ongoing, FALSE)
    )
}

make_patch_summary <- function(patch_df, taxa_group) {
  if (is.null(patch_df) || nrow(patch_df) == 0) return(NULL)
  
  patch_df %>%
    mutate(
      scientific_name = clean_names(scientific_name),
      category = as.character(category),
      order_name = clean_names(order_name),
      family_name = clean_names(family_name),
      description = clean_names(description),
      threats = as.character(threats),
      Synonyms = map(Synonyms, ~ unique(clean_names(.x))),
      taxa_group = taxa_group,
      has_logging_ongoing = str_detect(replace_na(threats, ""), THREAT_CODE_LOGGING),
      has_plantation_ongoing = str_detect(replace_na(threats, ""), THREAT_CODE_PLANTATION),
      has_wood_prod_ongoing = str_detect(replace_na(threats, ""), WOOD_PRODUCTION_PATTERN),
      has_fuelwood_charcoal_ongoing = FALSE,
      has_industrial_roundwood_ongoing = FALSE,
      has_small_scale_ongoing = str_detect(replace_na(threats, ""), "^(2_2_1|5_3_3)"),
      has_large_scale_ongoing = str_detect(replace_na(threats, ""), "^(2_2_2|5_3_4)")
    )
}

write_assessment_outputs <- function(prefix, data, synonyms_to_text = FALSE) {
  flatten_list_cell <- function(x) {
    if (is.null(x)) return(NA_character_)
    if (is.data.frame(x)) {
      vals <- unlist(x, recursive = TRUE, use.names = FALSE)
    } else if (is.list(x)) {
      vals <- unlist(x, recursive = TRUE, use.names = FALSE)
    } else {
      vals <- x
    }
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NA_character_)
    paste(as.character(vals), collapse = "; ")
  }
  
  csv_data <- data
  list_cols <- names(csv_data)[vapply(csv_data, is.list, logical(1))]
  if (length(list_cols) > 0) {
    csv_data <- csv_data %>%
      mutate(across(all_of(list_cols), ~ map_chr(.x, flatten_list_cell)))
  }
  
  if (synonyms_to_text && "Synonyms" %in% names(csv_data)) {
    csv_data <- csv_data %>%
      mutate(Synonyms = as.character(Synonyms))
  }
  
  write.csv(csv_data, file = file.path(output_dir, paste0(prefix, ".csv")), row.names = FALSE)
  saveRDS(data, file = file.path(output_dir, paste0(prefix, ".rds")))
}

prepare_taxon_outputs <- function(taxon_key, taxa_group, cache_file, keep_systems, missing_patch = NULL) {
  assessment_list <- readRDS(cache_file)
  keep_systems <- as.character(unlist(keep_systems, use.names = FALSE))
  
  threat_rows <- rl_assessment_extract(
    assessment_list,
    c("taxon", "red_list_category__code", "threats", "documentation", "systems"),
    format = "df",
    flatten = TRUE
  ) %>%
    standardize_threat_system_extract() %>%
    filter(system_code %in% keep_systems) %>%
    transmute(
      scientific_name = clean_names(scientific_name),
      family_name = clean_names(family_name),
      order_name = clean_names(order_name),
      category = as.character(category),
      threat = as.character(threat),
      timing = as.character(timing),
      description = clean_names(description),
      synonyms = synonyms
    )
  
  kept_species <- threat_rows %>%
    distinct(scientific_name, category)
  
  habitat_rows <- rl_assessment_extract(
    assessment_list,
    c("taxon", "red_list_category__code", "habitats"),
    format = "df",
    flatten = TRUE
  ) %>%
    standardize_habitat_extract() %>%
    semi_join(kept_species, by = c("scientific_name", "category")) %>%
    transmute(
      scientific_name = clean_names(scientific_name),
      family_name = clean_names(family_name),
      order_name = clean_names(order_name),
      category = as.character(category),
      habitat = as.character(habitat),
      suitability = as.character(suitability)
    )
  
  habitat_summary <- make_habitat_summary(habitat_rows)
  species_summary <- build_species_summary(threat_rows, taxa_group, habitat_summary)
  
  patch_summary <- make_patch_summary(missing_patch, taxa_group)
  if (!is.null(patch_summary)) {
    species_summary <- bind_rows(species_summary, patch_summary) %>%
      group_by(scientific_name, category) %>%
      arrange(desc(has_wood_prod_ongoing), scientific_name) %>%
      slice(1) %>%
      ungroup()
  }
  
  taxonomy_lookup <- species_summary %>%
    transmute(
      scientific_name,
      category,
      IUCN_Name = scientific_name,
      RedList_Cat = category,
      Family = family_name,
      Order = order_name,
      Threat_Logging = has_logging_ongoing,
      Threat_Plantation = has_plantation_ongoing,
      Threat_Wood_Production = has_wood_prod_ongoing,
      Synonyms
    )
  
  species_summary <- species_summary %>%
    mutate(
      taxon_key = taxon_key,
      .before = taxa_group
    )
  threat_rows <- threat_rows %>%
    mutate(
      taxa_group = taxa_group,
      taxon_key = taxon_key,
      .before = scientific_name
    )
  habitat_summary <- habitat_summary %>%
    mutate(
      taxa_group = taxa_group,
      taxon_key = taxon_key,
      .before = scientific_name
    )
  taxonomy_lookup <- taxonomy_lookup %>%
    mutate(
      taxa_group = taxa_group,
      taxon_key = taxon_key,
      .before = IUCN_Name
    )
  
  write_assessment_outputs(paste0(taxon_key, "_assessment_threat_rows"), threat_rows)
  write_assessment_outputs(paste0(taxon_key, "_assessment_habitat_summary"), habitat_summary)
  write_assessment_outputs(paste0(taxon_key, "_assessment_species_summary"), species_summary, synonyms_to_text = TRUE)
  write_assessment_outputs(paste0(taxon_key, "_assessment_taxonomy_lookup"), taxonomy_lookup, synonyms_to_text = TRUE)
  
  # Compatibility alias retained for scripts or notes that still reference
  # "assessment_consolidated". This now points to the richer species summary.
  write_assessment_outputs(paste0(taxon_key, "_assessment_consolidated"), species_summary, synonyms_to_text = TRUE)
  
  list(
    threat_rows = threat_rows,
    habitat_summary = habitat_summary,
    species_summary = species_summary,
    taxonomy_lookup = taxonomy_lookup
  )
}

bird_missing_patch <- tibble::tibble(
  scientific_name = c("Sterna dougallii", "Prinia fluviatilis", "Mitu mitu"),
  category = c("LC", "LC", "EW"),
  order_name = NA_character_,
  family_name = NA_character_,
  threats = c("5_1_1; 5_4_4", "", "2_1_4; 5_1_1"),
  description = NA_character_,
  Synonyms = list(character(0), character(0), character(0))
)

taxon_config <- tibble::tribble(
  ~taxon_key,    ~taxa_group,         ~cache_file,                                                ~keep_systems,
  "bird",        "Aves",              here("data", "raw", "iucn", "iucn_assessment_list.rds"),   list(c("0", "1")),
  "mammal",      "Mammalia",          here("data", "raw", "iucn", "mam_iucn_assessment_list.rds"), list(c("0", "1")),
  "amphibian",   "Amphibia",          here("data", "raw", "iucn", "amphib_iucn_assessment_list.rds"), list(c("0", "1")),
  "reptile",     "Reptilia",          here("data", "raw", "iucn", "rept_iucn_assessment_list.rds"), list(c("0", "1")),
  "fish",        "Freshwater Fish",   here("data", "raw", "iucn", "fish_iucn_assessment_list.rds"), list("1")
)

prepared <- vector("list", nrow(taxon_config))
names(prepared) <- taxon_config$taxon_key

for (i in seq_len(nrow(taxon_config))) {
  patch <- if (taxon_config$taxon_key[[i]] == "bird") bird_missing_patch else NULL
  prepared[[i]] <- prepare_taxon_outputs(
    taxon_key = taxon_config$taxon_key[[i]],
    taxa_group = taxon_config$taxa_group[[i]],
    cache_file = taxon_config$cache_file[[i]],
    keep_systems = taxon_config$keep_systems[[i]],
    missing_patch = patch
  )
}

all_threat_rows <- bind_rows(purrr::map(prepared, "threat_rows"))
all_habitat_summary <- bind_rows(purrr::map(prepared, "habitat_summary"))
all_species_summary <- bind_rows(purrr::map(prepared, "species_summary"))
all_taxonomy_lookup <- bind_rows(purrr::map(prepared, "taxonomy_lookup"))

write_assessment_outputs("all_taxa_assessment_threat_rows", all_threat_rows)
write_assessment_outputs("all_taxa_assessment_habitat_summary", all_habitat_summary)
write_assessment_outputs("all_taxa_assessment_species_summary", all_species_summary, synonyms_to_text = TRUE)
write_assessment_outputs("all_taxa_assessment_taxonomy_lookup", all_taxonomy_lookup, synonyms_to_text = TRUE)
write_assessment_outputs("all_taxa_assessment_consolidated", all_species_summary, synonyms_to_text = TRUE)

message("Canonical assessment tables written to ", output_dir)
