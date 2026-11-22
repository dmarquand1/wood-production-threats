rm(list = ls())
gc()

# =============================================================================
# SUMMARISE RED LIST THREAT-ASSESSMENT COMPLETENESS
#
# Creates summary tables for:
#   1. LC/DD species with listed threats, logging/plantation threats, no listed
#      threats, and no listed threats but a written threats description.
#   2. Assessments with any listed threat, split by whether listed threat rows
#      include scope/severity information.
#
# RStudio:
#   Click Source. Outputs are written to outputs/assessment_completeness.
# =============================================================================

find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "redListAnalysis.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = TRUE)
}

project_root <- find_project_root()
if (!identical(normalizePath(getwd(), winslash = "/", mustWork = TRUE), project_root)) {
  setwd(project_root)
  cat("Working directory set to: ", project_root, "\n", sep = "")
}

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)

output_dir <- file.path("outputs", "assessment_completeness")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse the species/threat-row detail files if they already exist. Set this to
# FALSE when the raw IUCN caches have changed and you want to rebuild everything.
REUSE_EXISTING_DETAIL_FILES <- TRUE
WRITE_DETAIL_FILES <- FALSE

taxon_config <- tibble::tribble(
  ~taxon_key,    ~taxa_group,       ~raw_file,
  "bird",        "Aves",            file.path("data", "raw", "iucn", "iucn_assessment_list.rds"),
  "mammal",      "Mammalia",        file.path("data", "raw", "iucn", "mam_iucn_assessment_list.rds"),
  "amphibian",   "Amphibia",        file.path("data", "raw", "iucn", "amphib_iucn_assessment_list.rds"),
  "reptile",     "Reptilia",        file.path("data", "raw", "iucn", "rept_iucn_assessment_list.rds"),
  "fish",        "Freshwater Fish", file.path("data", "raw", "iucn", "fish_iucn_assessment_list.rds")
)

category_lc_dd <- c("LC", "DD")
logging_pattern <- "^5_3"
plantation_pattern <- "^2_2"

clean_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\uFFFD|\u00A0", " ")
  x <- str_replace_all(x, "<[^>]+>", " ")
  x <- str_squish(x)
  na_if(x, "")
}

is_present <- function(x) {
  !is.na(x) & nzchar(as.character(x))
}

extract_description_en <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.data.frame(x) && "en" %in% names(x)) return(as.character(x$en))
  if (is.list(x) && "en" %in% names(x)) return(as.character(x$en))
  as.character(x)
}

extract_assessment_metadata <- function(assessment, taxon_key, taxa_group) {
  tibble(
    taxon_key = taxon_key,
    taxa_group = taxa_group,
    scientific_name = clean_text(assessment$taxon$scientific_name),
    category = as.character(assessment$red_list_category$code),
    assessment_id = as.character(assessment$assessment_id),
    threat_description_raw = clean_text(assessment$documentation$threats)
  )
}

extract_assessment_threats <- function(assessment, taxon_key, taxa_group) {
  threat_df <- assessment$threats
  if (!is.data.frame(threat_df) || nrow(threat_df) == 0) {
    return(tibble(
      taxon_key = character(),
      taxa_group = character(),
      scientific_name = character(),
      category = character(),
      assessment_id = character(),
      threat_code = character(),
      timing = character(),
      scope = character(),
      severity = character(),
      score = character(),
      threat_description = character()
    ))
  }

  get_col <- function(name) {
    if (name %in% names(threat_df)) as.character(threat_df[[name]]) else rep(NA_character_, nrow(threat_df))
  }

  tibble(
    taxon_key = taxon_key,
    taxa_group = taxa_group,
    scientific_name = clean_text(assessment$taxon$scientific_name),
    category = as.character(assessment$red_list_category$code),
    assessment_id = as.character(assessment$assessment_id),
    threat_code = get_col("code"),
    timing = get_col("timing"),
    scope = get_col("scope"),
    severity = get_col("severity"),
    score = get_col("score"),
    threat_description = clean_text(unlist(lapply(threat_df$description, extract_description_en), use.names = FALSE))
  )
}

extract_taxon_raw <- function(taxon_key, taxa_group, raw_file) {
  cat("Reading ", taxon_key, "\n", sep = "")
  assessments <- readRDS(raw_file)
  list(
    metadata = map_dfr(assessments, extract_assessment_metadata, taxon_key = taxon_key, taxa_group = taxa_group),
    threats = map_dfr(assessments, extract_assessment_threats, taxon_key = taxon_key, taxa_group = taxa_group)
  )
}

species_flags_file <- file.path(output_dir, "assessment_threat_completeness_species_flags.csv")
raw_threat_rows_file <- file.path(output_dir, "assessment_threat_completeness_raw_threat_rows.csv")

if (
  REUSE_EXISTING_DETAIL_FILES &&
    file.exists(species_flags_file) &&
    file.exists(raw_threat_rows_file)
) {
  cat("Reading existing detail files from: ", output_dir, "\n", sep = "")
  species_flags <- read.csv(species_flags_file, stringsAsFactors = FALSE) %>%
    mutate(
      across(
        c(
          processed_has_logging_ongoing,
          processed_has_plantation_ongoing,
          has_logging_raw,
          has_plantation_raw,
          has_scope,
          has_severity,
          has_scope_or_severity,
          all_rows_have_scope_and_severity,
          has_any_threat,
          has_logging,
          has_plantation,
          has_threat_description
        ),
        as.logical
      )
    )
  raw_threat_rows <- read.csv(raw_threat_rows_file, stringsAsFactors = FALSE)
} else {
  raw_products <- pmap(taxon_config, extract_taxon_raw)
  assessment_metadata <- bind_rows(map(raw_products, "metadata"))
  raw_threat_rows <- bind_rows(map(raw_products, "threats")) %>%
    filter(is_present(threat_code))

  processed_species <- readRDS(file.path("data", "processed", "assessments", "all_taxa_assessment_species_summary.rds")) %>%
    transmute(
      taxon_key,
      taxa_group,
      scientific_name = clean_text(scientific_name),
      category = as.character(category),
      processed_threats = clean_text(threats),
      processed_threat_description = clean_text(description),
      processed_has_logging_ongoing = as.logical(has_logging_ongoing),
      processed_has_plantation_ongoing = as.logical(has_plantation_ongoing)
    )

  species_flags <- processed_species %>%
    left_join(
      assessment_metadata %>%
        select(taxon_key, scientific_name, category, assessment_id, threat_description_raw),
      by = c("taxon_key", "scientific_name", "category")
    ) %>%
    left_join(
      raw_threat_rows %>%
        group_by(taxon_key, scientific_name, category) %>%
        summarise(
          n_threat_rows_raw = n(),
          n_logging_rows_raw = sum(str_detect(threat_code, logging_pattern), na.rm = TRUE),
          n_plantation_rows_raw = sum(str_detect(threat_code, plantation_pattern), na.rm = TRUE),
          has_logging_raw = any(str_detect(threat_code, logging_pattern), na.rm = TRUE),
          has_plantation_raw = any(str_detect(threat_code, plantation_pattern), na.rm = TRUE),
          has_scope = any(is_present(scope), na.rm = TRUE),
          has_severity = any(is_present(severity), na.rm = TRUE),
          has_scope_or_severity = any(is_present(scope) | is_present(severity), na.rm = TRUE),
          all_rows_have_scope_and_severity = all(is_present(scope) & is_present(severity)),
          .groups = "drop"
        ),
      by = c("taxon_key", "scientific_name", "category")
    ) %>%
    mutate(
      n_threat_rows_raw = replace_na(n_threat_rows_raw, 0L),
      n_logging_rows_raw = replace_na(n_logging_rows_raw, 0L),
      n_plantation_rows_raw = replace_na(n_plantation_rows_raw, 0L),
      has_any_threat = is_present(processed_threats) | n_threat_rows_raw > 0,
      has_logging = processed_has_logging_ongoing | replace_na(has_logging_raw, FALSE),
      has_plantation = processed_has_plantation_ongoing | replace_na(has_plantation_raw, FALSE),
      has_threat_description = is_present(processed_threat_description) | is_present(threat_description_raw),
      has_scope = replace_na(has_scope, FALSE),
      has_severity = replace_na(has_severity, FALSE),
      has_scope_or_severity = replace_na(has_scope_or_severity, FALSE),
      all_rows_have_scope_and_severity = replace_na(all_rows_have_scope_and_severity, FALSE),
      threat_detail_class = case_when(
        !has_any_threat ~ "no_listed_threats",
        has_scope_or_severity ~ "listed_threats_with_scope_or_severity",
        TRUE ~ "listed_threats_no_scope_or_severity"
      )
    )
}

summarise_lc_dd <- function(df) {
  lcdd <- df %>% filter(category %in% category_lc_dd)
  tibble(
    n_lc_dd = nrow(lcdd),
    n_lc = sum(lcdd$category == "LC", na.rm = TRUE),
    n_dd = sum(lcdd$category == "DD", na.rm = TRUE),
    n_with_any_threat = sum(lcdd$has_any_threat, na.rm = TRUE),
    prop_with_any_threat = n_with_any_threat / n_lc_dd,
    n_with_logging = sum(lcdd$has_logging, na.rm = TRUE),
    prop_lc_dd_with_logging = n_with_logging / n_lc_dd,
    prop_threatened_lc_dd_with_logging = n_with_logging / n_with_any_threat,
    n_with_plantation = sum(lcdd$has_plantation, na.rm = TRUE),
    prop_lc_dd_with_plantation = n_with_plantation / n_lc_dd,
    prop_threatened_lc_dd_with_plantation = n_with_plantation / n_with_any_threat,
    n_with_logging_or_plantation = sum(lcdd$has_logging | lcdd$has_plantation, na.rm = TRUE),
    prop_threatened_lc_dd_with_logging_or_plantation = n_with_logging_or_plantation / n_with_any_threat,
    n_with_no_listed_threats = sum(!lcdd$has_any_threat, na.rm = TRUE),
    prop_with_no_listed_threats = n_with_no_listed_threats / n_lc_dd,
    n_no_listed_threats_but_threat_description = sum(!lcdd$has_any_threat & lcdd$has_threat_description, na.rm = TRUE),
    prop_no_listed_threats_with_threat_description = n_no_listed_threats_but_threat_description / n_with_no_listed_threats
  )
}

lc_dd_species_summary <- bind_rows(
  species_flags %>%
    group_by(taxon_key, taxa_group) %>%
    group_modify(~ summarise_lc_dd(.x)) %>%
    ungroup(),
  species_flags %>%
    summarise_lc_dd() %>%
    mutate(taxon_key = "overall", taxa_group = "Overall", .before = 1)
)

lc_dd_threat_row_summary <- bind_rows(
  raw_threat_rows %>%
    filter(category %in% category_lc_dd) %>%
    group_by(taxon_key, taxa_group) %>%
    summarise(
      n_lc_dd_threat_rows = n(),
      n_logging_threat_rows = sum(str_detect(threat_code, logging_pattern), na.rm = TRUE),
      prop_threat_rows_logging = n_logging_threat_rows / n_lc_dd_threat_rows,
      n_plantation_threat_rows = sum(str_detect(threat_code, plantation_pattern), na.rm = TRUE),
      prop_threat_rows_plantation = n_plantation_threat_rows / n_lc_dd_threat_rows,
      n_logging_or_plantation_threat_rows = sum(str_detect(threat_code, logging_pattern) | str_detect(threat_code, plantation_pattern), na.rm = TRUE),
      prop_threat_rows_logging_or_plantation = n_logging_or_plantation_threat_rows / n_lc_dd_threat_rows,
      .groups = "drop"
    ),
  raw_threat_rows %>%
    filter(category %in% category_lc_dd) %>%
    summarise(
      taxon_key = "overall",
      taxa_group = "Overall",
      n_lc_dd_threat_rows = n(),
      n_logging_threat_rows = sum(str_detect(threat_code, logging_pattern), na.rm = TRUE),
      prop_threat_rows_logging = n_logging_threat_rows / n_lc_dd_threat_rows,
      n_plantation_threat_rows = sum(str_detect(threat_code, plantation_pattern), na.rm = TRUE),
      prop_threat_rows_plantation = n_plantation_threat_rows / n_lc_dd_threat_rows,
      n_logging_or_plantation_threat_rows = sum(str_detect(threat_code, logging_pattern) | str_detect(threat_code, plantation_pattern), na.rm = TRUE),
      prop_threat_rows_logging_or_plantation = n_logging_or_plantation_threat_rows / n_lc_dd_threat_rows
    )
)

threat_detail_summary <- bind_rows(
  species_flags %>%
    filter(has_any_threat) %>%
    group_by(taxon_key, taxa_group) %>%
    summarise(
      n_assessments_with_threats = n(),
      n_no_scope_or_severity = sum(!has_scope_or_severity, na.rm = TRUE),
      prop_no_scope_or_severity = n_no_scope_or_severity / n_assessments_with_threats,
      n_with_scope_or_severity = sum(has_scope_or_severity, na.rm = TRUE),
      prop_with_scope_or_severity = n_with_scope_or_severity / n_assessments_with_threats,
      n_all_threat_rows_have_scope_and_severity = sum(all_rows_have_scope_and_severity, na.rm = TRUE),
      prop_all_threat_rows_have_scope_and_severity = n_all_threat_rows_have_scope_and_severity / n_assessments_with_threats,
      .groups = "drop"
    ),
  species_flags %>%
    filter(has_any_threat) %>%
    summarise(
      taxon_key = "overall",
      taxa_group = "Overall",
      n_assessments_with_threats = n(),
      n_no_scope_or_severity = sum(!has_scope_or_severity, na.rm = TRUE),
      prop_no_scope_or_severity = n_no_scope_or_severity / n_assessments_with_threats,
      n_with_scope_or_severity = sum(has_scope_or_severity, na.rm = TRUE),
      prop_with_scope_or_severity = n_with_scope_or_severity / n_assessments_with_threats,
      n_all_threat_rows_have_scope_and_severity = sum(all_rows_have_scope_and_severity, na.rm = TRUE),
      prop_all_threat_rows_have_scope_and_severity = n_all_threat_rows_have_scope_and_severity / n_assessments_with_threats
    )
)

category_threat_summary <- bind_rows(
  species_flags %>%
    group_by(taxon_key, taxa_group, category) %>%
    summarise(
      n_assessments = n(),
      n_with_any_threat = sum(has_any_threat, na.rm = TRUE),
      prop_with_any_threat = n_with_any_threat / n_assessments,
      n_no_listed_threats_but_threat_description = sum(!has_any_threat & has_threat_description, na.rm = TRUE),
      prop_no_listed_threats_but_threat_description = n_no_listed_threats_but_threat_description / n_assessments,
      .groups = "drop"
    ),
  species_flags %>%
    group_by(category) %>%
    summarise(
      taxon_key = "overall",
      taxa_group = "Overall",
      n_assessments = n(),
      n_with_any_threat = sum(has_any_threat, na.rm = TRUE),
      prop_with_any_threat = n_with_any_threat / n_assessments,
      n_no_listed_threats_but_threat_description = sum(!has_any_threat & has_threat_description, na.rm = TRUE),
      prop_no_listed_threats_but_threat_description = n_no_listed_threats_but_threat_description / n_assessments,
      .groups = "drop"
    ) %>%
    select(taxon_key, taxa_group, category, everything())
)

write.csv(lc_dd_species_summary, file.path(output_dir, "lc_dd_threat_assessment_species_summary.csv"), row.names = FALSE)
write.csv(lc_dd_threat_row_summary, file.path(output_dir, "lc_dd_threat_row_type_summary.csv"), row.names = FALSE)
write.csv(threat_detail_summary, file.path(output_dir, "threat_scope_severity_summary.csv"), row.names = FALSE)
write.csv(category_threat_summary, file.path(output_dir, "red_list_category_threat_summary.csv"), row.names = FALSE)

if (WRITE_DETAIL_FILES || !file.exists(species_flags_file)) {
  write.csv(species_flags, species_flags_file, row.names = FALSE)
}
if (WRITE_DETAIL_FILES || !file.exists(raw_threat_rows_file)) {
  write.csv(raw_threat_rows, raw_threat_rows_file, row.names = FALSE)
}

cat("\nWrote summaries to: ", output_dir, "\n", sep = "")
cat("\nLC/DD species summary:\n")
print(lc_dd_species_summary)
cat("\nLC/DD threat-row type summary:\n")
print(lc_dd_threat_row_summary)
cat("\nScope/severity summary among assessments with any listed threats:\n")
print(threat_detail_summary)
