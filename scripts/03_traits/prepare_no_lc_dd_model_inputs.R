rm(list = ls())
gc()

# =============================================================================
# PREPARE TRAIT MODEL INPUTS WITHOUT LC/DD SPECIES
#
# Creates CR/EN/VU/NT-only model inputs from existing full model inputs and
# existing covariance matrices. This script does not fit models.
#
# RStudio:
#   Edit RSTUDIO_TAXA_TO_RUN below, then click Source.
#
# Command line:
#   Rscript scripts/03_traits/prepare_no_lc_dd_model_inputs.R all
#   Rscript scripts/03_traits/prepare_no_lc_dd_model_inputs.R amphibian,fish
# =============================================================================

# ---- RSTUDIO SETTINGS --------------------------------------------------------

RSTUDIO_TAXA_TO_RUN <- c("amphibian", "fish", "mammal", "reptile")
# Options: "all", "amphibian", "fish", "mammal", "reptile", "bird",
# or a vector such as c("amphibian", "fish").

# -----------------------------------------------------------------------------

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

default_taxa <- c("amphibian", "fish", "mammal", "reptile")
all_taxa <- c(default_taxa, "bird")
target_categories <- c("CR", "EN", "VU", "NT")

normalise_taxa <- function(x, include_bird_for_all = FALSE) {
  x <- unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE))
  x <- tolower(trimws(x))
  x <- x[x != ""]
  if (length(x) == 0 || any(x == "all")) {
    return(if (include_bird_for_all) all_taxa else default_taxa)
  }
  x
}

args <- commandArgs(trailingOnly = TRUE)
taxa_to_run <- if (interactive() || length(args) == 0) {
  normalise_taxa(RSTUDIO_TAXA_TO_RUN)
} else {
  normalise_taxa(args)
}

unknown_taxa <- setdiff(taxa_to_run, all_taxa)
if (length(unknown_taxa) > 0) {
  stop("Unknown taxon/taxa: ", paste(unknown_taxa, collapse = ", "), call. = FALSE)
}

clean_species_name <- function(x) {
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

standardise_phylo_name <- function(x) {
  gsub(" ", "_", clean_species_name(x))
}

read_csv_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
}

read_rds_required <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
  readRDS(path)
}

trait_dir <- function(taxon) file.path("data", "processed", "traits", taxon)

final_data_file <- function(taxon) {
  file.path(trait_dir(taxon), sprintf("%s_trait_model_final_data.csv", taxon))
}

assessment_lookup_file <- function(taxon) {
  file.path(
    "data", "processed", "assessments",
    sprintf("%s_assessment_taxonomy_lookup.rds", taxon)
  )
}

matched_trait_file <- function(taxon) {
  file.path(trait_dir(taxon), sprintf("%s_trait_data_matched.csv", taxon))
}

phylo_matrix_file <- function(taxon) {
  file.path(trait_dir(taxon), sprintf("%s_phylogeny_matrix.rds", taxon))
}

spatial_matrix_file <- function(taxon) {
  file.path(trait_dir(taxon), sprintf("%s_correlation_matrix.rds", taxon))
}

find_fallback_categories <- function(taxon) {
  path <- matched_trait_file(taxon)
  if (!file.exists(path)) {
    return(tibble(spatial_name = character(), RedList_Cat_fallback = character()))
  }
  
  matched <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  name_col <- intersect(
    c("Trait_Name", "IUCN_Name", "scientific_name", "spatial_name"),
    names(matched)
  )[1]
  category_col <- intersect(
    c(
      "RedList_Cat",
      "category",
      "X2024.IUCN.Red.List.category",
      "red_list_category",
      "iucn_category"
    ),
    names(matched)
  )[1]
  
  if (is.na(name_col) || is.na(category_col)) {
    return(tibble(spatial_name = character(), RedList_Cat_fallback = character()))
  }
  
  matched %>%
    transmute(
      spatial_name = clean_species_name(.data[[name_col]]),
      RedList_Cat_fallback = as.character(.data[[category_col]])
    ) %>%
    filter(!is.na(spatial_name), spatial_name != "") %>%
    distinct(spatial_name, .keep_all = TRUE)
}

subset_covariance <- function(mat, keep, label) {
  keep <- unique(keep)
  missing <- setdiff(keep, rownames(mat))
  if (length(missing) > 0) {
    stop(
      label,
      " is missing ",
      length(missing),
      " level(s), first few: ",
      paste(head(missing, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- mat[keep, keep, drop = FALSE]
  if (any(!is.finite(out))) stop(label, " contains non-finite values.", call. = FALSE)
  if (!isTRUE(all.equal(out, t(out), tolerance = 1e-8))) {
    warning(label, " is not exactly symmetric; forcing symmetry.", call. = FALSE)
    out[lower.tri(out)] <- t(out)[lower.tri(out)]
  }
  diag(out) <- 1
  out
}

prepare_one_taxon <- function(taxon) {
  cat("\n========================================\n")
  cat("Preparing ", taxon, " CR/EN/VU/NT inputs\n", sep = "")
  cat("========================================\n")
  
  out_dir <- file.path(trait_dir(taxon), "cr_en_vu_nt")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  final_data <- read_csv_required(final_data_file(taxon)) %>%
    mutate(
      phylo_name = standardise_phylo_name(phylo_name),
      spatial_name = clean_species_name(spatial_name)
    ) %>%
    select(-any_of(c("RedList_Cat", "category")))
  
  required_cols <- c("phylo_name", "spatial_name", "Threat_Logging", "Threat_Plantation")
  missing_cols <- setdiff(required_cols, names(final_data))
  if (length(missing_cols) > 0) {
    stop(taxon, " final data missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  
  lookup <- read_rds_required(assessment_lookup_file(taxon)) %>%
    transmute(
      spatial_name = clean_species_name(IUCN_Name),
      RedList_Cat = as.character(RedList_Cat)
    ) %>%
    filter(!is.na(spatial_name), spatial_name != "") %>%
    distinct(spatial_name, .keep_all = TRUE)
  
  category_map <- lookup %>%
    full_join(find_fallback_categories(taxon), by = "spatial_name") %>%
    mutate(RedList_Cat = coalesce(RedList_Cat, RedList_Cat_fallback)) %>%
    select(spatial_name, RedList_Cat) %>%
    distinct(spatial_name, .keep_all = TRUE)
  
  final_with_category <- final_data %>%
    left_join(category_map, by = "spatial_name")
  
  missing_category <- final_with_category %>%
    filter(is.na(RedList_Cat) | RedList_Cat == "") %>%
    transmute(
      spatial_name,
      phylo_name,
      RedList_Cat = NA_character_,
      drop_reason = "missing_red_list_category"
    )
  
  filtered_data <- final_with_category %>%
    filter(RedList_Cat %in% target_categories)
  
  phylo_cor <- read_rds_required(phylo_matrix_file(taxon))
  spatial_cor <- read_rds_required(spatial_matrix_file(taxon))
  
  has_cov <- filtered_data %>%
    mutate(
      has_phylo_cov = phylo_name %in% rownames(phylo_cor),
      has_spatial_cov = spatial_name %in% rownames(spatial_cor),
      drop_reason = case_when(
        !has_phylo_cov & !has_spatial_cov ~ "missing_phylo_and_spatial_covariance",
        !has_phylo_cov ~ "missing_phylo_covariance",
        !has_spatial_cov ~ "missing_spatial_covariance",
        TRUE ~ NA_character_
      )
    )
  
  missing_covariance <- has_cov %>%
    filter(!is.na(drop_reason)) %>%
    transmute(spatial_name, phylo_name, RedList_Cat, drop_reason)
  
  filtered_data <- has_cov %>%
    filter(is.na(drop_reason)) %>%
    select(-has_phylo_cov, -has_spatial_cov, -drop_reason)
  
  if (nrow(filtered_data) == 0) {
    stop(taxon, ": no rows remain after CR/EN/VU/NT and matrix filtering.", call. = FALSE)
  }
  
  phylo_keep <- unique(filtered_data$phylo_name)
  spatial_keep <- unique(filtered_data$spatial_name)
  
  phylo_filtered <- subset_covariance(phylo_cor, phylo_keep, paste(taxon, "phylogeny matrix"))
  spatial_filtered <- subset_covariance(spatial_cor, spatial_keep, paste(taxon, "spatial matrix"))
  
  rm(phylo_cor, spatial_cor)
  gc()
  
  out_data <- file.path(out_dir, sprintf("%s_trait_model_final_data_cr_en_vu_nt.csv", taxon))
  out_phylo <- file.path(out_dir, sprintf("%s_phylogeny_matrix_cr_en_vu_nt.rds", taxon))
  out_spatial <- file.path(out_dir, sprintf("%s_correlation_matrix_cr_en_vu_nt.rds", taxon))
  out_summary <- file.path(out_dir, sprintf("%s_cr_en_vu_nt_model_input_summary.csv", taxon))
  out_dropped <- file.path(out_dir, sprintf("%s_cr_en_vu_nt_dropped_species.csv", taxon))
  
  write.csv(filtered_data, out_data, row.names = FALSE)
  saveRDS(phylo_filtered, out_phylo)
  saveRDS(spatial_filtered, out_spatial)
  
  dropped_species <- bind_rows(missing_category, missing_covariance)
  write.csv(dropped_species, out_dropped, row.names = FALSE)
  
  category_counts <- filtered_data %>%
    count(RedList_Cat, name = "n") %>%
    transmute(metric = paste0("category_", RedList_Cat), value = as.character(n))
  
  summary_df <- bind_rows(
    tibble::tibble(
      metric = c(
        "taxon",
        "included_red_list_categories",
        "model_data_file",
        "phylogeny_matrix_file",
        "spatial_matrix_file",
        "n_original_rows",
        "n_rows_with_category",
        "n_after_category_filter",
        "n_final_rows",
        "n_unique_phylo_levels",
        "n_unique_spatial_levels",
        "n_missing_category",
        "n_missing_covariance",
        "n_logging_threatened",
        "n_logging_not_threatened",
        "n_plantation_threatened",
        "n_plantation_not_threatened"
      ),
      value = as.character(c(
        taxon,
        paste(target_categories, collapse = ";"),
        out_data,
        out_phylo,
        out_spatial,
        nrow(final_data),
        sum(!is.na(final_with_category$RedList_Cat) & final_with_category$RedList_Cat != ""),
        nrow(final_with_category %>% filter(RedList_Cat %in% target_categories)),
        nrow(filtered_data),
        length(phylo_keep),
        length(spatial_keep),
        nrow(missing_category),
        nrow(missing_covariance),
        sum(filtered_data$Threat_Logging %in% c(TRUE, "TRUE", "True", "true", 1, "1")),
        sum(filtered_data$Threat_Logging %in% c(FALSE, "FALSE", "False", "false", 0, "0")),
        sum(filtered_data$Threat_Plantation %in% c(TRUE, "TRUE", "True", "true", 1, "1")),
        sum(filtered_data$Threat_Plantation %in% c(FALSE, "FALSE", "False", "false", 0, "0"))
      ))
    ),
    category_counts
  )
  write.csv(summary_df, out_summary, row.names = FALSE)
  
  cat(sprintf(
    "Rows: %d -> %d | phylo levels: %d | spatial levels: %d\n",
    nrow(final_data), nrow(filtered_data), length(phylo_keep), length(spatial_keep)
  ))
  cat("Wrote: ", out_dir, "\n", sep = "")
  
  rm(filtered_data, phylo_filtered, spatial_filtered)
  gc()
  invisible(out_summary)
}

for (taxon in taxa_to_run) {
  prepare_one_taxon(taxon)
}

cat("\nPrepared CR/EN/VU/NT inputs for: ", paste(taxa_to_run, collapse = ", "), "\n", sep = "")
