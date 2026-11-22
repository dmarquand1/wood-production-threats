rm(list = ls())
gc()

# =============================================================================
# PREPARE BIRD TRAIT MODEL INPUTS FOR CR/EN/VU/NT SPECIES ONLY
#
# This script reuses the existing bird trait-model products and creates a
# parallel filtered input set for models where threat absence is interpretable:
# Critically Endangered, Endangered, Vulnerable, and Near Threatened species.
#
# It does not refit any model.
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(dplyr)

target_categories <- c("CR", "EN", "VU", "NT")

trait_dir <- file.path("data", "processed", "traits", "bird")
assessment_dir <- file.path("data", "processed", "assessments")

input_data_file <- file.path(trait_dir, "bird_trait_model_final_data.csv")
assessment_lookup_file <- file.path(
  assessment_dir,
  "bird_assessment_taxonomy_lookup.rds"
)
matched_trait_file <- file.path(trait_dir, "bird_trait_data_matched.csv")

phylo_matrix_file <- file.path(trait_dir, "bird_phylogeny_matrix.rds")
spatial_matrix_file <- file.path(trait_dir, "bird_correlation_matrix.rds")

out_dir <- file.path(trait_dir, "cr_en_vu_nt")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_data_file <- file.path(
  out_dir,
  "bird_trait_model_final_data_cr_en_vu_nt.csv"
)
out_phylo_matrix_file <- file.path(
  out_dir,
  "bird_phylogeny_matrix_cr_en_vu_nt.rds"
)
out_spatial_matrix_file <- file.path(
  out_dir,
  "bird_correlation_matrix_cr_en_vu_nt.rds"
)
out_summary_file <- file.path(
  out_dir,
  "bird_cr_en_vu_nt_model_input_summary.csv"
)
out_dropped_file <- file.path(
  out_dir,
  "bird_cr_en_vu_nt_dropped_species.csv"
)

clean_species_name <- function(x) {
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x <- gsub("[[:space:]]+", " ", x)
  trimws(x)
}

standardise_phylo_name <- function(x) {
  gsub(" ", "_", clean_species_name(x))
}

read_required <- function(path, reader) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path, call. = FALSE)
  }
  reader(path)
}

cat("\n========================================\n")
cat("Preparing bird CR/EN/VU/NT model inputs\n")
cat("========================================\n")

cat("\n1. Loading existing final bird model data...\n")
final_data <- read_required(
  input_data_file,
  function(path) read.csv(path, stringsAsFactors = FALSE, check.names = TRUE)
)

required_cols <- c(
  "phylo_name",
  "spatial_name",
  "Threat_Logging",
  "Threat_Plantation"
)
missing_cols <- setdiff(required_cols, names(final_data))
if (length(missing_cols) > 0) {
  stop(
    "Final model data is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

final_data <- final_data %>%
  mutate(
    phylo_name = standardise_phylo_name(phylo_name),
    spatial_name = clean_species_name(spatial_name)
  ) %>%
  select(-any_of(c("RedList_Cat", "category")))

cat(sprintf("   Loaded %d species from %s\n", nrow(final_data), input_data_file))

cat("\n2. Joining Red List category from processed bird assessments...\n")
assessment_lookup <- read_required(assessment_lookup_file, readRDS) %>%
  transmute(
    spatial_name = clean_species_name(IUCN_Name),
    RedList_Cat = as.character(RedList_Cat)
  ) %>%
  filter(!is.na(spatial_name), spatial_name != "") %>%
  distinct(spatial_name, .keep_all = TRUE)

category_map <- assessment_lookup

# Fallback for the small number of species that may be absent from the
# processed assessment lookup but still carry the BirdBase category column.
if (file.exists(matched_trait_file)) {
  matched_traits <- read.csv(
    matched_trait_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fallback_category_col <- intersect(
    c("RedList_Cat", "category", "X2024.IUCN.Red.List.category"),
    names(matched_traits)
  )[1]
  
  if (!is.na(fallback_category_col) && "Trait_Name" %in% names(matched_traits)) {
    fallback_categories <- matched_traits %>%
      transmute(
        spatial_name = clean_species_name(Trait_Name),
        RedList_Cat_fallback = as.character(.data[[fallback_category_col]])
      ) %>%
      filter(!is.na(spatial_name), spatial_name != "") %>%
      distinct(spatial_name, .keep_all = TRUE)
    
    category_map <- category_map %>%
      full_join(fallback_categories, by = "spatial_name") %>%
      mutate(RedList_Cat = coalesce(RedList_Cat, RedList_Cat_fallback)) %>%
      select(spatial_name, RedList_Cat) %>%
      distinct(spatial_name, .keep_all = TRUE)
  }
}

final_data_with_category <- final_data %>%
  left_join(category_map, by = "spatial_name")

missing_category <- final_data_with_category %>%
  filter(is.na(RedList_Cat) | RedList_Cat == "") %>%
  transmute(
    spatial_name,
    phylo_name,
    drop_reason = "missing_red_list_category"
  )

filtered_data <- final_data_with_category %>%
  filter(RedList_Cat %in% target_categories)

cat(sprintf(
  "   Retained %d species in categories: %s\n",
  nrow(filtered_data),
  paste(target_categories, collapse = ", ")
))

if (nrow(filtered_data) == 0) {
  stop("No species remained after category filtering.", call. = FALSE)
}

cat("\n3. Loading and subsetting existing covariance matrices...\n")
phylo_cor <- read_required(phylo_matrix_file, readRDS)
spatial_cor <- read_required(spatial_matrix_file, readRDS)

matrix_species_available <- filtered_data %>%
  mutate(
    has_phylo_cov = phylo_name %in% rownames(phylo_cor) &
      phylo_name %in% colnames(phylo_cor),
    has_spatial_cov = spatial_name %in% rownames(spatial_cor) &
      spatial_name %in% colnames(spatial_cor),
    drop_reason = case_when(
      !has_phylo_cov & !has_spatial_cov ~ "missing_phylo_and_spatial_covariance",
      !has_phylo_cov ~ "missing_phylo_covariance",
      !has_spatial_cov ~ "missing_spatial_covariance",
      TRUE ~ NA_character_
    )
  )

missing_covariance <- matrix_species_available %>%
  filter(!is.na(drop_reason)) %>%
  transmute(spatial_name, phylo_name, RedList_Cat, drop_reason)

filtered_data <- matrix_species_available %>%
  filter(is.na(drop_reason)) %>%
  select(-has_phylo_cov, -has_spatial_cov, -drop_reason)

if (nrow(filtered_data) == 0) {
  stop("No species remained after covariance matrix checks.", call. = FALSE)
}

if (anyDuplicated(filtered_data$phylo_name) > 0) {
  stop("Duplicate phylogenetic names remain after filtering.", call. = FALSE)
}

if (anyDuplicated(filtered_data$spatial_name) > 0) {
  stop("Duplicate spatial names remain after filtering.", call. = FALSE)
}

phylo_keep <- filtered_data$phylo_name
spatial_keep <- filtered_data$spatial_name

phylo_cor_filtered <- phylo_cor[phylo_keep, phylo_keep, drop = FALSE]
spatial_cor_filtered <- spatial_cor[spatial_keep, spatial_keep, drop = FALSE]

rm(phylo_cor, spatial_cor)
gc()

validate_covariance <- function(mat, matrix_name) {
  if (!identical(rownames(mat), colnames(mat))) {
    stop(matrix_name, " row and column names do not match.", call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop(matrix_name, " contains non-finite values.", call. = FALSE)
  }
  if (!isTRUE(all.equal(mat, t(mat), tolerance = 1e-8))) {
    warning(matrix_name, " is not exactly symmetric; forcing symmetry.", call. = FALSE)
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  }
  diag(mat) <- 1
  mat
}

phylo_cor_filtered <- validate_covariance(phylo_cor_filtered, "phylo_cor_filtered")
spatial_cor_filtered <- validate_covariance(spatial_cor_filtered, "spatial_cor_filtered")

cat(sprintf(
  "   Final input size: %d species | phylo matrix %d x %d | spatial matrix %d x %d\n",
  nrow(filtered_data),
  nrow(phylo_cor_filtered),
  ncol(phylo_cor_filtered),
  nrow(spatial_cor_filtered),
  ncol(spatial_cor_filtered)
))

cat("\n4. Writing filtered inputs...\n")
write.csv(filtered_data, out_data_file, row.names = FALSE)
saveRDS(phylo_cor_filtered, out_phylo_matrix_file)
saveRDS(spatial_cor_filtered, out_spatial_matrix_file)

dropped_species <- bind_rows(
  missing_category,
  missing_covariance
)
write.csv(dropped_species, out_dropped_file, row.names = FALSE)

category_counts <- filtered_data %>%
  count(RedList_Cat, name = "n_species") %>%
  mutate(metric = paste0("category_", RedList_Cat)) %>%
  transmute(metric, value = as.character(n_species))

threat_counts <- tibble::tibble(
  metric = c(
    "n_original_species",
    "n_species_with_category",
    "n_after_category_filter",
    "n_final_species",
    "n_missing_category",
    "n_missing_covariance",
    "n_logging_threatened",
    "n_logging_not_threatened",
    "n_plantation_threatened",
    "n_plantation_not_threatened"
  ),
  value = as.character(c(
    nrow(final_data),
    sum(!is.na(final_data_with_category$RedList_Cat) &
          final_data_with_category$RedList_Cat != ""),
    nrow(final_data_with_category %>% filter(RedList_Cat %in% target_categories)),
    nrow(filtered_data),
    nrow(missing_category),
    nrow(missing_covariance),
    sum(filtered_data$Threat_Logging %in% TRUE),
    sum(filtered_data$Threat_Logging %in% FALSE),
    sum(filtered_data$Threat_Plantation %in% TRUE),
    sum(filtered_data$Threat_Plantation %in% FALSE)
  ))
)

summary_df <- bind_rows(
  tibble::tibble(
    metric = c(
      "included_red_list_categories",
      "model_data_file",
      "phylogeny_matrix_file",
      "spatial_matrix_file"
    ),
    value = c(
      paste(target_categories, collapse = ";"),
      out_data_file,
      out_phylo_matrix_file,
      out_spatial_matrix_file
    )
  ),
  threat_counts,
  category_counts
)

write.csv(summary_df, out_summary_file, row.names = FALSE)

cat("\nDone. Wrote:\n")
cat("  ", out_data_file, "\n", sep = "")
cat("  ", out_phylo_matrix_file, "\n", sep = "")
cat("  ", out_spatial_matrix_file, "\n", sep = "")
cat("  ", out_summary_file, "\n", sep = "")
cat("  ", out_dropped_file, "\n", sep = "")
