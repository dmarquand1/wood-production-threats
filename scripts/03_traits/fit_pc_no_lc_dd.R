rm(list = ls())
gc()

# =============================================================================
# FIT CR/EN/VU/NT-ONLY TRAIT MODELS ON A REGULAR PC
#
# Fits logging and/or plantation models using filtered inputs prepared by:
#   scripts/03_traits/prepare_no_lc_dd_model_inputs.R
#
# Model settings:
#   chains = 4
#   warmup = 1000
#   iter = 3000  # 2000 post-warmup iterations
#
# RStudio usage:
#   1. Open this script.
#   2. Edit the RSTUDIO SETTINGS block below.
#   3. Click Source.
#
# Command-line usage from project root:
#   Rscript scripts/03_traits/fit_pc_no_lc_dd.R all both 4
#   Rscript scripts/03_traits/fit_pc_no_lc_dd.R amphibian both 4
#   Rscript scripts/03_traits/fit_pc_no_lc_dd.R fish logging 4
#
# Arguments:
#   1. taxa: all, one taxon, or a comma-separated list such as amphibian,fish
#   2. model: logging, plantation, or both. Default: both
#   3. cores: cores passed to brms. Default: min(4, detected cores)
# =============================================================================

# ---- RSTUDIO SETTINGS --------------------------------------------------------
# These are used when you click Source in RStudio, or when the script is run
# without command-line arguments.

RSTUDIO_TAXA_TO_RUN <- "all"
# Options: "all", "amphibian", "fish", "mammal", "reptile", "bird",
# or a vector such as c("amphibian", "fish").

RSTUDIO_MODEL_TO_RUN <- "both"
# Options: "both", "logging", "plantation".

RSTUDIO_CORES <- 4

# -----------------------------------------------------------------------------

find_project_root <- function(start_dir = getwd()) {
  current <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "redListAnalysis.Rproj"))) {
      return(current)
    }
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
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)

rstan::rstan_options(auto_write = TRUE)

allowed_taxa <- c("amphibian", "fish", "mammal", "reptile", "bird")
default_taxa <- c("amphibian", "fish", "mammal", "reptile", "bird")
valid_models <- c("logging", "plantation", "both")

normalise_taxa_arg <- function(x) {
  x <- unlist(strsplit(paste(x, collapse = ","), ",", fixed = TRUE))
  x <- tolower(trimws(x))
  x <- x[x != ""]
  if (length(x) == 0 || any(x == "all")) {
    return(default_taxa)
  }
  x
}

args <- commandArgs(trailingOnly = TRUE)
use_script_settings <- interactive() || length(args) == 0

taxa_to_run <- if (use_script_settings) {
  normalise_taxa_arg(RSTUDIO_TAXA_TO_RUN)
} else {
  normalise_taxa_arg(args[[1]])
}

model_to_run <- if (use_script_settings) {
  tolower(RSTUDIO_MODEL_TO_RUN)
} else if (length(args) >= 2) {
  tolower(args[[2]])
} else {
  "both"
}

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
default_cores <- min(4L, detected_cores)

requested_cores <- if (use_script_settings) {
  as.integer(RSTUDIO_CORES)
} else if (length(args) >= 3) {
  as.integer(args[[3]])
} else {
  default_cores
}
if (is.na(requested_cores) || requested_cores < 1L) {
  requested_cores <- default_cores
}
requested_cores <- min(requested_cores, detected_cores)
options(mc.cores = requested_cores)

unknown_taxa <- setdiff(taxa_to_run, allowed_taxa)
if (length(unknown_taxa) > 0) {
  stop(
    "Unknown taxon/taxa: ",
    paste(unknown_taxa, collapse = ", "),
    "\nAllowed taxa: all, ",
    paste(allowed_taxa, collapse = ", "),
    call. = FALSE
  )
}

if (!model_to_run %in% valid_models) {
  stop("Model must be one of: ", paste(valid_models, collapse = ", "), call. = FALSE)
}

models_to_run <- if (model_to_run == "both") {
  c("logging", "plantation")
} else {
  model_to_run
}

to_binary_int <- function(x) {
  as.integer(x %in% c(TRUE, "TRUE", "True", "true", 1, "1"))
}

set_reference <- function(x, ref, taxon, column) {
  x <- as.factor(x)
  if (ref %in% levels(x)) {
    return(relevel(x, ref = ref))
  }
  warning(
    taxon,
    ": reference level '",
    ref,
    "' not present in ",
    column,
    "; using default factor baseline.",
    call. = FALSE
  )
  x
}

taxon_formula_rhs <- function(taxon) {
  switch(
    taxon,
    amphibian = paste(
      "Body_Size_log * log_range_size",
      "habitat_breadth_log",
      "Breeding_Cat",
      "Forest + Savanna + Shrubland + Grassland + Wetlands",
      "(1 | gr(phylo_name, cov = phylo_cor))",
      "(1 | gr(spatial_name, cov = spatial_cor))",
      sep = " + "
    ),
    fish = paste(
      "MaxLength_Log * log_range_size",
      "Salinity_Breadth_sc",
      "Trophic_Level_sc",
      "PredPreyRatioMax_log_sc",
      "DemersPelag_Cat",
      "AirBreathing_Bin",
      "(1 | gr(phylo_name, cov = phylo_cor))",
      "(1 | gr(spatial_name, cov = spatial_cor))",
      sep = " + "
    ),
    mammal = paste(
      "Body_Mass_Log * log_range_size",
      "Habitat_Breadth",
      "Litter_Size_Log",
      "activity_cycle",
      "Forest + Savanna + Shrubland + Grassland + Wetlands",
      "binary_det_inv + binary_det_vend + binary_det_vect",
      "binary_det_vunk + binary_det_scav + binary_det_fruit",
      "binary_det_nect + binary_det_seed + binary_det_plantother",
      "(1 | gr(phylo_name, cov = phylo_cor))",
      "(1 | gr(spatial_name, cov = spatial_cor))",
      sep = " + "
    ),
    reptile = paste(
      "Max_Body_Mass_log * log_range_size",
      "habitat_breadth_log",
      "Active.time",
      "Reproductive.mode",
      "Fangs",
      "Forest + Savanna + Shrubland + Grassland + Wetlands",
      "(1 | gr(phylo_name, cov = phylo_cor))",
      "(1 | gr(spatial_name, cov = spatial_cor))",
      sep = " + "
    ),
    bird = paste(
      "Body_Mass_log * log_range_size",
      "Habitat_Breadth_log",
      "Clutch_Size_log",
      "HWI_log",
      "Forest + Savanna + Shrubland + Grassland + Wetlands",
      "Diet_Invertebrate + Diet_Fruit + Diet_Nectar + Diet_Seed",
      "Diet_Vertebrate + Diet_Fish + Diet_Scavenger + Diet_PlantMaterial",
      "(1 | gr(phylo_name, cov = phylo_cor))",
      "(1 | gr(spatial_name, cov = spatial_cor))",
      sep = " + "
    )
  )
}

required_predictors <- function(taxon) {
  switch(
    taxon,
    amphibian = c(
      "Body_Size_log", "log_range_size", "habitat_breadth_log",
      "Breeding_Cat", "Forest", "Savanna", "Shrubland", "Grassland", "Wetlands"
    ),
    fish = c(
      "MaxLength_Log", "log_range_size", "Salinity_Breadth_sc",
      "Trophic_Level_sc", "PredPreyRatioMax_log_sc",
      "DemersPelag_Cat", "AirBreathing_Bin"
    ),
    mammal = c(
      "Body_Mass_Log", "log_range_size", "Habitat_Breadth",
      "Litter_Size_Log", "activity_cycle", "Forest", "Savanna",
      "Shrubland", "Grassland", "Wetlands", "binary_det_inv",
      "binary_det_vend", "binary_det_vect", "binary_det_vunk",
      "binary_det_scav", "binary_det_fruit", "binary_det_nect",
      "binary_det_seed", "binary_det_plantother"
    ),
    reptile = c(
      "Max_Body_Mass_log", "log_range_size", "habitat_breadth_log",
      "Active.time", "Reproductive.mode", "Fangs", "Forest", "Savanna",
      "Shrubland", "Grassland", "Wetlands"
    ),
    bird = c(
      "Body_Mass_log", "log_range_size", "Habitat_Breadth_log",
      "Clutch_Size_log", "HWI_log", "Forest", "Savanna", "Shrubland",
      "Grassland", "Wetlands", "Diet_Invertebrate", "Diet_Fruit",
      "Diet_Nectar", "Diet_Seed", "Diet_Vertebrate", "Diet_Fish",
      "Diet_Scavenger", "Diet_PlantMaterial"
    )
  )
}

apply_reference_levels <- function(final_data, taxon) {
  if (taxon == "amphibian") {
    final_data$Breeding_Cat <- set_reference(
      final_data$Breeding_Cat,
      ref = "Lar",
      taxon = taxon,
      column = "Breeding_Cat"
    )
  }
  if (taxon == "fish") {
    final_data$DemersPelag_Cat <- set_reference(
      final_data$DemersPelag_Cat,
      ref = "benthopelagic",
      taxon = taxon,
      column = "DemersPelag_Cat"
    )
  }
  if (taxon == "mammal") {
    final_data$activity_cycle <- set_reference(
      final_data$activity_cycle,
      ref = "Nocturnal",
      taxon = taxon,
      column = "activity_cycle"
    )
  }
  if (taxon == "reptile") {
    final_data$Active.time <- set_reference(
      final_data$Active.time,
      ref = "Diurnal",
      taxon = taxon,
      column = "Active.time"
    )
    final_data$Reproductive.mode <- set_reference(
      final_data$Reproductive.mode,
      ref = "oviparous",
      taxon = taxon,
      column = "Reproductive.mode"
    )
    final_data$Fangs <- set_reference(
      final_data$Fangs,
      ref = "non-fanged",
      taxon = taxon,
      column = "Fangs"
    )
  }
  final_data
}

make_formula <- function(response, taxon) {
  bf(as.formula(paste(response, "~", taxon_formula_rhs(taxon))))
}

priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1), class = "sd", group = "phylo_name"),
  prior(normal(0, 1), class = "sd", group = "spatial_name")
)

fit_one_taxon <- function(taxon) {
  input_dir <- file.path("data", "processed", "traits", taxon, "cr_en_vu_nt")
  model_dir <- file.path("outputs", "traits", "models", paste0(taxon, "_cr_en_vu_nt"))
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  data_file <- file.path(
    input_dir,
    sprintf("%s_trait_model_final_data_cr_en_vu_nt.csv", taxon)
  )
  phylo_file <- file.path(
    input_dir,
    sprintf("%s_phylogeny_matrix_cr_en_vu_nt.rds", taxon)
  )
  spatial_file <- file.path(
    input_dir,
    sprintf("%s_correlation_matrix_cr_en_vu_nt.rds", taxon)
  )
  
  if (!file.exists(data_file) || !file.exists(phylo_file) || !file.exists(spatial_file)) {
    stop(
      "Filtered inputs are missing for ",
      taxon,
      ". Run scripts/03_traits/prepare_no_lc_dd_model_inputs.R ",
      taxon,
      " first.",
      call. = FALSE
    )
  }
  
  cat("\n========================================\n")
  cat("Fitting ", taxon, " CR/EN/VU/NT trait model(s)\n", sep = "")
  cat("========================================\n")
  cat("Model request: ", model_to_run, "\n", sep = "")
  cat("Chains: 4 | Warmup: 500 | Post-warmup: 1000 | Cores: ",
      requested_cores, "\n", sep = "")
  
  phylo_cor <- readRDS(phylo_file)
  spatial_cor <- readRDS(spatial_file)
  final_data <- read.csv(data_file, stringsAsFactors = FALSE)
  
  required_cols <- c(
    "phylo_name",
    "spatial_name",
    "Threat_Logging",
    "Threat_Plantation",
    required_predictors(taxon)
  )
  missing_cols <- setdiff(required_cols, names(final_data))
  if (length(missing_cols) > 0) {
    stop(
      taxon,
      " filtered model data is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!all(final_data$phylo_name %in% rownames(phylo_cor))) {
    stop(taxon, ": some phylo_name values are missing from the phylogeny matrix.")
  }
  if (!all(final_data$spatial_name %in% rownames(spatial_cor))) {
    stop(taxon, ": some spatial_name values are missing from the spatial matrix.")
  }
  
  phylo_keep <- unique(final_data$phylo_name)
  spatial_keep <- unique(final_data$spatial_name)
  
  phylo_cor <- phylo_cor[phylo_keep, phylo_keep, drop = FALSE]
  spatial_cor <- spatial_cor[spatial_keep, spatial_keep, drop = FALSE]
  
  final_data <- apply_reference_levels(final_data, taxon)
  final_data$Threat_Logging <- to_binary_int(final_data$Threat_Logging)
  final_data$Threat_Plantation <- to_binary_int(final_data$Threat_Plantation)
  
  cat(sprintf("Input rows: %d\n", nrow(final_data)))
  if ("RedList_Cat" %in% names(final_data)) {
    cat("Red List categories:\n")
    print(table(final_data$RedList_Cat))
  }
  cat("Threat counts:\n")
  print(table(logging = final_data$Threat_Logging))
  print(table(plantation = final_data$Threat_Plantation))
  
  fit_one_model <- function(model_name, response, seed_offset) {
    cat("\n----------------------------------------\n")
    cat("Fitting ", taxon, " ", model_name, " model\n", sep = "")
    cat("----------------------------------------\n")
    
    fit <- brm(
      formula = make_formula(response, taxon),
      data = final_data,
      data2 = list(phylo_cor = phylo_cor, spatial_cor = spatial_cor),
      family = bernoulli(link = "logit"),
      init = 0,
      refresh = 10,
      chains = 4,
      iter = 3000,
      warmup = 1000,
      backend = "rstan",
      cores = requested_cores,
      prior = priors,
      seed = 20260429L + seed_offset
    )
    
    combined_path <- file.path(model_dir, sprintf("%s_combined.rds", model_name))
    pc_path <- file.path(
      model_dir,
      sprintf("%s_cr_en_vu_nt_%s_pc_4chains.rds", taxon, model_name)
    )
    
    saveRDS(fit, combined_path)
    saveRDS(fit, pc_path)
    
    cat("Saved:\n")
    cat("  ", combined_path, "\n", sep = "")
    cat("  ", pc_path, "\n", sep = "")
    
    rm(fit)
    gc()
    invisible(NULL)
  }
  
  if ("logging" %in% models_to_run) {
    fit_one_model("logging", "Threat_Logging", seed_offset = 100L)
  }
  
  if ("plantation" %in% models_to_run) {
    fit_one_model("plantation", "Threat_Plantation", seed_offset = 200L)
  }
  
  rm(phylo_cor, spatial_cor, final_data)
  gc()
  invisible(NULL)
}

cat("\nRequested taxa: ", paste(taxa_to_run, collapse = ", "), "\n", sep = "")
cat("Requested model(s): ", paste(models_to_run, collapse = ", "), "\n", sep = "")

for (taxon in taxa_to_run) {
  fit_one_taxon(taxon)
}

cat("\nAll requested PC model fits complete.\n")
