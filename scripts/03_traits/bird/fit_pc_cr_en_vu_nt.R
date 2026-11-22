rm(list = ls())
gc()

# =============================================================================
# FIT BIRD CR/EN/VU/NT TRAIT MODELS ON A REGULAR PC
#
# Runs the filtered bird trait models in one local brms call per threat type:
#   - 4 chains
#   - 500 warmup iterations
#   - 1000 post-warmup iterations
#
# Usage from the project root:
#   Rscript scripts/03_traits/bird/fit_pc_cr_en_vu_nt.R
#   Rscript scripts/03_traits/bird/fit_pc_cr_en_vu_nt.R logging
#   Rscript scripts/03_traits/bird/fit_pc_cr_en_vu_nt.R plantation
#   Rscript scripts/03_traits/bird/fit_pc_cr_en_vu_nt.R both 4
#
# Optional arguments:
#   1. model to run: logging, plantation, or both. Default: both
#   2. cores to use. Default: min(4, detected cores)
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)

rstan::rstan_options(auto_write = TRUE)

args <- commandArgs(trailingOnly = TRUE)
model_to_run <- if (length(args) >= 1) tolower(args[[1]]) else "both"

valid_models <- c("logging", "plantation", "both")
if (!model_to_run %in% valid_models) {
  stop(
    "First argument must be one of: ",
    paste(valid_models, collapse = ", "),
    call. = FALSE
  )
}

detected_cores <- parallel::detectCores(logical = TRUE)
if (is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
default_cores <- min(4L, detected_cores)
requested_cores <- if (length(args) >= 2) as.integer(args[[2]]) else default_cores
if (is.na(requested_cores) || requested_cores < 1L) {
  requested_cores <- default_cores
}
requested_cores <- min(requested_cores, detected_cores)

options(mc.cores = requested_cores)

input_dir <- file.path("data", "processed", "traits", "bird", "cr_en_vu_nt")
model_dir <- file.path("outputs", "traits", "models", "bird_cr_en_vu_nt")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- file.path(input_dir, "bird_trait_model_final_data_cr_en_vu_nt.csv")
phylo_file <- file.path(input_dir, "bird_phylogeny_matrix_cr_en_vu_nt.rds")
spatial_file <- file.path(input_dir, "bird_correlation_matrix_cr_en_vu_nt.rds")

if (!file.exists(data_file) || !file.exists(phylo_file) || !file.exists(spatial_file)) {
  stop(
    "Filtered model inputs are missing. Run ",
    "scripts/03_traits/bird/06_prepare_cr_en_vu_nt_model_inputs.R first.",
    call. = FALSE
  )
}

cat("\n========================================\n")
cat("Fitting bird CR/EN/VU/NT trait model(s)\n")
cat("========================================\n")
cat("Model request: ", model_to_run, "\n", sep = "")
cat("Chains: 4 | Warmup: 500 | Post-warmup: 1000 | Cores: ",
    requested_cores, "\n", sep = "")

phylo_cor <- readRDS(phylo_file)
spatial_cor <- readRDS(spatial_file)
final_data <- read.csv(data_file, stringsAsFactors = FALSE)

allowed_categories <- c("CR", "EN", "VU", "NT")
if (!"RedList_Cat" %in% names(final_data)) {
  stop("Filtered model data is missing RedList_Cat.", call. = FALSE)
}
if (any(!final_data$RedList_Cat %in% allowed_categories)) {
  stop("Filtered model data contains species outside CR/EN/VU/NT.", call. = FALSE)
}

required_cols <- c(
  "phylo_name", "spatial_name",
  "Threat_Logging", "Threat_Plantation",
  "Body_Mass_log", "log_range_size", "Habitat_Breadth_log",
  "Clutch_Size_log", "HWI_log",
  "Forest", "Savanna", "Shrubland", "Grassland", "Wetlands",
  "Diet_Invertebrate", "Diet_Fruit", "Diet_Nectar", "Diet_Seed",
  "Diet_Vertebrate", "Diet_Fish", "Diet_Scavenger", "Diet_PlantMaterial"
)
missing_cols <- setdiff(required_cols, names(final_data))
if (length(missing_cols) > 0) {
  stop(
    "Filtered model data is missing required column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!all(final_data$phylo_name %in% rownames(phylo_cor))) {
  stop("Some phylo_name values are missing from the phylogeny matrix.")
}
if (!all(final_data$spatial_name %in% rownames(spatial_cor))) {
  stop("Some spatial_name values are missing from the spatial matrix.")
}

phylo_cor <- phylo_cor[final_data$phylo_name, final_data$phylo_name, drop = FALSE]
spatial_cor <- spatial_cor[
  final_data$spatial_name,
  final_data$spatial_name,
  drop = FALSE
]

final_data$Threat_Logging <- as.integer(final_data$Threat_Logging %in% TRUE)
final_data$Threat_Plantation <- as.integer(final_data$Threat_Plantation %in% TRUE)

cat(sprintf("Input species: %d\n", nrow(final_data)))
cat("Red List categories:\n")
print(table(final_data$RedList_Cat))
cat("Threat counts:\n")
print(table(logging = final_data$Threat_Logging))
print(table(plantation = final_data$Threat_Plantation))

model_formula_logging <- bf(
  Threat_Logging ~
    # Continuous
    Body_Mass_log * log_range_size +
    Habitat_Breadth_log +
    Clutch_Size_log +
    HWI_log +
    
    # Binary - habitat
    Forest + Savanna + Shrubland + Grassland + Wetlands +
    
    # Binary - diet
    Diet_Invertebrate + Diet_Fruit + Diet_Nectar + Diet_Seed +
    Diet_Vertebrate + Diet_Fish + Diet_Scavenger + Diet_PlantMaterial +
    
    # Random effects
    (1 | gr(phylo_name, cov = phylo_cor)) +
    (1 | gr(spatial_name, cov = spatial_cor))
)

model_formula_plantation <- bf(
  Threat_Plantation ~
    # Continuous
    Body_Mass_log * log_range_size +
    Habitat_Breadth_log +
    Clutch_Size_log +
    HWI_log +
    
    # Binary - habitat
    Forest + Savanna + Shrubland + Grassland + Wetlands +
    
    # Binary - diet
    Diet_Invertebrate + Diet_Fruit + Diet_Nectar + Diet_Seed +
    Diet_Vertebrate + Diet_Fish + Diet_Scavenger + Diet_PlantMaterial +
    
    # Random effects
    (1 | gr(phylo_name, cov = phylo_cor)) +
    (1 | gr(spatial_name, cov = spatial_cor))
)

priors <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1), class = "sd", group = "phylo_name"),
  prior(normal(0, 1), class = "sd", group = "spatial_name")
)

fit_one_model <- function(model_name, model_formula, seed_offset) {
  cat("\n----------------------------------------\n")
  cat("Fitting ", model_name, " model\n", sep = "")
  cat("----------------------------------------\n")
  
  fit <- brm(
    formula = model_formula,
    data = final_data,
    data2 = list(phylo_cor = phylo_cor, spatial_cor = spatial_cor),
    family = bernoulli(link = "logit"),
    init = 0,
    refresh = 10,
    chains = 4,
    iter = 1500,
    warmup = 500,
    backend = "rstan",
    cores = requested_cores,
    prior = priors,
    seed = 20260429L + seed_offset
  )
  
  combined_path <- file.path(model_dir, sprintf("%s_combined.rds", model_name))
  pc_path <- file.path(
    model_dir,
    sprintf("bird_cr_en_vu_nt_%s_pc_4chains.rds", model_name)
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

if (model_to_run %in% c("logging", "both")) {
  fit_one_model("logging", model_formula_logging, seed_offset = 1L)
}

if (model_to_run %in% c("plantation", "both")) {
  fit_one_model("plantation", model_formula_plantation, seed_offset = 2L)
}

cat("\nAll requested PC model fits complete.\n")
