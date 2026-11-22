local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(brms)

rstan::rstan_options(auto_write = TRUE)

input_dir <- file.path("data", "processed", "traits", "bird", "cr_en_vu_nt")
model_dir <- file.path("outputs", "traits", "models", "bird_cr_en_vu_nt")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

phylo_cor <- readRDS(
  file.path(input_dir, "bird_phylogeny_matrix_cr_en_vu_nt.rds")
)
spatial_cor <- readRDS(
  file.path(input_dir, "bird_correlation_matrix_cr_en_vu_nt.rds")
)
final_data <- read.csv(
  file.path(input_dir, "bird_trait_model_final_data_cr_en_vu_nt.csv"),
  stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
chain_id <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
if (is.na(chain_id) || chain_id < 1L) {
  stop("Pass a positive chain id, e.g. Rscript fit_plantation_hpc_cr_en_vu_nt.R 1")
}

allowed_categories <- c("CR", "EN", "VU", "NT")
if (!"RedList_Cat" %in% names(final_data)) {
  stop("Filtered model data is missing RedList_Cat. Run 06_prepare_cr_en_vu_nt_model_inputs.R first.")
}
if (any(!final_data$RedList_Cat %in% allowed_categories)) {
  stop("Filtered model data contains species outside CR/EN/VU/NT.")
}

model_formula <- bf(
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

fit <- brm(
  formula = model_formula,
  data = final_data,
  data2 = list(phylo_cor = phylo_cor, spatial_cor = spatial_cor),
  family = bernoulli(link = "logit"),
  init = 0,
  refresh = 10,
  chains = 1,
  iter = 1500,
  warmup = 1000,
  backend = "rstan",
  cores = 1,
  prior = priors,
  seed = chain_id
)

saveRDS(
  fit,
  file = file.path(
    model_dir,
    sprintf("bird_cr_en_vu_nt_plantation_chain%d.rds", chain_id)
  )
)
