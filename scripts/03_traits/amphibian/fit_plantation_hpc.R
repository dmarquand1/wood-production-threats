library(brms)
library(tidyverse)

rstan::rstan_options(auto_write = TRUE)

model_dir <- file.path("outputs", "traits", "models", "amphibian")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

phylo_cor = readRDS("data/processed/traits/amphibian/amphibian_phylogeny_matrix.rds")
spatial_cor = readRDS("data/processed/traits/amphibian/amphibian_correlation_matrix.rds")
final_data = read.csv("data/processed/traits/amphibian/amphibian_trait_model_final_data.csv")

chain_id <- as.integer(commandArgs(trailingOnly = TRUE))

# 1. Force "Lar" to be the reference (Baseline)
final_data$Breeding_Cat <- relevel(as.factor(final_data$Breeding_Cat), ref = "Lar")

# ---
# 2. DEFINE MODEL FORMULA
# ---
model_formula <- bf(
  Threat_Plantation ~ 
    # Continuous
    Body_Size_log * log_range_size +
    habitat_breadth_log +
    
    # Multinomial
    Breeding_Cat +       # Reference: Lar
    
    # Binary
    Forest +
    Savanna +
    Shrubland +
    Grassland +
    Wetlands +
    
    # Random effects
    (1 | gr(phylo_name, cov = phylo_cor)) + 
    (1 | gr(spatial_name, cov = spatial_cor)) 
)

# 3. Define Priors
priors <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1),   class = "sd", group = "phylo_name"),
  prior(normal(0, 1),   class = "sd", group = "spatial_name")
)

# ---
# 4. RUN MODEL
# ---
fit <- brm(
  formula  = model_formula,
  data     = final_data,
  data2    = list(phylo_cor = phylo_cor, spatial_cor = spatial_cor),
  family   = bernoulli(link = "logit"),
  init     = 0,
  refresh  = 100,
  chains   = 1,
  iter     = 1500,
  warmup   = 500,
  backend  = "rstan",
  cores    = 1,
  #control  = list(adapt_delta = 0.95, max_treedepth = 12),
  prior    = priors,
  seed     = chain_id
)

saveRDS(
  fit,
  file = file.path(model_dir, sprintf("amphibian_plantation_chain%d.rds", chain_id))
)
