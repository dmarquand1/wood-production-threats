library(brms)
library(tidyverse)

rstan::rstan_options(auto_write = TRUE)

model_dir <- file.path("outputs", "traits", "models", "fish")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

phylo_cor   <- readRDS("data/processed/traits/fish/fish_phylogeny_matrix.rds")
spatial_cor <- readRDS("data/processed/traits/fish/fish_correlation_matrix.rds")
final_data  <- read.csv("data/processed/traits/fish/fish_trait_model_final_data.csv")

chain_id <- as.integer(commandArgs(trailingOnly = TRUE))

# Set reference levels for categorical traits — only chain 1 writes to avoid race condition
final_data$DemersPelag_Cat <- relevel(as.factor(final_data$DemersPelag_Cat), ref = "benthopelagic")
if (chain_id == 1) {
write.csv(final_data, "data/processed/traits/fish/fish_trait_model_final_data.csv", row.names = FALSE)
}

model_formula <- bf(
  Threat_Plantation ~
    # Continuous
    MaxLength_Log * log_range_size +
    Salinity_Breadth_sc +
    Trophic_Level_sc +
    PredPreyRatioMax_log_sc +
    
    # Categorical
    DemersPelag_Cat +
    
    # Binary
    AirBreathing_Bin +
    
    # Random effects
    (1 | gr(phylo_name,   cov = phylo_cor)) +
    (1 | gr(spatial_name, cov = spatial_cor))
)

priors <- c(
  prior(normal(0, 1),   class = "Intercept"),
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1),   class = "sd", group = "phylo_name"),
  prior(normal(0, 1),   class = "sd", group = "spatial_name")
)

fit <- brm(
  formula  = model_formula,
  data     = final_data,
  data2    = list(phylo_cor = phylo_cor, spatial_cor = spatial_cor),
  family   = bernoulli(link = "logit"),
  init     = 0,
  refresh  = 10,
  chains   = 1,
  iter     = 1500,
  warmup   = 1000,
  backend  = "rstan",
  cores    = 1,
  #control  = list(adapt_delta = 0.95, max_treedepth = 12),
  prior    = priors,
  seed     = chain_id
)

saveRDS(
  fit,
  file = file.path(model_dir, sprintf("fish_plantation_chain%d.rds", chain_id))
)
