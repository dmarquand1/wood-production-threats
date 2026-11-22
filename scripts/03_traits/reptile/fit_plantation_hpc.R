library(brms)
library(tidyverse)

rstan::rstan_options(auto_write = TRUE)

model_dir <- file.path("outputs", "traits", "models", "reptile")
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

phylo_cor   <- readRDS("data/processed/traits/reptile/reptile_phylogeny_matrix.rds")
spatial_cor <- readRDS("data/processed/traits/reptile/reptile_correlation_matrix.rds")
final_data  <- read.csv("data/processed/traits/reptile/reptile_trait_model_final_data.csv")

chain_id <- as.integer(commandArgs(trailingOnly = TRUE))

# Set reference levels for categorical traits — only chain 1 writes to avoid race condition
final_data$Active.time       <- relevel(as.factor(final_data$Active.time),       ref = "Diurnal")
final_data$Reproductive.mode <- relevel(as.factor(final_data$Reproductive.mode), ref = "oviparous")
final_data$Fangs             <- relevel(as.factor(final_data$Fangs),             ref = "non-fanged")
if (chain_id == 1) {
write.csv(final_data, "data/processed/traits/reptile/reptile_trait_model_final_data.csv", row.names = FALSE)
}

model_formula <- bf(
  Threat_Plantation ~
    # Continuous
    Max_Body_Mass_log * log_range_size +
    habitat_breadth_log +
    
    # Categorical
    Active.time +
    Reproductive.mode +
    Fangs +
    
    # Binary - habitat
    Forest + Savanna + Shrubland + Grassland + Wetlands +
    
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
  iter     = 2000,
  warmup   = 1000,
  backend  = "rstan",
  cores    = 1,
  control  = list(adapt_delta = 0.95, max_treedepth = 12),
  prior    = priors,
  seed     = chain_id
)

saveRDS(
  fit,
  file = file.path(model_dir, sprintf("reptile_plantation_chain%d.rds", chain_id))
)
