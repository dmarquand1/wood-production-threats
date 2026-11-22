rm(list = ls())
gc()

# ==============================================================================
# BAYESIAN R2 VARIANCE PARTITION for a single taxon/model combination
#
# Usage: Rscript 10_bayesian_r2_partition.R <taxon> <logging|plantation>
#
# Partitions latent-scale (logit) variance into:
#   - Fixed effects   (marginal R2)
#   - Phylogeny       (semi-partial R2 for the phylo_name random effect)
#   - Spatial         (semi-partial R2 for the spatial_name random effect)
#   - Residual        (distribution-specific variance for bernoulli-logit, pi^2/3)
# Conditional R2 = Fixed + Phylogeny + Spatial.
#
# Method: de Villemereuil, Schielzeth, Nakagawa & Morrissey (2016, Genetics) /
# Nakagawa, Johnson & Schielzeth (2017, J. R. Soc. Interface), computed per
# posterior draw so every quantity carries a full Bayesian credible interval.
# phylo_name and spatial_name use correlation matrices (unit diagonal), so the
# per-observation variance contributed by each random effect is sd^2.
# ==============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 10_bayesian_r2_partition.R <taxon> <logging|plantation>")
}
taxon <- args[[1]]
model_type <- args[[2]]

model_dir <- file.path("outputs", "traits", "models", taxon)
diagnostic_dir <- file.path("outputs", "traits", "diagnostics", taxon)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)

model_path <- file.path(model_dir, sprintf("%s_combined.rds", model_type))
if (!file.exists(model_path)) {
  stop("Missing combined model file: ", model_path)
}

message(sprintf("[%s / %s] Loading %s ...", taxon, model_type, model_path))
fit <- readRDS(model_path)

stopifnot(fit$family$family == "bernoulli", fit$family$link == "logit")

message(sprintf("[%s / %s] Computing fixed-effects-only linear predictor draws...", taxon, model_type))
lp_fixed <- posterior_linpred(fit, re_formula = NA, scale = "linear")
var_fixed <- apply(lp_fixed, 1, var)
rm(lp_fixed)
gc()

message(sprintf("[%s / %s] Extracting random-effect SD draws...", taxon, model_type))
sd_draws <- as_draws_df(fit, variable = c("sd_phylo_name__Intercept", "sd_spatial_name__Intercept"))
var_phylo   <- sd_draws$sd_phylo_name__Intercept^2
var_spatial <- sd_draws$sd_spatial_name__Intercept^2

stopifnot(length(var_fixed) == length(var_phylo))

# Distribution-specific residual variance, bernoulli-logit (Nakagawa & Schielzeth 2013)
var_resid <- (pi^2) / 3

var_total <- var_fixed + var_phylo + var_spatial + var_resid

r2_fixed       <- var_fixed / var_total
r2_phylo       <- var_phylo / var_total
r2_spatial     <- var_spatial / var_total
r2_residual    <- var_resid / var_total
r2_conditional <- (var_fixed + var_phylo + var_spatial) / var_total

summarise_draws <- function(x) {
  c(median = median(x),
    q2.5   = unname(quantile(x, 0.025)),
    q97.5  = unname(quantile(x, 0.975)))
}

result <- data.frame(
  taxon = taxon,
  model = model_type,
  component = c("Fixed effects", "Phylogeny", "Spatial", "Residual", "Conditional (total)"),
  rbind(
    summarise_draws(r2_fixed),
    summarise_draws(r2_phylo),
    summarise_draws(r2_spatial),
    summarise_draws(r2_residual),
    summarise_draws(r2_conditional)
  ),
  row.names = NULL
)

out_path <- file.path(diagnostic_dir, sprintf("%s_r2_bayes_partition.csv", model_type))
write.csv(result, out_path, row.names = FALSE)
message(sprintf("[%s / %s] Saved summary: %s", taxon, model_type, out_path))

draws_out <- data.frame(
  taxon = taxon,
  model = model_type,
  r2_fixed = r2_fixed,
  r2_phylo = r2_phylo,
  r2_spatial = r2_spatial,
  r2_residual = r2_residual,
  r2_conditional = r2_conditional
)
draws_path <- file.path(diagnostic_dir, sprintf("%s_r2_bayes_partition_draws.rds", model_type))
saveRDS(draws_out, draws_path)
message(sprintf("[%s / %s] Saved posterior draws: %s", taxon, model_type, draws_path))

print(result)
