# =============================================================================
# FISH TRAIT SPACE ORDINATION
#
# Runs Gower PCoA, logging/plantation hulls, PERMANOVA, and envfit for the fish
# trait-model dataset.
# =============================================================================

trait_space_config <- list(
  taxon = "fish",
  taxon_label = "Freshwater fish",
  continuous_traits = c(
    "MaxLength_Log",
    "Trophic_Level_sc",
    "PredPreyRatioMax_log_sc",
    "Salinity_Breadth_sc"
  ),
  binary_traits = c(
    "Fresh",
    "Brack",
    "Salt",
    "AirBreathing_Bin",
    "Electrogenic_Bin",
    "Dangerous_Bin",
    "Wetlands",
    "Marine",
    "Artificial"
  ),
  categorical_traits = c(
    "DemersPelag_Cat",
    "EnvTemp_Cat"
  ),
  symmetric_binary_traits = c(
    "AirBreathing_Bin",
    "Electrogenic_Bin",
    "Dangerous_Bin"
  )
)

source(file.path("scripts", "03_traits", "trait_space_pcoa_common.R"), local = TRUE)
