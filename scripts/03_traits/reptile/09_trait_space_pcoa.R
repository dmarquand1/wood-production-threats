# =============================================================================
# REPTILE TRAIT SPACE ORDINATION
#
# Runs Gower PCoA, logging/plantation hulls, PERMANOVA, and envfit for the
# reptile trait-model dataset.
# =============================================================================

trait_space_config <- list(
  taxon = "reptile",
  taxon_label = "Reptile",
  continuous_traits = c(
    "Max_SVL_log",
    "Max_Body_Mass_log",
    "habitat_breadth_log"
  ),
  binary_traits = c(
    "Venomous",
    "Microhabitat_Arboreal",
    "Microhabitat_Fossorial",
    "Insular_or_Endemic",
    "Forest",
    "Savanna",
    "Shrubland",
    "Grassland",
    "Wetlands",
    "Rocky",
    "Caves",
    "Desert",
    "Marine",
    "Artificial",
    "Unknown"
  ),
  categorical_traits = c(
    "Active.time",
    "Reproductive.mode",
    "Fangs"
  ),
  symmetric_binary_traits = c(
    "Venomous",
    "Insular_or_Endemic"
  )
)

source(file.path("scripts", "03_traits", "trait_space_pcoa_common.R"), local = TRUE)
