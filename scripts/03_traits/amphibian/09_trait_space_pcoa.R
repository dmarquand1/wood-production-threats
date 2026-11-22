# =============================================================================
# AMPHIBIAN TRAIT SPACE ORDINATION
#
# Runs Gower PCoA, logging/plantation hulls, PERMANOVA, and envfit for the
# amphibian trait-model dataset.
# =============================================================================

trait_space_config <- list(
  taxon = "amphibian",
  taxon_label = "Amphibian",
  continuous_traits = c(
    "Body_Size_log",
    "habitat_breadth_log"
  ),
  binary_traits = c(
    "Forest",
    "Savanna",
    "Shrubland",
    "Grassland",
    "Wetlands",
    "Rocky",
    "Caves",
    "Artificial"
  ),
  categorical_traits = c(
    "Breeding_Cat"
  ),
  symmetric_binary_traits = character()
)

source(file.path("scripts", "03_traits", "trait_space_pcoa_common.R"), local = TRUE)
