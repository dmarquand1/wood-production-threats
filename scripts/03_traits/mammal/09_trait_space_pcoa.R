# =============================================================================
# MAMMAL TRAIT SPACE ORDINATION
#
# Runs Gower PCoA, logging/plantation hulls, PERMANOVA, and envfit for the
# mammal trait-model dataset.
# =============================================================================

trait_space_config <- list(
  taxon = "mammal",
  taxon_label = "Mammal",
  continuous_traits = c(
    "Litter_Size_Log",
    "Dispersal_Log",
    "Diet_Breadth",
    "Body_Mass_Log",
    "Habitat_Breadth"
  ),
  binary_traits = c(
    "binary_det_inv",
    "binary_det_vend",
    "binary_det_vect",
    "binary_det_vunk",
    "binary_det_scav",
    "binary_det_fruit",
    "binary_det_nect",
    "binary_det_seed",
    "binary_det_plantother",
    "binary_dphy_invertebrate",
    "binary_dphy_vertebrate",
    "binary_dphy_plant",
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
    "activity_cycle"
  ),
  symmetric_binary_traits = character()
)

source(file.path("scripts", "03_traits", "trait_space_pcoa_common.R"), local = TRUE)
