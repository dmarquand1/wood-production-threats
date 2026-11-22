library(dplyr)
library(tidyr)
library(purrr)
library(rredlist)
library(here)
library(sf)
library(ape)
library(stringr)

# ==============================================================================
# BIRD DATA PREPARATION
#
# Goal:
# Build the canonical matched bird dataset used by the active trait pipeline.
# The output pattern now mirrors the other taxa:
#   - bird_trait_data_matched.csv
#   - final_bird_tree_pruned.tre
#   - Master_Species_Crosswalk.csv
#
# Compatibility aliases are also refreshed:
#   - bird_model_data.csv
#   - Bird_Analysis_Tree.tre
#   - Analysis_Tree.tre
# ==============================================================================

THREAT_CODE_LOGGING <- "5_3_"
THREAT_CODE_PLANTATION <- "2_2_"

clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  
  x %>%
    str_replace_all("\uFFFD", " ") %>%
    str_replace_all("\u00A0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

tree_label <- function(x) {
  clean_names(x) %>%
    str_replace_all(" ", "_")
}

bird_dir <- here("data", "processed", "traits", "bird")
dir.create(bird_dir, recursive = TRUE, showWarnings = FALSE)

cat("--- Step 1: Loading Bird Trait Data ---\n")
trait_data <- read.csv(
  here("data", "raw", "traits", "bird", "BIRDBASE.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    scientific_name = clean_names(Latin..BirdLife...IOC...Clements.AviList.),
    Trait_Name = scientific_name,
    eBird_Name = coalesce(
      tree_label(eBird.Clements..V2024b.),
      tree_label(scientific_name)
    )
  )

# HWI is not present in BirdBase. The active pipeline now reconstructs it
# explicitly from AVONET1_BirdLife::Hand.Wing.Index.
avonet_hwi <- read.csv(
  here("data", "raw", "traits", "bird", "AVONET1_BirdLife.csv"),
  stringsAsFactors = FALSE
) %>%
  transmute(
    Trait_Name = clean_names(Species1),
    HWI = Hand.Wing.Index
  ) %>%
  distinct(Trait_Name, .keep_all = TRUE)

trait_data <- trait_data %>%
  left_join(avonet_hwi, by = "Trait_Name")

cat(sprintf(
  "BirdBase rows: %d | AVONET HWI matches: %d\n",
  nrow(trait_data),
  sum(!is.na(trait_data$HWI))
))

cat("--- Step 2: Loading prepared IUCN assessment backbone ---\n")
iucn_backbone <- readRDS(
  here("data", "processed", "assessments", "bird_assessment_taxonomy_lookup.rds")
) %>%
  filter(!RedList_Cat %in% c("EX", "EW")) %>%
  mutate(IUCN_Name = clean_names(IUCN_Name))

iucn_synonyms <- iucn_backbone %>%
  select(IUCN_Name, Synonyms) %>%
  tidyr::unnest_longer(Synonyms, values_to = "synonym_clean") %>%
  filter(!is.na(synonym_clean), synonym_clean != "") %>%
  mutate(synonym_clean = clean_names(synonym_clean)) %>%
  distinct(IUCN_Name, synonym_clean)

cat("--- Step 3: Loading Bird Ranges ---\n")
bird_ranges <- readRDS(here("data", "processed", "ranges", "all_bird_ranges_dissolved.rds")) %>%
  rename(scientific_name = sci_name) %>%
  mutate(scientific_name = clean_names(scientific_name))

cat("--- Step 4: Loading Bird Tree and Crosswalk ---\n")
bird_tree <- read.tree(here("data", "raw", "trees", "bird", "Bird_Consensus_Tree.tre"))

crosswalk <- read.csv(
  here("data", "raw", "traits", "bird", "BirdLife-BirdTree crosswalk.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Species1 = clean_names(Species1),
    Species3 = clean_names(Species3),
    Species3_formatted = str_replace_all(Species3, " ", "_")
  )

cat("--- Step 5: Defining Bird Data Universe ---\n")
data_species <- Reduce(intersect, list(
  unique(iucn_backbone$IUCN_Name),
  unique(trait_data$scientific_name),
  unique(bird_ranges$scientific_name)
))

cat(sprintf("Total species with IUCN + traits + ranges: %d\n", length(data_species)))

cat("--- Step 6: Crosswalk Matching ---\n")
stage1_matches <- crosswalk %>%
  filter(Species1 %in% data_species) %>%
  filter(Species3_formatted %in% bird_tree$tip.label) %>%
  transmute(
    data_name = Species1,
    tree_name = Species3_formatted,
    match_type = "Standard Crosswalk"
  )

missing_species <- setdiff(data_species, unique(stage1_matches$data_name))
cat(sprintf("Stage 1 matches: %d | still missing: %d\n",
            nrow(stage1_matches), length(missing_species)))

cat("--- Step 7: IUCN Synonym Rescue ---\n")
stage2_matches <- if (length(missing_species) > 0) {
  iucn_synonyms %>%
    filter(IUCN_Name %in% missing_species) %>%
    mutate(tree_name = str_replace_all(synonym_clean, " ", "_")) %>%
    filter(tree_name %in% bird_tree$tip.label) %>%
    transmute(
      data_name = IUCN_Name,
      tree_name = tree_name,
      match_type = "IUCN Synonym Rescue"
    ) %>%
    distinct(data_name, .keep_all = TRUE)
} else {
  tibble(data_name = character(), tree_name = character(), match_type = character())
}

cat(sprintf("Stage 2 rescues: %d\n", nrow(stage2_matches)))

cat("--- Step 8: Manual Patch Rescue ---\n")
patch_file <- here("data", "raw", "traits", "bird", "Manual_Patch_REQUIRED2.csv")

valid_patches <- if (file.exists(patch_file)) {
  read.csv(patch_file, stringsAsFactors = FALSE) %>%
    mutate(
      scientific_name = clean_names(scientific_name),
      manual_tree_name_clean = clean_names(manual_tree_name),
      tree_name = str_replace_all(manual_tree_name_clean, " ", "_")
    ) %>%
    filter(tree_name %in% bird_tree$tip.label) %>%
    transmute(
      data_name = scientific_name,
      tree_name = tree_name,
      match_type = "Manual Patch"
    ) %>%
    distinct(data_name, .keep_all = TRUE)
} else {
  tibble(data_name = character(), tree_name = character(), match_type = character())
}

cat(sprintf("Manual rescues: %d\n", nrow(valid_patches)))

cat("--- Step 9: Finalising Bird Matched Data ---\n")
final_map <- bind_rows(stage1_matches, stage2_matches, valid_patches) %>%
  distinct(data_name, .keep_all = TRUE)

final_analysis_data <- trait_data %>%
  filter(scientific_name %in% final_map$data_name) %>%
  left_join(final_map, by = c("scientific_name" = "data_name")) %>%
  left_join(iucn_backbone, by = c("scientific_name" = "IUCN_Name")) %>%
  rename(
    IUCN_Name = scientific_name,
    phylo_match = tree_name,
    match_method = match_type
  ) %>%
  mutate(
    Tree_Name = phylo_match,
    in_phylogeny = phylo_match %in% bird_tree$tip.label
  )

final_tree <- keep.tip(bird_tree, unique(final_analysis_data$phylo_match))

final_analysis_data <- final_analysis_data %>%
  arrange(match(phylo_match, final_tree$tip.label))

master_table <- final_analysis_data %>%
  transmute(
    Analysis_ID = row_number(),
    IUCN_Name = IUCN_Name,
    Trait_Name = Trait_Name,
    Polygon_Name = IUCN_Name,
    Tree_Name = phylo_match,
    eBird_Name = eBird_Name,
    Match_Method = match_method
  )

cat(sprintf("Final mapped bird species: %d\n", nrow(final_analysis_data)))

write.csv(
  final_analysis_data,
  file = here("data", "processed", "traits", "bird", "bird_trait_data_matched.csv"),
  row.names = FALSE
)

write.csv(
  final_analysis_data,
  file = here("data", "processed", "traits", "bird", "bird_model_data.csv"),
  row.names = FALSE
)

write.csv(
  master_table,
  file = here("data", "processed", "traits", "bird", "Master_Species_Crosswalk.csv"),
  row.names = FALSE
)

write.tree(final_tree, here("data", "processed", "traits", "bird", "final_bird_tree_pruned.tre"))
write.tree(final_tree, here("data", "processed", "traits", "bird", "Bird_Analysis_Tree.tre"))
write.tree(final_tree, here("data", "processed", "traits", "bird", "Analysis_Tree.tre"))

cat("Saved canonical bird matched data, pruned tree, and compatibility aliases.\n")
