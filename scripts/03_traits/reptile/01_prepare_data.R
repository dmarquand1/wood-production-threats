rm(list=ls())
gc()

library(tidyverse)
library(rredlist)
library(sf)
library(ape)
library(here)
library(stringr)
library(phangorn)
library(stringdist)

# ==============================================================================
# 0. CONFIGURATION & API KEYS
# ==============================================================================
# Ensure your IUCN Key is set. If not, uncomment and set it here:
# Sys.setenv(IUCN_REDLIST_KEY = "your_key_here")

# Define Threat Codes for Logging (5.3) and Wood Harvesting (2.2)
THREAT_CODE_LOGGING <- "5_3_"
THREAT_CODE_PLANTATION  <- "2_2_"

# Helper: Clean Names (Standardize spaces/chars)
clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  x %>%
    str_replace_all("\uFFFD", " ") %>%
    str_replace_all("\u00A0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

# ==============================================================================
# 1. LOAD PREPARED IUCN RED LIST BACKBONE
# ==============================================================================
cat("--- Step 1: Loading prepared IUCN assessment backbone ---\n")
iucn_clean <- readRDS(
  here::here(
    "data", "processed", "assessments", "reptile_assessment_taxonomy_lookup.rds"
  )
) %>%
  filter(!RedList_Cat %in% c("EX", "EW")) %>%
  mutate(
    IUCN_Name = clean_names(IUCN_Name),
    Synonyms = purrr::map(Synonyms, clean_names)
  )

cat(sprintf("   > IUCN Backbone Ready: %d species.\n", nrow(iucn_clean)))

# ==============================================================================
# 2. LOAD TRAIT DATA (ReptTraits)
# ==============================================================================
cat("--- Step 2: Loading Trait Data ---\n")
trait_df <- read.csv("data/raw/traits/reptile/ReptTraits dataset v1-2.csv", stringsAsFactors = FALSE)

# ==============================================================================
# 3. LOAD RANGE POLYGONS
# ==============================================================================
cat("--- Step 3: Loading Range Polygons ---\n")
# We only need the list of names for now to check availability
if(file.exists(here("data/processed/ranges/all_reptile_ranges_dissolved.rds"))) {
  ranges <- readRDS(here("data/processed/ranges/all_reptile_ranges_dissolved.rds"))
  # If sf object, drop geometry to speed up 'names check'
  range_names <- ranges$sci_name # Ensure column name matches
  rm(ranges) # Free memory
  cat("   > Range data found for", length(range_names), "species.\n")
} else {
  stop("Range file missing!")
}

# ==============================================================================
# 4. LOAD PHYLOGENY (Tonini et al 2016)
# ==============================================================================
cat("--- Step 4: Fetching Jetz & Pyron Tree ---\n")
# 1. Read trees
my_trees <- read.nexus("data/raw/trees/reptiletree/output.nex")

# 2. Identify the MCC tree
tree <- maxCladeCred(my_trees)
tree_tips = tree$tip.label

# 3. Plot tree
plot(tree, type = "fan", show.tip.label = FALSE, lwd = 0.5)

# Add a legend for the timeline/scale
add.scale.bar()
axisPhylo()

zoom(tree, focus = 1:50, col = "red")

# ==============================================================================
# 5. TAXONOMIC OVERLAP CHECK (Traits vs. Phylogeny)
# ==============================================================================
cat("--- Step 5: Checking Name Overlap ---\n")

# 1. Standardize formatting
# Most trees use underscores; reptIO uses spaces. 
# We'll create a cleaned version of the trait names to match the tree style.
trait_species_cleaned <- trait_df$Species %>%
  clean_names() %>%
  str_replace_all(" ", "_")

tree_tips_cleaned <- tree_tips %>% 
  clean_names()

# 2. Identify matches and mismatches
matches    <- intersect(trait_species_cleaned, tree_tips_cleaned)
only_trait <- setdiff(trait_species_cleaned, tree_tips_cleaned)
only_tree  <- setdiff(tree_tips_cleaned, trait_species_cleaned)

# 3. Calculate statistics
n_trait <- length(trait_species_cleaned)
n_tree  <- length(tree_tips_cleaned)
n_match <- length(matches)
perc_match <- (n_match / n_trait) * 100

# 4. Print Summary Report
cat("\n--- Overlap Summary ---\n")
cat(sprintf("Total Species in Trait Dataset:  %d\n", n_trait))
cat(sprintf("Total Tips in Phylogeny:        %d\n", n_tree))
cat(sprintf("Number of Matches:              %d\n", n_match))
cat(sprintf("Percentage of Traits Matched:   %.2f%%\n", perc_match))

if(length(only_trait) > 0) {
  cat(sprintf("\nNote: %d species in trait data were NOT found in the tree.\n", length(only_trait)))
  cat("First 5 mismatches:", head(only_trait, 5), "...\n")
}

# ==============================================================================
# 6. ADVANCED TAXONOMIC RESOLUTION (Direct -> Synonym -> Fuzzy)
# ==============================================================================
cat("--- Step 6: Resolving Taxonomy via Synonyms & Fuzzy Matching ---\n")

# Ensure required package is loaded
if (!require("stringdist")) install.packages("stringdist")
library(stringdist) 

# ------------------------------------------------------------------------------
# A. PREPARE THE DICTIONARY (Synonyms)
# ------------------------------------------------------------------------------
# 1. Clean Tree Tips
tree_tips_cleaned <- clean_names(tree_tips) %>% str_replace_all(" ", "_")
valid_tips_df <- data.frame(Tree_Tip = tree_tips_cleaned, stringsAsFactors = FALSE)

# 2. Build IUCN Synonym Dictionary
#    (Matches IUCN Synonyms to Valid Tree Tips)
synonyms_long <- iucn_clean %>%
  dplyr::select(IUCN_Name, Synonyms) %>%
  unnest(Synonyms) %>%
  mutate(Search_Name = clean_names(Synonyms) %>% str_replace_all(" ", "_")) %>%
  dplyr::select(IUCN_Name, Search_Name)

self_ref <- iucn_clean %>%
  mutate(Search_Name = clean_names(IUCN_Name) %>% str_replace_all(" ", "_")) %>%
  dplyr::select(IUCN_Name, Search_Name)

master_name_dictionary <- bind_rows(synonyms_long, self_ref) %>%
  distinct() %>%
  filter(!is.na(Search_Name) & Search_Name != "")

# Create the link: Any IUCN Name -> A Valid Tree Tip
iucn_to_tree_link <- master_name_dictionary %>%
  rename(Possible_Tree_Tip = Search_Name) %>%
  inner_join(valid_tips_df, by = c("Possible_Tree_Tip" = "Tree_Tip")) %>%
  group_by(IUCN_Name) %>%
  slice(1) %>%
  ungroup() %>%
  rename(Valid_Tree_Tip = Possible_Tree_Tip) %>%
  dplyr::select(IUCN_Name, Valid_Tree_Tip)

# ------------------------------------------------------------------------------
# B. PREPARE FUZZY MATCHING (Direct trait → tree)
# ------------------------------------------------------------------------------
trait_prep <- trait_df %>%
  mutate(Trait_Name_Clean = clean_names(Species) %>% str_replace_all(" ", "_"))

potential_fuzzy_candidates <- trait_prep %>%
  filter(!Trait_Name_Clean %in% tree_tips_cleaned) %>%
  pull(Trait_Name_Clean) %>%
  unique()

cat(sprintf("   > Attempting direct fuzzy match on %d unmatched trait names...\n", 
            length(potential_fuzzy_candidates)))

matches_index <- amatch(potential_fuzzy_candidates, tree_tips_cleaned, maxDist = 2, nomatch = NA)

fuzzy_lookup <- data.frame(
  Trait_Name_Clean = potential_fuzzy_candidates,
  Fuzzy_Tree_Tip   = tree_tips_cleaned[matches_index],
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Fuzzy_Tree_Tip))

# ------------------------------------------------------------------------------
# B2. SYNONYM + FUZZY MATCHING (IUCN synonym → fuzzy → tree tip)
#
# This catches double-difference cases where a species needs BOTH a synonym
# lookup AND a spelling correction to land on a valid tree tip.
# e.g. Trait: "Ablepharus sikimmensis"
#      → IUCN synonym: "Asymblepharus sikkimensis"  (name change)
#      → Tree tip:     "Asymblepharus_sikimmensis"  (spelling differs by 1 char)
# ------------------------------------------------------------------------------
synonym_fuzzy_candidates <- master_name_dictionary %>%
  filter(!Search_Name %in% tree_tips_cleaned) %>%   # exclude those already caught by exact synonym match
  pull(Search_Name) %>%
  unique()

cat(sprintf("   > Attempting synonym+fuzzy match on %d unmatched synonym names...\n",
            length(synonym_fuzzy_candidates)))

matches_syn_fuzzy_idx <- amatch(synonym_fuzzy_candidates, tree_tips_cleaned, maxDist = 2, nomatch = NA)

synonym_fuzzy_lookup <- data.frame(
  Search_Name        = synonym_fuzzy_candidates,
  Syn_Fuzzy_Tree_Tip = tree_tips_cleaned[matches_syn_fuzzy_idx],
  stringsAsFactors   = FALSE
) %>%
  filter(!is.na(Syn_Fuzzy_Tree_Tip))

# Link back to IUCN names: IUCN_Name → synonym → fuzzy tree tip
iucn_to_tree_link_fuzzy <- master_name_dictionary %>%
  inner_join(synonym_fuzzy_lookup, by = "Search_Name") %>%
  group_by(IUCN_Name) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::select(IUCN_Name, Syn_Fuzzy_Tree_Tip)

cat(sprintf("   > Synonym+Fuzzy links resolved for %d IUCN names.\n", nrow(iucn_to_tree_link_fuzzy)))

# ------------------------------------------------------------------------------
# C. MAP TRAIT DATA TO TREE (The Hierarchy)
# ------------------------------------------------------------------------------
trait_mapped <- trait_prep %>%
  mutate(Direct_Match_Exists = Trait_Name_Clean %in% tree_tips_cleaned) %>%
  left_join(master_name_dictionary, by = c("Trait_Name_Clean" = "Search_Name")) %>%
  left_join(iucn_to_tree_link,       by = "IUCN_Name") %>%
  left_join(iucn_to_tree_link_fuzzy, by = "IUCN_Name") %>%   # <-- NEW
  left_join(fuzzy_lookup,            by = "Trait_Name_Clean") %>%
  mutate(
    Final_Tree_Tip = case_when(
      # PRIORITY 1: Direct match (exact spelling)
      Direct_Match_Exists        ~ Trait_Name_Clean,
      # PRIORITY 2: Synonym match (biologically verified via IUCN, exact)
      !is.na(Valid_Tree_Tip)     ~ Valid_Tree_Tip,
      # PRIORITY 3: Synonym + Fuzzy (IUCN synonym + spelling correction)  <-- NEW
      !is.na(Syn_Fuzzy_Tree_Tip) ~ Syn_Fuzzy_Tree_Tip,
      # PRIORITY 4: Direct fuzzy (same name, spelling correction only)
      !is.na(Fuzzy_Tree_Tip)     ~ Fuzzy_Tree_Tip,
      TRUE ~ NA_character_
    ),
    Match_Type = case_when(
      Direct_Match_Exists        ~ "Direct",
      !is.na(Valid_Tree_Tip)     ~ "Synonym",
      !is.na(Syn_Fuzzy_Tree_Tip) ~ "Synonym_Fuzzy",            # <-- NEW
      !is.na(Fuzzy_Tree_Tip)     ~ "Fuzzy",
      TRUE ~ "None"
    )
  )

# ------------------------------------------------------------------------------
# D. DEDUPLICATION & CLEANING
# ------------------------------------------------------------------------------
# Deduplicate by Final_Tree_Tip, coalescing across duplicate rows to minimise
# NAs rather than simply dropping rows. The arrange() ensures we preferentially
# draw non-NA values from the best-quality match (Direct > Synonym > Fuzzy).
dupes_before <- nrow(trait_mapped %>% filter(!is.na(Final_Tree_Tip)))

final_trait_data <- trait_mapped %>%
  filter(!is.na(Final_Tree_Tip)) %>%
  group_by(Final_Tree_Tip) %>%
  arrange(factor(Match_Type, levels = c("Direct", "Synonym", "Synonym_Fuzzy", "Fuzzy")), Trait_Name_Clean) %>%
  summarise(
    across(everything(), ~ {
      non_na <- .x[!is.na(.x)]
      if (length(non_na) == 0) NA else non_na[1]
    }),
    .groups = "drop"
  )

n_dupes_removed <- dupes_before - nrow(final_trait_data)

cat("\n--- Matching Results ---\n")
print(table(final_trait_data$Match_Type))

if (n_dupes_removed > 0) {
  cat(sprintf("   > Merged %d duplicate Final_Tree_Tip entries, coalescing non-NA values.\n", n_dupes_removed))
} else {
  cat("   > No duplicate tree tips found.\n")
}

# ==============================================================================
# 7. UNMATCHED SPECIES REPORT (True Orphans)
# ==============================================================================
cat("\n--- Step 7: Analyzing Unmatched Species ---\n")

# Identify species that failed ALL matching attempts (Direct, Synonym, AND Fuzzy)
unmatched_species <- trait_mapped %>%
  filter(is.na(Final_Tree_Tip)) %>%
  dplyr::select(Species, IUCN_Name) %>%
  distinct()

n_unmatched <- nrow(unmatched_species)

cat(sprintf("Total True Orphans (No tree tip found): %d\n", n_unmatched))

if(n_unmatched > 0) {
  cat("\nExamples of species with NO phylogenetic match:\n")
  print(head(unmatched_species$Species, 50))
  
  # Diagnostic: Check if these missing species are in the IUCN list
  # If TRUE: The species is known to IUCN but the specific name used in the tree is missing or very different.
  # If FALSE: The species might be spelled totally wrong, or is Data Deficient/Not Evaluated in your IUCN set.
  cat("\nOf these unmatched species, how many exist in your IUCN list?\n")
  n_in_iucn <- sum(!is.na(unmatched_species$IUCN_Name))
  cat(sprintf("%d are in IUCN list but still didn't link to the tree.\n", n_in_iucn))
}

# Final Clean Dataset
final_data <- final_trait_data %>%
  dplyr::select(Species, Final_Tree_Tip, Match_Type, everything())

# ==============================================================================
# 8. FILTER TRAIT DATA TO MATCH RANGE AVAILABILITY
# ==============================================================================
cat("--- Step 8: Filtering by Range Data Availability ---\n")

current_count <- nrow(final_data)

# Check which of our current species actually have range data
species_with_ranges <- final_data$IUCN_Name[final_data$IUCN_Name %in% range_names]
missing_range_data  <- setdiff(final_data$IUCN_Name, range_names)

# Filter: Overwrite final_data to only include those with ranges
final_data <- final_data %>%
  filter(IUCN_Name %in% range_names)

cat(sprintf("   > Filter complete.\n"))
if (length(missing_range_data) > 0) {
  cat(sprintf("   > Removed %d species that were missing range data (e.g., %s...).\n", 
              length(missing_range_data), head(missing_range_data, 1)))
} else {
  cat("   > No species removed; all had range data.\n")
}

# Merge duplicate IUCN_Name rows (each range polygon must be unique).
# Rather than dropping rows, we coalesce across duplicates to minimise NAs.
# The arrange() ensures we preferentially draw non-NA values from the 
# best-quality match (Direct > Synonym > Fuzzy) when multiple are available.
dupes_before <- nrow(final_data)

final_data <- final_data %>%
  group_by(IUCN_Name) %>%
  arrange(factor(Match_Type, levels = c("Direct", "Synonym", "Fuzzy")), Species) %>%
  summarise(
    across(everything(), ~ {
      non_na <- .x[!is.na(.x)]
      if (length(non_na) == 0) NA else non_na[1]
    }),
    .groups = "drop"
  )

n_dupes_removed <- dupes_before - nrow(final_data)

if (n_dupes_removed > 0) {
  cat(sprintf("   > Merged %d duplicate IUCN_Name entries, coalescing non-NA values.\n", n_dupes_removed))
} else {
  cat("   > No duplicate IUCN_Names found.\n")
}

final_count <- nrow(final_data)

cat(sprintf("   > FINAL DATASET READY: %d species.\n", final_count))
cat("     (These species exist in IUCN Red List, Trait Database, and Range Polygons)\n")

# Prune Tree
tips_to_keep <- final_data$Final_Tree_Tip
pruned_tree <- keep.tip(tree, tips_to_keep)

# Sort Data
final_data_sorted <- final_data %>%
  arrange(match(Final_Tree_Tip, pruned_tree$tip.label))

# Show breakdown of how matches were found
cat("\n--- Match Methods Used ---\n")
table(final_data_sorted$Match_Type)

# ==============================================================================
# 9. ADD THREAT DATA (LOGGING & PLANTATION)
# ==============================================================================
cat("--- Step 9: Attaching Threat Columns (Logging & Plantation) ---\n")

# We join the final dataset with the specific threat columns from iucn_clean.
# We use 'IUCN_Name' as the key because that is the standardized name 
# shared by both dataframes.

final_data_with_threats <- final_data_sorted %>%
  left_join(
    iucn_clean %>% dplyr::select(IUCN_Name, Threat_Logging, Threat_Plantation),
    by = "IUCN_Name"
  )

# Quick check to ensure columns were added
if("Threat_Logging" %in% names(final_data_with_threats)) {
  cat(sprintf("   > Threat columns successfully added.\n"))
  cat(sprintf("   > Species threatened by logging:    %d\n", sum(final_data_with_threats$Threat_Logging, na.rm=TRUE)))
  cat(sprintf("   > Species threatened by plantation: %d\n", sum(final_data_with_threats$Threat_Plantation, na.rm=TRUE)))
} else {
  warning("Threat columns FAILED to attach!")
}

# Update the main variable name for the final save step
final_data_sorted <- final_data_with_threats

# ==============================================================================
# 10. FINALIZE & SAVE
# ==============================================================================
cat("--- Step 11: Finalizing Output ---\n")

# A. Save the Data
write.csv(final_data_sorted, "data/processed/traits/reptile/reptile_trait_data_matched.csv", row.names = FALSE)

# B. Save the Tree
write.tree(pruned_tree, "data/processed/traits/reptile/final_reptile_tree_pruned.tre")

cat("--- SCRIPT COMPLETE ---\n")

