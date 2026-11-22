library(tidyverse)
library(rredlist)
library(sf)
library(ape)
library(here)
library(stringr)
library(phangorn)

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
    "data", "processed", "assessments", "mammal_assessment_taxonomy_lookup.rds"
  )
) %>%
  filter(!RedList_Cat %in% c("EX", "EW")) %>%
  mutate(IUCN_Name = clean_names(IUCN_Name))

cat(sprintf("   > IUCN Backbone Ready: %d species.\n", nrow(iucn_clean)))

# ==============================================================================
# 2. LOAD TRAIT DATA (COMBINE)
# ==============================================================================
cat("--- Step 2: Loading Trait Data ---\n")
trait_df <- read.csv("data/raw/traits/mammal/combine_mammal_traits.csv", stringsAsFactors = FALSE)

# ==============================================================================
# 3. LOAD RANGE POLYGONS
# ==============================================================================
cat("--- Step 3: Loading Range Polygons ---\n")
# We only need the list of names for now to check availability
if(file.exists(here("data/processed/ranges/all_mammal_ranges_dissolved.rds"))) {
  ranges <- readRDS(here("data/processed/ranges/all_mammal_ranges_dissolved.rds"))
  # If sf object, drop geometry to speed up 'names check'
  range_names <- ranges$sci_name # Ensure column name matches
  rm(ranges) # Free memory
  cat("   > Range data found for", length(range_names), "species.\n")
} else {
  stop("Range file missing!")
}

# ==============================================================================
# 4. LOAD PHYLOGENY (Upham 2019)
# ==============================================================================
cat("--- Step 4: Fetching Upham Tree ---\n")
# 1. Read trees
my_trees <- read.nexus("data/raw/trees/mammaltree/mammaltree.nex")

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
# 5. DIAGNOSE & FIX DUPLICATES (STRICT "INPUT_NAME" UNIQUE)
# ==============================================================================
cat("--- Step 5: Bridging Taxonomy & Merging Data ---\n")

# 1. Load Data
crosswalk <- read.csv("data/raw/traits/mammal/taxonomy_crosswalk.csv", stringsAsFactors = FALSE)
# Ensure trait_df is loaded from Step 2

# 2. Prepare Crosswalk + Traits
# We merge them first, THEN deduplicate based on the column that joins to the Red List.
crosswalk_clean <- crosswalk %>%
  mutate(input_name = clean_names(input_name)) %>%
  select(iucn2020_binomial, input_name) %>%
  distinct()

# Merge traits to crosswalk
# This creates a "Master Trait Lookup" table
trait_lookup <- trait_df %>%
  inner_join(crosswalk_clean, by = "iucn2020_binomial")

# 3. IDENTIFY & REMOVE DUPLICATES
# This is the critical fix. We check for duplicates in 'input_name'.
# If 'input_name' is duplicated, it will break the Red List join.

# Optional: View the culprits (for your own knowledge)
duplicates <- trait_lookup %>%
  add_count(input_name) %>%
  filter(n > 1) %>%
  arrange(input_name)
 View(duplicates) # Uncomment to see exactly which 191 species are causing issues

# FORCE UNIQUE on 'input_name'
trait_lookup_unique <- trait_lookup %>%
  # We group by input_name and take the first entry.
  # This guarantees 1 row per input_name.
  group_by(input_name) %>%
  slice(1) %>%
  ungroup()

cat(sprintf("   > Trait Lookup Size: %d rows\n", nrow(trait_lookup_unique)))

# 4. Join to the IUCN Backbone
final_data <- iucn_clean %>%
  left_join(trait_lookup_unique, by = c("IUCN_Name" = "input_name"))

# ==============================================================================
# 6. VERIFICATION
# ==============================================================================

n_backbone <- nrow(iucn_clean)
n_final <- nrow(final_data)

cat("--- Quality Control Checks ---\n")
cat(sprintf("1. Original IUCN Backbone count: %d\n", n_backbone))
cat(sprintf("2. Final Data count:             %d\n", n_final))

if (n_backbone == n_final) {
  cat("   [SUCCESS] Row counts match! No duplicates created.\n")
} else {
  cat(sprintf("   [WARNING] Row counts differ by %d rows. Duplicates persist.\n", n_final - n_backbone))
}

# Check Match Rate
matched_count <- sum(!is.na(final_data$iucn2020_binomial))
cat(sprintf("3. Species successfully linked to traits: %d (%.1f%%)\n", 
            matched_count, (matched_count / n_backbone) * 100))

# ==============================================================================
# 7. FILTER TRAIT DATA TO MATCH RANGE AVAILABILITY
# ==============================================================================
cat("--- Step 6: Filtering by Range Data Availability ---\n")

current_count <- nrow(final_data)

# Check which of our current species actually have range data
species_with_ranges <- final_data$IUCN_Name[final_data$IUCN_Name %in% range_names]
missing_range_data  <- setdiff(final_data$IUCN_Name, range_names)

# Filter: Overwrite final_data to only include those with ranges
final_data <- final_data %>%
  filter(IUCN_Name %in% range_names)

final_count <- nrow(final_data)

cat(sprintf("   > Filter complete.\n"))
if (length(missing_range_data) > 0) {
  cat(sprintf("   > Removed %d species that were missing range data (e.g., %s...).\n", 
              length(missing_range_data), head(missing_range_data, 1)))
} else {
  cat("   > No species removed; all had range data.\n")
}

cat(sprintf("   > FINAL DATASET READY: %d species.\n", final_count))
cat("     (These species exist in IUCN Red List, Trait Database, and Range Polygons)\n")

# ==============================================================================
# 8. MATCH DATA TO PHYLOGENY (EXACT + SYNONYM + FUZZY)
# ==============================================================================
cat("--- Step 7: Synchronizing Data and Phylogeny (Maximized) ---\n")

# 1. Clean Tree Tips
tree_tips_clean <- tree$tip.label %>%
  str_replace_all("_", " ") %>%
  str_trim()

# Update tree immediately
tree$tip.label <- tree_tips_clean

# ------------------------------------------------------------------------------
# STAGE 1: EXACT & CROSSWALK MATCHING
# ------------------------------------------------------------------------------

# Create a lookup table of ALL synonyms in the crosswalk that appear in the tree
# This bridges: Data (iucn2020) -> Crosswalk (input_name) -> Tree Tip
valid_synonyms <- crosswalk %>%
  mutate(input_name = clean_names(input_name)) %>%
  # Only look at synonyms for species we actually have data for
  filter(iucn2020_binomial %in% final_data$iucn2020_binomial) %>%
  # Only keep synonyms that ARE in the tree
  filter(input_name %in% tree_tips_clean) %>%
  select(iucn2020_binomial, matched_tip = input_name) %>%
  distinct()

# Apply Exact and Synonym Matches
data_stage1 <- final_data %>%
  left_join(valid_synonyms, by = "iucn2020_binomial") %>%
  mutate(
    phylo_match = case_when(
      # 1. Try Standard Name
      iucn2020_binomial %in% tree_tips_clean ~ iucn2020_binomial,
      # 2. Try IUCN Input Name
      IUCN_Name %in% tree_tips_clean ~ IUCN_Name,
      # 3. Try Crosswalk Synonym
      !is.na(matched_tip) ~ matched_tip,
      TRUE ~ NA_character_
    ),
    match_method = case_when(
      !is.na(phylo_match) ~ "Exact/Synonym",
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------------------------
# STAGE 2: FUZZY MATCHING (Recovering Spelling/Suffix diffs)
# ------------------------------------------------------------------------------
# We only try to fuzzy match species that FAILED Stage 1
unmatched_species <- data_stage1 %>% filter(is.na(phylo_match))
matched_species   <- data_stage1 %>% filter(!is.na(phylo_match))

if(nrow(unmatched_species) > 0){
  cat(sprintf("   > Attempting fuzzy match for %d unmatched species...\n", nrow(unmatched_species)))
  
  # Calculate string distances between unmatched IUCN names and unused Tree Tips
  # We use 'agrep' logic or stringdist. Here is a simplified base R approach.
  # We allow a max distance of 1 or 2 characters (conservative).
  
  # Identify tips in tree that are NOT yet used
  used_tips <- unique(matched_species$phylo_match)
  available_tips <- setdiff(tree_tips_clean, used_tips)
  
  # Function to find closest match
  find_fuzzy <- function(name, candidates) {
    # aginst: approximate grep integer
    # max.distance = 0.1 means ~1-2 chars tolerance for typical species names
    matches <- agrep(name, candidates, max.distance = 0.05, value = TRUE)
    if(length(matches) == 1) return(matches) # Only accept unique close matches
    return(NA_character_)
  }
  
  # Apply fuzzy matching (progress bar helps if list is long)
  # This might take 30-60 seconds
  fuzzy_results <- unmatched_species %>%
    rowwise() %>%
    mutate(
      fuzzy_tip = find_fuzzy(iucn2020_binomial, available_tips)
    ) %>%
    ungroup()
  
  # Merge Fuzzy Results back
  data_stage2 <- matched_species %>%
    bind_rows(
      fuzzy_results %>%
        mutate(
          phylo_match = fuzzy_tip,
          match_method = ifelse(!is.na(fuzzy_tip), "Fuzzy", NA_character_)
        ) %>%
        select(-fuzzy_tip)
    )
} else {
  data_stage2 <- data_stage1
}

# ------------------------------------------------------------------------------
# STAGE 3: CLEANUP & PRUNING
# ------------------------------------------------------------------------------

# Filter down to matched data
final_data_phylo <- data_stage2 %>%
  filter(!is.na(phylo_match))

# Resolve Duplicates (Many-to-One)
# If multiple data rows match the SAME tree tip, keep the one with the most data 
# or just the first one.
final_data_unique <- final_data_phylo %>%
  add_count(phylo_match) %>%
  arrange(phylo_match, desc(match_method)) %>% # Prefer Exact over Fuzzy
  group_by(phylo_match) %>%
  slice(1) %>%
  ungroup()

# Prune Tree
tips_to_keep <- final_data_unique$phylo_match
pruned_tree <- keep.tip(tree, tips_to_keep)

# Sort Data
final_data_sorted <- final_data_unique %>%
  arrange(match(phylo_match, pruned_tree$tip.label))

# ==============================================================================
# 9. REPORTING
# ==============================================================================

n_total_input <- nrow(final_data)
n_matched     <- nrow(final_data_sorted)
n_lost        <- n_total_input - n_matched

cat("--- Final Match Report ---\n")
cat(sprintf("   > Input Species:   %d\n", n_total_input))
cat(sprintf("   > Matched Species: %d (%.1f%%)\n", n_matched, (n_matched/n_total_input)*100))
cat(sprintf("   > Lost Species:    %d\n", n_lost))

# Show breakdown of how matches were found
cat("\n--- Match Methods Used ---\n")
table(final_data_sorted$match_method)

# SAVE THE MISSING ONES FOR INSPECTION
missing_species <- data_stage2 %>%
  filter(is.na(phylo_match)) %>%
  select(IUCN_Name, iucn2020_binomial)

if(nrow(missing_species) > 0) {
write.csv(missing_species, "data/processed/traits/mammal/missing_species_report.csv", row.names = FALSE)
cat("\n   > [ACTION] 'data/processed/traits/mammal/missing_species_report.csv' has been saved.\n")
  cat("     Open this file to see which species are missing. \n")
  cat("     Most are likely recent splits (e.g., 'Bovidae' splits) absent from the 2019 tree.\n")
}

# Final Check
if (all(final_data_sorted$phylo_match == pruned_tree$tip.label)) {
  cat("\n   [SUCCESS] Data and Tree are perfectly aligned.\n")
} else {
  stop("   [ERROR] Mismatch in sorting.")
}

# ==============================================================================
# DIAGNOSING THE MISSING SPECIES
# ==============================================================================

# 1. Identify the 'Family' column name dynamically
# We look for a column in trait_df that contains "family" (case-insensitive)
fam_col <- grep("family", names(trait_df), ignore.case = TRUE, value = TRUE)[1]

if (is.na(fam_col)) {
  stop("   [ERROR] Could not find a column named 'Family' or 'family' in trait_df. Please check names(trait_df).")
} else {
  cat(sprintf("   > Using column '%s' for taxonomic grouping.\n", fam_col))
}

# 2. Run the analysis using the detected column
missing_analysis <- data_stage2 %>%
  filter(is.na(phylo_match)) %>%
  # Join to trait_df to get the Family info
  # We use !!sym(fam_col) to tell R to use the string we found above
  left_join(trait_df %>% select(iucn2020_binomial, family_name = !!sym(fam_col)), 
            by = "iucn2020_binomial") %>%
  # Count by the new standardized 'family_name' column
  count(family_name, sort = TRUE)

# 3. View the results
print(head(missing_analysis, 15))

# Optional: Check if we are losing any "important" groups (e.g., Primates, Carnivora)
# Adjust the filter list based on your study focus
cat("\n--- Losses in Key Groups ---\n")
key_groups <- c("Felidae", "Canidae", "Hominidae", "Cercopithecidae")
print(missing_analysis %>% filter(family_name %in% key_groups))

# ==============================================================================
# 10. FINALIZE & SAVE
# ==============================================================================
cat("--- Step 8: Finalizing Output ---\n")

# 1. Clean up columns
# We can drop the temporary matching columns now that sorting is done
final_analysis_data <- final_data_sorted %>%
  select(
    IUCN_Name,              # Your main ID
    iucn2020_binomial,      # The scientific name link
    RedList_Cat,            # Response Variable
    Threat_Logging,         # Predictor 1
    Threat_Plantation,      # Predictor 2
    # Add your specific trait columns here, e.g.:
    # Mass_g,
    # Diet_Type,
    # Family,
    everything()            # Keep the rest just in case
  )
  
# A. Save the Data
write.csv(final_analysis_data, "data/processed/traits/mammal/mammal_trait_data_matched.csv", row.names = FALSE)
  
# B. Save the Tree
# We use write.tree (Newick format) or write.nexus
write.tree(pruned_tree, "data/processed/traits/mammal/final_mammal_tree_pruned.tre")

