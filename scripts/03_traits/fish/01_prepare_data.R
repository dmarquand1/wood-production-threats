# =============================================================================
# Freshwater Fish Trait × Phylogeny Pipeline
#
# Purpose : Build a matched species × trait dataset and a pruned phylogeny for
#           freshwater Actinopterygii that have (a) an IUCN assessment,
#           (b) range polygon data, and (c) a place on the Fish Tree of Life.
#           Only traits with ≥ 60 % species coverage are retained.
#
# Outputs
#   data/processed/traits/fish/fish_high_coverage_traits.csv – wide trait matrix
#   data/processed/traits/fish/fish_trait_data_matched.csv   – matched + threat-annotated
#   outputs/traits/diagnostics/field_coverage_all_tables.csv – per-field coverage report
#   outputs/traits/diagnostics/table_coverage_summary.csv    – per-table coverage summary
#   data/processed/traits/fish/final_fish_tree_pruned.tre    – pruned MCC tree (Newick)
#   data/processed/traits/fish/pruned_trees_100.rds          – all 100 pruned trees
#
# Intermediate caches
#   data/raw/iucn/fish_iucn_assessment_list.rds
#   data/raw/traits/fish/fb_name_to_speccode.rds
#   data/raw/traits/fish/raw_rfishbase_tables.rds
# =============================================================================

rm(list = ls())
gc()

# =============================================================================
# 0. LIBRARIES & CONFIGURATION
# =============================================================================

library(rfishbase)
library(fishtree)
library(rredlist)
library(ape)
library(phangorn)
library(stringdist)
library(tidyverse)
library(sf)
library(here)

# ── Output directories ────────────────────────────────────────────────────────
dir.create(here("outputs/traits"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("data", "processed", "traits", "fish"), recursive = TRUE, showWarnings = FALSE)

# ── Threat codes ──────────────────────────────────────────────────────────────
THREAT_CODE_LOGGING    <- "5_3_"
THREAT_CODE_PLANTATION <- "2_2_"

# ── Coverage threshold ────────────────────────────────────────────────────────
COVERAGE_THRESHOLD <- 60   # percent

# ── Identifier columns to exclude from trait analysis ─────────────────────────
IDENTIFIER_COLS <- c(
  "SpecCode", "Species", "sciname", "Genus", "Species2",
  "StockCode", "AutoCtr", "TSN", "GenCode", "SubGenCode", "BodyShapeI"
)

# =============================================================================
# 1. HELPER FUNCTIONS
# =============================================================================

# Standardise Unicode/whitespace in species name strings
clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  x |>
    str_replace_all("\uFFFD", " ") |>
    str_replace_all("\u00A0", " ") |>
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") |>
    str_squish() |>
    str_trim()
}

# Safely call an rfishbase function; return NULL on error
safe_fetch <- function(fn, spp_vec, ...) {
  tryCatch(
    fn(spp_vec, ...),
    error = function(e) {
      message("  ERROR in ", deparse(substitute(fn)), ": ", conditionMessage(e))
      NULL
    }
  )
}

# Compute per-field species coverage for one rfishbase table
compute_coverage <- function(df, table_label, spp_vec) {
  empty <- tibble(table = table_label, field = NA_character_,
                  n_species = 0L, pct_coverage = 0)
  if (is.null(df) || nrow(df) == 0) return(empty)
  
  trait_cols <- setdiff(names(df), IDENTIFIER_COLS)
  trait_cols  <- trait_cols[map_lgl(df[trait_cols], ~ !all(is.na(.x)))]
  if (length(trait_cols) == 0) return(empty)
  
  sp_col <- intersect(c("Species", "sciname", "Species2"), names(df))[1]
  
  if (is.na(sp_col)) {
    map_dfr(trait_cols, function(col) {
      n <- sum(!is.na(df[[col]]))
      tibble(table = table_label, field = col, n_species = n,
             pct_coverage = round(100 * n / nrow(df), 1))
    })
  } else {
    map_dfr(trait_cols, function(col) {
      n_spp <- df |>
        filter(!is.na(.data[[col]])) |>
        pull(.data[[sp_col]]) |>
        unique() |>
        intersect(spp_vec) |>
        length()
      tibble(table = table_label, field = col, n_species = n_spp,
             pct_coverage = round(100 * n_spp / length(spp_vec), 1))
    })
  }
}

# Collapse a per-species table to one row per species:
#   numeric cols → mean; character cols → first non-NA
summarise_to_species <- function(df, spp_vec) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  sc <- intersect(c("Species", "sciname"), names(df))[1]
  if (is.na(sc)) return(NULL)
  
  df <- filter(df, .data[[sc]] %in% spp_vec)
  if (nrow(df) == 0) return(NULL)
  
  trait_cols <- setdiff(names(df), IDENTIFIER_COLS)
  num_cols   <- trait_cols[map_lgl(df[trait_cols], is.numeric)]
  chr_cols   <- trait_cols[!map_lgl(df[trait_cols], is.numeric)]
  
  df |>
    group_by(Species = .data[[sc]]) |>
    summarise(
      across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)),
      across(all_of(chr_cols), ~ na.omit(.x)[1]),
      .groups = "drop"
    )
}

# =============================================================================
# 2. IUCN RED LIST DATA
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 1 — Loading prepared IUCN assessment backbone\n")
cat(strrep("=", 70), "\n")

iucn_clean <- readRDS(
  here("data", "processed", "assessments", "fish_assessment_taxonomy_lookup.rds")
) |>
  filter(!RedList_Cat %in% c("EX", "EW")) |>
  mutate(
    IUCN_Name = clean_names(IUCN_Name),
    Synonyms = purrr::map(Synonyms, clean_names)
  )

cat(sprintf("IUCN backbone ready: %d freshwater extant species.\n", nrow(iucn_clean)))

# =============================================================================
# 3. VALIDATE SPECIES NAMES AGAINST FISHBASE
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 2 — Validating Species Names Against FishBase\n")
cat(strrep("=", 70), "\n")

# Use IUCN species names as the working species list
freshwater_spp_raw <- iucn_clean$IUCN_Name
validated_names    <- validate_names(freshwater_spp_raw)

unmatched_fb <- freshwater_spp_raw[!freshwater_spp_raw %in% validated_names]
if (length(unmatched_fb) > 0) {
  message(length(unmatched_fb),
          " name(s) not found in FishBase and will be excluded from trait fetch.")
}

spp <- validated_names[!is.na(validated_names)]
cat(sprintf("%d species names validated and retained for trait fetching.\n", length(spp)))

# =============================================================================
# 4. FETCH RFISHBASE TRAIT TABLES
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 3 — Fetching Trait Tables From FishBase\n")
cat(strrep("=", 70), "\n")

fb_tables_cache <- here("data/raw/traits/fish/raw_rfishbase_tables.rds")

if (file.exists(fb_tables_cache)) {
  cat("Loading rfishbase tables from local cache...\n")
  table_df_lookup <- readRDS(fb_tables_cache)
} else {
  cat("Fetching all trait tables from FishBase (this may take several minutes)...\n")
  
  fetch_msg <- function(label, fn, ...) {
    message("  Fetching: ", label, " ...")
    safe_fetch(fn, spp, ...)
  }
  
  table_df_lookup <- list(
    species       = fetch_msg("species",       species),
    ecology       = fetch_msg("ecology",       ecology),
    morphology    = fetch_msg("morphology",    morphology),
    morphometrics = fetch_msg("morphometrics", morphometrics),
    diet          = fetch_msg("diet",          diet),
    fooditems     = fetch_msg("fooditems",     fooditems),
    estimate      = fetch_msg("estimate",      estimate),
    popgrowth     = fetch_msg("popgrowth",     popgrowth),
    popchar       = fetch_msg("popchar",       popchar),
    length_weight = fetch_msg("length_weight", length_weight),
    length_length = fetch_msg("length_length", length_length),
    length_freq   = fetch_msg("length_freq",   length_freq),
    maturity      = fetch_msg("maturity",      maturity),
    reproduction  = fetch_msg("reproduction",  reproduction),
    spawning      = fetch_msg("spawning",      spawning),
    fecundity     = fetch_msg("fecundity",     fecundity),
    stocks        = fetch_msg("stocks",        stocks),
    oxygen        = fetch_msg("oxygen",        oxygen),
    swimming      = fetch_msg("swimming",      swimming),
    speed         = fetch_msg("speed",         speed),
    larvae        = fetch_msg("larvae",        larvae),
    genetics      = fetch_msg("genetics",      genetics),
    introductions = fetch_msg("introductions", introductions),
    predators     = fetch_msg("predators",     predators),
    brains        = fetch_msg("brains",        brains),
    country       = fetch_msg("country",       country),
    faoareas      = fetch_msg("faoareas",      faoareas),
    ecosystem     = fetch_msg("ecosystem",     ecosystem),
    popqb         = fetch_msg("popqb",         popqb),
    ration        = fetch_msg("ration",        ration)
  )
  
  saveRDS(table_df_lookup, fb_tables_cache)
  cat("rfishbase tables cached to", fb_tables_cache, "\n")
}

cat("All tables loaded.\n")

# =============================================================================
# 5. COMPUTE COVERAGE & BUILD HIGH-COVERAGE TRAIT MATRIX
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 4 — Computing Coverage & Building Trait Matrix\n")
cat(strrep("=", 70), "\n")

field_coverage <- imap_dfr(table_df_lookup, ~ compute_coverage(.x, .y, spp))

table_coverage_summary <- field_coverage |>
  filter(!is.na(field)) |>
  group_by(table) |>
  summarise(
    n_fields            = n(),
    max_spp_any_field   = max(n_species),
    max_pct_any_field   = max(pct_coverage),
    mean_pct_all_fields = round(mean(pct_coverage), 1),
    n_fields_gt10pct    = sum(pct_coverage > 10),
    n_fields_gt50pct    = sum(pct_coverage > 50),
    .groups = "drop"
  ) |>
  arrange(desc(max_pct_any_field))

cat("\nTable-level coverage summary (top tables by max field coverage):\n")
print(table_coverage_summary, n = Inf)

write.csv(field_coverage,
          here("outputs/traits/diagnostics/field_coverage_all_tables.csv"),
          row.names = FALSE)
write.csv(table_coverage_summary,
          here("outputs/traits/diagnostics/table_coverage_summary.csv"),
          row.names = FALSE)

# Identify field–table combinations at or above the threshold
high_cov_fields <- field_coverage |>
  filter(pct_coverage >= COVERAGE_THRESHOLD, !is.na(field)) |>
  arrange(desc(pct_coverage))

cat(sprintf("\n%d field-table combinations with >= %d%% coverage retained.\n",
            nrow(high_cov_fields), COVERAGE_THRESHOLD))

# Summarise each needed table to one row per species, then join
needed_tables     <- unique(high_cov_fields$table)
summarised_tables <- map(
  setNames(needed_tables, needed_tables),
  ~ summarise_to_species(table_df_lookup[[.x]], spp)
)

trait_values_wide <- tibble(Species = spp)

for (tbl in needed_tables) {
  fields_this_table <- high_cov_fields$field[high_cov_fields$table == tbl]
  summ              <- summarised_tables[[tbl]]
  if (is.null(summ)) next
  
  cols_present <- intersect(fields_this_table, names(summ))
  if (length(cols_present) == 0) next
  
  summ_prefixed <- summ |>
    select(Species, all_of(cols_present)) |>
    rename_with(~ paste(tbl, .x, sep = "_"), all_of(cols_present))
  
  trait_values_wide <- left_join(trait_values_wide, summ_prefixed, by = "Species")
}

trait_values_wide <- mutate(trait_values_wide,
                            across(where(is.numeric), ~ ifelse(is.na(.x), NaN, .x)))

write.csv(trait_values_wide,
          here("data/processed/traits/fish/fish_high_coverage_traits.csv"),
          row.names = FALSE, na = "NaN")

cat(sprintf("Trait matrix saved: %d species × %d columns (%d trait columns).\n",
            nrow(trait_values_wide), ncol(trait_values_wide), nrow(high_cov_fields)))

# =============================================================================
# 6. LOAD RANGE POLYGON DATA
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 5 — Loading Range Polygon Data\n")
cat(strrep("=", 70), "\n")

range_file <- here("data/processed/ranges/all_fish_ranges_dissolved.rds")
if (!file.exists(range_file)) stop("Range file missing: ", range_file)

ranges      <- readRDS(range_file)
range_names <- ranges$sci_name
rm(ranges); gc()
cat(sprintf("Range data found for %d species.\n", length(range_names)))

# =============================================================================
# 7. FETCH FISH TREE OF LIFE PHYLOGENY
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 6 — Fetching Fish Tree of Life Phylogeny\n")
cat(strrep("=", 70), "\n")

all_trees <- fishtree_complete_phylogeny()
tree      <- maxCladeCred(all_trees[[1]])          # use first tree for name matching only
tree_tips <- tree$tip.label
cat(sprintf("Full tree has %d tips across %d posterior trees.\n",
            length(tree_tips), length(all_trees)))

# =============================================================================
# 8. TAXONOMIC NAME RESOLUTION
#    Priority: Direct → IUCN synonym → FishBase SpecCode bridge → Fuzzy
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 7 — Taxonomic Name Resolution\n")
cat(strrep("=", 70), "\n")

# Working copies with cleaned, underscore-separated names
tree_tips_cleaned <- clean_names(tree_tips) |> str_replace_all(" ", "_")
valid_tips_df     <- data.frame(Tree_Tip = tree_tips_cleaned,
                                stringsAsFactors = FALSE)

trait_prep <- trait_values_wide |>
  mutate(Trait_Name_Clean = clean_names(Species) |> str_replace_all(" ", "_"))

# ── 7a. IUCN synonym dictionary ───────────────────────────────────────────────
synonyms_long <- iucn_clean |>
  select(IUCN_Name, Synonyms) |>
  unnest(Synonyms) |>
  mutate(Search_Name = clean_names(Synonyms) |> str_replace_all(" ", "_")) |>
  select(IUCN_Name, Search_Name)

self_ref <- iucn_clean |>
  mutate(Search_Name = clean_names(IUCN_Name) |> str_replace_all(" ", "_")) |>
  select(IUCN_Name, Search_Name)

master_name_dictionary <- bind_rows(synonyms_long, self_ref) |>
  distinct() |>
  filter(!is.na(Search_Name), Search_Name != "")

iucn_to_tree_link <- master_name_dictionary |>
  rename(Possible_Tree_Tip = Search_Name) |>
  inner_join(valid_tips_df, by = c("Possible_Tree_Tip" = "Tree_Tip")) |>
  group_by(IUCN_Name) |>
  slice(1) |>
  ungroup() |>
  rename(Valid_Tree_Tip = Possible_Tree_Tip) |>
  select(IUCN_Name, Valid_Tree_Tip)

cat(sprintf("IUCN synonym bridge: %d IUCN names linked to tree tips.\n",
            nrow(iucn_to_tree_link)))

# ── 7b. FishBase SpecCode bridge ──────────────────────────────────────────────
fb_cache_path <- here("data/raw/traits/fish/fb_name_to_speccode.rds")

if (file.exists(fb_cache_path)) {
  cat("Loading FishBase SpecCode map from cache...\n")
  fb_name_to_speccode <- readRDS(fb_cache_path)
} else {
  cat("Downloading complete FishBase synonym table (one-time, ~30 MB)...\n")
  fb_syn_full <- rfishbase::synonyms() |>
    dplyr::collect() |>
    filter(!is.na(SpecCode), !is.na(synonym)) |>
    filter(!str_detect(str_trim(synonym), "^NA( NA)?$"))
  
  fb_name_to_speccode <- fb_syn_full |>
    select(synonym, SpecCode) |>
    mutate(name_clean = clean_names(synonym) |> str_replace_all(" ", "_")) |>
    filter(!is.na(name_clean), name_clean != "",
           str_count(name_clean, "_") >= 1) |>
    select(name_clean, SpecCode) |>
    distinct()
  
  cat(sprintf("SpecCode map built: %d unique name entries.\n",
              nrow(fb_name_to_speccode)))
  saveRDS(fb_name_to_speccode, fb_cache_path)
}

# Map tree tips → SpecCodes
fb_tree_codes <- valid_tips_df |>
  left_join(fb_name_to_speccode, by = c("Tree_Tip" = "name_clean")) |>
  filter(!is.na(SpecCode)) |>
  rename(tree_name_clean = Tree_Tip) |>
  select(tree_name_clean, SpecCode) |>
  distinct()

cat(sprintf("SpecCodes found for %d / %d tree tips.\n",
            n_distinct(fb_tree_codes$tree_name_clean), length(tree_tips_cleaned)))

# Map unresolved trait & IUCN names → SpecCodes
unmatched_trait_clean <- trait_prep |>
  filter(!Trait_Name_Clean %in% tree_tips_cleaned) |>
  pull(Trait_Name_Clean) |>
  unique()

unmatched_iucn_clean <- iucn_clean |>
  mutate(name_clean = clean_names(IUCN_Name) |> str_replace_all(" ", "_")) |>
  filter(!name_clean %in% tree_tips_cleaned) |>
  pull(name_clean) |>
  unique()

query_names_clean <- unique(c(unmatched_trait_clean, unmatched_iucn_clean))

fb_query_codes <- tibble(query_name_clean = query_names_clean) |>
  left_join(fb_name_to_speccode, by = c("query_name_clean" = "name_clean")) |>
  filter(!is.na(SpecCode)) |>
  select(query_name_clean, SpecCode) |>
  distinct()

# Bridge: shared SpecCode → same biological species
fishbase_to_tree_link <- fb_query_codes |>
  inner_join(fb_tree_codes, by = "SpecCode", relationship = "many-to-many") |>
  filter(query_name_clean != tree_name_clean) |>
  group_by(query_name_clean) |>
  slice(1) |>
  ungroup() |>
  rename(Search_Name = query_name_clean, FishBase_Tree_Tip = tree_name_clean)

cat(sprintf("FishBase bridge resolved %d names to tree tips.\n",
            nrow(fishbase_to_tree_link)))

# ── 7c. Fuzzy matching (last resort) ─────────────────────────────────────────
fuzzy_candidates <- trait_prep |>
  filter(
    !Trait_Name_Clean %in% tree_tips_cleaned,
    !Trait_Name_Clean %in% fishbase_to_tree_link$Search_Name
  ) |>
  pull(Trait_Name_Clean) |>
  unique()

cat(sprintf("Attempting fuzzy match on %d still-unresolved names...\n",
            length(fuzzy_candidates)))

match_idx    <- amatch(fuzzy_candidates, tree_tips_cleaned, maxDist = 2, nomatch = NA)
fuzzy_lookup <- data.frame(
  Trait_Name_Clean = fuzzy_candidates,
  Fuzzy_Tree_Tip   = tree_tips_cleaned[match_idx],
  stringsAsFactors = FALSE
) |>
  filter(!is.na(Fuzzy_Tree_Tip))

cat(sprintf("Fuzzy matching resolved %d additional names.\n", nrow(fuzzy_lookup)))

# ── 7d. Assemble final mapping ────────────────────────────────────────────────
match_levels <- c("Direct", "Synonym_IUCN", "Synonym_FishBase",
                  "Synonym_FishBase_via_IUCN", "Fuzzy")

trait_mapped <- trait_prep |>
  mutate(Direct_Match_Exists = Trait_Name_Clean %in% tree_tips_cleaned) |>
  left_join(master_name_dictionary,   by = c("Trait_Name_Clean" = "Search_Name")) |>
  left_join(iucn_to_tree_link,        by = "IUCN_Name") |>
  left_join(fishbase_to_tree_link,    by = c("Trait_Name_Clean" = "Search_Name")) |>
  left_join(
    fishbase_to_tree_link |> rename(FishBase_Tree_Tip_via_IUCN = FishBase_Tree_Tip),
    by = c("IUCN_Name" = "Search_Name")
  ) |>
  left_join(fuzzy_lookup, by = "Trait_Name_Clean") |>
  mutate(
    Final_Tree_Tip = case_when(
      Direct_Match_Exists                ~ Trait_Name_Clean,
      !is.na(Valid_Tree_Tip)             ~ Valid_Tree_Tip,
      !is.na(FishBase_Tree_Tip)          ~ FishBase_Tree_Tip,
      !is.na(FishBase_Tree_Tip_via_IUCN) ~ FishBase_Tree_Tip_via_IUCN,
      !is.na(Fuzzy_Tree_Tip)             ~ Fuzzy_Tree_Tip,
      TRUE                               ~ NA_character_
    ),
    Match_Type = case_when(
      Direct_Match_Exists                ~ "Direct",
      !is.na(Valid_Tree_Tip)             ~ "Synonym_IUCN",
      !is.na(FishBase_Tree_Tip)          ~ "Synonym_FishBase",
      !is.na(FishBase_Tree_Tip_via_IUCN) ~ "Synonym_FishBase_via_IUCN",
      !is.na(Fuzzy_Tree_Tip)             ~ "Fuzzy",
      TRUE                               ~ "None"
    )
  )

# Deduplicate: keep highest-priority match when multiple species map to same tip
dupes_before <- nrow(filter(trait_mapped, !is.na(Final_Tree_Tip)))

final_trait_data <- trait_mapped |>
  filter(!is.na(Final_Tree_Tip)) |>
  group_by(Final_Tree_Tip) |>
  arrange(factor(Match_Type, levels = match_levels), Trait_Name_Clean) |>
  summarise(
    across(everything(), ~ { nna <- .x[!is.na(.x)]; if (length(nna) == 0) NA else nna[1] }),
    .groups = "drop"
  )

n_dupes <- dupes_before - nrow(final_trait_data)

cat("\nMatch type breakdown:\n")
print(table(final_trait_data$Match_Type))
if (n_dupes > 0)
  cat(sprintf("Merged %d duplicate Final_Tree_Tip entries (coalesced non-NA values).\n",
              n_dupes))

# =============================================================================
# 9. UNMATCHED SPECIES REPORT
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 8 — Unmatched Species Report\n")
cat(strrep("=", 70), "\n")

unmatched_species <- trait_mapped |>
  filter(is.na(Final_Tree_Tip)) |>
  select(Species, IUCN_Name) |>
  distinct()

cat(sprintf("True orphans (no tree tip found): %d\n", nrow(unmatched_species)))

if (nrow(unmatched_species) > 0) {
  cat("First 10 unmatched:\n")
  print(head(unmatched_species$Species, 10))
  
  orphan_names_clean <- clean_names(unmatched_species$Species) |>
    str_replace_all(" ", "_")
  n_in_fb  <- sum(orphan_names_clean %in% fb_name_to_speccode$name_clean)
  n_not_fb <- nrow(unmatched_species) - n_in_fb
  cat(sprintf("%d orphans in FishBase (genus absent from tree).\n", n_in_fb))
  cat(sprintf("%d orphans absent from FishBase entirely.\n",        n_not_fb))
}

# =============================================================================
# 10. FILTER TO SPECIES WITH RANGE DATA
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 9 — Filtering to Species With Range Data\n")
cat(strrep("=", 70), "\n")

final_data <- final_trait_data |>
  select(Species, Final_Tree_Tip, Match_Type, everything())

missing_ranges <- setdiff(final_data$IUCN_Name, range_names)

final_data <- filter(final_data, IUCN_Name %in% range_names)

if (length(missing_ranges) > 0) {
  cat(sprintf("Removed %d species lacking range data.\n", length(missing_ranges)))
} else {
  cat("All species had range data; none removed.\n")
}

# Deduplicate again by IUCN_Name after range filter
dupes_before <- nrow(final_data)

final_data <- final_data |>
  group_by(IUCN_Name) |>
  arrange(factor(Match_Type, levels = match_levels), Species) |>
  summarise(
    across(everything(), ~ { nna <- .x[!is.na(.x)]; if (length(nna) == 0) NA else nna[1] }),
    .groups = "drop"
  )

n_dupes <- dupes_before - nrow(final_data)
if (n_dupes > 0)
  cat(sprintf("Merged %d duplicate IUCN_Name entries.\n", n_dupes))

cat(sprintf("Species in final dataset: %d\n", nrow(final_data)))

# =============================================================================
# 11. PRUNE PHYLOGENY TO FINAL SPECIES SET
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 10 — Pruning Phylogeny\n")
cat(strrep("=", 70), "\n")

tips_to_keep <- final_data$Final_Tree_Tip

cat(sprintf("Pruning all %d trees to %d matched species...\n",
            length(all_trees), length(tips_to_keep)))

pruned_trees <- lapply(all_trees, function(t) keep.tip(t, tips_to_keep))
class(pruned_trees) <- "multiPhylo"

cat("Computing maximum clade credibility tree from pruned posterior...\n")
pruned_tree <- maxCladeCred(pruned_trees)

# Sort trait data rows to match tip order of pruned MCC tree
final_data <- final_data |>
  arrange(match(Final_Tree_Tip, pruned_tree$tip.label))

cat(sprintf("Pruned MCC tree has %d tips.\n", length(pruned_tree$tip.label)))

# =============================================================================
# 12. ATTACH THREAT COLUMNS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 11 — Attaching Threat Columns\n")
cat(strrep("=", 70), "\n")

final_data <- final_data |>
  left_join(
    select(iucn_clean, IUCN_Name, Threat_Logging, Threat_Plantation),
    by = "IUCN_Name"
  )

cat(sprintf("Species threatened by logging:    %d\n",
            sum(final_data$Threat_Logging,    na.rm = TRUE)))
cat(sprintf("Species threatened by plantation: %d\n",
            sum(final_data$Threat_Plantation, na.rm = TRUE)))

# =============================================================================
# 13. SAVE FINAL OUTPUTS
# =============================================================================
cat("\n", strrep("=", 70), "\n")
cat("STEP 12 — Saving Outputs\n")
cat(strrep("=", 70), "\n")

# Trait data
write.csv(final_data,
          here("data/processed/traits/fish/fish_trait_data_matched.csv"),
          row.names = FALSE)
cat("Saved: data/processed/traits/fish/fish_trait_data_matched.csv\n")

# MCC pruned tree (Newick)
write.tree(pruned_tree, here("data/processed/traits/fish/final_fish_tree_pruned.tre"))
cat("Saved: data/processed/traits/fish/final_fish_tree_pruned.tre\n")

# All 100 pruned posterior trees
saveRDS(pruned_trees, here("data/processed/traits/fish/pruned_trees_100.rds"))
cat("Saved: data/processed/traits/fish/pruned_trees_100.rds\n")

cat("\n", strrep("=", 70), "\n")
cat("PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Final species:      %d\n", nrow(final_data)))
cat(sprintf("  Trait columns:      %d\n", ncol(final_data) - 3))   # exc. Species/tip/match
cat(sprintf("  Pruned tree tips:   %d\n", length(pruned_tree$tip.label)))
cat(sprintf("  Posterior trees:    %d\n", length(pruned_trees)))
cat("\nOutputs written to:\n")
cat("  outputs/traits/\n")
cat("  data/processed/traits/fish/\n")

