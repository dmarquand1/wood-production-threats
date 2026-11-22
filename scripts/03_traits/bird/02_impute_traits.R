rm(list=ls())
gc()
# ==============================================================================
# PHYLOGENETIC TRAIT IMPUTATION PIPELINE — BIRDS
# ------------------------------------------------------------------------------
# OBJECTIVE: Fill missing trait data for birds using evolutionary
#            relationships encoded in a time-calibrated phylogenetic tree.
#
# METHODOLOGY: Mirrors the reptile imputation pipeline exactly.
#
# WORKFLOW:
#   1.  Load data and tree
#   1b. Fetch IUCN Red List habitat data → create binary IUCN_* habitat vars
#   2.  Data preparation (log-transform, encode binary/dummy traits)
#   3.  Define trait families (continuous / binary / multinomial)
#   4.  Filter: drop traits with <60% coverage
#   4b. Filter: drop binary traits with <1% positive prevalence
#   4c. Flag multinomial traits with any class <5% prevalence
#   4d. Recode: merge rare multinomial classes
#   5.  Validate: 10-fold repeated CV per trait family
#        Continuous  -> Rphylopars (lambda model) | metric: R²    | threshold: >0.70
#        Binary      -> Rphylopars LPM + prev.    | metric: Kappa | threshold: >0.60
#        Multinomial -> Rphylopars OvR + argmax    | metric: Kappa | threshold: >0.60
#   6.  Filter: drop traits that fail CV thresholds
#   7.  Final imputation using phylopars for all trait families
#   8.  Re-scale continuous traits and save output
# ==============================================================================


# ==============================================================================
# 0. LIBRARIES
# ==============================================================================

library(tidyverse)   # Data manipulation and pipe syntax
library(Rphylopars)  # Fast phylogenetic imputation (all trait families)
library(ape)         # Core phylogenetics (read trees, prune tips)
library(here)        # Clean relative file paths
library(caret)       # confusionMatrix and Kappa
library(doSNOW)      # Parallel backend with progress bar support
library(foreach)     # Parallel loop syntax
library(parallel)    # detectCores()
library(rredlist)    # IUCN Red List API client (>= 1.0.0 required)
library(dplyr)
library(stringr)

clean_join_name <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_replace_all("\u00A0", " ") %>%
    str_squish()
}


# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("--- Step 1: Loading Data ---\n")

analysis_df <- read.csv(
  "data/processed/traits/bird/bird_trait_data_matched.csv",
  stringsAsFactors = FALSE
)
tree        <- read.tree("data/processed/traits/bird/final_bird_tree_pruned.tre")

if ("phylo_match" %in% names(analysis_df)) {
  analysis_df <- analysis_df %>% rename(phylo_name = phylo_match)
} else if ("Tree_Name" %in% names(analysis_df)) {
  analysis_df <- analysis_df %>% rename(phylo_name = Tree_Name)
} else {
  stop("bird_trait_data_matched.csv is missing the matched tree-tip column ('phylo_match'/'Tree_Name'). Run 01_prepare_data.R first.")
}

# Fix tree: resolve polytomies and replace zero-length branches
tree <- multi2di(tree)
tree$edge.length[tree$edge.length == 0] <- 1e-6

# Keep only species present in the tree — phylogenetic models require an exact
# match between dataframe rows and tree tip labels.
analysis_df <- analysis_df %>%
  mutate(
    phylo_name = str_replace_all(phylo_name, " ", "_"),
    spatial_name = Trait_Name
  ) %>%
  filter(phylo_name %in% tree$tip.label)

cat(sprintf("Loaded %d species matched to tree.\n", nrow(analysis_df)))


# ==============================================================================
# 1b. IUCN RED LIST HABITAT DATA
# ------------------------------------------------------------------------------
# OBJECTIVE: Fetch habitat data from the IUCN Red List API for every species in
#   the tree, then create broad binary habitat variables (Forest, Savanna, …)
#   that replace the Primary.Habitat column previously used as a multinomial.
# ==============================================================================
cat("\n--- Step 1b: IUCN Red List Habitat Data ---\n")

assessment_df_consolidated <- readRDS(
  here("data", "processed", "assessments", "bird_assessment_habitat_summary.rds")
) %>%
  mutate(scientific_name = clean_join_name(scientific_name))
if (FALSE) {
assessment_df   <- rl_assessment_extract(
  assessment_list,
  c("taxon", "red_list_category__code", "threats", "habitats"),
  format  = "df",
  flatten = TRUE
) %>%
  rename(
    category = code.x,
    habitat  = code
  ) %>%
  mutate(
    scientific_name = clean_join_name(scientific_name)
  )

# Create binary habitat flags per row (one row per species × habitat entry)
df_with_flags <- assessment_df %>%
  mutate(
    is_forest     = str_detect(habitat, "^1(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_savanna    = str_detect(habitat, "^2(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_shrub      = str_detect(habitat, "^3(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_grass      = str_detect(habitat, "^4(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_wetland    = str_detect(habitat, "^5(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_rocky      = str_detect(habitat, "^6(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_cave       = str_detect(habitat, "^7(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_desert     = str_detect(habitat, "^8(\\_|$)")              & (is.na(suitability) | suitability != "Marginal"),
    is_marine     = str_detect(habitat, "^(9|10|11|12|13)(\\_|$)") & (is.na(suitability) | suitability != "Marginal"),
    is_artificial = str_detect(habitat, "^(14|15|16)(\\_|$)")     & (is.na(suitability) | suitability != "Marginal"),
    is_other      = str_detect(habitat, "^17(\\_|$)")             & (is.na(suitability) | suitability != "Marginal"),
    is_unknown    = str_detect(habitat, "^18(\\_|$)")             & (is.na(suitability) | suitability != "Marginal")
  )

# Collapse to one row per species: TRUE if the species uses that habitat in
# at least one non-marginal entry
assessment_df_consolidated <- df_with_flags %>%
  group_by(scientific_name, category) %>%
  summarise(
    Forest     = any(is_forest,     na.rm = TRUE),
    Savanna    = any(is_savanna,    na.rm = TRUE),
    Shrubland  = any(is_shrub,      na.rm = TRUE),
    Grassland  = any(is_grass,      na.rm = TRUE),
    Wetlands   = any(is_wetland,    na.rm = TRUE),
    Rocky      = any(is_rocky,      na.rm = TRUE),
    Caves      = any(is_cave,       na.rm = TRUE),
    Desert     = any(is_desert,     na.rm = TRUE),
    Marine     = any(is_marine,     na.rm = TRUE),
    Artificial = any(is_artificial, na.rm = TRUE),
    Other      = any(is_other,      na.rm = TRUE),
    Unknown    = any(is_unknown,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.logical), as.integer))%>% 
  mutate(habitat_breadth = rowSums(select(., Forest:Unknown)))
}

iucn_habitat_cols <- c(
  "Forest", "Savanna", "Shrubland", "Grassland", "Wetlands",
  "Rocky", "Caves", "Desert", "Marine", "Artificial", "Other", "Unknown", "habitat_breadth"
)

cat(sprintf("IUCN habitat data consolidated for %d species.\n",
            nrow(assessment_df_consolidated)))


# Join IUCN binary habitat columns.
# NOTE: adjust the left-hand join key ("IUCN_Name") if your bird dataframe
#       stores the IUCN-matched species name under a different column.
analysis_df <- analysis_df %>%
  mutate(join_name = clean_join_name(Trait_Name)) %>%
  left_join(
    assessment_df_consolidated %>% select(scientific_name, all_of(iucn_habitat_cols)),
    by = c("join_name" = "scientific_name")
  ) %>%
  select(-join_name)

# ==============================================================================
# 2. DATA PREPARATION
# ==============================================================================
cat("\n--- Step 2: Data Preparation ---\n")

analysis_df <- analysis_df %>%
  mutate(
    # --- Continuous traits (log-transform + scale) ---
    Body_Mass_log          = scale(log(Average.Mass)),
    Elevational_Range_log  = scale(log(Elevational.Range + 1)),
    ESI_scaled             = scale(ESI),
    Clutch_Size_log        = scale(log((Clutch_Min + Clutch_Max) / 2)),
    HWI_log = scale(log(HWI)),
    Habitat_Breadth_log = scale(log(habitat_breadth)),
    
    # --- Simple binary traits ---
    # Island status (already 0/1)
    Island_Status = ISL,
    
    # Flightlessness (1 = flightless or partial, 0 = volant)
    Is_Flightless = ifelse(Flightlessness %in% c("Yes", "Partial"), 1, 0),
    
    # Migration: sedentary species coded 0, all others 1
    Is_Migratory = case_when(
      is.na(Mig) & is.na(Alt) & is.na(Irreg) & is.na(Disp) & is.na(Sed) ~ NA_real_,
      Sed == 1 ~ 0,
      TRUE     ~ 1
    ),
    
    # Nesting substrates
    Is_Tree_Nester = case_when(
      is.na(Nest_SBS) ~ NA_real_,
      grepl("T|S", Nest_SBS) ~ 1,
      TRUE ~ 0
    ),
    Is_Stump_Nester = case_when(
      is.na(Nest_SBS) ~ NA_real_,
      grepl("C", Nest_SBS) ~ 1,
      TRUE ~ 0
    ),
    Is_Ground_Nester = case_when(
      is.na(Nest_SBS) ~ NA_real_,
      grepl("G", Nest_SBS) ~ 1,
      TRUE ~ 0
    )
  )

# Binary diet variables derived from weighted diet scores.
# Logic: 1 if the species has any non-zero weight for that diet category,
#        0 if the weight is explicitly 0,
#        NA if ALL diet columns are NA (no dietary information recorded).
diet_wt_cols    <- c("IN.Wt", "FR.Wt", "NE.Wt", "SE.Wt", "VE.Wt",
                     "FI.Wt", "SC.Wt", "PL.Wt", "MS.Wt")
diet_binary_cols <- c("Diet_Invertebrate", "Diet_Fruit", "Diet_Nectar",
                      "Diet_Seed", "Diet_Vertebrate", "Diet_Fish",
                      "Diet_Scavenger", "Diet_PlantMaterial", "Diet_Misc")

all_diet_na <- rowSums(!is.na(analysis_df[, diet_wt_cols])) == 0

for (i in seq_along(diet_wt_cols)) {
  wt_col  <- diet_wt_cols[i]
  bin_col <- diet_binary_cols[i]
  analysis_df[[bin_col]] <- case_when(
    all_diet_na                          ~ NA_real_,   # No diet info at all
    is.na(analysis_df[[wt_col]])         ~ 0,          # Other diets recorded but not this one
    analysis_df[[wt_col]] > 0            ~ 1,          # Diet present
    TRUE                                 ~ 0           # Weight explicitly 0
  )
}

# Handle split species (duplicates): average numeric columns
n_dups <- analysis_df %>% group_by(phylo_name) %>% filter(n() > 1) %>% nrow()
if (n_dups > 0) {
  cat(sprintf("   [INFO] Aggregating %d rows of split species...\n", n_dups))
  analysis_df <- analysis_df %>%
    group_by(phylo_name, spatial_name) %>%
    summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = "drop")
}


# ==============================================================================
# 3. DEFINE TRAIT FAMILIES
# ==============================================================================
cat("\n--- Step 3: Defining Trait Families ---\n")

continuous_traits <- c(
  "Body_Mass_log", "Elevational_Range_log", "ESI_scaled", "Clutch_Size_log", 
  "HWI_log", "Habitat_Breadth_log", ""
)

binary_traits <- c(
  "Island_Status", "Is_Flightless", "Is_Migratory",
  "Is_Tree_Nester", "Is_Stump_Nester", "Is_Ground_Nester",
  # IUCN Red List habitat binary variables (replace Primary.Habitat)
  iucn_habitat_cols,
  # Diet binary variables (derived from weighted diet scores; replace Primary.Diet)
  diet_binary_cols
)

multinomial_traits <- c()

# Keep only traits that actually exist in the dataframe after transformations
candidate_traits <- intersect(
  c(continuous_traits, binary_traits),
  names(analysis_df)
)


cat(sprintf("Candidate traits: %d continuous | %d binary | %d multinomial\n",
            sum(candidate_traits %in% continuous_traits),
            sum(candidate_traits %in% binary_traits),
            sum(candidate_traits %in% multinomial_traits)
            ))


# ==============================================================================
# 4. FILTER: COVERAGE (>60%)
# ==============================================================================
cat("\n--- Step 4: Coverage Check (threshold: >60%) ---\n")

coverage_results <- data.frame(
  Trait    = candidate_traits,
  Coverage = sapply(candidate_traits, function(t) mean(!is.na(analysis_df[[t]]))),
  stringsAsFactors = FALSE
) %>%
  mutate(Status = ifelse(Coverage > 0.60, "PASS", "FAIL")) %>%
  arrange(Coverage)

cat(sprintf(
  "Traits passing 60%% coverage: %d / %d\n",
  sum(coverage_results$Status == "PASS"),
  nrow(coverage_results)
))
print(coverage_results, row.names = FALSE)

traits_passing_coverage <- coverage_results %>%
  filter(Status == "PASS") %>%
  pull(Trait)

traits_100_coverage <- coverage_results %>%
  filter(Coverage == 1) %>%
  pull(Trait)

if (length(traits_100_coverage) > 0) {
  cat(sprintf("\n[INFO] %d trait(s) have 100%% coverage and will bypass CV:\n",
              length(traits_100_coverage)))
  cat(paste(" -", traits_100_coverage, collapse = "\n"), "\n")
}

# Subset each family to coverage-passing traits only
eval_continuous  <- intersect(continuous_traits,  traits_passing_coverage)
eval_binary      <- intersect(binary_traits,      traits_passing_coverage)


# ==============================================================================
# 4b. FILTER: MINIMUM POSITIVE PREVALENCE FOR BINARY TRAITS (threshold: >= 1%)
# ==============================================================================
cat("\n--- Step 4b: Minimum Positive Prevalence Check for Binary Traits (threshold: >= 1%) ---\n")

prev_threshold_bin <- 0.01

prevalence_results <- data.frame(
  Trait      = eval_binary,
  N_Observed = sapply(eval_binary, function(t) {
    sum(!is.na(as.numeric(analysis_df[[t]])))
  }),
  N_Positive = sapply(eval_binary, function(t) {
    sum(as.numeric(analysis_df[[t]]) == 1, na.rm = TRUE)
  }),
  N_Negative = sapply(eval_binary, function(t) {
    sum(as.numeric(analysis_df[[t]]) == 0, na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Prevalence = N_Positive / N_Observed,
    Status     = ifelse(
      Prevalence >= prev_threshold_bin,
      "PASS",
      sprintf("FAIL (%.1f%% < 1%% threshold)", Prevalence * 100)
    )
  ) %>%
  arrange(Prevalence)

n_pass <- sum(prevalence_results$Status == "PASS")
n_fail <- sum(prevalence_results$Status != "PASS")

cat(sprintf(
  "Binary traits passing minimum prevalence filter: %d / %d\n",
  n_pass, nrow(prevalence_results)
))
cat(sprintf("  Dropped (positive prevalence < 1%%): %d\n\n", n_fail))
print(prevalence_results, row.names = FALSE)

failed_bin <- prevalence_results %>% filter(Status != "PASS")
if (nrow(failed_bin) > 0) {
  cat("\nDropped binary traits:\n")
  for (i in seq_len(nrow(failed_bin))) {
    cat(sprintf(
      "  [DROP] %-40s Prevalence = %.1f%%  (N positive = %d / N observed = %d)\n",
      failed_bin$Trait[i],
      failed_bin$Prevalence[i] * 100,
      failed_bin$N_Positive[i],
      failed_bin$N_Observed[i]
    ))
  }
}

eval_binary <- prevalence_results %>%
  filter(Status == "PASS") %>%
  pull(Trait)

cat(sprintf("\n%d binary traits carried forward to cross-validation.\n", length(eval_binary)))

# ==============================================================================
# 5. CROSS-VALIDATION
# ==============================================================================
cat("\n\n=======================================================\n")
cat("   SECTION 5: 10-FOLD REPEATED CROSS-VALIDATION\n")
cat("=======================================================\n")

# Working dataframe for validation — species column must be named 'species'
df_valid <- analysis_df %>%
  dplyr::select(species = phylo_name, all_of(traits_passing_coverage))

# Parallel cluster — leave one core free for the OS
n_cores <- detectCores() - 1
cl      <- makeCluster(n_cores)
registerDoSNOW(cl)

# Stores primary validation metric (R² or Kappa) for each trait
validation_scores <- list()


# -------------------------------------------------------
# A. CONTINUOUS TRAITS
#    Model:   Rphylopars (Pagel's lambda)
#    Metric:  R² (coefficient of determination)
#    Repeats: 5 x 10-fold = 50 CV estimates per trait
# -------------------------------------------------------
cat("\n--- A. Continuous Traits (Rphylopars, R²) ---\n")

n_folds_cont   <- 10
n_repeats_cont <- 5

for (trait in eval_continuous) {
  
  if (trait %in% traits_100_coverage) {
    cat(sprintf("\n[%d/%d] %s | Skipping CV (100%% coverage already)\n",
                which(eval_continuous == trait), length(eval_continuous), trait))
    validation_scores[[trait]] <- 1.0
    next
  }
  
  valid_species <- df_valid %>% filter(!is.na(!!sym(trait))) %>% pull(species)
  if (length(valid_species) < n_folds_cont * 2) next
  
  cat(sprintf("\n[%d/%d] %s | %d species\n",
              which(eval_continuous == trait), length(eval_continuous),
              trait, length(valid_species)))
  
  repeat_r2 <- c()
  
  for (r in 1:n_repeats_cont) {
    
    set.seed(100 + r)
    shuffled_sp <- sample(valid_species)
    folds       <- cut(seq_along(shuffled_sp), breaks = n_folds_cont, labels = FALSE)
    
    pb   <- txtProgressBar(max = n_folds_cont, style = 3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    
    fold_results <- foreach(
      f = 1:n_folds_cont, .combine = rbind,
      .packages     = "Rphylopars",
      .export       = c("tree", "df_valid", "trait", "shuffled_sp", "folds"),
      .options.snow = opts
    ) %dopar% {
      
      test_species <- shuffled_sp[which(folds == f)]
      
      train_data <- df_valid[, c("species", trait)]
      train_data[train_data$species %in% test_species, trait] <- NA
      
      tryCatch({
        p      <- phylopars(trait_data = train_data, tree = tree,
                            model = "lambda", pheno_error = FALSE)
        preds  <- p$anc_recon[test_species, trait]
        truths <- df_valid[[trait]][match(test_species, df_valid$species)]
        cbind(preds, truths)
      }, error = function(e) {
        message("Fold ", f, " error (", trait, "): ", e$message); NULL
      })
    }
    close(pb)
    
    if (!is.null(fold_results)) {
      r2        <- cor(fold_results[, 1], fold_results[, 2], use = "complete.obs")^2
      repeat_r2 <- c(repeat_r2, r2)
    }
  }
  
  if (length(repeat_r2) > 0) {
    mean_r2 <- mean(repeat_r2)
    sd_r2   <- sd(repeat_r2)
    validation_scores[[trait]] <- mean_r2
    cat(sprintf("   -> Mean R²: %.3f (±%.3f SD) over %d repeats\n",
                mean_r2, sd_r2, n_repeats_cont))
  }
}


# -------------------------------------------------------
# B. BINARY TRAITS
#    Model:   Rphylopars LPM
#    Metric:  Cohen's Kappa
#    Thresholds:
#      (1) Naive fixed threshold of 0.5
#      (2) Prevalence-based threshold: proportion of
#          positives among training species per fold.
#          Adapts automatically to class imbalance with
#          no leakage risk.
#    Final score uses prevalence-based threshold.
#    Repeats: 5 x 10-fold
# -------------------------------------------------------
cat("\n--- B. Binary Traits (Rphylopars LPM, Kappa, prevalence-based threshold) ---\n")

n_folds_bin   <- 10
n_repeats_bin <- 5

oob_predictions    <- list()
optimal_thresholds <- list()

for (trait in eval_binary) {
  
  if (trait %in% traits_100_coverage) {
    cat(sprintf("\n[%d/%d] %s | Skipping CV (100%% coverage already)\n",
                which(eval_binary == trait), length(eval_binary), trait))
    validation_scores[[trait]]  <- 1.0
    optimal_thresholds[[trait]] <- 0.5
    next
  }
  
  valid_species <- df_valid %>% filter(!is.na(!!sym(trait))) %>% pull(species)
  
  if (length(valid_species) < n_folds_bin * 2) {
    cat(sprintf("\nSkipping %s: too few observations (%d)\n", trait, length(valid_species)))
    next
  }
  
  cat(sprintf("\n[%d/%d] %s | %d species\n",
              which(eval_binary == trait), length(eval_binary),
              trait, length(valid_species)))
  
  repeat_kappa_naive <- c()
  repeat_kappa_prev  <- c()
  repeat_thresholds  <- c()
  all_oob            <- NULL
  
  for (r in 1:n_repeats_bin) {
    
    set.seed(123 + r)
    shuffled_sp <- sample(valid_species)
    folds       <- cut(seq_along(shuffled_sp), breaks = n_folds_bin, labels = FALSE)
    
    pb   <- txtProgressBar(max = n_folds_bin, style = 3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    
    fold_results <- foreach(
      f = 1:n_folds_bin, .combine = rbind,
      .packages     = c("Rphylopars"),
      .export       = c("tree", "df_valid", "trait", "shuffled_sp", "folds"),
      .options.snow = opts
    ) %dopar% {
      
      test_species  <- shuffled_sp[which(folds == f)]
      train_species <- shuffled_sp[which(folds != f)]
      
      train_data          <- df_valid[, c("species", trait)]
      train_data[[trait]] <- as.numeric(as.character(train_data[[trait]]))
      train_data[train_data$species %in% test_species, trait] <- NA
      
      tryCatch({
        p <- phylopars(trait_data = train_data, tree = tree,
                       model = "lambda", pheno_error = FALSE)
        
        test_preds  <- pmax(0, pmin(1, p$anc_recon[test_species, trait]))
        test_truths <- as.numeric(as.character(
          df_valid[[trait]][match(test_species, df_valid$species)]))
        
        # Prevalence-based threshold from training fold
        train_truths     <- as.numeric(as.character(
          df_valid[[trait]][match(train_species, df_valid$species)]))
        train_truths_obs <- train_truths[!is.na(train_truths)]
        
        prev_thresh <- if (length(train_truths_obs) > 0 && sum(train_truths_obs) > 0)
          mean(train_truths_obs)
        else
          0.5
        
        data.frame(
          Pred_raw    = test_preds,
          Truth       = test_truths,
          Fold_thresh = prev_thresh
        )
        
      }, error = function(e) {
        message("Fold ", f, " error (", trait, "): ", e$message); NULL
      })
    }
    close(pb)
    
    if (!is.null(fold_results) && nrow(fold_results) > 0) {
      
      fold_results <- fold_results[
        !is.na(fold_results$Pred_raw) & !is.na(fold_results$Truth), ]
      
      if (nrow(fold_results) == 0) next
      
      pred_naive <- as.integer(fold_results$Pred_raw > 0.5)
      pred_prev  <- as.integer(fold_results$Pred_raw > fold_results$Fold_thresh)
      
      mean_thresh_r     <- mean(fold_results$Fold_thresh)
      repeat_thresholds <- c(repeat_thresholds, mean_thresh_r)
      
      if (length(unique(pred_naive)) == 1)
        message("WARNING: naive threshold gives constant predictions for ", trait)
      if (length(unique(pred_prev)) == 1)
        message("WARNING: prevalence threshold gives constant predictions for ", trait)
      
      tryCatch({
        k_naive <- confusionMatrix(
          factor(pred_naive,         levels = c(0, 1)),
          factor(fold_results$Truth, levels = c(0, 1))
        )$overall["Kappa"]
        if (!is.na(k_naive)) repeat_kappa_naive <- c(repeat_kappa_naive, k_naive)
      }, error = function(e) message("Kappa (naive) error for ", trait, ": ", e$message))
      
      tryCatch({
        k_prev <- confusionMatrix(
          factor(pred_prev,          levels = c(0, 1)),
          factor(fold_results$Truth, levels = c(0, 1))
        )$overall["Kappa"]
        if (!is.na(k_prev)) repeat_kappa_prev <- c(repeat_kappa_prev, k_prev)
      }, error = function(e) message("Kappa (prev) error for ", trait, ": ", e$message))
      
      if (r == 1) all_oob <- fold_results[, c("Pred_raw", "Truth", "Fold_thresh")]
    }
  }
  
  if (length(repeat_kappa_prev) > 0) {
    
    mean_kappa_naive <- if (length(repeat_kappa_naive) > 0) mean(repeat_kappa_naive) else NA_real_
    mean_kappa_prev  <- mean(repeat_kappa_prev)
    sd_kappa_prev    <- sd(repeat_kappa_prev)
    mean_thresh      <- mean(repeat_thresholds)
    sd_thresh        <- sd(repeat_thresholds)
    delta            <- mean_kappa_prev - mean_kappa_naive
    
    validation_scores[[trait]]  <- mean_kappa_prev
    optimal_thresholds[[trait]] <- mean_thresh
    oob_predictions[[trait]]    <- all_oob
    
    cat(sprintf(
      "   -> Naive (0.5)                   Kappa: %.3f\n",
      mean_kappa_naive))
    cat(sprintf(
      "   -> Prevalence (%.2f ± %.2f) Kappa: %.3f (± %.3f SD)  [delta: %+.3f]\n",
      mean_thresh, sd_thresh, mean_kappa_prev, sd_kappa_prev, delta))
    
  } else {
    cat(sprintf("   -> %s: no valid folds returned.\n", trait))
  }
}

# ==============================================================================
# 6. FILTER TRAITS BY VALIDATION THRESHOLDS
# ==============================================================================
cat("\n\n=======================================================\n")
cat("   SECTION 6: FILTERING BY VALIDATION SCORE\n")
cat("=======================================================\n")
cat("Thresholds: Continuous R² > 0.70 | Categorical Kappa > 0.60\n\n")

passed_traits <- c()

for (trait in traits_passing_coverage) {
  score <- validation_scores[[trait]]
  
  if (trait %in% eval_continuous) {
    threshold <- 0.70
    passes    <- !is.null(score) && !is.na(score) && score > threshold
    label     <- sprintf("R²: %.3f", ifelse(is.null(score), NA, score))
    
  } else if (trait %in% c(eval_binary)) {
    threshold <- 0.60
    passes    <- !is.null(score) && !is.na(score) && score > threshold
    label     <- sprintf("Kappa: %.3f", ifelse(is.null(score), NA, score))
    
  } else {
    next
  }
  
  if (passes) {
    passed_traits <- c(passed_traits, trait)
    cat(sprintf("[PASS] %-40s %s\n", trait, label))
  } else {
    cat(sprintf("[DROP] %-40s %s\n", trait, label))
  }
}

cat(sprintf("\n%d / %d traits passed validation.\n",
            length(passed_traits), length(traits_passing_coverage)))


# ==============================================================================
# 7. FINAL IMPUTATION
# ==============================================================================
cat("\n\n=======================================================\n")
cat("   SECTION 7: FINAL IMPUTATION\n")
cat("=======================================================\n")
cat("Using phylopars for all trait families for methodological consistency.\n\n")

imputed_df <- data.frame(species = tree$tip.label)

final_cont <- intersect(passed_traits, eval_continuous)
final_bin  <- intersect(passed_traits, eval_binary)

# Split each family into traits that need gap-filling vs traits that are already
# complete. 100%-coverage traits bypass phylopars entirely — their values are
# copied straight from df_valid, avoiding a model fit whose output would be
# discarded anyway.
final_cont_impute <- setdiff(final_cont, traits_100_coverage)
final_cont_copy   <- intersect(final_cont, traits_100_coverage)

final_bin_impute  <- setdiff(final_bin,  traits_100_coverage)
final_bin_copy    <- intersect(final_bin,  traits_100_coverage)


cat(sprintf(
  "Traits to impute  : %d continuous | %d binary\n",
  length(final_cont_impute), length(final_bin_impute)
))
cat(sprintf(
  "Traits to copy    : %d continuous | %d binary\n\n",
  length(final_cont_copy), length(final_bin_copy)
))


# --- Pass-through: copy 100%-coverage traits directly ---
all_copy_traits <- c(final_cont_copy, final_bin_copy)
if (length(all_copy_traits) > 0) {
  cat(sprintf("Copying %d fully-observed trait(s) directly...\n", length(all_copy_traits)))
  for (trait in all_copy_traits) {
    imputed_df[[trait]] <- df_valid[[trait]][match(tree$tip.label, df_valid$species)]
    cat(sprintf("  -> %s  [copied]\n", trait))
  }
  cat("   -> Done.\n\n")
}


# --- A. Continuous Imputation ---
# Only traits with gaps are passed to phylopars. All are imputed jointly in a
# single call so the model can exploit cross-trait correlations via the
# multivariate Brownian motion covariance structure.
if (length(final_cont_impute) > 0) {
  cat(sprintf("Imputing %d continuous traits (joint phylopars)...\n",
              length(final_cont_impute)))
  
  train_cont <- df_valid[, c("species", final_cont_impute)]
  p_cont     <- phylopars(trait_data = train_cont, tree = tree,
                          model = "lambda", pheno_error = FALSE)
  
  cont_recon         <- as.data.frame(
    p_cont$anc_recon[tree$tip.label, final_cont_impute, drop = FALSE])
  cont_recon$species <- rownames(cont_recon)
  imputed_df         <- left_join(imputed_df, cont_recon, by = "species")
  cat("   -> Done.\n\n")
} else {
  cat("No continuous traits require imputation (all at 100% coverage).\n\n")
}


# --- B. Binary Imputation ---
# Each trait with gaps is imputed separately via phylopars LPM.
# Threshold fixed at 0.5; observed values always preserved.
if (length(final_bin_impute) > 0) {
  cat(sprintf("Imputing %d binary traits (phylopars LPM)...\n", length(final_bin_impute)))
  
  for (trait in final_bin_impute) {
    cat(sprintf("  -> %s\n", trait))
    
    train_bin          <- df_valid[, c("species", trait)]
    train_bin[[trait]] <- as.numeric(as.character(train_bin[[trait]]))
    
    n_obs <- sum(!is.na(train_bin[[trait]]))
    if (n_obs == 0) {
      cat(sprintf("     [SKIP] %s has 0 observed values — dropping from imputation.\n", trait))
      next
    }
    
    p_bin <- phylopars(trait_data = train_bin, tree = tree,
                       model = "lambda", pheno_error = FALSE)
    preds <- p_bin$anc_recon[tree$tip.label, trait]
    preds <- pmax(0, pmin(1, preds))  # Clamp to [0,1]
    
    original <- as.numeric(as.character(
      df_valid[[trait]][match(tree$tip.label, df_valid$species)]))
    
    imputed_df[[trait]] <- ifelse(is.na(original), as.integer(preds > 0.5), original)
  }
  cat("   -> Done.\n\n")
} else {
  cat("No binary traits require imputation (all at 100% coverage).\n\n")
}


# ==============================================================================
# 8. POST-IMPUTATION RE-SCALING & SAVE
# ==============================================================================
cat("\n--- Step 8: Re-scaling Imputed Continuous Traits ---\n")

final_analysis_df <- imputed_df %>%
  mutate(
    across(
      all_of(final_cont_impute),
      ~ as.numeric(scale(.))
    )
  )

cat(sprintf("Re-scaled %d continuous traits to mean=0, sd=1.\n", length(final_cont_impute)))

write.csv(final_analysis_df, "data/processed/traits/bird/bird_traits_imputed_ready.csv", row.names = FALSE)
write.csv(final_analysis_df, "data/processed/traits/bird/Bird_Traits_Imputed_Final.csv", row.names = FALSE)
cat("\nSuccess! Final data saved to 'data/processed/traits/bird/bird_traits_imputed_ready.csv'\n")

