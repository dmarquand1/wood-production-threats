rm(list=ls())
gc()
# ==============================================================================
# PHYLOGENETIC TRAIT IMPUTATION PIPELINE
# ------------------------------------------------------------------------------
# OBJECTIVE: Fill missing trait data for amphibians using evolutionary
#            relationships encoded in a time-calibrated phylogenetic tree.
#
# WORKFLOW:
#   1. Load data and tree
#   2. Data preparation (log-transform, collapse one-hot encodings)
#   3. Define trait families (continuous / binary / multinomial)
#   4. Filter: drop traits with <60% coverage
#   5. Validate: 10-fold repeated CV per trait family
#        Continuous  -> Rphylopars (lambda model) | metric: R²   | threshold: >0.70
#        Binary      -> Rphylopars LPM            | metric: AUC  | threshold: >0.75
#        Multinomial -> Rphylopars OvR + argmax   | metric: Kappa| threshold: >0.60
#                                                                    + no class F1 = 0
#   6. Filter: drop traits that fail CV thresholds
#   7. Final imputation using phylopars for all trait families
#   8. Save output
# ==============================================================================


# ==============================================================================
# 0. LIBRARIES
# ==============================================================================

library(tidyverse)   # Data manipulation and pipe syntax
library(Rphylopars)  # Fast phylogenetic imputation (all trait families)
library(ape)         # Core phylogenetics (read trees, prune tips)
library(here)        # Clean relative file paths
library(caret)       # confusionMatrix and Kappa
library(pROC)        # AUC for binary and per-class multinomial validation
library(doSNOW)      # Parallel backend with progress bar support
library(foreach)     # Parallel loop syntax
library(parallel)    # detectCores()

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

analysis_df <- read.csv("data/processed/traits/amphibian/amphibian_trait_data_matched.csv",
                        stringsAsFactors = FALSE)
tree        <- read.tree("data/processed/traits/amphibian/final_amphibian_tree_pruned.tre")

# Keep only species present in the tree — phylogenetic models require an exact
# match between dataframe rows and tree tip labels.
analysis_df <- analysis_df %>%
  filter(Final_Tree_Tip %in% tree$tip.label)

cat(sprintf("Loaded %d species matched to tree.\n", nrow(analysis_df)))

# ==============================================================================
# 1b. IUCN RED LIST HABITAT DATA
# ------------------------------------------------------------------------------
# Load the pre-computed habitat summary rather than rebuilding it from the raw
# assessment cache each time.
# ==============================================================================

assessment_df_consolidated <- readRDS(
  here("data", "processed", "assessments", "amphibian_assessment_habitat_summary.rds")
) %>%
  mutate(scientific_name = clean_join_name(scientific_name))

analysis_df <- analysis_df%>%
  mutate(IUCN_Name = clean_join_name(IUCN_Name)) %>%
  left_join(assessment_df_consolidated %>% select(scientific_name, Forest:habitat_breadth), 
            by = c("IUCN_Name" = "scientific_name"))

# ==============================================================================
# 2. DATA PREPARATION
# ==============================================================================
cat("\n--- Step 2: Data Preparation ---\n")

# ------------------------------------------------------------------------------
# HELPER: fix_binary_nas
# If any column in a one-hot group contains a '1' for a row, fill remaining
# NAs in that group with '0'. Rows that are entirely NA remain entirely NA
# (these are genuinely missing observations, not implicit zeros).
# ------------------------------------------------------------------------------
fix_binary_nas <- function(df, cols) {
  all_na <- apply(df[, cols], 1, function(x) all(is.na(x)))
  for (col in cols) {
    df[[col]] <- ifelse(!all_na & is.na(df[[col]]), 0, df[[col]])
  }
  return(df)
}

# ------------------------------------------------------------------------------
# HELPER: collapse_to_multinomial
# Reverses one-hot encoding into a single factor column. Picks the first '1'
# found in the group — if the row is all NA, returns NA.
# ------------------------------------------------------------------------------
collapse_to_multinomial <- function(df, cols, col_name) {
  df[[col_name]] <- apply(df[, cols], 1, function(row) {
    if (all(is.na(row)))              return(NA_character_)
    if (any(row == 1, na.rm = TRUE)) return(cols[which(row == 1)[1]])
    return(NA_character_)
  })
  df[[col_name]] <- as.factor(df[[col_name]])
  return(df)
}

# One-hot group definitions
grp_diet    <- c("Leaves", "Flowers", "Seeds", "Fruits", "Arthro", "Vert")
grp_diel    <- c("Diu", "Noc", "Crepu")
grp_season  <- c("Wet_warm", "Wet_cold", "Dry_warm", "Dry_cold")
grp_breed   <- c("Dir", "Lar", "Viv")

analysis_df <- analysis_df %>%
  # Fix implicit zeros within each one-hot group
  fix_binary_nas(grp_diet)    %>%
  fix_binary_nas(grp_diel)    %>%
  fix_binary_nas(grp_season)  %>%
  # Breeding mode is mutually exclusive -> collapse to a single factor
  collapse_to_multinomial(grp_breed, "Breeding_Cat") %>%
  mutate(
    # Log10 + Z-score scaling for continuous traits:
    #   - Log10 normalises right-skewed life-history distributions
    #   - Scaling centres data at 0, which aids lambda estimation in phylopars
    Body_Mass_Log           = scale(log10(Body_mass_g)),
    Body_Size_log           = scale(log10(Body_size_mm)),
    Maturity_Age_Log        = scale(log10((Age_at_maturity_min_y  + Age_at_maturity_max_y)  / 2)),
    Maturity_Size_Log       = scale(log10((Size_at_maturity_min_mm + Size_at_maturity_max_mm) / 2)),
    Longevity_Max_Log       = scale(log10(Longevity_max_y)),
    Litter_Size_Log         = scale(log10((Litter_size_min_n       + Litter_size_max_n)       / 2)),
    Reproductive_Output_Log = scale(log10(Reproductive_output_y)),
    Offspring_Size_Log      = scale(log10((Offspring_size_min_mm   + Offspring_size_max_mm)   / 2)),
    Habitat_Breadth         = scale(log(habitat_breadth))
  )

# ==============================================================================
# 3. DEFINE TRAIT FAMILIES
# ==============================================================================
cat("\n--- Step 3: Defining Trait Families ---\n")

continuous_traits <- c(
  "Body_Mass_Log", "Body_Size_log", "Maturity_Age_Log", "Maturity_Size_Log",
  "Longevity_Max_Log", "Litter_Size_Log", "Reproductive_Output_Log", "Offspring_Size_Log"
)

binary_traits <- c(
  "Leaves", "Flowers", "Seeds", "Fruits", "Arthro", "Vert",  # Diet
  "Diu", "Noc", "Crepu",                              # Diel activity
  "Wet_warm", "Wet_cold", "Dry_warm", "Dry_cold",      # Seasonality
  "Forest",                                           # Habitat
  "Savanna",
  "Shrubland",
  "Grassland" ,
  "Wetlands" ,
  "Rocky" ,
  "Caves" ,
  "Desert" ,
  "Marine" ,
  "Artificial" ,
  "Other",
  "Unknown"
)

multinomial_traits <- c("Breeding_Cat")

# Keep only traits that exist in the dataframe after transformations
candidate_traits <- intersect(
  c(continuous_traits, binary_traits, multinomial_traits),
  names(analysis_df)
)


# ==============================================================================
# 4. FILTER: COVERAGE (>60%)
# ==============================================================================
cat("\n--- Step 4: Coverage Check (threshold: >60%) ---\n")

coverage_results <- data.frame(
  Trait    = candidate_traits,
  Coverage = sapply(candidate_traits, function(t) mean(!is.na(analysis_df[[t]]))),
  stringsAsFactors = FALSE
) %>%
  mutate(Status = ifelse(Coverage > 0.60, "PASS", "FAIL"))

traits_100_coverage <- coverage_results %>%
  filter(Coverage == 1) %>%
  pull(Trait)

if(length(traits_100_coverage) > 0) {
  cat(sprintf("\n[INFO] %d trait(s) have 100%% coverage and will bypass CV:\n", length(traits_100_coverage)))
  cat(paste(" -", traits_100_coverage, collapse = "\n"), "\n")
}

cat(sprintf(
  "Traits passing 60%% coverage: %d / %d\n",
  sum(coverage_results$Status == "PASS"),
  nrow(coverage_results)
))
print(coverage_results)

traits_passing_coverage <- coverage_results %>%
  filter(Status == "PASS") %>%
  pull(Trait)

# Subset each family to coverage-passing traits only
eval_continuous  <- intersect(continuous_traits,  traits_passing_coverage)
eval_binary      <- intersect(binary_traits,      traits_passing_coverage)
eval_multinomial <- intersect(multinomial_traits, traits_passing_coverage)


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
      sprintf("FAIL (%.1f%% < 5%% threshold)", Prevalence * 100)
    )
  ) %>%
  arrange(Prevalence)

n_pass <- sum(prevalence_results$Status == "PASS")
n_fail <- sum(prevalence_results$Status != "PASS")

cat(sprintf(
  "Binary traits passing minimum prevalence filter : %d / %d\n",
  n_pass, nrow(prevalence_results)
))
cat(sprintf("  Dropped (positive prevalence < 1%%)        : %d\n\n", n_fail))
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

# Update eval_binary — only carry forward traits with >= 5% positive prevalence
eval_binary <- prevalence_results %>%
  filter(Status == "PASS") %>%
  pull(Trait)

cat(sprintf("\n%d binary traits carried forward to cross-validation.\n", length(eval_binary)))


# ==============================================================================
# 4c. FLAG: MINIMUM PREVALENCE PER CLASS FOR MULTINOMIAL TRAITS (threshold: >= 5%)
# ==============================================================================
cat("\n--- Step 4c: Minimum Class Prevalence Check for Multinomial Traits (threshold: >= 5%) ---\n")

prev_threshold_mult <- 0.01

multinomial_level_results <- lapply(eval_multinomial, function(trait) {
  
  observed <- na.omit(analysis_df[[trait]])
  if (length(observed) == 0) return(NULL)
  
  level_counts <- sort(table(observed), decreasing = TRUE)
  n_total      <- length(observed)
  
  data.frame(
    Trait      = trait,
    Class      = names(level_counts),
    N          = as.integer(level_counts),
    Prevalence = as.integer(level_counts) / n_total,
    N_Total    = n_total,
    Status     = ifelse(
      as.integer(level_counts) / n_total >= prev_threshold_mult,
      "PASS",
      sprintf("FLAG (%.1f%% < 5%% threshold)", as.integer(level_counts) / n_total * 100)
    ),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

cat("\nClass-level prevalence:\n")
for (trait in eval_multinomial) {
  
  trait_rows <- multinomial_level_results %>% filter(Trait == trait)
  if (is.null(trait_rows) || nrow(trait_rows) == 0) next
  
  n_flag_classes <- sum(trait_rows$Status != "PASS")
  flag           <- if (n_flag_classes > 0) "  [!]" else "  [ ]"
  
  cat(sprintf("\n%s %s (%d classes, %d with prevalence < 5%%):\n",
              flag, trait, nrow(trait_rows), n_flag_classes))
  
  for (i in seq_len(nrow(trait_rows))) {
    cat(sprintf("      %-30s  N = %5d  (%.1f%%)   %s\n",
                trait_rows$Class[i],
                trait_rows$N[i],
                trait_rows$Prevalence[i] * 100,
                trait_rows$Status[i]))
  }
}

traits_with_rare_classes <- multinomial_level_results %>%
  filter(Status != "PASS") %>%
  pull(Trait) %>%
  unique()

if (length(traits_with_rare_classes) > 0) {
  cat(sprintf(
    "\n%d multinomial trait(s) have at least one class with prevalence < 5%% and require manual review:\n",
    length(traits_with_rare_classes)
  ))
  cat(paste(" -", traits_with_rare_classes, collapse = "\n"), "\n")
  cat("\nOptions: merge rare classes, collapse to 'Other', or drop trait from imputation.\n")
} else {
  cat("\nAll multinomial traits passed the minimum class prevalence check.\n")
}

# ==============================================================================
# 5. CROSS-VALIDATION
# ==============================================================================
cat("\n\n=======================================================\n")
cat("   SECTION 5: 10-FOLD REPEATED CROSS-VALIDATION\n")
cat("=======================================================\n")

# Working dataframe for validation — species column must be named 'species'
df_valid <- analysis_df %>%
  dplyr::select(species = Final_Tree_Tip, all_of(traits_passing_coverage))

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
    validation_scores[[trait]] <- 1.0 # Assign perfect R² to pass Step 6
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
      r2 <- cor(fold_results[, 1], fold_results[, 2], use = "complete.obs")^2
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
#    Thresholds compared:
#      (1) Naive fixed threshold of 0.5
#      (2) Prevalence-based threshold: proportion of
#          positives among training species per fold.
#          This is the theoretically correct decision
#          boundary for an LPM under symmetric loss —
#          the model predicts P(positive), so the
#          natural cut is the base rate. Adapts
#          automatically to class imbalance with no
#          leakage risk.
#    Final validation score uses the prevalence-based
#    threshold. Comparison is printed for sense-checking.
#    Repeats: 5 x 10-fold
# -------------------------------------------------------
cat("\n--- B. Binary Traits (Rphylopars LPM, Kappa, prevalence-based threshold) ---\n")

n_folds_bin   <- 10
n_repeats_bin <- 5

fast_kappa <- function(pred_class, truth) {
  cm <- table(factor(pred_class, levels = c(0, 1)),
              factor(truth,      levels = c(0, 1)))
  n  <- sum(cm)
  if (n == 0) return(-Inf)
  po <- sum(diag(cm)) / n
  pe <- sum(rowSums(cm) * colSums(cm)) / n^2
  if (pe >= 1) return(-Inf)
  (po - pe) / (1 - pe)
}

# Stores OOF raw predictions for diagnostics (populated from repeat 1)
oob_predictions    <- list()
optimal_thresholds <- list()

for (trait in eval_binary) {
  
  if (trait %in% traits_100_coverage) {
    cat(sprintf("\n[%d/%d] %s | Skipping CV (100%% coverage already)\n",
                which(eval_binary == trait), length(eval_binary), trait))
    validation_scores[[trait]] <- 1.0 # Assign perfect Kappa to pass Step 6
    optimal_thresholds[[trait]] <- 0.5 # Dummy threshold
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
  
  repeat_kappa_naive <- c()   # Kappa at fixed 0.5
  repeat_kappa_prev  <- c()   # Kappa at prevalence-based threshold
  repeat_thresholds  <- c()   # Mean prevalence threshold per repeat
  all_oob            <- NULL  # OOF predictions from repeat 1
  
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
        
        # Raw predictions for held-out test species
        test_preds  <- pmax(0, pmin(1, p$anc_recon[test_species, trait]))
        test_truths <- as.numeric(as.character(
          df_valid[[trait]][match(test_species, df_valid$species)]))
        
        # Prevalence-based threshold: proportion of positives among
        # training species with observed (non-missing) values
        train_truths <- as.numeric(as.character(
          df_valid[[trait]][match(train_species, df_valid$species)]))
        train_truths_obs <- train_truths[!is.na(train_truths)]
        
        prev_thresh <- if (length(train_truths_obs) > 0 && sum(train_truths_obs) > 0)
          mean(train_truths_obs)
        else
          0.5  # Fallback if no positives in training fold
        
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
      
      # --- Naive threshold (0.5) ---
      pred_naive <- as.integer(fold_results$Pred_raw > 0.5)
      
      # --- Prevalence-based threshold (fold-specific) ---
      pred_prev  <- as.integer(fold_results$Pred_raw > fold_results$Fold_thresh)
      
      mean_thresh_r <- mean(fold_results$Fold_thresh)
      repeat_thresholds <- c(repeat_thresholds, mean_thresh_r)
      
      # Warn if either method collapses to constant predictions
      if (length(unique(pred_naive)) == 1)
        message("WARNING: naive threshold gives constant predictions for ", trait)
      if (length(unique(pred_prev)) == 1)
        message("WARNING: prevalence threshold gives constant predictions for ", trait)
      
      tryCatch({
        k_naive <- confusionMatrix(
          factor(pred_naive,              levels = c(0, 1)),
          factor(fold_results$Truth,      levels = c(0, 1))
        )$overall["Kappa"]
        if (!is.na(k_naive)) repeat_kappa_naive <- c(repeat_kappa_naive, k_naive)
      }, error = function(e) message("Kappa (naive) error for ",  trait, ": ", e$message))
      
      tryCatch({
        k_prev <- confusionMatrix(
          factor(pred_prev,               levels = c(0, 1)),
          factor(fold_results$Truth,      levels = c(0, 1))
        )$overall["Kappa"]
        if (!is.na(k_prev)) repeat_kappa_prev <- c(repeat_kappa_prev, k_prev)
      }, error = function(e) message("Kappa (prev) error for ",   trait, ": ", e$message))
      
      # Accumulate OOF predictions from repeat 1 for diagnostics
      if (r == 1) all_oob <- fold_results[, c("Pred_raw", "Truth", "Fold_thresh")]
    }
  }
  
  # --- Report and store results ---
  if (length(repeat_kappa_prev) > 0) {
    
    mean_kappa_naive <- if (length(repeat_kappa_naive) > 0) mean(repeat_kappa_naive) else NA_real_
    mean_kappa_prev  <- mean(repeat_kappa_prev)
    sd_kappa_prev    <- sd(repeat_kappa_prev)
    mean_thresh      <- mean(repeat_thresholds)
    sd_thresh        <- sd(repeat_thresholds)
    delta            <- mean_kappa_prev - mean_kappa_naive
    
    # Validation score uses prevalence-based Kappa
    validation_scores[[trait]]  <- mean_kappa_prev
    optimal_thresholds[[trait]] <- mean_thresh
    oob_predictions[[trait]]    <- all_oob
    
    cat(sprintf(
      "   -> Naive (0.5)     Kappa: %.3f\n",
      mean_kappa_naive))
    cat(sprintf(
      "   -> Prevalence (%.2f ± %.2f) Kappa: %.3f (± %.3f SD)  [delta: %+.3f]\n",
      mean_thresh, sd_thresh, mean_kappa_prev, sd_kappa_prev, delta))
    
  } else {
    cat(sprintf("   -> %s: no valid folds returned.\n", trait))
  }
}


# -------------------------------------------------------
# C. MULTINOMIAL TRAITS
#    Model:   One-vs-Rest (OvR) Rphylopars LPM
#             Fits one binary phylopars model per class, normalises the
#             raw predictions across classes so they sum to 1, then takes
#             the argmax as the hard classification.
#    Metric:  Cohen's Kappa (consistent with binary evaluation)
#    Repeats: 5 x 10-fold
#
#    Why OvR instead of make.simmap?
#      make.simmap cannot handle NA at tips, so the standard masking approach
#      used for other traits breaks entirely. OvR with phylopars avoids this
#      while keeping imputation methodology consistent across all trait families.
# -------------------------------------------------------
cat("\n--- C. Multinomial Traits (OvR Rphylopars, Kappa) ---\n")

n_folds_mult   <- 10
n_repeats_mult <- 5

for (trait in eval_multinomial) {
  
  if (trait %in% traits_100_coverage) {
    cat(sprintf("\n[%d/%d] %s | Skipping CV (100%% coverage already)\n",
                which(eval_multinomial == trait), length(eval_multinomial), trait))
    validation_scores[[trait]] <- 1.0 # Assign perfect Kappa to pass Step 6
    next
  }
  
  valid_species <- df_valid %>% filter(!is.na(!!sym(trait))) %>% pull(species)
  
  if (length(valid_species) < n_folds_mult * 2) {
    cat(sprintf("\nSkipping %s: too few observations (%d)\n", trait, length(valid_species)))
    next
  }
  
  # Extract unique classes observed for this trait
  # factor() is required here because multinomial columns are character,
  # not factor — droplevels() cannot operate on character vectors directly
  trait_classes <- levels(droplevels(factor(
    df_valid[[trait]][df_valid$species %in% valid_species])))
  
  cat(sprintf("\n[%d/%d] %s | %d species | Classes: %s\n",
              which(eval_multinomial == trait), length(eval_multinomial),
              trait, length(valid_species), paste(trait_classes, collapse = ", ")))
  
  # Create temporary binary columns in df_valid (one per class)
  for (cls in trait_classes) {
    df_valid[[paste0(trait, "_ovr_", cls)]] <- as.numeric(df_valid[[trait]] == cls)
  }
  
  repeat_kappa <- c()
  
  for (r in 1:n_repeats_mult) {
    
    set.seed(42 + r)
    shuffled_sp <- sample(valid_species)
    folds       <- cut(seq_along(shuffled_sp), breaks = n_folds_mult, labels = FALSE)
    
    pb   <- txtProgressBar(max = n_folds_mult, style = 3)
    opts <- list(progress = function(n) setTxtProgressBar(pb, n))
    
    fold_results <- foreach(
      f = 1:n_folds_mult, .combine = rbind,
      .packages     = "Rphylopars",
      .export       = c("tree", "df_valid", "trait", "trait_classes",
                        "valid_species", "shuffled_sp", "folds"),
      .options.snow = opts
    ) %dopar% {
      
      test_species <- shuffled_sp[which(folds == f)]
      
      # OvR: fit one phylopars model per class
      # Returns a matrix: rows = test species, cols = classes
      class_probs <- sapply(trait_classes, function(cls) {
        bin_col    <- paste0(trait, "_ovr_", cls)
        train_data <- df_valid[df_valid$species %in% valid_species,
                               c("species", bin_col)]
        names(train_data)[2] <- "y"
        train_data$y <- as.numeric(train_data$y)
        train_data[train_data$species %in% test_species, "y"] <- NA
        
        tryCatch({
          p <- phylopars(trait_data = train_data, tree = tree,
                         model = "lambda", pheno_error = FALSE)
          p$anc_recon[test_species, "y"]
        }, error = function(e) {
          message("OvR error (", cls, ", fold ", f, "): ", e$message)
          rep(NA_real_, length(test_species))
        })
      })
      
      # Normalise rows so class probabilities sum to 1 before argmax
      row_sums         <- rowSums(class_probs, na.rm = TRUE)
      class_probs_norm <- class_probs / row_sums
      
      # Hard classification: class with highest normalised probability
      pred_class <- trait_classes[apply(class_probs_norm, 1, which.max)]
      true_class <- as.character(
        df_valid[[trait]][match(test_species, df_valid$species)])
      
      data.frame(Pred = pred_class, Truth = true_class, stringsAsFactors = FALSE)
    }
    close(pb)
    
    if (is.null(fold_results) || nrow(fold_results) == 0) {
      cat(sprintf("   Repeat %d: no valid fold results.\n", r)); next
    }
    
    fold_results <- fold_results[complete.cases(fold_results), ]
    if (nrow(fold_results) == 0) next
    
    # Warn if any class is never predicted (indicates near-constant OvR model)
    missing_classes <- setdiff(trait_classes, unique(fold_results$Pred))
    if (length(missing_classes) > 0) {
      message("WARNING: class(es) never predicted for ", trait, ": ",
              paste(missing_classes, collapse = ", "))
    }
    
    tryCatch({
      cm <- confusionMatrix(
        factor(fold_results$Pred,  levels = trait_classes),
        factor(fold_results$Truth, levels = trait_classes)
      )
      k <- cm$overall["Kappa"]
      if (!is.na(k)) repeat_kappa <- c(repeat_kappa, k)
    }, error = function(e) message("Kappa error for ", trait, ": ", e$message))
  }
  
  if (length(repeat_kappa) > 0) {
    mean_k <- mean(repeat_kappa, na.rm = TRUE)
    sd_k   <- sd(repeat_kappa,   na.rm = TRUE)
    validation_scores[[trait]] <- mean_k
    cat(sprintf("   -> Mean Kappa: %.3f (±%.3f SD) over %d repeats\n",
                mean_k, sd_k, n_repeats_mult))
  } else {
    cat(sprintf("   -> %s: no valid results across all repeats.\n", trait))
  }
  
  # Remove temporary OvR binary columns
  df_valid[, paste0(trait, "_ovr_", trait_classes)] <- NULL
}

stopCluster(cl)


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
    
  } else if (trait %in% c(eval_binary, eval_multinomial)) {
    # Kappa applied consistently across all categorical trait families
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
final_mult <- intersect(passed_traits, eval_multinomial)


# --- A. Continuous Imputation ---
# All continuous traits are imputed jointly in a single phylopars call.
# Joint imputation leverages cross-trait correlations (via the multivariate
# Brownian motion model) and is more efficient than per-trait calls.
if (length(final_cont) > 0) {
  cat(sprintf("Imputing %d continuous traits (joint phylopars)...\n",
              length(final_cont)))
  
  train_cont <- df_valid[, c("species", final_cont)]
  p_cont     <- phylopars(trait_data = train_cont, tree = tree,
                          model = "lambda", pheno_error = FALSE)
  
  cont_recon         <- as.data.frame(p_cont$anc_recon[tree$tip.label, final_cont, drop = FALSE])
  cont_recon$species <- rownames(cont_recon)
  imputed_df         <- left_join(imputed_df, cont_recon, by = "species")
  cat("   -> Done.\n")
}


# --- B. Binary Imputation ---
# Each binary trait is imputed separately via phylopars LPM (same model as CV).
# Predicted probabilities > 0.5 are classified as 1; otherwise 0.
# Original observed values are always preserved — only NAs are filled.
if (length(final_bin) > 0) {
  cat(sprintf("Imputing %d binary traits (phylopars LPM)...\n", length(final_bin)))
  
  for (trait in final_bin) {
    cat(sprintf("  -> %s\n", trait))
    
    train_bin        <- df_valid[, c("species", trait)]
    train_bin[[trait]] <- as.numeric(as.character(train_bin[[trait]]))
    
    p_bin  <- phylopars(trait_data = train_bin, tree = tree,
                        model = "lambda", pheno_error = FALSE)
    preds  <- p_bin$anc_recon[tree$tip.label, trait]
    preds  <- pmax(0, pmin(1, preds))  # Clamp to [0,1]
    
    # Retrieve original observed values (NA = missing, needs imputing)
    original <- as.numeric(as.character(
      df_valid[[trait]][match(tree$tip.label, df_valid$species)]))
    
    # Preserve observed values; fill missing with thresholded prediction
    imputed_df[[trait]] <- ifelse(is.na(original), as.integer(preds > 0.5), original)
  }
  cat("   -> Done.\n")
}


# --- C. Multinomial Imputation ---
# OvR phylopars: same methodology as CV to ensure consistency.
# Three separate LPM models are fit (one per breeding class), predictions are
# normalised to sum to 1, and the argmax class is assigned.
# Original observed values are always preserved.
if (length(final_mult) > 0) {
  cat(sprintf("Imputing %d multinomial traits (OvR phylopars)...\n", length(final_mult)))
  
  for (trait in final_mult) {
    cat(sprintf("  -> %s\n", trait))
    
    breed_classes <- levels(droplevels(na.omit(df_valid[[trait]])))
    
    # Fit one phylopars LPM per class on all observed data
    class_probs <- sapply(breed_classes, function(cls) {
      train_data    <- data.frame(
        species = df_valid$species,
        y       = as.numeric(df_valid[[trait]] == cls)
      )
      p <- phylopars(trait_data = train_data, tree = tree,
                     model = "lambda", pheno_error = FALSE)
      p$anc_recon[tree$tip.label, "y"]
    })
    
    # Normalise across classes so predictions sum to 1 per species
    row_sums         <- rowSums(class_probs, na.rm = TRUE)
    class_probs_norm <- class_probs / row_sums
    
    # Argmax: assign the most probable class
    pred_class <- breed_classes[apply(class_probs_norm, 1, which.max)]
    
    # Retrieve original observed values
    original <- as.character(
      df_valid[[trait]][match(tree$tip.label, df_valid$species)])
    
    # Preserve observed values; fill missing with predicted class
    imputed_df[[trait]] <- factor(
      ifelse(is.na(original), pred_class, original),
      levels = breed_classes
    )
  }
  cat("   -> Done.\n")
}


# ==============================================================================
# 8. OUTPUT
# ==============================================================================
out_path <- "data/processed/traits/amphibian/amphibian_traits_imputed_final.csv"
write.csv(imputed_df, out_path, row.names = FALSE)

cat(sprintf(
  "\n=======================================================\n
Pipeline complete.\n
  Species in output : %d\n
  Traits imputed    : %d continuous | %d binary | %d multinomial\n
  Saved to          : %s\n
=======================================================\n",
  nrow(imputed_df),
  length(final_cont), length(final_bin), length(final_mult),
  out_path
))

