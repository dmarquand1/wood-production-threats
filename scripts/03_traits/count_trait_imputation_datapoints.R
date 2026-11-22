.libPaths(c(normalizePath(".r-lib", winslash = "/", mustWork = FALSE), .libPaths()))

library(dplyr)
library(stringr)
library(ape)

clean_join_name <- function(x) {
  x %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_replace_all("\u00A0", " ") %>%
    str_squish()
}

zscale <- function(x) as.numeric(scale(x))

fix_binary_nas <- function(df, cols) {
  cols <- intersect(cols, names(df))
  if (!length(cols)) return(df)

  all_na <- apply(df[, cols, drop = FALSE], 1, function(x) all(is.na(x)))
  for (col in cols) {
    df[[col]] <- ifelse(!all_na & is.na(df[[col]]), 0, df[[col]])
  }
  df
}

collapse_to_multinomial <- function(df, cols, col_name) {
  cols <- intersect(cols, names(df))
  df[[col_name]] <- apply(df[, cols, drop = FALSE], 1, function(row) {
    if (all(is.na(row))) return(NA_character_)
    if (any(row == 1, na.rm = TRUE)) return(cols[which(row == 1)[1]])
    NA_character_
  })
  df[[col_name]] <- as.factor(df[[col_name]])
  df
}

make_binary_dummies <- function(df, col, prefix) {
  vals <- sort(unique(trimws(unlist(strsplit(na.omit(df[[col]]), "/", fixed = TRUE)))))
  vals <- vals[vals != ""]

  for (val in vals) {
    nm <- paste0(prefix, "_", gsub("[^A-Za-z0-9]", "_", trimws(val)))
    df[[nm]] <- ifelse(
      is.na(df[[col]]),
      NA_real_,
      ifelse(grepl(val, df[[col]], fixed = TRUE), 1, 0)
    )
  }
  df
}

read_hab <- function(tx) {
  readRDS(file.path("data/processed/assessments", paste0(tx, "_assessment_habitat_summary.rds"))) %>%
    mutate(scientific_name = clean_join_name(scientific_name))
}

get_model_traits <- function(tx, imp_file) {
  model <- read.csv(
    file.path("data/processed/traits", tx, paste0(tx, "_trait_model_final_data.csv")),
    check.names = FALSE
  )
  imp <- read.csv(file.path("data/processed/traits", tx, imp_file), check.names = FALSE)
  intersect(setdiff(names(imp), "species"), names(model))
}

count_missing <- function(tx, pre, key_col, imp_file) {
  model <- read.csv(
    file.path("data/processed/traits", tx, paste0(tx, "_trait_model_final_data.csv")),
    check.names = FALSE
  )
  traits <- get_model_traits(tx, imp_file)

  miss <- pre %>% select(species = all_of(key_col), all_of(traits))
  miss_map <- miss[match(model$phylo_name, miss$species), traits, drop = FALSE]

  by_trait <- data.frame(
    taxon = tx,
    trait = traits,
    model_rows = nrow(model),
    imputed_datapoints = colSums(is.na(miss_map)),
    stringsAsFactors = FALSE
  )

  list(
    summary = data.frame(
      taxon = tx,
      model_rows = nrow(model),
      distinct_species = length(unique(model$phylo_name)),
      model_trait_columns = length(traits),
      total_model_trait_cells = nrow(model) * length(traits),
      imputed_datapoints = sum(by_trait$imputed_datapoints),
      pct_imputed = round(100 * sum(by_trait$imputed_datapoints) / (nrow(model) * length(traits)), 2)
    ),
    by_trait = by_trait
  )
}

prep_amphibian <- function() {
  tree <- read.tree("data/processed/traits/amphibian/final_amphibian_tree_pruned.tre")
  df <- read.csv("data/processed/traits/amphibian/amphibian_trait_data_matched.csv", stringsAsFactors = FALSE) %>%
    filter(Final_Tree_Tip %in% tree$tip.label)

  df <- df %>%
    mutate(IUCN_Name = clean_join_name(IUCN_Name)) %>%
    left_join(
      read_hab("amphibian") %>% select(scientific_name, Forest:habitat_breadth),
      by = c("IUCN_Name" = "scientific_name")
    ) %>%
    fix_binary_nas(c("Leaves", "Flowers", "Seeds", "Fruits", "Arthro", "Vert")) %>%
    fix_binary_nas(c("Diu", "Noc", "Crepu")) %>%
    fix_binary_nas(c("Wet_warm", "Wet_cold", "Dry_warm", "Dry_cold")) %>%
    collapse_to_multinomial(c("Dir", "Lar", "Viv"), "Breeding_Cat") %>%
    mutate(Body_Size_log = zscale(log10(Body_size_mm)))

  df
}

prep_bird <- function() {
  tree <- read.tree("data/processed/traits/bird/final_bird_tree_pruned.tre")
  df <- read.csv("data/processed/traits/bird/bird_trait_data_matched.csv", stringsAsFactors = FALSE, check.names = FALSE)

  if ("phylo_match" %in% names(df)) {
    df <- df %>% rename(phylo_name = phylo_match)
  } else if ("Tree_Name" %in% names(df)) {
    df <- df %>% rename(phylo_name = Tree_Name)
  } else if ("IOC.World.Bird.List..v15.1." %in% names(df)) {
    df <- df %>% rename(phylo_name = `IOC.World.Bird.List..v15.1.`)
  }

  iucn_habitat_cols <- c(
    "Forest", "Savanna", "Shrubland", "Grassland", "Wetlands",
    "Rocky", "Caves", "Desert", "Marine", "Artificial", "Other", "Unknown",
    "habitat_breadth"
  )

  df <- df %>%
    mutate(
      phylo_name = str_replace_all(phylo_name, " ", "_"),
      spatial_name = Trait_Name
    ) %>%
    filter(phylo_name %in% tree$tip.label) %>%
    mutate(join_name = clean_join_name(Trait_Name)) %>%
    left_join(
      read_hab("bird") %>% select(scientific_name, all_of(iucn_habitat_cols)),
      by = c("join_name" = "scientific_name")
    ) %>%
    select(-join_name) %>%
    mutate(
      Body_Mass_log = zscale(log(Average.Mass)),
      Elevational_Range_log = zscale(log(Elevational.Range + 1)),
      ESI_scaled = zscale(ESI),
      Clutch_Size_log = zscale(log((Clutch_Min + Clutch_Max) / 2)),
      HWI_log = zscale(log(HWI)),
      Habitat_Breadth_log = zscale(log(habitat_breadth)),
      Island_Status = ISL,
      Is_Flightless = ifelse(Flightlessness %in% c("Yes", "Partial"), 1, 0),
      Is_Migratory = case_when(
        is.na(Mig) & is.na(Alt) & is.na(Irreg) & is.na(Disp) & is.na(Sed) ~ NA_real_,
        Sed == 1 ~ 0,
        TRUE ~ 1
      ),
      Is_Tree_Nester = case_when(is.na(Nest_SBS) ~ NA_real_, grepl("T|S", Nest_SBS) ~ 1, TRUE ~ 0),
      Is_Stump_Nester = case_when(is.na(Nest_SBS) ~ NA_real_, grepl("C", Nest_SBS) ~ 1, TRUE ~ 0),
      Is_Ground_Nester = case_when(is.na(Nest_SBS) ~ NA_real_, grepl("G", Nest_SBS) ~ 1, TRUE ~ 0)
    )

  diet_wt_cols <- c("IN.Wt", "FR.Wt", "NE.Wt", "SE.Wt", "VE.Wt", "FI.Wt", "SC.Wt", "PL.Wt", "MS.Wt")
  diet_binary_cols <- c(
    "Diet_Invertebrate", "Diet_Fruit", "Diet_Nectar", "Diet_Seed",
    "Diet_Vertebrate", "Diet_Fish", "Diet_Scavenger", "Diet_PlantMaterial",
    "Diet_Misc"
  )
  all_diet_na <- rowSums(!is.na(df[, diet_wt_cols])) == 0

  for (i in seq_along(diet_wt_cols)) {
    wt_col <- diet_wt_cols[i]
    bin_col <- diet_binary_cols[i]
    df[[bin_col]] <- case_when(
      all_diet_na ~ NA_real_,
      is.na(df[[wt_col]]) ~ 0,
      df[[wt_col]] > 0 ~ 1,
      TRUE ~ 0
    )
  }

  df %>%
    group_by(phylo_name, spatial_name) %>%
    summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE)), .groups = "drop")
}

prep_fish <- function() {
  tree <- read.tree("data/processed/traits/fish/final_fish_tree_pruned.tre")
  num <- function(x) suppressWarnings(as.numeric(x))

  read.csv("data/processed/traits/fish/fish_trait_data_matched.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    filter(Final_Tree_Tip %in% tree$tip.label) %>%
    mutate(IUCN_Name = clean_join_name(IUCN_Name)) %>%
    left_join(
      read_hab("fish") %>% select(scientific_name, Forest:habitat_breadth),
      by = c("IUCN_Name" = "scientific_name")
    ) %>%
    mutate(
      species_Length_num = num(species_Length),
      estimate_MaxLengthTL_num = num(estimate_MaxLengthTL),
      estimate_MaxLengthSL_num = num(estimate_MaxLengthSL),
      estimate_Troph_num = num(estimate_Troph),
      estimate_PredPreyRatioMax_num = num(estimate_PredPreyRatioMax),
      MaxLength_cm = coalesce(estimate_MaxLengthTL_num, estimate_MaxLengthSL_num, species_Length_num),
      MaxLength_Log = ifelse(!is.na(MaxLength_cm) & MaxLength_cm > 0, zscale(log10(MaxLength_cm)), NA_real_),
      Trophic_Level_sc = ifelse(!is.na(estimate_Troph_num), zscale(estimate_Troph_num), NA_real_),
      PredPreyRatioMax_log_sc = ifelse(
        !is.na(estimate_PredPreyRatioMax_num) & estimate_PredPreyRatioMax_num > 0,
        zscale(log10(estimate_PredPreyRatioMax_num)),
        NA_real_
      ),
      Fresh = case_when(species_Fresh %in% c(0, 1) ~ as.integer(species_Fresh), TRUE ~ NA_integer_),
      Brack = case_when(species_Brack %in% c(0, 1) ~ as.integer(species_Brack), TRUE ~ NA_integer_),
      Salt = case_when(species_Saltwater %in% c(0, 1) ~ as.integer(species_Saltwater), TRUE ~ NA_integer_),
      Salinity_Breadth = case_when(!is.na(Fresh) & !is.na(Brack) & !is.na(Salt) ~ Fresh + Brack + Salt, TRUE ~ NA_real_),
      Salinity_Breadth_sc = ifelse(!is.na(Salinity_Breadth) & Salinity_Breadth > 0, zscale(Salinity_Breadth), NA_real_),
      AirBreathing_Bin = case_when(
        species_AirBreathing %in% c("Facultative", "FacultativeGenus", "FacultativeObligate", "Obligate", "ObligateGenus") ~ 1L,
        species_AirBreathing == "WaterAssumed" ~ 0L,
        TRUE ~ NA_integer_
      ),
      Electrogenic_Bin = case_when(
        tolower(species_Electrogenic) %in% c("weakly discharging", "strongly discharging", "electrosensing only") ~ 1L,
        tolower(species_Electrogenic) == "no special ability" ~ 0L,
        TRUE ~ NA_integer_
      ),
      Dangerous_Bin = case_when(
        tolower(species_Dangerous) %in% c("harmless", "harmlessgenus") ~ 0L,
        is.na(species_Dangerous) | species_Dangerous == "" ~ NA_integer_,
        TRUE ~ 1L
      ),
      DemersPelag_Cat = as.factor(species_DemersPelag),
      EnvTemp_Cat = as.factor(stocks_EnvTemp),
      AirBreathing_Cat = as.factor(species_AirBreathing),
      Electrogenic_Cat = as.factor(species_Electrogenic)
    )
}

prep_mammal <- function() {
  tree <- read.tree("data/processed/traits/mammal/final_mammal_tree_pruned.tre")
  df <- read.csv("data/processed/traits/mammal/mammal_trait_data_matched.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    rename(phylo_name = phylo_match, spatial_name = IUCN_Name) %>%
    mutate(phylo_name = str_replace_all(phylo_name, " ", "_")) %>%
    filter(phylo_name %in% tree$tip.label)

  for (cc in grep("^det_", names(df), value = TRUE)) {
    df[[paste0("binary_", cc)]] <- as.integer(df[[cc]] > 0)
  }
  for (cc in grep("^dphy_", names(df), value = TRUE)) {
    df[[paste0("binary_", cc)]] <- as.integer(df[[cc]] > 0)
  }

  df %>%
    mutate(
      activity_cycle = case_when(
        activity_cycle == 1 ~ "Nocturnal",
        activity_cycle == 2 ~ "Cathemeral_Crepuscular",
        activity_cycle == 3 ~ "Diurnal",
        TRUE ~ NA_character_
      ),
      spatial_name = clean_join_name(spatial_name)
    ) %>%
    left_join(
      read_hab("mammal") %>% select(scientific_name, Forest:habitat_breadth),
      by = c("spatial_name" = "scientific_name")
    ) %>%
    mutate(
      Body_Mass_Log = zscale(log10(adult_mass_g)),
      Litter_Size_Log = zscale(log10(litter_size_n)),
      Dispersal_Log = zscale(log10(dispersal_km)),
      Habitat_Breadth = zscale(habitat_breadth),
      Diet_Breadth = zscale(det_diet_breadth_n),
      ESI = zscale(log10(100 / (habitat_breadth_n * habitat_breadth)))
    )
}

prep_reptile <- function() {
  tree <- read.tree("data/processed/traits/reptile/final_reptile_tree_pruned.tre")

  read.csv("data/processed/traits/reptile/reptile_trait_data_matched.csv", stringsAsFactors = FALSE, check.names = FALSE) %>%
    filter(Final_Tree_Tip %in% tree$tip.label) %>%
    mutate(IUCN_Name = clean_join_name(IUCN_Name)) %>%
    left_join(
      read_hab("reptile") %>% select(scientific_name, Forest:habitat_breadth),
      by = c("IUCN_Name" = "scientific_name")
    ) %>%
    make_binary_dummies("Dorsal.colour", "DorsalColour") %>%
    make_binary_dummies("Dorsal.pattern", "DorsalPattern") %>%
    make_binary_dummies("Microhabitat", "Microhabitat") %>%
    mutate(
      Max_Body_Mass_log = zscale(log10(Maximum.body.mass..g.)),
      Max_SVL_log = zscale(log10(Maximum.length...SVL...mm..straight.carapace.length.for.turtles...SCL...mm.)),
      Habitat_Breadth = zscale(log(habitat_breadth)),
      Insular_or_Endemic = case_when(Insular.endemic..yes.or.no. == "Yes" ~ 1, Insular.endemic..yes.or.no. == "No" ~ 0, TRUE ~ NA_real_),
      Venomous = case_when(Venomous..yes.or.no. == "Yes" ~ 1, Venomous..yes.or.no. == "No" ~ 0, TRUE ~ NA_real_)
    )
}

results <- list(
  count_missing("amphibian", prep_amphibian(), "Final_Tree_Tip", "amphibian_traits_imputed_final.csv"),
  count_missing("bird", prep_bird(), "phylo_name", "bird_traits_imputed_ready.csv"),
  count_missing("fish", prep_fish(), "Final_Tree_Tip", "fish_traits_imputed_final.csv"),
  count_missing("mammal", prep_mammal(), "phylo_name", "imputed_mammal_traits_ready.csv"),
  count_missing("reptile", prep_reptile(), "Final_Tree_Tip", "imputed_reptile_traits_ready.csv")
)

summary <- bind_rows(lapply(results, `[[`, "summary"))
by_trait <- bind_rows(lapply(results, `[[`, "by_trait"))

dir.create("outputs/traits/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(summary, "outputs/traits/tables/trait_imputation_datapoints_summary.csv", row.names = FALSE)
write.csv(by_trait, "outputs/traits/tables/trait_imputation_datapoints_by_trait.csv", row.names = FALSE)

print(summary)
cat("\nTOTAL\n")
print(summary %>%
  summarise(
    model_rows = sum(model_rows),
    distinct_species = sum(distinct_species),
    model_trait_columns = sum(model_trait_columns),
    total_model_trait_cells = sum(total_model_trait_cells),
    imputed_datapoints = sum(imputed_datapoints),
    pct_imputed = round(100 * imputed_datapoints / total_model_trait_cells, 2)
  ))

cat("\nBY TRAIT, NONZERO\n")
print.data.frame(by_trait %>% filter(imputed_datapoints > 0) %>% arrange(taxon, desc(imputed_datapoints)))
