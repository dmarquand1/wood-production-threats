# =========================================================
# Multi-taxon radial phylogram for wood-production threat
#
# This script builds the same circular tree plot for each of the five
# vertebrate taxa. It can be run on one taxon or on all taxa at once.
#
# Outputs:
#   - data/processed/traits/<taxon>/<taxon>_full_tree.tre
#   - data/processed/traits/<taxon>/<taxon>_iucn_tree_matched.csv
#   - outputs/extent/figures/<taxon>_radial_phylogeny.png
# =========================================================

rm(list = ls())
gc()

# =========================================================
# Package setup
# =========================================================
required_cran <- c(
  "ape", "ggplot2", "dplyr", "readr", "tibble", "stringr",
  "tidyr", "purrr", "rredlist", "here", "phangorn", "stringdist",
  "ggnewscale", "ggrepel", "scales", "patchwork"
)
local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

installed <- rownames(installed.packages(lib.loc = .libPaths()))
missing_required <- setdiff(c(required_cran, "ggtree"), installed)
if (length(missing_required) > 0) {
  stop(
    "Missing required R packages for the phylogeny script: ",
    paste(missing_required, collapse = ", "),
    ". Install them once in your project library before running this script."
  )
}

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(rredlist)
  library(here)
  library(phangorn)
  library(stringdist)
  library(ggnewscale)
  library(ggrepel)
  library(ggtree)
  library(patchwork)
  library(scales)
})

# =========================================================
# Configuration
# =========================================================
taxa_configs <- list(
  amphibian = list(
    taxon = "amphibian",
    iucn_class = "Amphibia",
    raw_tree_file = here::here("data", "raw", "trees", "amphibiantree", "output.nex"),
    tree_format = "nexus",
    raw_iucn_cache_file = here::here("data", "raw", "iucn", "amphib_iucn_assessment_list.rds"),
    processed_dir = here::here("data", "processed", "traits", "amphibian"),
    output_file = here::here("outputs", "extent", "figures", "amphibian_radial_phylogeny.png"),
    manual_patch_file = NULL
  ),
  bird = list(
    taxon = "bird",
    iucn_class = "Aves",
    raw_tree_file = here::here("data", "raw", "trees", "bird", "Bird_Consensus_Tree.tre"),
    tree_format = "tree",
    raw_iucn_cache_file = here::here("data", "raw", "iucn", "iucn_assessment_list.rds"),
    processed_dir = here::here("data", "processed", "traits", "bird"),
    output_file = here::here("outputs", "extent", "figures", "bird_radial_phylogeny.png"),
    manual_patch_file = here::here("data", "raw", "traits", "bird", "Manual_Patch_REQUIRED2.csv")
  ),
  mammal = list(
    taxon = "mammal",
    iucn_class = "Mammalia",
    raw_tree_file = here::here("data", "raw", "trees", "mammaltree", "mammaltree.nex"),
    tree_format = "nexus",
    raw_iucn_cache_file = here::here("data", "raw", "iucn", "mam_iucn_assessment_list.rds"),
    processed_dir = here::here("data", "processed", "traits", "mammal"),
    output_file = here::here("outputs", "extent", "figures", "mammal_radial_phylogeny.png"),
    manual_patch_file = NULL
  ),
  reptile = list(
    taxon = "reptile",
    iucn_class = "Reptilia",
    raw_tree_file = here::here("data", "raw", "trees", "reptiletree", "output.nex"),
    tree_format = "nexus",
    raw_iucn_cache_file = here::here("data", "raw", "iucn", "rept_iucn_assessment_list.rds"),
    processed_dir = here::here("data", "processed", "traits", "reptile"),
    output_file = here::here("outputs", "extent", "figures", "reptile_radial_phylogeny.png"),
    manual_patch_file = NULL
  ),
  fish = list(
    taxon = "fish",
    iucn_class = "Actinopterygii",
    raw_tree_file = here::here("data", "processed", "traits", "fish", "final_fish_tree_pruned.tre"),
    tree_format = "tree",
    raw_iucn_cache_file = here::here("data", "raw", "iucn", "fish_iucn_assessment_list.rds"),
    taxonomy_file = here::here("data", "raw", "trees", "fish", "fishtree_taxonomy.csv"),
    processed_dir = here::here("data", "processed", "traits", "fish"),
    output_file = here::here("outputs", "extent", "figures", "fish_radial_phylogeny.png"),
    manual_patch_file = NULL,
    label_monophyletic_families_only = TRUE
  )
)

args <- commandArgs(trailingOnly = TRUE)
requested_taxa <- if (length(args) == 0) names(taxa_configs) else tolower(args)
if (length(requested_taxa) == 1 && requested_taxa[[1]] == "all") {
  requested_taxa <- names(taxa_configs)
}
unknown_taxa <- setdiff(requested_taxa, names(taxa_configs))
if (length(unknown_taxa) > 0) {
  stop("Unknown taxon name(s): ", paste(unknown_taxa, collapse = ", "))
}

# =========================================================
# Global figure settings
# =========================================================
THREAT_CODE_LOGGING <- "5_3_"
THREAT_CODE_PLANTATION <- "2_2_"

col_none            <- "#d3d3d3"
col_logging_only    <- "#d9be00"
col_plantation_only <- "#0088d9"
col_both            <- "#007b00"

branch_low_col  <- "#2B6CB0"
branch_mid_col  <- "#F3EFCF"
branch_high_col <- "#D73027"
phylo_ring_interval_myr <- 20
phylo_ring_fill_col <- "grey96"

branch_line_size  <- 0.35
ring_offset       <- 0.00
ring_width        <- 0.045
family_ring_gap   <- 1.2
family_ring_width <- 0.018
figure_width      <- 23
figure_height     <- 17
figure_dpi        <- 600
combined_figure_width  <- 36
combined_figure_height <- 24

force_rebuild_full_tree    <- FALSE
force_rebuild_matched_data <- FALSE
force_rebuild_figure       <- FALSE
force_rebuild_d_statistics <- FALSE

# Fritz & Purvis D-statistic settings.
# Reduce permutations while iterating locally if you need faster turnaround.
d_stat_permutations <- 1000L
d_stat_seed         <- 42L

label_families_only_if_monophyletic <- FALSE
family_label_gap <- 15
family_label_stagger_bands <- 4L
family_label_stagger_step <- family_label_gap * 0.7
family_label_barsize <- 0.4
family_label_fontsize <- 3.4

col_family_highlight <- "black"
col_family_other     <- "grey85"

phylo_progress_updates <- 10L

taxon_plot_titles <- c(
  amphibian = "Amphibians",
  bird = "Birds",
  mammal = "Mammals",
  reptile = "Reptiles",
  fish = "Fish"
)

# =========================================================
# Helpers
# =========================================================
clean_names <- function(x) {
  if (is.null(x)) return(NULL)
  x %>%
    str_replace_all("\uFFFD", " ") %>%
    str_replace_all("\u00A0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

to_phylo_style <- function(x) {
  clean_names(x) %>% str_replace_all(" ", "_")
}

to_logical_safe <- function(x) {
  dplyr::case_when(
    x %in% c(TRUE, "TRUE", "True", "T", 1, "1") ~ TRUE,
    x %in% c(FALSE, "FALSE", "False", "F", 0, "0") ~ FALSE,
    is.na(x) ~ NA,
    TRUE ~ NA
  )
}

ensure_parent_dir <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

append_file_stem_suffix <- function(path, suffix) {
  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  if (nzchar(ext)) {
    paste0(stem, suffix, ".", ext)
  } else {
    paste0(path, suffix)
  }
}

format_elapsed_time <- function(start_time) {
  elapsed_secs <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (!is.finite(elapsed_secs) || elapsed_secs < 0) {
    elapsed_secs <- 0
  }
  
  if (elapsed_secs < 60) {
    return(sprintf("%.1fs", elapsed_secs))
  }
  
  elapsed_mins <- floor(elapsed_secs / 60)
  remaining_secs <- round(elapsed_secs %% 60)
  sprintf("%dm %02ds", elapsed_mins, remaining_secs)
}

emit_progress_message <- function(label, current, total, start_time) {
  pct_complete <- if (total > 0) (100 * current / total) else 100
  cat(sprintf(
    "      %s: %d/%d (%.0f%%) elapsed %s\n",
    label,
    current,
    total,
    pct_complete,
    format_elapsed_time(start_time)
  ))
}

should_emit_progress <- function(current, total, updates = phylo_progress_updates) {
  if (total <= 0) {
    return(FALSE)
  }
  
  progress_every <- max(1L, ceiling(total / updates))
  current %% progress_every == 0L || current == total
}

safe_write_csv <- function(x, path) {
  ensure_parent_dir(path)
  
  tryCatch(
    {
      write_csv(x, path)
      path
    },
    error = function(e) {
      fallback_path <- file.path(
        dirname(path),
        paste0(
          tools::file_path_sans_ext(basename(path)),
          "_fallback_",
          format(Sys.time(), "%Y%m%d_%H%M%S"),
          ".csv"
        )
      )
      write_csv(x, fallback_path)
      warning(
        "Could not write to ", path, ". Wrote fallback CSV instead: ",
        fallback_path, ". Original error: ", conditionMessage(e)
      )
      fallback_path
    }
  )
}

read_tree_file <- function(path, format) {
  if (!file.exists(path)) {
    stop("Missing raw tree file: ", path)
  }
  if (identical(tolower(format), "nexus")) {
    trees <- read.nexus(path)
    if (inherits(trees, "multiPhylo")) {
      if (length(trees) == 1) {
        return(trees[[1]])
      }
      return(maxCladeCred(trees))
    }
    if (inherits(trees, "phylo")) {
      return(trees)
    }
    stop("Unexpected object returned while reading nexus tree file: ", path)
  }
  read.tree(path)
}

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA else x[[1]]
}

first_field <- function(x, field) {
  if (is.null(x) || is.null(x[[field]])) return(character(0))
  as.character(x[[field]])
}

extract_synonyms <- function(taxon_info) {
  syn <- taxon_info$synonyms
  if (is.null(syn) || length(syn) == 0) {
    return(character(0))
  }
  
  if (is.data.frame(syn) && all(c("genus_name", "species_name") %in% names(syn))) {
    return(unique(na.omit(paste(syn$genus_name, syn$species_name))))
  }
  
  if (is.list(syn) && !is.data.frame(syn)) {
    genus <- syn$genus_name
    species <- syn$species_name
    if (!is.null(genus) && !is.null(species)) {
      return(unique(na.omit(paste(genus, species))))
    }
  }
  
  character(0)
}

build_iucn_taxonomy_lookup <- function(iucn_raw) {
  rows <- vector("list", length(iucn_raw))
  for (i in seq_along(iucn_raw)) {
    rows[[i]] <- tryCatch(
      {
        rec <- iucn_raw[[i]]
        taxon_info <- rec$taxon
        sci <- clean_names(first_non_na(first_field(taxon_info, "scientific_name")))
        if (is.na(sci) || sci == "") return(NULL)
        tibble(
          IUCN_Name = sci,
          Family_lookup = clean_names(first_non_na(first_field(taxon_info, "family_name"))),
          Order_lookup = clean_names(first_non_na(first_field(taxon_info, "order_name")))
        )
      },
      error = function(e) NULL
    )
  }
  
  bind_rows(purrr::compact(rows)) %>%
    filter(!is.na(IUCN_Name), IUCN_Name != "") %>%
    group_by(IUCN_Name) %>%
    summarise(
      Family_lookup = first_non_na(Family_lookup),
      Order_lookup = first_non_na(Order_lookup),
      .groups = "drop"
    )
}

standardize_match_table <- function(df, iucn_lookup = NULL) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  tip_candidates <- c(
    "phylo_name", "phylo_match", "Tree_Name", "Final_Tree_Tip",
    "Final.Tree.Tip", "matched_tip", "eBird_Name"
  )
  iucn_name_candidates <- c("IUCN_Name", "Trait_Name", "scientific_name", "Species", "Species1")
  family_candidates <- c(
    "Family", "family", "Family.IOC.15.1", "Family.AviList.v1.2025",
    "Family.Clements.v2024b", "Family.HBW.BirdLife.v9.1..2024."
  )
  order_candidates <- c("Order", "order")
  redlist_candidates <- c("RedList_Cat", "RedList_Category", "red_list_category", "category", "X2024.IUCN.Red.List.category")
  
  tip_col <- intersect(tip_candidates, names(df))[1]
  if (is.na(tip_col)) return(NULL)
  if (!all(c("Threat_Logging", "Threat_Plantation") %in% names(df))) return(NULL)
  
  iucn_name_col <- intersect(iucn_name_candidates, names(df))[1]
  family_col <- intersect(family_candidates, names(df))[1]
  order_col <- intersect(order_candidates, names(df))[1]
  redlist_col <- intersect(redlist_candidates, names(df))[1]
  
  match_type_vec <- if ("Match_Type" %in% names(df)) df[["Match_Type"]] else rep("Processed Match", nrow(df))
  iucn_name_vec <- if (!is.na(iucn_name_col)) df[[iucn_name_col]] else rep(NA_character_, nrow(df))
  redlist_vec <- if (!is.na(redlist_col)) df[[redlist_col]] else rep(NA_character_, nrow(df))
  
  out <- tibble(
    phylo_name = to_phylo_style(df[[tip_col]]),
    IUCN_Name = clean_names(iucn_name_vec),
    Family = if (!is.na(family_col)) clean_names(df[[family_col]]) else rep(NA_character_, nrow(df)),
    Order = if (!is.na(order_col)) clean_names(df[[order_col]]) else rep(NA_character_, nrow(df)),
    RedList_Cat = redlist_vec,
    Threat_Logging = to_logical_safe(df[["Threat_Logging"]]),
    Threat_Plantation = to_logical_safe(df[["Threat_Plantation"]]),
    Match_Type = match_type_vec
  ) %>%
    filter(!is.na(phylo_name), phylo_name != "") %>%
    group_by(phylo_name) %>%
    arrange(
      factor(
        Match_Type,
        levels = c(
          "Direct", "Synonym", "Manual Patch", "Synonym_IUCN",
          "Synonym_FishBase", "Synonym_FishBase_via_IUCN",
          "Processed Match", "Fuzzy"
        )
      ),
      IUCN_Name
    ) %>%
    summarise(across(everything(), first_non_na), .groups = "drop") %>%
    mutate(
      Threat_Logging = replace_na(Threat_Logging, FALSE),
      Threat_Plantation = replace_na(Threat_Plantation, FALSE)
    )
  
  if (!is.null(iucn_lookup) && nrow(iucn_lookup) > 0) {
    out <- out %>%
      left_join(iucn_lookup, by = "IUCN_Name") %>%
      mutate(
        Family = coalesce(Family, Family_lookup),
        Order = coalesce(Order, Order_lookup)
      ) %>%
      select(-Family_lookup, -Order_lookup)
  }
  
  out
}

extract_assessment_row <- function(rec, taxon) {
  taxon_info <- rec$taxon
  threats_df <- rec$threats
  systems_df <- rec$systems
  red_list_info <- rec$red_list_category
  assessment_id <- if (!is.null(rec$assessment_id)) rec$assessment_id else NA_integer_
  
  scientific_name <- first_non_na(first_field(taxon_info, "scientific_name"))
  family_name <- first_non_na(first_field(taxon_info, "family_name"))
  order_name <- first_non_na(first_field(taxon_info, "order_name"))
  category <- first_non_na(first_field(red_list_info, "code"))
  
  threat_codes <- first_field(threats_df, "code")
  threat_timing <- first_field(threats_df, "timing")
  system_codes <- first_field(systems_df, "code")
  
  keep_assessment <- if (taxon == "fish") {
    any(system_codes == "1", na.rm = TRUE)
  } else {
    !any(system_codes %in% "2", na.rm = TRUE)
  }
  
  if (!isTRUE(keep_assessment)) {
    return(NULL)
  }
  
  tibble(
    assessment_id = assessment_id,
    scientific_name = scientific_name,
    Family = family_name,
    Order = order_name,
    category = category,
    Threat_Logging = if (length(threat_codes) == 0) {
      FALSE
    } else {
      any(
        str_detect(threat_codes, paste0("^", THREAT_CODE_LOGGING)) & threat_timing == "Ongoing",
        na.rm = TRUE
      )
    },
    Threat_Plantation = if (length(threat_codes) == 0) {
      FALSE
    } else {
      any(
        str_detect(threat_codes, paste0("^", THREAT_CODE_PLANTATION)) & threat_timing == "Ongoing",
        na.rm = TRUE
      )
    },
    Synonyms = list(extract_synonyms(taxon_info))
  )
}

get_descendant_tips <- function(node, children_list, n_tips) {
  kids <- children_list[[as.character(node)]]
  if (is.null(kids)) return(integer(0))
  
  out <- integer(0)
  for (k in kids) {
    if (k <= n_tips) {
      out <- c(out, k)
    } else {
      out <- c(out, get_descendant_tips(k, children_list, n_tips))
    }
  }
  out
}

get_clade_label_node <- function(tree, tips) {
  tips <- unique(tips[!is.na(tips)])
  if (length(tips) == 0) {
    return(NA_integer_)
  }
  
  if (length(tips) == 1) {
    tip_idx <- match(tips[[1]], tree$tip.label)
    if (is.na(tip_idx)) {
      return(NA_integer_)
    }
    return(as.integer(tip_idx))
  }
  
  mrca_node <- getMRCA(tree, tips)
  if (length(mrca_node) == 0 || is.null(mrca_node) || is.na(mrca_node)) {
    return(NA_integer_)
  }
  
  as.integer(mrca_node)
}

get_family_anchor_tip_from_positions <- function(tips, tip_positions) {
  tips <- unique(tips[!is.na(tips)])
  if (length(tips) == 0) {
    return(NA_integer_)
  }
  
  family_pos <- tip_positions %>%
    filter(label %in% tips) %>%
    distinct(label, tip_node, anchor_angle)
  
  if (nrow(family_pos) == 0) {
    return(NA_integer_)
  }
  
  theta <- family_pos$anchor_angle * pi / 180
  mean_theta <- atan2(mean(sin(theta)), mean(cos(theta))) * 180 / pi
  if (mean_theta < 0) {
    mean_theta <- mean_theta + 360
  }
  
  circular_distance <- abs(((family_pos$anchor_angle - mean_theta + 180) %% 360) - 180)
  family_pos$tip_node[[which.min(circular_distance)]]
}

# Using plotmath via parse = TRUE.
# Safer than underline() for ggplot parse, which can be unreliable.
# We therefore encode:
# - top_count families: bold
# - top_prop families: italic
# - both: bolditalic
family_name_plotmath <- function(family, top_count, top_prop) {
  base_txt <- paste0('"', family, '"')
  
  if (isTRUE(top_count) && isTRUE(top_prop)) {
    return(paste0("bold(underline(", base_txt, "))"))
  }
  if (isTRUE(top_count)) {
    return(paste0("bold(", base_txt, ")"))
  }
  if (isTRUE(top_prop)) {
    return(paste0("underline(", base_txt, ")"))
  }
  base_txt
}

simulate_brownian_threshold_matrix <- function(tree, n_positive, n_sims, progress_label = NULL) {
  n_tips <- length(tree$tip.label)
  out <- matrix(
    0L,
    nrow = n_tips,
    ncol = n_sims,
    dimnames = list(tree$tip.label, paste0("sim_", seq_len(n_sims)))
  )
  
  if (n_positive <= 0L) {
    return(out)
  }
  if (n_positive >= n_tips) {
    out[] <- 1L
    return(out)
  }
  
  progress_start <- if (!is.null(progress_label)) Sys.time() else NULL
  for (i in seq_len(n_sims)) {
    sim_vals <- ape::rTraitCont(tree, model = "BM", sigma = 1)
    if (!is.null(names(sim_vals)) && all(tree$tip.label %in% names(sim_vals))) {
      sim_vals <- sim_vals[tree$tip.label]
    }
    ranked_idx <- order(sim_vals, decreasing = TRUE)
    out[ranked_idx[seq_len(n_positive)], i] <- 1L
    
    if (!is.null(progress_label) && should_emit_progress(i, n_sims)) {
      emit_progress_message(progress_label, i, n_sims, progress_start)
    }
  }
  
  out
}

load_fishtree_taxonomy <- function(path, tree_tips = NULL) {
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  
  taxonomy <- read_csv(path, show_col_types = FALSE)
  required_cols <- c("genus.species", "family", "order")
  missing_cols <- setdiff(required_cols, names(taxonomy))
  if (length(missing_cols) > 0) {
    stop(
      "FishTree taxonomy file is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  taxonomy <- taxonomy %>%
    transmute(
      phylo_name = to_phylo_style(.data[["genus.species"]]),
      Family_fishtree = clean_names(.data[["family"]]),
      Order_fishtree = clean_names(.data[["order"]])
    ) %>%
    filter(!is.na(phylo_name), phylo_name != "") %>%
    group_by(phylo_name) %>%
    summarise(
      Family_fishtree = first_non_na(Family_fishtree),
      Order_fishtree = first_non_na(Order_fishtree),
      .groups = "drop"
    )
  
  if (!is.null(tree_tips)) {
    tree_tips <- to_phylo_style(tree_tips)
    taxonomy <- taxonomy %>% filter(phylo_name %in% tree_tips)
  }
  
  taxonomy
}

apply_fishtree_taxonomy <- function(matched_dat, taxonomy, label = "fish") {
  if (is.null(taxonomy) || nrow(taxonomy) == 0) {
    warning("No FishTree taxonomy rows available for ", label, "; retaining existing family/order labels.")
    return(matched_dat)
  }
  
  matched_out <- matched_dat %>%
    left_join(taxonomy, by = "phylo_name") %>%
    mutate(
      Family = coalesce(Family_fishtree, Family),
      Order = coalesce(Order_fishtree, Order)
    ) %>%
    select(-Family_fishtree, -Order_fishtree)
  
  n_overridden <- matched_out %>%
    left_join(taxonomy, by = "phylo_name") %>%
    filter(!is.na(Family_fishtree), Family == Family_fishtree) %>%
    nrow()
  
  cat(sprintf(
    "   > Applied FishTree taxonomy to %d %s matched species.\n",
    n_overridden,
    label
  ))
  
  matched_out
}

phylo_d_change_sums <- function(tree, state_matrix) {
  if (is.null(dim(state_matrix))) {
    state_matrix <- matrix(state_matrix, ncol = 1)
  }
  if (is.null(colnames(state_matrix))) {
    colnames(state_matrix) <- paste0("V", seq_len(ncol(state_matrix)))
  }
  
  tree_pw <- ape::reorder.phylo(tree, order = "pruningwise")
  n_tips <- length(tree_pw$tip.label)
  
  if (nrow(state_matrix) != n_tips) {
    stop("State matrix row count does not match the number of tree tips.")
  }
  
  if (!is.null(rownames(state_matrix))) {
    if (!all(tree_pw$tip.label %in% rownames(state_matrix))) {
      stop("State matrix row names do not cover every tree tip.")
    }
    state_matrix <- state_matrix[tree_pw$tip.label, , drop = FALSE]
  }
  
  total_nodes <- n_tips + tree_pw$Nnode
  nodal_values <- matrix(
    NA_real_,
    nrow = total_nodes,
    ncol = ncol(state_matrix),
    dimnames = list(as.character(seq_len(total_nodes)), colnames(state_matrix))
  )
  nodal_values[seq_len(n_tips), ] <- state_matrix
  
  children_list <- split(tree_pw$edge[, 2], tree_pw$edge[, 1])
  edge_length_lookup <- setNames(tree_pw$edge.length, tree_pw$edge[, 2])
  parent_nodes <- unique(tree_pw$edge[, 1])
  change_sums <- setNames(numeric(ncol(state_matrix)), colnames(state_matrix))
  
  for (parent in parent_nodes) {
    children <- as.integer(children_list[[as.character(parent)]])
    bl <- as.numeric(edge_length_lookup[as.character(children)])
    vals <- nodal_values[children, , drop = FALSE]
    
    if (length(children) == 1L) {
      curr_nv <- vals[1, ]
      curr_contr <- rep(0, ncol(state_matrix))
    } else if (length(children) == 2L) {
      inv_bl <- 1 / bl
      curr_nv <- colSums(sweep(vals, 1, inv_bl, `*`)) / sum(inv_bl)
      curr_contr <- colSums(abs(sweep(vals, 2, curr_nv, `-`)))
    } else {
      curr_nv <- numeric(ncol(state_matrix))
      curr_contr <- numeric(ncol(state_matrix))
      
      for (j in seq_len(ncol(state_matrix))) {
        rv <- vals[, j]
        group <- if (stats::var(rv) < .Machine$double.eps) {
          c(TRUE, rep(FALSE, length(rv) - 1L))
        } else {
          rv > mean(rv)
        }
        
        inv_bl <- 1 / bl
        prod_val_bl <- c(
          sum(vals[!group, j] * inv_bl[!group]),
          sum(vals[group, j] * inv_bl[group])
        )
        sum_wght <- c(
          sum(inv_bl[!group]),
          sum(inv_bl[group])
        )
        sub_nv <- prod_val_bl / sum_wght
        sub_bl <- 1 / sum_wght
        
        curr_nv[j] <- sum(sub_nv * (1 / sub_bl)) / sum(1 / sub_bl)
        curr_contr[j] <- sum(abs(vals[, j] - curr_nv[j]))
      }
    }
    
    nodal_values[parent, ] <- curr_nv
    change_sums <- change_sums + curr_contr
  }
  
  change_sums
}

interpret_phylo_d <- function(d_estimate) {
  if (is.na(d_estimate)) {
    return(NA_character_)
  }
  if (d_estimate < 0) {
    return("More conserved than Brownian")
  }
  if (d_estimate < 1) {
    return("Between Brownian and random")
  }
  if (d_estimate > 1) {
    return("Overdispersed relative to random")
  }
  "Random expectation"
}

build_phylo_distance_context <- function(tree) {
  tree_pw <- reorder.phylo(tree, order = "postorder")
  n_tips <- length(tree_pw$tip.label)
  total_nodes <- n_tips + tree_pw$Nnode
  
  if (is.null(tree_pw$edge.length) || anyNA(tree_pw$edge.length)) {
    stop("SES.MPD and SES.MNTD require branch lengths for every edge.")
  }
  
  root <- setdiff(tree_pw$edge[, 1], tree_pw$edge[, 2])
  if (length(root) != 1L) {
    stop("Expected a single rooted tree when building the phylogenetic distance context.")
  }
  root <- as.integer(root)
  
  edge_length_to_parent <- numeric(total_nodes)
  edge_length_to_parent[tree_pw$edge[, 2]] <- tree_pw$edge.length
  
  children_by_node <- vector("list", total_nodes)
  split_children <- split(tree_pw$edge[, 2], tree_pw$edge[, 1])
  for (parent in names(split_children)) {
    children_by_node[[as.integer(parent)]] <- as.integer(split_children[[parent]])
  }
  
  preorder_nodes <- integer(0)
  postorder_nodes <- integer(0)
  stack_nodes <- root
  stack_expanded <- FALSE
  
  while (length(stack_nodes) > 0L) {
    current_idx <- length(stack_nodes)
    node <- stack_nodes[[current_idx]]
    expanded <- stack_expanded[[current_idx]]
    stack_nodes <- stack_nodes[-current_idx]
    stack_expanded <- stack_expanded[-current_idx]
    
    if (expanded) {
      postorder_nodes <- c(postorder_nodes, node)
      next
    }
    
    preorder_nodes <- c(preorder_nodes, node)
    stack_nodes <- c(stack_nodes, node)
    stack_expanded <- c(stack_expanded, TRUE)
    
    children <- children_by_node[[node]]
    if (length(children) > 0L) {
      stack_nodes <- c(stack_nodes, rev(children))
      stack_expanded <- c(stack_expanded, rep(FALSE, length(children)))
    }
  }
  
  list(
    n_tips = n_tips,
    total_nodes = total_nodes,
    root = root,
    children_by_node = children_by_node,
    edge_length_to_parent = edge_length_to_parent,
    non_root_nodes = setdiff(seq_len(total_nodes), root),
    preorder_internal = preorder_nodes[preorder_nodes > n_tips],
    postorder_internal = postorder_nodes[postorder_nodes > n_tips]
  )
}

calculate_pairwise_phylo_metrics <- function(phylo_context, selected_tip_indices) {
  selected_tip_indices <- sort(unique(as.integer(selected_tip_indices)))
  n_selected <- length(selected_tip_indices)
  
  if (n_selected < 2L) {
    return(list(mpd = NA_real_, mntd = NA_real_))
  }
  
  counts <- integer(phylo_context$total_nodes)
  counts[selected_tip_indices] <- 1L
  
  down_ids <- matrix(0L, nrow = phylo_context$total_nodes, ncol = 2L)
  down_distances <- matrix(Inf, nrow = phylo_context$total_nodes, ncol = 2L)
  down_ids[selected_tip_indices, 1] <- selected_tip_indices
  down_distances[selected_tip_indices, 1] <- 0
  
  for (node in phylo_context$postorder_internal) {
    children <- phylo_context$children_by_node[[node]]
    node_count <- 0L
    candidate_ids <- integer(0)
    candidate_distances <- numeric(0)
    
    for (child in children) {
      node_count <- node_count + counts[child]
      child_edge_length <- phylo_context$edge_length_to_parent[[child]]
      
      if (down_ids[child, 1] > 0L) {
        candidate_ids <- c(candidate_ids, down_ids[child, 1])
        candidate_distances <- c(candidate_distances, down_distances[child, 1] + child_edge_length)
      }
      if (down_ids[child, 2] > 0L) {
        candidate_ids <- c(candidate_ids, down_ids[child, 2])
        candidate_distances <- c(candidate_distances, down_distances[child, 2] + child_edge_length)
      }
    }
    
    counts[node] <- node_count
    
    if (length(candidate_ids) > 0L) {
      ord <- order(candidate_distances, candidate_ids)
      candidate_ids <- candidate_ids[ord]
      candidate_distances <- candidate_distances[ord]
      keep <- !duplicated(candidate_ids)
      candidate_ids <- candidate_ids[keep]
      candidate_distances <- candidate_distances[keep]
      n_keep <- min(2L, length(candidate_ids))
      down_ids[node, seq_len(n_keep)] <- candidate_ids[seq_len(n_keep)]
      down_distances[node, seq_len(n_keep)] <- candidate_distances[seq_len(n_keep)]
    }
  }
  
  pairwise_distance_sum <- sum(
    phylo_context$edge_length_to_parent[phylo_context$non_root_nodes] *
      counts[phylo_context$non_root_nodes] *
      (n_selected - counts[phylo_context$non_root_nodes])
  )
  mpd_observed <- pairwise_distance_sum / choose(n_selected, 2)
  
  up_ids <- matrix(0L, nrow = phylo_context$total_nodes, ncol = 2L)
  up_distances <- matrix(Inf, nrow = phylo_context$total_nodes, ncol = 2L)
  
  for (parent in phylo_context$preorder_internal) {
    children <- phylo_context$children_by_node[[parent]]
    if (length(children) == 0L) {
      next
    }
    
    for (child in children) {
      child_edge_length <- phylo_context$edge_length_to_parent[[child]]
      candidate_ids <- integer(0)
      candidate_distances <- numeric(0)
      
      if (up_ids[parent, 1] > 0L) {
        candidate_ids <- c(candidate_ids, up_ids[parent, 1])
        candidate_distances <- c(candidate_distances, up_distances[parent, 1] + child_edge_length)
      }
      if (up_ids[parent, 2] > 0L) {
        candidate_ids <- c(candidate_ids, up_ids[parent, 2])
        candidate_distances <- c(candidate_distances, up_distances[parent, 2] + child_edge_length)
      }
      
      sibling_nodes <- children[children != child]
      if (length(sibling_nodes) > 0L) {
        for (sibling in sibling_nodes) {
          sibling_edge_length <- phylo_context$edge_length_to_parent[[sibling]]
          
          if (down_ids[sibling, 1] > 0L) {
            candidate_ids <- c(candidate_ids, down_ids[sibling, 1])
            candidate_distances <- c(
              candidate_distances,
              down_distances[sibling, 1] + sibling_edge_length + child_edge_length
            )
          }
          if (down_ids[sibling, 2] > 0L) {
            candidate_ids <- c(candidate_ids, down_ids[sibling, 2])
            candidate_distances <- c(
              candidate_distances,
              down_distances[sibling, 2] + sibling_edge_length + child_edge_length
            )
          }
        }
      }
      
      if (length(candidate_ids) > 0L) {
        ord <- order(candidate_distances, candidate_ids)
        candidate_ids <- candidate_ids[ord]
        candidate_distances <- candidate_distances[ord]
        keep <- !duplicated(candidate_ids)
        candidate_ids <- candidate_ids[keep]
        candidate_distances <- candidate_distances[keep]
        n_keep <- min(2L, length(candidate_ids))
        up_ids[child, seq_len(n_keep)] <- candidate_ids[seq_len(n_keep)]
        up_distances[child, seq_len(n_keep)] <- candidate_distances[seq_len(n_keep)]
      }
    }
  }
  
  nearest_distances <- up_distances[selected_tip_indices, 1]
  nearest_distances[!is.finite(nearest_distances)] <- NA_real_
  mntd_observed <- if (all(is.na(nearest_distances))) {
    NA_real_
  } else {
    mean(nearest_distances, na.rm = TRUE)
  }
  
  list(
    mpd = mpd_observed,
    mntd = mntd_observed
  )
}

summarise_null_distribution <- function(observed_value, random_values) {
  random_mean <- mean(random_values, na.rm = TRUE)
  random_sd <- sd(random_values, na.rm = TRUE)
  z_score <- if (isTRUE(all.equal(random_sd, 0, tolerance = 1e-12))) {
    NA_real_
  } else {
    (observed_value - random_mean) / random_sd
  }
  
  observed_rank <- if (is.na(observed_value) || anyNA(random_values)) {
    NA_real_
  } else {
    rank(c(observed_value, random_values))[1]
  }
  
  list(
    observed = observed_value,
    random_mean = random_mean,
    random_sd = random_sd,
    z_score = z_score,
    p_value = if (is.na(observed_rank)) NA_real_ else observed_rank / (length(random_values) + 1)
  )
}

calculate_ses_pairwise_metrics <- function(
  phylo_context,
  states,
  permutations,
  seed = NULL,
  progress_label = NULL
) {
  states <- as.integer(states)
  
  if (length(states) != phylo_context$n_tips) {
    stop("SES.MPD and SES.MNTD state vector length does not match the number of tree tips.")
  }
  if (anyNA(states)) {
    stop("SES.MPD and SES.MNTD input contains missing states.")
  }
  if (!all(states %in% c(0L, 1L))) {
    stop("SES.MPD and SES.MNTD input must be binary (0/1).")
  }
  
  positive_tip_indices <- which(states == 1L)
  n_positive <- length(positive_tip_indices)
  if (n_positive < 2L) {
    return(list(
      n_positive = n_positive,
      ses_mpd = NA_real_,
      ses_mpd_obs = NA_real_,
      ses_mpd_rand_mean = NA_real_,
      ses_mpd_rand_sd = NA_real_,
      ses_mpd_p = NA_real_,
      ses_mntd = NA_real_,
      ses_mntd_obs = NA_real_,
      ses_mntd_rand_mean = NA_real_,
      ses_mntd_rand_sd = NA_real_,
      ses_mntd_p = NA_real_
    ))
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  observed_metrics <- calculate_pairwise_phylo_metrics(phylo_context, positive_tip_indices)
  mpd_random <- numeric(permutations)
  mntd_random <- numeric(permutations)
  progress_start <- if (!is.null(progress_label)) Sys.time() else NULL
  
  for (perm_idx in seq_len(permutations)) {
    random_tip_indices <- sample.int(phylo_context$n_tips, size = n_positive, replace = FALSE)
    random_metrics <- calculate_pairwise_phylo_metrics(phylo_context, random_tip_indices)
    mpd_random[[perm_idx]] <- random_metrics$mpd
    mntd_random[[perm_idx]] <- random_metrics$mntd
    
    if (!is.null(progress_label) && should_emit_progress(perm_idx, permutations)) {
      emit_progress_message(progress_label, perm_idx, permutations, progress_start)
    }
  }
  
  mpd_summary <- summarise_null_distribution(observed_metrics$mpd, mpd_random)
  mntd_summary <- summarise_null_distribution(observed_metrics$mntd, mntd_random)
  
  list(
    n_positive = n_positive,
    ses_mpd = mpd_summary$z_score,
    ses_mpd_obs = mpd_summary$observed,
    ses_mpd_rand_mean = mpd_summary$random_mean,
    ses_mpd_rand_sd = mpd_summary$random_sd,
    ses_mpd_p = mpd_summary$p_value,
    ses_mntd = mntd_summary$z_score,
    ses_mntd_obs = mntd_summary$observed,
    ses_mntd_rand_mean = mntd_summary$random_mean,
    ses_mntd_rand_sd = mntd_summary$random_sd,
    ses_mntd_p = mntd_summary$p_value
  )
}

calculate_phylo_d_statistic <- function(
  tree,
  states,
  permutations,
  seed = NULL,
  progress_label = NULL
) {
  states <- as.integer(states)
  
  if (length(states) != length(tree$tip.label)) {
    stop("State vector length does not match the number of tree tips.")
  }
  if (anyNA(states)) {
    stop("D-statistic input contains missing states.")
  }
  if (!all(states %in% c(0L, 1L))) {
    stop("D-statistic input must be binary (0/1).")
  }
  if (length(unique(states)) < 2L) {
    stop("D-statistic input must contain both binary states.")
  }
  if (any(tree$edge.length <= 0, na.rm = TRUE)) {
    stop("D-statistic requires strictly positive branch lengths.")
  }
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_positive <- sum(states == 1L)
  if (!is.null(progress_label)) {
    cat(sprintf("      %s: building random null permutations (%d)\n", progress_label, permutations))
  }
  random_matrix <- replicate(
    permutations,
    sample(states, replace = FALSE),
    simplify = "matrix"
  )
  colnames(random_matrix) <- paste0("random_", seq_len(permutations))
  
  if (!is.null(progress_label)) {
    cat(sprintf("      %s: simulating Brownian null permutations (%d)\n", progress_label, permutations))
  }
  brownian_matrix <- simulate_brownian_threshold_matrix(
    tree = tree,
    n_positive = n_positive,
    n_sims = permutations,
    progress_label = if (is.null(progress_label)) NULL else paste0(progress_label, " Brownian")
  )
  colnames(brownian_matrix) <- paste0("brownian_", seq_len(permutations))
  
  if (!is.null(progress_label)) {
    cat(sprintf("      %s: summarising observed vs null change sums\n", progress_label))
  }
  random_sums <- phylo_d_change_sums(
    tree,
    cbind(Observed = states, random_matrix)
  )
  brownian_sums <- phylo_d_change_sums(
    tree,
    cbind(Observed = states, brownian_matrix)
  )
  
  observed_sum <- unname(random_sums[["Observed"]])
  brownian_observed_sum <- unname(brownian_sums[["Observed"]])
  if (!isTRUE(all.equal(observed_sum, brownian_observed_sum, tolerance = 1e-8))) {
    stop("Observed D-statistic change sum did not match between random and Brownian runs.")
  }
  
  random_sums <- unname(random_sums[names(random_sums) != "Observed"])
  brownian_sums <- unname(brownian_sums[names(brownian_sums) != "Observed"])
  
  mean_random <- mean(random_sums)
  mean_brownian <- mean(brownian_sums)
  denominator <- mean_random - mean_brownian
  d_estimate <- if (isTRUE(all.equal(denominator, 0, tolerance = 1e-12))) {
    NA_real_
  } else {
    (observed_sum - mean_brownian) / denominator
  }
  
  list(
    d_estimate = d_estimate,
    p_random = mean(random_sums < observed_sum),
    p_brownian = mean(brownian_sums > observed_sum),
    observed_sum = observed_sum,
    mean_random_sum = mean_random,
    mean_brownian_sum = mean_brownian,
    n_positive = n_positive
  )
}

summarise_taxon_phylo_d <- function(
  tree,
  tip_data,
  taxon,
  permutations,
  seed_base,
  phylo_context = build_phylo_distance_context(tree)
) {
  trait_specs <- tibble(
    trait = c("Wood_Threat", "Threat_Logging", "Threat_Plantation"),
    trait_label = c("Wood production", "Logging", "Plantation")
  )
  
  purrr::pmap_dfr(
    trait_specs,
    function(trait, trait_label) {
      states <- as.integer(replace_na(tip_data[[trait]], FALSE))
      n_positive <- sum(states == 1L)
      n_tips <- length(states)
      trait_seed_offset <- match(trait, trait_specs$trait)
      cat(sprintf("   > %s: starting phylogenetic signal calculations\n", trait_label))
      
      if (length(unique(states)) < 2L) {
        return(tibble(
          taxon = taxon,
          trait = trait,
          trait_label = trait_label,
          n_tips = n_tips,
          n_threatened = n_positive,
          threat_prevalence = n_positive / n_tips,
          d_estimate = NA_real_,
          p_random = NA_real_,
          p_brownian = NA_real_,
          observed_sum = NA_real_,
          mean_random_sum = NA_real_,
          mean_brownian_sum = NA_real_,
          ses_mpd = NA_real_,
          ses_mpd_obs = NA_real_,
          ses_mpd_rand_mean = NA_real_,
          ses_mpd_rand_sd = NA_real_,
          ses_mpd_p = NA_real_,
          ses_mntd = NA_real_,
          ses_mntd_obs = NA_real_,
          ses_mntd_rand_mean = NA_real_,
          ses_mntd_rand_sd = NA_real_,
          ses_mntd_p = NA_real_,
          permutations = permutations,
          signal_class = NA_character_,
          clustered_vs_random = NA,
          note = "Only one observed state on the matched phylogeny tips."
        ))
      }
      
      result <- calculate_phylo_d_statistic(
        tree = tree,
        states = states,
        permutations = permutations,
        seed = seed_base + trait_seed_offset,
        progress_label = paste0(trait_label, " D")
      )
      cat(sprintf("      %s: D statistic complete\n", trait_label))
      ses_result <- calculate_ses_pairwise_metrics(
        phylo_context = phylo_context,
        states = states,
        permutations = permutations,
        seed = seed_base + 100L + trait_seed_offset,
        progress_label = paste0(trait_label, " SES")
      )
      cat(sprintf("      %s: SES.MPD and SES.MNTD complete\n", trait_label))
      
      tibble(
        taxon = taxon,
        trait = trait,
        trait_label = trait_label,
        n_tips = n_tips,
        n_threatened = result$n_positive,
        threat_prevalence = result$n_positive / n_tips,
        d_estimate = result$d_estimate,
        p_random = result$p_random,
        p_brownian = result$p_brownian,
        observed_sum = result$observed_sum,
        mean_random_sum = result$mean_random_sum,
        mean_brownian_sum = result$mean_brownian_sum,
        ses_mpd = ses_result$ses_mpd,
        ses_mpd_obs = ses_result$ses_mpd_obs,
        ses_mpd_rand_mean = ses_result$ses_mpd_rand_mean,
        ses_mpd_rand_sd = ses_result$ses_mpd_rand_sd,
        ses_mpd_p = ses_result$ses_mpd_p,
        ses_mntd = ses_result$ses_mntd,
        ses_mntd_obs = ses_result$ses_mntd_obs,
        ses_mntd_rand_mean = ses_result$ses_mntd_rand_mean,
        ses_mntd_rand_sd = ses_result$ses_mntd_rand_sd,
        ses_mntd_p = ses_result$ses_mntd_p,
        permutations = permutations,
        signal_class = interpret_phylo_d(result$d_estimate),
        clustered_vs_random = !is.na(result$d_estimate) &&
          result$d_estimate < 1 &&
          result$p_random <= 0.05,
        note = if (ses_result$n_positive < 2L) {
          "Fewer than two threatened species on the matched phylogeny tips; SES.MPD and SES.MNTD were not computed."
        } else {
          NA_character_
        }
      )
    }
  )
}

get_tree_depth <- function(tree) {
  max(ape::node.depth.edgelength(tree), na.rm = TRUE)
}

scale_tree_to_depth <- function(tree, target_depth = NULL) {
  if (is.null(target_depth) || !is.finite(target_depth) || target_depth <= 0) {
    return(tree)
  }
  if (is.null(tree$edge.length)) {
    return(tree)
  }
  
  current_depth <- get_tree_depth(tree)
  if (!is.finite(current_depth) || current_depth <= 0) {
    return(tree)
  }
  
  scaled_tree <- tree
  scaled_tree$edge.length <- scaled_tree$edge.length * (target_depth / current_depth)
  scaled_tree
}

build_phylogram_plot_set <- function(
  tree,
  tip_data,
  cfg,
  panel_title = NULL,
  target_tree_depth = NULL,
  plot_margin = margin(20, 180, 20, 120),
  print_family_tables = TRUE
) {
  source_tree_depth <- get_tree_depth(tree)
  plot_tree <- scale_tree_to_depth(tree, target_tree_depth)
  plot_tree_depth <- get_tree_depth(plot_tree)
  ring_scale_factor <- if (is.finite(source_tree_depth) && source_tree_depth > 0) {
    plot_tree_depth / source_tree_depth
  } else {
    1
  }
  plot_ring_interval <- phylo_ring_interval_myr * ring_scale_factor
  n_tips <- length(plot_tree$tip.label)
  children_list <- split(plot_tree$edge[, 2], plot_tree$edge[, 1])
  all_nodes <- seq_len(n_tips + plot_tree$Nnode)
  phylo_ring_band_starts <- seq(0, plot_tree_depth, by = plot_ring_interval * 2)
  phylo_ring_band_df <- tibble(
    xmin = phylo_ring_band_starts,
    xmax = pmin(phylo_ring_band_starts + plot_ring_interval, plot_tree_depth)
  ) %>%
    filter(xmax > xmin)
  
  node_summary <- lapply(all_nodes, function(node_id) {
    desc_tips <- if (node_id <= n_tips) node_id else get_descendant_tips(node_id, children_list, n_tips)
    tibble(
      node = node_id,
      prop_wood_threat = mean(tip_data$Wood_Threat[desc_tips], na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  p <- ggtree(plot_tree, layout = "circular", size = branch_line_size) %<+% node_summary
  if (nrow(phylo_ring_band_df) > 0) {
    p <- p +
      geom_rect(
        data = phylo_ring_band_df,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = phylo_ring_fill_col,
        colour = NA,
        show.legend = FALSE
      )
  }
  p <- p +
    geom_tree(aes(colour = prop_wood_threat), size = branch_line_size, lineend = "butt") +
    scale_colour_gradient2(
      low = branch_low_col,
      mid = branch_mid_col,
      high = branch_high_col,
      midpoint = 0.5,
      limits = c(0, 1),
      name = "Proportion of descendants threatened"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5),
      plot.margin = plot_margin
    )
  
  ring_df <- tip_data %>%
    select(label, tip_category) %>%
    as.data.frame()
  rownames(ring_df) <- ring_df$label
  ring_df$label <- NULL
  
  p2 <- gheatmap(
    p,
    ring_df,
    offset = ring_offset,
    width = ring_width,
    colnames = FALSE,
    color = NA
  ) +
    scale_fill_manual(
      values = c(
        "None" = col_none,
        "Logging only" = col_logging_only,
        "Plantation only" = col_plantation_only,
        "Both" = col_both
      ),
      limits = c("None", "Plantation only", "Both", "Logging only"),
      labels = c("None", "Plantation", "Both", "Natural Forestry"),
      name = "Production Pathway"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5)
    )
  family_ring_offset <- ring_offset + (plot_tree_depth * ring_width) + family_ring_gap
  tip_positions <- p$data %>%
    filter(isTip) %>%
    transmute(
      label,
      tip_node = node,
      anchor_y = y,
      anchor_angle = angle
    )
  
  family_summary <- tip_data %>%
    filter(!is.na(Family), Family != "") %>%
    group_by(Family) %>%
    summarise(
      total_species = n(),
      n_wood_threat = sum(Wood_Threat, na.rm = TRUE),
      prop_wood_threat = n_wood_threat / total_species,
      family_tips = list(label),
      .groups = "drop"
    ) %>%
    mutate(
      is_monophyletic = purrr::map_lgl(
        family_tips,
        ~ if (length(.x) > 1) ape::is.monophyletic(plot_tree, .x) else TRUE
      )
    )

  family_label_candidates <- if (isTRUE(cfg$label_monophyletic_families_only)) {
    family_summary %>% filter(is_monophyletic)
  } else {
    family_summary
  }
  
  top_family_count <- family_label_candidates %>%
    arrange(desc(n_wood_threat), Family) %>%
    slice_head(n = 10) %>%
    transmute(Family, top_count = TRUE)
  
  top_family_prop <- family_label_candidates %>%
    filter(n_wood_threat > 0) %>%
    arrange(desc(prop_wood_threat), desc(n_wood_threat), Family) %>%
    slice_head(n = 10) %>%
    transmute(Family, top_prop = TRUE)
  
  family_label_table <- family_summary %>%
    left_join(top_family_count, by = "Family") %>%
    left_join(top_family_prop, by = "Family") %>%
    mutate(
      top_count = replace_na(top_count, FALSE),
      top_prop = replace_na(top_prop, FALSE)
    ) %>%
    filter(top_count | top_prop) %>%
    mutate(
      clade_node = purrr::map_int(
        family_tips,
        ~ get_clade_label_node(plot_tree, .x)
      ),
      anchor_tip_node = purrr::map_int(
        family_tips,
        ~ get_family_anchor_tip_from_positions(.x, tip_positions)
      ),
      label_stats = paste0(
        '"',
        format(n_wood_threat, big.mark = ","),
        ' spp.; ',
        scales::percent(prop_wood_threat, accuracy = 0.1),
        '"'
      ),
      family_name_expr = purrr::pmap_chr(
        list(Family, top_count, top_prop),
        family_name_plotmath
      ),
      label_full_expr = paste0("atop(", family_name_expr, ", ", label_stats, ")")
    ) %>%
    left_join(tip_positions, by = c("anchor_tip_node" = "tip_node")) %>%
    mutate(
      raw_angle = dplyr::coalesce(anchor_angle, 0),
      label_hjust = if_else(raw_angle > 90 & raw_angle < 270, 1, 0)
    ) %>%
    arrange(
      anchor_tip_node,
      desc(top_count),
      desc(n_wood_threat),
      desc(top_prop),
      desc(prop_wood_threat),
      Family
    )
  
  labelled_families <- family_label_table %>%
    pull(Family) %>%
    unique()
  
  family_ring_df <- tip_data %>%
    mutate(
      family_highlight = if_else(
        Family %in% labelled_families,
        "Labelled family",
        "Other family"
      )
    ) %>%
    select(label, family_highlight) %>%
    as.data.frame()
  rownames(family_ring_df) <- family_ring_df$label
  family_ring_df$label <- NULL
  family_ring_df$family_highlight <- factor(
    family_ring_df$family_highlight,
    levels = c("Other family", "Labelled family")
  )
  
  p3 <- p2 +
    ggnewscale::new_scale_fill()
  
  p3 <- gheatmap(
    p3,
    family_ring_df,
    offset = family_ring_offset,
    width = family_ring_width,
    colnames = FALSE,
    color = NA
  ) +
    scale_fill_manual(
      values = c(
        "Other family" = col_family_other,
        "Labelled family" = col_family_highlight
      ),
      name = "Labelled family\nmembership"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10.5)
    )
  p_no_family_labels <- p3
  
  p3_build <- ggplot_build(p3)
  ring_layer_idx <- which(vapply(
    p3_build$data,
    function(layer_df) {
      all(c("xmin", "xmax", "ymin", "ymax") %in% names(layer_df)) &&
        nrow(layer_df) == nrow(family_ring_df)
    },
    logical(1)
  ))
  ring_outer_x <- if (length(ring_layer_idx) > 0) {
    max(p3_build$data[[tail(ring_layer_idx, 1)]]$xmax, na.rm = TRUE)
  } else {
    plot_tree_depth * (1.02 + family_ring_offset + family_ring_width)
  }
  
  if (nrow(family_label_table) == 0) {
    warning("No threatened families found for ", cfg$taxon)
  } else {
    if (isTRUE(print_family_tables)) {
      print(
        family_label_table %>%
          select(
            Family, total_species, n_wood_threat, prop_wood_threat,
            top_count, top_prop, is_monophyletic, clade_node, anchor_tip_node
          )
      )
    }
    
    family_label_table <- family_label_table %>%
      filter(
        !is.na(anchor_tip_node),
        !is.na(anchor_y),
        !is.na(label_hjust)
      ) %>%
      filter(
        !label_families_only_if_monophyletic | is_monophyletic
      ) %>%
      mutate(label_side = if_else(label_hjust == 1, "left", "right")) %>%
      group_by(label_side) %>%
      arrange(anchor_y, .by_group = TRUE) %>%
      mutate(
        stagger_band = (row_number() - 1L) %% max(1L, as.integer(family_label_stagger_bands)),
        anchor_x = ring_outer_x,
        label_nudge_x = family_label_gap + (stagger_band * family_label_stagger_step)
      ) %>%
      ungroup()
    
    p3 <- p3 +
      geom_text_repel(
        data        = family_label_table,
        aes(
          x     = anchor_x,
          y     = anchor_y,
          label = label_full_expr,
          hjust = label_hjust
        ),
        size        = family_label_fontsize,
        color       = "black",
        parse       = TRUE,
        inherit.aes = FALSE,
        nudge_x     = family_label_table$label_nudge_x,
        direction   = "y",
        force       = 1.5,
        force_pull  = 0.5,
        segment.color = "grey40",
        segment.size = 0.3,
        segment.linetype = "dotted",
        min.segment.length = 0,
        box.padding = 0.35,
        point.padding = 0.1,
        max.overlaps = Inf,
        seed = 42
      )
  }
  
  if (!is.null(panel_title)) {
    title_theme <- theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
    p_no_family_labels <- p_no_family_labels + ggtitle(panel_title) + title_theme
    p3 <- p3 + ggtitle(panel_title) + title_theme
  }
  
  list(
    no_family_labels_plot = p_no_family_labels,
    family_labels_plot = p3,
    plot_tree_depth = plot_tree_depth
  )
}

build_taxon_phylogram <- function(cfg) {
  cat("\n=========================================================\n")
  cat("Taxon: ", toupper(cfg$taxon), "\n", sep = "")
  cat("=========================================================\n")
  
  ensure_parent_dir(cfg$output_file)
  dir.create(cfg$processed_dir, recursive = TRUE, showWarnings = FALSE)
  
  full_tree_file <- if (cfg$taxon == "fish") {
    file.path(cfg$processed_dir, "final_fish_tree_pruned.tre")
  } else {
    file.path(cfg$processed_dir, paste0(cfg$taxon, "_full_tree.tre"))
  }
  matched_data_file <- file.path(cfg$processed_dir, paste0(cfg$taxon, "_iucn_tree_matched.csv"))
  
  # Step 1: load the cached full tree when available.
  cat("--- Step 1: Loading/building full tree ---\n")
  if (!force_rebuild_full_tree && file.exists(full_tree_file)) {
    tree <- read.tree(full_tree_file)
    cat("Loading full tree from cache...\n")
  } else {
    if (identical(tolower(cfg$tree_format), "nexus")) {
      cat("Computing single MCC tree from raw nexus posterior and caching as .tre...\n")
    } else {
      cat("Building full .tre cache from raw tree file...\n")
    }
    tree <- read_tree_file(cfg$raw_tree_file, cfg$tree_format)
    write.tree(tree, full_tree_file)
    cat("Saved full tree to: ", full_tree_file, "\n", sep = "")
  }
  cat(sprintf("   > Tree tips: %d\n", length(tree$tip.label)))
  
  selected_tree_file <- full_tree_file
  tree_tips_cleaned <- to_phylo_style(tree$tip.label)
  valid_tips_df <- tibble(Tree_Tip = tree_tips_cleaned)
  
  # Step 2: match IUCN assessments to tree tips
  cat("--- Step 2: Loading/building matched IUCN table ---\n")
  processed_iucn_lookup_file <- here::here(
    "data", "processed", "assessments", paste0(cfg$taxon, "_assessment_taxonomy_lookup.rds")
  )
  
  matched_dat <- NULL
  if (!force_rebuild_matched_data && file.exists(matched_data_file)) {
    if (file.exists(processed_iucn_lookup_file)) {
      iucn_lookup <- readRDS(processed_iucn_lookup_file) %>%
        transmute(
          IUCN_Name = clean_names(IUCN_Name),
          Family_lookup = clean_names(Family),
          Order_lookup = clean_names(Order)
        )
    } else if (file.exists(cfg$raw_iucn_cache_file)) {
      iucn_raw <- readRDS(cfg$raw_iucn_cache_file)
      iucn_lookup <- build_iucn_taxonomy_lookup(iucn_raw)
    } else {
      class_assessments <- rl_class(class = cfg$iucn_class, latest = TRUE, quiet = TRUE, scope_code = 1)
      iucn_raw <- rl_assessment_list(class_assessments$assessments$assessment_id)
      saveRDS(iucn_raw, cfg$raw_iucn_cache_file)
      iucn_lookup <- build_iucn_taxonomy_lookup(iucn_raw)
    }
    
    cached_candidate <- standardize_match_table(
      read_csv(matched_data_file, show_col_types = FALSE),
      iucn_lookup = iucn_lookup
    )
    
    if (!is.null(cached_candidate) && nrow(cached_candidate) > 0) {
      cached_overlap <- sum(unique(cached_candidate$phylo_name) %in% valid_tips_df$Tree_Tip)
      
      if (cached_overlap == nrow(cached_candidate)) {
        matched_dat <- cached_candidate
        cat("Loaded matched table from: ", matched_data_file, "\n", sep = "")
        cat(sprintf("   > Matched species available for plotting: %d\n", nrow(matched_dat)))
      } else {
        cat(
          "Cached matched table is out of sync with the current tree (",
          cached_overlap, " of ", nrow(cached_candidate),
          " rows overlap); rebuilding from IUCN assessments.\n",
          sep = ""
        )
      }
    }
  }
  
  if (is.null(matched_dat)) {
    if (file.exists(processed_iucn_lookup_file)) {
      iucn_clean <- readRDS(processed_iucn_lookup_file) %>%
        mutate(
          IUCN_Name = clean_names(IUCN_Name),
          Family = clean_names(Family),
          Order = clean_names(Order),
          Synonyms = purrr::map(Synonyms, clean_names)
        ) %>%
        filter(!RedList_Cat %in% c("EX", "EW"))
    } else if (file.exists(cfg$raw_iucn_cache_file)) {
      iucn_raw <- readRDS(cfg$raw_iucn_cache_file)
    } else {
      class_assessments <- rl_class(class = cfg$iucn_class, latest = TRUE, quiet = TRUE, scope_code = 1)
      iucn_raw <- rl_assessment_list(class_assessments$assessments$assessment_id)
      saveRDS(iucn_raw, cfg$raw_iucn_cache_file)
    }
    
    if (!exists("iucn_clean", inherits = FALSE)) {
      raw_rows <- vector("list", length(iucn_raw))
      for (i in seq_along(iucn_raw)) {
        raw_rows[[i]] <- tryCatch(
          extract_assessment_row(iucn_raw[[i]], taxon = cfg$taxon),
          error = function(e) {
            assessment_id <- if (!is.null(iucn_raw[[i]]$assessment_id)) iucn_raw[[i]]$assessment_id else i
            warning("Skipping malformed assessment ", assessment_id, " for ", cfg$taxon, ": ", conditionMessage(e))
            NULL
          }
        )
      }
      iucn_raw_df <- dplyr::bind_rows(purrr::compact(raw_rows))
      if (nrow(iucn_raw_df) == 0) {
        stop("No extant IUCN assessments were extracted for ", cfg$taxon)
      }
      
      iucn_clean <- iucn_raw_df %>%
        filter(!category %in% c("EX", "EW")) %>%
        group_by(assessment_id, scientific_name, Family, Order, category) %>%
        summarise(
          Threat_Logging = any(Threat_Logging, na.rm = TRUE),
          Threat_Plantation = any(Threat_Plantation, na.rm = TRUE),
          Synonyms = list(unique(unlist(Synonyms))),
          .groups = "drop"
        ) %>%
        rename(
          IUCN_Name = scientific_name,
          RedList_Cat = category
        ) %>%
        mutate(
          IUCN_Name = clean_names(IUCN_Name),
          Family = clean_names(Family),
          Order = clean_names(Order)
        ) %>%
        select(-assessment_id)
    }
    
    synonyms_long <- iucn_clean %>%
      select(IUCN_Name, Synonyms) %>%
      unnest(Synonyms) %>%
      transmute(
        IUCN_Name,
        Search_Name = to_phylo_style(Synonyms)
      ) %>%
      filter(!is.na(Search_Name), Search_Name != "")
    
    self_ref <- iucn_clean %>%
      transmute(
        IUCN_Name,
        Search_Name = to_phylo_style(IUCN_Name)
      )
    
    master_name_dictionary <- bind_rows(synonyms_long, self_ref) %>%
      distinct() %>%
      filter(!is.na(Search_Name), Search_Name != "")
    
    iucn_to_tree_link <- master_name_dictionary %>%
      inner_join(valid_tips_df, by = c("Search_Name" = "Tree_Tip")) %>%
      group_by(IUCN_Name) %>%
      slice(1) %>%
      ungroup() %>%
      transmute(IUCN_Name, Valid_Tree_Tip = Search_Name)
    
    iucn_names_clean <- iucn_clean %>%
      transmute(IUCN_Name, IUCN_Name_Clean = to_phylo_style(IUCN_Name))
    
    fuzzy_candidates <- iucn_names_clean %>%
      filter(!IUCN_Name_Clean %in% valid_tips_df$Tree_Tip) %>%
      anti_join(iucn_to_tree_link, by = "IUCN_Name") %>%
      pull(IUCN_Name_Clean) %>%
      unique()
    
    matches_index <- amatch(fuzzy_candidates, valid_tips_df$Tree_Tip, maxDist = 2, nomatch = NA)
    fuzzy_lookup <- tibble(
      IUCN_Name_Clean = fuzzy_candidates,
      Fuzzy_Tree_Tip = valid_tips_df$Tree_Tip[matches_index]
    ) %>%
      filter(!is.na(Fuzzy_Tree_Tip))
    
    matched_dat <- iucn_clean %>%
      left_join(iucn_names_clean, by = "IUCN_Name") %>%
      mutate(Direct_Match_Exists = IUCN_Name_Clean %in% valid_tips_df$Tree_Tip) %>%
      left_join(iucn_to_tree_link, by = "IUCN_Name") %>%
      left_join(fuzzy_lookup, by = "IUCN_Name_Clean") %>%
      mutate(
        phylo_name = case_when(
          Direct_Match_Exists ~ IUCN_Name_Clean,
          !is.na(Valid_Tree_Tip) ~ Valid_Tree_Tip,
          !is.na(Fuzzy_Tree_Tip) ~ Fuzzy_Tree_Tip,
          TRUE ~ NA_character_
        ),
        Match_Type = case_when(
          Direct_Match_Exists ~ "Direct",
          !is.na(Valid_Tree_Tip) ~ "Synonym",
          !is.na(Fuzzy_Tree_Tip) ~ "Fuzzy",
          TRUE ~ "None"
        )
      ) %>%
      filter(!is.na(phylo_name))
    
    if (!is.null(cfg$manual_patch_file) && file.exists(cfg$manual_patch_file)) {
      manual_patches <- read_csv(cfg$manual_patch_file, show_col_types = FALSE) %>%
        transmute(
          IUCN_Name = clean_names(scientific_name),
          phylo_name = to_phylo_style(manual_tree_name),
          Match_Type = "Manual Patch"
        ) %>%
        filter(!is.na(phylo_name), phylo_name %in% valid_tips_df$Tree_Tip)
      
      if (nrow(manual_patches) > 0) {
        manual_patches <- manual_patches %>%
          left_join(
            iucn_clean %>%
              select(IUCN_Name, Family, Order, RedList_Cat, Threat_Logging, Threat_Plantation),
            by = "IUCN_Name"
          )
        matched_dat <- bind_rows(
          matched_dat %>%
            select(phylo_name, IUCN_Name, Family, Order, RedList_Cat, Threat_Logging, Threat_Plantation, Match_Type),
          manual_patches %>%
            select(phylo_name, IUCN_Name, Family, Order, RedList_Cat, Threat_Logging, Threat_Plantation, Match_Type)
        )
      }
    }
    
    matched_dat <- matched_dat %>%
      group_by(phylo_name) %>%
      arrange(
        factor(Match_Type, levels = c("Direct", "Synonym", "Manual Patch", "Fuzzy")),
        IUCN_Name
      ) %>%
      summarise(across(everything(), first_non_na), .groups = "drop") %>%
      mutate(
        Threat_Logging = replace_na(Threat_Logging, FALSE),
        Threat_Plantation = replace_na(Threat_Plantation, FALSE)
      ) %>%
      select(
        phylo_name, IUCN_Name, Family, Order, RedList_Cat,
        Threat_Logging, Threat_Plantation, Match_Type
      )
    
    write_csv(matched_dat, matched_data_file)
    cat("Saved matched assessment table to: ", matched_data_file, "\n", sep = "")
  }
  
  if (identical(cfg$taxon, "fish")) {
    fish_taxonomy <- load_fishtree_taxonomy(
      path = cfg$taxonomy_file,
      tree_tips = tree$tip.label
    )
    matched_dat <- apply_fishtree_taxonomy(
      matched_dat = matched_dat,
      taxonomy = fish_taxonomy,
      label = cfg$taxon
    )
    write_csv(matched_dat, matched_data_file)
    cat("Updated fish matched table with FishTree family/order taxonomy: ", matched_data_file, "\n", sep = "")
  }
  
  cat(sprintf("   > Matched species available for plotting: %d\n", nrow(matched_dat)))
  
  # Step 3: build the tip-level threat table
  cat("--- Step 3: Joining matched data to all tree tips ---\n")
  tip_data <- tibble(label = tree$tip.label) %>%
    mutate(label = to_phylo_style(label)) %>%
    left_join(matched_dat, by = c("label" = "phylo_name")) %>%
    mutate(
      Threat_Logging = replace_na(Threat_Logging, FALSE),
      Threat_Plantation = replace_na(Threat_Plantation, FALSE),
      Wood_Threat = Threat_Logging | Threat_Plantation,
      tip_category = case_when(
        Threat_Logging & Threat_Plantation ~ "Both",
        Threat_Logging & !Threat_Plantation ~ "Logging only",
        !Threat_Logging & Threat_Plantation ~ "Plantation only",
        TRUE ~ "None"
      )
    )
  
  tip_data <- tip_data %>%
    filter(!is.na(IUCN_Name))
  
  if (nrow(tip_data) == 0) {
    stop("No tree tips with IUCN assessments were matched for ", cfg$taxon)
  }
  
  tree <- keep.tip(tree, tip_data$label)
  tip_data <- tip_data %>%
    arrange(match(label, tree$tip.label))
  
  if (!identical(tip_data$label, tree$tip.label)) {
    stop("Tip metadata are not aligned with the pruned tree for ", cfg$taxon)
  }
  
  tip_data$tip_category <- factor(
    tip_data$tip_category,
    levels = c("None", "Logging only", "Plantation only", "Both")
  )
  
  cat(sprintf("   > Tips with an IUCN match retained for plotting: %d\n", nrow(tip_data)))
  cat(sprintf("   > Logging-threatened tips: %d\n", sum(tip_data$Threat_Logging, na.rm = TRUE)))
  cat(sprintf("   > Plantation-threatened tips: %d\n", sum(tip_data$Threat_Plantation, na.rm = TRUE)))
  cat(sprintf("   > Branch line size used: %.2f\n", branch_line_size))
  
  # Step 3b: compute phylogenetic signal metrics on the matched, plotted tree
  d_statistics_file <- file.path(cfg$processed_dir, paste0(cfg$taxon, "_phylogenetic_d_statistics.csv"))
  phylo_context <- build_phylo_distance_context(tree)
  cat("--- Step 3b: Calculating D, SES.MPD, and SES.MNTD for threat conservatism ---\n")
  if (!force_rebuild_d_statistics && file.exists(d_statistics_file)) {
    cached_d_stat_results <- read_csv(d_statistics_file, show_col_types = FALSE)
    expected_traits <- c("Wood_Threat", "Threat_Logging", "Threat_Plantation")
    required_signal_columns <- c(
      "taxon", "trait", "trait_label", "n_tips", "n_threatened", "threat_prevalence",
      "d_estimate", "p_random", "p_brownian", "observed_sum", "mean_random_sum",
      "mean_brownian_sum", "ses_mpd", "ses_mpd_obs", "ses_mpd_rand_mean",
      "ses_mpd_rand_sd", "ses_mpd_p", "ses_mntd", "ses_mntd_obs",
      "ses_mntd_rand_mean", "ses_mntd_rand_sd", "ses_mntd_p", "permutations",
      "signal_class", "clustered_vs_random", "note"
    )
    cache_matches_current_tree <- all(cached_d_stat_results$n_tips == nrow(tip_data)) &&
      setequal(cached_d_stat_results$trait, expected_traits) &&
      all(required_signal_columns %in% names(cached_d_stat_results))
    
    if (cache_matches_current_tree) {
      d_stat_results <- cached_d_stat_results
      cat("Loading cached phylogenetic signal table from: ", d_statistics_file, "\n", sep = "")
    } else {
      cat(
        "Cached phylogenetic signal metrics are out of sync with the current matched tree or missing new SES columns (",
        paste(unique(cached_d_stat_results$n_tips), collapse = ", "),
        " cached tips vs ", nrow(tip_data),
        " current tips); recomputing.\n",
        sep = ""
      )
      d_stat_results <- summarise_taxon_phylo_d(
        tree = tree,
        tip_data = tip_data,
        taxon = cfg$taxon,
        permutations = d_stat_permutations,
        seed_base = d_stat_seed + (match(cfg$taxon, names(taxa_configs)) * 100L),
        phylo_context = phylo_context
      )
      saved_d_statistics_file <- safe_write_csv(d_stat_results, d_statistics_file)
      cat("Saved phylogenetic signal table to: ", saved_d_statistics_file, "\n", sep = "")
    }
  } else {
    d_stat_results <- summarise_taxon_phylo_d(
      tree = tree,
      tip_data = tip_data,
      taxon = cfg$taxon,
      permutations = d_stat_permutations,
      seed_base = d_stat_seed + (match(cfg$taxon, names(taxa_configs)) * 100L),
      phylo_context = phylo_context
    )
    saved_d_statistics_file <- safe_write_csv(d_stat_results, d_statistics_file)
    cat("Saved phylogenetic signal table to: ", saved_d_statistics_file, "\n", sep = "")
  }
  
  print(
    d_stat_results %>%
      select(
        taxon, trait_label, n_tips, n_threatened, threat_prevalence,
        d_estimate, p_random, p_brownian, ses_mpd, ses_mpd_p,
        ses_mntd, ses_mntd_p, signal_class, clustered_vs_random
      )
  )
  
  analysis_tree_depth <- get_tree_depth(tree)
  
  # Step 4: summarise threat across descendant clades
  cat("--- Step 4: Computing descendant-tip threat proportions ---\n")
  
  # Step 5: draw the circular tree and heatmap ring
  cat("--- Step 5: Building circular phylogram ---\n")
  
  # Step 6: label the top 10 families by threatened count and by threatened proportion.
  # Simpler implementation: family name + stats are combined into one plotmath label
  # with atop(), so they always move together under geom_text_repel().
  cat("--- Step 6: Adding top 10 family labels by count and proportion ---\n")
  plot_bundle <- build_phylogram_plot_set(
    tree = tree,
    tip_data = tip_data,
    cfg = cfg
  )
  p_no_family_labels <- plot_bundle$no_family_labels_plot
  p3 <- plot_bundle$family_labels_plot
  
  # Step 7: save the final figures
  cat("--- Step 7: Saving figures ---\n")
  no_family_label_output_file <- append_file_stem_suffix(cfg$output_file, "_no_family_labels")
  ensure_parent_dir(no_family_label_output_file)
  if (file.exists(no_family_label_output_file) && !force_rebuild_figure) {
    cat("Overwriting existing no-family-label figure: ", no_family_label_output_file, "\n", sep = "")
  }
  ggsave(
    filename = no_family_label_output_file,
    plot = p_no_family_labels,
    width = figure_width,
    height = figure_height,
    dpi = figure_dpi,
    bg = "white"
  )
  cat("Saved no-family-label figure to: ", no_family_label_output_file, "\n", sep = "")
  
  if (file.exists(cfg$output_file) && !force_rebuild_figure) {
    cat("Overwriting existing figure: ", cfg$output_file, "\n", sep = "")
  }
  ggsave(
    filename = cfg$output_file,
    plot = p3,
    width = figure_width,
    height = figure_height,
    dpi = figure_dpi,
    bg = "white"
  )
  
  cat("Saved figure to: ", cfg$output_file, "\n", sep = "")
  cat("--- ", toupper(cfg$taxon), " SCRIPT COMPLETE ---\n", sep = "")
  
  invisible(list(
    d_statistics = d_stat_results,
    tree_depth = analysis_tree_depth,
    plot_input = list(
      cfg = cfg,
      tree = tree,
      tip_data = tip_data
    )
  ))
}

build_combined_phylogram_figure <- function(
  taxon_results,
  plot_type = c("family_labels_plot", "no_family_labels_plot")
) {
  plot_type <- match.arg(plot_type)
  target_tree_depth <- max(purrr::map_dbl(taxon_results, "tree_depth"), na.rm = TRUE)
  plot_margin <- if (identical(plot_type, "family_labels_plot")) {
    margin(20, 180, 20, 120)
  } else {
    margin(20, 20, 20, 20)
  }
  
  plot_list <- purrr::map(
    taxon_results,
    function(result) {
      cfg <- result$plot_input$cfg
      plot_bundle <- build_phylogram_plot_set(
        tree = result$plot_input$tree,
        tip_data = result$plot_input$tip_data,
        cfg = cfg,
        panel_title = unname(taxon_plot_titles[[cfg$taxon]]),
        target_tree_depth = target_tree_depth,
        plot_margin = plot_margin,
        print_family_tables = FALSE
      )
      plot_bundle[[plot_type]]
    }
  )
  
  ncol_layout <- min(3L, length(plot_list))
  while (length(plot_list) %% ncol_layout != 0L) {
    plot_list[[length(plot_list) + 1L]] <- patchwork::plot_spacer()
  }
  
  subtitle_text <- if (identical(plot_type, "family_labels_plot")) {
    "All panels drawn to a shared tree depth, with family labels."
  } else {
    "All panels drawn to a shared tree depth, without family labels."
  }
  
  wrap_plots(plotlist = plot_list, ncol = ncol_layout, guides = "collect") +
    plot_annotation(
      title = "Wood-production phylogenies",
      subtitle = subtitle_text,
      theme = theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5)
      )
    ) &
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9.5)
    )
}

save_combined_phylogram_figures <- function(taxon_results) {
  combined_output_file <- here::here(
    "outputs", "extent", "figures", "combined_radial_phylogenies_same_size.png"
  )
  combined_no_family_output_file <- here::here(
    "outputs", "extent", "figures", "combined_radial_phylogenies_same_size_no_family_labels.png"
  )
  ensure_parent_dir(combined_output_file)
  
  combined_label_plot <- build_combined_phylogram_figure(
    taxon_results = taxon_results,
    plot_type = "family_labels_plot"
  )
  ggsave(
    filename = combined_output_file,
    plot = combined_label_plot,
    width = combined_figure_width,
    height = combined_figure_height,
    dpi = figure_dpi,
    bg = "white"
  )
  cat("Saved combined labelled phylogeny figure to: ", combined_output_file, "\n", sep = "")
  
  combined_no_family_plot <- build_combined_phylogram_figure(
    taxon_results = taxon_results,
    plot_type = "no_family_labels_plot"
  )
  ggsave(
    filename = combined_no_family_output_file,
    plot = combined_no_family_plot,
    width = combined_figure_width,
    height = combined_figure_height,
    dpi = figure_dpi,
    bg = "white"
  )
  cat(
    "Saved combined no-family-label phylogeny figure to: ",
    combined_no_family_output_file,
    "\n",
    sep = ""
  )
}

taxon_results <- purrr::map(
  requested_taxa,
  ~ build_taxon_phylogram(taxa_configs[[.x]])
)

if (length(taxon_results) > 1L) {
  save_combined_phylogram_figures(taxon_results)
} else {
  cat("Skipping combined phylogeny figure because only one taxon was requested.\n")
}

combined_d_statistics <- bind_rows(purrr::map(taxon_results, "d_statistics"))
combined_d_file <- here::here("outputs", "extent", "tables", "phylogenetic_d_statistics.csv")
saved_combined_d_file <- safe_write_csv(combined_d_statistics, combined_d_file)
cat("Saved combined phylogenetic signal table to: ", saved_combined_d_file, "\n", sep = "")
