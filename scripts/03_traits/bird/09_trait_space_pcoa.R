rm(list = ls())
gc()

# =============================================================================
# BIRD TRAIT SPACE ORDINATION
#
# Builds a mixed-trait Gower dissimilarity matrix, projects species into
# continuous trait space with PCoA, tests group separation with PERMANOVA, and
# fits traits back onto the ordination with vegan::envfit.
#
# Boundaries are drawn using kernel density estimation (KDE): a filled 95%
# region shows the broad trait envelope and a solid 50% contour shows the
# core trait space for each threat group. Both are robust to outlier species
# in a way that convex hulls are not.
#
# Run from the project root:
#   Rscript scripts/03_traits/bird/09_trait_space_pcoa.R
#
# Optional environment variables:
#   BIRD_TRAIT_PCOA_BOUNDARY_STYLE      kde, ellipse, both, none
#                                         Default: kde
#                                         "both" overlays KDE + t-ellipse
#   BIRD_TRAIT_PCOA_KDE_LEVELS          comma-separated probabilities for KDE
#                                         contours, default "0.50,0.95"
#   BIRD_TRAIT_PCOA_KDE_GRIDSIZE        integer grid resolution, default 200
#   BIRD_TRAIT_PCOA_INCLUDE_RANGE_SIZE  TRUE/FALSE
#   BIRD_TRAIT_PCOA_PERMUTATIONS        integer, default 999
#   BIRD_TRAIT_PCOA_MAX_SPECIES         integer subsample for quick checks
#   BIRD_TRAIT_PCOA_SAVE_DISTANCE       TRUE/FALSE, default FALSE
#   BIRD_TRAIT_PCOA_OUTPUT_TAG          optional suffix for test output folders
#   BIRD_TRAIT_PCOA_ENVFIT_MAX_FACTOR_LABELS integer, default 12
# =============================================================================

project_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib, winslash = "/", mustWork = FALSE), .libPaths()))
}

required_packages <- c("dplyr", "tidyr", "ggplot2", "cluster", "vegan", "tibble", "MASS")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this ordination script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(vegan)
  library(tibble)
  library(MASS)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

# =============================================================================
# Utility functions
# =============================================================================

to_bool <- function(x, default = FALSE) {
  if (length(x) == 0 || is.na(x) || x == "") {
    return(default)
  }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

env_number <- function(name, default = Inf) {
  value <- Sys.getenv(name, unset = "")
  if (value == "") {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    return(default)
  }
  parsed
}

pick_first_column <- function(df, candidates, path, required = TRUE) {
  matching <- candidates[candidates %in% names(df)]
  if (length(matching) > 0) {
    return(matching[[1]])
  }
  
  if (required) {
    stop(
      "Missing required column in ", path, ": ",
      paste(candidates, collapse = " / "),
      call. = FALSE
    )
  }
  
  NULL
}

clean_category <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

as_flag <- function(x) {
  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }
  
  if (is.numeric(x)) {
    return(!is.na(x) & x != 0)
  }
  
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

as_binary01 <- function(x) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_integer_, as.integer(x)))
  }
  
  if (is.numeric(x)) {
    return(ifelse(is.na(x), NA_integer_, as.integer(x != 0)))
  }
  
  x_chr <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    x_chr %in% c("true", "t", "1", "yes", "y", "present") ~ 1L,
    x_chr %in% c("false", "f", "0", "no", "n", "absent") ~ 0L,
    TRUE ~ NA_integer_
  )
}

has_variation <- function(x) {
  x <- x[!is.na(x)]
  length(unique(x)) > 1
}

# =============================================================================
# KDE boundary functions
# =============================================================================

# Compute kernel density contour polygons for each threat group.
#
# For each probability level p (e.g. 0.50, 0.95), the function finds the
# density threshold that encloses exactly p of the total probability mass and
# extracts the corresponding contour lines. This makes the boundaries robust
# to outlier species: an isolated outlier shifts the density surface only
# marginally and does not stretch the visible boundary.
#
# Returns a data frame with columns:
#   PCoA1, PCoA2, threat_group, prob_level, piece
# suitable for geom_polygon / geom_path.
make_kde_contours <- function(df,
                              prob_levels = c(0.50, 0.95),
                              gridsize   = 200L) {
  df <- df %>% filter(is.finite(PCoA1), is.finite(PCoA2))
  
  # Expand the grid limits slightly beyond the data range so that contours
  # at the edge of the point cloud are not clipped.
  x_pad <- diff(range(df$PCoA1)) * 0.15
  y_pad <- diff(range(df$PCoA2)) * 0.15
  global_lims <- c(
    range(df$PCoA1) + c(-x_pad, x_pad),
    range(df$PCoA2) + c(-y_pad, y_pad)
  )
  
  bind_rows(lapply(unique(as.character(df$threat_group)), function(grp) {
    sub <- df %>% filter(as.character(threat_group) == grp)
    
    if (nrow(sub) < 10) {
      warning(sprintf(
        "Group '%s' has only %d species; KDE contours may be unreliable. Skipping.",
        grp, nrow(sub)
      ))
      return(NULL)
    }
    
    # Bandwidth selection: use Sheather-Jones plug-in (MASS::bandwidth.SJ)
    # with a fallback to Silverman's rule-of-thumb if SJ fails (can happen
    # with low-n groups or degenerate distributions).
    bw_x <- tryCatch(
      MASS::bandwidth.SJ(sub$PCoA1),
      error = function(e) MASS::bandwidth.nrd(sub$PCoA1)
    )
    bw_y <- tryCatch(
      MASS::bandwidth.SJ(sub$PCoA2),
      error = function(e) MASS::bandwidth.nrd(sub$PCoA2)
    )
    bw_x <- max(bw_x, 1e-6)
    bw_y <- max(bw_y, 1e-6)
    
    kde <- MASS::kde2d(
      sub$PCoA1, sub$PCoA2,
      h    = c(bw_x, bw_y),
      n    = gridsize,
      lims = global_lims
    )
    
    # Flatten density grid for probability-mass thresholding.
    z_vec   <- as.vector(kde$z)
    z_total <- sum(z_vec)
    
    bind_rows(lapply(prob_levels, function(lvl) {
      # Sort densities descending; walk down until cumulative mass >= lvl.
      z_sorted   <- sort(z_vec, decreasing = TRUE)
      cumulative <- cumsum(z_sorted) / z_total
      threshold  <- z_sorted[which(cumulative >= lvl)[1]]
      
      # grDevices::contourLines returns a list of (x, y, level) paths.
      cl <- grDevices::contourLines(kde$x, kde$y, kde$z, levels = threshold)
      
      if (length(cl) == 0) {
        return(NULL)
      }
      
      bind_rows(lapply(seq_along(cl), function(i) {
        tibble(
          PCoA1        = cl[[i]]$x,
          PCoA2        = cl[[i]]$y,
          threat_group = grp,
          prob_level   = lvl,
          piece        = i
        )
      }))
    }))
  }))
}

# Build a tidy summary of KDE overlap between the two threat groups.
# Reports the proportion of each group's 50% and 95% density region that
# falls within the other group's corresponding region (estimated by point
# membership, which is a reasonable approximation for ordination data).
make_kde_overlap_summary <- function(kde_data, ordination_scores) {
  if (is.null(kde_data) || nrow(kde_data) == 0) {
    return(tibble())
  }
  
  groups     <- unique(kde_data$threat_group)
  prob_levels <- unique(kde_data$prob_level)
  
  bind_rows(lapply(prob_levels, function(lvl) {
    # For each species, check whether it falls inside the contour polygon of
    # each group at this probability level. Uses sp::point.in.polygon if
    # available, otherwise skips (informational only — not used in main
    # analysis).
    contour_sub <- kde_data %>% filter(prob_level == lvl)
    if (nrow(contour_sub) == 0) return(NULL)
    
    bind_rows(lapply(groups, function(focal) {
      focal_scores <- ordination_scores %>%
        filter(as.character(threat_group) == focal)
      if (nrow(focal_scores) == 0) return(NULL)
      
      bind_rows(lapply(groups, function(target) {
        poly_sub <- contour_sub %>%
          filter(as.character(threat_group) == target)
        if (nrow(poly_sub) == 0) return(NULL)
        
        # Aggregate all pieces into one bounding test (approximate).
        in_poly <- tryCatch({
          pieces <- split(poly_sub, poly_sub$piece)
          inside <- vapply(focal_scores$species_id, function(sp) {
            x <- focal_scores$PCoA1[focal_scores$species_id == sp]
            y <- focal_scores$PCoA2[focal_scores$species_id == sp]
            any(vapply(pieces, function(p) {
              # Simple winding-number test using grDevices
              pt <- sp::point.in.polygon(x, y, p$PCoA1, p$PCoA2)
              pt > 0
            }, logical(1)))
          }, logical(1))
          mean(inside)
        }, error = function(e) NA_real_)
        
        tibble(
          focal_group  = focal,
          target_group = target,
          prob_level   = lvl,
          prop_inside  = in_poly
        )
      }))
    }))
  }))
}

scale_arrow_endpoints <- function(arrows, site_scores, fraction = 0.28) {
  if (nrow(arrows) == 0 || all(arrows == 0)) {
    return(arrows)
  }
  
  x_range <- diff(range(site_scores$PCoA1, na.rm = TRUE))
  y_range <- diff(range(site_scores$PCoA2, na.rm = TRUE))
  
  if (!is.finite(x_range) || !is.finite(y_range) || x_range == 0 || y_range == 0) {
    return(arrows)
  }
  
  x_ratio    <- max(abs(arrows[, 1]), na.rm = TRUE) / x_range
  y_ratio    <- max(abs(arrows[, 2]), na.rm = TRUE) / y_range
  denominator <- max(x_ratio, y_ratio, na.rm = TRUE)
  
  if (!is.finite(denominator) || denominator == 0) {
    return(arrows)
  }
  
  arrows * (fraction / denominator)
}

make_group_palette <- function(groups) {
  groups <- as.character(groups)
  palette <- c(
    "Logging"            = "#D9BE00",
    "Timber plantations" = "#0088D9"
  )
  missing_groups <- setdiff(groups, names(palette))
  if (length(missing_groups) > 0) {
    fallback <- grDevices::hcl.colors(length(missing_groups), palette = "Dark 3")
    names(fallback) <- missing_groups
    palette <- c(palette, fallback)
  }
  
  palette[groups]
}

# =============================================================================
# Configuration
# =============================================================================

set.seed(20260524)

boundary_style <- Sys.getenv("BIRD_TRAIT_PCOA_BOUNDARY_STYLE", unset = "kde")
boundary_style <- match.arg(boundary_style, c("kde", "ellipse", "both", "none"))

# Parse KDE-specific settings
kde_levels_raw <- Sys.getenv("BIRD_TRAIT_PCOA_KDE_LEVELS", unset = "0.50,0.95")
kde_levels <- suppressWarnings(as.numeric(strsplit(kde_levels_raw, ",")[[1]]))
kde_levels <- sort(unique(kde_levels[is.finite(kde_levels) & kde_levels > 0 & kde_levels < 1]))
if (length(kde_levels) == 0) {
  kde_levels <- c(0.50, 0.95)
}

kde_gridsize <- as.integer(env_number("BIRD_TRAIT_PCOA_KDE_GRIDSIZE", 200))
if (is.na(kde_gridsize) || kde_gridsize < 50L) {
  kde_gridsize <- 200L
}

permutations <- as.integer(env_number("BIRD_TRAIT_PCOA_PERMUTATIONS", 999))
if (is.na(permutations) || permutations < 0) {
  permutations <- 999L
}

include_range_size      <- to_bool(Sys.getenv("BIRD_TRAIT_PCOA_INCLUDE_RANGE_SIZE", unset = "TRUE"))
save_distance           <- to_bool(Sys.getenv("BIRD_TRAIT_PCOA_SAVE_DISTANCE",       unset = "FALSE"))
max_species             <- env_number("BIRD_TRAIT_PCOA_MAX_SPECIES",           Inf)
envfit_plot_p_max       <- env_number("BIRD_TRAIT_PCOA_ENVFIT_P_MAX",          0.05)
envfit_max_factor_labels <- env_number("BIRD_TRAIT_PCOA_ENVFIT_MAX_FACTOR_LABELS", 12)
output_tag              <- Sys.getenv("BIRD_TRAIT_PCOA_OUTPUT_TAG", unset = "")
output_tag              <- gsub("[^A-Za-z0-9_-]+", "_", trimws(output_tag))
output_suffix           <- if (output_tag == "") "" else paste0("_", output_tag)

# =============================================================================
# File paths
# =============================================================================

model_path <- file.path(
  "data", "processed", "traits", "bird", "bird_trait_model_final_data.csv"
)
matched_path <- file.path(
  "data", "processed", "traits", "bird", "bird_trait_data_matched.csv"
)
assessment_lookup_path <- file.path(
  "data", "processed", "assessments", "bird_assessment_taxonomy_lookup.rds"
)

figure_dir    <- file.path("outputs", "traits", "figures",   paste0("bird_trait_space_pcoa", output_suffix))
table_dir     <- file.path("outputs", "traits", "tables",    paste0("bird_trait_space_pcoa", output_suffix))
processed_dir <- file.path("data", "processed", "traits", "bird", paste0("ordination", output_suffix))

dir.create(figure_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(model_path)) {
  stop("Missing model data file: ", model_path, call. = FALSE)
}
if (!file.exists(matched_path)) {
  stop("Missing matched bird trait file: ", matched_path, call. = FALSE)
}

# =============================================================================
# Load and prepare data
# =============================================================================

cat("--- Loading bird trait data ---\n")
bird_model_raw <- read.csv(model_path,   stringsAsFactors = FALSE, check.names = TRUE)
bird_matched   <- read.csv(matched_path, stringsAsFactors = FALSE, check.names = FALSE)

species_col <- pick_first_column(
  bird_model_raw,
  c("spatial_name", "Trait_Name", "species"),
  model_path
)
phylo_col <- pick_first_column(
  bird_model_raw,
  c("phylo_name", "species"),
  model_path,
  required = FALSE
)
logging_col <- pick_first_column(
  bird_model_raw,
  c("Threat_Logging", "Threat_Logging.x"),
  model_path
)
plantation_col <- pick_first_column(
  bird_model_raw,
  c("Threat_Plantation", "Threat_Plantation.x"),
  model_path
)

matched_species_col <- pick_first_column(
  bird_matched,
  c("Trait_Name", "IUCN_Name", "scientific_name"),
  matched_path
)
matched_category_col <- pick_first_column(
  bird_matched,
  c("RedList_Cat", "category", "X2024.IUCN.Red.List.category"),
  matched_path
)

status_lookup <- bird_matched %>%
  transmute(
    spatial_name      = .data[[matched_species_col]],
    redlist_category  = clean_category(.data[[matched_category_col]])
  ) %>%
  filter(!is.na(spatial_name), spatial_name != "") %>%
  distinct(spatial_name, .keep_all = TRUE)

bird_data <- bird_model_raw %>%
  mutate(
    species    = .data[[species_col]],
    species_id = if (!is.null(phylo_col)) .data[[phylo_col]] else .data[[species_col]],
    species_id = if_else(
      is.na(species_id) | species_id == "",
      species,
      species_id
    ),
    species_id       = make.unique(gsub(" ", "_", species_id)),
    Threat_Logging   = as_flag(.data[[logging_col]]),
    Threat_Plantation = as_flag(.data[[plantation_col]])
  ) %>%
  left_join(status_lookup, by = c("species" = "spatial_name"))

if (any(is.na(bird_data$redlist_category)) && file.exists(assessment_lookup_path)) {
  assessment_lookup <- readRDS(assessment_lookup_path) %>%
    transmute(
      species                     = scientific_name,
      redlist_category_assessment = clean_category(RedList_Cat)
    ) %>%
    distinct(species, .keep_all = TRUE)
  
  bird_data <- bird_data %>%
    left_join(assessment_lookup, by = "species") %>%
    mutate(redlist_category = coalesce(redlist_category, redlist_category_assessment)) %>%
    select(-redlist_category_assessment)
}

# =============================================================================
# Trait columns
# =============================================================================

continuous_trait_cols <- c(
  "ESI_scaled",
  "Body_Mass_log",
  "Habitat_Breadth_log",
  "Clutch_Size_log",
  "HWI_log"
)

if (include_range_size) {
  continuous_trait_cols <- c(continuous_trait_cols, "log_range_size")
}

binary_trait_cols <- c(
  "Island_Status",
  "Is_Tree_Nester",
  "Is_Stump_Nester",
  "Is_Ground_Nester",
  "Forest",
  "Savanna",
  "Shrubland",
  "Grassland",
  "Wetlands",
  "Rocky",
  "Desert",
  "Marine",
  "Artificial",
  "Diet_Invertebrate",
  "Diet_Fruit",
  "Diet_Nectar",
  "Diet_Seed",
  "Diet_Vertebrate",
  "Diet_Fish",
  "Diet_Scavenger",
  "Diet_PlantMaterial",
  "Diet_Misc"
)

continuous_trait_cols <- intersect(continuous_trait_cols, names(bird_data))
binary_trait_cols     <- intersect(binary_trait_cols,     names(bird_data))
all_trait_cols        <- c(continuous_trait_cols, binary_trait_cols)

if (length(all_trait_cols) < 2) {
  stop("Fewer than two usable trait columns were found in ", model_path, call. = FALSE)
}

# =============================================================================
# Build Gower dissimilarity matrix
# =============================================================================

cat("--- Preparing Gower trait matrix ---\n")
trait_data <- bird_data %>%
  mutate(
    across(all_of(continuous_trait_cols), ~ suppressWarnings(as.numeric(.x))),
    across(all_of(binary_trait_cols),     as_binary01)
  ) %>%
  filter(if_all(all_of(all_trait_cols), ~ !is.na(.x)))

if (is.finite(max_species) && nrow(trait_data) > max_species) {
  cat(sprintf(
    "Subsampling %d of %d species for this run.\n",
    as.integer(max_species),
    nrow(trait_data)
  ))
  trait_data <- trait_data %>%
    slice_sample(n = as.integer(max_species))
}

variation_flag        <- vapply(trait_data[all_trait_cols], has_variation, logical(1))
dropped_trait_cols    <- all_trait_cols[!variation_flag]
all_trait_cols        <- all_trait_cols[variation_flag]
continuous_trait_cols <- intersect(continuous_trait_cols, all_trait_cols)
binary_trait_cols     <- intersect(binary_trait_cols,     all_trait_cols)

if (length(all_trait_cols) < 2) {
  stop("Fewer than two variable trait columns remain after filtering.", call. = FALSE)
}

species_ids <- trait_data$species_id
group_counts <- tibble(
  threat_flag = c("Logging", "Timber plantations"),
  n = c(
    sum(trait_data$Threat_Logging),
    sum(trait_data$Threat_Plantation)
  ),
  note = "Groups are non-exclusive; species flagged for both are counted in both rows."
)

write.csv(
  group_counts,
  file.path(table_dir, "bird_trait_space_group_counts.csv"),
  row.names = FALSE
)

if (length(dropped_trait_cols) > 0) {
  write.csv(
    data.frame(dropped_trait = dropped_trait_cols),
    file.path(table_dir, "bird_trait_space_dropped_traits.csv"),
    row.names = FALSE
  )
}

gower_input          <- trait_data[, all_trait_cols, drop = FALSE]
rownames(gower_input) <- species_ids

# Habitat, diet, and nesting flags are treated as asymmetric binary variables:
# a shared presence increases similarity, but a shared absence does not.
asymmetric_binary_cols <- setdiff(binary_trait_cols, "Island_Status")
symmetrical_binary_cols <- intersect(binary_trait_cols, "Island_Status")

gower_type <- list()
if (length(asymmetric_binary_cols) > 0) gower_type$asymm <- asymmetric_binary_cols
if (length(symmetrical_binary_cols) > 0) gower_type$symm  <- symmetrical_binary_cols

gower_dist <- if (length(gower_type) > 0) {
  cluster::daisy(gower_input, metric = "gower", type = gower_type)
} else {
  cluster::daisy(gower_input, metric = "gower")
}

attr(gower_dist, "Labels") <- species_ids

if (save_distance) {
  saveRDS(
    gower_dist,
    file.path(processed_dir, "bird_trait_space_gower_distance.rds")
  )
}

# =============================================================================
# PCoA
# =============================================================================

cat("--- Running PCoA on Gower dissimilarities ---\n")
pcoa_fit <- vegan::wcmdscale(
  as.dist(gower_dist),
  k   = 2,
  eig = TRUE,
  add = "lingoes"
)

pcoa_points         <- as.data.frame(pcoa_fit$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")

positive_eigenvalues <- pcoa_fit$eig[pcoa_fit$eig > 0]
axis_variance        <- pcoa_fit$eig[seq_len(2)] / sum(positive_eigenvalues) * 100

ordination_scores <- bind_cols(
  trait_data %>%
    transmute(
      species_id,
      species,
      redlist_category,
      Threat_Logging,
      Threat_Plantation
    ),
  pcoa_points
)

write.csv(
  ordination_scores,
  file.path(processed_dir, "bird_trait_space_pcoa_scores.csv"),
  row.names = FALSE
)

eigen_table <- tibble(
  axis       = seq_along(pcoa_fit$eig),
  eigenvalue = pcoa_fit$eig,
  positive   = pcoa_fit$eig > 0,
  variance_explained = if_else(
    pcoa_fit$eig > 0,
    pcoa_fit$eig / sum(positive_eigenvalues),
    NA_real_
  )
) %>%
  mutate(cumulative_variance_explained = cumsum(replace_na(variance_explained, 0)))

write.csv(
  eigen_table,
  file.path(table_dir, "bird_trait_space_pcoa_eigenvalues.csv"),
  row.names = FALSE
)

# =============================================================================
# PERMANOVA
# =============================================================================

cat("--- Running PERMANOVA ---\n")
permanova_data <- ordination_scores %>%
  transmute(
    Threat_Logging    = factor(Threat_Logging,    levels = c(FALSE, TRUE)),
    Threat_Plantation = factor(Threat_Plantation, levels = c(FALSE, TRUE))
  )

if (
  nlevels(droplevels(permanova_data$Threat_Logging))    >= 2 ||
  nlevels(droplevels(permanova_data$Threat_Plantation)) >= 2
) {
  permanova_fit <- vegan::adonis2(
    as.dist(gower_dist) ~ Threat_Logging + Threat_Plantation,
    data         = permanova_data,
    permutations = permutations,
    by           = "margin"
  )
  
  permanova_table <- as.data.frame(permanova_fit) %>%
    rownames_to_column("term")
} else {
  warning("PERMANOVA skipped because neither threat flag has both TRUE and FALSE values.")
  permanova_table <- tibble(
    term      = c("Threat_Logging", "Threat_Plantation"),
    Df        = NA_real_,
    SumOfSqs  = NA_real_,
    R2        = NA_real_,
    F         = NA_real_,
    `Pr(>F)`  = NA_real_
  )
}

write.csv(
  permanova_table,
  file.path(table_dir, "bird_trait_space_permanova.csv"),
  row.names = FALSE
)

# =============================================================================
# Betadisper (homogeneity of group dispersions)
# =============================================================================

dispersion_table <- bind_rows(lapply(
  c("Threat_Logging", "Threat_Plantation"),
  function(flag_col) {
    tryCatch(
      {
        grouping <- droplevels(permanova_data[[flag_col]])
        if (nlevels(grouping) < 2) {
          return(tibble(
            predictor = flag_col,
            term      = "skipped",
            message   = "Fewer than two levels remain."
          ))
        }
        
        dispersion_fit  <- vegan::betadisper(as.dist(gower_dist), grouping, type = "centroid")
        dispersion_test <- vegan::permutest(dispersion_fit, permutations = permutations)
        as.data.frame(dispersion_test$tab) %>%
          rownames_to_column("term") %>%
          mutate(predictor = flag_col, .before = 1)
      },
      error = function(e) {
        tibble(predictor = flag_col, term = "betadisper_error", message = conditionMessage(e))
      }
    )
  }
))

write.csv(
  dispersion_table,
  file.path(table_dir, "bird_trait_space_betadisper.csv"),
  row.names = FALSE
)

# =============================================================================
# envfit
# =============================================================================

cat("--- Running vegan::envfit ---\n")
envfit_input <- trait_data[, all_trait_cols, drop = FALSE]
envfit_input <- envfit_input %>%
  mutate(across(
    all_of(binary_trait_cols),
    ~ factor(if_else(.x == 1, "Yes", "No"), levels = c("No", "Yes"))
  ))
rownames(envfit_input) <- species_ids

ordination_matrix          <- as.matrix(ordination_scores[, c("PCoA1", "PCoA2")])
rownames(ordination_matrix) <- species_ids

envfit_fit <- vegan::envfit(
  ordination_matrix,
  envfit_input,
  permutations = permutations,
  na.rm        = TRUE
)

saveRDS(
  envfit_fit,
  file.path(processed_dir, "bird_trait_space_envfit.rds")
)

vector_table <- tibble()
if (!is.null(envfit_fit$vectors) && nrow(envfit_fit$vectors$arrows) > 0) {
  vector_base      <- envfit_fit$vectors$arrows
  vector_base      <- sweep(vector_base, 1, sqrt(envfit_fit$vectors$r), "*")
  vector_endpoints <- scale_arrow_endpoints(vector_base, ordination_scores)
  
  vector_table <- as.data.frame(vector_endpoints) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("trait") %>%
    mutate(
      r2      = envfit_fit$vectors$r[trait],
      p_value = envfit_fit$vectors$pvals[trait]
    ) %>%
    arrange(p_value, desc(r2))
}

write.csv(
  vector_table,
  file.path(table_dir, "bird_trait_space_envfit_continuous_vectors.csv"),
  row.names = FALSE
)

factor_centroid_table <- bind_rows(lapply(binary_trait_cols, function(trait) {
  levels_here <- levels(envfit_input[[trait]])
  
  bind_rows(lapply(levels_here, function(level) {
    rows <- envfit_input[[trait]] == level
    tibble(
      trait    = trait,
      level    = level,
      n        = sum(rows, na.rm = TRUE),
      PCoA1    = mean(ordination_scores$PCoA1[rows], na.rm = TRUE),
      PCoA2    = mean(ordination_scores$PCoA2[rows], na.rm = TRUE),
      r2       = envfit_fit$factors$r[trait],
      p_value  = envfit_fit$factors$pvals[trait]
    )
  }))
})) %>%
  arrange(p_value, trait, level)

write.csv(
  factor_centroid_table,
  file.path(table_dir, "bird_trait_space_envfit_factor_centroids.csv"),
  row.names = FALSE
)

# =============================================================================
# Analysis summary
# =============================================================================

analysis_summary <- tibble(
  setting = c(
    "species_n",
    "continuous_traits_n",
    "categorical_traits_n",
    "boundary_style",
    "kde_prob_levels",
    "kde_gridsize",
    "permutations",
    "include_range_size",
    "permanova_model",
    "plot_group_membership",
    "gower_asymmetric_binary_traits",
    "pcoa_correction",
    "pcoa1_variance_percent",
    "pcoa2_variance_percent"
  ),
  value = c(
    nrow(trait_data),
    length(continuous_trait_cols),
    length(binary_trait_cols),
    boundary_style,
    paste(kde_levels, collapse = ", "),
    kde_gridsize,
    permutations,
    include_range_size,
    "Gower distance ~ Threat_Logging + Threat_Plantation; marginal tests",
    "Logging and timber plantation groups are non-exclusive; species with both flags are included in both KDE surfaces.",
    paste(asymmetric_binary_cols, collapse = "; "),
    "Lingoes additive correction",
    round(axis_variance[[1]], 3),
    round(axis_variance[[2]], 3)
  )
)

write.csv(
  analysis_summary,
  file.path(table_dir, "bird_trait_space_analysis_summary.csv"),
  row.names = FALSE
)

# =============================================================================
# Plotting: build threat overlay data
# =============================================================================

cat("--- Plotting PCoA trait space ---\n")
threat_overlay_scores <- bind_rows(
  ordination_scores %>%
    filter(Threat_Logging) %>%
    mutate(threat_group = "Logging"),
  ordination_scores %>%
    filter(Threat_Plantation) %>%
    mutate(threat_group = "Timber plantations")
) %>%
  mutate(threat_group = factor(threat_group, levels = c("Logging", "Timber plantations")))

group_palette <- make_group_palette(levels(threat_overlay_scores$threat_group))
group_label   <- "Threat flag"

# =============================================================================
# Build KDE and ellipse data
# =============================================================================

kde_data     <- NULL
ellipse_data <- NULL

if (boundary_style %in% c("kde", "both")) {
  cat(sprintf(
    "  Computing KDE contours (prob levels: %s; gridsize: %d)...\n",
    paste(kde_levels, collapse = ", "),
    kde_gridsize
  ))
  kde_data <- make_kde_contours(
    threat_overlay_scores,
    prob_levels = kde_levels,
    gridsize    = kde_gridsize
  )
  
  if (!is.null(kde_data) && nrow(kde_data) > 0) {
    # Save KDE contour coordinates for reproducibility
    write.csv(
      kde_data,
      file.path(processed_dir, "bird_trait_space_kde_contours.csv"),
      row.names = FALSE
    )
  }
}

if (boundary_style %in% c("ellipse", "both")) {
  ellipse_data <- threat_overlay_scores %>%
    group_by(threat_group) %>%
    filter(n() >= 3) %>%
    ungroup()
}

# For the KDE legend we need a label for each probability level.
# Inner (lower) level = core trait space; outer = broad envelope.
kde_level_labels <- setNames(
  sprintf("%d%% density region", as.integer(kde_levels * 100)),
  as.character(kde_levels)
)

# =============================================================================
# Assemble the ordination plot
# =============================================================================

ordination_plot <- ggplot(
  ordination_scores,
  aes(x = PCoA1, y = PCoA2)
) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "gray82") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "gray82")

# --- KDE layers -------------------------------------------------------------
# Outer (95% or highest level): translucent filled polygon per group.
# Inner (50% or lowest level): solid outlined path per group.
# Both use group colour for fill/colour so the legend is self-consistent.
if (boundary_style %in% c("kde", "both") && !is.null(kde_data) && nrow(kde_data) > 0) {
  kde_data <- kde_data %>%
    mutate(threat_group = factor(threat_group, levels = levels(threat_overlay_scores$threat_group)))
  
  outer_level <- max(kde_levels)
  inner_level <- min(kde_levels)
  
  # Filled outer region
  outer_kde <- kde_data %>% filter(prob_level == outer_level)
  if (nrow(outer_kde) > 0) {
    ordination_plot <- ordination_plot +
      geom_polygon(
        data = outer_kde,
        aes(
          x     = PCoA1,
          y     = PCoA2,
          fill  = threat_group,
          group = interaction(threat_group, piece)
        ),
        alpha       = 0.12,
        color       = NA,
        inherit.aes = FALSE
      )
  }
  
  # Dashed outline for outer region
  if (nrow(outer_kde) > 0) {
    ordination_plot <- ordination_plot +
      geom_path(
        data = outer_kde,
        aes(
          x     = PCoA1,
          y     = PCoA2,
          color = threat_group,
          group = interaction(threat_group, piece)
        ),
        linewidth   = 0.4,
        linetype    = "dashed",
        inherit.aes = FALSE
      )
  }
  
  # Solid outline for inner (core) region — drawn on top
  if (length(kde_levels) >= 2) {
    inner_kde <- kde_data %>% filter(prob_level == inner_level)
    if (nrow(inner_kde) > 0) {
      ordination_plot <- ordination_plot +
        geom_path(
          data = inner_kde,
          aes(
            x     = PCoA1,
            y     = PCoA2,
            color = threat_group,
            group = interaction(threat_group, piece)
          ),
          linewidth   = 0.75,
          linetype    = "solid",
          inherit.aes = FALSE
        )
    }
  }
}

# --- t-ellipse layer (used when boundary_style is "ellipse" or "both") ------
# Uses type = "t" (multivariate t-distribution) rather than "norm" so the
# ellipse is more robust to skewed or heavy-tailed distributions.
if (boundary_style %in% c("ellipse", "both") && !is.null(ellipse_data) && nrow(ellipse_data) > 0) {
  ordination_plot <- ordination_plot +
    stat_ellipse(
      data  = ellipse_data,
      aes(x = PCoA1, y = PCoA2, group = threat_group, color = threat_group),
      type        = "t",
      level       = 0.95,
      linewidth   = 0.65,
      inherit.aes = FALSE
    )
}

# --- Species points ----------------------------------------------------------
ordination_plot <- ordination_plot +
  geom_point(color = "gray74", alpha = 0.42, size = 0.9, stroke = 0) +
  geom_point(
    data = threat_overlay_scores,
    aes(x = PCoA1, y = PCoA2, color = threat_group),
    alpha       = 0.62,
    size        = 1.15,
    stroke      = 0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = group_palette, name = group_label, drop = FALSE) +
  scale_fill_manual( values = group_palette, name = group_label, drop = FALSE) +
  coord_equal() +
  labs(
    title    = "Bird trait space by logging and timber-plantation threat",
    subtitle = paste0(
      "Gower dissimilarity; PCoA uses Lingoes correction; n = ",
      nrow(ordination_scores),
      ".\n",
      if (boundary_style %in% c("kde", "both")) {
        paste0(
          "KDE boundaries: solid = ",
          as.integer(min(kde_levels) * 100),
          "% core; dashed = ",
          as.integer(max(kde_levels) * 100),
          "% envelope. ",
          "Sheather-Jones bandwidth; species with both threat flags contribute to both surfaces."
        )
      } else {
        "95% t-distribution ellipses. Species with both threat flags contribute to both ellipses."
      }
    ),
    x = sprintf("PCoA1 (%.1f%%)", axis_variance[[1]]),
    y = sprintf("PCoA2 (%.1f%%)", axis_variance[[2]])
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0),
    plot.subtitle = element_text(color = "gray30"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(8, 12, 8, 8)
  )

ggsave(
  file.path(figure_dir, "bird_trait_space_pcoa_logging_vs_plantations.png"),
  ordination_plot,
  width  = 9.2,
  height = 7.0,
  dpi    = 300,
  bg     = "white"
)
ggsave(
  file.path(figure_dir, "bird_trait_space_pcoa_logging_vs_plantations.pdf"),
  ordination_plot,
  width  = 9.2,
  height = 7.0
)

# =============================================================================
# envfit overlay plot
# =============================================================================

envfit_plot <- ordination_plot +
  guides(fill = "none") +
  labs(
    title    = "Bird trait space with fitted trait vectors and centroids",
    subtitle = paste0(
      "Continuous traits are arrows; categorical traits are fitted centroids.\n",
      "Full envfit tables include all fits; plotted factor labels are capped for readability."
    )
  )

plot_vectors <- vector_table %>%
  filter(!is.na(p_value), p_value <= envfit_plot_p_max)

if (nrow(plot_vectors) > 0) {
  envfit_plot <- envfit_plot +
    geom_segment(
      data = plot_vectors,
      aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
      inherit.aes = FALSE,
      color       = "gray12",
      linewidth   = 0.55,
      arrow       = grid::arrow(length = grid::unit(0.17, "cm"))
    )
  
  if (has_ggrepel) {
    envfit_plot <- envfit_plot +
      ggrepel::geom_text_repel(
        data        = plot_vectors,
        aes(x = PCoA1, y = PCoA2, label = trait),
        inherit.aes = FALSE,
        size        = 3.2,
        color       = "gray12",
        min.segment.length = 0,
        box.padding = 0.28,
        max.overlaps = Inf
      )
  } else {
    envfit_plot <- envfit_plot +
      geom_text(
        data        = plot_vectors,
        aes(x = PCoA1, y = PCoA2, label = trait),
        inherit.aes = FALSE,
        size        = 3.2,
        color       = "gray12"
      )
  }
}

plot_factor_centroids <- factor_centroid_table %>%
  filter(!is.na(p_value), p_value <= envfit_plot_p_max, level == "Yes") %>%
  arrange(p_value, desc(r2), trait) %>%
  mutate(label = paste0(trait, ": ", level))

if (is.finite(envfit_max_factor_labels)) {
  plot_factor_centroids <- plot_factor_centroids %>%
    slice_head(n = as.integer(envfit_max_factor_labels))
}

if (nrow(plot_factor_centroids) > 0) {
  envfit_plot <- envfit_plot +
    geom_point(
      data = plot_factor_centroids,
      aes(x = PCoA1, y = PCoA2),
      inherit.aes = FALSE,
      shape  = 21,
      size   = 2.1,
      stroke = 0.45,
      color  = "gray10",
      fill   = "white"
    )
  
  if (has_ggrepel) {
    envfit_plot <- envfit_plot +
      ggrepel::geom_text_repel(
        data        = plot_factor_centroids,
        aes(x = PCoA1, y = PCoA2, label = label),
        inherit.aes = FALSE,
        size        = 2.9,
        color       = "gray20",
        min.segment.length = 0,
        box.padding = 0.22,
        max.overlaps = Inf
      )
  } else {
    envfit_plot <- envfit_plot +
      geom_text(
        data        = plot_factor_centroids,
        aes(x = PCoA1, y = PCoA2, label = label),
        inherit.aes = FALSE,
        size        = 2.9,
        color       = "gray20"
      )
  }
}

ggsave(
  file.path(figure_dir, "bird_trait_space_pcoa_envfit.png"),
  envfit_plot,
  width  = 10.6,
  height = 7.6,
  dpi    = 300,
  bg     = "white"
)
ggsave(
  file.path(figure_dir, "bird_trait_space_pcoa_envfit.pdf"),
  envfit_plot,
  width  = 10.6,
  height = 7.6
)

cat("\nBird trait-space PCoA complete.\n")
cat("Figures:                  ", normalizePath(figure_dir,    winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Tables:                   ", normalizePath(table_dir,     winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("Processed ordination data:", normalizePath(processed_dir, winslash = "/", mustWork = FALSE), "\n", sep = "")