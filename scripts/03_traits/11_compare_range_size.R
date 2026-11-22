rm(list = ls())
gc()
# =============================================================================
# SCRIPT: CALCULATE RAW RANGE SIZES FOR ALL TAXA
#
# Description:
#   Loops over all five taxa, computes true (unscaled) range size in km2 from
#   the dissolved range polygons using a Mollweide equal-area projection, saves
#   the raw per-species km2 for each taxon, and builds a cross-taxon summary so
#   the taxa can be ranked by average range size.
#
#   The large spatial object for each taxon is dropped and memory is cleared
#   before the next taxon is loaded. A progress bar with elapsed time and ETA
#   is shown while each taxon is processed.
# =============================================================================

# --- 1. SETUP ----------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  dplyr,      # data manipulation
  sf,         # spatial vector data
  units,      # measurement units
  here,       # project-relative paths
  progress,   # progress bar with ETA
  ggplot2,    # plotting
  tidyr,      # reshaping (pivot_longer)
  scales      # axis label formatting
)

# --- 2. DEFINE PROJECTION ----------------------------------------------------
# Mollweide equal-area projection: standard Lat/Lon (WGS84) distorts area, so
# we reproject before measuring so areas are comparable across latitudes.
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# --- 3. CONFIG ---------------------------------------------------------------
taxa <- c("amphibian", "bird", "fish", "mammal", "reptile")

taxon_labels <- c(
  amphibian = "Amphibians",
  bird      = "Birds",
  fish      = "Fish",
  mammal    = "Mammals",
  reptile   = "Reptiles"
)

range_dir <- here("data", "processed", "ranges")
dir.create(range_dir, recursive = TRUE, showWarnings = FALSE)

# --- 4. PER-TAXON WORKER -----------------------------------------------------
# Heavy spatial objects live ONLY inside this function. When it returns, they
# go out of scope and become eligible for garbage collection; they are also
# explicitly removed + gc()'d before the (small) result is returned.
compute_taxon_ranges <- function(taxon) {
  label <- taxon_labels[[taxon]]
  path <- here("data", "processed", "ranges",
               sprintf("all_%s_ranges_dissolved.rds", taxon))
  
  if (!file.exists(path)) {
    warning("Missing range file for ", label, ": ", path, call. = FALSE)
    return(NULL)
  }
  
  message("\n=== ", label, " : loading ", basename(path), " ===")
  ranges <- readRDS(path)
  
  if (!"sci_name" %in% names(ranges)) {
    stop("No 'sci_name' column in ", path, call. = FALSE)
  }
  
  # One O(n) pass to index the rows belonging to each species. Subsetting by
  # integer index afterwards is far cheaper than filtering the whole sf object
  # once per species.
  row_index <- split(seq_len(nrow(ranges)), ranges$sci_name)
  species_list <- names(row_index)
  n <- length(species_list)
  
  pb <- progress::progress_bar$new(
    format = paste0(
      "  ", label,
      " [:bar] :current/:total (:percent) | elapsed: :elapsedfull | eta: :eta"
    ),
    total = n,
    clear = FALSE,
    width = 80
  )
  
  areas_km2 <- numeric(n)
  
  for (i in seq_len(n)) {
    idx <- row_index[[i]]
    
    areas_km2[i] <- tryCatch({
      g <- st_transform(ranges[idx, ], crs = moll_crs)
      # sum() handles species split across multiple polygons / rows
      as.numeric(set_units(sum(st_area(g)), km^2))
    }, error = function(e) NA_real_)
    
    pb$tick()
  }
  
  result <- data.frame(
    taxon = taxon,
    Taxon = label,
    species_name = species_list,
    range_size_km2 = areas_km2,
    stringsAsFactors = FALSE
  )
  
  # Drop the heavy objects and clear memory before returning the small result.
  rm(ranges, row_index)
  gc()
  
  result
}

# --- 5. LOOP OVER TAXA -------------------------------------------------------
overall_start <- Sys.time()
summary_list <- list()

for (taxon in taxa) {
  
  taxon_result <- compute_taxon_ranges(taxon)
  
  if (!is.null(taxon_result)) {
    
    # 5a. Save the raw, unscaled per-species range sizes for this taxon.
    out_csv <- file.path(range_dir, sprintf("%s_range_sizes_km2.csv", taxon))
    write.csv(taxon_result, out_csv, row.names = FALSE)
    message("  saved raw range sizes -> ", out_csv)
    
    # 5b. Summarise. Range-size distributions are heavily right-skewed, so the
    #     median and geometric mean are more representative of a "typical"
    #     range than the arithmetic mean (which a few huge ranges dominate).
    valid <- taxon_result$range_size_km2
    valid <- valid[is.finite(valid) & valid > 0]
    
    summary_list[[taxon]] <- data.frame(
      Taxon         = taxon_labels[[taxon]],
      n_species     = length(valid),
      mean_km2      = mean(valid),
      median_km2    = median(valid),
      geom_mean_km2 = 10^(mean(log10(valid))),
      min_km2       = min(valid),
      max_km2       = max(valid),
      stringsAsFactors = FALSE
    )
  }
  
  # Remove the per-taxon result and clear memory before the next taxon loads.
  rm(taxon_result)
  gc()
}

# --- 6. CROSS-TAXON SUMMARY --------------------------------------------------
range_summary <- do.call(rbind, summary_list) %>%
  arrange(desc(median_km2))

summary_csv <- file.path(range_dir, "range_size_summary_by_taxon.csv")
write.csv(range_summary, summary_csv, row.names = FALSE)

total_elapsed <- round(difftime(Sys.time(), overall_start, units = "mins"), 1)
message("\nTotal elapsed: ", total_elapsed, " mins")
message("Saved per-taxon summary -> ", summary_csv)
print(range_summary)

# =============================================================================
# 7. FIGURE: RAW RANGE SIZE BY THREAT STATUS, ACROSS TAXA
#
#   Joins the raw per-species range sizes saved above to the logging and
#   plantation threat flags, then compares the range size of threatened vs not
#   threatened species within each taxon, for each pressure. Range sizes are
#   raw equal-area km2 (directly comparable across taxa); the y-axis is log10
#   because range-size distributions are strongly right-skewed.
# =============================================================================

figure_dir <- here("outputs", "ranges", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# --- helpers (threat flags may be stored as logical / numeric / text) --------
as_logical_flag <- function(x) {
  if (is.logical(x)) { x[is.na(x)] <- FALSE; return(x) }
  if (is.numeric(x)) { y <- x != 0; y[is.na(y)] <- FALSE; return(y) }
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

pick_column <- function(df, candidates, path) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop("Missing required column in ", path, ": ",
         paste(candidates, collapse = " / "), call. = FALSE)
  }
  hit[[1]]
}

# --- read logging / plantation threat flags for one taxon --------------------
# Threat data is keyed by the same IUCN/spatial name used for the ranges, so it
# joins to the raw range table on species_name.
read_taxon_threats <- function(taxon) {
  path <- here("data", "processed", "traits", taxon,
               sprintf("%s_trait_model_final_data.csv", taxon))
  if (!file.exists(path)) {
    warning("Missing threat/model data for ", taxon, ": ", path, call. = FALSE)
    return(NULL)
  }
  
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  name_col  <- pick_column(df, c("spatial_name", "IUCN_Name"), path)
  log_col   <- pick_column(df, c("Threat_Logging", "Threat_Logging.x"), path)
  plant_col <- pick_column(df, c("Threat_Plantation", "Threat_Plantation.x"), path)
  
  out <- data.frame(
    species_name      = df[[name_col]],
    Threat_Logging    = as_logical_flag(df[[log_col]]),
    Threat_Plantation = as_logical_flag(df[[plant_col]]),
    stringsAsFactors  = FALSE
  )
  out[!duplicated(out$species_name), , drop = FALSE]
}

# --- assemble raw range size + threat status across taxa ---------------------
threat_range_data <- lapply(taxa, function(tx) {
  rs_path <- file.path(range_dir, sprintf("%s_range_sizes_km2.csv", tx))
  if (!file.exists(rs_path)) {
    warning("Missing raw range CSV for ", tx, ": ", rs_path, call. = FALSE)
    return(NULL)
  }
  
  rs <- read.csv(rs_path, stringsAsFactors = FALSE)
  th <- read_taxon_threats(tx)
  if (is.null(th)) return(NULL)
  
  # inner_join drops range-only species with no threat record
  inner_join(rs, th, by = "species_name")
}) %>%
  bind_rows()

# --- long format: one row per species per threat type ------------------------
plot_long <- threat_range_data %>%
  filter(is.finite(range_size_km2), range_size_km2 > 0) %>%
  pivot_longer(
    cols = c(Threat_Logging, Threat_Plantation),
    names_to = "threat_type",
    values_to = "threatened"
  ) %>%
  mutate(
    threat_type = recode(threat_type,
                         Threat_Logging    = "Logging threat",
                         Threat_Plantation = "Plantation threat"),
    threat_type = factor(threat_type,
                         levels = c("Logging threat", "Plantation threat")),
    status = factor(ifelse(threatened, "Threatened", "Not threatened"),
                    levels = c("Threatened", "Not threatened")),
    Taxon = factor(Taxon, levels = unname(taxon_labels))
  )

status_colors <- c("Threatened" = "#c0392b", "Not threatened" = "#9aa0a6")

threat_range_plot <- ggplot(
  plot_long,
  aes(x = Taxon, y = range_size_km2, fill = status)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7,
    outlier.alpha = 0.18,
    outlier.size = 0.7,
    linewidth = 0.4,
    colour = "gray25"
  ) +
  facet_wrap(~ threat_type, ncol = 1) +
  scale_y_log10(labels = scales::label_comma()) +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = status_colors, name = NULL) +
  labs(
    title = "Raw range size by threat status, across taxa",
    x = NULL,
    y = expression("Range size (km"^2 * ", log"[10] * " scale)"),
    caption = paste(
      "Within each taxon, species threatened vs not threatened by the pressure are compared.",
      "Equal-area range sizes in km2, directly comparable across taxa. Unknown status treated as not threatened.",
      sep = "\n"
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    legend.position = "top",
    strip.background = element_rect(fill = "gray92", colour = "gray45"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9, colour = "gray35"),
    plot.margin = margin(8, 12, 8, 8)
  )

threat_fig_png <- file.path(figure_dir, "all_taxa_raw_range_size_by_threat.png")
threat_fig_pdf <- file.path(figure_dir, "all_taxa_raw_range_size_by_threat.pdf")

ggsave(threat_fig_png, threat_range_plot, width = 11, height = 8, dpi = 300)
ggsave(threat_fig_pdf, threat_range_plot, width = 11, height = 8)

message("Saved raw range size vs threat figure -> ", threat_fig_png)
message("Saved raw range size vs threat figure -> ", threat_fig_pdf)

# =============================================================================
# 8. RED LIST CATEGORY COMPOSITION AND RANGE SIZE BY CATEGORY
#
#   For each taxon: the proportion of species in each IUCN Red List category
#   (DD, LC, NT, VU, EN, CR, EW, EX, ...) and the mean and median raw range
#   size (km2) of each category group.
#
#   Categories are taken from the processed IUCN assessment summary
#   (data/processed/assessments/all_taxa_assessment_species_summary.rds), the
#   same source used by 02_summarise_threat_assessment_completeness.R. It holds
#   taxon_key, scientific_name and the Red List category code for every assessed
#   species. That assessed set is the denominator for the proportions; range
#   means/medians use the species in each category that have a computed
#   (finite, positive) range from section 5.
# =============================================================================

table_dir <- here("outputs", "ranges", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

assessment_summary_path <- here("data", "processed", "assessments",
                                "all_taxa_assessment_species_summary.rds")

# Standard Red List ordering (most to least imperilled, then DD / NE).
std_levels <- c("EX", "EW", "CR", "EN", "VU", "NT", "LR/cd", "LC", "DD", "NE")

# Map codes / full names (any case) onto standard abbreviations.
normalise_category <- function(x) {
  raw <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    raw %in% c("ex", "extinct") ~ "EX",
    raw %in% c("ew", "extinct in the wild") ~ "EW",
    raw %in% c("cr", "critically endangered") ~ "CR",
    raw %in% c("en", "endangered") ~ "EN",
    raw %in% c("vu", "vulnerable") ~ "VU",
    raw %in% c("nt", "near threatened", "lr/nt",
               "lower risk/near threatened") ~ "NT",
    raw %in% c("lc", "least concern", "lr/lc",
               "lower risk/least concern") ~ "LC",
    raw %in% c("lr/cd", "cd", "conservation dependent",
               "lower risk/conservation dependent") ~ "LR/cd",
    raw %in% c("dd", "data deficient") ~ "DD",
    raw %in% c("ne", "not evaluated") ~ "NE",
    is.na(raw) | raw %in% c("", "na") ~ NA_character_,
    TRUE ~ toupper(raw)   # keep any unrecognised code, uppercased
  )
}

# Light name cleaning so assessment names line up with the range species names.
clean_name <- function(x) gsub("\\s+", " ", trimws(as.character(x)))

if (!file.exists(assessment_summary_path)) {
  
  message("Missing assessment summary: ", assessment_summary_path,
          " - section 8 skipped.")
  
} else {
  
  # taxon_key, scientific_name, category for every assessed species.
  category_lookup <- readRDS(assessment_summary_path) %>%
    transmute(
      taxon_key    = as.character(taxon_key),
      species_name = clean_name(scientific_name),
      category     = normalise_category(category)
    ) %>%
    filter(!is.na(species_name), species_name != "") %>%
    distinct(taxon_key, species_name, .keep_all = TRUE)
  
  # Join the raw km2 range sizes onto the assessed species, per taxon.
  category_range_data <- lapply(taxa, function(tx) {
    cats <- category_lookup %>% filter(taxon_key == tx)
    if (nrow(cats) == 0) return(NULL)
    
    rs_path <- file.path(range_dir, sprintf("%s_range_sizes_km2.csv", tx))
    if (file.exists(rs_path)) {
      rs <- read.csv(rs_path, stringsAsFactors = FALSE)
      rs$species_name <- clean_name(rs$species_name)
      cats <- left_join(cats, rs[, c("species_name", "range_size_km2")],
                        by = "species_name")
    } else {
      cats$range_size_km2 <- NA_real_
    }
    
    cats$Taxon <- taxon_labels[[tx]]
    cats
  }) %>%
    bind_rows()
  
  if (!is.null(category_range_data) && nrow(category_range_data) > 0) {
    
    cat_clean <- category_range_data %>%
      filter(!is.na(category)) %>%
      mutate(
        valid_range = ifelse(is.finite(range_size_km2) & range_size_km2 > 0,
                             range_size_km2, NA_real_),
        Taxon = factor(Taxon, levels = unname(taxon_labels))
      )
    
    taxon_totals  <- cat_clean %>% count(Taxon, name = "n_taxon")
    present_extra <- setdiff(unique(cat_clean$category), c(std_levels, NA))
    
    category_summary <- cat_clean %>%
      group_by(Taxon, category) %>%
      summarise(
        n_species        = n(),
        n_with_range     = sum(!is.na(valid_range)),
        mean_range_km2   = mean(valid_range, na.rm = TRUE),
        median_range_km2 = median(valid_range, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(taxon_totals, by = "Taxon") %>%
      mutate(
        prop_species = n_species / n_taxon,
        category = factor(category, levels = c(std_levels, present_extra))
      ) %>%
      arrange(Taxon, category) %>%
      select(Taxon, category, n_species, prop_species,
             n_with_range, mean_range_km2, median_range_km2)
    
    category_csv <- file.path(table_dir, "range_size_by_redlist_category.csv")
    write.csv(category_summary, category_csv, row.names = FALSE)
    
    message("Saved Red List category summary -> ", category_csv)
    print(as.data.frame(category_summary))
    
    # --- Figure: log range size by Red List category, per taxon -------------
    plot_cats <- cat_clean %>%
      filter(!is.na(valid_range)) %>%
      mutate(category = factor(category, levels = c(std_levels, present_extra)))
    
    # Approximate official IUCN Red List category colours.
    category_palette <- c(
      "EX"    = "#000000",
      "EW"    = "#5A1846",
      "CR"    = "#D81E05",
      "EN"    = "#FC7F3F",
      "VU"    = "#F9E814",
      "NT"    = "#CCE226",
      "LR/cd" = "#CCE226",
      "LC"    = "#60C659",
      "DD"    = "#D1D1C6",
      "NE"    = "#BCBCBC"
    )
    
    category_range_plot <- ggplot(
      plot_cats,
      aes(x = category, y = valid_range, fill = category)
    ) +
      geom_boxplot(
        width = 0.7,
        outlier.alpha = 0.15,
        outlier.size = 0.6,
        linewidth = 0.4,
        colour = "gray25",
        show.legend = FALSE
      ) +
      facet_wrap(~ Taxon, nrow = 2, scales = "free_x") +
      scale_y_log10(labels = scales::label_comma()) +
      annotation_logticks(sides = "l") +
      scale_fill_manual(values = category_palette, na.value = "grey70",
                        drop = TRUE, guide = "none") +
      labs(
        title = "Range size by Red List category, across taxa",
        x = "IUCN Red List category",
        y = expression("Range size (km"^2 * ", log"[10] * " scale)"),
        caption = "Raw equal-area range sizes (km2). Categories ordered most to least imperilled; only categories present in each taxon are shown."
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0),
        strip.background = element_rect(fill = "gray92", colour = "gray45"),
        strip.text = element_text(face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.caption = element_text(hjust = 0, size = 9, colour = "gray35"),
        plot.margin = margin(8, 12, 8, 8)
      )
    
    cat_fig_png <- file.path(figure_dir, "all_taxa_range_size_by_redlist_category.png")
    cat_fig_pdf <- file.path(figure_dir, "all_taxa_range_size_by_redlist_category.pdf")
    
    ggsave(cat_fig_png, category_range_plot, width = 12, height = 7, dpi = 300)
    ggsave(cat_fig_pdf, category_range_plot, width = 12, height = 7)
    
    message("Saved range size by category figure -> ", cat_fig_png)
    message("Saved range size by category figure -> ", cat_fig_pdf)
    
  } else {
    message("No Red List category data available; section 8 skipped.")
  }
}

