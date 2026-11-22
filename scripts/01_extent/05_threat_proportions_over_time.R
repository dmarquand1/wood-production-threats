rm(list = ls())
gc()

# =============================================================================
# THREAT PROPORTIONS OVER TIME: NATURAL FOREST LOGGING (5.3) VS PLANTATIONS (2.2)
#
# Question: among species with a CURRENT (latest) IUCN assessment, does the
# proportion flagged with an ongoing plantation threat (2.2) relative to an
# ongoing logging/natural-forest threat (5.3) differ depending on how long ago
# that assessment was published? If plantation production has grown steadily
# but plantation threats are rarer among species assessed longer ago, that is
# consistent with older assessments lagging behind the shift from natural
# forest logging towards plantation-based wood production.
#
# This deliberately looks only at the cross-section of current assessments,
# grouped by their publication year, rather than tracking re-assessment
# history for individual species (the threat classification scheme has
# changed over time, which would complicate a true longitudinal comparison).
#
# Outputs:
#   outputs/extent/tables/threat_proportions_over_time.csv
#   outputs/extent/figures/threat_proportions_over_time.png
# =============================================================================

library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)
library(here)
library(scales)
library(mgcv)

table_output_dir <- here("outputs", "extent", "tables")
figure_output_dir <- here("outputs", "extent", "figures")
dir.create(table_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

iucn_levels <- c("DD", "LC", "NT", "VU", "EN", "CR")
taxa_levels <- c("Amphibia", "Aves", "Freshwater Fish", "Mammalia", "Reptilia")
summary_taxa_levels <- c("All Vertebrates", taxa_levels)

`%||%` <- function(x, y) if (is.null(x)) y else x

# Same name-cleaning rule used when the canonical assessment tables were built
# (scripts/01_extent/01_prepare_assessments.R), so joins on scientific_name match.
clean_names <- function(x) {
  if (is.null(x)) return(NA_character_)
  x %>%
    str_replace_all("\uFFFD", " ") %>%
    str_replace_all("\u00A0", " ") %>%
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") %>%
    str_squish() %>%
    str_trim()
}

# --- 1. Pull publication year for each species' current assessment ---------
# The processed/canonical tables carry threat flags but not assessment year,
# so year_published is read directly from the raw rredlist assessment caches.
taxon_config <- tibble::tribble(
  ~taxon_key,    ~taxa_group,         ~raw_file,
  "bird",        "Aves",              here("data", "raw", "iucn", "iucn_assessment_list.rds"),
  "mammal",      "Mammalia",          here("data", "raw", "iucn", "mam_iucn_assessment_list.rds"),
  "amphibian",   "Amphibia",          here("data", "raw", "iucn", "amphib_iucn_assessment_list.rds"),
  "reptile",     "Reptilia",          here("data", "raw", "iucn", "rept_iucn_assessment_list.rds"),
  "fish",        "Freshwater Fish",   here("data", "raw", "iucn", "fish_iucn_assessment_list.rds")
)

extract_year_published <- function(taxon_key, raw_file) {
  cat("Reading publication years for:", taxon_key, "\n")
  assessments <- readRDS(raw_file)
  tibble(
    taxon_key = taxon_key,
    scientific_name = clean_names(map_chr(assessments, ~ as.character(.x$taxon$scientific_name %||% NA_character_))),
    year_published = as.integer(map_chr(assessments, ~ as.character(.x$year_published %||% NA_character_)))
  )
}

year_lookup <- pmap_dfr(taxon_config, function(taxon_key, taxa_group, raw_file) {
  extract_year_published(taxon_key, raw_file)
})

# --- 2. Join publication year onto the canonical species-level threat flags -
species_summary <- readRDS(
  here("data", "processed", "assessments", "all_taxa_assessment_consolidated.rds")
) %>%
  distinct(taxon_key, scientific_name, .keep_all = TRUE) %>%
  select(
    taxon_key, taxa_group, scientific_name, category,
    has_logging_ongoing, has_plantation_ongoing, has_wood_prod_ongoing
  )

dated_species <- species_summary %>%
  left_join(year_lookup, by = c("taxon_key", "scientific_name"))

n_unmatched <- sum(is.na(dated_species$year_published))
if (n_unmatched > 0) {
  cat(
    "Note:", n_unmatched, "of", nrow(dated_species),
    "species could not be matched to a publication year and are excluded.\n"
  )
}

# --- 3. Restrict to extant species with a known publication year -----------
# Converting to the modern 6-level factor first (then dropping NAs) excludes
# both EX/EW and the deprecated pre-2001 "Lower Risk" codes (LR/cd, LR/lc,
# LR/nt), consistent with how the other extent scripts (e.g. 02b) define the
# current species denominator.
dated_species <- dated_species %>%
  mutate(
    category = factor(category, levels = iucn_levels, ordered = TRUE),
    taxa_group = factor(taxa_group, levels = taxa_levels)
  ) %>%
  filter(!is.na(category), !is.na(year_published))

# --- 4. Bin publication year into assessment periods ------------------------
# Bin edges are chosen so that every taxon has a reasonable sample in every
# period (see scratch checks: the pre-2011 bin is thin for some taxa on a
# finer grid, but is well populated at 5-year resolution).
year_breaks <- c(-Inf, 2010, 2015, 2020, Inf)
year_bin_labels <- c("2010 or earlier", "2011-2015", "2016-2020", "2021-2025")

dated_species <- dated_species %>%
  mutate(
    assessment_period = cut(
      year_published,
      breaks = year_breaks,
      labels = year_bin_labels,
      right = TRUE,
      ordered_result = TRUE
    )
  )

# --- 5. Summarise proportions by period, for each taxon and pooled ---------
# Wilson score interval: more reliable than a normal approximation when a bin
# has a small sample or a proportion near 0 (both occur here, e.g. the small
# pre-2011 mammal bin).
wilson_ci <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  phat <- x / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  margin <- (z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2))) / denom
  tibble(lower = pmax(0, center - margin), upper = pmin(1, center + margin))
}

# Empirical (Haldane-Anscombe corrected) log-odds: adding 0.5 to each count
# keeps the estimate and its CI finite even in the zero-count cells present
# here (e.g. the pre-2011 mammal bin has zero plantation-threatened species),
# rather than the -Inf a plain logit(0) would give.
empirical_log_odds <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  log_odds <- log((x + 0.5) / (n - x + 0.5))
  se <- sqrt(1 / (x + 0.5) + 1 / (n - x + 0.5))
  tibble(log_odds = log_odds, lower = log_odds - z * se, upper = log_odds + z * se)
}

summarise_period <- function(df) {
  df %>%
    group_by(assessment_period) %>%
    summarise(
      n_species = n(),
      n_logging = sum(has_logging_ongoing, na.rm = TRUE),
      n_plantation = sum(has_plantation_ongoing, na.rm = TRUE),
      n_wood_prod = sum(has_wood_prod_ongoing, na.rm = TRUE),
      .groups = "drop"
    )
}

period_summary <- bind_rows(
  dated_species %>%
    group_by(taxa_group) %>%
    group_modify(~ summarise_period(.x)) %>%
    ungroup() %>%
    rename(summary_group = taxa_group) %>%
    mutate(summary_group = as.character(summary_group)),
  summarise_period(dated_species) %>%
    mutate(summary_group = "All Vertebrates", .before = 1)
) %>%
  mutate(
    prop_logging = n_logging / n_species,
    prop_plantation = n_plantation / n_species,
    prop_wood_prod = n_wood_prod / n_species,
    # Among species with a wood-production threat, the share specifically
    # attributable to plantations vs natural-forest logging (they are not
    # mutually exclusive, so these need not sum to 1).
    plantation_share_of_wood_prod = if_else(n_wood_prod > 0, n_plantation / n_wood_prod, NA_real_),
    logging_share_of_wood_prod = if_else(n_wood_prod > 0, n_logging / n_wood_prod, NA_real_),
    summary_group = factor(summary_group, levels = summary_taxa_levels)
  ) %>%
  bind_cols(
    wilson_ci(.$n_logging, .$n_species) %>% rename(logging_ci_lower = lower, logging_ci_upper = upper)
  ) %>%
  bind_cols(
    wilson_ci(.$n_plantation, .$n_species) %>% rename(plantation_ci_lower = lower, plantation_ci_upper = upper)
  ) %>%
  bind_cols(
    empirical_log_odds(.$n_logging, .$n_species) %>%
      rename(logging_log_odds = log_odds, logging_log_odds_lower = lower, logging_log_odds_upper = upper)
  ) %>%
  bind_cols(
    empirical_log_odds(.$n_plantation, .$n_species) %>%
      rename(plantation_log_odds = log_odds, plantation_log_odds_lower = lower, plantation_log_odds_upper = upper)
  ) %>%
  arrange(summary_group, assessment_period)

write.csv(
  period_summary %>% mutate(summary_group = as.character(summary_group)),
  file.path(table_output_dir, "threat_proportions_over_time.csv"),
  row.names = FALSE
)

cat("\nProportion of species with an ongoing logging (5.3) or plantation (2.2)\n")
cat("threat, by current-assessment publication period:\n")
print(as.data.frame(period_summary), row.names = FALSE)

# --- 6. Plot: proportion of species threatened, by assessment period -------
plot_data <- period_summary %>%
  filter(!is.na(assessment_period)) %>%
  select(
    summary_group, assessment_period, n_species,
    prop_logging, logging_ci_lower, logging_ci_upper,
    prop_plantation, plantation_ci_lower, plantation_ci_upper
  ) %>%
  pivot_longer(
    cols = c(prop_logging, prop_plantation),
    names_to = "threat_type",
    values_to = "proportion"
  ) %>%
  mutate(
    ci_lower = if_else(threat_type == "prop_logging", logging_ci_lower, plantation_ci_lower),
    ci_upper = if_else(threat_type == "prop_logging", logging_ci_upper, plantation_ci_upper),
    threat_type = recode(
      threat_type,
      prop_logging = "Natural Forestry",
      prop_plantation = "Plantation"
    ),
    threat_type = factor(threat_type, levels = c("Natural Forestry", "Plantation"))
  ) %>%
  select(-logging_ci_lower, -logging_ci_upper, -plantation_ci_lower, -plantation_ci_upper)

label_data <- plot_data %>%
  group_by(summary_group, threat_type) %>%
  filter(assessment_period == max(assessment_period)) %>%
  ungroup()

threat_colors <- c(
  "Natural Forestry" = "#eda100",
  "Plantation" = "#2a78d6"
)

threat_over_time_plot <- ggplot(
  plot_data,
  aes(x = assessment_period, y = proportion, color = threat_type, group = threat_type)
) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.5,
    position = position_dodge(width = 0.15)
  ) +
  geom_line(linewidth = 0.9, position = position_dodge(width = 0.15)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.15)) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label = percent(proportion, accuracy = 0.1)),
    size = 3,
    show.legend = FALSE,
    direction = "y",
    nudge_x = 0.35,
    hjust = 0,
    segment.size = 0.3,
    min.segment.length = 0
  ) +
  facet_wrap(~ summary_group, nrow = 2, scales = "free_y") +
  scale_color_manual(name = "Production Pathway", values = threat_colors) +
  scale_y_continuous(
    name = "Proportion of assessed species",
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  scale_x_discrete(name = "Publication year of current assessment") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )

print(threat_over_time_plot)

ggsave(
  filename = file.path(figure_output_dir, "threat_proportions_over_time.png"),
  plot = threat_over_time_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("\nTable saved to:", file.path(table_output_dir, "threat_proportions_over_time.csv"), "\n")
cat("Plot saved to:", file.path(figure_output_dir, "threat_proportions_over_time.png"), "\n")

# --- 6b. Plot: log-odds of the same threats, by assessment period ----------
# Same layout as the proportion plot above, but on the log-odds (logit)
# scale, which is the natural scale for a binomial trend (c.f. the binomial
# GLM discussion). Counts use the Haldane-Anscombe correction (+0.5) so that
# zero-count cells (e.g. plantation threats in the earliest mammal bin) give
# a finite log-odds and CI rather than -Inf.
log_odds_plot_data <- period_summary %>%
  filter(!is.na(assessment_period)) %>%
  select(
    summary_group, assessment_period, n_species,
    logging_log_odds, logging_log_odds_lower, logging_log_odds_upper,
    plantation_log_odds, plantation_log_odds_lower, plantation_log_odds_upper
  ) %>%
  pivot_longer(
    cols = c(logging_log_odds, plantation_log_odds),
    names_to = "threat_type",
    values_to = "log_odds"
  ) %>%
  mutate(
    ci_lower = if_else(
      threat_type == "logging_log_odds", logging_log_odds_lower, plantation_log_odds_lower
    ),
    ci_upper = if_else(
      threat_type == "logging_log_odds", logging_log_odds_upper, plantation_log_odds_upper
    ),
    threat_type = recode(
      threat_type,
      logging_log_odds = "Natural Forestry",
      plantation_log_odds = "Plantation"
    ),
    threat_type = factor(threat_type, levels = c("Natural Forestry", "Plantation"))
  ) %>%
  select(-logging_log_odds_lower, -logging_log_odds_upper, -plantation_log_odds_lower, -plantation_log_odds_upper)

log_odds_label_data <- log_odds_plot_data %>%
  group_by(summary_group, threat_type) %>%
  filter(assessment_period == max(assessment_period)) %>%
  ungroup()

log_odds_over_time_plot <- ggplot(
  log_odds_plot_data,
  aes(x = assessment_period, y = log_odds, color = threat_type, group = threat_type)
) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4, linetype = "dotted") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    linewidth = 0.5,
    position = position_dodge(width = 0.15)
  ) +
  geom_line(linewidth = 0.9, position = position_dodge(width = 0.15)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.15)) +
  ggrepel::geom_text_repel(
    data = log_odds_label_data,
    aes(label = sprintf("%.2f", log_odds)),
    size = 3,
    show.legend = FALSE,
    direction = "y",
    nudge_x = 0.35,
    hjust = 0,
    segment.size = 0.3,
    min.segment.length = 0
  ) +
  facet_wrap(~ summary_group, nrow = 2, scales = "free_y") +
  scale_color_manual(name = "Production Pathway", values = threat_colors) +
  scale_y_continuous(
    name = "Log odds of threat, ln(p / (1-p))",
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  scale_x_discrete(name = "Publication year of current assessment") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 10, face = "bold")
  )

print(log_odds_over_time_plot)

ggsave(
  filename = file.path(figure_output_dir, "threat_log_odds_over_time.png"),
  plot = log_odds_over_time_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("Plot saved to:", file.path(figure_output_dir, "threat_log_odds_over_time.png"), "\n")

# --- 6c. Plot: smoothed log-odds trend on the raw (unbinned) year ----------
# Alternative to the period-binned log-odds plot above: instead of four
# discrete periods, fit a binomial GAM smooth (logit link) directly against
# each species' assessment publication year, per taxon. This uses the full
# year-level resolution without the small-sample jitter a literal
# one-point-per-year version would show (many taxon-year cells have single-
# digit counts). predict(..., type = "link") gives the fitted curve and its
# SE directly on the log-odds scale, so no continuity correction is needed
# here (unlike the discrete-count empirical logit used for the binned plot).
fit_smooth_log_odds <- function(df, response_var, k_max = 8) {
  df$response <- as.numeric(df[[response_var]])
  n_years <- n_distinct(df$year_published)
  k <- max(3, min(k_max, n_years - 1))
  model <- gam(
    response ~ s(year_published, k = k),
    family = binomial(link = "logit"),
    method = "REML",
    data = df
  )
  year_seq <- seq(min(df$year_published), max(df$year_published), by = 1)
  pred <- predict(model, newdata = data.frame(year_published = year_seq), type = "link", se.fit = TRUE)
  tibble(
    year_published = year_seq,
    log_odds = pred$fit,
    ci_lower = pred$fit - 1.96 * pred$se.fit,
    ci_upper = pred$fit + 1.96 * pred$se.fit
  )
}

smooth_log_odds_data <- map_dfr(summary_taxa_levels, function(g) {
  group_data <- if (g == "All Vertebrates") dated_species else dated_species %>% filter(taxa_group == g)
  if (nrow(group_data) == 0) return(NULL)
  bind_rows(
    fit_smooth_log_odds(group_data, "has_logging_ongoing") %>% mutate(threat_type = "Natural Forestry"),
    fit_smooth_log_odds(group_data, "has_plantation_ongoing") %>% mutate(threat_type = "Plantation")
  ) %>%
    mutate(summary_group = g)
}) %>%
  mutate(
    summary_group = factor(summary_group, levels = summary_taxa_levels),
    threat_type = factor(threat_type, levels = c("Natural Forestry", "Plantation"))
  )

# Rug marks show where assessments actually fall along the year axis, so the
# curve isn't read as equally well-supported in sparsely-sampled years.
rug_data <- dated_species %>%
  transmute(summary_group = as.character(taxa_group), year_published) %>%
  bind_rows(dated_species %>% transmute(summary_group = "All Vertebrates", year_published)) %>%
  mutate(summary_group = factor(summary_group, levels = summary_taxa_levels))

smooth_log_odds_plot <- ggplot(
  smooth_log_odds_data,
  aes(x = year_published, y = log_odds, color = threat_type, fill = threat_type)
) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4, linetype = "dotted") +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  geom_rug(
    data = rug_data,
    aes(x = year_published),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.15,
    length = unit(0.02, "npc")
  ) +
  facet_wrap(~ summary_group, nrow = 2, scales = "free_y") +
  scale_color_manual(name = "Production Pathway", values = threat_colors) +
  scale_fill_manual(name = "Production Pathway", values = threat_colors) +
  scale_x_continuous(name = "Publication year of current assessment") +
  scale_y_continuous(
    name = "Log odds of threat, ln(p / (1-p))",
    expand = expansion(mult = c(0.08, 0.1))
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold")
  )

print(smooth_log_odds_plot)

ggsave(
  filename = file.path(figure_output_dir, "threat_log_odds_smooth_over_time.png"),
  plot = smooth_log_odds_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("Plot saved to:", file.path(figure_output_dir, "threat_log_odds_smooth_over_time.png"), "\n")

# =============================================================================
# --- 7. ASSESSMENT AGE: HOW OLD IS EACH TAXON'S CURRENT ASSESSMENT? ---------
# =============================================================================
# Age = years between the current assessment's publication and today. This is
# the distribution underlying the period bins above: it shows directly how
# much of each taxon's current Red List status rests on an old assessment.

current_year <- as.integer(format(Sys.Date(), "%Y"))

age_data <- dated_species %>%
  mutate(age_years = current_year - year_published)

age_summary <- age_data %>%
  group_by(taxa_group) %>%
  summarise(
    n_species = n(),
    min_age = min(age_years),
    median_age = median(age_years),
    mean_age = mean(age_years),
    max_age = max(age_years),
    .groups = "drop"
  )

write.csv(
  age_summary,
  file.path(table_output_dir, "assessment_age_summary.csv"),
  row.names = FALSE
)

cat("\nAssessment age (years since publication), by taxon:\n")
print(as.data.frame(age_summary), row.names = FALSE)

median_age_data <- age_summary %>%
  mutate(label = paste0("median = ", median_age, "y"))

assessment_age_plot <- ggplot(age_data, aes(x = age_years)) +
  geom_histogram(binwidth = 1, boundary = 0, fill = "#2a78d6", color = "white") +
  geom_vline(
    data = median_age_data,
    aes(xintercept = median_age),
    linetype = "dashed", color = "#0b0b0b", linewidth = 0.5
  ) +
  geom_text(
    data = median_age_data,
    aes(x = median_age, y = Inf, label = label),
    hjust = -0.1, vjust = 1.5, size = 3, color = "#0b0b0b"
  ) +
  facet_wrap(~ taxa_group, nrow = 2, scales = "free_y") +
  scale_x_continuous(
    name = "Assessment age (years since publication)",
    breaks = seq(0, 30, by = 5)
  ) +
  scale_y_continuous(name = "Number of species", labels = comma) +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"))

print(assessment_age_plot)

ggsave(
  filename = file.path(figure_output_dir, "assessment_age_histogram.png"),
  plot = assessment_age_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)

cat("\nTable saved to:", file.path(table_output_dir, "assessment_age_summary.csv"), "\n")
cat("Plot saved to:", file.path(figure_output_dir, "assessment_age_histogram.png"), "\n")
