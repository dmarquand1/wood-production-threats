# =============================================================================
#
# IUCN RED LIST ANALYSIS: ALL VERTEBRATE & FRESHWATER FISH THREATS
# FROM ONGOING WOOD PRODUCTION (PLOTS A, B, AND C)
#
# Date: 12 November 2025 (Updated)
#
# Description:
# This script fetches and combines threat and assessment data for five
# major vertebrate groups (Aves, Mammalia, Amphibia, Reptilia) and
# Freshwater Fish (Actinopterygii).

#
# It generates a three-part composite figure:
#   - Plot A: Horizontal stacked bars showing overall threat proportions
#             and a breakdown of wood production threat types.
#   - Plot B: A row of five faceted donut charts, one for each taxa,
#             visualizing the same data as Plot A.
#   - Plot C: A stacked radial bar chart breaking down threats by
#             taxa and IUCN Red List category.
#
# =============================================================================


# --- 1. SETUP: Load Libraries ---
# -----------------------------------------------------------------------------
# Note: The 'patchwork' and 'ggnewscale' packages are required.
# Please ensure they are installed via install.packages() if needed.

library(dplyr)       # For data manipulation (filter, group_by, etc.)
library(tidyr)       # For data tidying
library(purrr)       # For working with lists (used for API results)
library(ggplot2)     # For creating plots
library(ggnewscale)  # For multiple fill/color scales in one plot
library(here)        # For robust file path management
library(scales)      # For plot labels (e.g., percent_format)
library(patchwork)   # For combining multiple ggplot objects

figure_output_dir <- here::here("outputs", "extent", "figures")
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

load_assessment_summary <- function(taxon_key) {
  readRDS(
    here::here(
      "data",
      "processed",
      "assessments",
      paste0(taxon_key, "_assessment_species_summary.rds")
    )
  )
}

# --- 2. LOAD PREPARED ASSESSMENT SUMMARIES ---
# -----------------------------------------------------------------------------
cat("Loading pre-processed assessment summaries...\n")

all_assessment_df_consolidated <- bind_rows(
  load_assessment_summary("bird"),
  load_assessment_summary("mammal"),
  load_assessment_summary("amphibian"),
  load_assessment_summary("reptile"),
  load_assessment_summary("fish")
)

# Define standard factor levels for consistency across plots
iucn_levels <- c("DD", "LC", "NT", "VU", "EN", "CR")
taxa_levels <- c("Amphibia", "Aves", "Freshwater Fish", "Mammalia", "Reptilia")

all_assessment_df_consolidated <- all_assessment_df_consolidated %>%
  mutate(
    category = factor(category, levels = iucn_levels, ordered = TRUE),
    taxa_group = factor(taxa_group, levels = taxa_levels)
  ) %>%
  filter(!is.na(category))


# --- 8. COMBINED DATA ANALYSIS: ONGOING WOOD PRODUCTION THREATS ---
# -----------------------------------------------------------------------------
# Filter consolidated data to only extant species (remove EX, EW)
extant_species_df <- all_assessment_df_consolidated %>%
  filter(!category %in% c("EX", "EW"))

# Define the 4-category classification used in all plots
all_threat_labels_c_ordered <- c(
  "Wood production threat & extinction risk",    # Dark red  (VU/EN/CR + wood prod)
  "Wood production threat & no extinction risk", # Pink      (LC/NT + wood prod)
  "Wood production threat & data deficient",     # Dark grey (DD + wood prod)
  "No wood production threat"                    # Light blue
)

# Classify all extant species into one of the four groups
# This is the primary dataset used for all three plots
all_data_classified <- extant_species_df %>%
  mutate(
    has_wood_harvesting_threat = has_wood_prod_ongoing,
    is_threatened_iucn = category %in% c("VU", "EN", "CR")
  ) %>%
  mutate(threat_type_stacked = case_when(
    has_wood_harvesting_threat & is_threatened_iucn ~ all_threat_labels_c_ordered[1], # Yellow
    has_wood_harvesting_threat & !is_threatened_iucn & category %in% c("LC", "NT") ~ all_threat_labels_c_ordered[2], # Dark Blue
    has_wood_harvesting_threat & category == "DD" ~ all_threat_labels_c_ordered[3], # Black
    !has_wood_harvesting_threat ~ all_threat_labels_c_ordered[4] # Light Blue
  )) %>%
  filter(!is.na(threat_type_stacked))


# =============================================================================
# --- 9. DATA PREP FOR PLOT A (Horizontal Bars) ---
# =============================================================================

cat("Preparing data for Plot A (Horizontal Bars)...\n")

# --- 9a. Data for Top Bar (Overall Breakdown) ---
# Define the stacking order for the top bar (right-to-left)
plot_a1_factor_levels <- c(
  "No wood production threat",
  "Wood production threat & no extinction risk",
  "Wood production threat & extinction risk",
  "Wood production threat & data deficient"
)

# Summarize counts and proportions for the 4 categories
plot_a1_data <- all_data_classified %>%
  count(threat_type_stacked, name = "count") %>%
  mutate(
    proportion = count / sum(count),
    threat_type_stacked = factor(threat_type_stacked, levels = plot_a1_factor_levels),
    plot_group = "A" # A dummy variable for the x-axis
  )

# Calculate totals for the main label
total_threatened_count <- sum(plot_a1_data[plot_a1_data$threat_type_stacked != "Species without Wood Production Threat", "count"], na.rm = TRUE)
total_species_count <- sum(plot_a1_data$count)
total_threatened_prop <- total_threatened_count / total_species_count

plot_a1_label <- paste0(
  scales::percent(total_threatened_prop, accuracy = 0.1), "\n",
  "(", format(total_threatened_count, big.mark = ","), " spp.)"
)

# --- 9b. Data for Middle Bar (Logging vs Plantation Breakdown) ---
# Use the pre-computed species-level threat flags
species_threat_summary <- extant_species_df %>%
  filter(has_wood_prod_ongoing) %>%
  distinct(scientific_name, has_logging_ongoing, has_plantation_ongoing) %>%
  rename(
    has_logging = has_logging_ongoing,
    has_plantation = has_plantation_ongoing
  )

# Calculate counts and proportions for the middle bar
plot_a2_data <- species_threat_summary %>%
  mutate(
    breakdown_label = case_when(
      has_logging & has_plantation ~ "Both",
      has_logging & !has_plantation ~ "Logging",
      !has_logging & has_plantation ~ "Plantations"
    )
  ) %>%
  count(breakdown_label, name = "count") %>%
  mutate(
    proportion = count / sum(count),
    plot_group = "A",
    # Set the factor level order for stacking
    breakdown_label = factor(breakdown_label, levels = c("Plantations", "Both", "Logging"))
  ) %>%
  arrange(breakdown_label) %>%
  mutate(label_y_pos = cumsum(proportion) - 0.5 * proportion)


# =============================================================================
# --- 10. VISUALIZATION FOR PLOT A (Horizontal Bars with Legends) ---
# =============================================================================

cat("Generating Plot A...\n")

# --- 10a. Color Palettes ---
# Palette for the 4-category classification (used in A1)
plot_a1_colors <- c(
  "No wood production threat" = "#cce5ff", 
  "Wood production threat & extinction risk" = "#A30000", 
  "Wood production threat & no extinction risk" = "#FF8A8A", 
  "Wood production threat & data deficient" = "#3B3B3B"
)

# Palette for the wood production threat type (used in A2)
plot_a2_colors <- c(
  "Logging" = "#d9be00",       # Yellow
  "Plantations" = "#0088d9",   # Blue
  "Both" = "#007b00"           # Green
)

# --- Shared Theme Settings for Legends ---
# This ensures all legends look consistent in size
shared_legend_theme <- theme(
  legend.position = "right",
  legend.title = element_text(size = 9, face = "bold"),
  legend.text = element_text(size = 8),
  legend.key.size = unit(0.4, "cm"),
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(0,0,0,-10)
)

# --- 10b. Plot A1 (Top Bar) ---
p1 <- ggplot(plot_a1_data, aes(x = plot_group, y = proportion, fill = threat_type_stacked)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), color = "white", size = 1) +
  annotate("text", x = 1.265, y = 0.25, label = "25%", color = "black", size = 2.8, vjust = 0) +
  annotate("text", x = 1.265, y = 0.50, label = "50%", color = "black", size = 2.8, vjust = 0) +
  annotate("text", x = 1.265, y = 0.75, label = "75%", color = "black", size = 2.8, vjust = 0) +
  annotate(
    "text", x = 1, y = -0.03,
    label = plot_a1_label,
    color = "black", size = 3.5,
    hjust = 1
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    name = "Overall Status", # Title for the first legend
    values = plot_a1_colors, 
    limits = plot_a1_factor_levels
  ) +
  scale_y_continuous(
    labels = c("", "25%", "50%", "75%", ""),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    name = NULL,
    limits = c(-0.08, 1.0),
    expand = c(0, 0)
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 55)
  ) +
  shared_legend_theme

# --- 10c. Plot A2 (Middle Bar - Type) ---
p2 <- ggplot(plot_a2_data, aes(x = plot_group, y = proportion, fill = breakdown_label)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE), width = 0.35) +
  geom_text(
    aes(x = 0.8, y = label_y_pos, label = ifelse(proportion > 0.01, scales::percent(proportion, accuracy = 0.1), "")),
    vjust = 1, color = "black", size = 3.5, inherit.aes = FALSE
  ) +
  coord_flip(clip="off") +
  scale_fill_manual(
    name = "Production Pathway", # Title for the second legend
    values = plot_a2_colors, 
    limits = levels(plot_a2_data$breakdown_label),
    labels = c("Plantation", "Both", "Natural Forestry")
  ) +
  scale_y_continuous(labels = NULL, breaks = NULL, name = NULL, limits = c(-0.08, 1.0), expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 5, r = 0, b = -5, l = 55)
  ) +
  shared_legend_theme

# --- 10d. Combine Plots with Collected Legends ---
# guides = 'collect' puts all legends into a single column on the right
combined_plot_a <- p1 + p2 +
  plot_layout(
    ncol = 1,
    heights = c(1.8, 0.9),
    guides = 'collect'
  )

# Display Plot A
print("Displaying Plot A (Horizontal Bars with Legends)...")
print(combined_plot_a)

cat("Saving high-resolution image...\n")

ggsave(
  filename = file.path(figure_output_dir, "figure1a_bars.png"),
  plot = combined_plot_a,
  width = 18,       # Wide format width in inches
  height = 5,       # Height in inches
  dpi = 300,        # High resolution (standard for print)
  scale = 0.8,        # Increase to >1 (e.g., 1.2) if text appears too small
  bg = "white"      # Ensure background is white, not transparent
)

cat("Done.\n")

# =============================================================================
# --- 11. DATA PREP FOR PLOT B (Faceted Donut Plots) ---
# =============================================================================
# This plot uses geom_rect to manually build the donut/pie charts

cat("Preparing data for Plot B (Inner Pie)...\n")

# --- 11a. Data for Inner Pie (Status) ---
# Radius: 0 to 2.0
plot_b_pie_data_rect <- all_data_classified %>%
  count(taxa_group, threat_type_stacked, name = "count") %>%
  group_by(taxa_group) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  arrange(taxa_group, factor(threat_type_stacked, levels = all_threat_labels_c_ordered)) %>%
  group_by(taxa_group) %>%
  mutate(
    ymax = cumsum(proportion),
    ymin = c(0, head(ymax, n=-1)),
    # labelPosition = (ymax + ymin)/2 # Not needed if we aren't labeling the inner pie
  ) %>%
  ungroup() %>%
  mutate(
    xmin_plot = 0,   
    xmax_plot = 2.0  
  )

# --- 11b. Data for Middle Donut (Type: Logging vs Plantation) ---
cat("Preparing data for Plot B (Middle Donut - Type)...\n")
# Radius: 2.5 to 3.0

# 1. Get pre-computed threatened species by taxon
threatened_species_by_taxa <- extant_species_df %>%
  filter(has_wood_prod_ongoing)

# 2. Calculate proportions
plot_b_donut_data <- threatened_species_by_taxa %>%
  mutate(
    has_logging = has_logging_ongoing,
    has_plantation = has_plantation_ongoing,
    breakdown_label = case_when(
      has_logging & has_plantation ~ "Both",
      has_logging & !has_plantation ~ "Logging",
      !has_logging & has_plantation ~ "Plantations"
    )
  ) %>%
  count(taxa_group, breakdown_label, name = "count") %>%
  group_by(taxa_group) %>%
  mutate(proportion = count / sum(count)) %>%
  mutate(
    breakdown_label = factor(
      breakdown_label,
      levels = c("Plantations", "Both", "Logging")
    )
  ) %>%
  ungroup()

# 3. Rescale to 2/3 arc
donut_arc_fraction <- 2/3
plot_b_donut_data_scaled <- plot_b_donut_data %>%
  mutate(
    proportion_raw = proportion,
    proportion = proportion * donut_arc_fraction
  )

# 4. Create empty filler data (1/3 arc)
all_taxa <- unique(plot_b_donut_data_scaled$taxa_group)
empty_donut_data <- tibble(
  taxa_group = all_taxa,
  breakdown_label = NA_character_,
  proportion_raw = NA_real_,
  proportion = 1 - donut_arc_fraction
)

# 5. Bind and Calculate Coordinates
plot_b_donut_data_rect <- bind_rows(plot_b_donut_data_scaled, empty_donut_data) %>%
  mutate(breakdown_label_factor = factor(breakdown_label, levels = levels(plot_a2_data$breakdown_label))) %>%
  group_by(taxa_group) %>%
  # Arrange by REVERSE factor level to match Plot A stacking
  arrange(taxa_group, is.na(breakdown_label_factor), desc(breakdown_label_factor)) %>%
  mutate(
    ymax = cumsum(proportion),
    ymin = c(0, head(ymax, n=-1))
  ) %>%
  ungroup() %>%
  mutate(
    xmin_plot = 2.5, 
    xmax_plot = 3.0
  )

plot_b_donut_label_data <- plot_b_donut_data_rect %>%
  filter(!is.na(breakdown_label), !is.na(proportion_raw)) %>%
  mutate(
    label_y_base = (ymin + ymax) / 2,
    # Nudge the smaller neighbouring slices so labels do not stack on top of
    # each other; this intentionally moves labels slightly off the segment
    # midpoint when adjacent slices are both small.
    label_y = case_when(
      breakdown_label == "Plantations" & proportion_raw < 0.15 ~ label_y_base + 0.02,
      breakdown_label == "Both" & proportion_raw < 0.15 ~ label_y_base - 0.01,
      TRUE ~ label_y_base
    ),
    # Stagger labels radially by production type to reduce collisions.
    label_x = case_when(
      breakdown_label == "Logging" ~ 3.8,
      breakdown_label == "Plantations" ~ 3.54,
      breakdown_label == "Both" ~ 3.34,
      TRUE ~ 3.5
    ),
    label = scales::percent(proportion_raw, accuracy = 1)
  )


# =============================================================================
# --- 12. VISUALIZATION FOR PLOT B (Faceted Donut Plots) ---
# =============================================================================

cat("Generating Plot B (Faceted Donuts)...\n")

plot_b <- ggplot() +
  
  # --- 1. Draw Inner Pie (Status) ---
  geom_rect(
    data = plot_b_pie_data_rect,
    aes(ymin = ymin, ymax = ymax, xmin = xmin_plot, xmax = xmax_plot, fill = threat_type_stacked),
    color = "white", size = 0.1
  ) +
  scale_fill_manual(
    name = "Threat Category",
    values = plot_a1_colors,
    limits = all_threat_labels_c_ordered
  ) +
  
  # --- 2. Draw Middle Donut (Type) ---
  new_scale_fill() +
  geom_rect(
    data = plot_b_donut_data_rect,
    aes(ymin = ymin, ymax = ymax, xmin = xmin_plot, xmax = xmax_plot, fill = breakdown_label),
    color = "white", size = 0.1
  ) +
  geom_text(
    data = plot_b_donut_label_data,
    aes(x = label_x, y = label_y, label = label),
    size = 4,
    color = "black",
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    name = "Wood Production Type",
    values = plot_a2_colors,
    limits = levels(plot_a2_data$breakdown_label),
    na.value = "transparent" # Hides the 1/3 gap
  ) +
  
  # --- 3. Layout & Theme ---
  coord_polar(theta = "y", start = 0, clip ="off") +
  facet_wrap(~ taxa_group, nrow = 1, drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "none", # We rely on Plot A's legends
    strip.text = element_text(size = 15, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),
    panel.spacing = unit(0.6, "lines")
  )

# Display Plot B
print("Displaying Plot B (Faceted Donuts)...")
print(plot_b)

cat("Saving high-resolution image for Plot B...\n")

ggsave(
  filename = file.path(figure_output_dir, "figure1b_donuts.png"),
  plot = plot_b,
  width = 16,       # Very wide to accommodate 5 taxa side-by-side
  height = 5,       # Slightly taller to allow strip labels and percent labels
  dpi = 300,        # High resolution
  bg = "white"      # Ensure white background
)

cat("Done saving Plot B.\n")

# =============================================================================
# --- 13. DATA PREP FOR PLOT C (Stacked Radial Bar Plot - Flipped) ---
# =============================================================================

cat("Preparing data for Flipped Stacked Radial Bar Plot C...\n")

# A. Data is already prepared in `all_data_classified`
# We just need to count by taxa *and* category
plot_c_counts <- all_data_classified %>%
  # Rename to 'fill_group' for generic use
  rename(fill_group = threat_type_stacked) %>%
  count(taxa_group, category, fill_group, name = "species_count")

# B. Calculate proportions *within* each bar (Taxa x Category)
plot_c_proportions <- plot_c_counts %>%
  group_by(taxa_group, category) %>%
  mutate(total_in_bar = sum(species_count)) %>%
  ungroup() %>%
  mutate(total_in_bar = ifelse(total_in_bar == 0, 1, total_in_bar)) %>%
  mutate(proportion = species_count / total_in_bar) %>%
  # Ensure fill_group is a factor using the master list
  mutate(fill_group = factor(fill_group, levels = all_threat_labels_c_ordered)) %>%
  # Ensure correct stacking order
  arrange(taxa_group, category, factor(fill_group, levels = all_threat_labels_c_ordered)) %>%
  # Calculate cumulative proportion for stacking (ymin/ymax)
  group_by(taxa_group, category) %>%
  mutate(prop_end = cumsum(proportion)) %>%
  mutate(prop_start = prop_end - proportion) %>%
  ungroup()

# C. Add "empty bars" to create gaps between groups
empty_bar_small <- 1 # The small gap between most taxa
empty_bar_large <- 5 # The large gap at the end for labels
data_to_plot <- data.frame()
bar_id_counter <- 1

for (taxa in taxa_levels) {
  taxa_data <- plot_c_proportions %>%
    filter(taxa_group == taxa) %>%
    mutate(bar_id = as.numeric(category) + bar_id_counter - 1)
  
  data_to_plot <- bind_rows(data_to_plot, taxa_data)
  bar_id_counter <- bar_id_counter + length(iucn_levels)
  
  if (taxa == taxa_levels[length(taxa_levels)]) {
    current_empty_bars <- empty_bar_large
  } else {
    current_empty_bars <- empty_bar_small
  }
  
  # Create the empty bars
  empty_data <- data.frame(
    taxa_group = taxa,
    bar_id = seq(bar_id_counter, bar_id_counter + current_empty_bars - 1),
    category = NA,
    # Use a valid factor level for the fill_group
    fill_group = all_threat_labels_c_ordered[1],
    prop_start = 0, prop_end = 0
  )
  data_to_plot <- bind_rows(data_to_plot, empty_data)
  bar_id_counter <- bar_id_counter + current_empty_bars
}

# D. Define the bar height (radius)
bar_outer_radius <- 12 # 0% line (outermost)
bar_height <- 10       # Max height (grows inward)
bar_inner_radius <- bar_outer_radius - bar_height # 100% line (innermost)
bar_width_factor <- 0.8 # Factor to reduce bar width for internal spacing

# Scale the prop_start/end to the bar height
data_to_plot <- data_to_plot %>%
  mutate(
    ymin_plot = bar_outer_radius - (prop_start * bar_height),
    ymax_plot = bar_outer_radius - (prop_end * bar_height),
    xmin_bar = bar_id - 0.5 * bar_width_factor,
    xmax_bar = bar_id + 0.5 * bar_width_factor
  )

# E. Prepare data for bar labels (DD, LC, etc.)
label_data <- data_to_plot %>%
  filter(!is.na(category)) %>%
  distinct(taxa_group, category, bar_id)
label_data$y_pos <- bar_outer_radius + 0.5

# F. Prepare data for group labels (Aves, Mammals...)
base_data <- data_to_plot %>%
  filter(!is.na(category)) %>%
  group_by(taxa_group) %>%
  summarize(title = mean(c(min(bar_id), max(bar_id))))
taxa_label_y_pos <- bar_outer_radius + 4

# G. Define positions for concentric grid circles
grid_percentages <- c(0.25, 0.50, 0.75)
grid_radii <- bar_outer_radius - (grid_percentages * bar_height)


# =============================================================================
# --- 14. VISUALIZATION: PLOT C (Flipped Stacked Radial Bar Plot) ---
# =============================================================================

cat("Generating Flipped Stacked Radial Bar Plot C...\n")

plot_c <- ggplot() +
  
  # The stacked bars
  geom_rect(data = data_to_plot,
            aes(xmin = xmin_bar, xmax = xmax_bar,
                ymin = ymin_plot, ymax = ymax_plot,
                fill = fill_group)) +
  
  # Concentric gridlines and labels
  geom_hline(yintercept = grid_radii, color = "white", linetype = "solid", size = 1) +
  annotate("text", x = 39, y = grid_radii[1], label = "25%", hjust = 1.1, vjust = 0.5, size = 3.5, color = "black", fontface="bold") +
  annotate("text", x = 39, y = grid_radii[2], label = "50%", hjust = 1.1, vjust = 0.5, size = 3.5, color = "black", fontface="bold") +
  annotate("text", x = 39, y = grid_radii[3], label = "75%", hjust = 1.1, vjust = 0.5, size = 3.5, color = "black", fontface="bold") +
  
  # IUCN category labels (DD, LC, etc.)
  geom_text(data = label_data,
            aes(x = bar_id, y = y_pos+0.5, label = category),
            color = "black", alpha = 0.8, size = 3.5,
            hjust = 0.5, vjust = 0.5, angle = 0,
            inherit.aes = FALSE) +
  
  # Taxa group labels (Aves, Mammalia, etc.)
  geom_text(data = base_data,
            aes(x = title, y = taxa_label_y_pos, label = taxa_group),
            hjust = 0.5, colour = "black", alpha = 0.8, size = 4.5,
            fontface = "bold", inherit.aes = FALSE) +
  
  # Use the same fill scale as Plot A
  scale_fill_manual(
    name = "Threat Category",
    values = plot_a1_colors,
    limits = all_threat_labels_c_ordered
  ) +
  
  # Set y-limits to make room for labels
  ylim(bar_inner_radius - 2, taxa_label_y_pos + 1) +
  
  theme_minimal() +
  theme(
    legend.position = "none", # Removed, as it's collected by patchwork
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    # Use negative margins to crop empty space
    plot.margin = margin(t = -0, r = 0, b = -0, l = 0),
    aspect.ratio = 1
  ) +
  coord_polar() +
  labs(title = NULL, subtitle = NULL)

# Display Plot C
print("Displaying Plot C (Radial Plot)...")
print(plot_c)

cat("Saving high-resolution image for Plot C...\n")

ggsave(
  filename = file.path(figure_output_dir, "figure1c_category_radial.png"),
  plot = plot_c,
  width = 16,       # Very wide to accommodate 5 taxa side-by-side
  height = 16,       # Short height to keep donuts circular
  dpi = 300,        # High resolution
  bg = "white",      # Ensure white background
  scale = 0.5
)

cat("Done saving Plot C.\n")


# --- 15. COMBINE & SAVE FINAL FIGURE (A, B, C) ---
# =============================================================================

cat("Combining all three plots (A, B, C) into final figure...\n")

final_plot = patchwork::wrap_plots(
  combined_plot_a,
  plot_b,
  plot_c,
  ncol = 1,
  heights = c(1, 1.5, 3)
) +
  plot_annotation(tag_levels = 'A') &
  
  # Apply a consistent theme to the *collected* legends
  
  theme(
    
    legend.title = element_text(size = 10, face="bold"),
    
    legend.position = "none",
    
    legend.text = element_text(size = 8),
    
    legend.key.size = unit(0.4, "cm")
    
  )



# --- Display the final combined plot ---

print("Displaying final combined plot...")

print(final_plot)

# --- Save the final plot as high-resolution PNG ---
cat("Saving final plot as high-resolution PNG...\n")
ggsave(
  filename = file.path(figure_output_dir, "figure1_combined.png"),
  plot = final_plot,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

# --- Save the final plot as high-quality PDF ---
cat("Saving final plot as high-quality PDF...\n")
ggsave(
  filename = file.path(figure_output_dir, "figure1_combined.pdf"),
  plot = final_plot,
  width = 10,
  height = 20,
  units = "in",
  device = cairo_pdf # Use cairo_pdf for better font embedding
)

cat("All plots generated and saved.\n")

# =============================================================================
# --- 16. SUMMARY OUTPUT: % THREATENED BY WOOD PRODUCTION, BY TAXON & FAMILY ---
# =============================================================================

cat("Generating wood production threat summaries by taxon and family...\n")

# --- 16a. Use the prepared species-level assessment summary ---
species_threat_flags <- all_assessment_df_consolidated %>%
  filter(!category %in% c("EX", "EW")) %>%
  mutate(
    taxa_group = factor(taxa_group, levels = taxa_levels),
    category = factor(category, levels = iucn_levels, ordered = TRUE),
    has_wood_prod = has_wood_prod_ongoing,
    has_logging = has_logging_ongoing,
    has_plantation = has_plantation_ongoing
  ) %>%
  transmute(
    scientific_name,
    taxa_group,
    order_name,
    family_name,
    category,
    has_wood_prod,
    has_logging,
    has_plantation
  ) %>%
  filter(!is.na(category))

# --- 16c. TAXON-LEVEL SUMMARY (unchanged) ---
taxon_summary <- species_threat_flags %>%
  group_by(taxa_group) %>%
  summarise(
    total_species           = n_distinct(scientific_name),
    n_wood_prod_threatened  = sum(has_wood_prod),
    n_logging_threatened    = sum(has_logging),
    n_plantation_threatened = sum(has_plantation),
    pct_wood_prod           = round(100 * n_wood_prod_threatened  / total_species, 2),
    pct_logging             = round(100 * n_logging_threatened    / total_species, 2),
    pct_plantation          = round(100 * n_plantation_threatened / total_species, 2),
    .groups = "drop"
  ) %>%
  arrange(taxa_group)

cat("\n========== TAXON-LEVEL SUMMARY ==========\n")
print(taxon_summary)

# --- 16d. ORDER-LEVEL SUMMARY ---
order_summary <- species_threat_flags %>%
  group_by(taxa_group, order_name) %>%
  summarise(
    total_species           = n_distinct(scientific_name),
    n_wood_prod_threatened  = sum(has_wood_prod),
    n_logging_threatened    = sum(has_logging),
    n_plantation_threatened = sum(has_plantation),
    pct_wood_prod           = round(100 * n_wood_prod_threatened  / total_species, 2),
    pct_logging             = round(100 * n_logging_threatened    / total_species, 2),
    pct_plantation          = round(100 * n_plantation_threatened / total_species, 2),
    .groups = "drop"
  ) %>%
  arrange(taxa_group, desc(pct_wood_prod))

cat("\n========== ORDER-LEVEL SUMMARY (top 10 per taxon by % wood production) ==========\n")
order_summary %>%
  group_by(taxa_group) %>%
  slice_max(pct_wood_prod, n = 10) %>%
  print(n = Inf)

# --- 16e. FAMILY-LEVEL SUMMARY (now includes order_name) ---
family_summary <- species_threat_flags %>%
  group_by(taxa_group, order_name, family_name) %>%
  summarise(
    total_species           = n_distinct(scientific_name),
    n_wood_prod_threatened  = sum(has_wood_prod),
    n_logging_threatened    = sum(has_logging),
    n_plantation_threatened = sum(has_plantation),
    pct_wood_prod           = round(100 * n_wood_prod_threatened  / total_species, 2),
    pct_logging             = round(100 * n_logging_threatened    / total_species, 2),
    pct_plantation          = round(100 * n_plantation_threatened / total_species, 2),
    .groups = "drop"
  ) %>%
  arrange(taxa_group, order_name, desc(pct_wood_prod))

cat("\n========== FAMILY-LEVEL SUMMARY (top 10 per taxon by % wood production) ==========\n")
family_summary %>%
  group_by(taxa_group) %>%
  slice_max(pct_wood_prod, n = 10) %>%
  print(n = Inf)

# --- 16f. CATEGORY-LEVEL SUMMARY (by taxa and IUCN category) ---
category_summary <- species_threat_flags %>%
  group_by(taxa_group, category) %>%
  summarise(
    total_species           = n_distinct(scientific_name),
    n_wood_prod_threatened  = sum(has_wood_prod),
    n_logging_threatened    = sum(has_logging),
    n_plantation_threatened = sum(has_plantation),
    pct_wood_prod           = round(100 * n_wood_prod_threatened  / total_species, 2),
    pct_logging             = round(100 * n_logging_threatened    / total_species, 2),
    pct_plantation          = round(100 * n_plantation_threatened / total_species, 2),
    .groups = "drop"
  ) %>%
  arrange(taxa_group, category)

cat("\n========== CATEGORY-LEVEL SUMMARY (by taxa and IUCN Red List category) ==========\n")
print(category_summary, n = Inf)

# --- 16g. EXPORT TO CSV ---
write.csv(taxon_summary,    file = here::here("outputs", "extent", "tables", "taxon_wood_production_summary.csv"),    row.names = FALSE)
write.csv(order_summary,    file = here::here("outputs", "extent", "tables", "order_wood_production_summary.csv"),    row.names = FALSE)
write.csv(family_summary,   file = here::here("outputs", "extent", "tables", "family_wood_production_summary.csv"),   row.names = FALSE)
write.csv(category_summary, file = here::here("outputs", "extent", "tables", "category_wood_production_summary.csv"), row.names = FALSE)

cat("\nSummary CSVs saved to /outputs/extent/tables/\n")
cat("  - taxon_wood_production_summary.csv\n")
cat("  - order_wood_production_summary.csv\n")
cat("  - family_wood_production_summary.csv\n")
cat("  - category_wood_production_summary.csv\n")

# =============================================================================
# --- 17. PLOT D: ORDER-LEVEL SCATTER PLOT (Threat Ratio vs. Total Species) ---
# =============================================================================

cat("Generating Plot D (Order-level scatter plot)...\n")

# --- 17a. Data Prep ---
plot_d_data <- order_summary %>%
  mutate(taxa_group = factor(taxa_group, levels = taxa_levels))

# Slope = total threatened species in class / total species in class (species-weighted)
class_slopes <- taxon_summary %>%
  mutate(slope = n_wood_prod_threatened / total_species)

# Calculate the threat ratio: observed / expected
# Ratio > 1 = more threatened than class average
# Ratio < 1 = less threatened than class average
# Ratio = 1 = exactly at class average (the trend line)
# Small offset added to handle orders with 0 threatened species on log scale
plot_d_data <- plot_d_data %>%
  left_join(class_slopes %>% select(taxa_group, slope, pct_wood_prod), by = "taxa_group") %>%
  mutate(
    expected_threatened = total_species * slope,
    threat_ratio        = (n_wood_prod_threatened + 0.01) / expected_threatened
  )

# --- 17b. Label data with factor levels preserved ---
label_data_d <- bind_rows(
  plot_d_data %>% group_by(taxa_group) %>% slice_max(threat_ratio, n = 5, with_ties = FALSE),
  plot_d_data %>% group_by(taxa_group) %>% slice_min(threat_ratio, n = 5, with_ties = FALSE)
) %>%
  distinct(taxa_group, order_name, .keep_all = TRUE) %>%
  ungroup() %>%
  mutate(taxa_group = factor(taxa_group, levels = taxa_levels))

# --- 17c. Phantom points along the trend line (now y = 1 horizontal) ---
trend_phantom_points <- plot_d_data %>%
  group_by(taxa_group) %>%
  reframe(
    total_species = seq(min(total_species), max(total_species), length.out = 50),
    threat_ratio  = 1.0
  ) %>%
  mutate(taxa_group = factor(taxa_group, levels = taxa_levels))

# --- 17d. Color Palette ---
taxa_colors <- c(
  "Aves"            = "#4E79A7",
  "Mammalia"        = "#F28E2B",
  "Amphibia"        = "#59A14F",
  "Reptilia"        = "#E15759",
  "Freshwater Fish" = "#76B7B2"
)

# Legend labels showing class average %
slope_labels <- class_slopes %>%
  mutate(slope_label = paste0(taxa_group, " (", round(pct_wood_prod, 1), "% avg)")) %>%
  select(taxa_group, slope_label)

# --- 17e. Build the Plot ---
plot_d <- ggplot() +
  
  # Horizontal reference line at ratio = 1 (the class average)
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  
  # Per-class horizontal trend lines, coloured by taxa
  # These are all at y = 1 but coloured individually for visual association
  geom_segment(
    data = plot_d_data %>%
      group_by(taxa_group) %>%
      summarise(x_min = min(total_species), x_max = max(total_species), .groups = "drop"),
    aes(x = x_min, xend = x_max, y = 1, yend = 1, color = taxa_group),
    linetype = "dashed", linewidth = 0.8, alpha = 0.7
  ) +
  
  # Phantom points along trend line — ggrepel obstacles
  geom_point(
    data  = trend_phantom_points,
    aes(x = total_species, y = threat_ratio),
    size  = 0, alpha = 0, color = NA
  ) +
  
  # Order-level points
  geom_point(
    data = plot_d_data %>% mutate(taxa_group = factor(taxa_group, levels = taxa_levels)),
    aes(x = total_species, y = threat_ratio, color = taxa_group),
    size = 3, alpha = 0.85
  ) +
  
  # Labels for top/bottom 5 per class
  ggrepel::geom_text_repel(
    data               = label_data_d,
    aes(x = total_species, y = threat_ratio, label = order_name, color = taxa_group),
    size               = 2.8,
    max.overlaps       = 50,
    show.legend        = FALSE,
    box.padding        = 0.6,
    point.padding      = 0.8,
    force              = 3,
    force_pull         = 0.5,
    min.segment.length = 0.2,
    segment.alpha      = 0.6,
    segment.size       = 0.3
  ) +
  
  scale_color_manual(
    name   = "Taxonomic Class",
    values = taxa_colors,
    labels = setNames(slope_labels$slope_label, slope_labels$taxa_group)
  ) +
  scale_x_continuous(
    name   = "Total number of species in order",
    labels = scales::comma,
    trans  = scales::pseudo_log_trans(base = 10),
    breaks = c(0, 1, 10, 100, 1000, 10000),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  # scale_y_continuous(
  #   name   = "Threat ratio (observed / expected threatened species)",
  #   trans  = "log10",  # true log scale — symmetric above and below y = 1
  #   breaks = c(0.01, 0.1, 1, 10, 100),
  #   labels = function(x) ifelse(x == 1, "1 (avg)", as.character(x)),
  #   expand = expansion(mult = c(0.05, 0.05))
  # ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    axis.title      = element_text(size = 11),
    axis.text       = element_text(size = 9),
    plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40"),
    plot.margin     = margin(10, 15, 10, 10)
  ) +
  labs(
    title    = "Wood Production Threats Across Vertebrate Orders",
    subtitle = "Y-axis shows each order's threat ratio: observed threatened species / expected based on class average.\nRatio = 1 (dashed line) means the order matches its class average exactly.\nOn a log scale, distances above and below the line are now symmetric."
  )

print(plot_d)

# --- 17f. Save ---
ggsave(
  filename = file.path(figure_output_dir, "figure1d_order_scatter.png"),
  plot     = plot_d,
  width    = 12,
  height   = 9,
  dpi      = 300,
  bg       = "white"
)

cat("Done saving Plot D.\n")
