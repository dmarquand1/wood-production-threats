# =============================================================================
#
# FIGURE 1 — PUBLICATION VERSION (Biological Conservation)
#
#
# Description:
# Rebuilds the three-part Figure 1 from 02_make_figure1.R at final print size
# for Elsevier / Biological Conservation:
#   - A: Horizontal stacked bars — overall wood-production threat status and
#        the production pathway breakdown of wood-threatened species.
#   - B: Five donut charts, one per taxon, showing the same data as A.
#   - C: Stacked radial bar chart by taxon and IUCN Red List category.
# Taxon silhouettes are drawn in all three panels.
#
# Output (outputs/publication_figures/Figure_1/):
#   - Figure_1.pdf   vector, fonts embedded (cairo_pdf)   <- submit this
#   - Figure_1.eps   vector, fonts embedded (cairo_ps)
#   - Figure_1.tiff  1000-dpi LZW-compressed bitmap (fallback)
#   - Figure_1_preview.png  300-dpi preview for checking only
#
# Sizing follows the Elsevier artwork guide: double-column width = 190 mm;
# all text is 6–8 pt at final size (Arial). Panels A and B are drawn in
# millimetre-based coordinates so their printed size is known exactly and the
# layout can be packed with minimal white space.
#
# =============================================================================


# --- 1. SETUP ---
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(here)
library(scales)
library(patchwork)

figure_output_dir <- here::here("outputs", "publication_figures", "Figure_1")
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)

# Final print width (Elsevier double column). Height is derived from the
# panel heights in section 7.
fig_width_mm <- 190

# Typography: all sizes in pt at final print size. geom_text sizes are in mm,
# so convert with ggplot2::.pt.
font_family <- "Arial"
pt <- function(x) x / ggplot2::.pt
size_body   <- 7    # legend text, bar labels
size_small  <- 6    # IUCN category labels, gridline labels
size_title  <- 8    # taxon names, legend titles
size_tag    <- 10   # panel tags

gridline_width <- 0.35  # white separator lines (mm)


# --- 2. LOAD PREPARED ASSESSMENT SUMMARIES ---
# -----------------------------------------------------------------------------
load_assessment_summary <- function(taxon_key) {
  readRDS(
    here::here(
      "data", "processed", "assessments",
      paste0(taxon_key, "_assessment_species_summary.rds")
    )
  )
}

all_assessment_df_consolidated <- bind_rows(
  load_assessment_summary("bird"),
  load_assessment_summary("mammal"),
  load_assessment_summary("amphibian"),
  load_assessment_summary("reptile"),
  load_assessment_summary("fish")
)

iucn_levels <- c("DD", "LC", "NT", "VU", "EN", "CR")
taxa_levels <- c("Amphibia", "Aves", "Freshwater Fish", "Mammalia", "Reptilia")

all_assessment_df_consolidated <- all_assessment_df_consolidated %>%
  mutate(
    category = factor(category, levels = iucn_levels, ordered = TRUE),
    taxa_group = factor(taxa_group, levels = taxa_levels)
  ) %>%
  filter(!is.na(category))

extant_species_df <- all_assessment_df_consolidated %>%
  filter(!category %in% c("EX", "EW"))

all_threat_labels_c_ordered <- c(
  "Wood production threat & extinction risk",
  "Wood production threat & no extinction risk",
  "Wood production threat & data deficient",
  "No wood production threat"
)

all_data_classified <- extant_species_df %>%
  mutate(
    has_wood_harvesting_threat = has_wood_prod_ongoing,
    is_threatened_iucn = category %in% c("VU", "EN", "CR")
  ) %>%
  mutate(threat_type_stacked = case_when(
    has_wood_harvesting_threat & is_threatened_iucn ~ all_threat_labels_c_ordered[1],
    has_wood_harvesting_threat & !is_threatened_iucn & category %in% c("LC", "NT") ~ all_threat_labels_c_ordered[2],
    has_wood_harvesting_threat & category == "DD" ~ all_threat_labels_c_ordered[3],
    !has_wood_harvesting_threat ~ all_threat_labels_c_ordered[4]
  )) %>%
  filter(!is.na(threat_type_stacked))


# --- 3. COLOUR PALETTES (unchanged from 02_make_figure1.R) ---
# -----------------------------------------------------------------------------
plot_a1_colors <- c(
  "No wood production threat" = "#cce5ff",
  "Wood production threat & extinction risk" = "#A30000",
  "Wood production threat & no extinction risk" = "#FF8A8A",
  "Wood production threat & data deficient" = "#3B3B3B"
)

plot_a2_colors <- c(
  "Logging" = "#d9be00",
  "Plantations" = "#0088d9",
  "Both" = "#007b00"
)

# Legend order / wording shared by all panels
status_legend_levels <- c(
  "No wood production threat",
  "Wood production threat & no extinction risk",
  "Wood production threat & extinction risk",
  "Wood production threat & data deficient"
)
pathway_levels <- c("Plantations", "Both", "Logging")
pathway_labels <- c("Plantation", "Both", "Natural Forestry")


# --- 3b. TAXON SILHOUETTES ---
# -----------------------------------------------------------------------------
# Recoloured to solid black and cropped to their visible extent so that the
# requested sizes apply to the animal itself, not the transparent padding.
silhouette_dir <- here::here("archive", "reference_material", "species_silhouettes")
silhouette_files <- c(
  "Amphibia"        = "amphibian sil.png",
  "Aves"            = "bird sil.png",
  "Freshwater Fish" = "fish sil.png",
  "Mammalia"        = "mammal sil.png",
  "Reptilia"        = "reptile sil.png"
)

read_silhouette <- function(file) {
  img <- png::readPNG(file.path(silhouette_dir, file))
  alpha <- img[, , dim(img)[3]]
  keep_rows <- which(rowSums(alpha > 0.05) > 0)
  keep_cols <- which(colSums(alpha > 0.05) > 0)
  alpha <- alpha[range(keep_rows)[1]:range(keep_rows)[2], range(keep_cols)[1]:range(keep_cols)[2]]
  out <- array(0, dim = c(dim(alpha), 4))
  out[, , 4] <- alpha
  out
}
silhouettes <- lapply(silhouette_files, read_silhouette)
silhouette_aspect <- sapply(silhouettes, function(x) ncol(x) / nrow(x))  # width / height

# Width and height of a silhouette fitted inside a box, aspect ratio preserved.
# Box and result are in the same (isotropic) units.
fit_silhouette <- function(taxon, box_w, box_h) {
  a <- silhouette_aspect[[taxon]]
  if (a > box_w / box_h) c(w = box_w, h = box_w / a) else c(w = box_h * a, h = box_h)
}

# annotation_raster() layer centred on (x, y), for plots whose x and y units
# are both millimetres.
silhouette_layer <- function(taxon, x, y, box_w, box_h) {
  s <- fit_silhouette(taxon, box_w, box_h)
  annotation_raster(
    silhouettes[[taxon]],
    xmin = x - s[["w"]] / 2, xmax = x + s[["w"]] / 2,
    ymin = y - s[["h"]] / 2, ymax = y + s[["h"]] / 2,
    interpolate = TRUE
  )
}

# Raster grob centred on (x, y) in npc of a square viewport
silhouette_grob <- function(taxon, x, y, box_w, box_h) {
  s <- fit_silhouette(taxon, box_w, box_h)
  grid::rasterGrob(
    silhouettes[[taxon]], x = x, y = y,
    width = unit(s[["w"]], "npc"), height = unit(s[["h"]], "npc")
  )
}

# Shared theme pieces
base_theme <- theme(
  text = element_text(family = font_family, colour = "black"),
  plot.tag = element_text(size = size_tag, face = "bold", family = font_family, hjust = 0, vjust = 1),
  plot.tag.position = c(0, 1),
  legend.title = element_text(size = size_title, face = "bold", margin = margin(b = 1.5, unit = "mm")),
  legend.text = element_text(size = size_body, margin = margin(l = 1.5, unit = "mm")),
  legend.key.size = unit(3.2, "mm"),
  legend.key.spacing.y = unit(0.8, "mm"),
  legend.spacing.y = unit(3, "mm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(0, "mm"),
  legend.justification = c(0, 0.5)
)

# Printed width of the panel area (tags are overlaid, so panels span the page)
panel_width_mm <- fig_width_mm - 1


# =============================================================================
# --- 4. PANEL A: HORIZONTAL BARS ---
# =============================================================================
# x is in mm across the page; y is in mm up the panel.

plot_a1_data <- all_data_classified %>%
  count(threat_type_stacked, name = "count") %>%
  mutate(proportion = count / sum(count))

# Species threatened by wood production (any IUCN category)
total_threatened_count <- sum(plot_a1_data$count[plot_a1_data$threat_type_stacked != "No wood production threat"])
total_species_count    <- sum(plot_a1_data$count)
total_threatened_prop  <- total_threatened_count / total_species_count

plot_a1_label <- paste0(
  scales::percent(total_threatened_prop, accuracy = 0.1), "\n",
  "(", format(total_threatened_count, big.mark = ","), " spp.)"
)

plot_a2_data <- extant_species_df %>%
  filter(has_wood_prod_ongoing) %>%
  distinct(scientific_name, has_logging_ongoing, has_plantation_ongoing) %>%
  mutate(
    breakdown_label = case_when(
      has_logging_ongoing & has_plantation_ongoing ~ "Both",
      has_logging_ongoing & !has_plantation_ongoing ~ "Logging",
      !has_logging_ongoing & has_plantation_ongoing ~ "Plantations"
    )
  ) %>%
  count(breakdown_label, name = "count") %>%
  mutate(
    proportion = count / sum(count),
    breakdown_label = factor(breakdown_label, levels = c("Plantations", "Both", "Logging"))
  ) %>%
  arrange(breakdown_label)

# Horizontal geometry (mm): room for the summary label, then the bars
a_label_w <- 21
a_bar_x0  <- a_label_w
a_bar_w   <- panel_width_mm - a_label_w - 1
bar_x <- function(p) a_bar_x0 + p * a_bar_w

# Vertical geometry (mm, bottom to top)
a_pct_label_h  <- 3.2
bottom_bar     <- c(ymin = 3.8, ymax = 9.8)
bracket_low_y  <- 10.8   # bracket over the whole bottom bar
bracket_high_y <- 14.2   # bracket under the wood-threatened part of the top bar
bracket_tick   <- 0.8
top_bar        <- c(ymin = 15.2, ymax = 29.2)
a_height_mm    <- top_bar[["ymax"]] + 3.2

# Top bar, left to right: extinction risk, no extinction risk, DD, no threat
plot_a1_rect <- plot_a1_data %>%
  mutate(threat_type_stacked = factor(threat_type_stacked, levels = all_threat_labels_c_ordered)) %>%
  arrange(threat_type_stacked) %>%
  mutate(xmax = bar_x(cumsum(proportion)), xmin = bar_x(cumsum(proportion) - proportion))

plot_a2_rect <- plot_a2_data %>%
  mutate(xmax = bar_x(cumsum(proportion)), xmin = bar_x(cumsum(proportion) - proportion))

threat_end_x <- bar_x(total_threatened_prop)
bracket_lines <- tibble(
  x    = c(a_bar_x0, a_bar_x0, threat_end_x, (a_bar_x0 + threat_end_x) / 2, a_bar_x0, a_bar_x0, bar_x(1)),
  xend = c(a_bar_x0, threat_end_x, threat_end_x, (a_bar_x0 + threat_end_x) / 2, bar_x(1), a_bar_x0, bar_x(1)),
  y    = c(bracket_high_y - bracket_tick, bracket_high_y, bracket_high_y, bracket_high_y,
           bracket_low_y, bracket_low_y, bracket_low_y),
  yend = c(bracket_high_y, bracket_high_y, bracket_high_y - bracket_tick, bracket_low_y,
           bracket_low_y, bracket_low_y + bracket_tick, bracket_low_y + bracket_tick)
)

# Silhouettes in the "no threat" part of the top bar
a_sil_box_h <- 0.78 * (top_bar[["ymax"]] - top_bar[["ymin"]])
a_sil_box_w <- 15
a_sil_gap   <- 3
a_sil_sizes <- sapply(taxa_levels, fit_silhouette, box_w = a_sil_box_w, box_h = a_sil_box_h)
a_sil_total <- sum(a_sil_sizes["w", ]) + a_sil_gap * (length(taxa_levels) - 1)
a_sil_left  <- bar_x(0.61) - a_sil_total / 2 +
  c(0, cumsum(a_sil_sizes["w", ] + a_sil_gap))[seq_along(taxa_levels)]
a_sil_x     <- a_sil_left + a_sil_sizes["w", ] / 2

a_silhouette_layers <- lapply(seq_along(taxa_levels), function(i) {
  silhouette_layer(taxa_levels[i], a_sil_x[i], mean(top_bar), a_sil_box_w, a_sil_box_h)
})

p_a <- ggplot() +
  geom_rect(
    data = plot_a1_rect,
    aes(xmin = xmin, xmax = xmax, ymin = top_bar[["ymin"]], ymax = top_bar[["ymax"]], fill = threat_type_stacked)
  ) +
  scale_fill_manual(
    name = "Overall Status",
    values = plot_a1_colors,
    breaks = status_legend_levels,
    guide = guide_legend(order = 1)
  ) +
  new_scale_fill() +
  geom_rect(
    data = plot_a2_rect,
    aes(xmin = xmin, xmax = xmax, ymin = bottom_bar[["ymin"]], ymax = bottom_bar[["ymax"]], fill = breakdown_label)
  ) +
  scale_fill_manual(
    name = "Production Type",
    values = plot_a2_colors,
    breaks = pathway_levels,
    labels = pathway_labels,
    guide = guide_legend(order = 2)
  ) +
  annotate(
    "segment", x = bar_x(c(0.25, 0.5, 0.75)), xend = bar_x(c(0.25, 0.5, 0.75)),
    y = top_bar[["ymin"]], yend = top_bar[["ymax"]],
    colour = "white", linewidth = gridline_width
  ) +
  annotate(
    "text", x = bar_x(c(0.25, 0.5, 0.75)), y = top_bar[["ymax"]] + 0.5,
    label = c("25%", "50%", "75%"),
    size = pt(size_small), family = font_family, vjust = 0
  ) +
  annotate(
    "text", x = a_bar_x0 - 1.5, y = mean(top_bar), label = plot_a1_label,
    size = pt(size_body), family = font_family, hjust = 1, lineheight = 0.9
  ) +
  a_silhouette_layers +
  geom_segment(
    data = bracket_lines,
    aes(x = x, xend = xend, y = y, yend = yend),
    linewidth = 0.3, lineend = "square"
  ) +
  geom_text(
    data = plot_a2_rect,
    aes(x = (xmin + xmax) / 2, y = bottom_bar[["ymin"]] - 0.6, label = scales::percent(proportion, accuracy = 0.1)),
    vjust = 1, size = pt(size_body), family = font_family
  ) +
  scale_x_continuous(limits = c(0, panel_width_mm), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, a_height_mm), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  labs(tag = "A") +
  theme_void() +
  base_theme +
  theme(plot.margin = margin(0, 0, 0, 0))


# =============================================================================
# --- 5. PANEL B: DONUTS PER TAXON ---
# =============================================================================
# Donuts are drawn as polygons in mm so that they can be packed edge to edge
# (coord_polar() pads every panel). Angles are fractions of a full turn,
# clockwise from 12 o'clock, as in coord_polar(theta = "y").

plot_b_pie_data <- all_data_classified %>%
  count(taxa_group, threat_type_stacked, name = "count") %>%
  group_by(taxa_group) %>%
  mutate(proportion = count / sum(count)) %>%
  arrange(taxa_group, factor(threat_type_stacked, levels = all_threat_labels_c_ordered)) %>%
  mutate(t_end = cumsum(proportion), t_start = t_end - proportion) %>%
  ungroup()

plot_b_donut_data <- extant_species_df %>%
  filter(has_wood_prod_ongoing) %>%
  mutate(
    breakdown_label = case_when(
      has_logging_ongoing & has_plantation_ongoing ~ "Both",
      has_logging_ongoing & !has_plantation_ongoing ~ "Logging",
      !has_logging_ongoing & has_plantation_ongoing ~ "Plantations"
    )
  ) %>%
  count(taxa_group, breakdown_label, name = "count") %>%
  group_by(taxa_group) %>%
  mutate(proportion_raw = count / sum(count)) %>%
  ungroup()

# Outer ring covers 2/3 of the circle, running clockwise from 12 o'clock:
# Natural Forestry, Both, Plantation; the last 1/3 is left empty.
donut_arc_fraction <- 2 / 3
plot_b_ring_data <- plot_b_donut_data %>%
  mutate(breakdown_label = factor(breakdown_label, levels = c("Logging", "Both", "Plantations"))) %>%
  arrange(taxa_group, breakdown_label) %>%
  group_by(taxa_group) %>%
  mutate(
    t_end = cumsum(proportion_raw) * donut_arc_fraction,
    t_start = t_end - proportion_raw * donut_arc_fraction
  ) %>%
  ungroup()

# Geometry (mm). Radii in units of the original design (pie = 2, ring 2.5–3).
b_unit     <- 5.15                # mm per design unit
pie_r      <- 2.0 * b_unit
ring_r_in  <- 2.5 * b_unit
ring_r_out <- 3.0 * b_unit
b_label_r  <- ring_r_out + c(Logging = 2.3, Plantations = 2.1, Both = 1.8)  # text sits a fixed distance (mm) outside the ring
b_cell_w   <- panel_width_mm / length(taxa_levels)
b_centre_x <- b_cell_w * (seq_along(taxa_levels) - 0.5) - 0.06 * b_cell_w  # labels overhang more on the right
names(b_centre_x) <- taxa_levels
b_bottom   <- ring_r_out + 3.3    # lowest label below the centre
b_top      <- ring_r_out + 0.4
b_height_mm <- b_bottom + b_top
b_centre_y <- b_bottom

arc_polygon <- function(cx, cy, r_in, r_out, t0, t1, n = 120) {
  t <- seq(t0, t1, length.out = max(2, ceiling(n * (t1 - t0))))
  outer <- tibble(x = cx + r_out * sin(2 * pi * t), y = cy + r_out * cos(2 * pi * t))
  inner <- if (r_in > 0) {
    tibble(x = cx + r_in * sin(2 * pi * rev(t)), y = cy + r_in * cos(2 * pi * rev(t)))
  } else {
    tibble(x = cx, y = cy)
  }
  bind_rows(outer, inner)
}

make_polygons <- function(segments, fill_col, r_in, r_out) {
  segments %>%
    mutate(seg_id = paste(taxa_group, .data[[fill_col]], sep = "_")) %>%
    rowwise() %>%
    reframe(
      seg_id, fill = as.character(.data[[fill_col]]),
      arc_polygon(b_centre_x[[as.character(taxa_group)]], b_centre_y, r_in, r_out, t_start, t_end)
    )
}

pie_polygons  <- make_polygons(plot_b_pie_data, "threat_type_stacked", 0, pie_r)
ring_polygons <- make_polygons(plot_b_ring_data, "breakdown_label", ring_r_in, ring_r_out)

plot_b_label_data <- plot_b_ring_data %>%
  mutate(
    t_mid = (t_start + t_end) / 2,
    # Nudge small neighbouring slices apart so their labels do not overlap
    t_label = case_when(
      breakdown_label == "Plantations" & proportion_raw < 0.15 ~ t_mid + 0.02,
      breakdown_label == "Both" & proportion_raw < 0.15 ~ t_mid - 0.01,
      TRUE ~ t_mid
    ),
    r_label = b_label_r[as.character(breakdown_label)],
    x = b_centre_x[as.character(taxa_group)] + r_label * sin(2 * pi * t_label),
    y = b_centre_y + r_label * cos(2 * pi * t_label),
    label = scales::percent(proportion_raw, accuracy = 1)
  )

# Silhouettes sit in the left of each pie (the "no threat" side), clear of the
# threatened slices that start at 12 o'clock and run clockwise.
b_sil_box_w <- 0.85 * pie_r
b_sil_box_h <- 0.65 * pie_r
b_sil_dx    <- -0.5 * pie_r
b_sil_dy    <- 0.1 * pie_r
b_silhouette_layers <- lapply(taxa_levels, function(taxon) {
  silhouette_layer(taxon, b_centre_x[[taxon]] + b_sil_dx, b_centre_y + b_sil_dy, b_sil_box_w, b_sil_box_h)
})

p_b <- ggplot() +
  geom_polygon(
    data = pie_polygons, aes(x = x, y = y, group = seg_id, fill = fill),
    colour = "white", linewidth = 0.15
  ) +
  scale_fill_manual(values = plot_a1_colors, guide = "none") +
  new_scale_fill() +
  geom_polygon(
    data = ring_polygons, aes(x = x, y = y, group = seg_id, fill = fill),
    colour = "white", linewidth = 0.15
  ) +
  scale_fill_manual(values = plot_a2_colors, guide = "none") +
  b_silhouette_layers +
  geom_text(
    data = plot_b_label_data, aes(x = x, y = y, label = label),
    size = pt(size_body), family = font_family, fontface = "bold"
  ) +
  scale_x_continuous(limits = c(0, panel_width_mm), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, b_height_mm), expand = c(0, 0)) +
  coord_fixed(clip = "off") +
  labs(tag = "B") +
  theme_void() +
  base_theme +
  theme(plot.margin = margin(0, 0, 0, 0))


# =============================================================================
# --- 6. PANEL C: STACKED RADIAL BARS BY IUCN CATEGORY ---
# =============================================================================

plot_c_proportions <- all_data_classified %>%
  rename(fill_group = threat_type_stacked) %>%
  count(taxa_group, category, fill_group, name = "species_count") %>%
  group_by(taxa_group, category) %>%
  mutate(proportion = species_count / sum(species_count)) %>%
  ungroup() %>%
  mutate(fill_group = factor(fill_group, levels = all_threat_labels_c_ordered)) %>%
  arrange(taxa_group, category, fill_group) %>%
  group_by(taxa_group, category) %>%
  mutate(prop_end = cumsum(proportion), prop_start = prop_end - proportion) %>%
  ungroup()

# Gaps between taxa: one empty bar between groups, a wider gap at the end
empty_bar_small <- 1
empty_bar_large <- 5
data_to_plot <- data.frame()
bar_id_counter <- 1

for (taxa in taxa_levels) {
  taxa_data <- plot_c_proportions %>%
    filter(taxa_group == taxa) %>%
    mutate(bar_id = as.numeric(category) + bar_id_counter - 1)
  data_to_plot <- bind_rows(data_to_plot, taxa_data)
  bar_id_counter <- bar_id_counter + length(iucn_levels)

  current_empty_bars <- if (taxa == tail(taxa_levels, 1)) empty_bar_large else empty_bar_small
  bar_id_counter <- bar_id_counter + current_empty_bars
}
n_bar_slots <- bar_id_counter - 1

bar_outer_radius <- 12   # 0% line (outermost)
bar_height <- 10         # bars grow inward
bar_inner_radius <- bar_outer_radius - bar_height
bar_width_factor <- 0.8

data_to_plot <- data_to_plot %>%
  mutate(
    ymin_plot = bar_outer_radius - (prop_start * bar_height),
    ymax_plot = bar_outer_radius - (prop_end * bar_height),
    xmin_bar = bar_id - 0.5 * bar_width_factor,
    xmax_bar = bar_id + 0.5 * bar_width_factor
  )

label_data <- data_to_plot %>%
  distinct(taxa_group, category, bar_id) %>%
  mutate(y_pos = bar_outer_radius + 1)

taxa_label_y_pos <- bar_outer_radius + 3.4
base_data <- data_to_plot %>%
  group_by(taxa_group) %>%
  summarise(title = mean(range(bar_id)), .groups = "drop")

grid_percentages <- c(0.25, 0.50, 0.75)
grid_radii <- bar_outer_radius - (grid_percentages * bar_height)
c_r_max <- taxa_label_y_pos + 0.45

p_c <- ggplot() +
  geom_rect(
    data = data_to_plot,
    aes(xmin = xmin_bar, xmax = xmax_bar, ymin = ymin_plot, ymax = ymax_plot, fill = fill_group)
  ) +
  geom_hline(yintercept = grid_radii, colour = "white", linewidth = gridline_width) +
  annotate(
    "text", x = n_bar_slots, y = grid_radii, label = c("25%", "50%", "75%"),
    hjust = 1.1, vjust = 0.5, size = pt(size_small), family = font_family, fontface = "bold"
  ) +
  geom_text(
    data = label_data,
    aes(x = bar_id, y = y_pos, label = category),
    size = pt(size_small), family = font_family, colour = "grey20", inherit.aes = FALSE
  ) +
  geom_text(
    data = base_data,
    aes(x = title, y = taxa_label_y_pos, label = taxa_group),
    size = pt(size_title), family = font_family, fontface = "bold", colour = "grey20", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = plot_a1_colors) +
  scale_x_continuous(limits = c(0.5, n_bar_slots + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, c_r_max), expand = c(0, 0)) +
  coord_polar(clip = "off") +
  theme_void() +
  base_theme +
  theme(
    aspect.ratio = 1,
    legend.position = "none",  # shares the Overall Status legend from panel A
    plot.margin = margin(0, 0, 0, 0)
  )

# Silhouettes at the inner ends of each taxon's bars. coord_polar() maps the
# radius onto 0–0.4 npc around the panel centre, with angle 0 at the top
# running clockwise, so positions are computed directly in panel npc.
c_sil_radius <- 4.3
c_sil_positions <- base_data %>%
  mutate(
    theta = 2 * pi * (title - 0.5) / n_bar_slots,
    r_npc = 0.4 * c_sil_radius / c_r_max,
    x = 0.5 + r_npc * sin(theta),
    y = 0.5 + r_npc * cos(theta)
  )

p_c_gtable <- ggplotGrob(p_c)
c_panel <- which(grepl("^panel", p_c_gtable$layout$name))
for (i in seq_len(nrow(c_sil_positions))) {
  p_c_gtable$grobs[[c_panel]] <- grid::grobTree(
    p_c_gtable$grobs[[c_panel]],
    silhouette_grob(
      as.character(c_sil_positions$taxa_group[i]),
      x = c_sil_positions$x[i], y = c_sil_positions$y[i],
      box_w = 0.085, box_h = 0.065
    )
  )
}

# coord_polar() pads the circle to 80% of the panel (radius 0.4 npc). Draw the
# panel in a viewport enlarged by 1/0.8 so the labelled circle fills its cell.
radial_scale <- 1 / 0.8
p_c_filled <- wrap_elements(
  full = grid::grobTree(
    p_c_gtable,
    vp = grid::viewport(width = radial_scale, height = radial_scale)
  ),
  clip = FALSE
) +
  labs(tag = "C") +
  base_theme


# =============================================================================
# --- 7. ASSEMBLE ---
# =============================================================================
# Panels A and B span the full width; panel C sits bottom-left with the shared
# legends (from panel A) collected into the space on its right.

c_size_mm <- 96
gap_mm    <- 1   # vertical breathing room between panels

layout_design <- "
AA
##
BB
##
CE
"

figure_1 <- wrap_plots(
  A = p_a,
  B = p_b,
  C = p_c_filled,
  E = guide_area(),
  design = layout_design,
  heights = unit(c(a_height_mm, gap_mm, b_height_mm, gap_mm, c_size_mm), "mm"),
  widths = unit(c(c_size_mm + 8, fig_width_mm - c_size_mm - 8), "mm"),
  guides = "collect"
) &
  theme(legend.position = "right", legend.box.just = "left", plot.margin = margin(0, 0, 0, 0))

fig_height_mm <- a_height_mm + b_height_mm + c_size_mm + 2 * gap_mm + 1


# =============================================================================
# --- 8. SAVE ---
# =============================================================================

out_base <- file.path(figure_output_dir, "Figure_1")

# Vector PDF with embedded fonts (primary submission file)
ggsave(
  paste0(out_base, ".pdf"), figure_1,
  width = fig_width_mm, height = fig_height_mm, units = "mm",
  device = cairo_pdf
)

# Vector EPS with embedded fonts
ggsave(
  paste0(out_base, ".eps"), figure_1,
  width = fig_width_mm, height = fig_height_mm, units = "mm",
  device = cairo_ps, fallback_resolution = 1000
)

# 1000-dpi bitmap (line art requirement)
ggsave(
  paste0(out_base, ".tiff"), figure_1,
  width = fig_width_mm, height = fig_height_mm, units = "mm",
  dpi = 1000, device = ragg::agg_tiff, compression = "lzw", bg = "white"
)

# Low-resolution preview for checking the layout only
ggsave(
  paste0(out_base, "_preview.png"), figure_1,
  width = fig_width_mm, height = fig_height_mm, units = "mm",
  dpi = 300, device = ragg::agg_png, bg = "white"
)

cat("Figure 1 saved to", figure_output_dir, "\n")
cat(sprintf("Figure size: %.0f x %.0f mm\n", fig_width_mm, fig_height_mm))
cat("Species threatened by wood production:", total_threatened_count, "of", total_species_count,
    sprintf("(%.1f%%)\n", 100 * total_threatened_prop))
