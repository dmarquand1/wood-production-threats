# =============================================================================
#
# FIGURE S6 — PUBLICATION VERSION (Biological Conservation)
#
#
# Description:
# Circular phylogenies of the five taxa (A amphibians, B birds, C freshwater
# fishes, D mammals, E reptiles) with branch colour = proportion of
# descendant tips threatened by wood production, an inner ring of tip
# threat status by production pathway, and labelled families (black arcs +
# names). Families are bold if top 10 by number of threatened species,
# underlined if top 10 by proportion threatened, bold-underlined if both
# (same rules as build_phylogram_plot_set() in 03_make_phylogeny_plots.R).
#
# Data are loaded exactly as 03_make_phylogeny_plots.R does (cached full tree
# + cached matched IUCN table, standardised with the same helper functions,
# FishTree taxonomy for fish), so the trees and threat flags are identical.
#
# Layout. Each tree is drawn in Cartesian millimetre coordinates from
# ggtree's own circular layout (same ladderised tip order and node
# positions), so every label can be placed and checked exactly:
#   * labels sit in left/right columns hugging the circle, spaced vertically
#     (order preserved), with a leader line wherever a label moved away from
#     its family; families near the top/bottom may go on either side, and the
#     assignment with no conflicts is chosen;
#   * every candidate layout is checked for overlapping label boxes, crossing
#     leader lines, leaders passing through other labels, labels cutting into
#     the rings and leaders cutting across the circle;
#   * the largest tree radius for which all five panels pass every check and
#     fit their 95 mm column is used, and each row is cropped to its content,
#     so white space is minimal.
#
# Output (outputs/publication_figures/Figure_S6/):
#   - Figure_S6.pdf   vector, fonts embedded (cairo_pdf)   <- submit this
#   - Figure_S6.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_S6_preview.png  300-dpi preview for checking only
#   - Figure_S6_labelled_families.csv  families labelled in each panel
#
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(patchwork)
  library(here)
})

# --- 1. SHARED HELPERS FROM 03_make_phylogeny_plots.R ---
# -----------------------------------------------------------------------------
# Evaluate only the named function definitions from script 03, so the data
# loading here cannot drift from the analysis script.
phylo_script <- here::here("scripts", "01_extent", "03_make_phylogeny_plots.R")
needed_functions <- c(
  "clean_names", "to_phylo_style", "to_logical_safe", "first_non_na", "first_field",
  "standardize_match_table", "load_fishtree_taxonomy", "apply_fishtree_taxonomy"
)
for (expr in parse(phylo_script)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      as.character(expr[[2]]) %in% needed_functions) {
    eval(expr, envir = globalenv())
  }
}

taxa <- tibble::tribble(
  ~taxon,      ~tag, ~silhouette,
  "amphibian", "A",  "amphibian sil.png",
  "bird",      "B",  "bird sil.png",
  "fish",      "C",  "fish sil.png",
  "mammal",    "D",  "mammal sil.png",
  "reptile",   "E",  "reptile sil.png"
)

out_dir        <- here::here("outputs", "publication_figures", "Figure_S6")
silhouette_dir <- here::here("archive", "reference_material", "species_silhouettes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# --- 2. FIGURE SETTINGS (mm at print size) ---
# -----------------------------------------------------------------------------
fig_width_mm <- 190
col_w        <- 95              # each column
row_gap      <- 1.5             # vertical gap between rows
r_search     <- seq(27, 17, by = -0.25)   # tree radii tried, largest first

# Radial offsets beyond the tips (mm)
ring_in_off  <- 0.35            # tip threat-status ring
ring_out_off <- 1.85
fam_in_off   <- 2.3             # labelled-family arcs
fam_out_off  <- 2.75
anchor_off   <- 3.05            # leader lines start here
label_off    <- 4.3             # labels hug this circle

font_family   <- "Arial"
label_pt      <- 6
size_tag      <- 10
tag_box       <- c(w = 4.5, h = 4.5)   # space reserved for the panel letter
branch_lwd    <- 0.13           # ggplot linewidth (~0.28 pt)
label_pad_mm  <- 1.0            # vertical gap between stacked labels (clears underlines)
leader_min_mm <- 0.6            # draw a leader line if a label moved further
clearance_mm  <- 0.25           # minimum clearance used by the layout checks
flexible_cos  <- 0.35           # families this close to the top/bottom may go on either side
max_flexible  <- 8

phylo_ring_interval_myr <- 20
phylo_ring_fill_col     <- "grey96"

col_none            <- "#d3d3d3"
col_logging_only    <- "#d9be00"
col_plantation_only <- "#0088d9"
col_both            <- "#007b00"
branch_low_col  <- "#2B6CB0"
branch_mid_col  <- "#F3EFCF"
branch_high_col <- "#D73027"

tip_colours <- c(
  "None" = col_none, "Plantation only" = col_plantation_only,
  "Both" = col_both, "Logging only" = col_logging_only
)


# --- 3. DATA (as build_taxon_phylogram() in script 03, steps 1-3) ---
# -----------------------------------------------------------------------------
load_taxon_tree <- function(taxon) {
  processed_dir <- here::here("data", "processed", "traits", taxon)
  tree_file <- if (taxon == "fish") {
    file.path(processed_dir, "final_fish_tree_pruned.tre")
  } else {
    file.path(processed_dir, paste0(taxon, "_full_tree.tre"))
  }
  tree <- read.tree(tree_file)

  lookup_file <- here::here("data", "processed", "assessments", paste0(taxon, "_assessment_taxonomy_lookup.rds"))
  iucn_lookup <- readRDS(lookup_file) %>%
    transmute(
      IUCN_Name = clean_names(IUCN_Name),
      Family_lookup = clean_names(Family),
      Order_lookup = clean_names(Order)
    )
  matched <- standardize_match_table(
    read_csv(file.path(processed_dir, paste0(taxon, "_iucn_tree_matched.csv")), show_col_types = FALSE),
    iucn_lookup = iucn_lookup
  )
  if (taxon == "fish") {
    fish_taxonomy <- load_fishtree_taxonomy(
      here::here("data", "raw", "trees", "fish", "fishtree_taxonomy.csv"),
      tree_tips = tree$tip.label
    )
    matched <- apply_fishtree_taxonomy(matched, fish_taxonomy, label = taxon)
  }

  tip_data <- tibble(label = to_phylo_style(tree$tip.label)) %>%
    left_join(matched, by = c("label" = "phylo_name")) %>%
    mutate(
      Threat_Logging = replace_na(Threat_Logging, FALSE),
      Threat_Plantation = replace_na(Threat_Plantation, FALSE),
      Wood_Threat = Threat_Logging | Threat_Plantation,
      Family = str_to_title(Family),  # some sources store family names in capitals
      tip_category = case_when(
        Threat_Logging & Threat_Plantation ~ "Both",
        Threat_Logging ~ "Logging only",
        Threat_Plantation ~ "Plantation only",
        TRUE ~ "None"
      )
    ) %>%
    filter(!is.na(IUCN_Name))

  tree$tip.label <- to_phylo_style(tree$tip.label)
  tree <- keep.tip(tree, tip_data$label)
  tip_data <- tip_data %>% arrange(match(label, tree$tip.label))
  stopifnot(identical(tip_data$label, tree$tip.label))
  list(tree = tree, tip_data = tip_data)
}


# --- 4. FAMILY LABEL SELECTION (as build_phylogram_plot_set()) ---
# -----------------------------------------------------------------------------
select_label_families <- function(tree, tip_data, monophyletic_only) {
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
    mutate(is_monophyletic = map_lgl(
      family_tips, ~ if (length(.x) > 1) ape::is.monophyletic(tree, .x) else TRUE
    ))

  candidates <- if (monophyletic_only) filter(family_summary, is_monophyletic) else family_summary
  top_count <- candidates %>% arrange(desc(n_wood_threat), Family) %>% slice_head(n = 10) %>% pull(Family)
  top_prop  <- candidates %>% filter(n_wood_threat > 0) %>%
    arrange(desc(prop_wood_threat), desc(n_wood_threat), Family) %>% slice_head(n = 10) %>% pull(Family)

  family_summary %>%
    mutate(top_count = Family %in% top_count, top_prop = Family %in% top_prop) %>%
    filter(top_count | top_prop) %>%
    select(-family_tips)
}


# --- 5. GEOMETRY HELPERS ---
# -----------------------------------------------------------------------------
deg2rad <- function(d) d * pi / 180

# Many arcs at once (vectorised): one group per element of r / a0 / a1
arcs_df <- function(r, a0, a1, group, step_deg = 0.5) {
  n <- pmax(2, ceiling(abs(a1 - a0) / step_deg) + 1)
  idx <- rep(seq_along(r), n)
  frac <- (sequence(n) - 1) / rep(n - 1, n)
  a <- deg2rad(a0[idx] + (a1 - a0)[idx] * frac)
  tibble(x = r[idx] * cos(a), y = r[idx] * sin(a), group = group[idx], i = idx)
}

# Annular sector polygons between radii r0-r1 and angles a0-a1 (vectorised)
sectors_df <- function(r0, r1, a0, a1, group, step_deg = 0.5) {
  k <- length(a0)
  outer <- arcs_df(rep(r1, k), a0, a1, group, step_deg)
  inner <- arcs_df(rep(r0, k), a1, a0, group, step_deg)
  bind_rows(mutate(outer, part = 1), mutate(inner, part = 2)) %>%
    arrange(i, part) %>%
    select(x, y, group)
}

# Width/height (mm) of each label exactly as drawn (the parsed plotmath
# expression, so bold and underline are included), measured on a Cairo
# device so it matches the PDF output
measure_labels <- function(exprs) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  map_dfr(exprs, function(e) {
    g <- grid::textGrob(parse(text = e)[[1]], gp = grid::gpar(fontfamily = font_family, fontsize = label_pt))
    tibble(
      w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
      h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE)
    )
  })
}

# Size (mm) of any grob, measured on a Cairo device
grob_size_mm <- function(g) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(tmp)
  on.exit({ grDevices::dev.off(); unlink(tmp) })
  c(
    w = grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE),
    h = grid::convertHeight(grid::grobHeight(g), "mm", valueOnly = TRUE)
  )
}

plotmath_label <- function(family, bold, underline) {
  txt <- paste0('"', family, '"')
  if (underline) txt <- paste0("underline(", txt, ")")
  if (bold) txt <- paste0("bold(", txt, ")")
  txt
}

read_silhouette <- function(file) {
  img <- png::readPNG(file.path(silhouette_dir, file))
  alpha <- img[, , dim(img)[3]]
  keep_rows <- range(which(rowSums(alpha > 0.05) > 0))
  keep_cols <- range(which(colSums(alpha > 0.05) > 0))
  alpha <- alpha[keep_rows[1]:keep_rows[2], keep_cols[1]:keep_cols[2]]
  out <- array(0, dim = c(dim(alpha), 4))
  out[, , 4] <- alpha
  out
}


# --- 6. LABEL LAYOUT ---
# -----------------------------------------------------------------------------
# 1-D packing: move positions as little as possible so that consecutive items
# (kept in order) are at least (size_i + size_j)/2 + pad apart and lie within
# [lo, hi]. Overlapping neighbours merge into a block laid out edge to edge,
# centred on its members' preferred positions and clamped to the bounds.
pack_1d <- function(desired, size, lo, hi, pad) {
  n <- length(desired)
  if (n == 0) return(numeric(0))
  ord <- order(desired)
  d <- desired[ord]; s <- size[ord]
  first <- seq_len(n); last <- seq_len(n)
  place <- function(f, l) {
    offs <- c(0, cumsum((s[f:l][-1] + s[f:l][-(l - f + 1)]) / 2 + pad))
    p <- mean(d[f:l] - offs) + offs
    lo_b <- lo + s[f] / 2; hi_b <- hi - s[l] / 2
    if (p[1] < lo_b) p <- p + (lo_b - p[1])
    if (tail(p, 1) > hi_b) p <- p - (tail(p, 1) - hi_b)
    p
  }
  repeat {
    merged <- FALSE
    b <- 1
    while (b < length(first)) {
      p1 <- place(first[b], last[b]); p2 <- place(first[b + 1], last[b + 1])
      if (p2[1] - tail(p1, 1) < (s[last[b]] + s[first[b + 1]]) / 2 + pad) {
        last[b] <- last[b + 1]
        first <- first[-(b + 1)]; last <- last[-(b + 1)]
        merged <- TRUE
      } else {
        b <- b + 1
      }
    }
    if (!merged) break
  }
  out <- numeric(n)
  out[ord] <- unlist(map2(first, last, place))
  out
}

# Place labels for one side assignment (side = +1 right, -1 left)
place_labels <- function(labels, side, R) {
  label_r <- R + label_off
  y_lim <- label_r + 12
  labels %>%
    mutate(
      side = side,
      desired_y = label_r * sin(deg2rad(anchor_angle))
    ) %>%
    group_by(side) %>%
    mutate(ly = pack_1d(desired_y, h, -y_lim, y_lim, label_pad_mm)) %>%
    ungroup() %>%
    mutate(
      hug_x = side * sqrt(pmax(label_r^2 - ly^2, 1.5^2)),
      anchor_x = (R + anchor_off) * cos(deg2rad(anchor_angle)),
      # labels never start inward of their anchor, so leaders point away from the text
      lx = if_else(side > 0, pmax(hug_x, anchor_x + 0.6), pmin(hug_x, anchor_x - 0.6)),
      hjust = if_else(side > 0, 0, 1),
      xmin = if_else(side > 0, lx, lx - w), xmax = if_else(side > 0, lx + w, lx),
      ymin = ly - h / 2, ymax = ly + h / 2,
      ax = (R + anchor_off) * cos(deg2rad(anchor_angle)),
      ay = (R + anchor_off) * sin(deg2rad(anchor_angle)),
      bx = lx - side * 0.3, by = ly,
      moved_mm = sqrt((lx - label_r * cos(deg2rad(anchor_angle)))^2 + (ly - desired_y)^2),
      has_leader = moved_mm > leader_min_mm
    )
}

# Geometry tests
segments_cross <- function(p1, p2, q1, q2) {
  orient <- function(a, b, c) sign((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
  o1 <- orient(p1, p2, q1); o2 <- orient(p1, p2, q2)
  o3 <- orient(q1, q2, p1); o4 <- orient(q1, q2, p2)
  o1 != o2 && o3 != o4 && o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0
}
segment_hits_box <- function(p1, p2, box) {  # box = c(xmin, xmax, ymin, ymax)
  inside <- function(p) p[1] > box[1] && p[1] < box[2] && p[2] > box[3] && p[2] < box[4]
  if (inside(p1) || inside(p2)) return(TRUE)
  corners <- list(c(box[1], box[3]), c(box[2], box[3]), c(box[2], box[4]), c(box[1], box[4]))
  any(map_lgl(1:4, function(k) segments_cross(p1, p2, corners[[k]], corners[[k %% 4 + 1]])))
}
dist_origin_segment <- function(p1, p2) {
  d <- p2 - p1
  t <- max(0, min(1, -sum(p1 * d) / max(sum(d * d), 1e-12)))
  sqrt(sum((p1 + t * d)^2))
}
dist_origin_box <- function(b) {
  dx <- max(b[1], 0, -b[2]); dy <- max(b[3], 0, -b[4])
  sqrt(dx^2 + dy^2)
}

# Count every kind of conflict in a placed layout
layout_conflicts <- function(L, R) {
  n <- nrow(L)
  boxes <- lapply(seq_len(n), function(i) c(L$xmin[i], L$xmax[i], L$ymin[i], L$ymax[i]))
  grow <- function(b, m) b + c(-m, m, -m, m)
  lead <- which(L$has_leader)
  seg <- function(i) list(c(L$ax[i], L$ay[i]), c(L$bx[i], L$by[i]))

  text <- 0; cross <- 0; through <- 0
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    a <- grow(boxes[[i]], clearance_mm / 2); b <- grow(boxes[[j]], clearance_mm / 2)
    if (a[1] < b[2] && b[1] < a[2] && a[3] < b[4] && b[3] < a[4]) text <- text + 1
  }
  for (i in lead) {
    s <- seg(i)
    for (j in lead[lead > i]) {
      t <- seg(j)
      if (segments_cross(s[[1]], s[[2]], t[[1]], t[[2]])) cross <- cross + 1
    }
    for (j in setdiff(seq_len(n), i)) {
      if (segment_hits_box(s[[1]], s[[2]], grow(boxes[[j]], clearance_mm / 2))) through <- through + 1
    }
  }
  into_rings <- sum(map_dbl(boxes, dist_origin_box) < R + fam_out_off + clearance_mm)
  across <- sum(map_dbl(lead, function(i) { s <- seg(i); dist_origin_segment(s[[1]], s[[2]]) }) <
                  R + anchor_off - 0.05)
  c(text = text, cross = cross, through = through, into_rings = into_rings, across = across)
}

# Best layout at radius R: families near the top/bottom may go on either
# side; try the side assignments and keep the one with no conflicts that is
# narrowest (then moves labels least)
best_layout <- function(labels, R) {
  base_side <- if_else(cos(deg2rad(labels$anchor_angle)) >= 0, 1, -1)
  flexible <- which(abs(cos(deg2rad(labels$anchor_angle))) < flexible_cos)
  flexible <- head(flexible[order(abs(cos(deg2rad(labels$anchor_angle[flexible]))))], max_flexible)
  combos <- if (length(flexible) == 0) list(integer(0)) else
    unlist(lapply(0:length(flexible), function(k) combn(flexible, k, simplify = FALSE)), recursive = FALSE)

  best <- NULL
  for (flip in combos) {
    side <- base_side
    side[flip] <- -side[flip]
    L <- place_labels(labels, side, R)
    conf <- layout_conflicts(L, R)
    score <- c(sum(conf), max(L$xmax, R + fam_out_off) - min(L$xmin, -(R + fam_out_off)), sum(L$moved_mm))
    if (is.null(best) || score[1] < best$score[1] ||
        (score[1] == best$score[1] && (score[2] < best$score[2] - 1e-6 ||
          (abs(score[2] - best$score[2]) < 1e-6 && score[3] < best$score[3])))) {
      best <- list(labels = L, conflicts = conf, score = score)
    }
  }
  best
}

# Content bounding box of a panel (circle + labels + leaders)
content_bbox <- function(L, R) {
  r_out <- R + fam_out_off
  c(
    xmin = min(-r_out, L$xmin, L$ax), xmax = max(r_out, L$xmax, L$ax),
    ymin = min(-r_out, L$ymin, L$ay), ymax = max(r_out, L$ymax, L$ay)
  )
}


# --- 7. PREPARE EACH TAXON (radius-independent parts) ---
# -----------------------------------------------------------------------------
prepare_taxon <- function(taxon, tag, silhouette_file) {
  cat("Preparing", taxon, "\n")
  dat <- load_taxon_tree(taxon)
  tree <- dat$tree
  tip_data <- dat$tip_data
  n_tips <- length(tree$tip.label)

  # Proportion of descendant tips threatened, for every node (postorder sums)
  n_desc <- c(rep(1, n_tips), rep(0, tree$Nnode))
  n_thr  <- c(as.numeric(tip_data$Wood_Threat), rep(0, tree$Nnode))
  post <- reorder(tree, "postorder")$edge
  for (k in seq_len(nrow(post))) {
    n_desc[post[k, 1]] <- n_desc[post[k, 1]] + n_desc[post[k, 2]]
    n_thr[post[k, 1]]  <- n_thr[post[k, 1]] + n_thr[post[k, 2]]
  }
  prop_threat <- n_thr / n_desc

  # ggtree's circular layout: depth (x), angle per node; radii as fractions
  # of the tip radius so they scale with R
  layout <- ggtree(tree, layout = "circular")$data %>% select(node, parent, x, angle, isTip)
  depth <- max(layout$x)
  lay <- layout %>%
    mutate(u = x / depth, value = prop_threat[node]) %>%
    left_join(select(layout, parent = node, parent_u = x, parent_angle = angle), by = "parent") %>%
    mutate(parent_u = parent_u / depth) %>%
    filter(node != parent)

  arcs_unit <- arcs_df(lay$parent_u, lay$parent_angle, lay$angle, lay$node) %>% mutate(value = lay$value[i])
  radials_unit <- lay %>% transmute(
    x0 = parent_u * cos(deg2rad(angle)), y0 = parent_u * sin(deg2rad(angle)),
    x1 = u * cos(deg2rad(angle)), y1 = u * sin(deg2rad(angle)), value
  )
  band_starts <- seq(0, depth, by = phylo_ring_interval_myr * 2)
  bands_unit <- tibble(u0 = band_starts / depth, u1 = pmin(band_starts + phylo_ring_interval_myr, depth) / depth)

  tip_step <- 360 / n_tips
  tips <- lay %>% filter(isTip) %>%
    mutate(label = tree$tip.label[node]) %>%
    left_join(select(tip_data, label, tip_category, Family), by = "label") %>%
    arrange(angle)

  label_families <- select_label_families(tree, tip_data, monophyletic_only = taxon == "fish")
  runs <- tips %>%
    mutate(run_id = cumsum(c(TRUE, head(Family, -1) != tail(Family, -1)))) %>%
    filter(Family %in% label_families$Family) %>%
    group_by(Family, run_id) %>%
    summarise(a0 = min(angle) - tip_step / 2, a1 = max(angle) + tip_step / 2, n = n(), .groups = "drop")

  labels <- runs %>%
    group_by(Family) %>%
    slice_max(n, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(Family, anchor_angle = (a0 + a1) / 2) %>%
    left_join(label_families, by = "Family") %>%
    mutate(expr = pmap_chr(list(Family, top_count, top_prop), plotmath_label))
  labels <- bind_cols(labels, measure_labels(labels$expr))

  list(
    taxon = taxon, tag = tag, silhouette = read_silhouette(silhouette_file),
    arcs_unit = arcs_unit, radials_unit = radials_unit, bands_unit = bands_unit,
    tips = tips, tip_step = tip_step, runs = runs, labels = labels,
    families = label_families %>% mutate(taxon = taxon, .before = 1)
  )
}

prepared <- pmap(list(taxa$taxon, taxa$tag, taxa$silhouette), prepare_taxon)


# --- 8. CHOOSE THE TREE RADIUS ---
# -----------------------------------------------------------------------------
# Largest radius at which every panel is conflict-free and fits its column
cat("\nSearching for the largest conflict-free tree radius...\n")
chosen <- NULL
for (R in r_search) {
  layouts <- map(prepared, ~ best_layout(.x$labels, R))
  total_conflicts <- sum(map_dbl(layouts, ~ sum(.x$conflicts)))
  widths <- map_dbl(layouts, ~ { b <- content_bbox(.x$labels, R); b[["xmax"]] - b[["xmin"]] })
  if (total_conflicts == 0 && all(widths <= col_w - 1)) {
    chosen <- list(R = R, layouts = layouts)
    break
  }
}
if (is.null(chosen)) stop("No conflict-free layout found in the radius search range")
R <- chosen$R
cat(sprintf("Tree radius %.2f mm\n", R))
for (k in seq_along(prepared)) {
  cat(sprintf("  %-9s %2d labels | conflicts: %s\n", prepared[[k]]$taxon, nrow(prepared[[k]]$labels),
              paste(names(chosen$layouts[[k]]$conflicts), chosen$layouts[[k]]$conflicts, sep = "=", collapse = " ")))
}


# --- 9. PANELS ---
# -----------------------------------------------------------------------------
build_panel <- function(prep, layout, R, xlim, ylim) {
  L <- layout$labels
  arcs <- prep$arcs_unit %>% mutate(x = x * R, y = y * R)
  radials <- prep$radials_unit %>% mutate(across(c(x0, y0, x1, y1), ~ .x * R))
  bands <- pmap_dfr(list(prep$bands_unit$u0, prep$bands_unit$u1, seq_len(nrow(prep$bands_unit))), function(u0, u1, k) {
    outer <- arcs_df(u1 * R, 0, 360, k, 1) %>% mutate(hole = 1)
    inner <- if (u0 > 0) arcs_df(u0 * R, 360, 0, k, 1) %>% mutate(hole = 2)
    bind_rows(outer, inner) %>% mutate(band = k)
  })
  tip_ring <- sectors_df(R + ring_in_off, R + ring_out_off,
                         prep$tips$angle - prep$tip_step / 2, prep$tips$angle + prep$tip_step / 2,
                         prep$tips$node, 1) %>%
    left_join(select(prep$tips, group = node, cat = tip_category), by = "group")
  fam_arcs <- sectors_df(R + fam_in_off, R + fam_out_off, prep$runs$a0, prep$runs$a1,
                         paste(prep$runs$Family, prep$runs$run_id), 0.25)
  leaders <- L %>% filter(has_leader)

  img <- prep$silhouette
  a <- ncol(img) / nrow(img)
  box_w <- 0.62 * R; box_h <- 0.45 * R
  sil_w <- if (a > box_w / box_h) box_w else box_h * a
  sil_h <- if (a > box_w / box_h) box_w / a else box_h

  ggplot() +
    geom_polygon(data = bands, aes(x, y, group = band, subgroup = hole), fill = phylo_ring_fill_col, colour = NA) +
    geom_path(data = arcs, aes(x, y, group = group, colour = value), linewidth = branch_lwd, lineend = "butt") +
    geom_segment(data = radials, aes(x = x0, y = y0, xend = x1, yend = y1, colour = value),
                 linewidth = branch_lwd, lineend = "butt") +
    scale_colour_gradient2(low = branch_low_col, mid = branch_mid_col, high = branch_high_col,
                           midpoint = 0.5, limits = c(0, 1), guide = "none") +
    geom_polygon(data = tip_ring, aes(x, y, group = group, fill = cat), colour = NA) +
    scale_fill_manual(values = tip_colours, guide = "none") +
    geom_polygon(data = fam_arcs, aes(x, y, group = group), fill = "black", colour = NA) +
    geom_segment(data = leaders, aes(x = ax, y = ay, xend = bx, yend = by), linewidth = 0.2, colour = "black") +
    geom_text(data = L, aes(x = lx, y = ly, label = expr, hjust = hjust), vjust = 0.5,
              parse = TRUE, size = label_pt / .pt, family = font_family) +
    annotation_raster(img, xmin = -sil_w / 2, xmax = sil_w / 2, ymin = -sil_h / 2, ymax = sil_h / 2,
                      interpolate = TRUE) +
    annotate("text", x = xlim[1] + 0.3, y = ylim[2] - 0.3, label = prep$tag, hjust = 0, vjust = 1,
             size = size_tag / .pt, fontface = "bold", family = font_family) +
    coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

# Frame each panel: content centred horizontally in its column; the panel
# letter sits in the top-left corner, with extra height only if the content
# would reach into it
frames <- map2(prepared, chosen$layouts, function(prep, lay) {
  b <- content_bbox(lay$labels, R)
  cx <- (b[["xmin"]] + b[["xmax"]]) / 2
  xlim <- cx + c(-col_w / 2, col_w / 2)
  tag_x <- c(xlim[1], xlim[1] + tag_box[["w"]])
  tag_hits <- function(top) {
    tag <- c(tag_x, top - tag_box[["h"]], top)
    circle_hit <- dist_origin_box(tag) < R + fam_out_off
    label_hit <- any(lay$labels$xmin < tag[2] & lay$labels$xmax > tag[1] &
                       lay$labels$ymin < tag[4] & lay$labels$ymax > tag[3])
    circle_hit || label_hit
  }
  edge_pad <- 0.8   # text renders slightly beyond its measured box (underlines, ascenders)
  top <- b[["ymax"]] + edge_pad
  if (tag_hits(top)) top <- top + tag_box[["h"]]
  list(xlim = xlim, top = top, bottom = b[["ymin"]] - edge_pad)
})


# --- 10. LEGEND ---
# -----------------------------------------------------------------------------
legend_plot <- ggplot() +
  geom_tile(
    data = tibble(cat = factor(names(tip_colours), levels = names(tip_colours)), x = 1:4, y = 1),
    aes(x, y, fill = cat)
  ) +
  geom_point(data = tibble(x = 1:2, y = 1, v = c(0, 1)), aes(x, y, colour = v)) +
  scale_fill_manual(
    values = tip_colours,
    breaks = c("None", "Plantation only", "Both", "Logging only"),
    labels = c("None", "Plantation", "Both", "Natural Forestry"),
    name = "Production Pathway",
    guide = guide_legend(order = 1)
  ) +
  scale_colour_gradient2(
    low = branch_low_col, mid = branch_mid_col, high = branch_high_col,
    midpoint = 0.5, limits = c(0, 1), breaks = seq(0, 1, 0.25),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    name = "Proportion of descendants threatened",
    guide = guide_colourbar(order = 2, theme = theme(
      legend.key.height = unit(20, "mm"), legend.key.width = unit(3, "mm"), legend.ticks = element_blank()
    ))
  ) +
  theme_void(base_family = font_family) +
  theme(
    legend.position = "left",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.title = element_text(size = 7, margin = margin(b = 1.5, unit = "mm")),
    legend.text = element_text(size = 6.5, margin = margin(l = 1.2, unit = "mm")),
    legend.key.size = unit(3.2, "mm"),
    legend.spacing.y = unit(4, "mm")
  )
legend_gtable <- ggplotGrob(legend_plot)
legend_grob <- legend_gtable$grobs[[which(legend_gtable$layout$name == "guide-box-left")]]


rows <- list(c(1, 2), c(3, 4), c(5))
legend_h <- grob_size_mm(legend_grob)[["h"]] + 2   # measured legend height
row_h <- map_dbl(rows, function(r) {
  h <- max(map_dbl(frames[r], ~ .x$top - .x$bottom))
  if (length(r) == 1) h <- max(h, legend_h)
  h + row_gap
})

panel_plots <- map(seq_along(prepared), function(k) {
  r <- which(map_lgl(rows, ~ k %in% .x))
  fr <- frames[[k]]
  ylim <- c(fr$top - row_h[r], fr$top)
  build_panel(prepared[[k]], chosen$layouts[[k]], R, fr$xlim, ylim)
})



# --- 11. ASSEMBLE & SAVE ---
# -----------------------------------------------------------------------------
figure_s6 <- wrap_plots(c(panel_plots, list(wrap_elements(full = legend_grob))), ncol = 2) +
  plot_layout(widths = unit(rep(col_w, 2), "mm"), heights = unit(row_h, "mm")) &
  theme(plot.margin = margin(0, 0, 0, 0))   # incl. the legend cell, so the grid is exactly 190 mm wide

fig_height_mm <- sum(row_h)
out_base <- file.path(out_dir, "Figure_S6")
ggsave(paste0(out_base, ".pdf"), figure_s6,
       width = fig_width_mm, height = fig_height_mm, units = "mm", device = cairo_pdf)
ggsave(paste0(out_base, ".tiff"), figure_s6,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 1000, device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(paste0(out_base, "_preview.png"), figure_s6,
       width = fig_width_mm, height = fig_height_mm, units = "mm",
       dpi = 300, device = ragg::agg_png, bg = "white")

write_csv(bind_rows(map(prepared, "families")), file.path(out_dir, "Figure_S6_labelled_families.csv"))
cat(sprintf("Saved Figure_S6 (%.0f x %.0f mm) to %s\n", fig_width_mm, fig_height_mm, out_dir))

# Final check on the rendered image: nothing may touch the outer edge
# (anything there has been clipped)
preview <- png::readPNG(paste0(out_base, "_preview.png"))
ink <- apply(preview[, , 1:3], c(1, 2), min) < 0.95
edge_ink <- c(top = sum(ink[1, ]), bottom = sum(ink[nrow(ink), ]), left = sum(ink[, 1]), right = sum(ink[, ncol(ink)]))
cat("Pixels touching each edge of the rendered figure (should all be 0):\n")
print(edge_ink)
