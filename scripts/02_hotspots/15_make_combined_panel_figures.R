rm(list = ls())
gc()

# =============================================================================
# COMBINED PANEL FIGURES: BIVARIATE AND GRADIENT (5 taxa x raw/SES)
#
# Description:
# Composites the already-rendered per-taxon supplementary PNGs (from
# 14_make_transformed_supplementary_maps.R) into two combined panel figures,
# each a 5-row x 2-column grid (row order: amphibian, bird, fish, mammal,
# reptile; left column = raw score, right column = SES), tagged A-J, with no
# spacing between panels.
#   1. combined_bivariate_panel.png -- shared bivariate legend inset bottom-right
#   2. combined_gradient_panel.png  -- each panel keeps its own legend (taxa
#      have very different value ranges, so a shared scale would be misleading)
# =============================================================================

local_lib <- normalizePath(".r-lib", winslash = "/", mustWork = FALSE)
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

required_packages <- c("ggplot2", "patchwork", "png", "here")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))

out_dir <- here::here("outputs", "hotspots", "figures", "supplementary_transformed")
taxa_order <- c("amphibian", "bird", "fish", "mammal", "reptile")

img_panel <- function(path) {
  if (!file.exists(path)) stop("Missing image: ", path)
  img <- png::readPNG(path)
  patchwork::wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
}

# --- 1. Combined bivariate panel ---------------------------------------------
cat("Building combined bivariate panel...\n")

biv_panels <- list()
for (taxon in taxa_order) {
  biv_panels[[length(biv_panels) + 1]] <- img_panel(
    file.path(out_dir, paste0(taxon, "_transformed_bivariate_logging_vs_plantation.png"))
  )
  biv_panels[[length(biv_panels) + 1]] <- img_panel(
    file.path(out_dir, paste0(taxon, "_transformed_ses_bivariate_logging_vs_plantation.png"))
  )
}

biv_grid <- patchwork::wrap_plots(biv_panels, ncol = 2, nrow = 5) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 16, face = "bold")
  )

legend_panel <- img_panel(file.path(out_dir, "bivariate_legend.png"))

biv_combined <- biv_grid +
  patchwork::inset_element(
    legend_panel,
    left = 0.80, bottom = 0.0, right = 1.0, top = 0.16,
    align_to = "full", ignore_tag = TRUE
  )

ggsave(file.path(out_dir, "combined_bivariate_panel.png"),
       biv_combined, width = 12, height = 15, dpi = 300, bg = "white")

# --- 2. Combined gradient panel -----------------------------------------------
cat("Building combined gradient panel...\n")

grad_panels <- list()
for (taxon in taxa_order) {
  grad_panels[[length(grad_panels) + 1]] <- img_panel(
    file.path(out_dir, paste0(taxon, "_transformed_gradient_wood_production.png"))
  )
  grad_panels[[length(grad_panels) + 1]] <- img_panel(
    file.path(out_dir, paste0(taxon, "_transformed_ses_wood_production.png"))
  )
}

grad_grid <- patchwork::wrap_plots(grad_panels, ncol = 2, nrow = 5) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 16, face = "bold")
  )

ggsave(file.path(out_dir, "combined_gradient_panel.png"),
       grad_grid, width = 12, height = 15, dpi = 300, bg = "white")

cat("\nCombined panel figures saved to:", out_dir, "\n")
