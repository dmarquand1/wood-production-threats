rm(list = ls())
gc()

# =============================================================================
# COMBINED PANEL FIGURE: MAIN-TEXT RICHNESS MAPS (5 taxa x raw/SES)
#
# Description:
# Composites the existing main-text richness-based maps (built by
# 02_make_taxon_maps.R and 06_null_model_hotspots*.R) into a single 5-row x
# 2-column panel figure for the thesis main body (row order: amphibian, bird,
# fish, mammal, reptile; left column = raw threatened-richness gradient, right
# column = SES(threatened richness)), tagged A-J, no spacing between panels.
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

maps_dir          <- here::here("outputs", "hotspots", "figures", "maps")
null_model_dir    <- here::here("outputs", "hotspots", "figures", "null_model")
ses_diag_dir      <- here::here("outputs", "hotspots", "figures", "ses_diagnostics")
out_dir           <- here::here("outputs", "hotspots", "figures", "main_text_panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

taxa_order <- c("amphibian", "bird", "fish", "mammal", "reptile")

img_panel <- function(path) {
  if (!file.exists(path)) stop("Missing image: ", path)
  img <- png::readPNG(path)
  patchwork::wrap_elements(full = grid::rasterGrob(img, interpolate = TRUE))
}

cat("Building combined richness panel...\n")

panels <- list()
for (taxon in taxa_order) {
  panels[[length(panels) + 1]] <- img_panel(
    file.path(maps_dir, paste0(taxon, "_gradient_threatened.png"))
  )
  panels[[length(panels) + 1]] <- img_panel(
    file.path(null_model_dir, paste0(taxon, "_has_wood_prod_ongoing_ses_map.png"))
  )
}

richness_grid <- patchwork::wrap_plots(panels, ncol = 2, nrow = 5) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 16, face = "bold")
  )

ggsave(file.path(out_dir, "combined_richness_gradient_panel.png"),
       richness_grid, width = 12, height = 15, dpi = 300, bg = "white")

cat("\nSaved:", file.path(out_dir, "combined_richness_gradient_panel.png"), "\n")

# --- Combined bivariate panel (raw richness | SES, no legend) ----------------
cat("\nBuilding combined richness bivariate panel...\n")

biv_panels <- list()
for (taxon in taxa_order) {
  biv_panels[[length(biv_panels) + 1]] <- img_panel(
    file.path(ses_diag_dir, paste0(taxon, "_richness_bivariate_continuous_gradient.png"))
  )
  biv_panels[[length(biv_panels) + 1]] <- img_panel(
    file.path(ses_diag_dir, paste0(taxon, "_ses_bivariate_continuous_gradient.png"))
  )
}

biv_grid <- patchwork::wrap_plots(biv_panels, ncol = 2, nrow = 5) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 16, face = "bold")
  )

ggsave(file.path(out_dir, "combined_richness_bivariate_panel.png"),
       biv_grid, width = 12, height = 15, dpi = 300, bg = "white")

cat("\nSaved:", file.path(out_dir, "combined_richness_bivariate_panel.png"), "\n")
