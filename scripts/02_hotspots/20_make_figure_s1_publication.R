# =============================================================================
#
# FIGURE S1 — PUBLICATION VERSION (Biological Conservation)
#
#
# Description:
# Range-weighted version of Figure 2, using the same map engines
# (publication_map_helpers.R) so the supplementary figure matches the main
# text exactly. Rows are taxa (amphibian, bird, fish, mammal, reptile),
# tagged A–J.
#
#   left  = inverse cube-root range-weighted threat score,
#           T_c = sum over threatened species s of (a_{s,c} / A_s)^(1/3)
#   right = SES of that score from the weighted null model
#           (13_weighted_null_model_ses.R)
#
# Both columns come from the same cube-root weighted cell x species matrices,
# via weighted_threat_raster() and weighted_ses_raster() in the helpers
# (12_build_weighted_transformed_matrices.R, patched by
# 12d_recover_missing_species.R). The left column is therefore computed from
# the matrices rather than the cached *_combined_transformed_threat_raster.tif
# files, which were built with terra::rasterize(cover = TRUE) and silently
# omit species whose whole range is smaller than one grid cell.
#
# Output (outputs/publication_figures/Figure_S1/):
#   - Figure_S1.pdf   fonts embedded (cairo_pdf)   <- submit this
#   - Figure_S1.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_S1_preview.png  300-dpi preview for checking only
#
# =============================================================================

source(here::here("scripts", "02_hotspots", "publication_map_helpers.R"))

s1_left  <- list()
s1_right <- list()
for (taxon in taxa_order) {
  cat("Figure S1:", taxon, "\n")
  s1_left[[taxon]] <- plot_gradient_map(
    weighted_threat_raster(taxon), "Weighted Threat\n(Cube-Root)", label_digits = 2
  )
  s1_right[[taxon]] <- plot_ses_map(weighted_ses_raster(taxon))
}

# Wider legend title than Figure 2 ("Weighted Threat"), so slightly narrower maps
save_map_pair_figure(s1_left, s1_right, "Figure_S1", map_w = 80)
