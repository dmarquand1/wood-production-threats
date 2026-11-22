# =============================================================================
#
# FIGURE S3 — PUBLICATION VERSION (Biological Conservation)
#
#
# Description:
# Range-weighted version of Figure 3, using the same bivariate map engine
# (publication_map_helpers.R). Rows are taxa (amphibian, bird, fish, mammal,
# reptile), tagged A–J.
#
#   left  = bivariate cube-root range-weighted threat score, natural forestry
#           (logging) vs plantation
#   right = bivariate SES of those scores (13_weighted_null_model_ses.R)
#
# Each colour axis is scaled to its own panel maximum; SES values <= 0 are
# floored at 0 (grey). Scores come from the patched weighted matrices, as in
# Figures S1 and S2.
#
# Output (outputs/publication_figures/Figure_S3/):
#   - Figure_S3.pdf   fonts embedded (cairo_pdf)   <- submit this
#   - Figure_S3.tiff  1000-dpi LZW-compressed bitmap
#   - Figure_S3_preview.png  300-dpi preview for checking only
#
# =============================================================================

source(here::here("scripts", "02_hotspots", "publication_map_helpers.R"))

s3_raw <- list()
s3_ses <- list()
for (taxon in taxa_order) {
  cat("Figure S3:", taxon, "\n")
  weighted <- read_weighted_matrix(taxon)
  s3_raw[[taxon]] <- list(
    x = matrix_sum_raster(weighted, "has_logging_ongoing"),
    y = matrix_sum_raster(weighted, "has_plantation_ongoing")
  )
  rm(weighted)

  s3_ses[[taxon]] <- list(
    x = weighted_ses_raster(taxon, "has_logging_ongoing"),
    y = weighted_ses_raster(taxon, "has_plantation_ongoing")
  )
}

save_bivariate_figure(s3_raw, s3_ses, "Figure_S3", legend_title = "Inverse Cube-Root\nWeighted Threat")
