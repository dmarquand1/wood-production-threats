# =============================================================================
#
# FIGURES 2 & 3 — PUBLICATION VERSIONS (Biological Conservation)
#
#
# Description:
# Rebuilds the main-text hotspot map figures directly from the rasters (rather
# than compositing saved PNGs as 18_make_richness_combined_panel.R does), at
# final print size for Elsevier / Biological Conservation. Rows are taxa
# (amphibian, bird, fish, mammal, reptile), tagged A–J.
#
#   Figure 2: left  = threatened species richness (wood production),
#                     magma, sqrt-distributed colours (plot_gradient_sqrt()
#                     in 02_make_taxon_maps.R); cells with species but no
#                     threatened species are black
#             right = SES of threatened richness from the null model,
#                     symmetric blue–white–red (plot_continuous_ses_map() in
#                     08_bivariate_SES_mapping.R)
#   Figure 3: left  = bivariate threatened richness, natural forestry
#                     (logging) vs plantation
#             right = bivariate SES, natural forestry vs plantation
#                     (plot_true_continuous_bivariate() in
#                     08_bivariate_SES_mapping.R), one shared square legend
#
# Output (outputs/publication_figures/Figure_2/ and Figure_3/):
#   - Figure_N.pdf   fonts embedded (cairo_pdf); map cells are embedded
#                    images, all text/lines are vector   <- submit this
#   - Figure_N.tiff  1000-dpi LZW-compressed bitmap (meets every bitmap
#                    requirement in the artwork guide)
#   - Figure_N_preview.png  300-dpi preview for checking only
#
# Sizing follows the Elsevier artwork guide: full-page width = 190 mm; all
# text is 6–8 pt at final size (Arial).
#
# =============================================================================


source(here::here("scripts", "02_hotspots", "publication_map_helpers.R"))


# =============================================================================
# --- FIGURE 2: THREATENED RICHNESS AND SES ---
# =============================================================================
# Left: threatened species richness, with cells that hold no species of the
# taxon left as land background. Right: SES from the richness null model.

fig2_left  <- list()
fig2_right <- list()
for (taxon in taxa_order) {
  cat("Figure 2:", taxon, "\n")
  r_comb <- rast(file.path(raster_dir, paste0(taxon, "_combined_richness_raster.tif")))
  r_pot  <- rast(file.path(raster_dir, paste0(taxon, "_potential_richness_raster.tif")))
  fig2_left[[taxon]] <- plot_gradient_map(terra::ifel(r_pot > 0, r_comb, NA), "Threatened\nRichness")

  r_ses <- rast(file.path(ses_dir, paste0(taxon, "_has_wood_prod_ongoing_ses_raster.tif")))
  fig2_right[[taxon]] <- plot_ses_map(r_ses)
}

save_map_pair_figure(fig2_left, fig2_right, "Figure_2")


# =============================================================================
# --- FIGURE 3: BIVARIATE NATURAL FORESTRY vs PLANTATION ---
# =============================================================================
# Left: threatened richness for logging vs plantation, from the cached
# presence matrices exactly as 08_bivariate_SES_mapping.R computes it.
# Right: logging vs plantation SES.

fig3_raw <- list()
fig3_ses <- list()
for (taxon in taxa_order) {
  cat("Figure 3:", taxon, "\n")
  presence <- readRDS(file.path(matrix_dir, paste0(taxon, "_cell_species_presence_matrix.rds")))
  fig3_raw[[taxon]] <- list(
    x = matrix_sum_raster(presence, "has_logging_ongoing"),
    y = matrix_sum_raster(presence, "has_plantation_ongoing")
  )
  rm(presence)

  fig3_ses[[taxon]] <- list(
    x = rast(file.path(ses_pathway_dir, paste0(taxon, "_logging_ses.tif"))),
    y = rast(file.path(ses_pathway_dir, paste0(taxon, "_plantation_ses.tif")))
  )
}

save_bivariate_figure(fig3_raw, fig3_ses, "Figure_3", legend_title = "Threatened\nRichness/SES")
