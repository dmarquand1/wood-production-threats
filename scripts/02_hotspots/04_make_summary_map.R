rm(list = ls())
gc()

# ==============================================================================
# PREPARE THE HOTSPOT SUMMARY-FIGURE LAYER
# ==============================================================================

summary_dir <- file.path("outputs", "hotspots", "figures", "summary")
table_dir <- file.path("outputs", "hotspots", "tables")

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

summary_panels <- data.frame(
  taxon = c("all_vertebrates", "amphibian", "bird", "fish", "mammal", "reptile"),
  panel_file = c(
    file.path(summary_dir, "figure2_all_vertebrates_panel.png"),
    file.path(summary_dir, "figure2_amphibian_panel.png"),
    file.path(summary_dir, "figure2_bird_panel.png"),
    file.path(summary_dir, "figure2_fish_panel.png"),
    file.path(summary_dir, "figure2_mammal_panel.png"),
    file.path(summary_dir, "figure2_reptile_panel.png")
  ),
  stringsAsFactors = FALSE
)
summary_panels$exists <- file.exists(summary_panels$panel_file)

write.csv(
  summary_panels,
  file = file.path(table_dir, "figure2_summary_panel_manifest.csv"),
  row.names = FALSE
)

missing_panels <- summary_panels[!summary_panels$exists, , drop = FALSE]
if (nrow(missing_panels) > 0) {
  stop(
    "Missing expected summary panels: ",
    paste(missing_panels$taxon, collapse = ", ")
  )
}

source_panel <- file.path(summary_dir, "figure2_all_vertebrates_panel.png")
target_panel <- file.path(summary_dir, "figure2_summary_map.png")

file.copy(source_panel, target_panel, overwrite = TRUE)

message("Summary panel manifest written to ", file.path(table_dir, "figure2_summary_panel_manifest.csv"))
message("Canonical summary-map alias refreshed at ", target_panel)
