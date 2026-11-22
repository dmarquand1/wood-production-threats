rm(list = ls())

# ==============================================================================
# AGGREGATE BAYESIAN R2 PARTITION ACROSS ALL MAIN TAXA
# Reads outputs/traits/diagnostics/<taxon>/<model>_r2_bayes_partition.csv for
# each taxon/model that has one and combines into one file, replacing the
# previous hand-assembled r2_bayes_partition_all_taxa.csv.
# ==============================================================================

taxa <- c("amphibian", "bird", "fish", "mammal", "reptile")
models <- c("logging", "plantation")
diag_dir <- file.path("outputs", "traits", "diagnostics")

rows <- list()
missing <- character(0)
for (taxon in taxa) {
  for (model_name in models) {
    path <- file.path(diag_dir, taxon, sprintf("%s_r2_bayes_partition.csv", model_name))
    if (file.exists(path)) {
      rows[[length(rows) + 1]] <- read.csv(path, stringsAsFactors = FALSE)
    } else {
      missing <- c(missing, sprintf("%s/%s", taxon, model_name))
    }
  }
}

combined <- do.call(rbind, rows)
out_path <- file.path(diag_dir, "r2_bayes_partition_all_taxa.csv")
write.csv(combined, out_path, row.names = FALSE)

cat("Wrote", out_path, "\n")
if (length(missing) > 0) cat("Missing (not included):", paste(missing, collapse = ", "), "\n")
