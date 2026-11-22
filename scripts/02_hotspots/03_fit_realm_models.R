# =============================================================================
# SCRIPT: VERTEBRATE THREAT MODELING - INTERACTION ANALYSIS (BINOMIAL)
# =============================================================================

# --- 1. SETUP & LIBRARIES ---
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, tidyr, ggplot2, here, sf, terra, rnaturalearth, 
               broom, purrr, stringr, rredlist, cowplot) # Removed MASS

# --- 2. GLOBAL CONSTANTS ---
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
global_extent <- ext(-18040000, 18040000, -9020000, 9020000)
grid_template <- rast(resolution = 110000, crs = moll_crs, extent = global_extent, vals = 0)
realm_axis_title_size <- 25
realm_axis_text_size <- 20

# Realm Colors
my_colors <- c(
  "Afrotropical"   = "#FF6B6B", "Australian"     = "#D4A017",  
  "Madagascan"     = "#FFD54F", "Nearctic"       = "#55A630",
  "Neotropical"    = "#00C853", "Oceanina"       = "#00BFA5",  
  "Oriental"       = "#00B0FF", "Palearctic"     = "#2979FF",
  "Panamanian"     = "#6200EA", "Saharo-Arabian" = "#D500F9",  
  "Sino-Japanese"  = "#C51162"   
)

# Output Directories
table_dir <- here::here("outputs", "hotspots", "tables")
figure_dir <- here::here("outputs", "hotspots", "figures", "realm")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Canonical raster directory
raster_dir <- here::here("data", "processed", "hotspots", "rasters")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)
raster_path <- function(filename) file.path(raster_dir, filename)

# --- 3. LOAD DATA ---
cat("Loading World and Realm Data...\n")
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica") %>% st_transform(moll_crs)

path_shp <- here::here("data", "raw", "realms", "newRealms.shp")
path_rds <- here::here("data", "realms_global_sf.rds")

if (file.exists(path_shp)) {
  suppressMessages(sf::sf_use_s2(FALSE))
  realms_sf <- st_read(path_shp, quiet = TRUE)
} else if (file.exists(path_rds)) {
  realms_sf <- readRDS(path_rds)
} else {
  stop("Realm data not found.")
}

realms_sf <- st_make_valid(realms_sf) %>% st_wrap_dateline(options = c("WRAPDATELINE=YES")) %>% st_make_valid()
potential_cols <- names(realms_sf)[sapply(realms_sf, function(x) length(unique(x)) > 4 & length(unique(x)) < 20)]
realm_col_name <- grep("Realm", potential_cols, value = TRUE, ignore.case = TRUE)[1]
if (is.na(realm_col_name)) realm_col_name <- potential_cols[1]

# --- 4. HELPER FUNCTIONS ---

prepare_model_data <- function(r_total, r_threatened, realms_sf) {
  rx_masked <- terra::mask(r_total, vect(world_map))
  ry_masked <- terra::mask(r_threatened, vect(world_map))
  
  realms_proj <- st_transform(realms_sf, crs(rx_masked))
  realms_proj$rid <- as.numeric(as.factor(realms_proj[[realm_col_name]]))
  realm_r <- terra::rasterize(vect(realms_proj), rx_masked, field = "rid", touches = TRUE)
  
  target_mask <- rx_masked; target_mask[!is.na(target_mask)] <- 1 
  for(i in 1:15) { 
    missing_mask <- !is.na(target_mask) & is.na(realm_r)
    if (global(missing_mask, "sum", na.rm=TRUE)[1,1] == 0) break
    filled <- terra::focal(realm_r, w=5, fun="modal", na.policy="only", na.rm=TRUE)
    realm_r <- terra::cover(realm_r, filled)
  }
  
  stack_data <- c(rx_masked, ry_masked, realm_r)
  names(stack_data) <- c("Total", "Threatened", "realm_id")
  
  df <- as.data.frame(stack_data, xy = FALSE) %>% 
    filter(!is.na(Total) & !is.na(Threatened) & !is.na(realm_id)) %>%
    filter(Total > 0)
  
  lookup <- unique(st_drop_geometry(realms_proj)[, c("rid", realm_col_name)])
  colnames(lookup) <- c("realm_id", "Realm_Name")
  
  df <- df %>% left_join(lookup, by = "realm_id") %>% 
    filter(!is.na(Realm_Name)) %>%
    mutate(Realm_Name = as.factor(Realm_Name))
  
  return(df)
}

# --- 5. EXECUTION LOOP ---

taxa_list <- c("bird", "mammal", "amphibian", "reptile", "fish")
model_results <- list()

for (taxon in taxa_list) {
  cat("\n========================================\n")
  cat("MODELING:", toupper(taxon), "\n")
  
  # 1. Load Rasters
  path_pot   <- raster_path(paste0(taxon, "_potential_richness_raster.tif"))
  path_comb  <- raster_path(paste0(taxon, "_combined_richness_raster.tif"))
  
  if (!file.exists(path_pot)) { cat("  [!] Raster not found, skipping.\n"); next }
  
  r_pot <- rast(path_pot)
  r_comb <- rast(path_comb)
  
  # 2. Prepare Data
  cat("  -> Extracting spatial data...\n")
  mod_df <- prepare_model_data(r_pot, r_comb, realms_sf) 
  
  # --- CRITICAL STEP FOR BINOMIAL MODELS ---
  # We need the count of "Failures" (Non-Threatened species)
  mod_df$Not_Threatened <- mod_df$Total - mod_df$Threatened
  mod_df <- mod_df %>% filter(Not_Threatened >= 0)
  
  # 3. Run Binomial GLM
  cat("  -> Running Binomial GLM...\n")
  
  # Note: family = binomial(link = "logit")
  # Response is a matrix of cbind(Successes, Failures)
  mod <- glm(cbind(Threatened, Not_Threatened) ~ Total * Realm_Name, 
             family = binomial(link = "logit"), 
             data = mod_df)
  
  # 4. Save Summary
  tidy_mod <- tidy(mod) %>% mutate(Taxon = taxon)
  glance_mod <- glance(mod) %>% mutate(Taxon = taxon)
  
  # Store results
  # We store 'mod' so we can use it for predictions later
  model_results[[taxon]] <- list(tidy = tidy_mod, glance = glance_mod, model_obj = mod)
  
  # Save individual CSV
  write.csv(tidy_mod, file.path(table_dir, paste0(taxon, "_realm_coefficients.csv")), row.names = FALSE)
}

# --- 6. VISUALIZATION & SUMMARY ---

if(length(model_results) > 0) {
  
  # Combine all tidy results
  all_coefs <- bind_rows(lapply(model_results, function(x) x$tidy))
  
  cat("\n========================================\n")
  cat("GENERATING INTERACTION PLOT\n")
  
  # Filter for Interaction Terms
  interaction_plot_data <- all_coefs %>%
    filter(grepl(":", term)) %>%
    mutate(
      Realm = str_remove(term, "Total:Realm_Name"),
      # Clean up the term name for better plotting
      Realm = str_replace_all(Realm, "_", " ")
    )
  
  # Plot: Interaction Estimates
  p_effects <- ggplot(interaction_plot_data, aes(x = estimate, y = Realm, color = Realm)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = estimate - std.error, xmax = estimate + std.error), height = 0.2) +
    facet_wrap(~Taxon, scales = "free_x") +
    scale_color_manual(values = my_colors) +
    labs(
      title = "Difference in Threat Probability Slope by Realm",
      subtitle = "Interaction Terms (Total Richness * Realm)\nValues > 0 indicate the probability of being threatened increases FASTER with richness than in the baseline.",
      x = "Interaction Estimate (Difference in Log-Odds Slope)",
      y = "Realm"
    ) +
    theme_cowplot() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = realm_axis_title_size),
      axis.text = element_text(size = realm_axis_text_size),
      panel.grid.major.x = element_line(color = "grey90")
    )
  
  ggsave(file.path(figure_dir, "interaction_plot.png"), p_effects, width = 12, height = 8, bg = "white")
  
  cat("\nDone! Results saved to:", figure_dir, "and", table_dir, "\n")
}

# =============================================================================
# PART C: VISUALIZATION - BINOMIAL TREND LINES WITH 95% CI
# =============================================================================
cat("\n========================================\n")
cat("GENERATING TREND LINES WITH 95% CI RIBBONS\n")

plot_dir <- figure_dir

for (taxon in taxa_list) {
  
  cat("  -> Processing:", toupper(taxon), "\n")
  
  # 1. Load Data
  path_pot  <- raster_path(paste0(taxon, "_potential_richness_raster.tif"))
  path_comb <- raster_path(paste0(taxon, "_combined_richness_raster.tif"))
  
  if(!file.exists(path_pot)) next
  
  r_pot <- rast(path_pot); r_comb <- rast(path_comb)
  plot_df <- prepare_model_data(r_pot, r_comb, realms_sf)
  
  # 2. Run Binomial Model
  # We model the probability of a species being threatened
  plot_df$Not_Threatened <- plot_df$Total - plot_df$Threatened
  plot_df <- plot_df %>% filter(Not_Threatened >= 0)
  
  cat("     Running Binomial GLM...\n")
  mod_obj <- glm(cbind(Threatened, Not_Threatened) ~ Total * Realm_Name, 
                 family = binomial(link = "logit"), 
                 data = plot_df)
  
  # 3. Generate Predictions with Standard Error
  cat("     Calculating Confidence Intervals...\n")
  
  # Create a smooth range of Total richness for each Realm
  pred_grid <- plot_df %>%
    dplyr::group_by(Realm_Name) %>%
    dplyr::summarise(min_t = min(Total), max_t = max(Total)) %>%
    dplyr::mutate(data = purrr::map2(min_t, max_t, ~data.frame(Total = seq(.x, .y, length.out = 100)))) %>%
    tidyr::unnest(data) %>%
    dplyr::select(Realm_Name, Total)
  
  # Predict on the LINK scale (logit) to get proper SEs
  preds <- predict(mod_obj, newdata = pred_grid, type = "link", se.fit = TRUE)
  
  # Calculate CI on logit scale
  pred_grid$fit_link <- preds$fit
  pred_grid$se_link  <- preds$se.fit
  pred_grid$upper_link <- preds$fit + (1.96 * preds$se.fit)
  pred_grid$lower_link <- preds$fit - (1.96 * preds$se.fit)
  
  # Transform back to Probability (0-1) using plogis (inverse logit)
  pred_grid$Prob_Fit   <- plogis(pred_grid$fit_link)
  pred_grid$Prob_Upper <- plogis(pred_grid$upper_link)
  pred_grid$Prob_Lower <- plogis(pred_grid$lower_link)
  
  # Scale up to Counts (Probability * Total Species)
  pred_grid$Count_Fit   <- pred_grid$Prob_Fit * pred_grid$Total
  pred_grid$Count_Upper <- pred_grid$Prob_Upper * pred_grid$Total
  pred_grid$Count_Lower <- pred_grid$Prob_Lower * pred_grid$Total
  
  # 4. Plot
  cat("     Generating plot...\n")
  background_points <- plot_df %>% sample_n(min(nrow(plot_df), 5000))
  
  p_trend <- ggplot() +
    # A. Background Points (Faint)
    geom_point(data = background_points, aes(x = Total, y = Threatened, color = Realm_Name), 
               alpha = 0.25, size = 1.5, shape = 16) +
    
    # B. Confidence Ribbons (First layer so lines go on top)
    geom_ribbon(data = pred_grid, aes(x = Total, ymin = Count_Lower, ymax = Count_Upper, fill = Realm_Name), 
                alpha = 0.2) +
    
    # C. Trend Lines
    geom_line(data = pred_grid, aes(x = Total, y = Count_Fit, color = Realm_Name), 
              linewidth = 1.5) +
    
    # D. Styling
    scale_color_manual(values = my_colors, name = "Zoogeographic Realm") +
    scale_fill_manual(values = my_colors, name = "Zoogeographic Realm") +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    labs(
      title = paste0(tools::toTitleCase(taxon), ": Threat Accumulation Trends"),
      subtitle = "Model: Binomial GLM with 95% Confidence Intervals",
      x = "Total Species Richness",
      y = "Threatened Species Richness"
    ) +
    theme_cowplot() +
    theme(
      aspect.ratio = 1,
      legend.position = "right",
      panel.grid.major = element_blank(),
      axis.title = element_text(size = realm_axis_title_size),
      axis.text = element_text(size = realm_axis_text_size)
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, linewidth = 2)))
  
  # 5. Save
  ggsave(file.path(plot_dir, paste0(taxon, "_realm_model_trend.png")), p_trend, width = 10, height = 8, bg = "white")
}
cat("\nAll plots with ribbons saved.\n")

