rm(list = ls())
gc()

# ==============================================================================
# Combined coefficient forest plot for five vertebrate taxa
#
# Uses the saved fixed-effect summary tables from each taxon's summarise script.
# Panel heights are proportional to the number of plotted trait rows, with a
# small fixed allowance for titles and axes, so trait spacing stays consistent.
# ==============================================================================

project_lib <- file.path(getwd(), ".r-lib")
if (dir.exists(project_lib)) {
  .libPaths(c(normalizePath(project_lib, winslash = "/", mustWork = FALSE), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

diagnostic_root <- file.path("outputs", "traits", "diagnostics")
figure_dir <- file.path("outputs", "traits", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

taxa <- tibble::tribble(
  ~taxon,       ~label,
  "amphibian",  "Amphibians",
  "bird",       "Birds",
  "fish",       "Fish",
  "mammal",     "Mammals",
  "reptile",    "Reptiles"
)

threat_colors <- c(
  "ThreatLogging"    = "#d9be00",
  "ThreatPlantation" = "#0088d9"
)

threat_labels <- c(
  "ThreatLogging"    = "Natural Forestry",
  "ThreatPlantation" = "Plantation"
)

term_labels <- c(
  Body_Size_log = "Body Size",
  Body_Mass_log = "Body Mass",
  Body_Mass_Log = "Body Mass",
  Max_Body_Mass_log = "Body Mass",
  MaxLength_Log = "Maximum Length",
  log_range_size = "Range Size",
  habitat_breadth_log = "Habitat Breadth",
  Habitat_Breadth_log = "Habitat Breadth",
  Habitat_Breadth = "Habitat Breadth",
  Salinity_Breadth_sc = "Salinity Breadth",
  Trophic_Level_sc = "Trophic Level",
  PredPreyRatioMax_log_sc = "Predator:Prey Ratio",
  Clutch_Size_log = "Clutch Size",
  HWI_log = "Hand-Wing Index",
  Litter_Size_Log = "Litter Size",
  AirBreathing_Bin = "Air Breathing",
  Body_Size_log.log_range_size = "Body Size x Range Size",
  Body_Mass_log.log_range_size = "Body Mass x Range Size",
  Body_Mass_Log.log_range_size = "Body Mass x Range Size",
  Max_Body_Mass_log.log_range_size = "Body Mass x Range Size",
  MaxLength_Log.log_range_size = "Maximum Length x Range Size",
  Breeding_CatDir = "Direct Development",
  Breeding_CatViv = "Viviparous",
  Active.timeCrepuscular_Cathemeral = "Crepuscular/Cathemeral",
  Active.timeCrepuscular_Cathermal = "Crepuscular/Cathemeral",
  Active.timeDiurnal = "Diurnal",
  Active.timeNocturnal = "Nocturnal",
  Reproductive.modeviviparous_ovoviviparous = "Viviparous/Ovoviviparous",
  FangsfrontMfanged = "Front-fanged",
  FangsrearMfanged = "Rear-fanged",
  FangsnonMfanged = "Non-fanged",
  Forest = "Forest",
  Savanna = "Savanna",
  Shrubland = "Shrubland",
  Grassland = "Grassland",
  Wetlands = "Wetlands",
  Diet_Invertebrate = "Diet: Invertebrate",
  Diet_Fruit = "Diet: Fruit",
  Diet_Nectar = "Diet: Nectar",
  Diet_Seed = "Diet: Seed",
  Diet_Vertebrate = "Diet: Vertebrate",
  Diet_Fish = "Diet: Fish",
  Diet_Scavenger = "Diet: Scavenger",
  Diet_PlantMaterial = "Diet: Plant Material",
  binary_det_inv = "Invertebrates",
  binary_det_vend = "Endotherm Vertebrates",
  binary_det_vect = "Ectotherm Vertebrates",
  binary_det_vunk = "Unknown Vertebrates",
  binary_det_scav = "Scavenger",
  binary_det_fruit = "Fruit",
  binary_det_nect = "Nectar",
  binary_det_seed = "Seed",
  binary_det_plantother = "Other Plant Material",
  activity_cycleCathemeral_Crepuscular = "Cathemeral/Crepuscular",
  activity_cycleDiurnal = "Diurnal",
  DemersPelag_Catdemersal = "Demersal/Pelagic: demersal",
  DemersPelag_Catpelagic = "Demersal/Pelagic: pelagic"
)

clean_term_label <- function(term) {
  lookup_key <- gsub(":", ".", term, fixed = TRUE)
  labels <- unname(term_labels[lookup_key])

  missing_labels <- is.na(labels)
  if (any(missing_labels)) {
    fallback <- term[missing_labels]
    fallback <- gsub("_", " ", fallback, fixed = TRUE)
    fallback <- gsub("\\.", " ", fallback)
    fallback <- gsub("([a-z])([A-Z])", "\\1 \\2", fallback)
    fallback <- gsub(" M ", "-", fallback)
    labels[missing_labels] <- fallback
  }

  labels
}

harmonise_reptile_factor_draws <- function(plot_data) {
  factor_labels <- c(
    "Crepuscular/Cathemeral", "Diurnal", "Nocturnal",
    "Front-fanged", "Non-fanged", "Rear-fanged"
  )
  
  if (!any(plot_data$Clean_Label == "Non-fanged") &&
      !any(plot_data$Clean_Label == "Diurnal")) {
    return(plot_data)
  }
  
  indexed <- plot_data %>%
    group_by(Response, Clean_Label) %>%
    mutate(.draw = row_number()) %>%
    ungroup()
  
  get_draws <- function(response_label, clean_label) {
    indexed %>%
      filter(Response == response_label, Clean_Label == clean_label) %>%
      select(.draw, .value)
  }
  
  make_single_draw <- function(draws, response_label, clean_label, term_name) {
    draws %>%
      transmute(
        term = term_name,
        .value = .value,
        Clean_Label = clean_label,
        Response = response_label
      )
  }
  
  logging_activity <- indexed %>%
    filter(
      Response == "ThreatLogging",
      Clean_Label %in% c("Crepuscular/Cathemeral", "Nocturnal")
    ) %>%
    select(term, .value, Clean_Label, Response)
  
  logging_fangs <- indexed %>%
    filter(
      Response == "ThreatLogging",
      Clean_Label %in% c("Front-fanged", "Rear-fanged")
    ) %>%
    select(term, .value, Clean_Label, Response)
  
  plantation_diurnal <- get_draws("ThreatPlantation", "Diurnal")
  plantation_nocturnal <- get_draws("ThreatPlantation", "Nocturnal")
  plantation_non_fanged <- get_draws("ThreatPlantation", "Non-fanged")
  plantation_rear_fanged <- get_draws("ThreatPlantation", "Rear-fanged")
  
  plantation_activity <- bind_rows(
    make_single_draw(
      plantation_diurnal %>% mutate(.value = -.value),
      "ThreatPlantation",
      "Crepuscular/Cathemeral",
      "Active.timeCrepuscular_Cathemeral_common_ref_Diurnal"
    ),
    plantation_nocturnal %>%
      inner_join(plantation_diurnal, by = ".draw", suffix = c("_nocturnal", "_diurnal")) %>%
      transmute(.draw, .value = .value_nocturnal - .value_diurnal) %>%
      make_single_draw(
        "ThreatPlantation",
        "Nocturnal",
        "Active.timeNocturnal_common_ref_Diurnal"
      )
  )
  
  plantation_fangs <- bind_rows(
    make_single_draw(
      plantation_non_fanged %>% mutate(.value = -.value),
      "ThreatPlantation",
      "Front-fanged",
      "FangsfrontMfanged_common_ref_nonMfanged"
    ),
    plantation_rear_fanged %>%
      inner_join(plantation_non_fanged, by = ".draw", suffix = c("_rear", "_non")) %>%
      transmute(.draw, .value = .value_rear - .value_non) %>%
      make_single_draw(
        "ThreatPlantation",
        "Rear-fanged",
        "FangsrearMfanged_common_ref_nonMfanged"
      )
  )
  
  plot_data %>%
    filter(!Clean_Label %in% factor_labels) %>%
    bind_rows(logging_activity, plantation_activity, logging_fangs, plantation_fangs)
}

read_reptile_draw_coefficients <- function(label, model_name, response_label) {
  draw_path <- file.path("outputs", "traits", "models", "reptile", "hpc_coefficient_draws.rds")
  if (!file.exists(draw_path)) {
    return(NULL)
  }
  
  readRDS(draw_path) %>%
    harmonise_reptile_factor_draws() %>%
    filter(Response == response_label) %>%
    group_by(Response, term, Clean_Label) %>%
    summarise(
      estimate = median(.value),
      lower = quantile(.value, 0.025),
      upper = quantile(.value, 0.975),
      .groups = "drop"
    ) %>%
    transmute(
      Taxon = label,
      Model = model_name,
      Response = response_label,
      term = term,
      Clean_Label = Clean_Label,
      estimate = estimate,
      lower = lower,
      upper = upper
    )
}

read_taxon_coefficients <- function(taxon, label, model_name, response_label) {
  if (taxon == "reptile") {
    reptile_coefficients <- read_reptile_draw_coefficients(label, model_name, response_label)
    if (!is.null(reptile_coefficients)) {
      return(reptile_coefficients)
    }
  }
  
  coefficient_path <- file.path(
    diagnostic_root,
    taxon,
    sprintf("%s_fixed_effects_summary.csv", model_name)
  )
  if (!file.exists(coefficient_path)) {
    stop("Missing fixed-effect summary: ", coefficient_path)
  }

  read.csv(coefficient_path, stringsAsFactors = FALSE) %>%
    filter(term != "Intercept") %>%
    transmute(
      Taxon = label,
      Model = model_name,
      Response = response_label,
      term = term,
      Clean_Label = clean_term_label(term),
      estimate = estimate,
      lower = q2.5,
      upper = q97.5
    )
}

plot_summary <- bind_rows(lapply(seq_len(nrow(taxa)), function(i) {
  taxon <- taxa$taxon[[i]]
  label <- taxa$label[[i]]

  bind_rows(
    read_taxon_coefficients(taxon, label, "logging", "ThreatLogging"),
    read_taxon_coefficients(taxon, label, "plantation", "ThreatPlantation")
  )
}))

term_order <- plot_summary %>%
  group_by(Taxon, Clean_Label) %>%
  summarise(Sort_Value = mean(estimate), .groups = "drop") %>%
  arrange(Taxon, Sort_Value) %>%
  group_by(Taxon) %>%
  mutate(
    y = row_number(),
    trait_count = n()
  ) %>%
  ungroup()

plot_summary <- plot_summary %>%
  left_join(term_order, by = c("Taxon", "Clean_Label")) %>%
  mutate(
    Response = factor(Response, levels = names(threat_colors)),
    response_offset = case_when(
      Response == "ThreatLogging" ~ -0.18,
      Response == "ThreatPlantation" ~ 0.18,
      TRUE ~ 0
    ),
    y_plot = y + response_offset
  )

write.csv(
  plot_summary %>%
    arrange(Taxon, y, Response) %>%
    select(Taxon, Model, Response, term, Clean_Label, estimate, lower, upper, y),
  file.path(figure_dir, "all_taxa_coefficients.csv"),
  row.names = FALSE
)

make_taxon_plot <- function(taxon_label) {
  taxon_summary <- plot_summary %>%
    filter(Taxon == taxon_label)
  
  taxon_order <- term_order %>%
    filter(Taxon == taxon_label)
  
  x_limits <- range(c(taxon_summary$lower, taxon_summary$upper), na.rm = TRUE)
  if (!all(is.finite(x_limits)) || diff(x_limits) == 0) {
    x_limits <- x_limits + c(-0.5, 0.5)
  } else {
    x_padding <- diff(x_limits) * 0.05
    x_limits <- x_limits + c(-x_padding, x_padding)
  }
  
  ggplot(taxon_summary, aes(color = Response)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    geom_segment(
      aes(x = lower, xend = upper, y = y_plot, yend = y_plot),
      linewidth = 1.2,
      lineend = "round"
    ) +
    geom_point(
      aes(x = estimate, y = y_plot),
      size = 3
    ) +
    scale_color_manual(
      values = threat_colors,
      labels = threat_labels,
      name = "Production Pathway"
    ) +
    scale_y_continuous(
      breaks = taxon_order$y,
      labels = taxon_order$Clean_Label,
      limits = c(0.4, max(taxon_order$y) + 0.6),
      expand = expansion(mult = 0)
    ) +
    scale_x_continuous(expand = expansion(mult = 0)) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(
      title = NULL,
      x = "Effect on risk (log-odds)",
      y = NULL
    ) +
    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(size = 23, face = "bold", hjust = 0),
      axis.title.x = element_text(size = 20, margin = margin(t = 10)),
      axis.text.x = element_text(size = 17),
      axis.text.y = element_text(size = 17),
      legend.title = element_text(size = 19),
      legend.text = element_text(size = 18),
      legend.position = "bottom",
      panel.grid.major.y = element_line(color = "gray95", linetype = "dotted"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 18, 8, 8)
    )
}

plots <- setNames(
  lapply(taxa$label, make_taxon_plot),
  taxa$label
)

trait_counts <- term_order %>%
  distinct(Taxon, trait_count) %>%
  tibble::deframe()

fixed_plot_rows <- 6
area_rows <- trait_counts + fixed_plot_rows

legend_rows <- 5
panel_rows <- max(
  area_rows["Birds"] + area_rows["Mammals"],
  area_rows["Amphibians"] + area_rows["Fish"] + area_rows["Reptiles"]
)
layout_rows <- panel_rows + legend_rows

left_rows <- area_rows["Birds"] + area_rows["Mammals"]
right_rows <- area_rows["Amphibians"] + area_rows["Fish"] + area_rows["Reptiles"]

left_gap <- panel_rows - left_rows
right_gap <- panel_rows - right_rows

left_bird_top <- 1
left_bird_bottom <- area_rows["Birds"]
left_mammal_top <- left_bird_bottom + 1
left_mammal_bottom <- left_mammal_top + area_rows["Mammals"] + left_gap - 1

right_amphibian_top <- 1
right_amphibian_bottom <- area_rows["Amphibians"]
right_fish_top <- right_amphibian_bottom + 1
right_fish_bottom <- right_fish_top + area_rows["Fish"] - 1
right_reptile_top <- right_fish_bottom + 1
right_reptile_bottom <- right_reptile_top + area_rows["Reptiles"] + right_gap - 1
legend_top <- panel_rows + 1

layout_design <- c(
  patchwork::area(
    t = left_bird_top,
    l = 1,
    b = left_bird_bottom,
    r = 1
  ),
  patchwork::area(
    t = left_mammal_top,
    l = 1,
    b = left_mammal_bottom,
    r = 1
  ),
  patchwork::area(
    t = right_amphibian_top,
    l = 2,
    b = right_amphibian_bottom,
    r = 2
  ),
  patchwork::area(
    t = right_fish_top,
    l = 2,
    b = right_fish_bottom,
    r = 2
  ),
  patchwork::area(
    t = right_reptile_top,
    l = 2,
    b = right_reptile_bottom,
    r = 2
  ),
  patchwork::area(
    t = legend_top,
    l = 1,
    b = layout_rows,
    r = 2
  )
)

combined_plot <-
  plots$Birds + plots$Mammals +
  plots$Amphibians + plots$Fish + plots$Reptiles + patchwork::guide_area() +
  plot_layout(
    design = layout_design,
    widths = c(1, 1),
    heights = rep(1, layout_rows),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.key.width = unit(2.2, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.spacing.x = unit(1.0, "cm")
  )

output_png <- file.path(figure_dir, "all_taxa_coefficients.png")
output_pdf <- file.path(figure_dir, "all_taxa_coefficients.pdf")

figure_width <- 16
figure_height <- layout_rows * ((9 - 2) / 21) + 1

ggsave(
  output_png,
  combined_plot,
  width = figure_width,
  height = figure_height,
  dpi = 300
)

ggsave(
  output_pdf,
  combined_plot,
  width = figure_width,
  height = figure_height
)

message("Saved combined coefficient plot to ", output_png)
message("Saved combined coefficient plot to ", output_pdf)
