library(tidyverse)

# ==========================================
# 1. Load data
# ==========================================
df <- read.csv("data/processed/assessments/all_taxa_assessment_consolidated.csv")

# Ensure taxa_group doesn't have missing or empty strings
df <- df %>% filter(!is.na(taxa_group) & taxa_group != "")

# ==========================================
# 2. Define the Lookup Table for Threat Names
# ==========================================
lookup_raw <- tibble(text = c(
  "1_1 Housing & urban areas", "1_2 Commercial & industrial areas", "1_3 Tourism & recreation areas",
  "2_1 Annual & perennial non-timber crops", "2_2 Wood & pulp plantations", "2_3 Livestock farming & ranching", "2_4 Marine & freshwater aquaculture",
  "3_1 Oil & gas drilling", "3_2 Mining & quarrying", "3_3 Renewable energy",
  "4_1 Roads & railroads", "4_2 Utility & service lines", "4_3 Shipping lanes", "4_4 Flight paths",
  "5_1 Hunting & collecting terrestrial animals", "5_2 Gathering terrestrial plants", "5_3 Logging & wood harvesting", "5_4 Fishing & harvesting aquatic resources",
  "6_1 Recreational activities", "6_2 War, civil unrest & military exercises", "6_3 Work & other activities",
  "7_1 Fire & fire suppression", "7_2 Dams & water management/use", "7_3 Other ecosystem modifications",
  "8_1 Invasive non-native/alien species/diseases", "8_2 Problematic native species/diseases", "8_3 Introduced genetic material", "8_4 Problematic species/diseases of unknown origin", "8_5 Viral/prion-induced diseases", "8_6 Diseases of unknown cause",
  "9_1 Domestic & urban waste water", "9_2 Industrial & military effluents", "9_3 Agricultural & forestry effluents", "9_4 Garbage & solid waste", "9_5 Air-borne pollutants", "9_6 Excess energy",
  "10_1 Volcanoes", "10_2 Earthquakes/tsunamis", "10_3 Avalanches/landslides",
  "11_1 Habitat shifting & alteration", "11_2 Droughts", "11_3 Temperature extremes", "11_4 Storms & flooding", "11_5 Other impacts",
  "12_1 Other threat",
  "WOOD_PROD Custom: Wood production (2.2 and/or 5.3)" # Add our custom label to the lookup table
))

threat_metadata <- lookup_raw %>%
  extract(text, into = c("threats", "threat_name"), regex = "^([A-Z0-9_]+)\\s+(.*)$")

# ==========================================
# 3. Clean, Prepare & Inject Custom Threat
# ==========================================
# Standard cleaning and extraction step
clean_long_base <- df %>%
  select(scientific_name, taxa_group, threats) %>%
  filter(!is.na(threats) & threats != "") %>%
  separate_rows(threats, sep = ";\\s*") %>% 
  filter(!is.na(threats) & threats != "") %>%
  mutate(threats = str_extract(threats, "^\\d+[._]\\d+"),
         threats = str_replace(threats, "\\.", "_")) %>%
  filter(!is.na(threats)) %>%
  distinct(scientific_name, taxa_group, threats)

# Identify species that meet our "wood production" criteria (threatened by 5_3 and/or 2_2)
wood_prod_species <- clean_long_base %>%
  filter(threats %in% c("5_3", "2_2")) %>%
  distinct(scientific_name, taxa_group) %>%
  mutate(threats = "WOOD_PROD") # Tag them with our custom indicator

# Combine original dataset with our custom injected rows
clean_long_data <- bind_rows(clean_long_base, wood_prod_species)


# ==========================================
# 4. Generate the ALL VERTEBRATES Summary
# ==========================================
total_vert_species <- n_distinct(df$scientific_name)

all_vertebrates_summary <- clean_long_data %>%
  group_by(threats) %>%
  summarize(
    species_count = n(),
    proportion = species_count / total_vert_species,
    percentage = proportion * 100
  ) %>%
  left_join(threat_metadata, by = "threats") %>%
  mutate(threat_code_clean = str_replace(threats, "_", "."),
         threat_full_label = ifelse(threats == "WOOD_PROD", threat_name, paste0(threat_code_clean, ": ", threat_name))) %>%
  select(threat_full_label, species_count, proportion, percentage) %>%
  arrange(desc(species_count))


# ==========================================
# 5. Generate SEPARATE Summaries for Each Taxa
# ==========================================
taxa_denominators <- df %>%
  group_by(taxa_group) %>%
  summarize(total_group_species = n_distinct(scientific_name))

taxa_summaries <- clean_long_data %>%
  group_by(taxa_group, threats) %>%
  summarize(species_count = n(), .groups = "drop") %>%
  left_join(taxa_denominators, by = "taxa_group") %>%
  mutate(
    proportion = species_count / total_group_species,
    percentage = proportion * 100
  ) %>%
  left_join(threat_metadata, by = "threats") %>%
  mutate(threat_code_clean = str_replace(threats, "_", "."),
         threat_full_label = ifelse(threats == "WOOD_PROD", threat_name, paste0(threat_code_clean, ": ", threat_name))) %>%
  select(taxa_group, threat_full_label, species_count, total_group_species, proportion, percentage) %>%
  arrange(taxa_group, desc(species_count))


# ==========================================
# 6. Print Results (Safely bypassing AppLocker)
# ==========================================

# --- Print All Vertebrates Combined ---
cat("\n======================================================\n")
cat("SUMMARY FOR ALL VERTEBRATES COMBINED (WITH CUSTOM WOOD PRODUCTION)\n")
cat(paste("Total Unique Species Evaluated:", total_vert_species, "\n"))
cat("======================================================\n")
print(as.data.frame(head(all_vertebrates_summary, 20))) # Extended head to 20 to ensure custom category is visible

# --- Print Taxa-by-Taxa Breakdowns ---
unique_taxa <- unique(taxa_summaries$taxa_group)

for (taxon in unique_taxa) {
  # Filter out data for just this specific group
  taxon_data <- taxa_summaries %>% filter(taxa_group == taxon)
  sample_size <- taxon_data$total_group_species[1]
  
  cat("\n------------------------------------------------------\n")
  cat(paste("TAXA GROUP SUMMARY:", toupper(taxon), "\n"))
  cat(paste("Total Unique Species in this Group:", sample_size, "\n"))
  cat("------------------------------------------------------\n")
  
  display_table <- taxon_data %>% select(-taxa_group, -total_group_species)
  print(as.data.frame(head(display_table, 20)))
}