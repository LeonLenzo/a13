#!/usr/bin/env Rscript
# Split Defense category into Stress Response and Secondary Metabolites
# Generated: 2026-05-07

library(readr)
library(dplyr)
library(stringr)

# Read the completed manual curation
curation <- read_csv("go_hierarchy_analysis/manual_curation_COMPLETED.csv", show_col_types = FALSE)

# Split Defense into two subcategories
curation_split <- curation %>%
  mutate(
    Final_Category_Split = case_when(
      # Non-defense categories remain unchanged
      Final_Category != "Defense" ~ Final_Category,

      # Secondary Metabolites - specialized compounds for defense
      grepl("lignin|phenylpropanoid|flavonoid|terpenoid|diterpenoid|alkaloid",
            GO_Description, ignore.case = TRUE) ~ "Secondary Metabolites",

      # Cell wall modification (structural defense) - Secondary Metabolites
      grepl("cell wall", GO_Description, ignore.case = TRUE) ~ "Secondary Metabolites",

      # Chitin-related (anti-fungal defense) - could go either way, putting in Secondary Metabolites
      grepl("chitin", GO_Description, ignore.case = TRUE) ~ "Secondary Metabolites",

      # Ethylene biosynthesis (hormone for defense signaling) - Stress Response
      grepl("ethylene", GO_Description, ignore.case = TRUE) ~ "Stress Response",

      # Everything else in Defense goes to Stress Response
      # This includes: stress response, antioxidants, ROS, defense response to pathogen,
      # innate immune response, glutathione, detoxification, wound response, etc.
      TRUE ~ "Stress Response"
    )
  )

# Print summary
cat("\n=== DEFENSE CATEGORY SPLIT ===\n\n")
defense_split_summary <- curation_split %>%
  filter(Final_Category == "Defense") %>%
  group_by(Final_Category_Split) %>%
  summarise(N_Terms = n(), .groups = 'drop')

print(defense_split_summary)

cat("\n=== SAMPLE TERMS BY SUBCATEGORY ===\n\n")

cat("Stress Response terms (first 10):\n")
stress_terms <- curation_split %>%
  filter(Final_Category_Split == "Stress Response") %>%
  select(GO_Description) %>%
  head(10)
print(stress_terms, n = Inf)

cat("\nSecondary Metabolites terms (first 10):\n")
secondary_terms <- curation_split %>%
  filter(Final_Category_Split == "Secondary Metabolites") %>%
  select(GO_Description) %>%
  head(10)
print(secondary_terms, n = Inf)

# Save the split curation
write_csv(curation_split, "go_hierarchy_analysis/manual_curation_DEFENSE_SPLIT.csv")

cat("\n✅ Split curation saved to: go_hierarchy_analysis/manual_curation_DEFENSE_SPLIT.csv\n")
