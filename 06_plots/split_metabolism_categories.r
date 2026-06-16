#!/usr/bin/env Rscript
# Split Metabolism category into Photosynthesis and Biosynthesis
# Generated: 2026-05-07

library(readr)
library(dplyr)
library(stringr)

# Read the defense-split curation
curation <- read_csv("go_hierarchy_analysis/manual_curation_DEFENSE_SPLIT.csv", show_col_types = FALSE)

# Split Metabolism into two subcategories
curation_final <- curation %>%
  mutate(
    Final_Category_Full_Split = case_when(
      # Keep Defense subcategories as-is
      Final_Category_Split %in% c("Stress Response", "Secondary Metabolites", "Growth") ~ Final_Category_Split,

      # Split Metabolism into Photosynthesis vs Biosynthesis
      Final_Category_Split == "Metabolism" &
        grepl("photo|chloro|thylakoid|plastid", GO_Description, ignore.case = TRUE) ~ "Photosynthesis",

      # Everything else in Metabolism goes to Biosynthesis
      Final_Category_Split == "Metabolism" ~ "Biosynthesis",

      TRUE ~ Final_Category_Split
    )
  )

# Print summary
cat("\n=== METABOLISM CATEGORY SPLIT ===\n\n")
metabolism_split_summary <- curation_final %>%
  filter(Final_Category == "Metabolism") %>%
  group_by(Final_Category_Full_Split) %>%
  summarise(N_Terms = n(), .groups = 'drop')

print(metabolism_split_summary)

cat("\n=== ALL CATEGORY SUMMARY ===\n\n")
all_summary <- curation_final %>%
  group_by(Final_Category_Full_Split) %>%
  summarise(N_Terms = n(), .groups = 'drop') %>%
  arrange(desc(N_Terms))

print(all_summary)

cat("\n=== SAMPLE TERMS BY SUBCATEGORY ===\n\n")

cat("Photosynthesis terms (all):\n")
photo_terms <- curation_final %>%
  filter(Final_Category_Full_Split == "Photosynthesis") %>%
  select(GO_Description, Ontology)
print(photo_terms, n = Inf)

cat("\n\nBiosynthesis terms (first 20):\n")
biosyn_terms <- curation_final %>%
  filter(Final_Category_Full_Split == "Biosynthesis") %>%
  select(GO_Description) %>%
  head(20)
print(biosyn_terms, n = Inf)

# Save the full split curation
write_csv(curation_final, "go_hierarchy_analysis/manual_curation_FULL_SPLIT.csv")

cat("\n✅ Full split curation saved to: go_hierarchy_analysis/manual_curation_FULL_SPLIT.csv\n")
cat("\nFinal categories (5 total):\n")
cat("  1. Stress Response (Defense subcategory)\n")
cat("  2. Secondary Metabolites (Defense subcategory)\n")
cat("  3. Photosynthesis (Metabolism subcategory)\n")
cat("  4. Biosynthesis (Metabolism subcategory)\n")
cat("  5. Growth\n")
