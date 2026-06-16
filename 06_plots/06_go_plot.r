#!/usr/bin/env Rscript
# GO enrichment three-column plot
# Input:  05_clusterprofiler/results/GO_enrichment_*.csv
#         03_kallisto/reference/manual_curation_FULL_SPLIT.csv
# Output: 06_plots/GO_enrichment_three_column.{png,pdf,svg}
# Run from: toxa13-wheat-rna/
# Requires: 05_go_enrichment.r to have been run first

library(readr)
library(dplyr)
library(ggplot2)
library(ggtext)
library(ggh4x)
library(stringr)

cat("=== GO Enrichment Plot ===\n\n")

comparisons <- c("SN15_vs_Tween", "3KO_vs_Tween", "3KO_vs_SN15")

all_enrichment_results <- list()
for (comp_name in comparisons) {
  f <- paste0("05_clusterprofiler/results/GO_enrichment_", comp_name, ".csv")
  if (file.exists(f)) {
    all_enrichment_results[[comp_name]] <- read_csv(f, show_col_types = FALSE)
  }
}

if (length(all_enrichment_results) == 0) {
  stop("No GO enrichment results — run 05_clusterprofiler/05_go_enrichment.r first.")
}

curation <- read_csv("reference/manual_curation_FULL_SPLIT.csv",
                     show_col_types = FALSE) %>%
  select(GO_ID, Subcategory = Final_Category_Full_Split) %>%
  filter(!is.na(Subcategory), Subcategory != "Growth") %>%
  mutate(
    Category = ifelse(Subcategory %in% c("Stress Response", "Secondary Metabolites"),
                      "Defense", "Metabolism")
  )

cat("Curation lookup loaded:", nrow(curation), "GO terms\n\n")

ontology_colors <- c("BP" = "#112180", "MF" = "#990c0c", "CC" = "#801c6f")

dataset_labels <- c(
  "SN15_vs_Tween" = "SN15 vs Tween",
  "3KO_vs_Tween"  = "3KO vs Tween",
  "3KO_vs_SN15"   = "3KO vs SN15"
)

subcategory_levels <- c("Stress Response", "Secondary Metabolites",
                        "Photosynthesis", "Biosynthesis")
category_levels    <- c("Defense", "Metabolism")

combined_data <- bind_rows(all_enrichment_results) %>%
  mutate(Dataset = dataset_labels[Dataset])

plot_data <- combined_data %>%
  left_join(curation, by = "GO_ID") %>%
  filter(!is.na(Subcategory)) %>%
  group_by(Subcategory, Dataset) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    GO_Description_Short = ifelse(
      nchar(GO_Description) > 50,
      paste0(substr(GO_Description, 1, 50), "..."),
      GO_Description
    ),
    ont_color = ontology_colors[Ontology],
    GO_label  = paste0(
      "<span style='color:", ont_color, "'>", GO_Description_Short, "</span>"
    ),
    Category    = factor(Category,    levels = category_levels),
    Subcategory = factor(Subcategory, levels = subcategory_levels),
    Dataset     = factor(Dataset,     levels = unname(dataset_labels))
  )

p <- ggplot(plot_data,
            aes(x = RichFactor, y = reorder(GO_label, -p.adjust))) +
  geom_point(aes(fill = -log10(p.adjust), size = Count,
                 shape = Direction, color = Ontology),
             stroke = 0.3, alpha = 0.8) +
  scale_fill_gradient(low = "blue", high = "red",
                      name = "-log10(Adj. P-value)") +
  scale_color_manual(name = "GO Terms",
                     values = ontology_colors,
                     guide = guide_legend(
                       order = 4, title.position = "top", title.hjust = 0.5,
                       ncol = 3, override.aes = list(size = 3, shape = 16, alpha = 1)
                     )) +
  scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
  scale_shape_manual(values = c("Up" = 24, "Down" = 25), name = "Direction") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.1, 0.1)),
                     breaks = c(0, 0.5, 1)) +
  scale_y_discrete(position = "right", expand = expansion(add = c(1.2, 1.2))) +
  facet_nested(Category + Subcategory ~ Dataset,
               scales = "free_y", space = "free_y", switch = "y",
               strip = strip_nested(
                 background_y = list(
                   element_rect(fill = "gray85", color = "black", linewidth = 1),
                   element_blank()
                 ),
                 text_y = list(
                   element_text(size = 12, face = "bold", angle = 90),
                   element_text(size = 10, angle = 90)
                 ),
                 by_layer_y = TRUE
               )) +
  theme_bw() +
  theme(
    axis.text.y.right    = element_markdown(size = 8, face = "italic"),
    axis.text.x          = element_text(size = 10),
    axis.text.y.left     = element_blank(),
    axis.ticks           = element_blank(),
    axis.ticks.x.bottom  = element_line(color = "black", linewidth = 0.3),
    legend.position      = "bottom",
    legend.box           = "horizontal",
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    panel.grid.major.y   = element_line(color = "lightgray", linewidth = 0.1),
    strip.text.x         = element_text(size = 12, face = "bold"),
    strip.placement      = "outside",
    panel.border         = element_rect(color = "#000000", fill = NA, linewidth = 0.5),
    ggh4x.facet.nestline = element_line(color = "black", linewidth = 1)
  ) +
  labs(x = "Rich Factor", y = "") +
  guides(
    fill  = guide_colorbar(order = 1, title.position = "top",
                           title.hjust = 0.5, barwidth = 8, barheight = 0.5,
                           direction = "horizontal"),
    size  = guide_legend(order = 2, title.position = "top", title.hjust = 0.5,
                         override.aes = list(shape = 21), ncol = 3),
    shape = guide_legend(order = 3, title.position = "top",
                         title.hjust = 0.5, ncol = 2)
  )

ggsave("06_plots/GO_enrichment_three_column.png", p, width = 12, height = 14, dpi = 300, bg = "white")
ggsave("06_plots/GO_enrichment_three_column.pdf", p, width = 12, height = 14, bg = "white")
ggsave("06_plots/GO_enrichment_three_column.svg", p, width = 12, height = 14, bg = "white")

cat("Saved: 06_plots/GO_enrichment_three_column.{png,pdf,svg}\n")
cat("=== GO Plot Complete ===\n")
