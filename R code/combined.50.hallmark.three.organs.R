# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)
library (readr)
library(viridisLite)  # For colorblind-friendly color palette

names <- c("Immune & Inflammation", "Proliferation & Cell Cycle", "Metabolism", "Canonical Signaling", "Stress & Damage Response", "Development & Differentiation")
names_ordered <- c("Immune & Inflammation", "Stress & Damage Response", "Proliferation & Cell Cycle", "Metabolism", "Canonical Signaling", "Development & Differentiation")


category1 = c(
  "INFLAMMATORY RESPONSE", "INTERFERON ALPHA RESPONSE", "INTERFERON GAMMA RESPONSE",
  "IL2 STAT5 SIGNALING", "IL6 JAK STAT3 SIGNALING", "COMPLEMENT",
  "ALLOGRAFT REJECTION", "TNFA SIGNALING VIA NFKB", "COAGULATION"
)
category2 = c(
  "E2F TARGETS", "G2M CHECKPOINT", "MITOTIC SPINDLE",
  "MYC TARGETS V1", "MYC TARGETS V2", "DNA REPAIR"
)
category3 = c(
  "FATTY ACID METABOLISM", "BILE ACID METABOLISM", "ADIPOGENESIS",
  "OXIDATIVE PHOSPHORYLATION", "XENOBIOTIC METABOLISM", "HEME METABOLISM",
  "CHOLESTEROL HOMEOSTASIS", "HYPOXIA", "GLYCOLYSIS", "PROTEIN SECRETION", "PEROXISOME"
)
category4 = c(
  "MTORC1 SIGNALING", "PI3K AKT MTOR SIGNALING", "TGF BETA SIGNALING",
  "NOTCH SIGNALING", "WNT BETA CATENIN SIGNALING", "HEDGEHOG SIGNALING",
  "KRAS SIGNALING UP", "KRAS SIGNALING DN", "ESTROGEN RESPONSE EARLY",
  "ESTROGEN RESPONSE LATE", "ANDROGEN RESPONSE"
)
category5 = c(
  "UNFOLDED PROTEIN RESPONSE", "REACTIVE OXYGEN SPECIES PATHWAY",
  "P53 PATHWAY", "UV RESPONSE UP", "UV RESPONSE DN", "APOPTOSIS"
)
category6 = c(
  "EPITHELIAL MESENCHYMAL TRANSITION", "ANGIOGENESIS",
  "MYOGENESIS", "SPERMATOGENESIS", "APICAL JUNCTION", "APICAL SURFACE", "PANCREAS BETA CELLS"
)


# Read in GSEA results
lung   <- read_csv("sasha/LUNG.COMBINED.TABLE.csv")
# Assuming your dataframe is named df and the column is called "text"
lung$Pathways <- gsub("  +", " ", lung$Pathways)

kidney <- read_csv("sasha/KIDNEY.COMBINED.TABLE.csv")
kidney$Pathways <- gsub("  +", " ", kidney$Pathways)

liver  <- read_csv("sasha/LIVER.COMBINED.TABLE.csv")
liver$Pathways <- gsub("  +", " ", liver$Pathways)

# Label groups
lung$Organ <- "Lung"
kidney$Organ <- "Kidney"
liver$Organ <- "Liver"

# Standardize column names
colnames(lung) <- c("Pathway", "GeneCount", "NES", "FDR", "Organ")
colnames(kidney) <- c("Pathway", "GeneCount", "NES", "FDR", "Organ")
colnames(liver) <- c("Pathway", "GeneCount", "NES", "FDR", "Organ")

# Combine datasets
combined <- bind_rows(lung, kidney, liver)

# Set Organ order: Liver, Lung, Kidney
combined$Organ <- factor(combined$Organ, levels = c("Liver", "Kidney", "Lung"))


combined$PathwayGroup <- case_when(
  combined$Pathway %in% category1 ~ names[1],
  combined$Pathway %in% category2 ~ names[2],
  combined$Pathway %in% category3 ~ names[3],
  combined$Pathway %in% category4 ~ names[4],
  combined$Pathway %in% category5 ~ names[5],
  combined$Pathway %in% category6 ~ names[6],
  
)


# Reorder PathwayGroup by descending average NES
group_order <- combined %>%
  group_by(PathwayGroup) %>%
  summarise(avg_NES = mean(NES, na.rm = TRUE)) %>%
  arrange(desc(avg_NES)) %>%
  pull(PathwayGroup)

combined$PathwayGroup <- factor(combined$PathwayGroup, levels = names_ordered)

# Reorder Pathways within each group by average NES
combined <- combined %>%
  group_by(PathwayGroup, Pathway) %>%
  mutate(avg_pathway_NES = mean(NES, na.rm = TRUE)) %>%
  ungroup()

# Create combined label for sorting purposes
combined <- combined %>%
  mutate(Pathway_combined = paste(PathwayGroup, Pathway, sep = "___"))

# Set levels of the combined factor based on average NES
pathway_order <- combined %>%
  distinct(PathwayGroup, Pathway_combined, avg_pathway_NES) %>%
  arrange(PathwayGroup, avg_pathway_NES) %>%
  pull(Pathway_combined)

save(pathway_order, file = "pathway_order_hallmark.RData")
#load("pathway_order_hallmark.RData")


combined$Pathway_combined <- factor(combined$Pathway_combined, levels = pathway_order)

# Update the Pathway factor (used on y-axis) 
combined <- combined %>%
  arrange(PathwayGroup, avg_pathway_NES) %>%  
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)))

# Create faceted dot plot with jitter to avoid overlapping
p <- ggplot(combined, aes(x = NES, y = Pathway)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(
    aes(size = GeneCount, fill = FDR),
    shape = 21,                        # Circle with fill and border
    color = "black",                   # Black border
    stroke = 0.3,                      # Thin border
    alpha = 0.9,
    position = position_jitter(width = 0.2, height = 0)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = -1,
    name = "FDR q-val",
    trans = "reverse"
  ) +
  scale_size_continuous(range = c(3, 10), name = "Gene Count") +
  facet_grid(PathwayGroup ~ Organ, scales = "free_y", space = "free_y")+
  #  facet_grid(. ~ Organ) +
  labs(
    title = "GSEA Hallmark Pathways: Septic Response Comparison across Tissues",
    x = "Normalized Enrichment Score (NES)",
    y = ""
  ) +
  guides(
    fill = guide_colorbar(order = 1, title.position = "top", title.hjust = 0.5),
    size = guide_legend(order = 2, title.position = "top", title.hjust = 0.5)
  ) +
  theme_bw(base_size = 17) +
  theme(
    axis.text.y = element_text(size = 18),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.box = "vertical",
    legend.position = "right",
    legend.title.align = 0.5,
    legend.spacing.y = unit(0.8, "lines")
  )

# Show the plot
print(p)

# Save the plot as high-resolution PNG
ggsave("GSEA_dotplot_3organs_all_pathways.png", plot = p, width = 20, height = 24, dpi = 300)
