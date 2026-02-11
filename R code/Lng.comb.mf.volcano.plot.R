# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# Load data
df <- read_csv("~/Desktop/RNASEQ/DESeq2 data YMF/lung.full.ymf/gsa_resultc.csv")

# Rename and prepare relevant columns
df <- df %>%
  rename(
    FoldChange = `Fold change (slurry vs dextrose)`,
    FDR = `FDR step up (slurry vs dextrose)`
  ) %>%
  mutate(
    FoldChange = as.numeric(FoldChange),
    FDR = as.numeric(FDR),
    log2FC = log2(abs(FoldChange)) * sign(FoldChange)
  )

# Define thresholds
log2FC_thresh <- log2(1.5)  # ~0.585
FDR_thresh <- 0.05

# Categorize genes
df <- df %>%
  mutate(Significance = case_when(
    FDR < FDR_thresh & log2FC > log2FC_thresh ~ "Upregulated",
    FDR < FDR_thresh & log2FC < -log2FC_thresh ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

# Count per group and label
group_counts <- df %>%
  count(Significance) %>%
  mutate(label = paste0(Significance, " (", n, ")"))

# Create label mapping
label_map <- setNames(group_counts$label, group_counts$Significance)

# Add labeled group to main df
df <- df %>%
  mutate(Significance_labeled = label_map[Significance])

# Define color map
color_map <- c(
  setNames("red", label_map["Upregulated"]),
  setNames("blue", label_map["Downregulated"]),
  setNames("grey", label_map["Not significant"])
)

# Volcano plot with fixed axes
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = Significance_labeled)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-log2FC_thresh, log2FC_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_thresh), linetype = "dashed") +
  scale_color_manual(values = color_map) +
  labs(
    title = "Lung: Slurry vs Dextrose",
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("FDR")),
    color = "Gene Group"
  ) +
  scale_x_continuous(limits = c(-5, 10)) +        # fixed x-axis from -3 to 3
  scale_y_continuous(limits = c(0, 50)) +        # fixed y-axis from 0 to 10
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )