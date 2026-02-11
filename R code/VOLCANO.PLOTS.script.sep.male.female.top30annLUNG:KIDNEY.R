# -------------------------
# 1. Load libraries
# -------------------------
library(ggplot2)
library(dplyr)
library(readr)
library(biomaRt)
library(ggrepel)

# -------------------------
# 2. Load DESeq2 results
# -------------------------
# Female
df_female <- read_csv("~/Desktop/RNASEQ/DESeq2 data YMF/biol sex/female.slur.dext.lung.csv") %>%
  rename(
    FoldChange = `Fold change (slurry,female vs dextrose,female)`,
    FDR = `FDR step up (slurry,female vs dextrose,female)`,
    gene_id = `Gene ID`
  ) %>%
  mutate(
    FoldChange = as.numeric(FoldChange),
    FDR = as.numeric(FDR),
    log2FC = log2(abs(FoldChange)) * sign(FoldChange)
  )

# Male
df_male <- read_csv("~/Desktop/RNASEQ/DESeq2 data YMF/biol sex/gsa_resultmale.csv") %>%
  rename(
    FoldChange = `Fold change (slurry,male vs dextrose,male)`,
    FDR = `FDR step up (slurry,male vs dextrose,male)`,
    gene_id = `Gene ID`
  ) %>%
  mutate(
    FoldChange = as.numeric(FoldChange),
    FDR = as.numeric(FDR),
    log2FC = log2(abs(FoldChange)) * sign(FoldChange)
  )

# -------------------------
# 3. Map Ensembl IDs → gene symbols
# -------------------------
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = unique(c(df_female$gene_id, df_male$gene_id)),
  mart = mart
)

df_female <- df_female %>%
  left_join(gene_map, by = c("gene_id" = "ensembl_gene_id")) %>%
  mutate(mgi_symbol = ifelse(mgi_symbol == "" | is.na(mgi_symbol), gene_id, mgi_symbol))

df_male <- df_male %>%
  left_join(gene_map, by = c("gene_id" = "ensembl_gene_id")) %>%
  mutate(mgi_symbol = ifelse(mgi_symbol == "" | is.na(mgi_symbol), gene_id, mgi_symbol))

# -------------------------
# 4. Define thresholds
# -------------------------
log2FC_thresh <- log2(1.5)
FDR_thresh <- 0.05

# Categorize genes
df_female <- df_female %>%
  mutate(Significance = case_when(
    FDR < FDR_thresh & log2FC > log2FC_thresh ~ "Upregulated",
    FDR < FDR_thresh & log2FC < -log2FC_thresh ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

df_male <- df_male %>%
  mutate(Significance = case_when(
    FDR < FDR_thresh & log2FC > log2FC_thresh ~ "Upregulated",
    FDR < FDR_thresh & log2FC < -log2FC_thresh ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

# -------------------------
# 5. Top 30 upregulated inflammatory genes per sex
# -------------------------
# Load inflammatory panel (assuming you saved as .txt with one gene per line)
inflammatory_genes <- read.delim("~/Desktop/RNASEQ/inflammatory_120_list.txt", header = FALSE)$V1

df_inflam_female <- df_female %>%
  filter(FDR < FDR_thresh, log2FC > log2FC_thresh, mgi_symbol %in% inflammatory_genes) %>%
  arrange(desc(log2FC)) %>%      # order by strongest upreg
  slice_head(n = 30)

df_inflam_male <- df_male %>%
  filter(FDR < FDR_thresh, log2FC > log2FC_thresh, mgi_symbol %in% inflammatory_genes) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 30)


# -------------------------
# 7. Color mapping
# -------------------------
color_map <- c(
  "Upregulated" = "red",
  "Downregulated" = "blue",
  "Not significant" = "grey"
)

# -------------------------
# 8. Shared axis limits
# -------------------------
x_limits <- range(c(df_female$log2FC, df_male$log2FC), na.rm = TRUE)
y_limits <- range(c(-log10(df_female$FDR), -log10(df_male$FDR)), na.rm = TRUE)

# -------------------------
# 9. Female volcano (circles)
# -------------------------
ggplot(df_female, aes(x = log2FC, y = -log10(FDR), color = Significance)) +
  geom_point(shape = 17, alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-log2FC_thresh, log2FC_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_thresh), linetype = "dashed") +
  scale_color_manual(values = color_map) +
  geom_text_repel(
    data = df_inflam_female,
    aes(label = mgi_symbol),
    size = 4, fontface = "italic", color = "black",
    box.padding = 0.5, point.padding = 0.3,
    max.overlaps = Inf, bg.color = "white", bg.r = 0.15
  ) +
  labs(title = "Volcano Plot: Lung Female (Slurry vs Dextrose)",
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("FDR")),
       color = "Gene Group") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  xlim(x_limits) +
  ylim(y_limits)

# -------------------------
# 10. Male volcano (squares)
# -------------------------
ggplot(df_male, aes(x = log2FC, y = -log10(FDR), color = Significance)) +
  geom_point(shape = 15, alpha = 0.7, size = 2) +
  geom_vline(xintercept = c(-log2FC_thresh, log2FC_thresh), linetype = "dashed") +
  geom_hline(yintercept = -log10(FDR_thresh), linetype = "dashed") +
  scale_color_manual(values = color_map) +
  geom_text_repel(
    data = df_inflam_male,
    aes(label = mgi_symbol),
    size = 4, fontface = "italic", color = "black",
    box.padding = 0.5, point.padding = 0.3,
    max.overlaps = Inf, bg.color = "white", bg.r = 0.15
  ) +
  labs(title = "Volcano Plot: Lung Male (Slurry vs Dextrose)",
       x = expression(log[2]("Fold Change")),
       y = expression(-log[10]("FDR")),
       color = "Gene Group") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  ) +
  xlim(x_limits) +
  ylim(y_limits)

