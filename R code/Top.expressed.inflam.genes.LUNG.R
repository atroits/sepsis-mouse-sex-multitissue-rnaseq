# -------------------------
# 1. Load libraries
# -------------------------
library(readr)
library(dplyr)
library(tibble)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -------------------------
# 2. Load counts file
# -------------------------
counts <- read_tsv("~/Desktop/RNASEQ/DESeq2 data YMF/lung.full.ymf/counts.txt")
counts <- counts %>%
  column_to_rownames(var = colnames(counts)[1])  # first column = Ensembl IDs

# -------------------------
# 3. Map Ensembl IDs → gene symbols
# -------------------------
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(counts),
  mart = mart
)

counts_annot <- counts %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_map, by = "ensembl_gene_id") %>%
  mutate(gene_symbol = ifelse(mgi_symbol == "" | is.na(mgi_symbol), ensembl_gene_id, mgi_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)

counts <- counts_annot %>%
  column_to_rownames("gene_symbol")
counts <- counts[, !(colnames(counts) %in% c("ensembl_gene_id", "mgi_symbol"))]

# -------------------------
# 4. Load DE results
# -------------------------
de_female <- read_csv("~/Desktop/RNASEQ/DESeq2 data YMF/biol sex/female.slur.dext.lung.csv") %>%
  rename(
    FoldChange = `Fold change (slurry,female vs dextrose,female)`,
    FDR = `FDR step up (slurry,female vs dextrose,female)`,
    gene_id = `Gene ID`
  ) %>%
  mutate(
    log2FC = log2(abs(FoldChange)) * sign(FoldChange),
    FDR = as.numeric(FDR)
  ) %>%
  left_join(gene_map, by = c("gene_id" = "ensembl_gene_id")) %>%
  mutate(gene_symbol = ifelse(mgi_symbol == "" | is.na(mgi_symbol), gene_id, mgi_symbol))

de_male <- read_csv("~/Desktop/RNASEQ/DESeq2 data YMF/biol sex/gsa_resultmale.csv") %>%
  rename(
    FoldChange = `Fold change (slurry,male vs dextrose,male)`,
    FDR = `FDR step up (slurry,male vs dextrose,male)`,
    gene_id = `Gene ID`
  ) %>%
  mutate(
    log2FC = log2(abs(FoldChange)) * sign(FoldChange),
    FDR = as.numeric(FDR)
  ) %>%
  left_join(gene_map, by = c("gene_id" = "ensembl_gene_id")) %>%
  mutate(gene_symbol = ifelse(mgi_symbol == "" | is.na(mgi_symbol), gene_id, mgi_symbol))

# -------------------------
# 5. Top 30 upregulated inflammatory genes per sex
# -------------------------
FDR_thresh <- 0.05
log2FC_thresh <- log2(1.5)

inflammatory_genes <- read.delim("~/Desktop/RNASEQ/inflammatory_120_list.txt", header = FALSE)$V1

df_inflam_female <- de_female %>%
  filter(FDR < FDR_thresh, log2FC > log2FC_thresh, gene_symbol %in% inflammatory_genes) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 30)

df_inflam_male <- de_male %>%
  filter(FDR < FDR_thresh, log2FC > log2FC_thresh, gene_symbol %in% inflammatory_genes) %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 30)

# -------------------------
# 6. Combine & deduplicate top inflammatory genes
# -------------------------
combined_genes <- unique(c(df_inflam_female$gene_symbol, df_inflam_male$gene_symbol))

# -------------------------
# 7. Prepare expression matrix
# -------------------------
expr_mat <- counts[rownames(counts) %in% combined_genes, ]
expr_log <- log2(expr_mat + 1)
expr_scaled <- t(scale(t(expr_log)))

# -------------------------
# 8. Sample annotations and order
# -------------------------
sample_groups <- data.frame(
  Condition = factor(c(rep("Slurry_Male",4),
                       rep("Slurry_Female",5),
                       rep("Dextrose_Male",3),
                       rep("Dextrose_Female",3)))
)
rownames(sample_groups) <- c("S1","S5","S40","S41",
                             "S6","S7","S34","S36","S38",
                             "S2","S3","S39",
                             "S33","S35","S37")

expr_scaled <- expr_scaled[, rownames(sample_groups)]

# -------------------------
# 9. ComplexHeatmap plotting with fixed -2,0,2 scale
# -------------------------
ha_col <- HeatmapAnnotation(df = sample_groups, col = list(
  Condition = c("Slurry_Male"="#66C2A5","Slurry_Female"="#FC8D62",
                "Dextrose_Male"="#8DA0CB","Dextrose_Female"="#E78AC3"))
)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue","white","red"))

Heatmap(
  expr_scaled,
  name = "Expression",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontface = "italic", fontsize = 10, col = "black"),
  column_names_gp = gpar(fontsize = 10),
  top_annotation = ha_col,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_title = "Top Upregulated Inflammatory Genes (Male + Female, Deduplicated)"
)
