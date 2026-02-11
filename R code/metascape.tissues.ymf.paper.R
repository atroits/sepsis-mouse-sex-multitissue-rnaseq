library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

plot_metascape <- function(file_path, tissue_name) {
  
  # Define groups with pretty labels
  groups <- list(
    "Group 1: Shared in all three tissues" = c(
      "Positive regulation of response to external stimulus",
      "Positive regulation of cytokine production",
      "Inflammatory response",
      "Leukocyte activation",
      "Cytokine-mediated signaling pathway",
      "Regulation of cell activation",
      "Positive regulation of cell migration",
      "Positive regulation of programmed cell death",
      "Regulation of MAPK cascade",
      "Regulation of protein phosphorylation",
      "Enzyme-linked receptor protein signaling pathway",
      "Pathways in cancer - Mus musculus (house mouse)",
      "Tube morphogenesis"
    ),
    "Group 2: Shared in two tissues" = c(
      "Response to molecule of bacterial origin",
      "Regulation of immune effector process",
      "Negative regulation of immune system process",
      "Negative regulation of cell population proliferation",
      "Hemostasis"
    ),
    "Group 3: Tissue-specific" = c(
      "Negative regulation of cell adhesion",
      "Cellular response to lipid",
      "Response to peptide hormone",
      "Signaling by Rho GTPases, Miro GTPases and RHOBTB3",
      "Response to growth factor",
      "Response to virus",
      "Response to interferon-beta",
      "Cytokine signaling in immune system",
      "Myeloid leukocyte activation",
      "Locomotion",
      "Negative regulation of cell differentiation"
      
    )
  )
  
  pathway_order <- unlist(groups, use.names = FALSE)
  group_map <- rep(names(groups), lengths(groups))
  
  # Load and clean input
  go_data <- read.csv(file_path) %>%
    mutate(Description = trimws(Description),
           GO_Term_lower = tolower(Description))
  
  pathway_df <- data.frame(
    GO_Term_master = pathway_order,
    GO_Term_lower = tolower(pathway_order),
    Group = group_map,
    stringsAsFactors = FALSE
  )
  
  merged <- pathway_df %>%
    left_join(go_data, by = "GO_Term_lower") %>%
    dplyr::select(Group, GO_Term_master,
           D_LogP = X_LogP_Downregulated,
           U_LogP = X_LogP_Upregulated)
  
  data_long <- merged %>%
    pivot_longer(cols = c(D_LogP, U_LogP),
                 names_to = "Direction",
                 values_to = "LogP") %>%
    mutate(
      Direction = ifelse(Direction == "U_LogP", "Upregulated", "Downregulated"),
      NES_like = ifelse(is.na(LogP), 0,
                        ifelse(Direction == "Upregulated", abs(LogP), -abs(LogP))),
      GO_Term = factor(GO_Term_master, levels = rev(pathway_order)),
      Group = factor(Group, levels = names(groups))
    )
  
  # Compute divider lines
  divider_positions <- data_long %>%
    distinct(GO_Term, Group) %>%
    mutate(y = as.numeric(GO_Term)) %>%
    group_by(Group) %>%
    summarise(min_y = min(y), .groups = "drop") %>%
    mutate(divider = min_y - 0.5) %>%
    pull(divider)
  
  # Plot
  ggplot(data_long, aes(x = NES_like, y = GO_Term, fill = Direction)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
    scale_x_continuous(
      limits = c(-20, 60),
      breaks = c(-20, 0, 20, 40, 60),
      labels = c("20", "0", "20", "40", "60")
    ) +
    scale_y_discrete(position = "right") +
    geom_hline(yintercept = divider_positions, color = "grey40", linewidth = 0.6) +
    labs(
      x = expression(-log[10](P-value)),
      y = NULL,
      title = paste0("GO Enrichment: ", tissue_name, "\nSlurry vs Dextrose")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey70", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.ticks.y = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text.y = element_text(size = 12, hjust = 0)
    )
}


plot_metascape("~/Desktop/YMF paper/paper latest versions/KIDNEY.heatmap.1.5.csv", "Kidney")
plot_metascape("~/Desktop/YMF paper/paper latest versions/LUNG.heatmap.1.5.csv", "Lung")
plot_metascape("~/Desktop/YMF paper/paper latest versions/LIVER.heatmap.1.5.csv", "Liver")
