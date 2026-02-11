# Sex-Specific Murine Sepsis RNA-seq Analysis

This repository contains analysis inputs and scripts used in the study investigating sex-specific transcriptional responses to experimental murine sepsis across lung, kidney, and liver tissues.

## Differential Expression
Differential expression analysis was performed using DESeq2 within Partek Flow.

## GSEA Analysis
Gene Set Enrichment Analysis (GSEA Desktop) was conducted using:
- Expression matrices (.gct)
- Phenotype label files (.cls)
- Ranked metric files (.rnk)

GSEA was condicted on combiend septic responses (male + female) for lung, kidney, and liver, and on separate septic responses (male vs female) for lung and kidney.  

Hallmark gene sets were obtained from MSigDB (version mh.all.v2025.1.Mm).
GSEA was performed on expression datasets (normalized counts files generated with DESeq2 from Partek). 
All other GSEA analysis settings were used as default.

## Metascape Analysis
Gene lists used for Metascape enrichment analysis are included in the `Metascape input` directory.

## Data Availability
Raw RNA-sequencing data will be available in GEO under accession number GSEXXXXX (to be released upon publication).
