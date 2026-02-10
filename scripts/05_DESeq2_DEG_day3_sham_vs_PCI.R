# --------------------------------------------------- 
#                 library requirements
# ---------------------------------------------------

library(GEOquery)
library(DESeq2)
library(readr)
library(data.table)
library(ggrepel)
library(ggplot2)
#----------------------------------------------------
#                       INPUT 
#-------------------------------------------------------

# Loading the counts
count <-fread("GSE198862_137_HYC_counts.csv.gz", data.table = F)

#set row names as the first column 
row.names(count) <- count[[1]]
count<- count[,-1]

# Downloading the metadata from GEO
gse <- getGEO("GSE198862", GSEMatrix = TRUE)
metadata = pData(gse[[1]])

#-------------------------------------------------------------------
#                   Data wrangling
#-------------------------------------------------------------------
# creating a sample id column in metadata
metadata$sampleID <- sub("\\s*\\(.*", "",metadata$title)

# filtering low count genes to avoid nosiy genes which is not expressed in 80% samples
perc_keep <- 0.8
gene_keep <- rowSums(count > 0) >= ceiling(perc_keep * ncol(count))
count_tbl_low_rm <- count[gene_keep, ]

#-------------------------------------------------------------------
#              Data subset according to the study 
#------------------------------------------------------------------
# subset of metadata
D_3_meta <- metadata[metadata$characteristics_ch1.1 %in% c("treatment: sham","treatment: PCI")& metadata$characteristics_ch1.2 == "time: 3 days", ]
# subset for counts 
D_3_count <- count_tbl_low_rm[, D_3_meta$sampleID, drop = F]

# metadata for DEseq2
meta <- data.frame(
  sampleID = colnames(D_3_count),
  Treatment = D_3_meta$`treatment:ch1`
)
#------------------------------------------------------------------
#                     Experimental Design for this study
#-----------------------------------------------------------------
# Identifying the up & downregulated genes in 3 days sham vs PCI
# PCI == Peritoneal contamination and infection 
# sham == healthy group

#-----------------------------------------------------------------
#                              DEG testing using Deseq2
#-----------------------------------------------------------------
# DEseq2 DEG
# creating Deseq2 object 
obj_deseq2 <- DESeqDataSetFromMatrix(
  countData = D_3_count,
  colData = meta,
  design = ~ Treatment
)

# set reference level and perfrom DE analysis
obj_deseq2$Treatment <- relevel(obj_deseq2$Treatment, ref = "sham")
obj_deseq2 <- DESeq(obj_deseq2)

res <- results(obj_deseq2)
deg_results <- as.data.frame(res)

# To calculate overlap between 2 or more methods 
day3_deseq <- rownames(subset(deg_results, padj < 0.05 & abs(log2FoldChange) > 1))
#-----------------------------------------------------------------
#                                   NOTE
#-----------------------------------------------------------------
# The res from the deseq2 is used for further analysis like GSEA while the other
#  formatted outputs separated using cutoff are for visual representation 

# Reason in further analysis we will order genes based on the statistic test rather thresholds becuase thresholds
# leads to missing biology


#-----------------------------------------------------------------
#                                  VISUALIZATION OF RESULTS
#-----------------------------------------------------------------
# Refining the log2change values from Deseq2
resLFC <- lfcShrink(obj_deseq2, coef = "Treatment_PCI_vs_sham", type = "apeglm")

res_deseq2 <- resLFC[order(resLFC$padj),]

deg <- as.data.frame(res_deseq2)

# converting log2FC to log10 scale to make plots 
deg$logP <- -log10(deg$padj)

# Adding annotations according to the log2FC
deg$significance <- "Not Sig"

# separating up and down regualted genes based on  adj.p-value
deg$significance[deg$padj < 0.05 & deg$log2FoldChange >1 ] <- "up"

up_regualted_genes <- deg[deg$significance == "up",]

deg$significance[deg$padj<0.05 & deg$log2FoldChange < -1] <- "down"

down_regualted_genes <- deg[deg$significance == "down", ]

#-------------------------------------------------------------
#                             VOLCANO PLOT
#-------------------------------------------------------------

# plotting the DEG's 
ggplot(deg, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "up" = "red",
    "down" = "blue",
    "Not Sig" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(
    title = "Volcano Plot (Day-3 PCI vs sham)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  )



