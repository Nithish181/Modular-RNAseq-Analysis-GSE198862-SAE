# --------------------------------------------------- 
#              library requirements
# ---------------------------------------------------

library(GEOquery)
library(edgeR)
library(readr)
library(data.table)
library(ggrepel)
library(ggplot2)
#----------------------------------------------------
#                       INPUT 
#-------------------------------------------------------
# Loading the counts
count <-fread("GSE198862_137_HYC_counts.csv.gz", data.table = F)

#set rownames as the first column 
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

# filtering low count genes 
#perc_keep <- 0.8
#gene_keep <- rowSums(count > 0) >= ceiling(perc_keep * ncol(count))
#count_tbl_low_rm <- count[gene_keep, ]

#-------------------------------------------------------------------
#              Data subset according to the study 
#------------------------------------------------------------------

# subset of metadata
D_20_meta <- metadata[metadata$characteristics_ch1.1 %in% c("treatment: sham","treatment: PCI")& metadata$characteristics_ch1.2 == "time: 20 days", ]

D_20_count <- count[, D_20_meta$sampleID, drop = F]


#------------------------------------------------------------------
#                     Experimental Design for this study
#-----------------------------------------------------------------
# Identifying the up & downregulated genes in 3 days sham vs PCI
# PCI == Peritoneal contamination and infection 
# sham == healthy group

#-----------------------------------------------------------------
#                              DEG testing using edgeR
#-----------------------------------------------------------------

# metadata for edgeR
meta <- data.frame(
  sampleID = colnames(D_20_count),
  Treatment = D_20_meta$`treatment:ch1`
)

# create a DEGList object
dge <- DGEList(counts = D_20_count, samples = meta)

# Set factor + reference
dge$samples$Treatment <- factor(dge$samples$Treatment, levels = c("sham", "PCI"))


# Design matrix
design <- model.matrix(~ Treatment, data = dge$samples)

# 5. Filter using edgeR's method
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes=FALSE]


# Normalize
dge <- calcNormFactors(dge, method = "TMM")

# Estimate dispersion
dge <- estimateDisp(dge, design, robust = T)

# 8. Check dispersion plot
plotBCV(dge)

# Fit model
fit <- glmFit(dge, design, robust = T)

# Test PCI vs sham
qlf <- glmLRT(fit, coef = 2)

res_qlf <- topTags(qlf, n = Inf, adjust.method = "BH")
res_qlf_df <- res_qlf$table

# to caluclate the overlap between 2 or more methods 
day20_edger <- rownames(subset(res_qlf_df, FDR < 0.05 & abs(logFC) > 1))

# Load DEG table
dega <- res_qlf_df

# Convert FDR to -log10 scale
dega$logP <- -log10(dega$FDR)

# Create significance labels
dega$significance <- "Not Sig"

dega$significance[
  dega$FDR < 0.05 & dega$logFC > 1
] <- "Up"

dega$significance[
  dega$FDR < 0.05 & dega$logFC < -1
] <- "Down"

# Volcano plot
ggplot(dega, aes(x = logFC, y = logP, color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "Not Sig" = "grey"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(
    title = "EdgeR Volcano Plot 20 days (PCI vs Sham)",
    x = "log2 Fold Change",
    y = "-log10(FDR)"

  )
