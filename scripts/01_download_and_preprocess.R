source(file.path(DIR_SCRIPTS,"00_setup.R"))

# Load raw counts
counts <- fread(file.path(DIR_RAW,"GSE198862_137_HYC_counts.csv.gz"), data.table = FALSE)

# Set row names as the first column
rownames(counts) <- counts[[1]]

# removing the duplicate gene column
counts <- counts[, -1]
head(counts)

# Download Metadata
gse <- readRDS(file.path(DIR_RAW,"GSE198862_ExpressionSet.rds"))
metadata <- pData(gse[[1]])

# subsetting sample id from metadata to have same sample id in both metadata and raw count matrix
metadata$sample_id <- sub("\\s*\\(.*", "", metadata$title)
head(metadata)

# Create DGE list object to compute library size 
dge <- DGEList(counts = counts)

# Normalize to account for library size 
dge <- calcNormFactors(dge)

# Calculate CPM on LINEAR scale to filter low count reads
cpm_mat <- cpm(dge)

# Filter genes with low expressed counts and it should be atleast 3 samples 
keep <- rowSums(cpm_mat >= 1) >= 3

# Subset RAW counts (do NOT transform here)
counts_filtered <- counts[keep, ]

# Re-create DGE List using FILTERED raw counts
dge_filtered <- DGEList(counts = counts_filtered)

# Normalization done for filtered gene-set 
dge_filtered <- calcNormFactors(dge_filtered)

#logCPM for visualization
logCPM_counts <- cpm(dge_filtered, log = TRUE, prior.count = 1)

# Save processed objects
saveRDS(counts_filtered,file.path(DIR_PROCESSED,"counts_filtered.rds"))
saveRDS(logCPM_counts,file.path(DIR_PROCESSED,"logCPM_counts.rds"))
saveRDS(metadata, file.path(DIR_PROCESSED, "metadata_clean.rds"))
