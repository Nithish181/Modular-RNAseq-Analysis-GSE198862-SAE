source(file.path(DIR_SCRIPTS,"00_setup.R"))

subsets <- readRDS(file.path(DIR_PROCESSED,"subsets.rds"))

plot_corr <- function(counts, meta, title, out) {
  cor_mat <- cor(counts, method = "pearson")
  
  meta_ord <- meta[match(colnames(counts), meta$sample_id), ]
  ann <- data.frame(Condition = meta_ord$characteristics_ch1.1)
  rownames(ann) <- colnames(cor_mat)
  
  png(out, width = 900, height = 800)
  pheatmap(
    cor_mat,
    annotation_row = ann,
    annotation_col = ann,
    clustering_method = "complete",
    main = title
  )
  dev.off()
}

plot_corr(
  subsets$d3$counts,
  subsets$d3$meta,
  "Sample Correlation: 3d Sham vs PCI",
  "figures/correlation/day3_sham_vs_pci.png"
)

plot_corr(
  subsets$d20$counts,
  subsets$d20$meta,
  "Sample Correlation: 20d Sham vs PCI",
  "figures/correlation/day20_sham_vs_pci.png"
)

plot_corr(
  subsets$pci$counts,
  subsets$pci$meta,
  "Sample Correlation: PCI 3d vs 20d",
  "figures/correlation/pci_day3_vs_day20.png"
)
