source(file.path(DIR_SCRIPTS,"00_setup.R"))

subsets <- readRDS("subsets.rds")

run_pca <- function(counts, meta, out) {
  pca <- prcomp(t(counts), scale. = TRUE)
  
  df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Time = ifelse(grepl("3 days", meta$characteristics_ch1.2), "3days", "20days"),
    Condition = ifelse(grepl("PCI", meta$characteristics_ch1.1), "PCI", "sham"),
    Sample = rownames(pca$x)
  )
  
  p <- ggplot(df, aes(PC1, PC2, color = Time, shape = Condition)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = Sample), size = 3) +
    theme_minimal()
  
  ggsave(out, p, width = 7, height = 6)
}

run_pca(subsets$d3$counts, subsets$d3$meta,
        "figures/pca/pca_day3_sham_vs_pci.png")

run_pca(subsets$d20$counts, subsets$d20$meta,
        "figures/pca/pca_day20_sham_vs_pci.png")

run_pca(subsets$pci$counts, subsets$pci$meta,
        "figures/pca/pca_pci_day3_vs_day20.png")
