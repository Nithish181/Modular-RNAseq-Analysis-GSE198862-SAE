# Sepsis-Associated-Encephalopathy-SAE-temporal-transcriptomics
Benchmarking DESeq2 vs edgeR for temporal analysis of sepsis-associated encephalopathy (GSE198862)
# Exploratory Transcriptomic Analysis of Sepsis-Associated Encephalopathy (GSE198862)

This repository contains a modular exploratory data analysis (EDA) pipeline for bulk RNA-seq data generated from a murine model of sepsis-associated encephalopathy (SAE). The analysis is based on the publicly available GEO dataset **GSE198862**, which profiles brain transcriptomic changes following peritoneal contamination and infection (PCI) at multiple time points.

The primary goal of this project is to assess data quality, sample coherence, and global gene-expression structure prior to performing differential expression and downstream functional analyses.

---

## Dataset

- **GEO accession:** GSE198862  
- **Organism:** *Mus musculus*  
- **Experimental design:**  
  - Sham vs PCI treatment  
  - Two time points: **3 days** and **20 days post-PCI**  

Raw count data were downloaded separately, and sample metadata were obtained using the GEOquery package.

---

## What this pipeline does

This workflow focuses strictly on **pre-DE exploratory analysis**, which is a critical but often under-documented step in RNA-seq studies.

Specifically, the pipeline performs:

- Library size normalization using edgeR
- CPM-based filtering to remove lowly expressed genes
- Log-CPM transformation for visualization and exploratory analysis
- Metadata-aware subsetting of samples based on treatment and time point
- Sample–sample Pearson correlation analysis to assess global expression similarity
- Principal component analysis (PCA) to evaluate separation by treatment and time

The analysis compares:
- 3-day sham vs 3-day PCI samples  
- 20-day sham vs 20-day PCI samples  
- PCI samples across 3 days vs 20 days to examine temporal progression

No differential expression testing is performed in this repository; this project is intended to establish a clean and well-understood foundation before DEG and pathway analysis.

---



<h2>Repository structure</h2>

<pre>
Sepsis-Associated-Encephalopathy-SAE-temporal-transcriptomics/
├── data/
│   ├── raw/
│   └── processed/
├── scripts/
│   ├── 00_setup.R
│   ├── 01_download_and_preprocess.R
│   ├── 02_data_wrangling.R
│   ├── 03_sample_correlation.R
│   └── 04_pca_analysis.R
├── figures/
│   ├── correlation/
│   └── pca/
└── README.md
</pre>


## How to run the analysis

Each step is modular and can be run sequentially:

```r
source("scripts/01_download_and_preprocess.R")
source("scripts/02_data_wrangling.R")
source("scripts/03_sample_correlation.R")
source("scripts/04_pca_analysis.R")
````






