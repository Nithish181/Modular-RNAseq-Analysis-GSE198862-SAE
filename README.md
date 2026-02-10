# Sepsis-Associated-Encephalopathy (SAE): Temporal Transcriptomic Analysis  
**GEO dataset:** GSE198862

This repository contains a modular RNA-seq analysis workflow for studying temporal gene-expression changes in a murine model of sepsis-associated encephalopathy (SAE). The analysis is based on the publicly available GEO dataset **GSE198862**, which profiles brain transcriptomic responses following peritoneal contamination and infection (PCI).

The project is structured to progress from exploratory data analysis (EDA) and quality control to differential expression analysis (DEG), with an ongoing comparison of DESeq2 and edgeR at biologically relevant contrasts.

---

## Dataset

- **GEO accession:** GSE198862  
- **Organism:** *Mus musculus*  
- **Experimental design:**  
  - Sham vs PCI treatment  
  - Two time points: **3 days** and **20 days post-PCI**

Raw count data were downloaded separately, and sample metadata were retrieved using the GEOquery package.

---

## What this project does

This workflow emphasizes a **stepwise and biologically driven RNA-seq analysis**, beginning with exploratory analysis to validate data quality and experimental structure, followed by contrast-specific differential expression analysis.

### Exploratory data analysis (EDA)

The EDA component focuses on pre-DE validation and includes:

- Library size normalization using edgeR
- CPM-based filtering to remove lowly expressed genes
- Log-CPM transformation for visualization
- Metadata-aware subsetting by treatment and time point
- Sample–sample Pearson correlation analysis
- Principal component analysis (PCA) to assess global expression structure

Exploratory analyses compare:
- 3-day sham vs 3-day PCI samples  
- 20-day sham vs 20-day PCI samples  
- PCI samples across 3 days vs 20 days to assess temporal progression

---

### Differential expression analysis and benchmarking (in progress)

Following EDA, differential expression analysis is implemented as **contrast-specific case studies** using both **DESeq2 and edgeR**. Each method is applied independently to the same biological comparisons, enabling direct inspection of consistency and method-specific differences in detected signals.

Current DEG analyses focus on:
- Day-3 sham vs PCI  
- Day-20 sham vs PCI  
- Temporal comparison of PCI samples (Day-3 vs Day-20)

At this stage, the comparison between DESeq2 and edgeR is **qualitative and contrast-driven**, intended to validate biological trends and model behavior rather than to serve as a large-scale benchmarking study. Refactoring toward a generalized and reusable DEG framework is planned after these case studies are stabilized.

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
│   ├── 04_pca_analysis.R
│   ├── 05_DESeq2_DEG_day3_sham_vs_PCI.R
│   ├── 06_DESeq2_DEG_day20_sham_vs_PCI.R
│   ├── 07_DESeq2_DEG_PCI_day3_vs_day20.R # To study disease progression and it is ongoing
│   ├── 08_edgeR_DEG_day3_sham_vs_PCI.R
    ├── 09_edgeR_DEG_day20_sham_vs_PCI.R
    ├── 10_edgeR_DEG_PCI_day3_vs_day20.R  # To study disease progression and it is ongoing
│   └── 11_DEG_template_generalized.R   # planned
├── figures/
│   ├── correlation/
│   └── pca/
└── README.md
</pre>

---

## How to run the analysis

Each step is modular and can be run sequentially:

```r
source("scripts/01_download_and_preprocess.R")
source("scripts/02_data_wrangling.R")
source("scripts/03_sample_correlation.R")
source("scripts/04_pca_analysis.R")

