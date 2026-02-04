source(file.path(DIR_SCRIPTS,"00_setup.R"))
#-------------------------------------
#             DATA WRANGLING
#---------------------------------------

# Data wrangle according to experimental design
# Sham --- NO PCI treatement ( PCI = peritoneal contamination and infection)
# PCI --- PCI treatment
# Metadata for 3 days of sham and PCI

#logCPM for visualization
logCPM_counts <- readRDS(file.path(DIR_PROCESSED, "logCPM_counts.rds"))

metadata <- readRDS(file.path(DIR_PROCESSED,"metadata_clean.rds"))

#------------------------------------------
#      DAY-3 PCI VS SHAM
#------------------------------------------

# subsetting day-3 metadata from the whole metadata
Day_3_all_metadata <- metadata[
  metadata$characteristics_ch1.1 %in% c("treatment: sham", "treatment: PCI") &
    metadata$characteristics_ch1.2 == "time: 3 days",
]

# subsetting day-3 logcpm counts from the log transformed counts table 
Day_3_all_counts <- logCPM_counts[, Day_3_all_metadata$sample_id, drop = FALSE]

#------------------------------------------
#      DAY-20 PCI VS SHAM
#------------------------------------------

# subsetting day-20 metadata from the whole metadata
Day_20_all_metadata <- metadata[
  metadata$characteristics_ch1.1 %in% c("treatment: sham", "treatment: PCI") &
    metadata$characteristics_ch1.2 == "time: 20 days",
]
# subsetting day-3 logcpm counts from the log transformed counts table 
Day_20_all_counts <- logCPM_counts[, Day_20_all_metadata$sample_id, drop = FALSE]

#------------------------------------------
#      DAY-3 PCI VS DAY-20 PCI
#------------------------------------------

D3_D20_meta <- metadata[metadata$characteristics_ch1.1 == "treatment: PCI", ]
D3_D20_counts <- logCPM_counts[, D3_D20_meta$sample_id]

# saving the separated counts and metadata for visualization
saveRDS(list(
  d3 = list(meta = Day_3_all_metadata, counts = Day_3_all_counts),
  d20 = list(meta = Day_20_all_metadata, counts = Day_20_all_counts),
  pci = list(meta = D3_D20_meta, counts = D3_D20_counts)
), "data/processed/subsets.rds")
