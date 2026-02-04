#--------------------------------------
#  Creating directory structure to save files
#------------------------------------------

library(here)

# Project root 
PROJECT_ROOT <- here::here()

# Define all paths RELATIVE to root
DIR_RAW    <- here("data", "raw")
DIR_PROCESSED <- here("data", "processed")
DIR_SCRIPTS <- here("scripts")
DIR_RESULTS   <- here("results")
DIR_FIGURES  <- here("figures")

# Create directories if missing
dirs <- c(DIR_RAW, DIR_PROCESSED,DIR_SCRIPTS, DIR_RESULTS, DIR_FIGURES)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Optional sanity check
message("Project root: ", PROJECT_ROOT)

#----------------------------------------
#            setup
#----------------------------------------
library(GEOquery)
library(readr)
library(data.table)
library(edgeR)
library(pheatmap)
library(ggrepel)