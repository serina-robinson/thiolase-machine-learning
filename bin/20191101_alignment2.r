# Pipeline for files with # Install packages
pacman::p_load("tidyverse", "readxl", "randomcoloR", "RColorBrewer", "ggplot2", "broom", "maditr")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the data
dat <- readAAStringSet("~/Documents/University_of_Minnesota/Wackett_Lab/github/mBIO_paper_figs/data/73_OleA_aligned_trimmed.fasta")
kyto <- dat[grep("Kytococcus|Xanthomonas", names(dat))]
BrowseSeqs(AAStringSet(kyto))
