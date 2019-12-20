# Preparation of the dataset for machine learning

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Install packages
pacman::p_load("tidyverse", "readxl", "ggplot2", "RColorBrewer")

# Read in the activity data
activity <- read_csv("output/substrate_comparisons/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Read in the chemical PC data

