# Install packages
pacman::p_load("tidyverse", "DECIPHER", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the raw chemical data
dat <- read_csv("data/substrate_comparisons/15pNPs_159_selected_molecular_properties.csv") %>%
  dplyr::select(contains(c("log", "cmpnd_abbrev")))
head(dat)

# Combine with other desscriptors
subkey <- read_excel("shiny_app/data/pNP_substrate_key.xlsx")

# Combine all
xel <- subkey %>%
  inner_join(., dat, by = "cmpnd_abbrev") %>%
  dplyr::select(-ALogP_octanol_water_partition_coefficient)
head(xel)
write_csv(xel, "output/15_pNPs_chemical_properties.csv")
