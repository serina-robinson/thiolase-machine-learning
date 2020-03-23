# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret", "colorspace",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "RColorBrewer")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the principal componenets of molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

# Read in the sequence features 
seq_fts <- read_csv("data/machine_learning/73_12angstrom_4aa_features.csv") %>%
  dplyr::mutate(raw_org = word(enzyme, sep = "\\.1", 2)) %>%
  dplyr::mutate(org = paste0(word(raw_org, sep = "_", 2), " ", word(raw_org, sep = "_", 3))) %>%
  dplyr::select(-raw_org) # remove enzyme
seq_fts$org[seq_fts$enzyme == "4KU5_Xanthomonas_campestris"] <- "Xanthomonas campestris"
seq_fts$enzyme[seq_fts$enzyme == "4KU5_Xanthomonas_campestris"] <- "NP_635607.1_4KU5_Xanthomonas_campestris"
seq_fts <- seq_fts %>%
  dplyr::mutate(acc = word(enzyme, sep = "\\.1", 1))
