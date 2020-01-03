# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", 
               "randomForest","gbm","nnet","e1071","svmpath","lars",
               "glmnet","svmpath", "data.table")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Look at original channel residues
channel_a <- c(253, 258, 261, 284, 291, 292, 295, 343, 345, 349, 351, 353)
channel_b <- c(176, 173, 172, 242, 243, 112, 111, 171, 117, 316, 203, 246)
writeLines(as.character(channels), sep = "+")
channels <- sort(c(channel_a, channel_b))

# 8 angstrom overlap
ang8 <- fread("bin/pymol_scripts/resn_8angstroms.txt", header = F, data.table = F) %>%
  dplyr::mutate(resi = as.numeric(gsub("\\'", "", stringr::word(start = 1, sep = "\\)", V2)))) %>%
  dplyr::distinct()
ang8

intersect(ang8$resi, channels) # only 6 overlapping
setdiff(ang8$resi, channels)


# 10 angstrom overlap
ang10 <- fread("bin/pymol_scripts/resn_10angstroms.txt", header = F, data.table = F) %>%
  dplyr::mutate(resi = as.numeric(gsub("\\'", "", stringr::word(start = 1, sep = "\\)", V2)))) %>%
  dplyr::distinct()
ang10

length(intersect(ang10$resi, channels)) # 13 overlapping, about half
setdiff(ang10$resi, channels)

# 10 angstrom overlap
ang10 <- fread("bin/pymol_scripts/resn_10angstroms.txt", header = F, data.table = F) %>%
  dplyr::mutate(resi = as.numeric(gsub("\\'", "", stringr::word(start = 1, sep = "\\)", V2)))) %>%
  dplyr::distinct()
ang10

length(intersect(ang10$resi, channels)) # 13 overlapping, about half
setdiff(ang10$resi, channels)

# 12 angstrom overlap
ang12 <- fread("bin/pymol_scripts/resn_12angstroms.txt", header = F, data.table = F) %>%
  dplyr::mutate(resi = as.numeric(gsub("\\'", "", stringr::word(start = 1, sep = "\\)", V2)))) %>%
  dplyr::distinct()
ang12$resi
nrow(ang12) # 84 residues total
length(intersect(ang12$resi, channels)) # 19 overlapping
write_csv(data.frame(ang12$resi, stringsAsFactors = F), "data/machine_learning/84_residues_12_angstroms_4KU5_S143.csv")

# 14 angstrom overlap
ang14 <- fread("bin/pymol_scripts/resn_14angstroms.txt", header = F, data.table = F) %>%
  dplyr::mutate(resi = as.numeric(gsub("\\'", "", stringr::word(start = 1, sep = "\\)", V2)))) %>%
  dplyr::distinct(resi)
ang14
write_csv(data.frame(ang14$resi, stringsAsFactors = F), "data/machine_learning/358_residues_14_angstroms_4KU5_S143.csv")


length(intersect(ang14$resi, channels)) # all 24 overlapping
nrow(ang14) # 358 residues
