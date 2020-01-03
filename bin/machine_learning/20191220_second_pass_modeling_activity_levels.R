# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

# Read in the sequence features 
seq_fts <- read_csv("data/machine_learning/73_channels_a_b_features.csv") %>%
  dplyr::mutate(raw_org = word(enzyme, sep = "\\.1", 2)) %>%
  dplyr::mutate(org = paste0(word(raw_org, sep = "_", 2), " ", word(raw_org, sep = "_", 3))) %>%
  dplyr::select(-raw_org)
seq_fts$org[seq_fts$enzyme == "4KU5_Xanthomonas_campestris"] <- "Xanthomonas campestris"

# Read in the activity data
activity <- read_csv("data/machine_learning/20191228_enzyme_activity_level.csv") %>%
  dplyr::rename(substrate = cmpnd) %>%
  dplyr::rename(activity = log_slope) %>%
  dplyr::mutate(org = paste0(word(org, sep = " ", 1), " ", word(org, sep = " ", 2))) %>%
  dplyr::select(org, substrate, activity, activity_level)
head(activity)

# Fix discrepancies in merging names
activity$org[grep("XC", activity$org)] <- "Xanthomonas campestris"
pseudo1 <- seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)]
pseudo2 <- activity$org[grep("Pseudoxanthomonas", activity$org)][1]
seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)] <- pseudo2

leif1 <- seq_fts$org[grep("Leifsonia", seq_fts$org)]
leif2 <- activity$org[grep("Leifsonia", activity$org)][1]
seq_fts$org[grep("Leifsonia", seq_fts$org)] <- leif2

activity$org[!activity$org %in% seq_fts$org]
activity$substrate <- gsub(" ", "\\.", activity$substrate)
activity$substrate %in% molec_fts$substrate

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., seq_fts, by = "org")
comb[is.na(comb$PC2),]

# Now remove duplicate rows (hopefully there aren't any)
dedup <- comb[complete.cases(comb),] # no duplicates
# dedupr <- comb[duplicated(comb$id),]

write_csv(dedup, "data/machine_learning/20191228_1095_training_examples_pass2.csv")
dedup$activity_level
head(cbind(dedup$activity_level, dedup$activity))

