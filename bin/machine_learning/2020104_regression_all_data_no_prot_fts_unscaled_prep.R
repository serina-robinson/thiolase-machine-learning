# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "skimr", "caret",
               "cowplot", "tidymodels", "ranger", "tree", "rsample", "randomForest")

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

# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Fix discrepancies in merging names
pseudo1 <- seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)]
pseudo2 <- activity$org[grep("Pseudoxanthomonas", activity$org)][1]
seq_fts$org[grep("Pseudoxanthomonas", seq_fts$org)] <- pseudo2

leif1 <- seq_fts$org[grep("Leifsonia", seq_fts$org)]
leif2 <- activity$org[grep("Leifsonia", activity$org)][1]
seq_fts$org[grep("Leifsonia", seq_fts$org)] <- leif2

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., seq_fts, by = "org") 

# Now remove duplicate rows (hopefully there aren't any)
dedup <- comb[complete.cases(comb),] # no duplicates
dedup <- dedup[!duplicated(dedup),]
tow <- colnames(dedup)[grepl("4KU5pos_", colnames(dedup))] 

# Only keep variables with nonzero variance
nozdat <- nearZeroVar(dedup, saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]
# write_csv(data.frame(which_rem, stringsAsFactors = F), "shiny_app/data/vars_to_remove.csv")

dat <- dedup %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::filter(is_active == "Y") %>%
  dplyr::select(-which_rem, -org, -substrate, -enzyme, -is_active, -acc)


# Write unscaled dataset to file
write_csv(dat, "data/machine_learning/20200323_550_regression_training_examples_no_prot_feats_unscaled.csv")
