# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings", "fastDummies", "caret")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/thiolase-machine-learning/")

# Read in the reference
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})

sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
source("lib/extract_12angstrom_residues_one_hot.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_12angstrom_one_hot(query_fils[x]) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                            stringsAsFactors=FALSE)

# Write unencoded data frame to file
onehot <- fastDummies::dummy_cols(extract_84_df, remove_selected_columns = T, ignore_na = T)
onehot_df <- tibble(names(sqs)) %>%
  bind_cols(onehot)
colnames(onehot_df)[1] <- "enzyme"
onehot_fts <- onehot_df %>%
  dplyr::mutate(raw_org = word(enzyme, sep = "\\.1", 2)) %>%
  dplyr::mutate(org = paste0(word(raw_org, sep = "_", 2), " ", word(raw_org, sep = "_", 3))) %>%
  dplyr::mutate(acc = word(enzyme, sep = "\\.1", 1)) %>%
  dplyr::select(-raw_org) # remove enzyme
onehot_fts$org[1] <- "Xanthomonas campestris"

# Merge with the rest of the data
# Read in the principal componenets of molecular features
molec_fts <- read_csv("data/machine_learning/PC7_molecular_descriptors.csv") %>%
  dplyr::rename(substrate = V1)
molec_fts$substrate <- gsub(" ", "\\.", molec_fts$substrate)
molec_fts$substrate[molec_fts$substrate == "oxidazole"] <- "oxadiazole"

# Read in the activity data
activity <- read_csv("data/machine_learning/20191218_all_cmpnds_avg_log_slopes_for_modeling.csv")

# Fix discrepancies in merging names
pseudo1 <- onehot_fts$org[grep("Pseudoxanthomonas", onehot_fts$org)]
pseudo2 <- activity$org[grep("Pseudoxanthomonas", activity$org)][1]
onehot_fts$org[grep("Pseudoxanthomonas", seq_fts$org)] <- pseudo2

leif1 <- onehot_fts$org[grep("Leifsonia", onehot_fts$org)]
leif2 <- activity$org[grep("Leifsonia", activity$org)][1]
onehot_fts$org[grep("Leifsonia", onehot_fts$org)] <- leif2

# Read in the protein features
prot_fts <- read_csv("data/machine_learning/73_overall_calculated_protein_properties.csv") 

# Now merge everything...
comb <- activity %>%
  dplyr::left_join(., molec_fts, by = "substrate") %>%
  dplyr::left_join(., onehot_fts, by = "org") %>%
  dplyr::left_join(., prot_fts, by = "acc")

# Now remove duplicate rows
dedup <- comb[complete.cases(comb),] # no duplicates
dedup <- dedup[!duplicated(dedup),]

# write_csv(dedup, "data/machine_learning/20200111_1095_training_examples_12angstrom_features.csv")

# Create ID variable
dat <- dedup %>%
  dplyr::mutate(id = paste0(org, "_", substrate)) %>%
  dplyr::mutate(is_active = as.factor(case_when(activity == 0 ~ "N",
                                                TRUE ~ "Y"))) %>%
  dplyr::select(-org, -substrate, -enzyme, -activity, -sqs, -core, -acc, -nams) #-IUPAC, -SMILES, -cmpnd_abbrev, ) 

# Identify columns that have "NaNs"
cols <- colSums(mapply('==', 'NaN', dat))
new.df <- dat[,which(cols == 0)]

write_csv(new.df, "data/machine_learning/20200322_1095_training_examples_12angstrom_one_hot_encoded_unscaled.csv")
