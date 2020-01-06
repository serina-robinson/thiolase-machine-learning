# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the reference
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})


sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
source("lib/extract_14angstrom_aas.R")

extract_84_list <- lapply(1:length(sqs), function(x) { extract_14angstrom_aas(query_fils[x]) })
extract_84_df <- data.frame(matrix(unlist(extract_84_list), nrow = length(extract_84_list), byrow=T), 
                            stringsAsFactors=FALSE)
head(extract_84_df)


ftnams <- colnames(extract_84_list[[1]])
colnames(extract_84_df) <- ftnams

rownames(extract_84_df) <- names(sqs)
towrite <- data.frame(cbind(rownames(extract_84_df), extract_84_df), stringsAsFactors = F)
colnames(towrite)[1] <- "enzyme"
colnames(extract_84_df)
write_csv(towrite, "data/machine_learning/73_14angstrom_4aa_features.csv")
