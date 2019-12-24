# Install packages
pacman::p_load("tidyverse", "DECIPHER", "Biostrings")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the reference
sqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})


sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
source("lib/extract_channel_aas.r")

# extract_34_list <- lapply(1:length(sqs), function(x) { extract_channel_aas(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(extract_34_list), byrow=T), 
                            stringsAsFactors=FALSE)


ftnams <- c("WOLS870101", "WOLS870102", "WOLS870103", "FAUJ880109", "GRAR740102", "RADA880108",
         "ZIMJ680103", "TSAJ990101", "CHOP780201", "CHOP780202", "CHOP780203", "ZIMJ680104",
         "NEU1", "NEU2", "NEU3")
num_fts <- ncol(extract_34_df)/length(ftnams)
clnams <- paste0(ftnams, "_", rep(1:num_fts, each = length(ftnams)))
colnames(extract_34_df)[1:ncol(extract_34_df)] <- clnams
clnams
rownames(extract_34_df) <- names(sqs)
towrite <- data.frame(cbind(rownames(extract_34_df), extract_34_df), stringsAsFactors = F)
colnames(towrite)[1] <- "enzyme"
write_csv(towrite, "data/machine_learning/73_channels_a_b_features.csv")
