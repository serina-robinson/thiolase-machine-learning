# Install packages
pacman::p_load("tidyverse", "Biostrings", "DECIPHER", "readxl")

# Set working directory
setwd("~/Documents/University_of_Minnesota/Wackett_Lab/github/synbio-data-analysis/")

# Read in the 73 accession numbers
seqs <- readAAStringSet("data/alignments/73_OleA_JGI_unaligned.fasta")
acc <- word(names(seqs), sep = "\\.1", 1)

# Read in the freezer stock key
sheet1 <- read_excel("data/plasmid_sequences/JGI_freezer_stock_KEY.xlsx", sheet = 1) %>%
  janitor::clean_names()
sheet1_olea <- sheet1[grep(paste0(acc, collapse = "|"), sheet1$users_name),]
sheet1_olea$sequence
sheet1_seqs <- AAStringSet(sheet1_olea$sequence)
names(sheet1_seqs) <- sheet1_olea$users_name


sheet2 <- read_excel("data/plasmid_sequences/JGI_freezer_stock_KEY.xlsx", sheet = 2) %>%
  janitor::clean_names()
sheet2_olea <- sheet2[grep(paste0(acc, collapse = "|"), sheet2$user_id),]
sheet2_olea$sequence
sheet2_seqs <- AAStringSet(sheet2_olea$sequence)
names(sheet2_seqs) <- sheet2_olea$user_id

# Combine all sequences
allsqs <- c(sheet1_seqs, sheet2_seqs)
finalsqs <- AAStringSet(allsqs)

# Pull in the OleA sequences
olea <- readAAStringSet("data/plasmid_sequences/Xc_OleA_WT_plasmid_sequence.fa")

# Combine with 72 seqs
combsqs <- AAStringSet(c(olea, finalsqs))
combsqs
writeXStringSet(combsqs, "data/plasmid_sequences/73_oleA_plasmid_sequences.fa")


